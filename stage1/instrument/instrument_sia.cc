// Stage 1 instrumentation driver for the PISM SIA stress-balance boundary.
//
// Uses the PISM verification test F exact solution (Bueler, Brown, Lingle,
// J. Glaciol. 53(182), 2007) as the input state, so every dumped quantity has
// an exact reference value. The flow law is `arr` (Paterson-Budd, cold case),
// no sliding (test F is a no-slip SIA solution).
//
// Dumps (instantaneous, one stress-balance update at t = 0):
//
//   inputs.nc              - all stress-balance inputs (geometry, enthalpy, age)
//   parameters.csv         - config keys affecting the SIA boundary
//   columns.csv            - per cell-center column (i,j,k):
//                              geometry, enthalpy, u3/v3/w3, strain heating,
//                              and the exact-solution reference (T,u,v,w,sigma)
//   staggered_D.csv        - per staggered point: h_x, h_y, |grad h|, D
//   diffusive_flux.csv     - per staggered point: D, grad h -> q = -D grad h
//   delta.csv              - per staggered point and level: z, depth, pressure,
//                              E, stress, flow (recomputed), delta integrand,
//                              D (the trapezoidal integral of delta)
//   basal_frictional_heat.csv - F_b (zero here: no sliding)
//   theta.csv              - Schoof bed-smoothness factor (1 for the flat bed)
//   exact_solution.csv     - exactFG reference at cell centers (inside the sheet)
//   meta.txt               - PISM revision, config file, command line, grid
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/sia.json
//
// The SIAFD internals that are not public (m_delta_0/1, m_e_factor) are
// exposed through the InstrumentedSIAFD subclass below.

#include <petsc.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "pism/util/Config.hh"
#include "pism/util/Context.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Grid.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/Units.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/sia/BedSmoother.hh"
#include "pism/stressbalance/sia/SIAFD.hh"
#include "pism/stressbalance/ShallowStressBalance.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM SIA stress-balance component.\n"
  "Dumps SIA I/O traces (test-F exact-solution input state) to "
  "stage1/dumps/sia/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

using pism::rheology::FlowLaw;
using pism::rheology::FlowLawFactory;
using pism::stressbalance::SIAFD;

//! Subclass exposing SIAFD internals needed for the dumps.
class InstrumentedSIAFD : public SIAFD {
public:
  InstrumentedSIAFD(std::shared_ptr<const Grid> g) : SIAFD(g) {}

  const array::Array3D &delta_0() const { return m_delta_0; }
  const array::Array3D &delta_1() const { return m_delta_1; }
  double e_factor() const { return m_e_factor; }
};

// ---------------------------------------------------------------------------
// Output helpers (same as instrument_rheology.cc)
// ---------------------------------------------------------------------------

std::ofstream open_output(const std::string &path) {
  std::ofstream out(path);
  if (not out) {
    throw std::runtime_error("Cannot open output file '" + path + "'");
  }
  out << std::setprecision(17);
  return out;
}

void header(std::ofstream &out, const std::string &key, const std::string &value) {
  out << "# " << key << " = " << value << "\n";
}

//! Dump config parameters affecting the SIA boundary.
void dump_parameters(const Config &config, const std::string &dir) {
  std::ofstream out = open_output(dir + "/parameters.csv");
  header(out, "description",
         "Configuration parameters affecting the SIA stress-balance component "
         "(stress_balance.sia.*, flow_law.*, constants.*, grid.*, enthalpy_converter.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {
      "stress_balance.sia.", "flow_law.", "constants.", "grid.",
      "enthalpy_converter.", "energy.model"};

  for (const auto &key : config.keys()) {
    bool match = false;
    for (const auto &prefix : prefixes) {
      if (key.rfind(prefix, 0) == 0) {
        match = true;
        break;
      }
    }
    if (not match) {
      continue;
    }
    if (key.size() > 6 and
        (key.compare(key.size() - 4, 4, "_doc") == 0 or
         key.compare(key.size() - 5, 5, "_type") == 0 or
         key.compare(key.size() - 6, 6, "_units") == 0 or
         key.compare(key.size() - 4, 4, "_opt") == 0 or
         key.compare(key.size() - 7, 7, "_option") == 0)) {
      continue;
    }
    try {
      const std::string type = config.type(key);
      if (type == "number") {
        out << key << "," << config.get_number(key) << " " << config.units(key) << "\n";
      } else if (type == "boolean") {
        out << key << "," << (config.get_flag(key) ? "true" : "false") << "\n";
      } else {
        out << key << "," << config.get_string(key) << "\n";
      }
    } catch (const std::exception &) {
      out << key << ",<unreadable>\n";
    }
  }
  out.close();
}

} // namespace
} // namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::stressbalance;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "sia_test");
    auto config = ctx->config();

    // Deterministic, age-decoupled config (e_factor = flow_law.enhancement_factor)
    config->set_flag("stress_balance.sia.grain_size_age_coupling", false);
    config->set_flag("stress_balance.sia.e_age_coupling", false);
    config->set_string("stress_balance.sia.flow_law", "arr");

    set_config_from_options(*config);
    config->resolve_filenames();

    // output directory
    std::string out_dir = "stage1/dumps/sia";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    // grid: test-F domain
    grid::Parameters P(*config);
    P.Lx = 900.0e3;
    P.Ly = P.Lx;
    const double Lz = 4000.0;
    const unsigned int Mz = config->get_number("grid.Mz");
    P.z = grid::compute_vertical_levels(Lz, Mz, grid::EQUAL);
    P.ownership_ranges_from_options(*config, ctx->size());

    auto grid = std::make_shared<Grid>(ctx, P);

    // record run metadata
    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "test", "F (Bueler-Brown-Lingle exact solution, t=0)");
      header(meta, "flow_law", "arr (Paterson-Budd cold)");
      header(meta, "sliding", "none (ZeroSliding)");
      header(meta, "Mx", std::to_string(grid->Mx()));
      header(meta, "My", std::to_string(grid->My()));
      header(meta, "Mz", std::to_string(grid->Mz()));
      header(meta, "Lx", std::to_string(P.Lx));
      header(meta, "Lz", std::to_string(Lz));
      header(meta, "dx", std::to_string(grid->dx()));
      header(meta, "command", "");
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    std::shared_ptr<EnthalpyConverter> EC(new ColdEnthalpyConverter(*config));

    const int WIDE_STENCIL = config->get_number("grid.max_stencil_width");
    array::Array3D enthalpy(grid, "enthalpy", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL),
        age(grid, "age", array::WITHOUT_GHOSTS, grid->z());
    age.set(0.0);

    Geometry geometry(grid);
    geometry.sea_level_elevation.set(0.0);

    // fill the exact test-F state
    {
      const double LforFG = 750000.0;
      geometry.bed_elevation.set(0.0);
      geometry.cell_type.set(MASK_GROUNDED);
      {
        array::AccessScope list{&geometry.ice_thickness};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          const double r = std::max(grid::radius(*grid, i, j), 1.0);
          geometry.ice_thickness(i, j) =
              (r > LforFG - 1.0) ? 0.0 : exactFG(0.0, r, grid->z(), 0.0).H;
        }
      }
      geometry.ice_thickness.update_ghosts();
      geometry.ice_surface_elevation.copy_from(geometry.ice_thickness);
      geometry.ice_surface_elevation.update_ghosts();

      // enthalpy from the exact temperature profile
      {
        array::AccessScope list{&enthalpy, &geometry.ice_thickness};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          const double r = std::max(grid::radius(*grid, i, j), 1.0);
          const double H = geometry.ice_thickness(i, j);
          double *E_ij = enthalpy.get_column(i, j);
          if (r > LforFG - 1.0) {
            for (unsigned int k = 0; k < grid->Mz(); ++k) {
              E_ij[k] = EC->enthalpy(223.15, 0.0, EC->pressure(1000.0));
            }
          } else {
            TestFGParameters F = exactFG(0.0, r, grid->z(), 0.0);
            for (unsigned int k = 0; k < grid->Mz(); ++k) {
              const double depth = H - grid->z(k);
              E_ij[k] = EC->enthalpy_permissive(F.T[k], 0.0, EC->pressure(depth));
            }
          }
        }
        enthalpy.update_ghosts();
      }
    }

    geometry.ensure_consistency(config->get_number("geometry.ice_free_thickness_standard"));

    // SIA solver (no sliding)
    std::shared_ptr<InstrumentedSIAFD> sia(new InstrumentedSIAFD(grid));
    std::shared_ptr<ZeroSliding> no_sliding(new ZeroSliding(grid));
    StressBalance stress_balance(grid, no_sliding, sia);
    stress_balance.init();

    stressbalance::Inputs inputs;
    inputs.geometry = &geometry;
    inputs.water_column_pressure = nullptr;
    inputs.enthalpy = &enthalpy;
    inputs.age = &age;

    stress_balance.update(inputs, true);

    // ------------------------------------------------------------------
    // Dumps
    // ------------------------------------------------------------------
    const array::Array3D &u3 = stress_balance.velocity_u(),
                         &v3 = stress_balance.velocity_v(),
                         &w3 = stress_balance.velocity_w(),
                         &sigma = stress_balance.volumetric_strain_heating();
    const array::Staggered &h_x = sia->surface_gradient_x(),
                           &h_y = sia->surface_gradient_y(),
                           &D = sia->diffusivity(),
                           &q = stress_balance.diffusive_flux();
    const array::Scalar &F_b = stress_balance.basal_frictional_heating();
    const array::Array3D &delta0 = sia->delta_0(), &delta1 = sia->delta_1();

    const double e_factor = sia->e_factor();
    const double theta_min = config->get_number("stress_balance.sia.bed_smoother.theta_min");

    // recompute theta with the same bed smoother for the delta chain dump
    array::Scalar2 theta(grid, "theta");
    sia->bed_smoother().theta(geometry.ice_surface_elevation, theta);

    // flow law for recomputing the flow factor in the delta chain
    double n = config->get_number("stress_balance.sia.Glen_exponent");
    FlowLawFactory factory(config, EC);
    auto flow_law = factory.create("arr", n);
    const double grain_size = config->get_number("constants.ice.grain_size", "m");

    // exact-solution scaling (Sig is K/s -> J/(s m^3))
    const double ice_rho = config->get_number("constants.ice.density"),
                 ice_c = config->get_number("constants.ice.specific_heat_capacity"),
                 LforFG = 750000.0;

    {
      // all inputs + outputs as NetCDF (mirrors the PISM input dump pattern;
      // note: PISM's own StressBalance::Inputs::dump() is broken -- it calls
      // io::write_config() on an undefined variable, see bugs.md)
      {
        auto writer = std::make_shared<SynchronousOutputWriter>(grid->com, *config);
        writer->initialize({}, true);
        OutputFile file(writer, out_dir + "/inputs.nc");

        const array::Array *vecs[] = {
            &geometry.ice_surface_elevation, &geometry.ice_thickness,
            &geometry.bed_elevation,         &geometry.cell_type,
            &enthalpy,                       &age,
            &u3,                             &v3,
            &w3,                             &sigma,
            &D,                              &h_x,
            &h_y,                            &q,
            &F_b};

        auto time = ctx->time();
        file.define_variable(time->metadata());
        for (const auto *vec : vecs) {
          for (auto &var : vec->all_metadata()) {
            file.define_variable(var);
          }
        }
        file.append_time(time->current());
        for (const auto *vec : vecs) {
          vec->write(file);
        }
      }

      // cell-center columns: geometry + velocity + heating + exact reference
      std::ofstream out = open_output(out_dir + "/columns.csv");
      header(out, "description",
             "SIA output columns at cell centers (test-F exact state): "
             "per (i,j,k): geometry, enthalpy, computed u3/v3/w3/sigma, "
             "and exact reference T/u/v/w/sigma");
      header(out, "columns",
             "i,j,x,y,r,H,surface,bed,mask,k,z,E,T_exact,u3,v3,w3,sigma3,"
             "sigma_exact_Jpm3,u_exact,v_exact,w_exact");

      array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                              &geometry.bed_elevation, &geometry.cell_type, &enthalpy,
                              &u3, &v3, &w3, &sigma};

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double x = grid->x(i), y = grid->y(j);
        const double r = std::max(grid::radius(*grid, i, j), 1.0);

        const double H = geometry.ice_thickness(i, j);

        const double *E_ij = enthalpy.get_column(i, j);
        const double *u_ij = u3.get_column(i, j), *v_ij = v3.get_column(i, j),
                     *w_ij = w3.get_column(i, j), *sig_ij = sigma.get_column(i, j);

        // exact reference exists only inside the sheet (0 < r < L)
        const bool have_exact = (r <= LforFG - 1.0);
        TestFGParameters F(grid->Mz());
        if (have_exact) {
          F = exactFG(0.0, r, grid->z(), 0.0);
        }
        const double inv_r = (r > 1.0) ? (1.0 / r) : 0.0;

        for (unsigned int k = 0; k < grid->Mz(); ++k) {
          const double T_x = have_exact ? F.T[k] : std::numeric_limits<double>::quiet_NaN();
          const double sig_x =
              have_exact ? F.Sig[k] * ice_rho * ice_c : std::numeric_limits<double>::quiet_NaN();
          const double u_x =
              have_exact ? inv_r * x * F.U[k] : std::numeric_limits<double>::quiet_NaN();
          const double v_x =
              have_exact ? inv_r * y * F.U[k] : std::numeric_limits<double>::quiet_NaN();
          const double w_x =
              have_exact ? F.w[k] : std::numeric_limits<double>::quiet_NaN();
          out << i << "," << j << "," << x << "," << y << "," << r << "," << H << ","
              << geometry.ice_surface_elevation(i, j) << "," << geometry.bed_elevation(i, j)
              << "," << (int)geometry.cell_type(i, j) << "," << k << "," << grid->z(k) << ","
              << E_ij[k] << "," << T_x << "," << u_ij[k] << "," << v_ij[k] << ","
              << w_ij[k] << "," << sig_ij[k] << "," << sig_x << "," << u_x << "," << v_x
              << "," << w_x << "\n";
        }
      }
      out.close();
    }

    // staggered grid: h_x, h_y, alpha, D
    {
      std::ofstream out = open_output(out_dir + "/staggered_D.csv");
      header(out, "description",
             "SIA staggered points: surface gradient (h_x,h_y), |grad h|, "
             "and the flow-law diffusivity D (m^2 s^-1)");
      header(out, "columns", "i,j,o,x_stag,y_stag,h_x,h_y,alpha,D");

      array::AccessScope list{&h_x, &h_y, &D};

      for (int o = 0; o < 2; ++o) {
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          if (i >= (int)grid->Mx() - 1 or j >= (int)grid->My() - 1) {
            continue;
          }
          const double dx = 0.5 * grid->dx(), dy = 0.5 * grid->dy();
          const double x_stag = grid->x(i) + (o == 0 ? dx : 0.0);
          const double y_stag = grid->y(j) + (o == 1 ? dy : 0.0);
          const double alpha = std::sqrt(h_x(i, j, o) * h_x(i, j, o) +
                                         h_y(i, j, o) * h_y(i, j, o));
          out << i << "," << j << "," << o << "," << x_stag << "," << y_stag << ","
              << h_x(i, j, o) << "," << h_y(i, j, o) << "," << alpha << "," << D(i, j, o)
              << "\n";
        }
      }
      out.close();
    }

    // diffusive flux q = -D grad h
    {
      std::ofstream out = open_output(out_dir + "/diffusive_flux.csv");
      header(out, "description",
             "SIA diffusive flux at staggered points: q = -D grad h (m^2 s^-1)");
      header(out, "columns", "i,j,o,x_stag,y_stag,D,h_x,h_y,q_x,q_y");

      array::AccessScope list{&h_x, &h_y, &D, &q};

      for (int o = 0; o < 2; ++o) {
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          if (i >= (int)grid->Mx() - 1 or j >= (int)grid->My() - 1) {
            continue;
          }
          const double dx = 0.5 * grid->dx(), dy = 0.5 * grid->dy();
          const double x_stag = grid->x(i) + (o == 0 ? dx : 0.0);
          const double y_stag = grid->y(j) + (o == 1 ? dy : 0.0);
          out << i << "," << j << "," << o << "," << x_stag << "," << y_stag << ","
              << D(i, j, o) << "," << h_x(i, j, o) << "," << h_y(i, j, o) << ","
              << q(i, j, o) << "," << q(i, j, 1 - o) << "\n";
        }
      }
      out.close();
    }

    // delta chain: (z, depth, pressure, E, stress, flow) -> delta -> D
    {
      std::ofstream out = open_output(out_dir + "/delta.csv");
      header(out, "description",
             "SIA diffusivity chain at staggered points: per level k: z, depth, "
             "pressure p(z), enthalpy E, stress alpha*p, flow-law strain factor, "
             "delta = e*theta*2*p*flow; and D = trapezoid(delta) at the point");
      header(out, "columns",
             "i,j,o,k,z,depth,pressure,E,stress,flow,delta,D,theta");

      array::AccessScope list{&h_x, &h_y, &D, &delta0, &delta1, &enthalpy, &theta,
                              &geometry.ice_thickness, &geometry.ice_surface_elevation};

      const std::vector<double> &z = grid->z();
      const unsigned int Mz = grid->Mz();
      std::vector<double> depth(Mz), pressure(Mz), stress(Mz), Ecol(Mz);

      for (int o = 0; o < 2; ++o) {
        const array::Array3D &delta = (o == 0) ? delta0 : delta1;
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          if (i >= (int)grid->Mx() - 1 or j >= (int)grid->My() - 1) {
            continue;
          }
          const int oi = 1 - o, oj = o;

          // thickness at the staggered point (from the smoothed thickness used
          // by SIAFD; here geometry.ice_thickness is unsmoothed and the bed is
          // flat, so smoothed == raw)
          const double thk = 0.5 * (geometry.ice_thickness(i, j) +
                                    geometry.ice_thickness(i + oi, j + oj));
          if (thk == 0.0) {
            continue;
          }
          const int ks = grid->kBelowHeight(thk);

          for (unsigned int k = 0; k <= (unsigned int)ks; ++k) {
            depth[k] = thk - z[k];
          }
          EC->pressure(depth, ks, pressure);

          const double *E_ij = enthalpy.get_column(i, j),
                       *E_off = enthalpy.get_column(i + oi, j + oj);
          for (unsigned int k = 0; k <= (unsigned int)ks; ++k) {
            Ecol[k] = 0.5 * (E_ij[k] + E_off[k]);
          }

          const double alpha = std::sqrt(h_x(i, j, o) * h_x(i, j, o) +
                                         h_y(i, j, o) * h_y(i, j, o));
          for (unsigned int k = 0; k <= (unsigned int)ks; ++k) {
            stress[k] = alpha * pressure[k];
          }

          // recompute the flow factor with the same flow law (deterministic)
          std::vector<double> flow(ks + 1);
          flow_law->flow_n(stress.data(), Ecol.data(), pressure.data(),
                           &grain_size, ks + 1, flow.data());

          const double theta_local = 0.5 * (theta(i, j) + theta(i + oi, j + oj));

          for (unsigned int k = 0; k <= (unsigned int)ks; ++k) {
            out << i << "," << j << "," << o << "," << k << "," << z[k] << ","
                << depth[k] << "," << pressure[k] << "," << Ecol[k] << "," << stress[k]
                << "," << flow[k] << "," << delta.get_column(i, j)[k] << "," << D(i, j, o) << ","
                << theta_local << "\n";
          }
        }
      }
      out.close();
    }

    // basal frictional heating (zero: no sliding) and theta map
    {
      std::ofstream out = open_output(out_dir + "/basal_frictional_heat.csv");
      header(out, "description", "Basal frictional heating F_b (W m^-2); zero for no-sliding");
      header(out, "columns", "i,j,x,y,F_b");
      array::AccessScope list{&F_b};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << F_b(i, j)
            << "\n";
      }
      out.close();

      std::ofstream tout = open_output(out_dir + "/theta.csv");
      header(tout, "description",
             "Schoof bed-smoothness factor theta at cell centers (1 for the flat bed)");
      header(tout, "columns", "i,j,x,y,theta");
      array::AccessScope tlist{&theta};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        tout << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << ","
             << theta(i, j) << "\n";
      }
      tout.close();
    }

    // exact solution reference at cell centers (inside the sheet)
    {
      std::ofstream out = open_output(out_dir + "/exact_solution.csv");
      header(out, "description",
             "Exact test-F solution at cell centers (inside the sheet): H, and "
             "per level T, radial U, w, Sig (K s^-1)");
      header(out, "columns", "i,j,x,y,r,H,k,z,T,U,w,Sig");
      array::AccessScope list{&geometry.ice_thickness};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        if (r > LforFG - 1.0) {
          continue;
        }
        const double H = geometry.ice_thickness(i, j);
        if (H <= 0.0) {
          continue;
        }
        TestFGParameters F = exactFG(0.0, r, grid->z(), 0.0);
        for (unsigned int k = 0; k < grid->Mz(); ++k) {
          out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
              << H << "," << k << "," << grid->z(k) << "," << F.T[k] << "," << F.U[k]
              << "," << F.w[k] << "," << F.Sig[k] << "\n";
        }
      }
      out.close();
    }

    // scalar diagnostics
    std::cout << "SIA instrumentation complete.\n"
              << "  max_diffusivity D_max = " << stress_balance.max_diffusivity()
              << " m^2 s^-1\n"
              << "  e_factor = " << e_factor << ", theta_min = " << theta_min << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
