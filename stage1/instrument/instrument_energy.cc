// Stage 1 instrumentation driver for the PISM ice-energy and ice-age
// boundaries (the first *integrated*, time-advancing boundaries instrumented).
//
// In the coupled model the energy model consumes the 3-D velocity field
// (u3,v3,w3) and volumetric strain heating produced by the stress balance,
// and the age model consumes (u3,v3,w3).  To make the traces fully
// self-contained and analytic, we use the PISM verification test F exact
// solution (Bueler, Brown, Lingle, J. Glaciol. 53(182), 2007) *directly* as
// the forcing:
//
//   u3 = U(r,z) x/r,  v3 = U(r,z) y/r,  w3 = w(r,z)      (exactFG)
//   Sigma = Sig(r,z) * rho * c          (exact strain heating, K/s -> W/m^3)
//   surface_temp = Tmin + ST*r          (exact surface temperature)
//   basal_heat_flux = 0.042 W/m^2       (the Ggeo used by the exact solution)
//   shelf_base_temp, surface liquid frac, till water thickness = 0 / 0
//   F_b = 0                             (test F has no sliding)
//
// The exact temperature field is the steady state of the continuous energy
// equation for this forcing, so one enthalpy step from it produces only a
// *small* change (the discrete-vs-continuous gap).  Where the exact
// temperature exceeds the pressure-melting temperature (the deep interior
// near the dome, T > T_pmp) the enthalpy model clamps to the temperate state
// and produces a small basal melt rate -- both are validation targets the
// .esm reimplementation must reproduce.
//
// The age model starts from age = 0 everywhere; after one step of dt the age
// in the deep interior must equal dt exactly (ice ages one second per second)
// while the age-0 surface boundary keeps the top layers younger.
//
// Dumps (one step of dt years, both models):
//
//   energy/{meta.txt, parameters.csv, inputs.nc, columns.csv,
//           basal_melt_rate.csv}
//   age/{meta.txt, parameters.csv, inputs.nc, columns.csv}
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/energy_age.json

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
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/Units.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Array3D.hh"
#include "pism/energy/EnergyModel.hh"
#include "pism/energy/EnthalpyModel.hh"
#include "pism/age/AgeModel.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM ice-energy (EnthalpyModel) and\n"
  "ice-age (AgeModel) components. Dumps one-step I/O traces driven by the\n"
  "test-F exact solution to stage1/dumps/energy/ and stage1/dumps/age/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

//! Subclass exposing the (protected) enthalpy state so we can set the initial
//! condition directly from the exact test-F temperature field.
class InstrumentedEnthalpyModel : public energy::EnthalpyModel {
public:
  InstrumentedEnthalpyModel(std::shared_ptr<const Grid> g)
      : EnthalpyModel(g, nullptr) {}

  void set_enthalpy(const array::Array3D &E) {
    m_ice_enthalpy.copy_from(E);
    m_ice_enthalpy.update_ghosts();
  }
};

// ---------------------------------------------------------------------------
// Output helpers (same pattern as the other instrument drivers)
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

void write_meta(const std::shared_ptr<const Grid> &grid, int argc, char *argv[],
                const std::string &dir, const std::string &component, const std::string &desc,
                double dt) {
  std::ofstream meta = open_output(dir + "/meta.txt");
  header(meta, "pism_revision", pism::revision);
  header(meta, "pism_config_file", pism::config_file);
  header(meta, "component", component);
  header(meta, "description", desc);
  header(meta, "dt_seconds", std::to_string(dt));
  header(meta, "Mx", std::to_string(grid->Mx()));
  header(meta, "My", std::to_string(grid->My()));
  header(meta, "Mz", std::to_string(grid->Mz()));
  header(meta, "Lx", std::to_string(grid->Lx()));
  header(meta, "Ly", std::to_string(grid->Ly()));
  header(meta, "Lz", std::to_string(grid->Lz()));
  header(meta, "dx", std::to_string(grid->dx()));
  meta << "# command:";
  for (int i = 0; i < argc; ++i) {
    meta << " " << argv[i];
  }
  meta << "\n";
  meta.close();
}

void dump_parameters(const Config &config, const std::string &dir,
                     const std::set<std::string> &prefixes, const std::string &description) {
  std::ofstream out = open_output(dir + "/parameters.csv");
  header(out, "description", description);
  header(out, "columns", "key,value");

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

//! Write a set of arrays to a NetCDF file (same pattern as the other drivers).
void write_netcdf(const std::shared_ptr<const Context> &ctx, const std::string &path,
                  const std::vector<const array::Array *> &vecs) {
  auto writer = std::make_shared<SynchronousOutputWriter>(ctx->com(), *ctx->config());
  writer->initialize({}, true);
  OutputFile file(writer, path);

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

} // namespace
} // namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "energy_age_test");
    auto config = ctx->config();

    set_config_from_options(*config);
    config->resolve_filenames();

    // output directories (energy/ and age/ under the -dumps_dir parent)
    std::string dumps_dir = "stage1/dumps";
    options::String output_dir("-dumps_dir", "Parent output directory", dumps_dir);
    dumps_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + dumps_dir + "/energy " + dumps_dir + "/age";
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directories");
    }

    const double dt_years = options::Real(ctx->unit_system(), "-dt_years",
                                          "Time step length, in years", "years", 10.0);
    auto sys = ctx->unit_system();
    const double dt = units::convert(sys, dt_years, "years", "seconds");

    // grid: test-F domain (same as instrument_sia.cc)
    grid::Parameters P(*config);
    P.Lx = 900.0e3;
    P.Ly = P.Lx;
    const double Lz = 4000.0;
    const unsigned int Mz = config->get_number("grid.Mz");
    P.z = grid::compute_vertical_levels(Lz, Mz, grid::EQUAL);
    P.ownership_ranges_from_options(*config, ctx->size());

    auto grid = std::make_shared<Grid>(ctx, P);

    write_meta(grid, argc, argv, dumps_dir + "/energy", "ice_energy (EnthalpyModel)",
               "test-F exact solution as analytic forcing; one enthalpy step from "
               "the exact steady-state temperature", dt);
    write_meta(grid, argc, argv, dumps_dir + "/age", "ice_age (AgeModel)",
               "test-F exact velocity field; one age step from age = 0", dt);

    dump_parameters(*config, dumps_dir + "/energy",
                    {"energy.", "constants.", "grid.", "enthalpy_converter."},
                    "Configuration parameters affecting the ice-energy component");
    dump_parameters(*config, dumps_dir + "/age",
                    {"age.", "constants.", "grid."},
                    "Configuration parameters affecting the ice-age component");

    // Use the same enthalpy converter the energy model uses internally
    // (the context's; standard EnthalpyConverter with pressure-dependent
    // melting temperature), so the dumped initial state is in the model's
    // own enthalpy convention.
    auto EC = ctx->enthalpy_converter();

    const int WIDE_STENCIL = config->get_number("grid.max_stencil_width");
    const double ice_rho = config->get_number("constants.ice.density"),
                 ice_c = config->get_number("constants.ice.specific_heat_capacity"),
                 Ggeo = 0.042, // W/m^2, the exact test-F geothermal flux
                 Tmin = 223.15, ST = 1.67e-5; // exact test-F surface temperature T_s(r)=Tmin+ST*r

    // ---------------------------------------------------------------
    // Test-F state: geometry + enthalpy from the exact temperature field
    // ---------------------------------------------------------------
    Geometry geometry(grid);
    geometry.sea_level_elevation.set(0.0);
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
    }
    geometry.ensure_consistency(config->get_number("geometry.ice_free_thickness_standard"));

    // ---------------------------------------------------------------
    // Exact 3-D forcing: u3, v3, w3, strain heating Sigma
    // ---------------------------------------------------------------
    array::Array3D u3(grid, "u3", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL),
        v3(grid, "v3", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL),
        w3(grid, "w3", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL),
        sigma(grid, "strain_heating", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL),
        enthalpy(grid, "enthalpy", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL);

    {
      const double LforFG = 750000.0;
      array::AccessScope list{&geometry.ice_thickness, &u3, &v3, &w3, &sigma, &enthalpy};

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double x = grid->x(i), y = grid->y(j);
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        const double H = geometry.ice_thickness(i, j);
        const double inv_r = (r > 1.0) ? (1.0 / r) : 0.0;

        const bool in_sheet = (r <= LforFG - 1.0);
        TestFGParameters F(grid->Mz());
        if (in_sheet) {
          F = exactFG(0.0, r, grid->z(), 0.0);
        }

        double *u_ij = u3.get_column(i, j), *v_ij = v3.get_column(i, j),
               *w_ij = w3.get_column(i, j), *sig_ij = sigma.get_column(i, j),
               *E_ij = enthalpy.get_column(i, j);

        for (unsigned int k = 0; k < grid->Mz(); ++k) {
          if (in_sheet) {
            u_ij[k] = inv_r * x * F.U[k];
            v_ij[k] = inv_r * y * F.U[k];
            w_ij[k] = F.w[k];
            sig_ij[k] = F.Sig[k] * ice_rho * ice_c; // K/s -> W/m^3
          } else {
            u_ij[k] = v_ij[k] = w_ij[k] = sig_ij[k] = 0.0;
          }
          // initial enthalpy from the exact temperature field; ice-free
          // columns are set the way the model sets them (local surface
          // temperature enthalpy) so that E_in == E_out there.
          if (in_sheet) {
            const double depth = H - grid->z(k);
            E_ij[k] = EC->enthalpy_permissive(F.T[k], 0.0, EC->pressure(depth));
          } else {
            E_ij[k] = EC->enthalpy_permissive(Tmin + ST * r, 0.0, 0.0);
          }
        }
      }
      u3.update_ghosts();
      v3.update_ghosts();
      w3.update_ghosts();
      sigma.update_ghosts();
      enthalpy.update_ghosts();
    }

    // ---------------------------------------------------------------
    // 2-D energy forcing fields
    // ---------------------------------------------------------------
    array::Scalar surface_temp(grid, "surface_temp"), shelf_base_temp(grid, "shelf_base_temp"),
        surface_liquid_fraction(grid, "surface_liquid_fraction"),
        till_water_thickness(grid, "till_water_thickness"), basal_heat_flux(grid, "bheatflx"),
        F_b(grid, "basal_frictional_heating");
    shelf_base_temp.set(271.15);
    surface_liquid_fraction.set(0.0);
    till_water_thickness.set(0.0);
    basal_heat_flux.set(Ggeo);
    F_b.set(0.0);
    {
      array::AccessScope list{&surface_temp};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        surface_temp(i, j) = Tmin + ST * r;
      }
    }

    // ---------------------------------------------------------------
    // Enthalpy model: one step from the exact steady-state temperature
    // ---------------------------------------------------------------
    auto energy_model = std::make_shared<InstrumentedEnthalpyModel>(grid);
    energy_model->set_enthalpy(enthalpy);

    array::Array3D E_in(grid, "enthalpy_in", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL);
    E_in.copy_from(energy_model->enthalpy());

    energy::Inputs einputs;
    einputs.cell_type = &geometry.cell_type;
    einputs.basal_frictional_heating = &F_b;
    einputs.basal_heat_flux = &basal_heat_flux;
    einputs.ice_thickness = &geometry.ice_thickness;
    einputs.surface_liquid_fraction = &surface_liquid_fraction;
    einputs.shelf_base_temp = &shelf_base_temp;
    einputs.surface_temp = &surface_temp;
    einputs.till_water_thickness = &till_water_thickness;
    einputs.volumetric_heating_rate = &sigma;
    einputs.u3 = &u3;
    einputs.v3 = &v3;
    einputs.w3 = &w3;

    energy_model->update(0.0, dt, einputs);

    const array::Array3D &E_out = energy_model->enthalpy();
    const array::Scalar &bmelt = energy_model->basal_melt_rate();

    // ---------------------------------------------------------------
    // Age model: one step from age = 0
    // ---------------------------------------------------------------
    auto age_model = std::make_shared<AgeModel>(grid, nullptr);
    age_model->init(InputOptions(INIT_BOOTSTRAP, "", 0));

    array::Array3D age_in(grid, "age_in", array::WITH_GHOSTS, grid->z(), WIDE_STENCIL);
    age_in.copy_from(age_model->age());

    AgeModelInputs ainputs(&geometry.ice_thickness, &u3, &v3, &w3);
    age_model->update(0.0, dt, ainputs);

    const array::Array3D &age_out = age_model->age();

    // ------------------------------------------------------------------
    // Dumps
    // ------------------------------------------------------------------
    {
      // energy inputs + outputs
      const std::vector<const array::Array *> energy_vecs = {
          &geometry.ice_surface_elevation, &geometry.ice_thickness,
          &geometry.bed_elevation,         &geometry.cell_type,
          &surface_temp,                   &shelf_base_temp,
          &surface_liquid_fraction,        &till_water_thickness,
          &basal_heat_flux,                &F_b,
          &sigma,                          &u3,
          &v3,                             &w3,
          &E_in,                           &E_out,
          &bmelt};
      write_netcdf(ctx, dumps_dir + "/energy/inputs.nc", energy_vecs);

      // per-column energy trace
      {
        std::ofstream out = open_output(dumps_dir + "/energy/columns.csv");
        header(out, "description",
               "One enthalpy step from the exact test-F temperature: per (i,j,k): "
               "geometry, surface temperature, exact forcing (u3,v3,w3,sigma), "
               "enthalpy before/after, exact reference temperature");
        header(out, "columns",
               "i,j,x,y,r,H,k,z,depth,surface_temp,u3,v3,w3,sigma,E_in,E_out,T_exact");

        array::AccessScope list{&geometry.ice_thickness, &surface_temp, &u3, &v3, &w3,
                                &sigma, &E_in, &E_out};

        const double LforFG = 750000.0;
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          const double r = std::max(grid::radius(*grid, i, j), 1.0);
          const double H = geometry.ice_thickness(i, j);
          const bool in_sheet = (r <= LforFG - 1.0) and (H > 0.0);
          TestFGParameters F(grid->Mz());
          if (in_sheet) {
            F = exactFG(0.0, r, grid->z(), 0.0);
          }

          const double *u_ij = u3.get_column(i, j), *v_ij = v3.get_column(i, j),
                       *w_ij = w3.get_column(i, j), *sig_ij = sigma.get_column(i, j),
                       *ein_ij = E_in.get_column(i, j), *eout_ij = E_out.get_column(i, j);

          for (unsigned int k = 0; k < grid->Mz(); ++k) {
            out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
                << H << "," << k << "," << grid->z(k) << "," << H - grid->z(k) << ","
                << surface_temp(i, j) << "," << u_ij[k] << "," << v_ij[k] << "," << w_ij[k]
                << "," << sig_ij[k] << "," << ein_ij[k] << "," << eout_ij[k] << ","
                << (in_sheet ? F.T[k] : std::numeric_limits<double>::quiet_NaN()) << "\n";
          }
        }
        out.close();
      }

      // basal melt rate
      {
        std::ofstream out = open_output(dumps_dir + "/energy/basal_melt_rate.csv");
        header(out, "description",
               "Basal melt rate from the energy model (m s^-1; zero for the "
               "cold-based exact test-F state)");
        header(out, "columns", "i,j,x,y,basal_melt_rate");
        array::AccessScope list{&bmelt};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << bmelt(i, j)
              << "\n";
        }
        out.close();
      }
    }

    {
      // age inputs + outputs
      const std::vector<const array::Array *> age_vecs = {
          &geometry.ice_thickness, &u3, &v3, &w3, &age_in, &age_out};
      write_netcdf(ctx, dumps_dir + "/age/inputs.nc", age_vecs);

      // per-column age trace
      {
        std::ofstream out = open_output(dumps_dir + "/age/columns.csv");
        header(out, "description",
               "One age step from age = 0 with the exact test-F velocity field: "
               "per (i,j,k): geometry, u3/v3/w3, age before/after (s)");
        header(out, "columns", "i,j,x,y,r,H,k,z,u3,v3,w3,age_in,age_out");

        array::AccessScope list{&geometry.ice_thickness, &u3, &v3, &w3, &age_in, &age_out};

        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          const double r = std::max(grid::radius(*grid, i, j), 1.0);
          const double H = geometry.ice_thickness(i, j);

          const double *u_ij = u3.get_column(i, j), *v_ij = v3.get_column(i, j),
                       *w_ij = w3.get_column(i, j), *ain_ij = age_in.get_column(i, j),
                       *aout_ij = age_out.get_column(i, j);

          for (unsigned int k = 0; k < grid->Mz(); ++k) {
            out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
                << H << "," << k << "," << grid->z(k) << "," << u_ij[k] << "," << v_ij[k]
                << "," << w_ij[k] << "," << ain_ij[k] << "," << aout_ij[k] << "\n";
          }
        }
        out.close();
      }
    }

    // ------------------------------------------------------------------
    // Validation summaries
    // ------------------------------------------------------------------
    {
      // energy: max |dE|, max |bmelt|, count of temperate cells (should be 0)
      double max_dE = 0.0, max_bmelt = 0.0;
      long temperate = 0;
      const double LforFG = 750000.0;

      array::AccessScope list{&geometry.ice_thickness, &E_in, &E_out, &bmelt};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        const double H = geometry.ice_thickness(i, j);
        if (H <= 0.0) {
          continue;
        }
        max_bmelt = std::max(max_bmelt, std::fabs(bmelt(i, j)));
        const unsigned int ks = grid->kBelowHeight(H);
        const double *ein_ij = E_in.get_column(i, j), *eout_ij = E_out.get_column(i, j);
        for (unsigned int k = 0; k <= ks; ++k) {
          max_dE = std::max(max_dE, std::fabs(eout_ij[k] - ein_ij[k]));
          const double p = EC->pressure(H - grid->z(k));
          if (EC->is_temperate(eout_ij[k], p)) {
            ++temperate;
          }
        }
      }

      const double gmax_dE = GlobalMax(grid->com, max_dE);
      const double gmax_bmelt = GlobalMax(grid->com, max_bmelt);
      const long gtemperate = GlobalSum(grid->com, (int)temperate);

      std::cout << "ENERGY instrumentation complete.\n"
                << "  max |E_out - E_in| = " << gmax_dE << " J/kg"
                << " (dE/dT ~ " << gmax_dE / ice_c << " K)\n"
                << "  max basal melt rate = " << gmax_bmelt << " m/s"
                << " (= " << units::convert(sys, gmax_bmelt, "m second^-1", "m year^-1")
                << " m/yr)\n"
                << "  temperate cells = " << gtemperate
                << " (deep interior near the dome where T_exact > T_pmp)\n"
                << "  reduced_accuracy = " << energy_model->stats().reduced_accuracy_counter
                << ", bulge = " << energy_model->stats().bulge_counter << "\n";
    }

    {
      // age: max age after one step must be ~dt (<= dt, == dt in the interior)
      double max_age = 0.0;
      array::AccessScope list{&geometry.ice_thickness, &age_out};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double H = geometry.ice_thickness(i, j);
        if (H <= 0.0) {
          continue;
        }
        const unsigned int ks = grid->kBelowHeight(H);
        const double *aout_ij = age_out.get_column(i, j);
        for (unsigned int k = 0; k <= ks; ++k) {
          max_age = std::max(max_age, aout_ij[k]);
        }
      }
      const double gmax_age = GlobalMax(grid->com, max_age);
      std::cout << "AGE instrumentation complete.\n"
                << "  dt = " << dt << " s"
                << " (= " << units::convert(sys, dt, "seconds", "years") << " years)\n"
                << "  max age after one step = " << gmax_age << " s"
                << " (= " << units::convert(sys, gmax_age, "seconds", "years") << " years;"
                << " interior must equal dt)\n";
    }

    // CFL sanity for the explicit horizontal advection
    {
      double max_speed = 0.0;
      array::AccessScope list{&u3, &v3};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const unsigned int kt = grid->Mz() - 1;
        const double u = u3.get_column(i, j)[kt], v = v3.get_column(i, j)[kt];
        max_speed = std::max(max_speed, std::sqrt(u * u + v * v));
      }
      const double gmax_speed = GlobalMax(grid->com, max_speed);
      std::cout << "  max surface |(u,v)| = "
                << units::convert(sys, gmax_speed, "m second^-1", "m year^-1")
                << " m/yr; CFL = u*dt/dx = " << gmax_speed * dt / grid->dx() << "\n";
    }

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
