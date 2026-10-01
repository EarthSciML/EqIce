// Stage 1 instrumentation driver for the PISM calving and frontal-melt
// boundaries (boundaries.md section 9 "Front retreat / calving" and section
// 6.4 "Frontal melt").
//
// State: the verification test V "circular flat-bottomed continent" (CFBC)
// configuration -- the van der Veen flow-line floating shelf -- set up exactly
// as instrument_ssa.cc does, but on a *square* grid (Mx = My) with a few
// ocean columns beyond the calving front:
//
//   * the rate-based calving components (EigenCalving, vonMisesCalving,
//     HayhurstCalving) reject non-square grid cells in init()
//     (|dx - dy| / min(dx, dy) <= 1e-2 is required), so the SSA run's
//     61 x 3 flow-line grid cannot be used;
//   * the front cell (first ice-free-ocean cell, next to the floating shelf)
//     must be at least 2 interior cells away from the east domain edge because
//     the calving rate is sampled at +-2 cells from the front.
//
// The shelf state is uniform in y, so the exact solution H(x),
// u(x) = V0 H0 / H(x), v = 0 is unchanged. The actual SSA solve is run to
// produce the velocity the calving models act on (as in the model).
//
// Dumps (instantaneous):
//   inputs.nc                        - geometry, SSA velocity, enthalpy, and all
//                                      calving / frontal-melt rate outputs
//   eigen_calving.csv                - per-cell strain-rate invariants (eigen1,
//                                      eigen2) + rate (exact-velocity run);
//                                      with the shelf velocity eigen2 = 0
//                                      exactly and the rate is zero (tests the
//                                      compressive branch)
//   eigen_calving_divergent.csv      - same with a manufactured divergent
//                                      velocity (v = 0.05 u cos(2 pi y / L));
//                                      exercises the nonzero branch
//   vonmises_calving.csv             - per-cell hardness, eigen rates,
//                                      sigma_tilde, rate
//   hayhurst_calving.csv             - per-cell omega, sigma_0, rate
//   calving_at_thickness.csv         - mask/thickness in/out for the
//                                      thickness-threshold calving model
//   float_kill.csv                   - mask/thickness in/out for float-kill
//   frontal_melt_undercutting.csv    - FrontalMeltPhysics pointwise kernel
//   frontal_melt_ismip6.csv          - FrontalMeltPhysics pointwise kernel
//   frontal_melt_constant.csv        - Constant frontal-melt model melt +
//                                      retreat-rate maps
//   parameters.csv, meta.txt
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/calving.json

#include <petsc.h>
#include <netcdf.h>
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
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/CellType.hh"
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/ssa/SSA.hh"
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFDBase.hh"
#include "pism/frontretreat/calving/EigenCalving.hh"
#include "pism/frontretreat/calving/vonMisesCalving.hh"
#include "pism/frontretreat/calving/HayhurstCalving.hh"
#include "pism/frontretreat/calving/CalvingAtThickness.hh"
#include "pism/frontretreat/calving/FloatKill.hh"
#include "pism/coupler/FrontalMelt.hh"
#include "pism/coupler/frontalmelt/FrontalMeltPhysics.hh"
#include "pism/coupler/frontalmelt/Constant.hh"
#include "pism/rheology/FlowLawFactory.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM calving and frontal-melt\n"
  "components. Dumps calving-rate maps, in-place mask-removal traces, and\n"
  "frontal-melt kernel samples on the van der Veen CFBC shelf state to\n"
  "stage1/dumps/calving/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

using pism::stressbalance::SSAFD;

// ---------------------------------------------------------------------------
// Subclasses exposing protected internals
// ---------------------------------------------------------------------------

class InstrumentedSSAFD : public SSAFD {
public:
  InstrumentedSSAFD(std::shared_ptr<const Grid> g) : SSAFD(g, false) {}
};

class InstrumentedEigenCalving : public calving::EigenCalving {
public:
  InstrumentedEigenCalving(std::shared_ptr<const Grid> grid)
    : calving::EigenCalving(grid) {}
  const array::Array2D<stressbalance::PrincipalStrainRates> &strain_rates() const {
    return m_strain_rates;
  }
};

class InstrumentedVonMisesCalving : public calving::vonMisesCalving {
public:
  InstrumentedVonMisesCalving(std::shared_ptr<const Grid> grid,
                              std::shared_ptr<const rheology::FlowLaw> flow_law)
    : calving::vonMisesCalving(grid, flow_law) {}
  const array::Array2D<stressbalance::PrincipalStrainRates> &strain_rates() const {
    return m_strain_rates;
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

void dump_parameters(const Config &config, const std::string &dir) {
  std::ofstream out = open_output(dir + "/parameters.csv");
  header(out, "description",
         "Configuration parameters affecting the calving and frontal-melt "
         "components (calving.*, frontal_melt.*, constants.*, grid.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {
      "calving.", "frontal_melt.", "constants.", "grid."};

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

// ---------------------------------------------------------------------------
// Minimal NetCDF-3 input file for the thickness-calving model. The constructor
// of CalvingAtThickness requires an input file (it reads the
// "calving.thickness_calving.file" config key); because the file does not
// contain the "thickness_calving_threshold" variable, init() falls back to the
// constant threshold from "calving.thickness_calving.threshold".
// ---------------------------------------------------------------------------

void create_thickness_calving_input_file(const std::string &path) {
  int ncid;
  if (nc_create(path.c_str(), NC_CLOBBER, &ncid) != NC_NOERR) {
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Cannot create the thickness-calving input file '" + path + "'");
  }
  int dim_x = 0, dim_y = 0, dim_t = 0, var_t = 0;
  if (nc_def_dim(ncid, "x", 1, &dim_x) != NC_NOERR or
      nc_def_dim(ncid, "y", 1, &dim_y) != NC_NOERR or
      nc_def_dim(ncid, "time", NC_UNLIMITED, &dim_t) != NC_NOERR or
      nc_def_var(ncid, "time", NC_DOUBLE, 1, &dim_t, &var_t) != NC_NOERR or
      nc_enddef(ncid) != NC_NOERR) {
    nc_close(ncid);
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Cannot define variables in the thickness-calving input file '" + path + "'");
  }
  const double t0 = 0.0;
  if (nc_put_var_double(ncid, var_t, &t0) != NC_NOERR) {
    nc_close(ncid);
    throw RuntimeError(PISM_ERROR_LOCATION,
                       "Cannot write to the thickness-calving input file '" + path + "'");
  }
  nc_close(ncid);
}

// ---------------------------------------------------------------------------
// van der Veen exact solution (from pism_ssa_test_cfbc / instrument_ssa.cc)
// ---------------------------------------------------------------------------

double H_exact(double V0, double H0, double C, double x) {
  const double Q0 = V0 * H0;
  return pow(4.0 * C / Q0 * x + 1.0 / pow(H0, 4), -0.25);
}

double u_exact(double V0, double H0, double C, double x) {
  const double Q0 = V0 * H0;
  return Q0 / H_exact(V0, H0, C, x);
}

// analytic along-flow strain rate du/dx = C H^3 (derived from u = Q0/H)
double u_x_analytic(double C, double H) {
  return C * pow(H, 3);
}

// ---------------------------------------------------------------------------
// Reconstruct the EigenCalving rate from the dumped strain-rate invariants and
// the mask, exactly as EigenCalving::update does (the .esm test target).
// offset = m_stencil_width = 2; eigenCalvOffset = 0.
// ---------------------------------------------------------------------------

double eigen_calving_reconstructed(const array::CellType1 &mask,
                                   const array::Array2D<stressbalance::PrincipalStrainRates> &strain,
                                   double K, int offset, int i, int j, double *e1_out,
                                   double *e2_out, int *N_out,
                                   unsigned int Mx, unsigned int My) {
  if (not mask.ice_free_ocean(i, j) or not mask.next_to_floating_ice(i, j) or
      // the mask has a 1-cell-wide ghost region only; the +-offset sampling
      // must stay inside the interior (this excludes front cells at the
      // domain edges, which the model itself would sample out of bounds)
      i < (int)offset or i > (int)Mx - 1 - (int)offset) {
    *e1_out = std::numeric_limits<double>::quiet_NaN();
    *e2_out = std::numeric_limits<double>::quiet_NaN();
    *N_out = -1;
    return 0.0;
  }
  double eigen1 = 0.0, eigen2 = 0.0;
  int N = 0;
  for (int p = -1; p < 2; p += 2) {
    const int I = i + p * offset;
    if (mask.floating_ice(I, j) and not mask.ice_margin(I, j)) {
      eigen1 += strain(I, j).eigen1;
      eigen2 += strain(I, j).eigen2;
      N += 1;
    }
  }
  // the y-direction is periodic but the +-offset sampling is only valid away
  // from the y edges for the 1-cell-wide mask ghosts
  if (j >= (int)offset and j <= (int)My - 1 - (int)offset) {
    for (int q = -1; q < 2; q += 2) {
      const int J = j + q * offset;
      if (mask.floating_ice(i, J) and not mask.ice_margin(i, J)) {
        eigen1 += strain(i, J).eigen1;
        eigen2 += strain(i, J).eigen2;
        N += 1;
      }
    }
  }
  if (N > 0) {
    eigen1 /= N;
    eigen2 /= N;
  }
  const double eigenCalvOffset = 0.0;
  double rate = 0.0;
  if (eigen2 > eigenCalvOffset and eigen1 > 0.0) {
    rate = K * eigen1 * (eigen2 - eigenCalvOffset);
  }
  *e1_out = eigen1;
  *e2_out = eigen2;
  *N_out = N;
  return rate;
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
    std::shared_ptr<Context> ctx = context_from_options(com, "calving_test");
    auto config = ctx->config();

    // Same configuration as pism_ssa_test_cfbc (test V / van der Veen):
    config->set_number("flow_law.isothermal_Glen.ice_softness",
                       pow(1.9e8, -config->get_number("stress_balance.ssa.Glen_exponent")));
    config->set_flag("stress_balance.ssa.compute_surface_gradient_inward", false);
    config->set_flag("stress_balance.calving_front_stress_bc", true);
    config->set_flag("stress_balance.ssa.fd.flow_line_mode", true);
    config->set_flag("stress_balance.ssa.fd.extrapolate_at_margins", false);
    config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");

    // Calving / frontal-melt configuration:
    //   - a nonzero eigen-calving constant K so the formula branch is
    //     exercised (the config default is 0.0, which would zero the rate).
    //   - a thickness-calving threshold of 300 m so the front (H ~ 224 m)
    //     is calved; the input file is created below and contains no
    //     "thickness_calving_threshold" variable, so init() uses this constant.
    //   - frontal melt: keep the constant model's default melt rate (1 m/day).
    config->set_number("calving.eigen_calving.K", 3.0e16); // m s, ~1e9 m a
    config->set_number("calving.thickness_calving.threshold", 300.0); // m

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/calving";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    const std::string threshold_file = out_dir + "/thickness_threshold_input.nc";
    create_thickness_calving_input_file(threshold_file);
    config->set_string("calving.thickness_calving.file", threshold_file);

    const unsigned int Mx = config->get_number("grid.Mx");
    const unsigned int My = config->get_number("grid.My");

    // test V domain: 250 km half-width, cell-centered, periodic in Y, square
    // cells (required by the rate-based calving components)
    auto grid = Grid::Shallow(ctx, 250.0e3, 250.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::Y_PERIODIC);
    const double dx = grid->dx(), dy = grid->dy();

    // the shelf occupies columns i in [0, N_ice-1] and the calving front is at
    // the first ice-free-ocean column i = N_ice (H = 0 there and beyond); the
    // front cell must stay >= 2 cells from the east edge for the +-2 sampling.
    const int N_ice = (int)Mx - 5;

    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "test", "V (van der Veen flow-line shelf, cfbc)");
      header(meta, "flow_law", "isothermal_glen, B = 1.9e8 Pa s^(1/n)");
      header(meta, "sliding", "none (floating shelf, tauc = 0)");
      header(meta, "Mx", std::to_string(Mx));
      header(meta, "My", std::to_string(My));
      header(meta, "Lx", "250000");
      header(meta, "dx", std::to_string(dx));
      header(meta, "dy", std::to_string(dy));
      header(meta, "grid_note",
             "square grid (Mx = My) required by EigenCalving / vonMisesCalving / "
             "HayhurstCalving init(); 5 ocean columns beyond the calving front so "
             "the +-2-cell sampling stays in the interior");
      header(meta, "calving_front_column", std::to_string(N_ice));
      header(meta, "eigen_calving_K_m_s", "3e16 (config default is 0.0)");
      header(meta, "thickness_calving_threshold_m", "300 (config default is 50)");
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    EnthalpyConverter EC(*config);
    auto sys = ctx->unit_system();

    const double V0 = units::convert(sys, 300.0, "m year^-1", "m second^-1");
    const double H0 = 600.0;
    const double C = 2.45e-18;
    const double to_m_yr = units::convert(sys, 1.0, "m second^-1", "m year^-1");
    const double to_m_day = units::convert(sys, 1.0, "m second^-1", "m day^-1");

    array::Scalar1 tauc(grid, "tauc");
    array::Array3D enthalpy(grid, "enthalpy", array::WITH_GHOSTS, grid->z(), 1);
    array::Vector2 bc_values(grid, "_bc");
    array::Scalar2 bc_mask(grid, "bc_mask");
    bc_mask.set_interpolation_type(NEAREST);

    Geometry geometry(grid);

    // initialize coefficients (as in SSATestCaseCFBC / instrument_ssa.cc)
    tauc.set(0.0);
    geometry.bed_elevation.set(-1000.0);
    enthalpy.set(EC.enthalpy(273.15, 0.01, 0.0));

    {
      array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                              &bc_mask, &bc_values, &geometry.cell_type};

      const double x_min = grid->x(0);

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double x = grid->x(i);

        if (i < N_ice) {
          geometry.ice_thickness(i, j) = H_exact(V0, H0, C, x - x_min);
        } else {
          geometry.ice_thickness(i, j) = 0.0;
        }

        if (i == 0) {
          bc_mask(i, j) = 1;
          bc_values(i, j) = {V0, 0.0};
        } else {
          bc_mask(i, j) = 0;
          bc_values(i, j) = {0.0, 0.0};
        }
      }
    }
    geometry.ensure_consistency(0.0);
    bc_mask.update_ghosts();
    bc_values.update_ghosts();

    // SSA solve (FD Picard) -- the velocity the calving models act on
    auto ssa = std::make_shared<InstrumentedSSAFD>(grid);
    ssa->init();

    Inputs ssa_inputs;
    ssa_inputs.water_column_pressure = nullptr;
    ssa_inputs.geometry = &geometry;
    ssa_inputs.enthalpy = &enthalpy;
    ssa_inputs.basal_yield_stress = &tauc;
    ssa_inputs.bc_mask = &bc_mask;
    ssa_inputs.bc_values = &bc_values;

    ssa->update(ssa_inputs, true);

    const array::Vector1 &velocity_ssa = ssa->velocity();
    // local copy with valid ghosts for the strain-rate stencils
    array::Vector1 velocity(grid, "velocity");
    velocity.copy_from(velocity_ssa);
    velocity.update_ghosts();

    // manufactured divergent velocity to exercise the nonzero eigen-calving
    // branch: v = 0.05 u cos(2 pi y / L), periodic in y
    const double L = 2.0 * 250.0e3;
    array::Vector1 velocity_div(grid, "velocity_divergent");
    {
      velocity_div.copy_from(velocity);
      array::AccessScope list{&velocity_div};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double v_pert = 0.05 * velocity_div(i, j).u * std::cos(2.0 * M_PI * grid->y(j) / L);
        velocity_div(i, j) = {velocity_div(i, j).u, v_pert};
      }
      velocity_div.update_ghosts();
    }

    // ------------------------------------------------------------------
    // (a) Rate-based calving components
    // ------------------------------------------------------------------

    const double K_eigen = config->get_number("calving.eigen_calving.K");
    const double sigma_max = config->get_number("calving.vonmises_calving.sigma_max");

    auto eigen_calving = std::make_shared<InstrumentedEigenCalving>(grid);
    eigen_calving->init();
    eigen_calving->update(geometry.cell_type, velocity);

    const array::Array2D<stressbalance::PrincipalStrainRates> &eigen_strain =
        eigen_calving->strain_rates();

    auto eigen_calving_div = std::make_shared<InstrumentedEigenCalving>(grid);
    eigen_calving_div->init();
    eigen_calving_div->update(geometry.cell_type, velocity_div);

    // von Mises calving with the same flow law the SSA uses (isothermal_glen,
    // B = 1.9e8, n = 3)
    rheology::FlowLawFactory flow_law_factory(config, grid->ctx()->enthalpy_converter());
    auto flow_law = flow_law_factory.create(
        config->get_string("stress_balance.ssa.flow_law"),
        config->get_number("stress_balance.ssa.Glen_exponent"));

    auto vonmises_calving = std::make_shared<InstrumentedVonMisesCalving>(grid, flow_law);
    vonmises_calving->init();
    vonmises_calving->update(geometry.cell_type, geometry.ice_thickness, velocity, enthalpy);

    const array::Array2D<stressbalance::PrincipalStrainRates> &vonmises_strain =
        vonmises_calving->strain_rates();

    auto hayhurst_calving = std::make_shared<calving::HayhurstCalving>(grid);
    hayhurst_calving->init();
    hayhurst_calving->update(geometry.cell_type, geometry.ice_thickness,
                             geometry.sea_level_elevation, geometry.bed_elevation);

    const array::Scalar &eigen_rate = eigen_calving->calving_rate();
    const array::Scalar &vonmises_rate = vonmises_calving->calving_rate();
    const array::Scalar &hayhurst_rate = hayhurst_calving->calving_rate();

    // per-cell strain-rate / calving-rate dumps
    {
      auto dump_strain = [&](const std::string &path, const char *description,
                             const array::Array2D<stressbalance::PrincipalStrainRates> &strain,
                             const array::Vector1 &vel, const array::Scalar &rate,
                             bool divergent) {
        std::ofstream out = open_output(path);
        header(out, "description", description);
        header(out, "columns",
               "i,j,x,y,H,mask,u,v,eigen1,eigen2,du_dx_analytic,"
               "eigen1_avg,eigen2_avg,N_avg,rate_m_s,rate_m_yr");

        array::AccessScope list{&geometry.ice_thickness, &geometry.cell_type, &vel,
                                &strain, &rate};
        const double x_min = grid->x(0);
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          const double x = grid->x(i);
          const bool have_exact = (i < N_ice);
          const double H = geometry.ice_thickness(i, j);
          const double Hx = have_exact ? H_exact(V0, H0, C, x - x_min) : 0.0;

          double e1_avg = std::numeric_limits<double>::quiet_NaN();
          double e2_avg = std::numeric_limits<double>::quiet_NaN();
          int N_avg = -1;
          const double rate_recon =
              eigen_calving_reconstructed(geometry.cell_type, strain, K_eigen, 2,
                                          i, j, &e1_avg, &e2_avg, &N_avg, Mx, My);

          out << i << "," << j << "," << x << "," << grid->y(j) << "," << H << ","
              << (int)geometry.cell_type(i, j) << "," << vel(i, j).u << "," << vel(i, j).v
              << "," << strain(i, j).eigen1 << "," << strain(i, j).eigen2 << ","
              << (divergent ? std::numeric_limits<double>::quiet_NaN()
                            : u_x_analytic(C, Hx))
              << "," << e1_avg << "," << e2_avg << "," << N_avg << ","
              << rate(i, j) << "," << rate(i, j) * to_m_yr << "\n";

          if (divergent and rate(i, j) != rate_recon) {
            std::cout << "  [WARN] divergent eigen-calving reconstruction mismatch at ("
                      << i << "," << j << "): component " << rate(i, j)
                      << " vs reconstructed " << rate_recon << "\n";
          }
        }
        out.close();
      };

      dump_strain(out_dir + "/eigen_calving.csv",
                  "EigenCalving::update on the van der Veen shelf with the exact "
                  "SSA velocity: per-cell strain-rate invariants (eigen1 >= eigen2, "
                  "from compute_2D_principal_strain_rates) and the calving rate; "
                  "eigen1_avg/eigen2_avg/N_avg are the front-cell sampling average "
                  "reconstructed as in EigenCalving::update (offset 2); with the "
                  "shelf velocity eigen2 = 0 exactly so the rate is zero "
                  "(compressive branch). du_dx_analytic = C H^3.",
                  eigen_strain, velocity, eigen_rate, false);

      dump_strain(out_dir + "/eigen_calving_divergent.csv",
                  "EigenCalving::update with a manufactured divergent velocity "
                  "v = 0.05 u cos(2 pi y / L) (u = SSA velocity): exercises the "
                  "nonzero branch rate = K * eigen1_avg * (eigen2_avg - 0).",
                  eigen_calving_div->strain_rates(), velocity_div,
                  eigen_calving_div->calving_rate(), true);
    }

    // von Mises calving dump: per-cell strain rates and rate; at the front
    // cells also the reconstructed averaging (velocity magnitude, hardness,
    // effective tensile strain rate, sigma_tilde) used by the component.
    {
      std::ofstream out = open_output(out_dir + "/vonmises_calving.csv");
      header(out, "description",
             "vonMisesCalving::update on the van der Veen shelf (isothermal_glen "
             "B = 1.9e8, n = 3, cold ice): per-cell eigen rates and calving rate. "
             "At front cells (ice-free ocean next to ice) the reconstructed "
             "averages over icy +-offset cells are dumped (|v|_avg, hardness_avg, "
             "e_s, sigma_tilde); rate = |v|_avg sigma_tilde / sigma_max.");
      header(out, "columns",
             "i,j,x,y,H,mask,u,v,eigen1,eigen2,|v|_avg,hardness_avg,e_s,sigma_tilde,"
             "sigma_max,rate_m_s,rate_m_yr");

      const double glen_exponent = flow_law->exponent();
      const double *z = enthalpy.levels().data();

      array::AccessScope list{&geometry.ice_thickness, &geometry.cell_type, &velocity,
                              &vonmises_strain, &vonmises_rate, &enthalpy};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();

        double vel_avg = std::numeric_limits<double>::quiet_NaN();
        double hard_avg = std::numeric_limits<double>::quiet_NaN();
        double e_s = std::numeric_limits<double>::quiet_NaN();
        double sigma_tilde = std::numeric_limits<double>::quiet_NaN();

        if (geometry.cell_type.ice_free_ocean(i, j) and geometry.cell_type.next_to_ice(i, j) and
            i >= 2 and i <= (int)Mx - 3) {
          double vel_sum = 0.0, hard_sum = 0.0, e1_sum = 0.0, e2_sum = 0.0;
          int N = 0;
          for (int p = -1; p < 2; p += 2) {
            const int I = i + p * 2;
            if (geometry.cell_type.icy(I, j)) {
              vel_sum += velocity(I, j).magnitude();
              const double H = geometry.ice_thickness(I, j);
              auto k = grid->kBelowHeight(H);
              hard_sum += averaged_hardness(*flow_law, H, k, z, enthalpy.get_column(I, j));
              e1_sum += vonmises_strain(I, j).eigen1;
              e2_sum += vonmises_strain(I, j).eigen2;
              N += 1;
            }
          }
          if (j >= 2 and j <= (int)My - 3) {
            for (int q = -1; q < 2; q += 2) {
              const int J = j + q * 2;
              if (geometry.cell_type.icy(i, J)) {
                vel_sum += velocity(i, J).magnitude();
                const double H = geometry.ice_thickness(i, J);
                auto k = grid->kBelowHeight(H);
                hard_sum += averaged_hardness(*flow_law, H, k, z, enthalpy.get_column(i, J));
                e1_sum += vonmises_strain(i, J).eigen1;
                e2_sum += vonmises_strain(i, J).eigen2;
                N += 1;
              }
            }
          }
          if (N > 0) {
            vel_avg = vel_sum / N;
            hard_avg = hard_sum / N;
            const double e1 = e1_sum / N, e2 = e2_sum / N;
            e_s = std::sqrt(0.5 * (std::pow(std::max(0.0, e1), 2) +
                                   std::pow(std::max(0.0, e2), 2)));
            sigma_tilde = std::sqrt(3.0) * hard_avg * std::pow(e_s, 1.0 / glen_exponent);
          }
        }

        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << ","
            << geometry.ice_thickness(i, j) << "," << (int)geometry.cell_type(i, j) << ","
            << velocity(i, j).u << "," << velocity(i, j).v << ","
            << vonmises_strain(i, j).eigen1 << "," << vonmises_strain(i, j).eigen2 << ","
            << vel_avg << "," << hard_avg << "," << e_s << "," << sigma_tilde << ","
            << sigma_max << "," << vonmises_rate(i, j) << "," << vonmises_rate(i, j) * to_m_yr
            << "\n";
      }
      out.close();
    }

    // Hayhurst calving dump: pure pointwise rate.
    {
      std::ofstream out = open_output(out_dir + "/hayhurst_calving.csv");
      header(out, "description",
             "HayhurstCalving::update on the van der Veen shelf: per-cell inputs "
             "(H, sea_level, bed, water_depth, omega, sigma_0, threshold) and the "
             "calving rate. omega = water_depth / H (with the floating-shelf "
             "adjustment), sigma_0 = (0.4 - 0.45 (omega - 0.065)^2) rho_i g H.");
      header(out, "columns",
             "i,j,x,y,H,mask,sea_level,bed,water_depth,omega,sigma_0_Pa,"
             "sigma_threshold_Pa,rate_m_s,rate_m_yr");

      const double ice_rho = config->get_number("constants.ice.density");
      const double water_rho = config->get_number("constants.sea_water.density");
      const double gravity = config->get_number("constants.standard_gravity");
      const double B_tilde = config->get_number("calving.hayhurst_calving.B_tilde");
      const double exponent_r = config->get_number("calving.hayhurst_calving.exponent_r");
      const double sigma_threshold =
          config->get_number("calving.hayhurst_calving.sigma_threshold", "Pa");

      array::AccessScope list{&geometry.ice_thickness, &geometry.cell_type,
                              &geometry.sea_level_elevation, &geometry.bed_elevation,
                              &hayhurst_rate};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double H = geometry.ice_thickness(i, j);
        const double sea_level = geometry.sea_level_elevation(i, j);
        const double bed = geometry.bed_elevation(i, j);
        const double water_depth = sea_level - bed;

        double omega = std::numeric_limits<double>::quiet_NaN();
        double sigma_0 = std::numeric_limits<double>::quiet_NaN();
        if (geometry.cell_type.icy(i, j) and water_depth > 0.0) {
          omega = water_depth / H;
          if (omega > ice_rho / water_rho) {
            const double freeboard = (1.0 - ice_rho / water_rho) * H;
            const double H_fixed = water_depth + freeboard;
            omega = water_depth / H_fixed;
          }
          sigma_0 = std::max((0.4 - 0.45 * std::pow(omega - 0.065, 2.0)) * ice_rho * gravity * H,
                             sigma_threshold);
        }

        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << H << ","
            << (int)geometry.cell_type(i, j) << "," << sea_level << "," << bed << ","
            << water_depth << "," << omega << "," << sigma_0 << "," << sigma_threshold << ","
            << hayhurst_rate(i, j) << "," << hayhurst_rate(i, j) * to_m_yr << "\n";
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // (b) In-place mask-removal models (on copies of the state)
    // ------------------------------------------------------------------

    // CalvingAtThickness
    {
      array::CellType1 mask(grid, "cat_mask");
      array::Scalar1 thickness(grid, "cat_thickness");
      mask.copy_from(geometry.cell_type);
      thickness.copy_from(geometry.ice_thickness);

      auto cat = std::make_shared<calving::CalvingAtThickness>(grid);
      cat->init();
      const double dt = units::convert(sys, 1.0, "year", "second");
      cat->update(0.0, dt, mask, thickness);
      const array::Scalar &threshold = cat->threshold();

      std::ofstream out = open_output(out_dir + "/calving_at_thickness.csv");
      header(out, "description",
             "CalvingAtThickness::update(t=0, dt=1 yr, cell_type, thickness) on the "
             "van der Veen shelf with calving.thickness_calving.threshold = 300 m: "
             "cells that are floating AND next to ice-free ocean AND H < threshold "
             "are removed in place (thickness = 0, mask = ice-free ocean).");
      header(out, "columns",
             "i,j,x,y,H_in,H_out,mask_in,mask_out,threshold,removed");

      array::AccessScope list{&geometry.cell_type, &geometry.ice_thickness, &mask,
                              &thickness, &threshold};
      int n_removed = 0;
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const int mask_in = (int)geometry.cell_type(i, j);
        const int mask_out = (int)mask(i, j);
        const bool removed = (mask_out != mask_in);
        if (removed) {
          n_removed += 1;
        }
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << ","
            << geometry.ice_thickness(i, j) << "," << thickness(i, j) << ","
            << mask_in << "," << mask_out << "," << threshold(i, j) << ","
            << (removed ? 1 : 0) << "\n";
      }
      out.close();
      std::cout << "CALVING AT THICKNESS: removed " << n_removed << " cells\n";
    }

    // FloatKill
    {
      array::CellType1 mask(grid, "fk_mask");
      array::Scalar1 thickness(grid, "fk_thickness");
      mask.copy_from(geometry.cell_type);
      thickness.copy_from(geometry.ice_thickness);

      auto fk = std::make_shared<calving::FloatKill>(grid);
      fk->init();
      fk->update(mask, thickness);

      std::ofstream out = open_output(out_dir + "/float_kill.csv");
      header(out, "description",
             "FloatKill::update(cell_type, thickness) on the van der Veen shelf "
             "(defaults: margin_only = no, calve_near_grounding_line = yes): all "
             "floating cells are removed in place. h_f = (sea_level - bed) rho_w / "
             "rho_i is the analytic flotation thickness; a cell is floating iff "
             "H < h_f (and bed below sea level).");
      header(out, "columns",
             "i,j,x,y,H_in,H_out,mask_in,mask_out,h_f,floating,removed");

      const double ice_rho = config->get_number("constants.ice.density");
      const double water_rho = config->get_number("constants.sea_water.density");

      array::AccessScope list{&geometry.cell_type, &geometry.ice_thickness, &mask,
                              &thickness, &geometry.bed_elevation, &geometry.sea_level_elevation};
      int n_removed = 0, n_floating = 0;
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const int mask_in = (int)geometry.cell_type(i, j);
        const int mask_out = (int)mask(i, j);
        const bool removed = (mask_out != mask_in);
        const double h_f =
            (geometry.sea_level_elevation(i, j) - geometry.bed_elevation(i, j)) * water_rho / ice_rho;
        const bool floating =
            mask::floating_ice(mask_in) and geometry.ice_thickness(i, j) > 0.0;
        if (floating) {
          n_floating += 1;
        }
        if (removed) {
          n_removed += 1;
        }
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << ","
            << geometry.ice_thickness(i, j) << "," << thickness(i, j) << ","
            << mask_in << "," << mask_out << "," << h_f << "," << (floating ? 1 : 0) << ","
            << (removed ? 1 : 0) << "\n";
      }
      out.close();
      std::cout << "FLOAT KILL: removed " << n_removed << " of " << n_floating
                << " floating cells\n";
    }

    // ------------------------------------------------------------------
    // (c) Frontal-melt physics kernels + the Constant model
    // ------------------------------------------------------------------

    frontalmelt::FrontalMeltPhysics physics(*config);

    const std::vector<double> h_grid = {0.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0};
    const std::vector<double> q_grid = {-0.1, 0.0, 0.1, 0.5, 1.0, 2.0, 5.0};
    const std::vector<double> T_grid = {-1.0, 0.0, 0.5, 1.0, 2.0, 4.0};

    {
      std::ofstream out = open_output(out_dir + "/frontal_melt_undercutting.csv");
      header(out, "description",
             "FrontalMeltPhysics::frontal_melt_from_undercutting(h, q_sg, TF): "
             "q_m = (A h q_sg^alpha + B) TF^beta with the routing tuning parameters; "
             "returns 0 for h <= 0 or q_sg < 0 or TF < 0. Units: h m, q_sg m/day, "
             "TF Celsius, q_m m/day.");
      header(out, "columns", "h,q_sg,TF,q_m_m_day");

      for (double h : h_grid) {
        for (double q : q_grid) {
          for (double TF : T_grid) {
            out << h << "," << q << "," << TF << ","
                << physics.frontal_melt_from_undercutting(h, q, TF) << "\n";
          }
        }
      }
      out.close();
    }

    {
      std::ofstream out = open_output(out_dir + "/frontal_melt_ismip6.csv");
      header(out, "description",
             "FrontalMeltPhysics::frontal_melt_from_ismip6(h, q_sg, TF): the same "
             "formula but WITHOUT the nonnegativity clamps (applies the raw power "
             "law for negative inputs too).");
      header(out, "columns", "h,q_sg,TF,q_m_m_day");

      for (double h : h_grid) {
        for (double q : q_grid) {
          for (double TF : T_grid) {
            out << h << "," << q << "," << TF << ","
                << physics.frontal_melt_from_ismip6(h, q, TF) << "\n";
          }
        }
      }
      out.close();
    }

    // Constant frontal-melt model: default config (grounded only) and with
    // frontal_melt.include_floating_ice = true.
    {
      array::Scalar1 water_flux(grid, "subglacial_water_flux");
      water_flux.set(0.0);

      std::ofstream out = open_output(out_dir + "/frontal_melt_constant.csv");
      header(out, "description",
             "frontalmelt::Constant on the van der Veen shelf (fully floating): "
             "default config (frontal_melt.include_floating_ice = no) gives a zero "
             "melt rate because apply() only targets grounded ice; with "
             "include_floating_ice = yes the melt rate equals the constant "
             "frontal_melt.constant.melt_rate on icy cells and on ice-free-ocean "
             "cells next to ice. The retreat rate is frontal_melt_rate scaled by "
             "H_submerged / H_threshold (part-grid threshold thickness); for a "
             "floating front H_submerged = (rho_i/rho_w) H_threshold.");
      header(out, "columns",
             "i,j,x,y,H,mask,include_floating_ice,melt_rate_m_day,retreat_rate_m_day,"
             "H_threshold,H_submerged");

      const double melt_rate_config =
          config->get_number("frontal_melt.constant.melt_rate", "m day-1");

      for (bool include_floating : {false, true}) {
        config->set_flag("frontal_melt.include_floating_ice", include_floating);
        auto fm = std::make_shared<frontalmelt::Constant>(grid);
        fm->init(geometry);

        FrontalMeltInputs inputs;
        inputs.geometry = &geometry;
        inputs.subglacial_water_flux = &water_flux;
        fm->update(inputs, 0.0, units::convert(sys, 1.0, "year", "second"));

        const array::Scalar &melt_rate = fm->frontal_melt_rate();
        const array::Scalar &retreat_rate = fm->retreat_rate();

        GeometryCalculator gc(*config);
        const double ice_rho = config->get_number("constants.ice.density");
        const double water_rho = config->get_number("constants.sea_water.density");
        const double alpha = ice_rho / water_rho;

        array::AccessScope list{&geometry.cell_type, &geometry.ice_thickness,
                                &geometry.ice_surface_elevation, &geometry.bed_elevation,
                                &geometry.sea_level_elevation, &melt_rate, &retreat_rate};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();

          double H_threshold = std::numeric_limits<double>::quiet_NaN();
          double H_submerged = std::numeric_limits<double>::quiet_NaN();
          if (geometry.cell_type.ice_free_ocean(i, j) and geometry.cell_type.next_to_ice(i, j)) {
            const double bed = geometry.bed_elevation(i, j);
            const double sea_level = geometry.sea_level_elevation(i, j);
            auto H = geometry.ice_thickness.star(i, j);
            auto h = geometry.ice_surface_elevation.star(i, j);
            auto M = geometry.cell_type.star_int(i, j);
            H_threshold = part_grid_threshold_thickness(M, H, h, bed);
            const int m = gc.mask(sea_level, bed, H_threshold);
            H_submerged = (mask::grounded(m) ? std::max(sea_level - bed, 0.0)
                                             : alpha * H_threshold);
          }

          out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << ","
              << geometry.ice_thickness(i, j) << "," << (int)geometry.cell_type(i, j) << ","
              << (include_floating ? 1 : 0) << ","
              << melt_rate(i, j) * to_m_day << "," << retreat_rate(i, j) * to_m_day << ","
              << H_threshold << "," << H_submerged << "\n";
        }
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // NetCDF: geometry + velocity + enthalpy + all rate outputs
    // ------------------------------------------------------------------
    {
      const std::vector<const array::Array *> vecs = {
          &geometry.ice_surface_elevation, &geometry.ice_thickness,
          &geometry.bed_elevation,         &geometry.cell_type,
          &geometry.sea_level_elevation,   &velocity,
          &enthalpy,                       &eigen_rate,
          &vonmises_rate,                  &hayhurst_rate};
      write_netcdf(ctx, out_dir + "/inputs.nc", vecs);
    }

    // ------------------------------------------------------------------
    // Validation summaries (printed to stdout; recorded in runs/calving.json)
    // ------------------------------------------------------------------

    // (a) strain-rate verification: eigen1 must match the analytic du/dx = C H^3
    //     on interior floating cells (v = 0 shelf; eigen1 = u_x, eigen2 = 0).
    //     The comparison is restricted to i in [8, N_ice-4]: near the west
    //     Dirichlet boundary (i=0,1) the velocity ghost at i=-1 is stale
    //     (non-periodic boundary, not updated by DMLocalToLocal) and the
    //     centered difference of the sharply curved u(x) has large truncation
    //     error, exactly as in the model; the same holds near the front.
    {
      double max_rel_err = 0.0, max_abs_e2 = 0.0;
      double west_edge_rel_err = 0.0, front_rel_err = 0.0;
      int n_checked = 0;
      array::AccessScope list{&geometry.ice_thickness, &geometry.cell_type, &eigen_strain};
      const double x_min = grid->x(0);
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if (not geometry.cell_type.floating_ice(i, j)) {
          continue;
        }
        // skip margin cells (one-sided differences)
        if (geometry.cell_type.ice_margin(i, j)) {
          continue;
        }
        const double H = H_exact(V0, H0, C, grid->x(i) - x_min);
        const double du = u_x_analytic(C, H);
        const double rel_err = std::fabs(eigen_strain(i, j).eigen1 - du) / du;
        if (i == 1) {
          west_edge_rel_err = rel_err;
        }
        if (i == N_ice - 2) {
          front_rel_err = rel_err;
        }
        if (i < 8 or i > N_ice - 4) {
          continue;
        }
        max_rel_err = std::max(max_rel_err, rel_err);
        max_abs_e2 = std::max(max_abs_e2, std::fabs(eigen_strain(i, j).eigen2));
        n_checked += 1;
      }
      std::cout << "STRAIN RATE (exact-velocity run, interior floating cells "
                   "i in [8, N_ice-4]):\n"
                << "  cells checked = " << n_checked << "\n"
                << "  max |eigen1 - C H^3| / (C H^3) = " << max_rel_err << "\n"
                << "  max |eigen2| = " << max_abs_e2 << " (must be ~0; eigen2 = A - q"
                   " = 0 exactly for the 1-D shelf)\n"
                << "  for reference: rel err at i=1 (west boundary, stale velocity"
                   " ghost + FD truncation) = " << west_edge_rel_err
                << ", at i=" << (N_ice - 2) << " (front, steep gradient) = "
                << front_rel_err << "\n";
    }

    // (b) eigen-calving: exact run zero everywhere; divergent run matches the
    //     reconstructed formula
    {
      double max_rate = 0.0;
      array::AccessScope list{&eigen_rate};
      for (auto p : grid->points()) {
        max_rate = std::max(max_rate, std::fabs(eigen_rate(p.i(), p.j())));
      }
      std::cout << "EIGEN CALVING (exact-velocity run): max |rate| = " << max_rate
                << " m/s (must be 0: eigen2 = 0, condition eigen2 > 0 fails)\n";
    }
    {
      const array::Scalar &rate_div = eigen_calving_div->calving_rate();
      const auto &strain_div = eigen_calving_div->strain_rates();
      double max_rate = 0.0;
      double max_consistency_err = 0.0;
      int n_nonzero = 0;
      array::AccessScope list{&rate_div, &geometry.cell_type, &strain_div};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        max_rate = std::max(max_rate, std::fabs(rate_div(i, j)));
        if (rate_div(i, j) != 0.0) {
          n_nonzero += 1;
        }
        double e1 = 0, e2 = 0;
        int N = 0;
        const double r_recon =
            eigen_calving_reconstructed(geometry.cell_type, strain_div,
                                        K_eigen, 2, i, j, &e1, &e2, &N, Mx, My);
        max_consistency_err = std::max(max_consistency_err, std::fabs(rate_div(i, j) - r_recon));
      }
      std::cout << "EIGEN CALVING (divergent-velocity run):\n"
                << "  max |rate| = " << max_rate << " m/s = " << max_rate * to_m_yr << " m/yr\n"
                << "  cells with nonzero rate = " << n_nonzero << "\n"
                << "  max |component rate - reconstructed formula| = "
                << max_consistency_err << " (must be ~0)\n";
    }

    // (c) von Mises: front-cell rate and consistency of the reconstructed
    //     intermediates
    {
      double max_consistency_err = 0.0;
      double front_rate = 0.0, front_rate_yr = 0.0;
      int n_front = 0;
      const double glen_exponent = flow_law->exponent();
      const double *z = enthalpy.levels().data();
      array::AccessScope list{&geometry.cell_type, &geometry.ice_thickness, &velocity,
                              &vonmises_strain, &vonmises_rate, &enthalpy};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if (not (geometry.cell_type.ice_free_ocean(i, j) and geometry.cell_type.next_to_ice(i, j)) or
            i < 2 or i > (int)Mx - 3) {
          continue;
        }
        double vel_sum = 0.0, hard_sum = 0.0, e1_sum = 0.0, e2_sum = 0.0;
        int N = 0;
        for (int p = -1; p < 2; p += 2) {
          const int I = i + p * 2;
          if (geometry.cell_type.icy(I, j)) {
            vel_sum += velocity(I, j).magnitude();
            const double H = geometry.ice_thickness(I, j);
            auto k = grid->kBelowHeight(H);
            hard_sum += averaged_hardness(*flow_law, H, k, z, enthalpy.get_column(I, j));
            e1_sum += vonmises_strain(I, j).eigen1;
            e2_sum += vonmises_strain(I, j).eigen2;
            N += 1;
          }
        }
        if (j >= 2 and j <= (int)My - 3) {
          for (int q = -1; q < 2; q += 2) {
            const int J = j + q * 2;
            if (geometry.cell_type.icy(i, J)) {
              vel_sum += velocity(i, J).magnitude();
              const double H = geometry.ice_thickness(i, J);
              auto k = grid->kBelowHeight(H);
              hard_sum += averaged_hardness(*flow_law, H, k, z, enthalpy.get_column(i, J));
              e1_sum += vonmises_strain(i, J).eigen1;
              e2_sum += vonmises_strain(i, J).eigen2;
              N += 1;
            }
          }
        }
        if (N == 0) {
          continue;
        }
        const double vel_avg = vel_sum / N;
        const double hard_avg = hard_sum / N;
        const double e1 = e1_sum / N, e2 = e2_sum / N;
        const double e_s = std::sqrt(0.5 * (std::pow(std::max(0.0, e1), 2) +
                                            std::pow(std::max(0.0, e2), 2)));
        const double sigma_tilde = std::sqrt(3.0) * hard_avg * std::pow(e_s, 1.0 / glen_exponent);
        const double rate_recon = vel_avg * sigma_tilde / sigma_max;
        max_consistency_err = std::max(max_consistency_err,
                                       std::fabs(vonmises_rate(i, j) - rate_recon));
        front_rate = vonmises_rate(i, j);
        front_rate_yr = vonmises_rate(i, j) * to_m_yr;
        n_front += 1;
      }
      std::cout << "VON MISES CALVING: front cells = " << n_front << "\n"
                << "  front-cell rate = " << front_rate << " m/s = " << front_rate_yr
                << " m/yr\n"
                << "  max |component rate - reconstructed| = " << max_consistency_err
                << " (must be ~0)\n";
    }

    // (d) Hayhurst: a positive rate on icy cells (water depth > 0); on
    //     ice-free cells next to ice the component (intentionally) propagates
    //     the mean of the icy neighbors' rates (used at the front by the
    //     retreat-rate code and for visualization), so "zero elsewhere" is not
    //     the contract -- verify the neighbor average instead.
    {
      double min_rate = std::numeric_limits<double>::infinity();
      double max_rate = -std::numeric_limits<double>::infinity();
      double max_propagation_err = 0.0;
      int n_icy = 0, n_propagated = 0;
      array::AccessScope list{&geometry.cell_type, &hayhurst_rate};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if (geometry.cell_type.icy(i, j)) {
          min_rate = std::min(min_rate, hayhurst_rate(i, j));
          max_rate = std::max(max_rate, hayhurst_rate(i, j));
          n_icy += 1;
        } else if (geometry.cell_type.ice_free(i, j) and geometry.cell_type.next_to_ice(i, j) and
                   i >= 1 and i <= (int)Mx - 2 and j >= 1 and j <= (int)My - 2) {
          // average of the icy neighbors (the second loop of HayhurstCalving::update)
          double sum = 0.0;
          int N = 0;
          for (int di = -1; di <= 1; di += 2) {
            if (geometry.cell_type.icy(i + di, j)) {
              sum += hayhurst_rate(i + di, j);
              N += 1;
            }
          }
          for (int dj = -1; dj <= 1; dj += 2) {
            if (geometry.cell_type.icy(i, j + dj)) {
              sum += hayhurst_rate(i, j + dj);
              N += 1;
            }
          }
          if (N > 0) {
            max_propagation_err =
                std::max(max_propagation_err, std::fabs(hayhurst_rate(i, j) - sum / N));
            n_propagated += 1;
          }
        }
      }
      std::cout << "HAYHURST CALVING: icy cells = " << n_icy << "\n"
                << "  rate range over icy cells = [" << min_rate << ", " << max_rate
                << "] m/s = [" << min_rate * to_m_yr << ", " << max_rate * to_m_yr << "] m/yr\n"
                << "  front cells with the propagated neighbor-average rate = "
                << n_propagated << "\n"
                << "  max |rate - mean(icy neighbors)| on those = " << max_propagation_err
                << " (must be ~0)\n";
    }

    // (e) frontal-melt kernels: monotonicity + clamps
    {
      bool mono_h = true, mono_q = true, mono_T = true;
      for (double h : h_grid) {
        for (double q : q_grid) {
          double prev = -std::numeric_limits<double>::infinity();
          for (double TF : T_grid) {
            const double qm = physics.frontal_melt_from_undercutting(h, q, TF);
            if (qm < prev) {
              mono_T = false;
            }
            prev = qm;
          }
        }
      }
      for (double q : q_grid) {
        for (double TF : T_grid) {
          double prev = -std::numeric_limits<double>::infinity();
          for (double h : h_grid) {
            const double qm = physics.frontal_melt_from_undercutting(h, q, TF);
            if (qm < prev) {
              mono_h = false;
            }
            prev = qm;
          }
        }
      }
      for (double h : h_grid) {
        for (double TF : T_grid) {
          double prev = -std::numeric_limits<double>::infinity();
          for (double q : q_grid) {
            const double qm = physics.frontal_melt_from_undercutting(h, q, TF);
            if (qm < prev) {
              mono_q = false;
            }
            prev = qm;
          }
        }
      }

      const double qm_clamp_h = physics.frontal_melt_from_undercutting(-5.0, 1.0, 1.0);
      const double qm_clamp_q = physics.frontal_melt_from_undercutting(100.0, -0.5, 1.0);
      const double qm_clamp_T = physics.frontal_melt_from_undercutting(100.0, 1.0, -1.0);
      const double qm_ismip6_clamp_h = physics.frontal_melt_from_ismip6(-5.0, 1.0, 1.0);
      const double qm_ref = physics.frontal_melt_from_undercutting(100.0, 1.0, 1.0);

      std::cout << "FRONTAL-MELT PHYSICS (undercutting kernel):\n"
                << "  monotone non-decreasing in h: " << (mono_h ? "OK" : "FAIL") << "\n"
                << "  monotone non-decreasing in q_sg: " << (mono_q ? "OK" : "FAIL") << "\n"
                << "  monotone non-decreasing in TF: " << (mono_T ? "OK" : "FAIL") << "\n"
                << "  q_m(h=-5) = " << qm_clamp_h << ", q_m(q_sg=-0.5) = " << qm_clamp_q
                << ", q_m(TF=-1) = " << qm_clamp_T << " (must all be 0: clamp)\n"
                << "  q_m(100, 1, 1) = " << qm_ref << " m/day\n"
                << "  ISMIP6 kernel (no clamp): q_m(h=-5) = " << qm_ismip6_clamp_h
                << " (must NOT be 0)\n";
    }

    // (f) Constant frontal-melt model
    {
      // recompute the maps for the summary from the dumped CSV logic: run again
      // cheaply for the two flag settings and report front-cell values
      array::Scalar1 water_flux(grid, "subglacial_water_flux");
      water_flux.set(0.0);
      const double melt_rate_config =
          config->get_number("frontal_melt.constant.melt_rate", "m day-1");

      for (bool include_floating : {false, true}) {
        config->set_flag("frontal_melt.include_floating_ice", include_floating);
        auto fm = std::make_shared<frontalmelt::Constant>(grid);
        fm->init(geometry);
        FrontalMeltInputs inputs;
        inputs.geometry = &geometry;
        inputs.subglacial_water_flux = &water_flux;
        fm->update(inputs, 0.0, units::convert(sys, 1.0, "year", "second"));

        double max_melt = 0.0, max_retreat = 0.0, front_retreat = 0.0;
        int n_melt_nonzero = 0;
        array::AccessScope list{&fm->frontal_melt_rate(), &fm->retreat_rate()};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          max_melt = std::max(max_melt, std::fabs(fm->frontal_melt_rate()(i, j)));
          max_retreat = std::max(max_retreat, std::fabs(fm->retreat_rate()(i, j)));
          if (fm->frontal_melt_rate()(i, j) != 0.0) {
            n_melt_nonzero += 1;
          }
          if (fm->retreat_rate()(i, j) != 0.0) {
            front_retreat = fm->retreat_rate()(i, j);
          }
        }
        std::cout << "CONSTANT FRONTAL MELT (include_floating_ice = "
                  << (include_floating ? "yes" : "no") << "):\n"
                  << "  max |melt rate| = " << max_melt * to_m_day << " m/day"
                  << " (" << n_melt_nonzero << " nonzero cells)\n"
                  << "  max |retreat rate| = " << max_retreat * to_m_day << " m/day\n";
        if (max_melt != 0.0) {
          const double ice_rho = config->get_number("constants.ice.density");
          const double water_rho = config->get_number("constants.sea_water.density");
          std::cout << "  front-cell retreat rate = " << front_retreat * to_m_day
                    << " m/day (must equal (rho_i/rho_w) * " << melt_rate_config
                    << " = " << (ice_rho / water_rho) * melt_rate_config << " m/day)\n";
        }
      }
    }

    // velocity error report (same computation as SSATestCase::report)
    {
      double maxvecerr = 0.0, avvecerr = 0.0, avuerr = 0.0, avverr = 0.0;
      double maxuerr = 0.0, maxverr = 0.0, exactvelmax = 0.0;

      array::AccessScope list{&velocity};
      const double x_min = grid->x(0);

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if (i >= N_ice) {
          continue;
        }
        const double x = grid->x(i);
        const double uexact = u_exact(V0, H0, C, x - x_min);
        const double vexact = 0.0;

        exactvelmax = std::max(exactvelmax, std::sqrt(uexact * uexact + vexact * vexact));

        const double uerr = std::fabs(velocity(i, j).u - uexact);
        const double verr = std::fabs(velocity(i, j).v - vexact);
        avuerr += uerr;
        avverr += verr;
        maxuerr = std::max(maxuerr, uerr);
        maxverr = std::max(maxverr, verr);
        avvecerr += std::sqrt(uerr * uerr + verr * verr);
        maxvecerr = std::max(maxvecerr, std::sqrt(uerr * uerr + verr * verr));
      }

      const unsigned int N = N_ice * grid->My();
      const double gexactvelmax = GlobalMax(grid->com, exactvelmax);
      const double gmaxuerr = GlobalMax(grid->com, maxuerr);
      const double gmaxverr = GlobalMax(grid->com, maxverr);
      const double gavuerr = GlobalSum(grid->com, avuerr) / N;
      const double gavverr = GlobalSum(grid->com, avverr) / N;
      const double gmaxvecerr = GlobalMax(grid->com, maxvecerr);
      const double gavvecerr = GlobalSum(grid->com, avvecerr) / N;

      std::cout << "NUMERICAL ERRORS in velocity relative to exact solution:\n";
      std::cout << "  max vector err = "
                << units::convert(sys, gmaxvecerr, "m second^-1", "m year^-1")
                << " m/year, av vector err = "
                << units::convert(sys, gavvecerr, "m second^-1", "m year^-1")
                << " m/year\n";
      std::cout << "  max u err = " << units::convert(sys, gmaxuerr, "m second^-1", "m year^-1")
                << " m/year, max v err = "
                << units::convert(sys, gmaxverr, "m second^-1", "m year^-1") << " m/year\n";
      std::cout << "  average percent error = " << (gavvecerr / gexactvelmax) * 100.0 << "%\n";
    }

    std::cout << "Calving / frontal-melt instrumentation complete. Dumps in " << out_dir
              << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
