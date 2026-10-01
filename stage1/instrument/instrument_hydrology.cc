// Stage 1 instrumentation driver for the PISM subglacial-hydrology boundary
// (hydrology::Routing, PISM's default hydrology.model = "routing").
//
// Runs ONE update step of the Routing (Shreve) subglacial hydrology model on
// the verification test-F exact geometry (flat bed, thickness = exact H,
// Mx x My = 31 x 31, Lx = Ly = 900 km), with prescribed deterministic forcing:
//
//   surface_input_rate  = 0.2 m/yr water-equivalent, constant under ice
//                         (passed in the model's units, kg m^-2 s^-1)
//   basal_melt_rate     = 0.1 m/yr * exp(-(r/200 km)^2)  (dome-centered blob,
//                         ice-equivalent m/s)
//   ice_sliding_speed   = 100 m/yr, constant (Routing *ignores* it -- it is a
//                         Distributed-model input; dumped because it is part of
//                         the Hydrology Inputs contract)
//   no_model_mask       = 0 everywhere (whole domain modeled)
//   W_till (state-in)   = 0 everywhere (dry till; absorbs the input over dt)
//   W      (state-in)   = 0.5 m * exp(-(r/150 km)^2)  (transportable-water blob
//                         whose drainage drives the flux trace)
//   dt                  = 0.25 yr (3 months; the model takes internal substeps
//                         capped by max_timestep_W_cfl() at ~5 days -- dominated
//                         by its eps = 1e-6 regularization at these velocities,
//                         so 19 substeps over the step)
//
// All forcing is radially symmetric and the geometry is the radially symmetric
// test-F dome, so the flux field must be exactly centrally anti-symmetric and
// the water balance must close to round-off.
//
// Dumps (stage1/dumps/hydrology/):
//
//   meta.txt        - PISM revision, config file, grid, dt, command line
//   parameters.csv  - hydrology.*, constants.*, grid.* config keys
//   inputs.nc       - geometry + all inputs + state in/out + outputs
//   columns.csv     - per cell: (i,j,x,y,r,H,mask, surface_input_rate,
//                     basal_melt_rate, sliding_speed, W_till_in, W_till_out,
//                     W_in, W_out, overburden_pressure, flux_u, flux_v)
//   flux_law.csv    - pointwise Routing flux law: K = k W^(alpha-1) |grad R|^(beta-2)
//                     and |q| = K W |grad R| sampled over (W, |grad R|) grids
//                     (analytic .esm test target)
//   staggered.csv   - internal subassembly fields from the last internal substep:
//                     Wstag, Kstag, Vstag (u,v), Qstag (u,v) and the accumulated
//                     Qstag_average (u,v) at cell edges
//
// The dumped flux() is the *advective* flux q = V W_upwind (V = -K grad R with
// R = P + rho_w g b).  The diffusive part rho_w g K W grad W is applied
// separately inside the W update (W_change_due_to_flow) and is NOT part of the
// dumped flux field.
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/hydrology.json

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
#include "pism/util/Grid.hh"
#include "pism/util/Mask.hh"
#include "pism/util/Time.hh"
#include "pism/util/Units.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/util/array/Vector.hh"
#include "pism/util/array/Staggered.hh"
#include "pism/hydrology/Routing.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM subglacial-hydrology boundary.\n"
  "Runs one update step of hydrology::Routing on the test-F exact geometry and\n"
  "dumps the I/O trace to stage1/dumps/hydrology/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

//! Subclass exposing Routing internals so we can dump the state and the
//! internal (staggered) subassembly fields without touching the vendored source.
class InstrumentedRouting : public hydrology::Routing {
public:
  InstrumentedRouting(std::shared_ptr<const Grid> g) : hydrology::Routing(g) {}

  // state (m_W is protected in Hydrology; exposed here for the W_in copy)
  const array::Scalar &W() const { return m_W; }

  // the model's converted input fields (m/s water-equivalent)
  const array::Scalar &model_surface_input() const { return m_surface_input_rate; }
  const array::Scalar &model_basal_melt() const { return m_basal_melt_rate; }

  // internal subassembly fields (last internal substep)
  const array::Staggered &W_stag() const { return m_Wstag; }
  const array::Staggered &K_stag() const { return m_Kstag; }
  const array::Staggered &V_stag() const { return m_Vstag; }
  const array::Staggered &Q_stag() const { return m_Qstag; }
  const array::Staggered &Q_stag_average() const { return m_Qstag_average; }

  // the bed surface the model routes water over (ice_bottom_surface)
  const array::Scalar &bottom_surface() const { return m_bottom_surface; }

  // hydrology-recommended max time step (from the last-substep state)
  double max_dt_cfl() const { return max_timestep_W_cfl(); }
  double max_dt_diff(double KW_max) const { return max_timestep_W_diff(KW_max); }
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
         "Configuration parameters affecting the hydrology::Routing component "
         "(hydrology.*, constants.*, grid.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {"hydrology.", "constants.", "grid."};

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

} // namespace
} // namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::hydrology;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "hydrology_test");
    auto config = ctx->config();

    // the model we are instrumenting
    config->set_string("hydrology.model", "routing");

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/hydrology";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    const unsigned int Mx = config->get_number("grid.Mx");
    const unsigned int My = config->get_number("grid.My");

    const double dt_years = options::Real(ctx->unit_system(), "-dt_years",
                                          "Time step length, in years", "years", 1.0);
    const double dt = units::convert(ctx->unit_system(), dt_years, "years", "seconds");

    // grid: test-F domain (2-D, same as the basal-strength driver)
    auto grid = Grid::Shallow(ctx, 900.0e3, 900.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::NOT_PERIODIC);

    // ------------------------------------------------------------------
    // config keys that enter the Routing physics
    // ------------------------------------------------------------------
    const double ice_rho = config->get_number("constants.ice.density"),
                 rho_w = config->get_number("constants.fresh_water.density"),
                 g_std = config->get_number("constants.standard_gravity"),
                 k_cond = config->get_number("hydrology.hydraulic_conductivity"),
                 alpha = config->get_number("hydrology.thickness_power_in_flux"),
                 beta = config->get_number("hydrology.gradient_power_in_flux"),
                 W_till_max = config->get_number("hydrology.tillwat_max"),
                 C_decay = config->get_number("hydrology.tillwat_decay_rate", "m / second");
    const double rg = rho_w * g_std;

    // ------------------------------------------------------------------
    // meta + parameters
    // ------------------------------------------------------------------
    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "component", "hydrology (Routing, Shreve: q = -K grad psi)");
      header(meta, "description",
             "One dt step of hydrology::Routing on the test-F exact geometry "
             "(flat bed, H = exact FG) with prescribed deterministic forcing: "
             "constant surface input rate 0.2 m/yr, dome-centered basal melt "
             "blob 0.1 m/yr * exp(-(r/200km)^2), constant sliding speed 100 "
             "m/yr (ignored by Routing), W_till_in = 0, W_in = 0.5 m * "
             "exp(-(r/150km)^2) transportable-water blob");
      header(meta, "dt_seconds", std::to_string(dt));
      header(meta, "dt_years", std::to_string(dt_years));
      header(meta, "Mx", std::to_string(Mx));
      header(meta, "My", std::to_string(My));
      header(meta, "Lx", "900000");
      header(meta, "Ly", "900000");
      header(meta, "dx", std::to_string(grid->dx()));
      header(meta, "cell_area", std::to_string(grid->cell_area()));
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    // ------------------------------------------------------------------
    // Test-F exact geometry
    // ------------------------------------------------------------------
    const double LforFG = 750000.0;
    Geometry geometry(grid);
    geometry.sea_level_elevation.set(0.0);
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
    geometry.ensure_consistency(config->get_number("geometry.ice_free_thickness_standard"));

    // ------------------------------------------------------------------
    // Inputs: deterministic forcing fields
    // ------------------------------------------------------------------
    // no_model_mask: whole domain modeled
    array::Scalar1 no_model_mask(grid, "no_model_mask");
    no_model_mask.set(0.0);

    // surface_input_rate in the model's units (kg m^-2 s^-1): 0.2 m/yr of water
    const double surface_input_kgm2s =
        rho_w * units::convert(ctx->unit_system(), 0.2, "m / year", "m / second");
    array::Scalar1 surface_input_rate(grid, "surface_input_rate");
    surface_input_rate.set(surface_input_kgm2s);

    // basal melt rate in ice-equivalent m/s: dome-centered blob, 0.1 m/yr peak
    const double bmelt_scale = units::convert(ctx->unit_system(), 0.1, "m / year", "m / second");
    const double bmelt_r0 = 200.0e3;
    array::Scalar1 basal_melt_rate(grid, "basal_melt_rate");
    {
      array::AccessScope list{&basal_melt_rate};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        basal_melt_rate(i, j) = bmelt_scale * std::exp(-(r / bmelt_r0) * (r / bmelt_r0));
      }
    }

    // ice sliding speed in m/s: constant 100 m/yr (not used by Routing)
    const double sliding_speed = units::convert(ctx->unit_system(), 100.0, "m / year", "m / second");
    array::Scalar1 ice_sliding_speed(grid, "ice_sliding_speed");
    ice_sliding_speed.set(sliding_speed);

    // ------------------------------------------------------------------
    // Initial state: W_till = 0 (dry till), W = 0.5 m exp(-(r/150km)^2)
    // ------------------------------------------------------------------
    const double W0 = 0.5, W_r0 = 150.0e3;
    array::Scalar1 W_till(grid, "till_water_thickness");
    array::Scalar1 W(grid, "subglacial_water_thickness");
    W_till.set(0.0);
    {
      array::AccessScope list{&W};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        W(i, j) = W0 * std::exp(-(r / W_r0) * (r / W_r0));
      }
    }
    W.update_ghosts();

    // ------------------------------------------------------------------
    // Model: init state, copy state-in, run one step
    // ------------------------------------------------------------------
    auto model = std::make_shared<InstrumentedRouting>(grid);

    // overburden pressure field for init (ignored by init_impl, but part of
    // the Hydrology::init contract)
    array::Scalar1 P0(grid, "overburden_pressure");
    P0.set(0.0);

    model->init(W_till, W, P0);

    // state-in copies (self-contained trace)
    array::Scalar1 W_till_in(grid, "W_till_in"), W_in(grid, "W_in");
    W_till_in.copy_from(model->till_water_thickness());
    W_in.copy_from(model->subglacial_water_thickness());

    Inputs hinputs;
    hinputs.no_model_mask = &no_model_mask;
    hinputs.geometry = &geometry;
    hinputs.surface_input_rate = &surface_input_rate;
    hinputs.basal_melt_rate = &basal_melt_rate;
    hinputs.ice_sliding_speed = &ice_sliding_speed;

    model->update(0.0, dt, hinputs);

    const array::Scalar &W_till_out = model->till_water_thickness();
    const array::Scalar &W_out = model->subglacial_water_thickness();
    const array::Scalar &P_over = model->overburden_pressure();
    const array::Vector &Q = model->flux();
    const array::Scalar &surf_m = model->model_surface_input();
    const array::Scalar &basal_m = model->model_basal_melt();

    // ------------------------------------------------------------------
    // Dumps
    // ------------------------------------------------------------------
    {
      const std::vector<const array::Array *> vecs = {
          &geometry.ice_surface_elevation, &geometry.ice_thickness,
          &geometry.bed_elevation,         &geometry.cell_type,
          &geometry.sea_level_elevation,   &no_model_mask,
          &surface_input_rate,             &basal_melt_rate,
          &ice_sliding_speed,              &surf_m,
          &basal_m,                        &W_till_in,
          &W_till_out,                     &W_in,
          &W_out,                          &P_over,
          &Q};
      write_netcdf(ctx, out_dir + "/inputs.nc", vecs);
    }

    // per-cell columns
    {
      std::ofstream out = open_output(out_dir + "/columns.csv");
      header(out, "description",
             "One dt step of hydrology::Routing on the test-F exact geometry: per "
             "cell: geometry (H, mask), the raw forcing inputs (surface_input_rate "
             "in kg m^-2 s^-1, basal_melt_rate in ice-equivalent m/s, sliding "
             "speed in m/s), the model's converted input fields (m/s "
             "water-equivalent), the state before/after (W_till, W in m), the "
             "overburden pressure, and the dumped advective flux components "
             "(m^2 s^-1)");
      header(out, "columns",
             "i,j,x,y,r,H,mask,surface_input_rate,basal_melt_rate,sliding_speed,"
             "surface_input_m_s,basal_melt_m_s,W_till_in,W_till_out,W_in,W_out,"
             "overburden_pressure,flux_u,flux_v");

      array::AccessScope list{&geometry.ice_thickness, &geometry.cell_type,
                              &surface_input_rate, &basal_melt_rate, &ice_sliding_speed,
                              &surf_m, &basal_m, &W_till_in, &W_till_out, &W_in,
                              &W_out, &P_over, &Q};

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
            << geometry.ice_thickness(i, j) << "," << (int)geometry.cell_type(i, j) << ","
            << surface_input_rate(i, j) << "," << basal_melt_rate(i, j) << ","
            << ice_sliding_speed(i, j) << "," << surf_m(i, j) << "," << basal_m(i, j) << ","
            << W_till_in(i, j) << "," << W_till_out(i, j) << "," << W_in(i, j) << ","
            << W_out(i, j) << "," << P_over(i, j) << "," << Q(i, j).u << "," << Q(i, j).v
            << "\n";
      }
      out.close();
    }

    // pointwise Routing flux law:
    //   K    = k W^(alpha-1) (|grad R|^2 + eps^2)^((beta-2)/2),  eps = (beta<2) ? 1 : 0
    //   |q|  = K W |grad R|
    // (the exact staggered-grid formula from compute_conductivity())
    {
      std::ofstream out = open_output(out_dir + "/flux_law.csv");
      header(out, "description",
             "Pointwise Routing flux law (analytic .esm test target): "
             "K = k W^(alpha-1) (G^2 + eps^2)^((beta-2)/2) with eps = 1 for "
             "beta < 2 (the regularization used in compute_conductivity), "
             "G = |grad R|, R = P + rho_w g b; |q| = K W G is the advective "
             "flux magnitude of the continuous law q = -K W grad R.");
      header(out, "columns", "W,G,alpha,beta,k,K,q_mag");

      const std::vector<double> Ws = {0.0, 0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5, 2.0, 4.0};
      const std::vector<double> Gs = {0.0, 0.1, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0};
      const double eps = (beta < 2.0) ? 1.0 : 0.0;
      for (double Wv : Ws) {
        for (double G : Gs) {
          const double K =
              k_cond * std::pow(Wv, alpha - 1.0) * std::pow(G * G + eps * eps, (beta - 2.0) / 2.0);
          out << Wv << "," << G << "," << alpha << "," << beta << "," << k_cond << ","
              << K << "," << K * Wv * G << "\n";
        }
      }
      out.close();
    }

    // staggered subassembly fields (last internal substep)
    {
      std::ofstream out = open_output(out_dir + "/staggered.csv");
      header(out, "description",
             "Internal Routing subassembly fields from the last internal "
             "substep, on the staggered (edge-centered) grid: (i,j,o=0) is the "
             "east edge of cell (i,j) at x + dx/2, (i,j,o=1) the north edge at "
             "y + dy/2. Wstag = edge water thickness, Kstag = nonlinear "
             "conductivity, Vstag = edge velocity V = -K grad R, Qstag = "
             "advective flux V*W_upwind, Qstag_average = accumulated "
             "sum(hdt*Qstag) over the whole step.");
      header(out, "columns", "i,j,o,x,y,Wstag,Kstag,Vstag,Qstag,Qstag_avg");

      const array::Staggered &Ws = model->W_stag(), &Ks = model->K_stag(),
                             &Vs = model->V_stag(), &Qs = model->Q_stag(),
                             &Qa = model->Q_stag_average();

      array::AccessScope list{&Ws, &Ks, &Vs, &Qs, &Qa};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        for (int o = 0; o < 2; ++o) {
          const double xo = grid->x(i) + (o == 0 ? 0.5 * grid->dx() : 0.0);
          const double yo = grid->y(j) + (o == 1 ? 0.5 * grid->dy() : 0.0);
          out << i << "," << j << "," << o << "," << xo << "," << yo << "," << Ws(i, j, o)
              << "," << Ks(i, j, o) << "," << Vs(i, j, o) << "," << Qs(i, j, o) << ","
              << Qa(i, j, o) << "\n";
        }
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // Validation summaries
    // ------------------------------------------------------------------

    // Water-mass conservation:
    //   the model's own accounting (kg) must balance:
    //     m_total_change = m_input_change + m_flow_change
    //                     + (grounded_margin + grounding_line + conservation_error
    //                        + no_model_mask)_change
    //   and an independent recomputation from the state fields must match:
    //     dV = dt*(surface + basal input volume) - dt*sum(div Q) - V_boundary
    {
      const double cell_area = grid->cell_area();
      const double kg_per_m = rho_w * cell_area;

      double sum_total = 0.0, sum_input = 0.0, sum_flow = 0.0, sum_boundary = 0.0;
      double sum_Win = 0.0, sum_Wtill_in = 0.0, sum_Wout = 0.0, sum_Wtill_out = 0.0;
      double sum_surf = 0.0, sum_basal = 0.0;

      const array::Scalar &d_total = model->mass_change();
      const array::Scalar &d_input = model->mass_change_due_to_input();
      const array::Scalar &d_flow = model->mass_change_due_to_lateral_flow();
      const array::Scalar &d_margin = model->mass_change_at_grounded_margin();
      const array::Scalar &d_gl = model->mass_change_at_grounding_line();
      const array::Scalar &d_cons = model->mass_change_due_to_conservation_error();
      const array::Scalar &d_nmm = model->mass_change_at_domain_boundary();

      array::AccessScope list{&d_total, &d_input, &d_flow, &d_margin, &d_gl, &d_cons,
                              &d_nmm, &W_till_in, &W_in, &W_till_out, &W_out, &surf_m,
                              &basal_m};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        sum_total += d_total(i, j);
        sum_input += d_input(i, j);
        sum_flow += d_flow(i, j);
        sum_boundary += d_margin(i, j) + d_gl(i, j) + d_cons(i, j) + d_nmm(i, j);
        sum_Win += W_in(i, j);
        sum_Wtill_in += W_till_in(i, j);
        sum_Wout += W_out(i, j);
        sum_Wtill_out += W_till_out(i, j);
        sum_surf += surf_m(i, j);
        sum_basal += basal_m(i, j);
      }

      const double gsum_total = GlobalSum(grid->com, sum_total);
      const double gsum_input = GlobalSum(grid->com, sum_input);
      const double gsum_flow = GlobalSum(grid->com, sum_flow);
      const double gsum_boundary = GlobalSum(grid->com, sum_boundary);
      const double gsum_Win = GlobalSum(grid->com, sum_Win);
      const double gsum_Wtill_in = GlobalSum(grid->com, sum_Wtill_in);
      const double gsum_Wout = GlobalSum(grid->com, sum_Wout);
      const double gsum_Wtill_out = GlobalSum(grid->com, sum_Wtill_out);
      const double gsum_surf = GlobalSum(grid->com, sum_surf);
      const double gsum_basal = GlobalSum(grid->com, sum_basal);

      const double V_in = dt * (gsum_surf + gsum_basal) * cell_area;                 // m^3
      const double dV = (gsum_Wout + gsum_Wtill_out - gsum_Win - gsum_Wtill_in) *
                        cell_area;                                                   // m^3
      // V_boundary < 0 when water leaves the system (margin/grounding-line/
      // conservation-error/no-model corrections), > 0 when water is added
      const double V_boundary = gsum_boundary / rho_w;                               // m^3

      // residual of the model's own accounting
      const double residual_model = gsum_total - (gsum_input + gsum_flow + gsum_boundary);
      // independent residual from the state fields:
      //   dV = V_in + V_flow + V_boundary,  V_flow ~ 0 (no flux through the
      //   outer boundary; the flow is internal and conservative)
      const double residual_state = dV - V_in - gsum_flow / rho_w - V_boundary;

      std::cout << "HYDROLOGY (Routing) mass-conservation validation:\n"
                << "  dt = " << dt << " s = " << dt_years << " years\n"
                << "  water volume in  (W_in + W_till_in)  = "
                << gsum_Win + gsum_Wtill_in << " m (sum) -> " << (gsum_Win + gsum_Wtill_in) * cell_area
                << " m^3\n"
                << "  water volume out (W_out + W_till_out) = "
                << gsum_Wout + gsum_Wtill_out << " m (sum) -> "
                << (gsum_Wout + gsum_Wtill_out) * cell_area << " m^3\n"
                << "  dV(state) = " << dV << " m^3\n"
                << "  input volume dt*(surf+basal) = " << V_in << " m^3\n"
                << "  boundary-correction volume = " << V_boundary
                << " m^3 (< 0: water left via margin/grounding-line/conservation-error/no-model corrections)\n"
                << "  residual (model accounting)  = " << residual_model << " kg"
                << " (= d(total) - input - flow - boundary; must be ~0)\n"
                << "  residual (independent state) = " << residual_state << " m^3"
                << " (= dV - input - flow - boundary; must be ~0)\n"
                << "  sum(flow change) = " << gsum_flow << " kg"
                << " (= -dt*sum(div q) * rho_w; ~0 when no flux leaves the outer boundary)\n";
    }

    // W_till non-negativity, W bounds, and overburden pressure check
    {
      double min_Wtill = std::numeric_limits<double>::infinity(),
             max_Wtill = -std::numeric_limits<double>::infinity(),
             min_W = std::numeric_limits<double>::infinity(),
             max_W = -std::numeric_limits<double>::infinity(),
             max_Perr = 0.0;
      long grounded_cells = 0;

      array::AccessScope list{&W_till_out, &W_out, &P_over, &geometry.ice_thickness,
                              &geometry.cell_type};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        min_Wtill = std::min(min_Wtill, W_till_out(i, j));
        max_Wtill = std::max(max_Wtill, W_till_out(i, j));
        min_W = std::min(min_W, W_out(i, j));
        max_W = std::max(max_W, W_out(i, j));
        if (geometry.cell_type.grounded_ice(i, j)) {
          ++grounded_cells;
          const double P_expected = ice_rho * g_std * geometry.ice_thickness(i, j);
          max_Perr = std::max(max_Perr, std::fabs(P_over(i, j) - P_expected));
        }
      }
      const double gmin_Wtill = GlobalMin(grid->com, min_Wtill),
                   gmax_Wtill = GlobalMax(grid->com, max_Wtill),
                   gmin_W = GlobalMin(grid->com, min_W),
                   gmax_W = GlobalMax(grid->com, max_W),
                   gmax_Perr = GlobalMax(grid->com, max_Perr);
      const long ggrounded = GlobalSum(grid->com, (int)grounded_cells);

      std::cout << "  W_till in range [" << gmin_Wtill << ", " << gmax_Wtill
                << "] m (must be >= 0)\n"
                << "  W in range [" << gmin_W << ", " << gmax_W << "] m\n"
                << "  max |overburden_pressure - rho_ice*g*H| = " << gmax_Perr
                << " Pa over " << ggrounded << " grounded cells (must be ~0)\n";
    }

    // flux direction: q = V W_upwind with V = -K grad R, R = P + rho_w g b
    //   => q . grad R <= 0 where W > 0
    {
      array::Scalar1 R(grid, "R");
      {
        array::AccessScope list{&P_over, &model->bottom_surface(), &R};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          R(i, j) = P_over(i, j) + rg * model->bottom_surface()(i, j);
        }
      }
      R.update_ghosts();

      double max_positive = 0.0, max_qmag = 0.0;
      long checked = 0;
      array::AccessScope list{&R, &Q, &geometry.cell_type};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if (not geometry.cell_type.grounded_ice(i, j)) {
          continue;
        }
        const double qmag = std::sqrt(Q(i, j).u * Q(i, j).u + Q(i, j).v * Q(i, j).v);
        max_qmag = std::max(max_qmag, qmag);
        if (qmag <= 0.0) {
          continue;
        }
        ++checked;
        const double dRdx = (R(i + 1, j) - R(i - 1, j)) / (2.0 * grid->dx());
        const double dRdy = (R(i, j + 1) - R(i, j - 1)) / (2.0 * grid->dy());
        const double dot = Q(i, j).u * dRdx + Q(i, j).v * dRdy;
        max_positive = std::max(max_positive, dot);
      }
      const double gmax_positive = GlobalMax(grid->com, max_positive);
      const double gmax_qmag = GlobalMax(grid->com, max_qmag);
      const long gchecked = GlobalSum(grid->com, (int)checked);
      std::cout << "  flux-direction check (q . grad R <= 0): max q.gradR = "
                << gmax_positive << " Pa m s^-1 over " << gchecked
                << " grounded cells with |q| > 0 (must be <= 0)\n"
                << "  max |q| = " << gmax_qmag << " m^2 s^-1 (= "
                << units::convert(ctx->unit_system(), gmax_qmag, "m^2 second^-1", "m^2 day^-1")
                << " m^2/day)\n";
    }

    // central antisymmetry of the flux: q(i,j) = -q(Mx-1-i, My-1-j) for a
    // radially symmetric forcing + geometry
    {
      double max_antisym = 0.0, max_q = 0.0;
      array::AccessScope list{&Q};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const int im = (int)Mx - 1 - i, jm = (int)My - 1 - j;
        const double u_m = Q(im, jm).u, v_m = Q(im, jm).v;
        const double au = Q(i, j).u + u_m, av = Q(i, j).v + v_m;
        max_antisym = std::max(max_antisym, std::sqrt(au * au + av * av));
        max_q = std::max(max_q, std::sqrt(Q(i, j).u * Q(i, j).u + Q(i, j).v * Q(i, j).v));
      }
      const double gmax_antisym = GlobalMax(grid->com, max_antisym);
      const double gmax_q = GlobalMax(grid->com, max_q);
      std::cout << "  central antisymmetry: max |q(i,j) + q(mirror)| = "
                << gmax_antisym << " m^2 s^-1 (max |q| = " << gmax_q
                << "; relative = " << gmax_antisym / gmax_q << ")\n";
    }

    // hydrology-recommended max time step from the last-substep state
    {
      double KW_max = 0.0;
      const array::Staggered &Ws = model->W_stag(), &Ks = model->K_stag();
      array::AccessScope list{&Ws, &Ks};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        for (int o = 0; o < 2; ++o) {
          KW_max = std::max(KW_max, Ks(i, j, o) * Ws(i, j, o));
        }
      }
      KW_max = GlobalMax(grid->com, KW_max);
      const double dt_cfl = model->max_dt_cfl();
      const double dt_diff = model->max_dt_diff(KW_max);
      std::cout << "  internal-substep limits (last substep state): dt_cfl = "
                << dt_cfl << " s, dt_diff = " << dt_diff << " s (KW_max = " << KW_max << ")\n";
    }

    std::cout << "Hydrology instrumentation complete. Dumps in " << out_dir << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
