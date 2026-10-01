// Stage 1 instrumentation driver for the PISM bed-deformation boundary
// (bed::PointwiseIsostasy, PISM config bed_deformation.model = "iso").
//
// The pointwise isostasy model is the *instantaneous* local isostasy model:
//
//   update_impl(load, t, dt):
//     topg_{n+1} = topg_last - f * (load - load_last),   f = rho_ice / rho_mantle
//     load_last  = load
//
// (see stage1/pism/src/earth/PointwiseIsostasy.cc).  It has NO relaxation
// time: a load change is fully compensated in the single update that follows
// it.  The exponential approach to isostatic equilibrium belongs to the
// Lingle-Clark model (mantle_viscosity + elastic lithosphere), which is NOT
// instrumented here (deferred subassembly; see stage1/dumps/README.md).
//
// The load passed to update_impl() is the *time-averaged* load over the update
// interval: BedDef::update() accumulates load * dt into m_load_accumulator on
// every call and only when t_final reaches the next update time (m_t_last +
// bed_deformation.update_interval) computes
//
//   load = m_load_accumulator / dt_beddef,  dt_beddef = t_final - m_t_last,
//
// then calls update_impl() and sets
//
//   uplift = (topg - topg_last) / dt_beddef,  topg_last = topg.
//
// The driver runs TWO full update intervals of the model on the verification
// test-F exact geometry (flat bed, thickness = exact FG H, Mx x My = 31 x 31,
// Lx = Ly = 900 km), with a deterministic load history:
//
//   bed_deformation.update_interval  = 100 "365day years" (3.1536e9 s)
//   dt_call                          = 20 "365day years" per update() call
//                                      (interval/5), so each interval is
//                                      exactly tiled by 5 calls and the
//                                      accumulator sees a clean time average.
//   Interval 1 (t = 0 -> 100 yr): the ice load grows linearly in 4 equal
//       steps from the test-F load H0(r) to H0(r) + dH, dH = 300 m (uniform
//       thickening applied only where ice exists, r <= LforFG - 1).  The
//       loads used by the 5 calls are H0 + k*dH/4, k = 0..4, so the
//       time-averaged load is exactly H0 + dH/2 = H0 + 150 m.
//   Interval 2 (t = 100 -> 200 yr): the load is constant at H0 + dH, so the
//       time-averaged load is exactly H0 + 300 m = load_last of interval 1
//       plus another 150 m step.
//
// Because f = rho_ice/rho_mantle = 910/3300 and the response is linear, the
// analytic results per ice cell are (sea level = 0, bed = 0):
//
//   bed_out1 = -f * 150 m            = -41.363636... m     (after interval 1)
//   bed_out2 = -f * 300 m            = -82.727272... m     (after interval 2)
//   uplift1 = uplift2 = -41.3636 m / 100 "365day yr"  (~ -0.414 m/365day-yr)
//   ice-free far-field cells: load = 0 always, bed stays 0.
//
// Every dumped bed value must satisfy  bed_out = bed_in - f*(load - load_last)
// to machine precision (that is the pointwise I-boundary contract) and
// uplift == (bed_out - bed_in)/dt_beddef (the base-class contract).
//
// Dumps (stage1/dumps/bed/):
//
//   meta.txt        - PISM revision, config file, grid, dt, command line
//   parameters.csv  - bed_deformation.*, constants.*, grid.*,
//                     time_stepping.*, bootstrapping.* config keys
//   inputs.nc       - geometry + the full state-in/state-out trace for both
//                     update intervals (loads, bed in/out, uplift)
//   columns.csv     - per cell: (i,j,x,y,r,H0,H1, load_avg_1, load_avg_2,
//                     bed_in_1, bed_out_1, uplift_1, bed_in_2, bed_out_2,
//                     uplift_2, bed_analytic_1, bed_analytic_2)
//   pointwise.csv   - the two pointwise laws sampled over small grids
//                     (compute_load and the topg update law) -- pure .esm
//                     test targets
//   load_history.csv- the deterministic load sequence (call, interval, t, dt,
//                     dH) that produced the trace
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/bed.json

#include <petsc.h>
#include <array>
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
#include "pism/earth/BedDef.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM bed-deformation boundary.\n"
  "Runs two update intervals of bed::PointwiseIsostasy on the test-F exact\n"
  "geometry with a deterministic load history and dumps the I/O trace to\n"
  "stage1/dumps/bed/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

//! Subclass exposing the (protected) bed-deformation state so we can dump the
//! load accumulator, the time-averaged load, the old load, and the bed
//! elevation at the start of each interval, without touching the vendored source.
class InstrumentedPointwiseIsostasy : public bed::PointwiseIsostasy {
public:
  InstrumentedPointwiseIsostasy(std::shared_ptr<const Grid> g)
      : bed::PointwiseIsostasy(g) {}

  // time-averaged load that update_impl() just consumed (m_load after an update)
  const array::Scalar &load() const { return m_load; }
  // load at the time of the last update (the update_impl input m_load_last)
  const array::Scalar &load_last() const { return m_load_last; }
  // the running load accumulator (reset to 0 after each update)
  const array::Scalar &load_accumulator() const { return m_load_accumulator; }
  // bed elevation at the start of the current interval
  const array::Scalar &topg_last() const { return m_topg_last; }
  // time of the last update
  double t_last() const { return m_t_last; }
  double update_interval_seconds() const { return m_update_interval; }
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
         "Configuration parameters affecting the bed-deformation component "
         "(bed_deformation.*, constants.*, grid.*, time_stepping.*, bootstrapping.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {"bed_deformation.", "constants.", "grid.",
                                          "time_stepping.", "bootstrapping."};

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

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "bed_def_test");
    auto config = ctx->config();

    // the model we are instrumenting: pointwise isostasy
    config->set_string("bed_deformation.model", "iso");

    // one full update interval: 100 "365day years" (the default is 10).  Each
    // update() call below is dt = interval/5 = 20 "365day years", so 5 calls
    // exactly tile each interval and the accumulator is a clean time average.
    config->set_number("bed_deformation.update_interval", 100.0);

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/bed";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    const unsigned int Mx = config->get_number("grid.Mx");
    const unsigned int My = config->get_number("grid.My");

    // grid: test-F domain (2-D, same as the hydrology/basal-strength drivers)
    auto grid = Grid::Shallow(ctx, 900.0e3, 900.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::NOT_PERIODIC);

    // ------------------------------------------------------------------
    // Scenario parameters (deterministic, documented in stage1/runs/bed.json)
    // ------------------------------------------------------------------
    const double interval_s = config->get_number("bed_deformation.update_interval", "seconds");
    const unsigned int calls_per_interval = 5;
    const double dt_call = interval_s / calls_per_interval; // 20 "365day years"
    const double year365_s = 365.0 * 86400.0;               // one "365day year"
    const unsigned int n_intervals = 2;

    const double dH = 300.0; // m, total uniform thickening over interval 1
    const double rho_ice = config->get_number("constants.ice.density"),
                 rho_mantle = config->get_number("bed_deformation.mantle_density");
    const double f = rho_ice / rho_mantle; // isostatic compensation factor

    // fraction of dH applied at call k of interval 1: k = 0..4 -> 0, 1/4, 1/2, 3/4, 1
    auto dH_k_frac = [](unsigned int k) { return 0.25 * (double)k; };

    // ------------------------------------------------------------------
    // meta + parameters
    // ------------------------------------------------------------------
    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "component", "bed deformation (PointwiseIsostasy, model = 'iso')");
      header(meta, "description",
             "Two full update intervals of bed::PointwiseIsostasy on the test-F "
             "exact geometry (flat bed, H = exact FG). Deterministic load "
             "history: interval 1 (t = 0..100 '365day yr') grows the ice load "
             "uniformly from H0(r) to H0(r)+300 m in 5 equal steps (calls of "
             "dt = 20 '365day yr'), so the time-averaged load is exactly "
             "H0 + 150 m; interval 2 holds the load constant at H0 + 300 m. "
             "PointwiseIsostasy has NO relaxation time: it responds "
             "instantaneously to the load change with bed_out = bed_in - "
             "f*(load_avg - load_last), f = rho_ice/rho_mantle. Expected "
             "bed_out_1 = -f*150 m, bed_out_2 = -f*300 m (ice cells).");
      header(meta, "dt_seconds_per_call", std::to_string(dt_call));
      header(meta, "dt_365day_years_per_call", std::to_string(dt_call / year365_s));
      header(meta, "calls_per_interval", std::to_string(calls_per_interval));
      header(meta, "update_interval_seconds", std::to_string(interval_s));
      header(meta, "update_interval_365day_years",
             std::to_string(interval_s / year365_s));
      header(meta, "n_intervals", std::to_string(n_intervals));
      header(meta, "total_calls", std::to_string(calls_per_interval * n_intervals));
      header(meta, "Mx", std::to_string(Mx));
      header(meta, "My", std::to_string(My));
      header(meta, "Lx", "900000");
      header(meta, "Ly", "900000");
      header(meta, "dx", std::to_string(grid->dx()));
      header(meta, "cell_area", std::to_string(grid->cell_area()));
      header(meta, "rho_ice", std::to_string(rho_ice));
      header(meta, "rho_mantle", std::to_string(rho_mantle));
      header(meta, "f_rho_ice_over_rho_mantle", std::to_string(f));
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    // ------------------------------------------------------------------
    // Test-F exact geometry (initial state: bed = 0, sea level = 0,
    // H0 = exact FG thickness; the load on the bed is H0, ice-equivalent)
    // ------------------------------------------------------------------
    const double LforFG = 750000.0;
    array::Scalar1 H0(grid, "ice_thickness_initial");
    array::Scalar1 sea_level(grid, "sea_level_elevation");
    sea_level.set(0.0);
    {
      array::AccessScope list{&H0};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        H0(i, j) = (r > LforFG - 1.0) ? 0.0 : exactFG(0.0, r, grid->z(), 0.0).H;
      }
    }

    // load fields passed to update() at each call: H = H0 + dH_k on ice cells
    array::Scalar1 H(grid, "ice_thickness");
    auto set_load = [&](double dH_k) {
      array::AccessScope list{&H0, &H};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        H(i, j) = (H0(i, j) > 0.0) ? (H0(i, j) + dH_k) : 0.0;
      }
    };

    // ------------------------------------------------------------------
    // Model: init from the test-F state (bed = 0, load_last = H0), then run
    // two intervals of 5 calls each, capturing the state in/out of each
    // interval (the model state after an interval is also the state-in of
    // the next).
    // ------------------------------------------------------------------
    auto model = std::make_shared<InstrumentedPointwiseIsostasy>(grid);
    model->init(InputOptions(INIT_OTHER, "", 0), H0, sea_level);

    // state-in copies
    array::Scalar1 load_last_0(grid, "bed_def_load_last_initial");
    array::Scalar1 bed_in_1(grid, "bed_elevation_in_1");
    array::Scalar1 load_avg_1(grid, "bed_def_load_avg_1"), load_avg_2(grid, "bed_def_load_avg_2"),
        bed_out_1(grid, "bed_elevation_out_1"), bed_in_2(grid, "bed_elevation_in_2"),
        uplift_1(grid, "bed_uplift_1"), uplift_2(grid, "bed_uplift_2");
    load_last_0.copy_from(model->load_last());
    bed_in_1.copy_from(model->bed_elevation());

    double t = 0.0;
    std::vector<std::array<double, 4>> load_history; // {t_start, dt, dH_k, interval}

    for (unsigned int interval = 1; interval <= n_intervals; ++interval) {
      for (unsigned int k = 0; k < calls_per_interval; ++k) {
        const double dH_k = (interval == 1) ? dH * dH_k_frac(k) : dH;
        set_load(dH_k);
        load_history.push_back({t, dt_call, dH_k, (double)interval});
        model->update(H, sea_level, t, dt_call);
        t += dt_call;
      }
      if (interval == 1) {
        // state-out of interval 1 == state-in of interval 2
        load_avg_1.copy_from(model->load());
        bed_out_1.copy_from(model->bed_elevation());
        uplift_1.copy_from(model->uplift());
        bed_in_2.copy_from(model->topg_last());
      } else {
        load_avg_2.copy_from(model->load());
        uplift_2.copy_from(model->uplift());
      }
    }

    // final bed (== bed_out_2)
    const array::Scalar &topg_out = model->bed_elevation();

    // ------------------------------------------------------------------
    // Analytic reference fields (the pointwise law applied to our inputs)
    // ------------------------------------------------------------------
    array::Scalar1 bed_analytic_1(grid, "bed_analytic_1"), bed_analytic_2(grid, "bed_analytic_2");
    {
      array::AccessScope list{&H0, &bed_in_1, &load_avg_1, &load_avg_2, &bed_out_1,
                              &bed_analytic_1, &bed_analytic_2};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        // interval 1: bed_in = 0, load_last = H0 (initial load)
        bed_analytic_1(i, j) = bed_in_1(i, j) - f * (load_avg_1(i, j) - H0(i, j));
        // interval 2: bed_in = bed_out_1, load_last = load_avg_1
        bed_analytic_2(i, j) = bed_out_1(i, j) - f * (load_avg_2(i, j) - load_avg_1(i, j));
      }
    }

    // ------------------------------------------------------------------
    // Dumps
    // ------------------------------------------------------------------
    {
      const std::vector<const array::Array *> vecs = {
          &topg_out,       &H,           &sea_level,    &H0,
          &load_last_0,    &load_avg_1,  &load_avg_2,   &bed_in_1,
          &bed_out_1,      &bed_in_2,    &uplift_1,     &uplift_2,
          &bed_analytic_1, &bed_analytic_2};
      write_netcdf(ctx, out_dir + "/inputs.nc", vecs);
    }

    // per-cell columns
    {
      std::ofstream out = open_output(out_dir + "/columns.csv");
      header(out, "description",
             "Two update intervals of bed::PointwiseIsostasy on the test-F "
             "exact geometry: per cell: geometry (H0 = initial test-F "
             "thickness, H1 = H0 + dH = load after the thickening), the "
             "time-averaged loads fed to update_impl (ice-equivalent m), the "
             "bed elevation in/out and uplift for both intervals, and the "
             "analytic reference bed from the pointwise law "
             "bed_out = bed_in - f*(load - load_last). The load, bed, and "
             "sea level are ice-equivalent meters / meters / meters; uplift "
             "is m/s. Ice-free far-field cells have load = 0 and stay at 0.");
      header(out, "columns",
             "i,j,x,y,r,H0,H1,load_avg_1,load_avg_2,bed_in_1,bed_out_1,uplift_1,"
             "bed_in_2,bed_out_2,uplift_2,bed_analytic_1,bed_analytic_2");

      array::AccessScope list{&H0, &H, &load_avg_1, &load_avg_2, &bed_in_1, &bed_out_1,
                              &uplift_1, &bed_in_2, &topg_out, &uplift_2,
                              &bed_analytic_1, &bed_analytic_2};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
            << H0(i, j) << "," << H(i, j) << "," << load_avg_1(i, j) << ","
            << load_avg_2(i, j) << "," << bed_in_1(i, j) << "," << bed_out_1(i, j) << ","
            << uplift_1(i, j) << "," << bed_in_2(i, j) << "," << topg_out(i, j) << ","
            << uplift_2(i, j) << "," << bed_analytic_1(i, j) << ","
            << bed_analytic_2(i, j) << "\n";
      }
      out.close();
    }

    // pointwise laws (analytic .esm test targets)
    {
      std::ofstream out = open_output(out_dir + "/pointwise.csv");
      header(out, "description",
             "Pointwise laws of the bed-deformation boundary (analytic .esm "
             "test targets). (1) compute_load(bed, ice_thickness, "
             "sea_level_elevation): the ice-equivalent load transmitted to the "
             "bed: ice_load = H, ocean_load = (rho_w/rho_i)*max(sea_level-bed,0), "
             "returns ice_load if ice_load > ocean_load else 0 (the load of "
             "floating ice is excluded). (2) update_law: topg_out = bed_in - "
             "f*(load - load_last), f = rho_ice/rho_mantle; load_last <- load "
             "(PointwiseIsostasy has no relaxation time).");
      header(out, "columns",
             "function,bed,H,sea_level,rho_ice,rho_mantle,ocean_depth,load,"
             "bed_in,load_last,f,bed_out,new_load_last");

      const double rho_w = config->get_number("constants.sea_water.density");
      const std::vector<double> beds = {-200.0, 0.0, 100.0, 500.0};
      const std::vector<double> Hs = {0.0, 100.0, 500.0, 2000.0, 3000.0};
      const std::vector<double> seas = {0.0, 50.0, 100.0};
      for (double bed : beds) {
        for (double Hv : Hs) {
          for (double sl : seas) {
            const double ice_load = Hv;
            const double ocean_depth = std::max(sl - bed, 0.0);
            const double ocean_load = (rho_w / rho_ice) * ocean_depth;
            const double load = ice_load > ocean_load ? ice_load : 0.0;
            out << "compute_load," << bed << "," << Hv << "," << sl << "," << rho_ice
                << "," << rho_mantle << "," << ocean_depth << "," << load << ",,,,,\n";
          }
        }
      }
      const std::vector<double> bed_ins = {-200.0, -100.0, 0.0, 100.0};
      const std::vector<double> loads = {0.0, 150.0, 300.0, 1000.0, 2000.0};
      const std::vector<double> last = {0.0, 150.0, 300.0};
      for (double bed_in : bed_ins) {
        for (double load : loads) {
          for (double load_last : last) {
            const double bed_out = bed_in - f * (load - load_last);
            // function, bed, H, sea_level, rho_ice, rho_mantle, ocean_depth,
            // load, bed_in, load_last, f, bed_out, new_load_last
            out << "update_law," << "" << "," << "" << "," << "" << "," << rho_ice << ","
                << rho_mantle << "," << "" << "," << load << "," << bed_in << ","
                << load_last << "," << f << "," << bed_out << "," << load << "\n";
          }
        }
      }
      out.close();
    }

    // load history
    {
      std::ofstream out = open_output(out_dir + "/load_history.csv");
      header(out, "description",
             "The deterministic load sequence fed to BedDef::update() (the "
             "driver passes H = H0(r) + dH on ice cells at each call). "
             "BedDef::update accumulates load*dt into m_load_accumulator and "
             "only triggers update_impl() when t_final reaches the next "
             "update time (m_t_last + update_interval), feeding it the "
             "time-averaged load m_load_accumulator/dt_beddef. Times are in "
             "seconds and '365day years' (1 yr365 = 31536000 s).");
      header(out, "columns", "call,interval,t_start_s,t_start_yr365,dt_s,dt_yr365,dH_m,t_end_s");
      for (size_t c = 0; c < load_history.size(); ++c) {
        const double t0 = load_history[c][0], dt = load_history[c][1], dH_k = load_history[c][2];
        const double iv = load_history[c][3];
        out << c + 1 << "," << iv << "," << t0 << "," << t0 / year365_s << "," << dt
            << "," << dt / year365_s << "," << dH_k << "," << t0 + dt << "\n";
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // Validation summaries
    // ------------------------------------------------------------------

    // (1) Pointwise contract: bed_out = bed_in - f*(load - load_last) for
    // both intervals, to machine precision.
    {
      double max_err1 = 0.0, max_err2 = 0.0, max_bed_farfield = 0.0;
      long ice_cells = 0, far_cells = 0;
      array::AccessScope list{&H0, &bed_in_1, &bed_out_1, &bed_in_2, &topg_out,
                              &load_avg_1, &load_avg_2, &uplift_1, &uplift_2,
                              &bed_analytic_1, &bed_analytic_2};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        max_err1 = std::max(max_err1, std::fabs(bed_out_1(i, j) - bed_analytic_1(i, j)));
        max_err2 = std::max(max_err2, std::fabs(topg_out(i, j) - bed_analytic_2(i, j)));
        if (H0(i, j) > 0.0) {
          ++ice_cells;
        } else {
          ++far_cells;
          max_bed_farfield = std::max(max_bed_farfield,
                                      std::fabs(topg_out(i, j)) + std::fabs(bed_out_1(i, j)));
        }
      }
      const double gmax_err1 = GlobalMax(grid->com, max_err1);
      const double gmax_err2 = GlobalMax(grid->com, max_err2);
      const double gmax_far = GlobalMax(grid->com, max_bed_farfield);
      const long gice = GlobalSum(grid->com, (int)ice_cells);
      const long gfar = GlobalSum(grid->com, (int)far_cells);

      std::cout << "BED DEFORMATION (PointwiseIsostasy) validation:\n"
                << "  f = rho_ice/rho_mantle = " << f << "\n"
                << "  update_interval = " << interval_s / year365_s << " '365day yr'"
                << " (" << interval_s << " s); dt_call = " << dt_call / year365_s
                << " '365day yr'; " << calls_per_interval << " calls/interval, "
                << n_intervals << " intervals\n"
                << "  max |bed_out_1 - analytic_1| = " << gmax_err1
                << " m (must be ~0; analytic_1 = -f*150 = " << -f * 150.0
                << " m under the ice)\n"
                << "  max |bed_out_2 - analytic_2| = " << gmax_err2
                << " m (must be ~0; analytic_2 = -f*300 = " << -f * 300.0
                << " m under the ice)\n"
                << "  ice cells = " << gice << ", ice-free far-field cells = " << gfar
                << "; max |bed| there = " << gmax_far << " m (must be 0)\n";
    }

    // (2) Uplift identity: uplift = (bed_out - bed_in)/dt_beddef per cell.
    {
      double max_uplift_err = 0.0, max_uplift1 = 0.0, max_uplift2 = 0.0;
      array::AccessScope list{&bed_in_1, &bed_out_1, &uplift_1, &bed_in_2, &topg_out,
                              &uplift_2};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double u1_expected = (bed_out_1(i, j) - bed_in_1(i, j)) / interval_s;
        const double u2_expected = (topg_out(i, j) - bed_in_2(i, j)) / interval_s;
        max_uplift_err = std::max(max_uplift_err,
                                  std::fabs(uplift_1(i, j) - u1_expected) +
                                      std::fabs(uplift_2(i, j) - u2_expected));
        max_uplift1 = std::max(max_uplift1, std::fabs(uplift_1(i, j)));
        max_uplift2 = std::max(max_uplift2, std::fabs(uplift_2(i, j)));
      }
      const double gmax_uplift_err = GlobalMax(grid->com, max_uplift_err);
      const double gmax_uplift1 = GlobalMax(grid->com, max_uplift1);
      const double gmax_uplift2 = GlobalMax(grid->com, max_uplift2);
      std::cout << "  max |uplift - (bed_out-bed_in)/dt_beddef| = " << gmax_uplift_err
                << " m/s (must be ~0)\n"
                << "  max |uplift_1| = " << gmax_uplift1 << " m/s ("
                << gmax_uplift1 * year365_s << " m/365day-yr; expected "
                << -f * 150.0 / 100.0 << " m/365day-yr)\n"
                << "  max |uplift_2| = " << gmax_uplift2 << " m/s ("
                << gmax_uplift2 * year365_s << " m/365day-yr; expected "
                << -f * 150.0 / 100.0 << " m/365day-yr)\n";
    }

    // (3) Load-accumulator semantics: the load fed to update_impl is the
    // time average of the loads passed to the calls over the interval.
    {
      double max_acc_err = 0.0;
      // recompute the interval-1 average independently from H0 and the dH history
      array::Scalar1 acc(grid, "acc");
      acc.set(0.0);
      {
        array::AccessScope l2{&H0, &acc};
        for (unsigned int k = 0; k < calls_per_interval; ++k) {
          const double dH_k = dH * dH_k_frac(k);
          for (auto p : grid->points()) {
            const int i = p.i(), j = p.j();
            acc(i, j) += ((H0(i, j) > 0.0) ? (H0(i, j) + dH_k) : 0.0) * dt_call;
          }
        }
      }
      acc.scale(1.0 / interval_s);
      array::AccessScope l3{&acc, &load_avg_1};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        max_acc_err = std::max(max_acc_err, std::fabs(acc(i, j) - load_avg_1(i, j)));
      }
      const double gmax_acc_err = GlobalMax(grid->com, max_acc_err);
      std::cout << "  load-accumulator check: max |recomputed_avg_1 - load_avg_1| = "
                << gmax_acc_err << " m (must be ~0; the load is the time "
                << "average over the interval, here exactly H0 + dH/2 = H0 + 150 m)\n";
    }

    // (4) Time-averaged load and bed values at the dome center
    {
      array::AccessScope list{&H0, &load_avg_1, &load_avg_2, &bed_out_1, &topg_out};
      const int ic = (int)Mx / 2, jc = (int)My / 2;
      std::cout << "  dome cell (" << ic << "," << jc << "): H0 = " << H0(ic, jc)
                << " m, load_avg_1 = " << load_avg_1(ic, jc)
                << " m (expect H0 + 150 = " << H0(ic, jc) + 150.0
                << "), load_avg_2 = " << load_avg_2(ic, jc)
                << " m (expect H0 + 300 = " << H0(ic, jc) + 300.0 << ")\n"
                << "  bed_out_1 = " << bed_out_1(ic, jc) << " m (expect "
                << -f * 150.0 << "), bed_out_2 = " << topg_out(ic, jc) << " m (expect "
                << -f * 300.0 << ")\n";
    }

    std::cout << "Bed deformation instrumentation complete. Dumps in " << out_dir << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
