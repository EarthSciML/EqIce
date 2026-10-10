// Stage-2 analysis driver for the PISM subglacial-hydrology boundary
// (hydrology::Routing, PISM v2.3.2). Stage-2 companion to
// stage1/instrument/instrument_hydrology.cc: re-runs the SAME one-step test-F
// configuration but replicates the internal Routing::update_impl() substep
// loop in a subclass so that the per-substep state can be dumped.
//
// The replication is a verbatim copy of the loop in
// stage1/pism/src/hydrology/Routing.cc (Routing::update_impl, lines ~841-976),
// with dump hooks added; every model method called is the protected base
// method itself, so the numbers are the numbers the base class would compute.
// The driver validates the replication: its final state must be bit-for-bit
// identical to the stage-1 dump (columns.csv flux/W_out/W_till_out,
// staggered.csv Qstag_avg).
//
// Dumps (written next to the stage-1 dumps, stage1/dumps/hydrology/):
//
//   substeps.csv          - per internal substep: k (1-based), hdt (s),
//                           cumulative time ht at the START of the substep
//   W_per_substep.csv     - per substep: k, i, j, W, Wtill at the START of the
//                           substep (k = 1 is the initial state; the largest k
//                           is the input of the final update)
//   staggered_initial.csv - staggered subassembly fields at the start of
//                           substep 1 (the initial state): (i,j,o) with
//                           Wstag, Kstag, Vstag, Qstag -- the analytic-state
//                           pinning target for the .esm operator templates
//   staggered_last.csv    - the same fields at the start of the LAST substep
//                           (this is what stage1/dumps/hydrology/staggered.csv
//                           samples; the files must agree)
//
// Build (from the eqice repo root, after sourcing stage1/runs/env.sh):
//   mpicxx -std=c++17 -O2 $(pkg-config --cflags pism) $(pkg-config --cflags petsc) \
//     stage2/tools/instrument_hydrology_substeps.cc -o /tmp/instrument_hydrology_substeps \
//     $(pkg-config --libs pism) $(pkg-config --libs petsc) $(pkg-config --libs gsl) \
//     $(pkg-config --libs netcdf) $(pkg-config --libs fftw3) -L<build-local>/lib -ludunits2
// Run (from the repo root):
//   source stage1/runs/env.sh
//   /tmp/instrument_hydrology_substeps -Mx 31 -My 31 -dt_years 0.25 \
//     -dumps_dir stage1/dumps/hydrology

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
  "Stage-2 analysis driver: re-runs the stage-1 hydrology (Routing) one-step test-F\n"
  "configuration with per-substep dump hooks (see the file header).\n\n";

namespace pism {

namespace {

class SubstepDumpingRouting : public hydrology::Routing {
public:
  SubstepDumpingRouting(std::shared_ptr<const Grid> g) : hydrology::Routing(g) {}

  // dump gating: write the staggered fields and the W/Wtill state when
  // step_counter == dump_substep (0 = never)
  int dump_substep = 0;
  std::ofstream *stag_log = nullptr;  // k, i, j, o, Wstag, Kstag, Vstag, Qstag
  std::ofstream *W_log = nullptr;     // k, i, j, W, Wtill
  std::ofstream *substep_log = nullptr; // k, hdt, ht

  void dump_staggered(std::ofstream &out, int k) const {
    array::AccessScope list{&m_Wstag, &m_Kstag, &m_Vstag, &m_Qstag};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      for (int o = 0; o < 2; ++o) {
        out << k << "," << i << "," << j << "," << o << ","
            << m_Wstag(i, j, o) << "," << m_Kstag(i, j, o) << ","
            << m_Vstag(i, j, o) << "," << m_Qstag(i, j, o) << "\n";
      }
    }
  }

  void dump_W(int k) const {
    array::AccessScope list{&m_W, &m_Wtill};
    for (auto p : m_grid->points()) {
      const int i = p.i(), j = p.j();
      (*W_log) << k << "," << i << "," << j << "," << m_W(i, j) << "," << m_Wtill(i, j) << "\n";
    }
  }

  // Replication of Routing::update_impl() with per-substep dump hooks.
  void update_impl(double t, double dt, const hydrology::Inputs &inputs) override {
    ice_bottom_surface(*inputs.geometry, m_bottom_surface);

    double ht = t, hdt = 0.0;
    const double t_final = t + dt;
    const double dt_max = m_config->get_number("hydrology.maximum_time_step", "seconds");
    const double tillwat_max = m_config->get_number("hydrology.tillwat_max");

    m_Qstag_average.set(0.0);

    // make sure W has valid ghosts before starting hydrology steps
    m_W.update_ghosts();

    unsigned int step_counter = 0;
    for (; ht < t_final; ht += hdt) {
      step_counter++;

      // updates ghosts of m_Wstag
      water_thickness_staggered(m_W, inputs.geometry->cell_type, m_Wstag);

      double maxKW = 0.0;
      compute_conductivity(m_Wstag, subglacial_water_pressure(), m_bottom_surface,
                           m_Kstag, maxKW);

      compute_velocity(m_Wstag, subglacial_water_pressure(), m_bottom_surface,
                       m_Kstag, inputs.no_model_mask, m_Vstag);

      advective_fluxes(m_Vstag, m_W, m_Qstag);

      m_Qstag_average.add(hdt, m_Qstag);

      {
        const double dt_cfl = max_timestep_W_cfl();
        const double dt_diff_w = max_timestep_W_diff(maxKW);
        hdt = std::min(t_final - ht, dt_max);
        hdt = std::min(hdt, dt_cfl);
        hdt = std::min(hdt, dt_diff_w);
      }

      if (substep_log) {
        (*substep_log) << step_counter << "," << hdt << "," << ht << "\n";
      }
      if (W_log) {
        dump_W((int)step_counter);
      }
      if ((int)step_counter == dump_substep and stag_log) {
        dump_staggered(*stag_log, (int)step_counter);
      }

      {
        update_Wtill(hdt, m_Wtill, m_surface_input_rate, m_basal_melt_rate, m_Wtillnew);
        enforce_bounds(inputs.geometry->cell_type, inputs.no_model_mask,
                       0.0, tillwat_max, m_Wtillnew,
                       m_grounded_margin_change, m_grounding_line_change,
                       m_conservation_error_change, m_no_model_mask_change);
      }

      {
        update_W(hdt, m_surface_input_rate, m_basal_melt_rate,
                 m_W, m_Wstag, m_Wtill, m_Wtillnew, m_Kstag, m_Qstag, m_Wnew);
        enforce_bounds(inputs.geometry->cell_type, inputs.no_model_mask,
                       0.0, 0.0, m_Wnew,
                       m_grounded_margin_change, m_grounding_line_change,
                       m_conservation_error_change, m_no_model_mask_change);

        m_W.copy_from(m_Wnew);
      }

      m_Wtill.copy_from(m_Wtillnew);
    } // end of the time-stepping loop

    staggered_to_regular(inputs.geometry->cell_type, m_Qstag_average,
                         m_config->get_flag("hydrology.routing.include_floating_ice"), m_Q);
    m_Q.scale(1.0 / dt);
  }
};

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

    auto grid = Grid::Shallow(ctx, 900.0e3, 900.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::NOT_PERIODIC);

    // Test-F exact geometry (same as the stage-1 driver)
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

    // Forcing (same as the stage-1 driver)
    array::Scalar1 no_model_mask(grid, "no_model_mask");
    no_model_mask.set(0.0);

    array::Scalar1 surface_input_rate(grid, "surface_input_rate");
    surface_input_rate.set(config->get_number("constants.fresh_water.density") *
                           units::convert(ctx->unit_system(), 0.2, "m / year", "m / second"));

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

    array::Scalar1 ice_sliding_speed(grid, "ice_sliding_speed");
    ice_sliding_speed.set(units::convert(ctx->unit_system(), 100.0, "m / year", "m / second"));

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
    // Dump targets
    // ------------------------------------------------------------------
    std::ofstream substep_log = open_output(out_dir + "/substeps.csv");
    header(substep_log, "description",
           "Per internal substep of the Routing update: k (1-based), hdt (s), "
           "and the cumulative time ht at the START of the substep");
    header(substep_log, "columns", "k,hdt_s,ht_s");

    std::ofstream W_log = open_output(out_dir + "/W_per_substep.csv");
    header(W_log, "description",
           "W and Wtill at the start of each internal substep (k = 1 is the "
           "initial state; the largest k is the input of the final update)");
    header(W_log, "columns", "k,i,j,W,Wtill");

    std::ofstream stag_initial = open_output(out_dir + "/staggered_initial.csv");
    header(stag_initial, "description",
           "Staggered subassembly fields at the start of substep 1 (the initial "
           "analytic state): Wstag, Kstag, Vstag, Qstag");
    header(stag_initial, "columns", "k,i,j,o,Wstag,Kstag,Vstag,Qstag");

    std::ofstream stag_last = open_output(out_dir + "/staggered_last.csv");
    header(stag_last, "description",
           "Staggered subassembly fields at the start of the LAST internal substep "
           "(the state sampled by stage1/dumps/hydrology/staggered.csv)");
    header(stag_last, "columns", "k,i,j,o,Wstag,Kstag,Vstag,Qstag");

    Inputs hinputs;
    hinputs.no_model_mask = &no_model_mask;
    hinputs.geometry = &geometry;
    hinputs.surface_input_rate = &surface_input_rate;
    hinputs.basal_melt_rate = &basal_melt_rate;
    hinputs.ice_sliding_speed = &ice_sliding_speed;

    // Run 1: substep log + W per substep + initial-state staggered fields
    {
      auto model = std::make_shared<SubstepDumpingRouting>(grid);
      array::Scalar1 P0(grid, "overburden_pressure");
      P0.set(0.0);
      model->init(W_till, W, P0);

      model->substep_log = &substep_log;
      model->W_log = &W_log;
      model->stag_log = &stag_initial;
      model->dump_substep = 1;

      model->update(0.0, dt, hinputs);
    }
    substep_log.close();
    W_log.close();
    stag_initial.close();

    // count the substeps
    int n_substeps = 0;
    {
      std::ifstream in(out_dir + "/substeps.csv");
      std::string line;
      while (std::getline(in, line)) {
        if (not line.empty() and line[0] != '#') {
          n_substeps++;
        }
      }
    }
    std::cout << "took " << n_substeps << " substeps\n";

    // Run 2: last-substep staggered fields + final-state dump
    std::ofstream final_log;
    {
      auto model = std::make_shared<SubstepDumpingRouting>(grid);
      array::Scalar1 P0(grid, "overburden_pressure");
      P0.set(0.0);
      model->init(W_till, W, P0);

      model->substep_log = nullptr;
      model->W_log = nullptr;
      model->stag_log = &stag_last;
      model->dump_substep = n_substeps;

      model->update(0.0, dt, hinputs);

      final_log = open_output(out_dir + "/final_state.csv");
      header(final_log, "description",
             "Final state after the full dt step and the regular-grid time-averaged "
             "flux (the replication's output; must match stage1 columns.csv)");
      header(final_log, "columns", "i,j,W,Wtill,overburden_pressure,flux_u,flux_v");
      {
        const array::Scalar &Wf = model->subglacial_water_thickness();
        const array::Scalar &Wtf = model->till_water_thickness();
        const array::Scalar &Pf = model->overburden_pressure();
        const array::Vector &Qf = model->flux();
        array::AccessScope list{&Wf, &Wtf, &Pf, &Qf};
        for (auto p : grid->points()) {
          const int i = p.i(), j = p.j();
          final_log << i << "," << j << "," << Wf(i, j) << "," << Wtf(i, j) << ","
                    << Pf(i, j) << "," << Qf(i, j).u << "," << Qf(i, j).v << "\n";
        }
      }
      final_log.close();
    }
    stag_last.close();

    std::cout << "done; per-substep dumps in " << out_dir << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
