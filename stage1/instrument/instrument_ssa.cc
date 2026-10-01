// Stage 1 instrumentation driver for the PISM SSA stress-balance boundary.
//
// Uses the verification test V "circular flat-bottomed continent" (CFBC)
// configuration -- the van der Veen flow-line floating shelf -- exactly as
// PISM's own pism_ssa_test_cfbc does, with the calving-front stress boundary
// condition. The exact solution is known: H(x), u(x) = V0*H0/H(x), v = 0.
//
// The flow law is isothermal_glen with hardness B = 1.9e8 Pa s^(1/n) (the test
// sets flow_law.isothermal_Glen.ice_softness accordingly), which corresponds
// to the van der Veen rate constant C = 2.45e-18.
//
// Dumps (instantaneous, one SSA solve):
//
//   inputs.nc       - all SSA inputs + outputs (geometry, tauc, bc_mask,
//                     bc_values, enthalpy, velocity, F_b, driving stress,
//                     nuH, ice hardness)
//   parameters.csv  - config keys affecting the SSA boundary
//   columns.csv     - per cell: geometry, tauc, computed (u, v, F_b, taud_x/y,
//                     nuH, Bbar) and exact reference (H, u, v)
//   meta.txt        - PISM revision, config file, grid, command line
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/ssa.json

#include <petsc.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <set>
#include <string>

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
#include "pism/stressbalance/StressBalance.hh"
#include "pism/stressbalance/ssa/SSA.hh"
#include "pism/stressbalance/ssa/SSAFD.hh"
#include "pism/stressbalance/ssa/SSAFDBase.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM SSA stress-balance component.\n"
  "Dumps SSA I/O traces (test V / van der Veen CFBC shelf) to "
  "stage1/dumps/ssa/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

using pism::stressbalance::SSAFD;

//! Subclass exposing SSAFDBase internals needed for the dumps.
class InstrumentedSSAFD : public SSAFD {
public:
  InstrumentedSSAFD(std::shared_ptr<const Grid> g) : SSAFD(g, false) {}

  const array::Staggered1 &nuH() const { return m_nuH; }
  const array::Staggered &hardness() const { return m_hardness; }
  const array::Vector &taud() const { return m_taud; }
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

void dump_parameters(const Config &config, const std::string &dir) {
  std::ofstream out = open_output(dir + "/parameters.csv");
  header(out, "description",
         "Configuration parameters affecting the SSA stress-balance component "
         "(stress_balance.ssa.*, flow_law.*, constants.*, grid.*, "
         "basal_resistance.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {
      "stress_balance.ssa.", "stress_balance.calving_front_stress_bc",
      "flow_law.", "constants.", "grid.", "basal_resistance.",
      "enthalpy_converter."};

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

// van der Veen exact solution (from pism_ssa_test_cfbc)
double H_exact(double V0, double H0, double C, double x) {
  const double Q0 = V0 * H0;
  return pow(4.0 * C / Q0 * x + 1.0 / pow(H0, 4), -0.25);
}

double u_exact(double V0, double H0, double C, double x) {
  const double Q0 = V0 * H0;
  return Q0 / H_exact(V0, H0, C, x);
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
    std::shared_ptr<Context> ctx = context_from_options(com, "ssa_test");
    auto config = ctx->config();

    // Same configuration as pism_ssa_test_cfbc (test V / van der Veen):
    config->set_number("flow_law.isothermal_Glen.ice_softness",
                       pow(1.9e8, -config->get_number("stress_balance.ssa.Glen_exponent")));
    config->set_flag("stress_balance.ssa.compute_surface_gradient_inward", false);
    config->set_flag("stress_balance.calving_front_stress_bc", true);
    config->set_flag("stress_balance.ssa.fd.flow_line_mode", true);
    config->set_flag("stress_balance.ssa.fd.extrapolate_at_margins", false);
    config->set_string("stress_balance.ssa.flow_law", "isothermal_glen");

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/ssa";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    unsigned int Mx = config->get_number("grid.Mx");
    unsigned int My = config->get_number("grid.My");

    // test V domain: 250 km half-width, cell-centered, periodic in Y
    auto grid = Grid::Shallow(ctx, 250.0e3, 250.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::Y_PERIODIC);

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
      header(meta, "dx", std::to_string(grid->dx()));
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

    array::Scalar1 tauc(grid, "tauc");
    array::Array3D enthalpy(grid, "enthalpy", array::WITH_GHOSTS, grid->z(), 1);
    array::Vector2 bc_values(grid, "_bc");
    array::Scalar2 bc_mask(grid, "bc_mask");
    bc_mask.set_interpolation_type(NEAREST);

    Geometry geometry(grid);

    // initialize coefficients (as in SSATestCaseCFBC)
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

        if (i != (int)grid->Mx() - 1) {
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

    // SSA solve (FD Picard)
    auto ssa = std::make_shared<InstrumentedSSAFD>(grid);
    ssa->init();

    Inputs inputs;
    inputs.water_column_pressure = nullptr;
    inputs.geometry = &geometry;
    inputs.enthalpy = &enthalpy;
    inputs.basal_yield_stress = &tauc;
    inputs.bc_mask = &bc_mask;
    inputs.bc_values = &bc_values;

    ssa->update(inputs, true);

    // ------------------------------------------------------------------
    // Dumps
    // ------------------------------------------------------------------
    const array::Vector1 &velocity = ssa->velocity();
    const array::Scalar &F_b = ssa->basal_frictional_heating();
    const array::Vector &taud = ssa->taud();
    const array::Staggered1 &nuH = ssa->nuH();
    const array::Staggered &hardness = ssa->hardness();

    {
      // inputs + outputs as NetCDF
      auto writer = std::make_shared<SynchronousOutputWriter>(grid->com, *config);
      writer->initialize({}, true);
      OutputFile file(writer, out_dir + "/inputs.nc");

      const array::Array *vecs[] = {
          &geometry.ice_surface_elevation, &geometry.ice_thickness,
          &geometry.bed_elevation,         &geometry.cell_type,
          &tauc,                           &bc_mask,
          &bc_values,                      &enthalpy,
          &velocity,                       &F_b,
          &taud,                           &nuH,
          &hardness};

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

    {
      // per-cell columns: inputs + outputs + exact reference
      std::ofstream out = open_output(out_dir + "/columns.csv");
      header(out, "description",
             "SSA output columns (test V / van der Veen CFBC shelf): per cell: "
             "geometry, tauc, computed u/v/F_b/taud/nuH/Bbar, exact H/u/v");
      header(out, "columns",
             "i,j,x,y,H,surface,bed,mask,tauc,u,v,F_b,taud_x,taud_y,nuH,Bbar,"
             "H_exact,u_exact,v_exact");

      array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                              &geometry.bed_elevation, &geometry.cell_type, &tauc,
                              &velocity, &F_b, &taud, &nuH, &hardness};

      const double x_min = grid->x(0);

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double x = grid->x(i), y = grid->y(j);

        const bool have_exact = (i != (int)grid->Mx() - 1);
        const double Hx = have_exact ? H_exact(V0, H0, C, x - x_min) : 0.0;
        const double ux = have_exact ? u_exact(V0, H0, C, x - x_min) : 0.0;

        out << i << "," << j << "," << x << "," << y << ","
            << geometry.ice_thickness(i, j) << "," << geometry.ice_surface_elevation(i, j)
            << "," << geometry.bed_elevation(i, j) << "," << (int)geometry.cell_type(i, j)
            << "," << tauc(i, j) << "," << velocity(i, j).u << "," << velocity(i, j).v
            << "," << F_b(i, j) << "," << taud(i, j).u << "," << taud(i, j).v << ","
            << nuH(i, j, 0) << "," << hardness(i, j, 0) << "," << Hx << "," << ux << ",0"
            << "\n";
      }
      out.close();
    }

    // velocity error report (same computation as SSATestCase::report)
    {
      double maxvecerr = 0.0, avvecerr = 0.0, avuerr = 0.0, avverr = 0.0;
      double maxuerr = 0.0, maxverr = 0.0, exactvelmax = 0.0;

      array::AccessScope list{&velocity};
      const double x_min = grid->x(0);

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if (i == (int)grid->Mx() - 1) {
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

      const unsigned int N = grid->Mx() * grid->My();
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

    std::cout << "SSA instrumentation complete.\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
