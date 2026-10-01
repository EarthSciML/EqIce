// Stage 1 instrumentation driver for the PISM ice-rheology component boundary.
//
// Dumps instantaneous inputs/outputs of the ice flow law subassembly:
//
//   - FlowLaw::softness(E, p)                    -> A  (Pa^-n s^-1)
//   - FlowLaw::hardness(E, p)                    -> B  (Pa s^(1/n))
//   - FlowLaw::flow(stress, E, p, gs)            -> strain rate (s^-1)
//   - FlowLaw::flow_n (vectorized batch path)
//   - FlowLaw::effective_viscosity(B, gamma, eps)-> (nu, dnu)
//   - secondInvariant_2D(U_x, U_y)               -> gamma (s^-2)
//   - averaged_hardness (column integral)        -> \bar B
//   - EnthalpyConverter methods used to convert flow law inputs
//
// for every flow law registered in PISM's FlowLawFactory.
//
// The dump files are the authoritative reference traces that Stage-2 .esm
// component tests are built from (see PLAN.md section 4 and catalog.md).
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/rheology.json

#include <petsc.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "pism/util/Config.hh"
#include "pism/util/EnthalpyConverter.hh"
#include "pism/util/Logger.hh"
#include "pism/util/Vector2d.hh"
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/Units.hh"
#include "pism/rheology/FlowLaw.hh"
#include "pism/rheology/FlowLawFactory.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM ice-rheology component.\n"
  "Dumps flow law I/O traces to stage1/dumps/rheology/.\n\n";

namespace pism {
extern const char *revision; // defined in pism_config.cc (compiled in)
} // namespace pism

namespace {

using pism::rheology::FlowLaw;
using pism::rheology::FlowLawFactory;
using pism::EnthalpyConverter;
using pism::Config;

// ---------------------------------------------------------------------------
// Sample grids. These are the deterministic input distributions used for the
// dumps; they are recorded in stage1/runs/rheology.json.
//
// Enthalpy range: valid (E < E_liquid(p)) for every pressure in m_pressures.
//   E_liquid(p=4.5e7 Pa) = E_cts + L ~= 416000 J kg^-1 with default config.
// ---------------------------------------------------------------------------

const std::vector<double> m_enthalpies = {
    0.0,     2.0e4,  4.0e4,  6.0e4,  8.0e4,  9.5e4, 1.1e5,
    1.35e5,  1.7e5,  2.2e5,  2.8e5,  3.5e5,  4.05e5}; // J kg^-1

const std::vector<double> m_pressures = {
    0.0, 1.0e6, 5.0e6, 1.0e7, 2.0e7, 3.0e7, 4.5e7}; // Pa

const std::vector<double> m_stresses = {
    1.0e3, 3.1622776601683795e3, 1.0e4, 3.1622776601683795e4,
    1.0e5, 3.1622776601683795e5, 1.0e6}; // Pa

const std::vector<double> m_grain_sizes = {1.0e-4, 1.0e-3, 3.0e-3, 1.0e-2}; // m

const std::vector<double> m_gammas = {
    1.0e-14, 3.1622776601683795e-14, 1.0e-13, 3.1622776601683795e-13,
    1.0e-12, 3.1622776601683795e-12, 1.0e-11, 3.1622776601683795e-11,
    1.0e-10, 3.1622776601683795e-10, 1.0e-9, 3.1622776601683795e-9,
    1.0e-8, 3.1622776601683795e-8, 1.0e-7, 3.1622776601683795e-7,
    1.0e-6}; // s^-2 (second invariant gamma = 0.5 D_ij D_ij)

const std::vector<double> m_hardnesses = {
    1.0e4, 3.1622776601683795e4, 1.0e5, 3.1622776601683795e5,
    1.0e6, 3.1622776601683795e6, 1.0e7, 3.1622776601683795e7,
    1.0e8, 3.1622776601683795e8, 1.0e9}; // Pa s^(1/n)

const std::vector<std::string> m_flow_law_names = {
    "isothermal_glen", "pb", "arr", "arrwarm", "gpbld", "hooke", "gk"};

// ---------------------------------------------------------------------------
// Output helpers
// ---------------------------------------------------------------------------

std::ofstream open_output(const std::string &path) {
  std::ofstream out(path);
  if (not out) {
    throw std::runtime_error("Cannot open output file '" + path + "'");
  }
  out << std::setprecision(17);
  return out;
}

//! Write a `# key = value` header line.
void header(std::ofstream &out, const std::string &key, const std::string &value) {
  out << "# " << key << " = " << value << "\n";
}

//! Number of vertical levels below height H in a z-level array, replicating
//! Grid::kBelowHeight() semantics (GSL accel find: largest k with z[k] <= H,
//! i.e. the interval index).
unsigned int k_below_height(const std::vector<double> &z, double H) {
  unsigned int k = 0;
  for (size_t i = 0; i < z.size(); ++i) {
    if (z[i] < H) {
      k = static_cast<unsigned int>(i);
    }
  }
  return k;
}

// ---------------------------------------------------------------------------
// Dump sections
// ---------------------------------------------------------------------------

//! Dump the EnthalpyConverter boundary (law-independent).
void dump_enthalpy_converter(const Config &config, const EnthalpyConverter &EC,
                             const std::string &dir) {

  std::ofstream out = open_output(dir + "/enthalpy_converter.csv");

  header(out, "description",
         "EnthalpyConverter boundary: E,p -> T, T_pa, omega, is_temperate, E_cts, E_l, T_m");
  header(out, "columns", "E,p,T,T_pa,omega,is_temperate,E_cts,E_l,T_m");
  for (double E : m_enthalpies) {
    for (double p : m_pressures) {
      out << E << "," << p << ","
          << EC.temperature(E, p) << ","
          << EC.pressure_adjusted_temperature(E, p) << ","
          << EC.water_fraction(E, p) << ","
          << (EC.is_temperate(E, p) ? 1 : 0) << ","
          << EC.enthalpy_cts(p) << ","
          << EC.enthalpy_liquid(p) << ","
          << EC.melting_temperature(p) << "\n";
    }
  }
  out.close();

  // pressure(depth) boundary
  std::ofstream pout = open_output(dir + "/enthalpy_converter_pressure.csv");
  header(pout, "description", "EnthalpyConverter::pressure(depth): p_air + rho_i g depth");
  header(pout, "columns", "depth,pressure");
  const double depths[] = {-10.0, 0.0, 1.0, 500.0, 1000.0, 2000.0, 4000.0};
  for (double d : depths) {
    pout << d << "," << EC.pressure(d) << "\n";
  }
  pout.close();

  // scalar constants
  std::ofstream cout = open_output(dir + "/enthalpy_converter_constants.csv");
  header(cout, "description", "EnthalpyConverter scalar constants and derived values");
  header(cout, "columns", "name,value");
  cout << "c_i," << EC.c() << "\n";
  cout << "L(T_m(0))," << EC.L(EC.melting_temperature(0.0)) << "\n";
  cout << "L(T_m(4.5e7))," << EC.L(EC.melting_temperature(4.5e7)) << "\n";
  cout.close();
}

//! Dump secondInvariant_2D (law-independent; free function in FlowLaw.hh).
void dump_second_invariant(const std::string &dir) {

  std::ofstream out = open_output(dir + "/second_invariant.csv");
  header(out, "description",
         "secondInvariant_2D: velocity gradient vectors -> gamma = 0.5 D_ij D_ij (s^-2)");
  header(out, "columns", "u_x,u_y,v_x,v_y,gamma");

  const double rates[] = {1.0e-13, 1.0e-12, 1.0e-11, 1.0e-10, 1.0e-9};
  // canonical gradient tensors: pure shear, pure normal (incompressible),
  // simple shear, mixed
  const struct {
    double u_x, u_y, v_x, v_y;
  } tensors[] = {
      {0.0, 1.0, 1.0, 0.0},    // pure shear
      {1.0, 0.0, 0.0, -1.0},   // pure normal (w_z = 0)
      {1.0, 0.0, 0.0, 0.0},    // simple shear x
      {0.0, 0.0, 1.0, 0.0},    // simple shear y
      {1.0, 0.5, -0.5, -1.0},  // mixed
      {1.0, -1.0, 1.0, 1.0},   // mixed
  };

  for (double r : rates) {
    for (const auto &t : tensors) {
      pism::Vector2d U_x{t.u_x * r, t.u_y * r};
      pism::Vector2d U_y{t.v_x * r, t.v_y * r};
      out << U_x.u << "," << U_x.v << "," << U_y.u << "," << U_y.v << ","
          << pism::secondInvariant_2D(U_x, U_y) << "\n";
    }
  }
  out.close();
}

//! Dump effective_viscosity (law-independent; implemented in the FlowLaw base class).
void dump_effective_viscosity(const FlowLaw &flow_law, double schoof_reg,
                              const std::string &dir, const std::string &law_name) {

  std::ofstream out = open_output(dir + "/effective_viscosity_" + law_name + ".csv");
  header(out, "description",
         "FlowLaw::effective_viscosity(B, gamma, eps) -> nu, dnu; "
         "eps = 0 and eps = Schoof regularization");
  header(out, "columns", "B,gamma,eps,nu,dnu");
  for (double B : m_hardnesses) {
    for (double gamma : m_gammas) {
      for (double eps : {0.0, schoof_reg}) {
        double nu = 0.0, dnu = 0.0;
        flow_law.effective_viscosity(B, gamma, eps, &nu, &dnu);
        out << B << "," << gamma << "," << eps << "," << nu << "," << dnu << "\n";
      }
    }
  }
  out.close();
}

//! Dump the per-flow-law boundary: softness, hardness, flow (scalar path).
void dump_flow_law(const FlowLaw &flow_law, const EnthalpyConverter &EC,
                   const std::string &dir, const std::string &law_name) {

  std::ofstream out = open_output(dir + "/flowlaw_" + law_name + ".csv");
  header(out, "description",
         "FlowLaw boundary: stress,E,p,gs -> softness, hardness, flow; "
         "plus EC conversions used to map E,p -> T");
  header(out, "columns",
         "stress,E,p,gs,softness,hardness,flow,T,T_pa,omega,is_temperate");

  for (double stress : m_stresses) {
    for (double E : m_enthalpies) {
      for (double p : m_pressures) {
        for (double gs : m_grain_sizes) {
          double softness = 0.0, hardness = 0.0, flow = 0.0;
          try {
            softness = flow_law.softness(E, p);
          } catch (const std::exception &) {
            softness = std::numeric_limits<double>::quiet_NaN();
          }
          try {
            hardness = flow_law.hardness(E, p);
          } catch (const std::exception &) {
            hardness = std::numeric_limits<double>::quiet_NaN();
          }
          try {
            flow = flow_law.flow(stress, E, p, gs);
          } catch (const std::exception &) {
            flow = std::numeric_limits<double>::quiet_NaN();
          }

          out << stress << "," << E << "," << p << "," << gs << ","
              << softness << "," << hardness << "," << flow << ","
              << EC.temperature(E, p) << ","
              << EC.pressure_adjusted_temperature(E, p) << ","
              << EC.water_fraction(E, p) << ","
              << (EC.is_temperate(E, p) ? 1 : 0) << "\n";
        }
      }
    }
  }
  out.close();

  // vectorized (flow_n) path, used by the SIA stress balance
  std::ofstream bout = open_output(dir + "/flowlaw_" + law_name + "_batched.csv");
  header(bout, "description",
         "FlowLaw::flow_n: batched (stress,E,p,gs) -> flow (vectorized path)");
  header(bout, "columns", "index,stress,E,p,gs,flow_n");

  const size_t batch = 8;
  std::vector<double> stress(batch), E(batch), p(batch), gs(batch), result(batch);
  for (size_t k = 0; k < batch; ++k) {
    stress[k] = m_stresses[(3 * k) % m_stresses.size()];
    E[k]      = m_enthalpies[(5 * k) % m_enthalpies.size()];
    p[k]      = m_pressures[(2 * k) % m_pressures.size()];
    gs[k]     = m_grain_sizes[k % m_grain_sizes.size()];
  }
  flow_law.flow_n(stress.data(), E.data(), p.data(), gs.data(), batch, result.data());
  for (size_t k = 0; k < batch; ++k) {
    bout << k << "," << stress[k] << "," << E[k] << "," << p[k] << "," << gs[k] << ","
         << result[k] << "\n";
  }
  bout.close();

  // column-averaged hardness (subassembly: trapezoidal + rectangle integration)
  std::ofstream aout = open_output(dir + "/averaged_hardness_" + law_name + ".csv");
  header(aout, "description",
         "averaged_hardness: column of (zlevels, enthalpy) -> \\bar B (Pa s^(1/n)); "
         "trapezoidal rule below kbelowH, rectangle rule above");
  header(aout, "columns", "H,kbelowH,B_averaged");

  const std::vector<double> zlevels = {0.0,   200.0, 400.0,  600.0,  800.0, 1000.0,
                                       1200.0, 1400.0, 1600.0, 1800.0, 2000.0}; // m
  const std::vector<double> E_profile = {1.5e5, 1.42e5, 1.34e5, 1.26e5, 1.18e5,
                                         1.1e5, 1.0e5,  9.0e4,  8.0e4,  7.0e4,
                                         6.0e4}; // J kg^-1, warm at bed -> cold at surface

  for (double H : {2000.0, 850.0}) {
    unsigned int kbelowH = k_below_height(zlevels, H);
    double B = pism::rheology::averaged_hardness(flow_law, H, kbelowH,
                                                 zlevels.data(), E_profile.data());
    aout << H << "," << kbelowH << "," << B << "\n";
  }
  aout.close();
}

//! Dump all config parameters that affect the rheology component.
void dump_parameters(const Config &config, const std::string &dir) {
  std::ofstream out = open_output(dir + "/parameters.csv");
  header(out, "description",
         "Configuration parameters affecting the ice-rheology component "
         "(flow_law.*, constants.*, enthalpy_converter.*, energy.model, "
         "surface.pressure, stress_balance.*.Glen_exponent)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {
      "flow_law.", "constants.", "enthalpy_converter.", "energy.model",
      "surface.pressure", "stress_balance.sia.Glen_exponent",
      "stress_balance.ssa.Glen_exponent"};

  auto keys = config.keys();
  for (const auto &key : keys) {
    bool match = false;
    for (const auto &prefix : prefixes) {
      if (key.rfind(prefix, 0) == 0) {
        match = true;
        break;
      }
    }
    if (match) {
      // Skip metadata keys (key_doc, key_type, key_units, ...)
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
  }
  out.close();
}

} // namespace

int main(int argc, char *argv[]) {

  using namespace pism;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;

  try {
    units::System::Ptr sys(new units::System);
    auto logger = std::make_shared<Logger>(com, 6);
    auto config = config_from_options(com, sys);
    set_config_from_options(*config);
    config->resolve_filenames();

    auto EC = std::make_shared<EnthalpyConverter>(*config);

    std::string out_dir = "stage1/dumps/rheology";
    options::String output_dir("-dumps_dir", "Output directory for the dumps",
                               out_dir);
    out_dir = output_dir;

    // make sure the output directory exists
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    std::ofstream meta = open_output(out_dir + "/meta.txt");
    header(meta, "pism_revision", pism::revision);
    header(meta, "pism_config_file", pism::config_file);
    header(meta, "mpi_rank0", "1");
    meta << "# command:";
    for (int i = 0; i < argc; ++i) {
      meta << " " << argv[i];
    }
    meta << "\n";
    meta.close();

    // Schoof regularization parameter, exactly as computed in FlowLaw's ctor:
    double schoof_len = config->get_number("flow_law.Schoof_regularizing_length", "m"),
           schoof_vel = config->get_number("flow_law.Schoof_regularizing_velocity",
                                           "m second-1");
    double schoof_reg = PetscSqr(schoof_vel / schoof_len);

    // parameters (before flow law construction, same values are used by all)
    dump_parameters(*config, out_dir);

    // law-independent boundaries
    dump_enthalpy_converter(*config, *EC, out_dir);
    dump_second_invariant(out_dir);

    // per-law boundaries
    FlowLawFactory factory(config, EC);
    double n = config->get_number("stress_balance.sia.Glen_exponent");
    for (const auto &law_name : m_flow_law_names) {
      auto flow_law = factory.create(law_name, n);
      dump_effective_viscosity(*flow_law, schoof_reg, out_dir, law_name);
      dump_flow_law(*flow_law, *EC, out_dir, law_name);
      std::cout << "Dumped: " << law_name << " (" << flow_law->name() << ")\n";
    }

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
