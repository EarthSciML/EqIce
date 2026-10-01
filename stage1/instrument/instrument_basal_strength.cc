// Stage 1 instrumentation driver for the PISM basal-strength boundary.
//
// Dumps instantaneous inputs/outputs of the basal yield stress and basal
// resistance subassemblies:
//
//   (a) MohrCoulombPointwise (pure pointwise functions, ideal .esm test targets):
//         - effective_pressure(delta, P_overburden, water_thickness) -> N_till
//         - yield_stress(delta, P_overburden, water_thickness, phi)  -> tauc
//         - till_friction_angle(delta, P_overburden, water_thickness,
//                               yield_stress) -> phi  (inverse of yield_stress)
//       sampled over deterministic grids recorded in stage1/runs/basal_strength.json.
//   (b) the basal resistance laws (basal_resistance.{hh,cc}):
//         - IceBasalResistancePlasticLaw::drag / drag_with_derivative
//         - IceBasalResistancePseudoPlasticLaw::drag / drag_with_derivative
//         - IceBasalResistanceRegularizedLaw::drag / drag_with_derivative
//       over a grid of (tauc, vx, vy) spanning zero / slow / fast sliding.
//   (c) a tauc map computed by the actual MohrCoulombYieldStress component on the
//       verification test-F exact geometry (flat bed, thickness = exact H),
//       with a prescribed smooth till-water-thickness field and a prescribed
//       till friction angle field; the effective pressure field is recomputed
//       pointwise from the same inputs (update_impl does not store it).
//
// The dump files are the authoritative reference traces that Stage-2 .esm
// component tests are built from (see PLAN.md section 4 and boundaries.md
// section 5).
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/basal_strength.json

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
#include "pism/util/error_handling.hh"
#include "pism/util/petscwrappers/PetscInitializer.hh"
#include "pism/util/pism_options.hh"
#include "pism/util/pism_utilities.hh"
#include "pism/util/array/Scalar.hh"
#include "pism/basalstrength/MohrCoulombPointwise.hh"
#include "pism/basalstrength/MohrCoulombYieldStress.hh"
#include "pism/basalstrength/basal_resistance.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM basal-strength component.\n"
  "Dumps Mohr-Coulomb yield stress, till friction angle, effective pressure,\n"
  "and basal resistance law traces to stage1/dumps/basal_strength/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

// ---------------------------------------------------------------------------
// Sample grids. These are the deterministic input distributions used for the
// dumps; they are recorded in stage1/runs/basal_strength.json.
//
// Physical regime (with default config):
//   hydrology.tillwat_max                        = 2.0 m
//   basal_yield_stress.mohr_coulomb.delta.file   = "" -> delta = config value
//   P_overburden up to ~3.6e7 Pa (4000 m of ice).
// Water thickness is sampled across [0, 2*Wmax] plus one negative value (the
// pointwise formulas accept it; hydrology normally keeps it non-negative).
// ---------------------------------------------------------------------------

const std::vector<double> m_deltas = {0.01, 0.02, 0.05, 0.2};   // config default 0.02
const std::vector<double> m_P_overburden = {
    0.0, 1.0e6, 5.0e6, 1.0e7, 2.0e7, 3.0e7, 5.0e7};            // Pa
const std::vector<double> m_water_thickness = {
    -1.0, 0.0, 0.1, 0.5, 1.0, 1.5, 2.0, 2.5, 4.0};             // m (Wmax = 2.0)
const std::vector<double> m_phis = {0.0, 5.0, 10.0, 20.0, 30.0, 45.0}; // degrees

// basal resistance law sample grids
const std::vector<double> m_taucs = {0.0, 1.0e3, 1.0e4, 5.0e4, 1.0e5, 2.0e5, 1.0e6}; // Pa
// speeds in m/s: 0, ~0.3 m/yr, ~3 m/yr, ~30 m/yr, 100 m/yr (pseudo threshold),
// ~300 m/yr, ~3000 m/yr
const std::vector<double> m_speeds = {
    0.0, 1.0e-8, 1.0e-7, 1.0e-6, 3.1689e-6, 1.0e-5, 1.0e-4};

// direction unit vectors (beta depends only on |v|; these show that)
const std::vector<std::pair<double, double>> m_directions = {
    {1.0, 0.0}, {0.0, 1.0}, {1.0 / std::sqrt(2.0), 1.0 / std::sqrt(2.0)},
    {1.0 / std::sqrt(2.0), -1.0 / std::sqrt(2.0)}};

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
         "Configuration parameters affecting the basal-strength component "
         "(basal_yield_stress.*, basal_resistance.*, hydrology.tillwat_max, "
         "constants.*, grid.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {
      "basal_yield_stress.", "basal_resistance.", "hydrology.tillwat_max",
      "constants.", "grid."};

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
    std::shared_ptr<Context> ctx = context_from_options(com, "basal_strength_test");
    auto config = ctx->config();

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/basal_strength";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    const unsigned int Mx = config->get_number("grid.Mx");
    const unsigned int My = config->get_number("grid.My");

    // config keys that enter the Mohr-Coulomb pointwise formulas
    const double W_till_max    = config->get_number("hydrology.tillwat_max");
    const double till_cohesion = config->get_number("basal_yield_stress.mohr_coulomb.till_cohesion");
    const double delta_config  = config->get_number("basal_yield_stress.mohr_coulomb.till_effective_fraction_overburden");
    const double N0            = config->get_number("basal_yield_stress.mohr_coulomb.till_reference_effective_pressure");
    const double e0            = config->get_number("basal_yield_stress.mohr_coulomb.till_reference_void_ratio");
    const double Cc            = config->get_number("basal_yield_stress.mohr_coulomb.till_compressibility_coefficient");

    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "component", "basal_strength (yield stress + basal resistance)");
      header(meta, "description",
             "MohrCoulombPointwise (N_till, tauc, phi) sample grids; basal "
             "resistance law drag/derivative grids; MohrCoulombYieldStress tauc "
             "map on the test-F exact geometry (Mx x My, flat bed)");
      header(meta, "Mx", std::to_string(Mx));
      header(meta, "My", std::to_string(My));
      header(meta, "Lx", "900000");
      header(meta, "dx", std::to_string(2.0 * 900.0e3 / Mx));
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    MohrCoulombPointwise mc(config);

    // ------------------------------------------------------------------
    // (a) Mohr-Coulomb pointwise functions
    // ------------------------------------------------------------------

    // effective pressure N_till (with the pre-saturation interior value, so a
    // .esm test can separate the formula from the min(P, .) clamp)
    {
      std::ofstream out = open_output(out_dir + "/effective_pressure.csv");
      header(out, "description",
             "MohrCoulombPointwise::effective_pressure: "
             "N_till = min(P_overburden, N0 (delta P / N0)^s 10^((e0/Cc)(1-s))) "
             "with s = water_thickness / W_till_max; N_interior is the "
             "un-clamped formula value");
      header(out, "columns",
             "delta,P_overburden,water_thickness,N_interior,N_till");

      for (double delta : m_deltas) {
        for (double P : m_P_overburden) {
          for (double W : m_water_thickness) {
            const double s = W / W_till_max;
            const double N_interior =
                N0 * std::pow(delta * P / N0, s) *
                std::pow(10.0, (e0 / Cc) * (1.0 - s));
            const double N_till = mc.effective_pressure(delta, P, W);
            out << delta << "," << P << "," << W << "," << N_interior << ","
                << N_till << "\n";
          }
        }
      }
      out.close();
    }

    // yield stress tauc = c0 + tan(phi) N_till
    {
      std::ofstream out = open_output(out_dir + "/yield_stress.csv");
      header(out, "description",
             "MohrCoulombPointwise::yield_stress: tauc = c0 + tan(phi_rad) N_till "
             "(phi in degrees in the input column)");
      header(out, "columns",
             "delta,P_overburden,water_thickness,phi_deg,N_till,tauc");

      for (double delta : m_deltas) {
        for (double P : m_P_overburden) {
          for (double W : m_water_thickness) {
            const double N_till = mc.effective_pressure(delta, P, W);
            for (double phi : m_phis) {
              const double tauc = mc.yield_stress(delta, P, W, phi);
              out << delta << "," << P << "," << W << "," << phi << "," << N_till
                  << "," << tauc << "\n";
            }
          }
        }
      }
      out.close();
    }

    // till_friction_angle: the inverse of yield_stress. Dump the round trip
    // phi -> tauc -> phi_inverse (must be exact).
    {
      std::ofstream out = open_output(out_dir + "/till_friction_angle.csv");
      header(out, "description",
             "MohrCoulombPointwise::till_friction_angle (inverse of "
             "yield_stress): for each (delta, P, W, phi) we dump tauc = "
             "yield_stress(delta, P, W, phi) and phi_inverse = "
             "till_friction_angle(delta, P, W, tauc); the round trip must be exact");
      header(out, "columns",
             "delta,P_overburden,water_thickness,phi_deg,tauc,phi_inverse_deg");

      for (double delta : m_deltas) {
        for (double P : m_P_overburden) {
          for (double W : m_water_thickness) {
            for (double phi : m_phis) {
              const double tauc = mc.yield_stress(delta, P, W, phi);
              const double phi_inv = mc.till_friction_angle(delta, P, W, tauc);
              out << delta << "," << P << "," << W << "," << phi << "," << tauc
                  << "," << phi_inv << "\n";
            }
          }
        }
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // (b) basal resistance laws
    // ------------------------------------------------------------------
    struct DragLawDesc {
      const char *name;
      std::unique_ptr<IceBasalResistancePlasticLaw> law;
    };
    std::vector<DragLawDesc> drag_laws;
    drag_laws.push_back({"plastic", std::make_unique<IceBasalResistancePlasticLaw>(*config)});
    drag_laws.push_back({"pseudo_plastic", std::make_unique<IceBasalResistancePseudoPlasticLaw>(*config)});
    drag_laws.push_back({"regularized", std::make_unique<IceBasalResistanceRegularizedLaw>(*config)});

    for (const auto &dl : drag_laws) {
      std::ofstream out = open_output(out_dir + "/drag_" + std::string(dl.name) + ".csv");
      header(out, "description",
             std::string("IceBasalResistance") + (std::string(dl.name) == "plastic"
                          ? "PlasticLaw" : std::string(dl.name) == "pseudo_plastic"
                          ? "PseudoPlasticLaw" : "RegularizedLaw") +
             ": beta = drag(tauc, vx, vy), dbeta = drag_with_derivative(...) "
             "w.r.t. alpha = 0.5 |v|^2; tau_b = -beta v. dbeta_fd is a central "
             "finite-difference estimate of dbeta (magnitude perturbation); "
             "dbeta_rel_err = |dbeta_fd - dbeta| / |dbeta| (NaN at |v| = 0).");
      header(out, "columns",
             "tauc,vx,vy,|v|,beta,dbeta,dbeta_fd,dbeta_rel_err,|tau_b|");

      for (double tauc : m_taucs) {
        for (double speed : m_speeds) {
          for (const auto &dir : m_directions) {
            const double vx = speed * dir.first, vy = speed * dir.second;
            double beta = dl.law->drag(tauc, vx, vy);
            double dbeta = 0.0;
            dl.law->drag_with_derivative(tauc, vx, vy, &beta, &dbeta);

            // central difference of dbeta/d(alpha), alpha = 0.5 |v|^2, by
            // perturbing the speed magnitude
            double beta_fd = std::numeric_limits<double>::quiet_NaN();
            double rel_err = std::numeric_limits<double>::quiet_NaN();
            if (speed > 0.0) {
              const double eps = 1.0e-5;
              const double v_p = speed * (1.0 + eps), v_m = speed * (1.0 - eps);
              const double bp = dl.law->drag(tauc, v_p * dir.first, v_p * dir.second);
              const double bm = dl.law->drag(tauc, v_m * dir.first, v_m * dir.second);
              const double alpha_p = 0.5 * v_p * v_p, alpha_m = 0.5 * v_m * v_m;
              beta_fd = (bp - bm) / (alpha_p - alpha_m);
              rel_err = std::fabs(beta_fd - dbeta) / std::max(1.0, std::fabs(dbeta));
            }

            const double mag_v = std::sqrt(vx * vx + vy * vy);
            out << tauc << "," << vx << "," << vy << "," << mag_v << "," << beta
                << "," << dbeta << "," << beta_fd << "," << rel_err << ","
                << beta * mag_v << "\n";
          }
        }
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // (c) tauc map on the test-F exact geometry (MohrCoulombYieldStress)
    // ------------------------------------------------------------------
    auto grid = Grid::Shallow(ctx, 900.0e3, 900.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::NOT_PERIODIC);
    const double LforFG = 750000.0;

    Geometry geometry(grid);
    geometry.sea_level_elevation.set(0.0);
    geometry.bed_elevation.set(0.0);
    geometry.cell_type.set(MASK_GROUNDED);

    array::Scalar1 W_till(grid, "till_water_thickness");
    array::Scalar1 W_subglacial(grid, "subglacial_water_thickness");
    array::Scalar1 no_model(grid, "no_model_mask");
    array::Scalar1 phi(grid, "tillphi");

    // thickness = exact test-F H (flat bed, surface = H); deterministic smooth
    // fields: W_till = Wmax (r/L)^2 (0 at the dome -> saturated at the margin)
    // and phi = 20 + 10 cos(pi r/L) deg (30 at the dome -> 10 at the margin).
    {
      array::AccessScope list{&geometry.ice_thickness, &W_till, &phi};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        geometry.ice_thickness(i, j) =
            (r > LforFG - 1.0) ? 0.0 : exactFG(0.0, r, grid->z(), 0.0).H;
        const double sr = std::min(1.0, r / LforFG);
        W_till(i, j) = W_till_max * sr * sr;
        phi(i, j) = 20.0 + 10.0 * std::cos(M_PI * sr);
      }
      geometry.ice_thickness.update_ghosts();
      geometry.ice_surface_elevation.copy_from(geometry.ice_thickness);
      geometry.ice_surface_elevation.update_ghosts();
      W_till.update_ghosts();
      phi.update_ghosts();
    }
    geometry.ensure_consistency(config->get_number("geometry.ice_free_thickness_standard"));
    W_subglacial.set(0.0);
    no_model.set(0.0);

    YieldStressInputs ys_inputs;
    ys_inputs.geometry = &geometry;
    ys_inputs.till_water_thickness = &W_till;
    ys_inputs.subglacial_water_thickness = &W_subglacial;
    ys_inputs.no_model_mask = &no_model;

    auto mcys = std::make_shared<MohrCoulombYieldStress>(grid);
    mcys->init(ys_inputs); // initializes tauc with the default tillphi (30 deg)
    mcys->set_till_friction_angle(phi);
    mcys->update(ys_inputs, 0.0, 1.0);

    const array::Scalar &tauc = mcys->basal_material_yield_stress();

    // effective pressure field, recomputed pointwise with the same inputs and
    // config (update_impl does not store it; no slippery_grounding_lines /
    // add_transportable_water here, so water = W_till)
    const double ice_rho = config->get_number("constants.ice.density"),
                 g_std = config->get_number("constants.standard_gravity");
    array::Scalar1 N_till_field(grid, "N_till");
    {
      array::AccessScope list{&geometry.ice_thickness, &W_till, &N_till_field};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double P_overburden = ice_rho * g_std * geometry.ice_thickness(i, j);
        N_till_field(i, j) = mc.effective_pressure(delta_config, P_overburden, W_till(i, j));
      }
    }

    // NetCDF: inputs + outputs
    {
      const std::vector<const array::Array *> vecs = {
          &geometry.ice_surface_elevation, &geometry.ice_thickness,
          &geometry.bed_elevation,         &geometry.cell_type,
          &geometry.sea_level_elevation,   &W_till,
          &W_subglacial,                   &no_model,
          &phi,                            &tauc,
          &N_till_field};
      write_netcdf(ctx, out_dir + "/inputs.nc", vecs);
    }

    // per-cell tauc map CSV
    {
      std::ofstream out = open_output(out_dir + "/tauc_map.csv");
      header(out, "description",
             "MohrCoulombYieldStress::update on the test-F exact geometry: per "
             "cell: geometry, prescribed W_till and phi, P_overburden, "
             "recomputed N_till, and the computed yield stress tauc");
      header(out, "columns",
             "i,j,x,y,r,H,surface,bed,mask,W_till,phi_deg,delta,P_overburden,N_till,tauc");

      array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                              &geometry.bed_elevation, &geometry.cell_type, &W_till,
                              &phi, &tauc, &N_till_field};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        const double P_overburden = ice_rho * g_std * geometry.ice_thickness(i, j);
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
            << geometry.ice_thickness(i, j) << "," << geometry.ice_surface_elevation(i, j)
            << "," << geometry.bed_elevation(i, j) << "," << (int)geometry.cell_type(i, j)
            << "," << W_till(i, j) << "," << phi(i, j) << "," << delta_config << ","
            << P_overburden << "," << N_till_field(i, j) << "," << tauc(i, j) << "\n";
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // Validation summaries
    // ------------------------------------------------------------------

    // (a) Mohr-Coulomb pointwise checks
    {
      const double N_W0 = mc.effective_pressure(0.02, 1.0e7, 0.0);
      const double N_Wmax = mc.effective_pressure(0.02, 1.0e7, W_till_max);
      const double N_Whalf = mc.effective_pressure(0.02, 1.0e7, 0.5 * W_till_max);
      const double tauc_W0 = mc.yield_stress(0.02, 1.0e7, 0.0, 30.0);
      const double tauc_Wmax = mc.yield_stress(0.02, 1.0e7, W_till_max, 30.0);

      // monotonicity: N_till non-increasing in water_thickness
      bool monotone = true;
      double max_roundtrip_err = 0.0;
      for (double delta : m_deltas) {
        for (double P : m_P_overburden) {
          double prev = std::numeric_limits<double>::infinity();
          for (double W : m_water_thickness) {
            const double N = mc.effective_pressure(delta, P, W);
            if (N > prev + 1.0e-9 * std::max(1.0, std::fabs(N))) {
              monotone = false;
            }
            prev = N;
            for (double phi : m_phis) {
              const double tauc_r = mc.yield_stress(delta, P, W, phi);
              const double phi_inv = mc.till_friction_angle(delta, P, W, tauc_r);
              max_roundtrip_err =
                  std::max(max_roundtrip_err, std::fabs(phi_inv - phi));
            }
          }
        }
      }

      std::cout << "BASAL STRENGTH (Mohr-Coulomb pointwise) validation:\n"
                << "  N_till(0.02, 1e7 Pa, W=0)       = " << N_W0
                << " Pa (must equal P_overburden)\n"
                << "  N_till(0.02, 1e7 Pa, W=Wmax)    = " << N_Wmax
                << " Pa (must equal delta*P = " << 0.02 * 1.0e7 << ")\n"
                << "  N_till(0.02, 1e7 Pa, W=Wmax/2)  = " << N_Whalf
                << " Pa (still saturated at P = 1e7)\n"
                << "  tauc(0.02, 1e7 Pa, W=0, phi=30) = " << tauc_W0
                << " Pa (= tan(30 deg) * 1e7 + c0)\n"
                << "  tauc(0.02, 1e7 Pa, W=Wmax, phi=30) = " << tauc_Wmax
                << " Pa (= tan(30 deg) * 2e5 + c0)\n"
                << "  monotone non-increasing in W: " << (monotone ? "OK" : "FAIL") << "\n"
                << "  max |phi_inverse(tauc(phi)) - phi| = " << max_roundtrip_err
                << " deg (round trip)\n";
    }

    // (b) basal resistance law checks: plastic limit |tau_b| -> tauc, and
    //     drag_with_derivative vs finite differences
    {
      const double reg = config->get_number("basal_resistance.plastic.regularization", "m second-1");
      for (const auto &dl : drag_laws) {
        double max_rel_err = 0.0, max_beta_spread = 0.0;
        // plastic limit at a fast speed
        double beta_fast = dl.law->drag(1.0e5, 1.0e-4, 0.0);
        for (double tauc : m_taucs) {
          for (double speed : m_speeds) {
            if (speed == 0.0) {
              continue;
            }
            for (const auto &dir : m_directions) {
              const double vx = speed * dir.first, vy = speed * dir.second;
              double beta = dl.law->drag(tauc, vx, vy);
              double dbeta = 0.0;
              dl.law->drag_with_derivative(tauc, vx, vy, &beta, &dbeta);
              const double eps = 1.0e-5;
              const double v_p = speed * (1.0 + eps), v_m = speed * (1.0 - eps);
              const double bp = dl.law->drag(tauc, v_p * dir.first, v_p * dir.second);
              const double bm = dl.law->drag(tauc, v_m * dir.first, v_m * dir.second);
              const double alpha_p = 0.5 * v_p * v_p, alpha_m = 0.5 * v_m * v_m;
              const double fd = (bp - bm) / (alpha_p - alpha_m);
              max_rel_err = std::max(max_rel_err, std::fabs(fd - dbeta) /
                                                      std::max(1.0, std::fabs(dbeta)));
            }
            // direction invariance: beta must not depend on the direction
            double b0 = dl.law->drag(tauc, speed, 0.0);
            for (const auto &dir : m_directions) {
              const double vx = speed * dir.first, vy = speed * dir.second;
              const double b = dl.law->drag(tauc, vx, vy);
              max_beta_spread = std::max(max_beta_spread, std::fabs(b - b0));
            }
          }
        }
        const double beta_v0 = dl.law->drag(1.0e5, 0.0, 0.0);
        std::cout << "  [" << dl.name << "] beta(1e5, 0, 0) = " << beta_v0
                  << " Pa s m^-1 (regularized at v=0; plastic: tauc/eps = "
                  << 1.0e5 / reg << ")\n"
                  << "  [" << dl.name << "] beta(1e5, 1e-4 m/s, 0) = " << beta_fast
                  << " Pa s m^-1 (|tau_b| = beta*v -> tauc at high speed)\n"
                  << "  [" << dl.name << "] max |dbeta_fd - dbeta|/|dbeta| = "
                  << max_rel_err << "\n"
                  << "  [" << dl.name << "] max beta spread across directions = "
                  << max_beta_spread << " (must be ~0)\n";
      }
    }

    // (c) tauc map checks
    {
      const double high_tauc = config->get_number("basal_yield_stress.ice_free_bedrock");
      double min_tauc = std::numeric_limits<double>::infinity(),
             max_tauc = -std::numeric_limits<double>::infinity(),
             max_N = -std::numeric_limits<double>::infinity(),
             min_N = std::numeric_limits<double>::infinity();
      bool ice_free_ok = true;
      double dome_tauc = 0.0, dome_N = 0.0, dome_P = 0.0;

      array::AccessScope list{&geometry.ice_thickness, &geometry.cell_type, &W_till,
                              &phi, &tauc, &N_till_field};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double H = geometry.ice_thickness(i, j);
        if (geometry.cell_type.ice_free(i, j)) {
          if (tauc(i, j) != high_tauc) {
            ice_free_ok = false;
          }
          continue;
        }
        min_tauc = std::min(min_tauc, tauc(i, j));
        max_tauc = std::max(max_tauc, tauc(i, j));
        min_N = std::min(min_N, N_till_field(i, j));
        max_N = std::max(max_N, N_till_field(i, j));
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        if (r < 1000.0) { // dome-center cell
          dome_tauc = tauc(i, j);
          dome_N = N_till_field(i, j);
          dome_P = ice_rho * g_std * H;
        }
      }
      // the dome-center tauc must equal c0 + tan(phi) P (W_till ~ 0 there)
      std::cout << "TAUC MAP (MohrCoulombYieldStress, test-F geometry) validation:\n"
                << "  ice-free cells all = " << high_tauc << " Pa: "
                << (ice_free_ok ? "OK" : "FAIL") << "\n"
                << "  grounded tauc range: [" << min_tauc << ", " << max_tauc << "] Pa\n"
                << "  N_till range: [" << min_N << ", " << max_N << "] Pa\n"
                << "  dome-center: P = " << dome_P << " Pa, N_till = " << dome_N
                << " Pa (must equal P: W_till ~ 0), tauc = " << dome_tauc
                << " Pa (must equal c0 + tan(phi) N_till)\n";
    }

    std::cout << "Basal-strength instrumentation complete. Dumps in " << out_dir
              << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
