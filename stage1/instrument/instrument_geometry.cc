// Stage 1 instrumentation driver for the PISM geometry (mass continuity)
// boundary: GeometryEvolution::flow_step, GeometryEvolution::source_term_step,
// and the pointwise part_grid_threshold_thickness function.
//
// Runs the two integrated (time-stepping) steps on the verification test-F
// exact geometry (flat bed, thickness = exact FG H, sea level = 0, Mx x My =
// 31 x 31, Lx = Ly = 900 km -- all cells grounded or ice-free bedrock), each
// over the same dt = 1 year, from the SAME initial geometry (as the two steps
// are independent in the model: the real coupled model runs flow_step first
// and then source_term_step on the *updated* geometry, but the boundary
// contract of each step is a trace from an arbitrary input geometry).
//
// Prescribed deterministic forcing (all exactly specified):
//
//   advective_velocity = the exact test-F surface velocity
//                        (u, v) = U_surf(r) * (x/r, y/r),  U_surf(r) =
//                        exactFG(0, r, z, 0).U[top level]  -- the analytic
//                        surface value of the exact SIA velocity field (the
//                        "SSA-like" horizontal velocity of the mass continuity
//                        scheme).  Zero at/outside the ice margin r >= L.
//                        Note: test F is a cold no-slide SIA solution, so
//                        these speeds are small (max ~2.5 m/yr); the advective
//                        part of the flow_step trace is therefore a small
//                        correction, while the diffusive part carries the
//                        visible signal.
//   diffusive_flux     = a prescribed, exactly-specified Gaussian radial flux
//                        on the staggered grid:
//                          q(r) = q0 * (r/Rq) * exp(-(r/Rq)^2) * (x/r, y/r),
//                        q0 = 2e-2 m^2 s^-1, Rq = 300 km, sampled at the edge
//                        midpoints (q_x at the east edge, q_y at the north
//                        edge of each cell, matching the SIA diffusive-flux
//                        layout).  This is a documented surrogate for the real
//                        test-F SIA flux (max |q_SIA| ~ 2e-2 m^2 s^-1 at this
//                        resolution, so the magnitude is calibrated); its
//                        divergence is analytic:
//                          div q = (2 q0 / Rq) exp(-s) (1 - s),  s = (r/Rq)^2
//                        so dH/dt = -dt div q is a clean .esm test target.
//   smb_rate           = the exact test-F surface mass balance M(r) * rho_ice
//                        (kg m^-2 s^-1; exactFG returns M in m s^-1).
//   basal_melt_rate    = 0.05 m/yr * exp(-(r/200 km)^2)  (ice-equivalent m/s)
//   thickness_bc_mask  = 0 everywhere (no Dirichlet thickness B.C.s)
//   dt                 = 1 year
//
// flow_step contract: H_change + Href_change = dt*(SMB + BMB - flux_divergence)
// + conservation_error (here SMB = BMB = 0).  Because the outermost ring of the
// domain is ice-free, the flux limiters zero every flux on the outer boundary,
// so sum(div Q) telescopes to 0 exactly and the domain-sum of H_change equals
// sum(conservation_error) to round-off -- the "closed domain" mass-balance
// check.  source_term_step contract: H_change = dH_SMB + dH_BMB with
// dH_SMB = dt*smb/rho_ice, dH_BMB = -dt*bmr (clamped by effective_change()).
//
// Dumps (stage1/dumps/geometry/):
//
//   meta.txt          - PISM revision, config file, grid, dt, forcing, command
//   parameters.csv    - geometry.*, constants.*, grid.* config keys
//   flow_inputs.nc    - flow_step: geometry in, velocity, diffusive flux, bc
//                       mask, flux (after limiting), flux divergence, thickness
//                       change, Href change, conservation error, geometry out
//   source_inputs.nc  - source_term_step: geometry in, bc mask, smb rate, basal
//                       melt rate, effective SMB/BMB (m), geometry out
//   columns.csv       - per cell for BOTH steps:
//                       (i,j,x,y,r,H_in,surface_in,bed,mask_in,
//                        u,v,qx,qy,smb_rate,basal_melt_rate,
//                        H_out_flow,surface_out_flow,mask_out_flow,dH_flow,
//                        dH_dt_flow_myr,H_out_source,surface_out_source,
//                        mask_out_source,dH_source,dH_dt_source_myr)
//                       with qx = diffusive flux x-component on the east edge,
//                       qy = y-component on the north edge (staggered layout)
//   pointwise.csv     - part_grid_threshold_thickness sampled over
//                       (cell_type star, thickness star, surface star, bed)
//                       covering the branches: all icy (grounded / floating /
//                       mixed), 3 icy + 1 ice-free, 1 icy + 3 ice-free, all
//                       ice-free (ocean and bedrock); thin/thick and
//                       bed-above/below-surface cases.  The function averages
//                       H and h over the icy N/S/E/W neighbors and returns
//                       max(min(h_avg - bed, H_avg), 0) (0 if no icy
//                       neighbors) -- an analytic .esm test target.
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/geometry.json

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
#include "pism/geometry/Geometry.hh"
#include "pism/geometry/GeometryEvolution.hh"
#include "pism/geometry/part_grid_threshold_thickness.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM geometry (mass continuity) boundary.\n"
  "Runs GeometryEvolution::flow_step and source_term_step plus the pointwise\n"
  "part_grid_threshold_thickness on the test-F exact geometry and dumps the I/O\n"
  "traces to stage1/dumps/geometry/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

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
         "Configuration parameters affecting the GeometryEvolution component "
         "(geometry.*, constants.*, grid.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {"geometry.", "constants.", "grid."};

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
// Test-F exact-geometry setup (same as the hydrology / basal-strength drivers)
// ---------------------------------------------------------------------------
void set_testF_geometry(const std::shared_ptr<const Grid> &grid, Geometry &geometry,
                        double ice_free_thickness_standard) {
  const double LforFG = 750000.0;
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
  geometry.ensure_consistency(ice_free_thickness_standard);
}

// Copy one Geometry into another (field by field; the array classes do not
// define a public deep copy).
void copy_geometry(const Geometry &from, Geometry &to) {
  to.sea_level_elevation.copy_from(from.sea_level_elevation);
  to.bed_elevation.copy_from(from.bed_elevation);
  to.ice_thickness.copy_from(from.ice_thickness);
  to.ice_area_specific_volume.copy_from(from.ice_area_specific_volume);
  to.cell_type.copy_from(from.cell_type);
  to.ice_surface_elevation.copy_from(from.ice_surface_elevation);
}

} // namespace
} // namespace pism

int main(int argc, char *argv[]) {

  using namespace pism;
  using namespace pism::mask;

  MPI_Comm com = MPI_COMM_WORLD;
  petsc::Initializer petsc(argc, argv, help);
  com = MPI_COMM_WORLD;

  try {
    std::shared_ptr<Context> ctx = context_from_options(com, "geometry_test");
    auto config = ctx->config();

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/geometry";
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

    // grid: test-F domain (2-D, same as the hydrology driver)
    auto grid = Grid::Shallow(ctx, 900.0e3, 900.0e3, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::NOT_PERIODIC);

    // ------------------------------------------------------------------
    // config keys that enter the geometry physics
    // ------------------------------------------------------------------
    const double ice_rho = config->get_number("constants.ice.density");
    const double ice_free_thickness =
        config->get_number("geometry.ice_free_thickness_standard");
    const double myr = 31556926.0; // seconds per year (365.2422 days, matches exactFG)

    // forcing constants (exactly specified; documented in the run config)
    const double LforFG = 750000.0;
    const double q0 = 2.0e-2;  // m^2 s^-1, diffusive flux amplitude
    const double Rq = 300.0e3; // m, diffusive flux length scale
    const double bmr0 = units::convert(ctx->unit_system(), 0.05, "m / year", "m / second");
    const double Rb = 200.0e3; // m, basal melt blob radius

    // ------------------------------------------------------------------
    // meta + parameters
    // ------------------------------------------------------------------
    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "component",
             "geometry (GeometryEvolution: flow_step + source_term_step + "
             "part_grid_threshold_thickness)");
      header(meta, "description",
             "One dt = 1 yr step each of GeometryEvolution::flow_step (advective + "
             "diffusive mass transport) and source_term_step (surface + basal mass "
             "balance) from the SAME test-F exact geometry (flat bed, H = exact FG, "
             "sea level 0, grounded), plus a pointwise sample of "
             "part_grid_threshold_thickness. Forcing: advective_velocity = exact "
             "test-F surface velocity U_surf(r) (x/r, y/r); diffusive_flux = "
             "prescribed Gaussian radial flux q = q0 (r/Rq) exp(-(r/Rq)^2) e_r "
             "(q0 = 2e-2 m^2/s, Rq = 300 km, staggered edge midpoints); smb_rate = "
             "exact test-F M(r) * rho_ice; basal_melt_rate = 0.05 m/yr "
             "exp(-(r/200 km)^2); thickness_bc_mask = 0 everywhere");
      header(meta, "dt_seconds", std::to_string(dt));
      header(meta, "dt_years", std::to_string(dt_years));
      header(meta, "Mx", std::to_string(Mx));
      header(meta, "My", std::to_string(My));
      header(meta, "Lx", "900000");
      header(meta, "Ly", "900000");
      header(meta, "dx", std::to_string(grid->dx()));
      header(meta, "cell_area", std::to_string(grid->cell_area()));
      meta << "# forcing: advective_velocity = exact test-F surface velocity "
              "U_surf(r)*(x/r,y/r)\n";
      meta << "# forcing: diffusive_flux = q0*(r/Rq)*exp(-(r/Rq)^2)*(x/r,y/r), "
              "q0 = 2e-2 m^2/s, Rq = 300 km (staggered)\n";
      meta << "# forcing: smb_rate = rho_ice * exactFG(...).M (kg m^-2 s^-1)\n";
      meta << "# forcing: basal_melt_rate = 0.05 m/yr * exp(-(r/200 km)^2) (m/s)\n";
      meta << "# forcing: thickness_bc_mask = 0 everywhere\n";
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    // ------------------------------------------------------------------
    // Test-F exact geometry (shared by both steps)
    // ------------------------------------------------------------------
    Geometry geometry(grid);
    set_testF_geometry(grid, geometry, ice_free_thickness);

    // ------------------------------------------------------------------
    // Inputs: deterministic forcing fields
    // ------------------------------------------------------------------
    // thickness Dirichlet B.C. mask: no B.C.s
    array::Scalar1 bc_mask(grid, "thk_bc_mask");
    bc_mask.set(0.0);

    // advective velocity: exact test-F surface velocity (radial)
    array::Vector1 velocity(grid, "advective_velocity");
    {
      array::AccessScope list{&velocity};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double x = grid->x(i), y = grid->y(j);
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        double U_surf = 0.0;
        if (r < LforFG) {
          U_surf = exactFG(0.0, r, grid->z(), 0.0).U.back(); // top level = surface
        }
        const double inv_r = (r > 1.0) ? 1.0 / r : 0.0;
        velocity(i, j).u = U_surf * inv_r * x;
        velocity(i, j).v = U_surf * inv_r * y;
      }
    }
    velocity.update_ghosts();

    // diffusive flux: prescribed Gaussian radial field on the staggered grid
    // (q_x at the east edge midpoint, q_y at the north edge midpoint of each cell)
    array::Staggered1 diffusive_flux(grid, "diffusive_flux");
    {
      array::AccessScope list{&diffusive_flux};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        for (int o = 0; o < 2; ++o) {
          const double xs = grid->x(i) + (o == 0 ? 0.5 * grid->dx() : 0.0);
          const double ys = grid->y(j) + (o == 1 ? 0.5 * grid->dy() : 0.0);
          const double rs = std::max(std::hypot(xs, ys), 1.0);
          const double s = (rs / Rq) * (rs / Rq);
          const double qmag = q0 * (rs / Rq) * std::exp(-s);
          diffusive_flux(i, j, o) = (o == 0) ? qmag * (xs / rs) : qmag * (ys / rs);
        }
      }
    }
    diffusive_flux.update_ghosts();

    // surface mass balance: exact test-F M(r) in kg m^-2 s^-1
    array::Scalar1 smb_rate(grid, "surface_mass_balance_rate");
    {
      array::AccessScope list{&smb_rate};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        smb_rate(i, j) =
            (r < LforFG) ? exactFG(0.0, r, grid->z(), 0.0).M * ice_rho : 0.0;
      }
    }

    // basal melt rate: dome-centered blob, 0.05 m/yr peak (ice-equivalent m/s)
    array::Scalar1 basal_melt_rate(grid, "basal_melt_rate");
    {
      array::AccessScope list{&basal_melt_rate};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        basal_melt_rate(i, j) = bmr0 * std::exp(-(r / Rb) * (r / Rb));
      }
    }

    // ------------------------------------------------------------------
    // The component
    // ------------------------------------------------------------------
    GeometryEvolution geom_evo(grid);
    geom_evo.init(InputOptions(INIT_OTHER, "", 0));

    // ==================================================================
    // (a) flow_step trace
    // ==================================================================
    geom_evo.reset();
    geom_evo.flow_step(geometry, dt, velocity, diffusive_flux, bc_mask);

    const array::Scalar &dH_flow = geom_evo.thickness_change_due_to_flow();
    const array::Scalar &dHref_flow = geom_evo.area_specific_volume_change_due_to_flow();
    const array::Scalar &cons_err = geom_evo.conservation_error();
    const array::Scalar &flux_div = geom_evo.flux_divergence();
    const array::Staggered1 &flux_stag = geom_evo.flux_staggered();

    // geometry out = geometry in + changes, then consistent mask/surface
    Geometry geometry_flow_out(grid);
    copy_geometry(geometry, geometry_flow_out);
    geometry_flow_out.ice_thickness.add(1.0, dH_flow);
    geometry_flow_out.ice_area_specific_volume.add(1.0, dHref_flow);
    geometry_flow_out.ensure_consistency(ice_free_thickness);
    // distinct NetCDF names (the *_out arrays carry the geometry-in names otherwise)
    geometry_flow_out.ice_thickness.metadata().set_name("thk_flow_out");
    geometry_flow_out.ice_surface_elevation.metadata().set_name("usurf_flow_out");
    geometry_flow_out.cell_type.metadata().set_name("mask_flow_out");

    {
      const std::vector<const array::Array *> vecs = {
          &geometry.ice_surface_elevation,  &geometry.ice_thickness,
          &geometry.bed_elevation,          &geometry.cell_type,
          &geometry.sea_level_elevation,    &velocity,
          &diffusive_flux,                  &bc_mask,
          &flux_stag,                       &flux_div,
          &dH_flow,                         &dHref_flow,
          &cons_err,                        &geometry_flow_out.ice_thickness,
          &geometry_flow_out.ice_surface_elevation, &geometry_flow_out.cell_type};
      write_netcdf(ctx, out_dir + "/flow_inputs.nc", vecs);
    }

    // ==================================================================
    // (b) source_term_step trace
    // ==================================================================
    geom_evo.source_term_step(geometry, dt, bc_mask, smb_rate, basal_melt_rate);

    const array::Scalar &dH_SMB = geom_evo.top_surface_mass_balance();
    const array::Scalar &dH_BMB = geom_evo.bottom_surface_mass_balance();

    Geometry geometry_source_out(grid);
    copy_geometry(geometry, geometry_source_out);
    geom_evo.apply_mass_fluxes(geometry_source_out);
    geometry_source_out.ensure_consistency(ice_free_thickness);
    geometry_source_out.ice_thickness.metadata().set_name("thk_source_out");
    geometry_source_out.ice_surface_elevation.metadata().set_name("usurf_source_out");
    geometry_source_out.cell_type.metadata().set_name("mask_source_out");

    {
      const std::vector<const array::Array *> vecs = {
          &geometry.ice_surface_elevation,  &geometry.ice_thickness,
          &geometry.bed_elevation,          &geometry.cell_type,
          &geometry.sea_level_elevation,    &bc_mask,
          &smb_rate,                        &basal_melt_rate,
          &dH_SMB,                          &dH_BMB,
          &geometry_source_out.ice_thickness, &geometry_source_out.ice_surface_elevation,
          &geometry_source_out.cell_type};
      write_netcdf(ctx, out_dir + "/source_inputs.nc", vecs);
    }

    // ==================================================================
    // per-cell columns (both steps)
    // ==================================================================
    {
      std::ofstream out = open_output(out_dir + "/columns.csv");
      header(out, "description",
             "One dt step each of flow_step and source_term_step on the test-F "
             "exact geometry: per cell: geometry in (H, surface, bed, mask), the "
             "prescribed forcing (u, v, qx = diffusive flux x-comp on the east "
             "edge, qy = y-comp on the north edge, smb_rate in kg m^-2 s^-1, "
             "basal_melt_rate in ice-equivalent m/s), and the geometry + "
             "thickness change after each step. dH/dt is in m/yr (positive = "
             "thickening).");
      header(out, "columns",
             "i,j,x,y,r,H_in,surface_in,bed,mask_in,u,v,qx,qy,smb_rate,"
             "basal_melt_rate,H_out_flow,surface_out_flow,mask_out_flow,dH_flow,"
             "dH_dt_flow_myr,H_out_source,surface_out_source,mask_out_source,"
             "dH_source,dH_dt_source_myr");

      array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                              &geometry.bed_elevation, &geometry.cell_type,
                              &velocity, &diffusive_flux, &smb_rate, &basal_melt_rate,
                              &dH_flow, &geometry_flow_out.ice_thickness,
                              &geometry_flow_out.ice_surface_elevation,
                              &geometry_flow_out.cell_type,
                              &dH_SMB, &dH_BMB,
                              &geometry_source_out.ice_thickness,
                              &geometry_source_out.ice_surface_elevation,
                              &geometry_source_out.cell_type};

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        const double dH_source = dH_SMB(i, j) + dH_BMB(i, j);
        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
            << geometry.ice_thickness(i, j) << ","
            << geometry.ice_surface_elevation(i, j) << ","
            << geometry.bed_elevation(i, j) << "," << (int)geometry.cell_type(i, j) << ","
            << velocity(i, j).u << "," << velocity(i, j).v << ","
            << diffusive_flux(i, j, 0) << "," << diffusive_flux(i, j, 1) << ","
            << smb_rate(i, j) << "," << basal_melt_rate(i, j) << ","
            << geometry_flow_out.ice_thickness(i, j) << ","
            << geometry_flow_out.ice_surface_elevation(i, j) << ","
            << (int)geometry_flow_out.cell_type(i, j) << ","
            << dH_flow(i, j) << "," << dH_flow(i, j) * myr / dt << ","
            << geometry_source_out.ice_thickness(i, j) << ","
            << geometry_source_out.ice_surface_elevation(i, j) << ","
            << (int)geometry_source_out.cell_type(i, j) << ","
            << dH_source << "," << dH_source * myr / dt << "\n";
      }
      out.close();
    }

    // ==================================================================
    // (c) part_grid_threshold_thickness pointwise sampling
    // ==================================================================
    {
      std::ofstream out = open_output(out_dir + "/pointwise.csv");
      header(out, "description",
             "part_grid_threshold_thickness(cell_type, thickness, surface, bed): "
             "a pointwise sample over deterministic (mask star, thickness star, "
             "surface star, bed) inputs. The function averages H and h over the "
             "icy N/S/E/W neighbors (mask 2 or 3) and returns "
             "max(min(h_avg - bed, H_avg), 0), or 0 if there are no icy "
             "neighbors. Columns: pattern name, number of icy neighbors N, the "
             "5-point mask/thickness/surface stars (c,n,e,s,w), bed, H_avg, "
             "h_avg, and the computed threshold (m).");
      header(out, "columns",
             "pattern,N,mask_c,mask_n,mask_e,mask_s,mask_w,H_c,H_n,H_e,H_s,H_w,"
             "h_c,h_n,h_e,h_s,h_w,bed,H_avg,h_avg,threshold");

      struct Pattern {
        std::string name;
        int m[5]; // c, n, e, s, w
      };
      const std::vector<Pattern> patterns = {
          {"all_grounded", {MASK_GROUNDED, MASK_GROUNDED, MASK_GROUNDED, MASK_GROUNDED,
                            MASK_GROUNDED}},
          {"all_floating", {MASK_FLOATING, MASK_FLOATING, MASK_FLOATING, MASK_FLOATING,
                            MASK_FLOATING}},
          {"grounded_floating_ocean",
           {MASK_GROUNDED, MASK_FLOATING, MASK_GROUNDED, MASK_GROUNDED, MASK_ICE_FREE_OCEAN}},
          {"3_icy_1_ocean", {MASK_GROUNDED, MASK_GROUNDED, MASK_GROUNDED, MASK_GROUNDED,
                             MASK_ICE_FREE_OCEAN}},
          {"1_icy_3_icefree", {MASK_GROUNDED, MASK_GROUNDED, MASK_ICE_FREE_BEDROCK,
                               MASK_ICE_FREE_OCEAN, MASK_ICE_FREE_BEDROCK}},
          {"all_ice_free_ocean", {MASK_ICE_FREE_OCEAN, MASK_ICE_FREE_OCEAN, MASK_ICE_FREE_OCEAN,
                                  MASK_ICE_FREE_OCEAN, MASK_ICE_FREE_OCEAN}},
          {"all_ice_free_land", {MASK_ICE_FREE_BEDROCK, MASK_ICE_FREE_BEDROCK, MASK_ICE_FREE_BEDROCK,
                                 MASK_ICE_FREE_BEDROCK, MASK_ICE_FREE_BEDROCK}}};

      const std::vector<double> Hs = {0.1, 100.0, 1000.0, 3000.0};
      const std::vector<double> hs = {-500.0, 0.0, 500.0, 2000.0, 4000.0};
      const std::vector<double> beds = {-1000.0, 0.0, 500.0, 2000.0};

      for (const auto &pat : patterns) {
        for (double H : Hs) {
          for (double h : hs) {
            for (double bed : beds) {
              // set all five star points; the function only uses the N/S/E/W
              auto set_star = [](auto &star, const auto *vals) {
                star.c = vals[0];
                star.n = vals[1];
                star.e = vals[2];
                star.s = vals[3];
                star.w = vals[4];
              };
              stencils::Star<int> m;
              stencils::Star<double> thk, surf;
              int mi[5];
              double Hv[5], hv[5];
              for (int k = 0; k < 5; ++k) {
                mi[k] = pat.m[k];
                const bool icy_cell = (mi[k] == MASK_GROUNDED or mi[k] == MASK_FLOATING);
                Hv[k] = icy_cell ? H : 0.0;
                hv[k] = icy_cell ? h : 0.0;
              }
              set_star(m, mi);
              set_star(thk, Hv);
              set_star(surf, hv);

              int N = 0;
              double H_avg = 0.0, h_avg = 0.0;
              for (auto d : {North, East, South, West}) {
                if (mask::icy(m[d])) {
                  H_avg += thk[d];
                  h_avg += surf[d];
                  N++;
                }
              }
              if (N > 0) {
                H_avg /= N;
                h_avg /= N;
              }

              const double threshold = part_grid_threshold_thickness(m, thk, surf, bed);

              out << pat.name << "," << N << "," << m.c << "," << m.n << "," << m.e << ","
                  << m.s << "," << m.w << "," << thk.c << "," << thk.n << "," << thk.e << ","
                  << thk.s << "," << thk.w << "," << surf.c << "," << surf.n << "," << surf.e
                  << "," << surf.s << "," << surf.w << "," << bed << "," << H_avg << ","
                  << h_avg << "," << threshold << "\n";
            }
          }
        }
      }
      out.close();
    }

    // ==================================================================
    // Validation
    // ==================================================================
    const double cell_area = grid->cell_area();

    // --- flow_step mass balance ---
    double sum_dH_flow = 0.0, sum_dHref_flow = 0.0, sum_flux_div = 0.0, sum_cons = 0.0;
    {
      array::AccessScope list{&dH_flow, &dHref_flow, &flux_div, &cons_err};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        sum_dH_flow += dH_flow(i, j);
        sum_dHref_flow += dHref_flow(i, j);
        sum_flux_div += flux_div(i, j);
        sum_cons += cons_err(i, j);
      }
      sum_dH_flow = GlobalSum(grid->com, sum_dH_flow);
      sum_dHref_flow = GlobalSum(grid->com, sum_dHref_flow);
      sum_flux_div = GlobalSum(grid->com, sum_flux_div);
      sum_cons = GlobalSum(grid->com, sum_cons);
    }
    const double sum_dt_divQ = -dt * sum_flux_div;
    const double flow_residual = sum_dH_flow + sum_dHref_flow - sum_dt_divQ - sum_cons;

    // --- source_term_step mass balance ---
    double sum_dSMB = 0.0, sum_dBMB = 0.0;
    double sum_analytic_smb = 0.0, sum_analytic_bmb = 0.0;
    long n_smb_clamped = 0, n_bmb_clamped = 0, n_src_bc = 0;
    {
      array::AccessScope list{&dH_SMB, &dH_BMB, &smb_rate, &basal_melt_rate, &bc_mask,
                              &geometry.cell_type, &geometry.ice_thickness};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        sum_dSMB += dH_SMB(i, j);
        sum_dBMB += dH_BMB(i, j);
        const double dsmb_analytic = dt * smb_rate(i, j) / ice_rho;
        const double dbmb_analytic = -dt * basal_melt_rate(i, j);
        sum_analytic_smb += dsmb_analytic;
        sum_analytic_bmb += dbmb_analytic;
        // clamp detection (same rule as effective_change):
        if (bc_mask(i, j) > 0.5 or geometry.cell_type.ice_free_ocean(i, j)) {
          n_src_bc++;
        } else {
          const double H = geometry.ice_thickness(i, j);
          if (H + dsmb_analytic < 0.0) {
            n_smb_clamped++;
          }
          if (H + dsmb_analytic + dbmb_analytic < 0.0) {
            n_bmb_clamped++;
          }
        }
      }
      sum_dSMB = GlobalSum(grid->com, sum_dSMB);
      sum_dBMB = GlobalSum(grid->com, sum_dBMB);
      sum_analytic_smb = GlobalSum(grid->com, sum_analytic_smb);
      sum_analytic_bmb = GlobalSum(grid->com, sum_analytic_bmb);
      n_smb_clamped = GlobalSum(grid->com, (int)n_smb_clamped);
      n_bmb_clamped = GlobalSum(grid->com, (int)n_bmb_clamped);
      n_src_bc = GlobalSum(grid->com, (int)n_src_bc);
    }
    const double sum_dH_source = sum_dSMB + sum_dBMB;
    const double source_truncation =
        (sum_analytic_smb - sum_dSMB) + (sum_analytic_bmb - sum_dBMB);

    // --- mask changes ---
    long n_mask_flow = 0, n_mask_source = 0;
    {
      array::AccessScope list{&geometry.cell_type, &geometry_flow_out.cell_type,
                              &geometry_source_out.cell_type};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        if ((int)geometry.cell_type(i, j) != (int)geometry_flow_out.cell_type(i, j)) {
          n_mask_flow++;
        }
        if ((int)geometry.cell_type(i, j) != (int)geometry_source_out.cell_type(i, j)) {
          n_mask_source++;
        }
      }
      n_mask_flow = GlobalSum(grid->com, (int)n_mask_flow);
      n_mask_source = GlobalSum(grid->com, (int)n_mask_source);
    }

    // --- field extremes ---
    double max_vel = 0.0, max_q_diff = 0.0, max_q_total = 0.0, min_smb = 1e300,
           max_smb = -1e300, max_dH_flow_abs = 0.0, max_dH_source_abs = 0.0;
    double min_H = 1e300, max_H = -1e300;
    {
      array::AccessScope list{&velocity, &diffusive_flux, &flux_stag, &smb_rate,
                              &dH_flow, &dH_SMB, &dH_BMB, &geometry.ice_thickness};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        max_vel = std::max(max_vel, std::hypot(velocity(i, j).u, velocity(i, j).v));
        max_q_diff = std::max(max_q_diff, std::fabs(diffusive_flux(i, j, 0)));
        max_q_diff = std::max(max_q_diff, std::fabs(diffusive_flux(i, j, 1)));
        max_q_total = std::max(max_q_total, std::fabs(flux_stag(i, j, 0)));
        max_q_total = std::max(max_q_total, std::fabs(flux_stag(i, j, 1)));
        min_smb = std::min(min_smb, smb_rate(i, j));
        max_smb = std::max(max_smb, smb_rate(i, j));
        max_dH_flow_abs = std::max(max_dH_flow_abs, std::fabs(dH_flow(i, j)));
        max_dH_source_abs =
            std::max(max_dH_source_abs, std::fabs(dH_SMB(i, j) + dH_BMB(i, j)));
        min_H = std::min(min_H, geometry.ice_thickness(i, j));
        max_H = std::max(max_H, geometry.ice_thickness(i, j));
      }
      max_vel = GlobalMax(grid->com, max_vel);
      max_q_diff = GlobalMax(grid->com, max_q_diff);
      max_q_total = GlobalMax(grid->com, max_q_total);
      min_smb = GlobalMin(grid->com, min_smb);
      max_smb = GlobalMax(grid->com, max_smb);
      max_dH_flow_abs = GlobalMax(grid->com, max_dH_flow_abs);
      max_dH_source_abs = GlobalMax(grid->com, max_dH_source_abs);
      min_H = GlobalMin(grid->com, min_H);
      max_H = GlobalMax(grid->com, max_H);
    }

    // --- center cell (dome) values ---
    const int ic = Mx / 2, jc = My / 2;
    double dome_dH_flow = 0.0, dome_dH_smb = 0.0, dome_dH_bmb = 0.0;
    double dome_divq_analytic = 0.0, dome_H = 0.0, dome_M_myr = 0.0;
    {
      array::AccessScope list{&dH_flow, &dH_SMB, &dH_BMB, &geometry.ice_thickness, &smb_rate};
      dome_dH_flow = dH_flow(ic, jc);
      dome_dH_smb = dH_SMB(ic, jc);
      dome_dH_bmb = dH_BMB(ic, jc);
      dome_H = geometry.ice_thickness(ic, jc);
      // analytic diffusive divergence at the dome (r -> 0): div q = 2 q0 / Rq
      dome_divq_analytic = 2.0 * q0 / Rq;
      dome_M_myr = smb_rate(ic, jc) / ice_rho * myr;
    }

    std::cout << "GEOMETRY (GeometryEvolution) validation:\n"
              << "  dt = " << dt << " s = " << dt_years << " year(s)\n"
              << "  grid " << Mx << " x " << My << ", dx = " << grid->dx() << " m, "
              << "cell_area = " << cell_area << " m^2\n"
              << "  H_in range [" << min_H << ", " << max_H << "] m\n"
              << "  forcing extremes: max |u| = " << max_vel << " m/s"
              << " (= " << max_vel * myr << " m/yr); max |q_diffusive| = " << max_q_diff
              << " m^2/s; max |q_total(limited)| = " << max_q_total
              << " m^2/s; smb_rate range [" << min_smb << ", " << max_smb
              << "] kg m^-2 s^-1\n"
              << "  dome cell (" << ic << "," << jc << "): H = " << dome_H << " m\n"
              << "  --- flow_step ---\n"
              << "    sum(H_change) = " << sum_dH_flow << " m (sum over cells)\n"
              << "    sum(Href_change) = " << sum_dHref_flow << " m\n"
              << "    sum(-dt*div Q) = " << sum_dt_divQ << " m\n"
              << "    sum(conservation_error) = " << sum_cons << " m\n"
              << "    contract residual = sum(dH + dHref - (-dt divQ) - cons) = "
              << flow_residual << " m (must be ~0)\n"
              << "    closed-domain mass balance: sum(-dt*div Q)*cell_area = "
              << sum_dt_divQ * cell_area << " m^3 (must be ~0: boundary fluxes are "
              << "limited to zero on the ice-free outer ring)\n"
              << "    dome cell dH_flow = " << dome_dH_flow
              << " m (analytic diffusive only: -dt*2*q0/Rq = " << -dt * dome_divq_analytic
              << " m; the advective part is ~0 at r = 0)\n"
              << "    mask changes: " << n_mask_flow << " cells\n"
              << "    max |dH_flow| = " << max_dH_flow_abs << " m\n"
              << "  --- source_term_step ---\n"
              << "    sum(dH_SMB) = " << sum_dSMB << " m; analytic dt*sum(smb/rho) = "
              << sum_analytic_smb << " m\n"
              << "    sum(dH_BMB) = " << sum_dBMB << " m; analytic -dt*sum(bmr) = "
              << sum_analytic_bmb << " m\n"
              << "    sum(H_change) = " << sum_dH_source
              << " m; truncation (clamped cells) = " << source_truncation
              << " m (" << n_smb_clamped << " SMB-clamped, " << n_bmb_clamped
              << " BMB-clamped, " << n_src_bc << " BC/ocean cells skipped)\n"
              << "    dome cell dH_SMB = " << dome_dH_smb << " m (= dt*M(r~0) = "
              << dt * dome_M_myr / myr << " m, M(0) = " << dome_M_myr << " m/yr); "
              << "dH_BMB = " << dome_dH_bmb << " m\n"
              << "    mask changes: " << n_mask_source << " cells\n"
              << "    max |dH_source| = " << max_dH_source_abs << " m\n";

    // part-grid pointwise cross-checks (a few hand-computed values)
    {
      // N = 0 (all ice-free): threshold must be 0
      stencils::Star<int> m0;
      m0.set(MASK_ICE_FREE_OCEAN);
      stencils::Star<double> thk0, surf0;
      thk0.set(1000.0);
      surf0.set(500.0);
      const double t0 = part_grid_threshold_thickness(m0, thk0, surf0, 0.0);
      // one grounded icy neighbor, H_avg = 1000, h_avg = 500, bed = -1000:
      // bed + H_avg = 0 > h_avg = 500? no -> threshold = H_avg = 1000
      stencils::Star<int> m1;
      m1.set(MASK_ICE_FREE_BEDROCK);
      m1.n = MASK_GROUNDED;
      stencils::Star<double> thk1, surf1;
      thk1.set(0.0);
      thk1.n = 1000.0;
      surf1.set(0.0);
      surf1.n = 500.0;
      const double t1 = part_grid_threshold_thickness(m1, thk1, surf1, -1000.0);
      // bed + H_avg = 0 > h_avg = 500? no -> H_avg = 1000
      // bed = 2000: bed + H_avg = 3000 > h_avg = 500 -> h_avg - bed = -1500 -> clamp 0
      const double t2 = part_grid_threshold_thickness(m1, thk1, surf1, 2000.0);
      std::cout << "  --- part_grid_threshold_thickness pointwise cross-checks ---\n"
                << "    all-ice-free-ocean neighbors: threshold = " << t0 << " (must be 0)\n"
                << "    one icy neighbor (H=1000,h=500), bed=-1000: threshold = " << t1
                << " (must be 1000 = min(h_avg-bed, H_avg) = min(1500,1000))\n"
                << "    one icy neighbor (H=1000,h=500), bed=2000: threshold = " << t2
                << " (must be 0 = max(min(-1500,1000),0))\n";
    }

    std::cout << "Geometry instrumentation complete. Dumps in " << out_dir << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
