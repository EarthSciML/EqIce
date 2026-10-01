// Stage 1 instrumentation driver for the atmosphere / ocean / surface forcing
// boundaries (the "coupler" components feeding the ice-sheet model).
//
// Constructs the default deterministic parameterization of each family and
// dumps its instantaneous outputs:
//
//   (a) the PIK atmosphere (atmosphere::PIK, "atmosphere.pik.parameterization
//       = martin"): air_temperature() (= the mean annual T_ma of the Martin et
//       al. 2011 parameterization, T_ma = 273.15 + 30 - 0.0075*usurf -
//       0.68775*lat*(-1)) and precipitation() (a constant-in-time field the
//       model normally reads from a file; the driver prescribes it) over the
//       domain, plus pointwise kernel samples (T_ma as a function of latitude
//       and surface elevation, the precipitation field, and the yearly-cycle
//       time series -- flat for the martin parameterization because
//       T_ms = T_ma, verified through the component's own temp_time_series).
//
//   (b) the PIK ocean (ocean::PIK, the Beckmann-Goosse-style sub-shelf
//       parameterization of Martin et al. 2011): shelf_base_temperature()
//       (= the pressure melting point T0 - beta_CC * rho_i * g * H), the
//       shelf_base_mass_flux() (proportional to the ocean-to-ice heat flux
//       Q = melt_factor * rho_w * c_p * gamma_T * (T_ocean - T_f), positive =
//       melting, linear in the shelf-base depth z_b = -(rho_i/rho_w) H), and
//       average_water_column_pressure() (= 0.5 rho_w g h_w^2 / H over the ice
//       column at the margin).  Dumps the maps and pointwise kernels over a
//       grid of depths, and the Constant ocean (trivial cross-check) maps.
//
//   (c) the PIK surface (surface::PIK): accumulation(), melt(), runoff(),
//       temperature() and mass_flux() maps.  surface::PIK holds the SMB
//       constant (from a file in real runs; prescribed here as
//       atmosphere.precipitation() minus a deterministic ablation blob) and
//       partitions it with dummy_accumulation/dummy_melt/dummy_runoff
//       (accumulation = max(smb, 0), melt = runoff = max(-smb, 0)).  Its
//       temperature uses the same Martin et al. 2011 formula as the PIK
//       atmosphere, so surface temperature == atmosphere air_temperature
//       field-for-field.
//
// Geometry: a radial ice dome (verification test-F exact thickness, extent
// 750 km) on a prescribed Gaussian ocean trough/cavity bed
// bed(r) = -1500 - 1500 exp(-((r - 450 km)/200 km)^2) m, sea level 0, so the
// dome center is grounded, the mid-radii are a floating shelf over a cavity
// whose depth varies 1500..3000 m, and the outer ring is ice-free ocean.  A
// prescribed latitude field spanning -80..-60 deg (linear in y) gives a clear
// latitudinal temperature gradient for the atmosphere/surface.
//
// Dumps (stage1/dumps/surface_ocean/):
//   meta.txt, parameters.csv, inputs.nc, columns.csv,
//   atmosphere_temperature.csv, atmosphere_precipitation.csv,
//   atmosphere_temperature_timeseries.csv, shelf_base_temperature.csv,
//   beckmann_goosse_mass_flux.csv, water_column_pressure.csv,
//   surface_partition.csv
//
// Build: see stage1/build/build_instrument.sh
// Run:   see stage1/runs/surface_ocean.json
//
// NOTE on the sign of shelf_base_mass_flux: the code path is
// combine_basal_melt_rate() -> geometry evolution with dH_BMB = -dt*mf/rho_i,
// so POSITIVE mass flux = melting.  The PIK mass flux is positive everywhere
// under the shelf; the comment in ConstantPIK.cc claiming it is "always
// negative" and "positive if ice is freezing on" is wrong (see bugs.md).

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
#include "pism/coupler/AtmosphereModel.hh"
#include "pism/coupler/OceanModel.hh"
#include "pism/coupler/SurfaceModel.hh"
#include "pism/coupler/atmosphere/PIK.hh"
#include "pism/coupler/ocean/ConstantPIK.hh"
#include "pism/coupler/ocean/Constant.hh"
#include "pism/coupler/surface/ConstantPIK.hh"
#include "pism/geometry/Geometry.hh"
#include "pism/verification/tests/exactTestsFG.hh"
#include "pism/util/io/SynchronousOutputWriter.hh"
#include "pism/util/io/io_helpers.hh"

static char help[] =
  "Stage-1 instrumentation driver for the PISM atmosphere/ocean/surface\n"
  "forcing boundaries. Dumps PIK atmosphere, PIK ocean (Beckmann-Goosse) and\n"
  "PIK surface traces plus the Constant ocean cross-check to\n"
  "stage1/dumps/surface_ocean/.\n\n";

namespace pism {

extern const char *revision; // defined in pism_config.cc (compiled in)

namespace {

// ---------------------------------------------------------------------------
// Subclasses exposing protected state so the driver can prescribe the fields
// the real models would read from an input file (no file needed; everything
// deterministic).
// ---------------------------------------------------------------------------

//! PIK atmosphere with the precipitation field prescribed by the driver
//! instead of read from a file.
class InstrumentedAtmospherePIK : public atmosphere::PIK {
public:
  InstrumentedAtmospherePIK(std::shared_ptr<const Grid> g) : atmosphere::PIK(g) {}

  array::Scalar &precip_field() { return m_precipitation; }

protected:
  void init_impl(const Geometry &geometry) override {
    (void) geometry; // precipitation is set by the driver; nothing to read
  }
};

//! PIK surface with the SMB field prescribed by the driver instead of read
//! from a file.
class InstrumentedSurfacePIK : public surface::PIK {
public:
  InstrumentedSurfacePIK(std::shared_ptr<const Grid> g,
                         std::shared_ptr<atmosphere::AtmosphereModel> atmosphere)
    : surface::PIK(g, atmosphere) {}

  array::Scalar &smb_field() { return *m_mass_flux; }

protected:
  void init_impl(const Geometry &geometry) override {
    (void) geometry; // SMB is set by the driver; nothing to read
  }
};

// ---------------------------------------------------------------------------
// Analytic mirror functions for the pointwise kernels (the .esm test targets).
// These reproduce the C++ formulas in coupler/{atmosphere/PIK.cc,
// ocean/ConstantPIK.cc, ocean/Constant.cc, ocean/OceanModel.cc,
// surface/ConstantPIK.cc} exactly; the component outputs are verified against
// them on the map.
// ---------------------------------------------------------------------------

// Martin et al. (2011), eq. (1): T_ma (K) as a function of surface elevation
// (m) and latitude (degrees north).
double martin2011_mean_annual(double elevation, double latitude) {
  return 273.15 + 30.0 - 0.0075 * elevation - 0.68775 * latitude * (-1.0);
}

// Pressure melting temperature of ice at depth H (m), from
// ocean::PIK::melting_point_temperature (also ocean::Constant).
double shelf_base_temperature(double H, double T0, double beta_CC, double rho_i, double g) {
  return T0 - beta_CC * rho_i * g * H;
}

// Beckmann-Goosse-style shelf base mass flux (kg m^-2 s^-1), from
// ocean::PIK::mass_flux.  Positive = melting (see file header).
double beckmann_goosse_mass_flux(double H, double melt_factor, double rho_i, double rho_w,
                                 double c_p_ocean, double gamma_T, double T_ocean,
                                 double salinity, double L) {
  const double shelfbaseelev = -(rho_i / rho_w) * H; // depth of the shelf base
  const double T_f = 273.15 + (0.0939 - 0.057 * salinity + 7.64e-4 * shelfbaseelev);
  const double ocean_heat_flux = melt_factor * rho_w * c_p_ocean * gamma_T * (T_ocean - T_f);
  return ocean_heat_flux / L;
}

// average_water_column_pressure from util/pism_utilities.cc (vertically
// averaged over the ice column at a margin).
double water_column_pressure(double H, double bed, double z_s, double rho_i, double rho_w,
                             double g) {
  const double ice_bottom = std::max(bed, z_s - rho_i / rho_w * H);
  const double h_w = std::max(z_s - ice_bottom, 0.0);
  if (H > 0.0) {
    return 0.5 * rho_w * g * pow(h_w, 2.0) / H;
  }
  return 0.0;
}

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
         "Configuration parameters affecting the atmosphere/ocean/surface "
         "forcing boundaries (atmosphere.*, ocean.*, surface.*, constants.*, "
         "grid.*)");
  header(out, "columns", "key,value");

  const std::set<std::string> prefixes = {
      "atmosphere.", "ocean.", "surface.", "constants.", "grid."};

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

void write_netcdf(const std::shared_ptr<const Grid> &grid, const std::string &path,
                  const std::vector<const array::Array *> &vecs) {
  auto writer = std::make_shared<SynchronousOutputWriter>(grid->com, *grid->ctx()->config());
  writer->initialize({}, true);
  OutputFile file(writer, path);

  auto time = grid->ctx()->time();
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
    std::shared_ptr<Context> ctx = context_from_options(com, "surface_ocean_test");
    auto config = ctx->config();
    auto sys = ctx->unit_system();

    set_config_from_options(*config);
    config->resolve_filenames();

    std::string out_dir = "stage1/dumps/surface_ocean";
    options::String output_dir("-dumps_dir", "Output directory for the dumps", out_dir);
    out_dir = output_dir;
    std::string mkdir_cmd = "mkdir -p " + out_dir;
    if (std::system(mkdir_cmd.c_str()) != 0) {
      throw RuntimeError(PISM_ERROR_LOCATION, "Cannot create the output directory");
    }

    const unsigned int Mx = config->get_number("grid.Mx");
    const unsigned int My = config->get_number("grid.My");

    // ------------------------------------------------------------------
    // Config values entering the kernels (recorded in parameters.csv too)
    // ------------------------------------------------------------------
    const double rho_i   = config->get_number("constants.ice.density");
    const double rho_w   = config->get_number("constants.sea_water.density");
    const double g       = config->get_number("constants.standard_gravity");
    const double T0      = config->get_number("constants.fresh_water.melting_point_temperature");
    const double beta_CC = config->get_number("constants.ice.beta_Clausius_Clapeyron");
    const double L       = config->get_number("constants.fresh_water.latent_heat_of_fusion");
    const double melt_factor = config->get_number("ocean.pik_melt_factor");
    const double c_p_ocean   = 3974.0; // J/(K kg), hard-coded in ConstantPIK.cc
    const double gamma_T     = 1e-4;   // m/s, hard-coded in ConstantPIK.cc
    const double salinity    = 35.0;   // g/kg, hard-coded in ConstantPIK.cc
    const double T_ocean     = units::convert(sys, -1.7, "degree_Celsius", "kelvin");

    // ------------------------------------------------------------------
    // Grid and geometry
    // ------------------------------------------------------------------
    const double Lx = 900.0e3, Ly = 900.0e3;   // domain half-widths
    const double L_for = 750.0e3;              // test-F ice extent
    auto grid = Grid::Shallow(ctx, Lx, Ly, 0.0, 0.0, Mx, My,
                              grid::CELL_CENTER, grid::NOT_PERIODIC);

    Geometry geometry(grid);
    geometry.sea_level_elevation.set(0.0);

    {
      array::AccessScope list{&geometry.ice_thickness, &geometry.bed_elevation,
                              &geometry.latitude, &geometry.longitude};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        const double y = grid->y(j);

        // latitude: linear in y, -80 (south edge) .. -60 (north edge)
        geometry.latitude(i, j) = -70.0 + 10.0 * (y / Ly);
        geometry.longitude(i, j) = 0.0;
        // ocean cavity/trough: deepens from -1500 m to -3000 m at r = 450 km
        geometry.bed_elevation(i, j) =
            -1500.0 - 1500.0 * std::exp(-std::pow((r - 450.0e3) / 200.0e3, 2.0));
        // test-F exact dome thickness (radial, extent 750 km)
        geometry.ice_thickness(i, j) =
            (r > L_for - 1.0) ? 0.0 : exactFG(0.0, r, grid->z(), 0.0).H;
      }
      geometry.ice_thickness.update_ghosts();
      geometry.bed_elevation.update_ghosts();
      geometry.latitude.update_ghosts();
      geometry.longitude.update_ghosts();
    }
    geometry.ensure_consistency(config->get_number("geometry.ice_free_thickness_standard"));

    // ------------------------------------------------------------------
    // (a) PIK atmosphere
    // ------------------------------------------------------------------
    auto atmos = std::make_shared<InstrumentedAtmospherePIK>(grid);

    // Prescribed constant-in-time precipitation field, in kg m^-2 s^-1:
    //   P(r) = (200 + 100 cos(pi r/Lx)) kg m^-2 year^-1
    {
      const double P_ref = units::convert(sys, 100.0, "kg m^-2 year^-1", "kg m^-2 second^-1");
      array::AccessScope list{&atmos->precip_field()};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        atmos->precip_field()(i, j) = P_ref * (2.0 + std::cos(M_PI * r / Lx));
      }
    }

    atmos->init(geometry);
    atmos->update(geometry, 0.0, units::convert(sys, 1.0, "year", "second"));

    // ------------------------------------------------------------------
    // (b) PIK ocean + Constant ocean (cross-check)
    // ------------------------------------------------------------------
    ocean::Inputs ocean_inputs;
    ocean_inputs.geometry = &geometry;

    auto ocean_pik = std::make_shared<ocean::PIK>(grid);
    ocean_pik->init(geometry);
    ocean_pik->update(ocean_inputs, 0.0, units::convert(sys, 1.0, "year", "second"));

    auto ocean_const = std::make_shared<ocean::Constant>(grid);
    ocean_const->init(geometry);
    ocean_const->update(ocean_inputs, 0.0, units::convert(sys, 1.0, "year", "second"));

    // ------------------------------------------------------------------
    // (c) PIK surface, "driven by" the PIK atmosphere: SMB = precipitation -
    //     a deterministic ablation blob (surface::PIK itself ignores the
    //     atmosphere pointer and holds the SMB constant, as in real runs
    //     where it is read from a file).
    // ------------------------------------------------------------------
    auto surf = std::make_shared<InstrumentedSurfacePIK>(grid, atmos);

    {
      const double abl_peak =
          units::convert(sys, 300.0, "kg m^-2 year^-1", "kg m^-2 second^-1");
      const double sigma = 250.0e3;
      array::AccessScope list{&surf->smb_field(), &atmos->precip_field()};
      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double x = grid->x(i), y = grid->y(j);
        // ablation blob centered on the warm-north part of the floating shelf
        // (x = 0, y = +450 km, the lat = -60 margin): peak 300 kg m^-2 yr^-1,
        // sigma 250 km, so a band of north-shelf cells has negative SMB
        // (melt > 0) while the rest of the ice stays in accumulation.
        const double d2 = ((x * x + (y - 0.5 * Ly) * (y - 0.5 * Ly)) / (sigma * sigma));
        surf->smb_field()(i, j) = atmos->precip_field()(i, j) - abl_peak * std::exp(-d2);
      }
    }

    surf->init(geometry);
    surf->update(geometry, 0.0, units::convert(sys, 1.0, "year", "second"));

    // ------------------------------------------------------------------
    // Dumps
    // ------------------------------------------------------------------

    {
      std::ofstream meta = open_output(out_dir + "/meta.txt");
      header(meta, "pism_revision", pism::revision);
      header(meta, "pism_config_file", pism::config_file);
      header(meta, "component",
             "atmosphere PIK + ocean PIK (Beckmann-Goosse) + ocean Constant + surface PIK");
      header(meta, "description",
             "PIK atmosphere (martin parameterization) air_temperature and "
             "precipitation maps + pointwise kernels; PIK ocean (Beckmann-Goosse) "
             "shelf_base_temperature / shelf_base_mass_flux / "
             "average_water_column_pressure maps + pointwise depth kernels + "
             "Constant ocean cross-check; PIK surface accumulation / melt / "
             "runoff / mass_flux / temperature maps (SMB = atmosphere "
             "precipitation minus a deterministic ablation blob)");
      header(meta, "geometry",
             "test-F exact dome (extent 750 km) on a Gaussian ocean trough "
             "bed(r) = -1500 - 1500 exp(-((r-450 km)/200 km)^2) m, sea level 0; "
             "latitude = -70 + 10 y/Ly deg (spanning -80..-60)");
      header(meta, "Mx", std::to_string(Mx));
      header(meta, "My", std::to_string(My));
      header(meta, "Lx", std::to_string(Lx));
      header(meta, "dx", std::to_string(grid->dx()));
      header(meta, "time", "t = 0, dt = 1 year (all update_impl ignore t/dt)");
      meta << "# command:";
      for (int i = 0; i < argc; ++i) {
        meta << " " << argv[i];
      }
      meta << "\n";
      meta.close();
    }

    dump_parameters(*config, out_dir);

    // --- inputs.nc: geometry + atmosphere + both oceans + surface ---
    // The Constant ocean outputs are copied to uniquely-named fields (both
    // ocean models use the same internal variable names).
    array::Scalar1 shelfbtemp_const(grid, "shelfbtemp_constant");
    shelfbtemp_const.copy_from(ocean_const->shelf_base_temperature());
    shelfbtemp_const.metadata(0).units("kelvin");

    array::Scalar1 shelfbmassflux_const(grid, "shelfbmassflux_constant");
    shelfbmassflux_const.copy_from(ocean_const->shelf_base_mass_flux());
    shelfbmassflux_const.metadata(0).units("kg m^-2 second^-1");

    array::Scalar1 awcp_const(grid, "average_water_column_pressure_constant");
    awcp_const.copy_from(ocean_const->average_water_column_pressure());
    awcp_const.metadata(0).units("Pa");

    {
      const std::vector<const array::Array *> vecs = {
          &geometry.latitude,              &geometry.longitude,
          &geometry.bed_elevation,         &geometry.sea_level_elevation,
          &geometry.ice_thickness,         &geometry.ice_surface_elevation,
          &geometry.cell_type,             &geometry.cell_grounded_fraction,
          &atmos->precipitation(),         &atmos->air_temperature(),
          &atmos->mean_summer_temp(),      &ocean_pik->shelf_base_temperature(),
          &ocean_pik->shelf_base_mass_flux(), &ocean_pik->average_water_column_pressure(),
          &shelfbtemp_const,               &shelfbmassflux_const,
          &awcp_const,                     &surf->mass_flux(),
          &surf->temperature(),            &surf->accumulation(),
          &surf->melt(),                   &surf->runoff()};
      write_netcdf(grid, out_dir + "/inputs.nc", vecs);
    }

    // --- columns.csv: per-cell map with analytic references ---
    {
      std::ofstream out = open_output(out_dir + "/columns.csv");
      header(out, "description",
             "Per-cell trace of the atmosphere/ocean/surface forcing maps: "
             "geometry, the PIK atmosphere outputs, the PIK + Constant ocean "
             "outputs, the PIK surface outputs, and the analytic reference "
             "values (martin T_ma, pressure melting temperature, Beckmann-Goosse "
             "mass flux, water column pressure) the component outputs must match");
      header(out, "columns",
             "i,j,x,y,r,lat,lon,H,surface,bed,mask,sea_level,"
             "precip_atm,air_temp_atm,T_ms_atm,"
             "shelfbtemp,shelfbmassflux,awcp,"
             "shelfbtemp_const,shelfbmassflux_const,awcp_const,"
             "smb_surface,surface_accumulation,surface_melt,surface_runoff,"
             "surface_mass_flux,surface_temperature,"
             "T_ma_analytic,T_shelf_analytic,mf_analytic,P_analytic");

      array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                              &geometry.bed_elevation, &geometry.cell_type,
                              &geometry.sea_level_elevation, &geometry.latitude,
                              &geometry.longitude, &atmos->precipitation(),
                              &atmos->air_temperature(), &atmos->mean_summer_temp(),
                              &ocean_pik->shelf_base_temperature(),
                              &ocean_pik->shelf_base_mass_flux(),
                              &ocean_pik->average_water_column_pressure(),
                              &ocean_const->shelf_base_temperature(),
                              &ocean_const->shelf_base_mass_flux(),
                              &ocean_const->average_water_column_pressure(),
                              &surf->mass_flux(), &surf->temperature(),
                              &surf->accumulation(), &surf->melt(), &surf->runoff()};

      for (auto p : grid->points()) {
        const int i = p.i(), j = p.j();
        const double r = std::max(grid::radius(*grid, i, j), 1.0);
        const double H = geometry.ice_thickness(i, j);
        const double usurf = geometry.ice_surface_elevation(i, j);
        const double lat = geometry.latitude(i, j);

        out << i << "," << j << "," << grid->x(i) << "," << grid->y(j) << "," << r << ","
            << lat << "," << geometry.longitude(i, j) << "," << H << "," << usurf << ","
            << geometry.bed_elevation(i, j) << "," << (int)geometry.cell_type(i, j) << ","
            << geometry.sea_level_elevation(i, j) << "," << atmos->precipitation()(i, j) << ","
            << atmos->air_temperature()(i, j) << "," << atmos->mean_summer_temp()(i, j) << ","
            << ocean_pik->shelf_base_temperature()(i, j) << ","
            << ocean_pik->shelf_base_mass_flux()(i, j) << ","
            << ocean_pik->average_water_column_pressure()(i, j) << ","
            << ocean_const->shelf_base_temperature()(i, j) << ","
            << ocean_const->shelf_base_mass_flux()(i, j) << ","
            << ocean_const->average_water_column_pressure()(i, j) << ","
            << surf->mass_flux()(i, j) << "," << surf->accumulation()(i, j) << ","
            << surf->melt()(i, j) << "," << surf->runoff()(i, j) << ","
            << surf->mass_flux()(i, j) << "," << surf->temperature()(i, j) << ","
            << martin2011_mean_annual(usurf, lat) << ","
            << shelf_base_temperature(H, T0, beta_CC, rho_i, g) << ","
            << beckmann_goosse_mass_flux(H, melt_factor, rho_i, rho_w, c_p_ocean, gamma_T,
                                         T_ocean, salinity, L)
            << "," << water_column_pressure(H, geometry.bed_elevation(i, j),
                                            geometry.sea_level_elevation(i, j), rho_i, rho_w, g)
            << "\n";
      }
      out.close();
    }

    // --- pointwise kernels (pure .esm test targets) ---

    // atmosphere T_ma vs (latitude, elevation); martin: T_ms = T_ma
    {
      std::ofstream out = open_output(out_dir + "/atmosphere_temperature.csv");
      header(out, "description",
             "atmosphere::PIK (martin parameterization) mean annual / mean "
             "summer temperature: T_ma = 303.15 + 0.68775*lat - 0.0075*usurf "
             "(K, lat in degrees north), T_ms = T_ma (no seasonal cycle)");
      header(out, "columns", "latitude_deg,surface_elevation_m,T_ma_K,T_ms_K");
      for (double lat : {-60.0, -65.0, -70.0, -75.0, -80.0}) {
        for (double h : {0.0, 500.0, 1000.0, 1500.0, 2000.0, 2500.0, 3000.0}) {
          const double T_ma = martin2011_mean_annual(h, lat);
          out << lat << "," << h << "," << T_ma << "," << T_ma << "\n";
        }
      }
      out.close();
    }

    // atmosphere precipitation field (time-independent): r -> P
    {
      std::ofstream out = open_output(out_dir + "/atmosphere_precipitation.csv");
      header(out, "description",
             "atmosphere::PIK precipitation field prescribed by the driver: "
             "P(r) = (200 + 100 cos(pi r/Lx)) kg m^-2 year^-1 (converted to "
             "kg m^-2 s^-1); time-independent (precip_time_series is constant)");
      header(out, "columns", "r_m,P_kg_m-2_s-1,P_kg_m-2_yr-1");
      for (double r = 0.0; r <= Lx + 1.0; r += Lx / 20.0) {
        const double P_s = units::convert(sys, 100.0 * (2.0 + std::cos(M_PI * r / Lx)),
                                          "kg m^-2 year^-1", "kg m^-2 second^-1");
        out << r << "," << P_s << "," << 100.0 * (2.0 + std::cos(M_PI * r / Lx)) << "\n";
      }
      out.close();
    }

    // atmosphere temperature time series over one year at three sample points
    // (uses the component's own init_timeseries / temp_time_series path).
    {
      const double sec_per_year = units::convert(sys, 1.0, "year", "second");
      std::vector<double> ts;
      for (int k = 0; k <= 12; ++k) {
        ts.push_back(k * sec_per_year / 12.0);
      }

      std::vector<double> T_center, T_north, T_south, P_center;
      atmos->init_timeseries(ts);
      atmos->begin_pointwise_access();
      atmos->temp_time_series(15, 15, T_center);
      atmos->temp_time_series(15, 30, T_north);
      atmos->temp_time_series(15, 0, T_south);
      atmos->precip_time_series(15, 15, P_center);
      atmos->end_pointwise_access();

      std::ofstream out = open_output(out_dir + "/atmosphere_temperature_timeseries.csv");
      header(out, "description",
             "Yearly-cycle temperature time series from the component's "
             "temp_time_series (T(t) = T_ma + (T_ms - T_ma) cos(2 pi (year_frac "
             "- summerday_frac))); for the martin parameterization T_ms = T_ma, "
             "so the cycle is flat (zero amplitude). Sample points: (15,15) "
             "near-dome, (15,30) north edge (lat=-60), (15,0) south edge "
             "(lat=-80).");
      header(out, "columns", "k,time_s,T_center_K,T_north_K,T_south_K,P_center_kg_m-2_s-1");
      for (size_t k = 0; k < ts.size(); ++k) {
        out << k << "," << ts[k] << "," << T_center[k] << "," << T_north[k] << ","
            << T_south[k] << "," << P_center[k] << "\n";
      }
      out.close();
    }

    // shelf base temperature vs depth (pressure melting point)
    {
      std::ofstream out = open_output(out_dir + "/shelf_base_temperature.csv");
      header(out, "description",
             "ocean::PIK::melting_point_temperature (same in ocean::Constant): "
             "T_shelf = T0 - beta_CC * rho_i * g * H (linear pressure melting "
             "point of fresh ice at the shelf base, decreasing with H)");
      header(out, "columns", "H_m,pressure_Pa,T_shelf_K");
      for (double H = 0.0; H <= 3000.0 + 1.0; H += 100.0) {
        out << H << "," << rho_i * g * H << "," << shelf_base_temperature(H, T0, beta_CC, rho_i, g)
            << "\n";
      }
      out.close();
    }

    // Beckmann-Goosse mass flux vs depth
    {
      std::ofstream out = open_output(out_dir + "/beckmann_goosse_mass_flux.csv");
      header(out, "description",
             "ocean::PIK::mass_flux: z_b = -(rho_i/rho_w) H; "
             "T_f = 273.15 + 0.0939 - 0.057*S + 7.64e-4*z_b; "
             "Q = melt_factor*rho_w*c_p_ocean*gamma_T*(T_ocean - T_f); "
             "mass_flux = Q/L. Positive = melting (dH_BMB = -dt*mf/rho_i in the "
             "geometry update). Linear in H; always melting with the default "
             "T_ocean = -1.7 C.");
      header(out, "columns", "H_m,z_b_m,T_f_K,T_ocean_K,ocean_heat_flux_W_m2,mass_flux_kg_m-2_s-1");
      for (double H = 0.0; H <= 3000.0 + 1.0; H += 100.0) {
        const double z_b = -(rho_i / rho_w) * H;
        const double T_f = 273.15 + (0.0939 - 0.057 * salinity + 7.64e-4 * z_b);
        const double Q = melt_factor * rho_w * c_p_ocean * gamma_T * (T_ocean - T_f);
        out << H << "," << z_b << "," << T_f << "," << T_ocean << "," << Q << "," << Q / L
            << "\n";
      }
      out.close();
    }

    // average water column pressure vs (H, bed)
    {
      std::ofstream out = open_output(out_dir + "/water_column_pressure.csv");
      header(out, "description",
             "average_water_column_pressure(H, bed, z_s=0): "
             "ice_bottom = max(bed, -rho_i/rho_w H); h_w = max(z_s - ice_bottom, 0); "
             "P = 0.5 rho_w g h_w^2 / H for H > 0, else 0 (depth-averaged water "
             "pressure over the ice column at a margin)");
      header(out, "columns", "H_m,bed_m,z_s_m,ice_bottom_m,h_w_m,P_Pa");
      for (double H : {0.0, 500.0, 1000.0, 2000.0, 3000.0}) {
        for (double bed : {0.0, -500.0, -1500.0, -3000.0}) {
          const double ice_bottom = std::max(bed, -rho_i / rho_w * H);
          const double h_w = std::max(0.0 - ice_bottom, 0.0);
          const double P = water_column_pressure(H, bed, 0.0, rho_i, rho_w, g);
          out << H << "," << bed << ",0," << ice_bottom << "," << h_w << "," << P << "\n";
        }
      }
      out.close();
    }

    // surface SMB partition (dummy_accumulation/dummy_melt/dummy_runoff)
    {
      std::ofstream out = open_output(out_dir + "/surface_partition.csv");
      header(out, "description",
             "surface::PIK SMB partition (SurfaceModel::dummy_*): "
             "accumulation = max(smb, 0), melt = runoff = max(-smb, 0), "
             "so mass_flux = accumulation - runoff always");
      header(out, "columns", "smb_kg_m-2_s-1,accumulation,melt,runoff,mass_flux_identity");
      for (double smb : {-5.0e-5, -3.0e-5, -1.0e-5, -1.0e-8, 0.0, 1.0e-8, 1.0e-5, 3.0e-5, 5.0e-5}) {
        const double acc = std::max(smb, 0.0);
        const double run = std::max(-smb, 0.0);
        out << smb << "," << acc << "," << run << "," << run << "," << acc - run << "\n";
      }
      out.close();
    }

    // ------------------------------------------------------------------
    // Validation
    // ------------------------------------------------------------------
    double max_err_T_ma = 0.0, max_err_T_shelf = 0.0, max_err_mf = 0.0, max_err_P = 0.0;
    double max_err_smb_identity = 0.0, max_err_surf_T = 0.0;
    double max_err_T_shelf_const = 0.0, max_err_mf_const = 0.0, max_err_P_const = 0.0;
    double min_air_T = std::numeric_limits<double>::infinity(),
           max_air_T = -std::numeric_limits<double>::infinity();
    double min_mf = std::numeric_limits<double>::infinity(),
           max_mf = -std::numeric_limits<double>::infinity();
    double max_awcp = -std::numeric_limits<double>::infinity(),
           min_awcp = std::numeric_limits<double>::infinity();
    double dome_T_ma = 0.0, dome_lat = 0.0, dome_usurf = 0.0, min_r = std::numeric_limits<double>::infinity();
    long n_grounded = 0, n_floating = 0, n_ocean = 0, n_melt_cells = 0;

    // Constant ocean: mass flux is uniform = melt_rate * ice_density (the
    // field itself is compared inside the AccessScope below).
    const double mf_const_value =
        config->get_number("ocean.constant.melt_rate", "m second-1") * rho_i;

    array::AccessScope list{&geometry.ice_thickness, &geometry.ice_surface_elevation,
                            &geometry.bed_elevation, &geometry.cell_type,
                            &geometry.sea_level_elevation, &geometry.latitude,
                            &atmos->precipitation(), &atmos->air_temperature(),
                            &ocean_pik->shelf_base_temperature(),
                            &ocean_pik->shelf_base_mass_flux(),
                            &ocean_pik->average_water_column_pressure(),
                            &ocean_const->shelf_base_temperature(),
                            &ocean_const->shelf_base_mass_flux(),
                            &ocean_const->average_water_column_pressure(),
                            &surf->mass_flux(), &surf->temperature(),
                            &surf->accumulation(), &surf->melt(), &surf->runoff()};

    for (auto p : grid->points()) {
      const int i = p.i(), j = p.j();
      const double r = std::max(grid::radius(*grid, i, j), 1.0);
      const double H = geometry.ice_thickness(i, j);
      const double usurf = geometry.ice_surface_elevation(i, j);
      const double lat = geometry.latitude(i, j);
      const double bed = geometry.bed_elevation(i, j);

      max_err_T_ma =
          std::max(max_err_T_ma, std::fabs(atmos->air_temperature()(i, j) -
                                           martin2011_mean_annual(usurf, lat)));
      max_err_T_shelf =
          std::max(max_err_T_shelf,
                   std::fabs(ocean_pik->shelf_base_temperature()(i, j) -
                             shelf_base_temperature(H, T0, beta_CC, rho_i, g)));
      max_err_mf = std::max(max_err_mf, std::fabs(ocean_pik->shelf_base_mass_flux()(i, j) -
                                                  beckmann_goosse_mass_flux(
                                                      H, melt_factor, rho_i, rho_w, c_p_ocean,
                                                      gamma_T, T_ocean, salinity, L)));
      max_err_P = std::max(max_err_P,
                           std::fabs(ocean_pik->average_water_column_pressure()(i, j) -
                                     water_column_pressure(H, bed, 0.0, rho_i, rho_w, g)));
      max_err_T_shelf_const =
          std::max(max_err_T_shelf_const,
                   std::fabs(ocean_const->shelf_base_temperature()(i, j) -
                             shelf_base_temperature(H, T0, beta_CC, rho_i, g)));
      max_err_mf_const =
          std::max(max_err_mf_const,
                   std::fabs(ocean_const->shelf_base_mass_flux()(i, j) - mf_const_value));
      max_err_P_const =
          std::max(max_err_P_const,
                   std::fabs(ocean_const->average_water_column_pressure()(i, j) -
                             water_column_pressure(H, bed, 0.0, rho_i, rho_w, g)));
      max_err_smb_identity = std::max(
          max_err_smb_identity,
          std::fabs(surf->mass_flux()(i, j) -
                    (surf->accumulation()(i, j) - surf->runoff()(i, j))));
      max_err_surf_T = std::max(max_err_surf_T, std::fabs(surf->temperature()(i, j) -
                                                          atmos->air_temperature()(i, j)));

      min_air_T = std::min(min_air_T, atmos->air_temperature()(i, j));
      max_air_T = std::max(max_air_T, atmos->air_temperature()(i, j));
      min_mf = std::min(min_mf, ocean_pik->shelf_base_mass_flux()(i, j));
      max_mf = std::max(max_mf, ocean_pik->shelf_base_mass_flux()(i, j));
      min_awcp = std::min(min_awcp, ocean_pik->average_water_column_pressure()(i, j));
      max_awcp = std::max(max_awcp, ocean_pik->average_water_column_pressure()(i, j));

      if (geometry.cell_type.grounded(i, j)) {
        ++n_grounded;
      } else if (geometry.cell_type.floating_ice(i, j)) {
        ++n_floating;
      } else if (geometry.cell_type.ice_free_ocean(i, j)) {
        ++n_ocean;
      }
      if (surf->melt()(i, j) > 0.0) {
        ++n_melt_cells;
      }

      if (r < min_r) { // dome-center cell (minimum radius)
        min_r = r;
        dome_T_ma = atmos->air_temperature()(i, j);
        dome_lat = lat;
        dome_usurf = usurf;
      }
    }

    // Beckmann-Goosse analytic checks at sample depths
    const double mf_H0 = beckmann_goosse_mass_flux(0.0, melt_factor, rho_i, rho_w, c_p_ocean,
                                                   gamma_T, T_ocean, salinity, L);
    const double mf_H2000 = beckmann_goosse_mass_flux(2000.0, melt_factor, rho_i, rho_w,
                                                      c_p_ocean, gamma_T, T_ocean, salinity, L);
    const double T_shelf_H0 = shelf_base_temperature(0.0, T0, beta_CC, rho_i, g);
    const double T_shelf_H2000 = shelf_base_temperature(2000.0, T0, beta_CC, rho_i, g);
    const double P_sample = water_column_pressure(2000.0, -3000.0, 0.0, rho_i, rho_w, g);
    const double P_grounded = water_column_pressure(2500.0, -1509.5, 0.0, rho_i, rho_w, g);

    const double lat_lo = -80.0, lat_hi = -60.0;
    std::cout
        << "SURFACE/OCEAN (PIK atmosphere + PIK ocean + Constant ocean + PIK surface) "
           "instrumentation validation:\n"
        << "  grid: " << Mx << "x" << My << ", dx = " << grid->dx()
        << " m; cell types: grounded = " << n_grounded << ", floating = " << n_floating
        << ", ice-free ocean = " << n_ocean << "\n"
        << "  domain air_temperature range: [" << min_air_T << ", " << max_air_T << "] K\n"
        << "  domain PIK shelf_base_mass_flux range: [" << min_mf << ", " << max_mf
        << "] kg m^-2 s^-1 (positive = melting; always positive under the shelf)\n"
        << "  domain average_water_column_pressure range: [" << min_awcp << ", " << max_awcp
        << "] Pa\n"
        << "  surface melt > 0 cells: " << n_melt_cells << "\n"
        << "  ATMOSPHERE PIK (martin):\n"
        << "    T_ma(lat=" << lat_lo << ", h=0) analytic = "
        << martin2011_mean_annual(0.0, lat_lo) << " K\n"
        << "    T_ma(lat=" << lat_hi << ", h=0) analytic = "
        << martin2011_mean_annual(0.0, lat_hi) << " K\n"
        << "    T_ma(lat=-70, h=1500) analytic = " << martin2011_mean_annual(1500.0, -70.0)
        << " K\n"
        << "    near-dome cell: lat = " << dome_lat << ", usurf = " << dome_usurf
        << ", air_temp = " << dome_T_ma << " K\n"
        << "    max |air_temp - T_ma(usurf, lat)| over the map = " << max_err_T_ma << " K\n"
        << "  OCEAN PIK (Beckmann-Goosse):\n"
        << "    T_shelf(H=0) analytic = " << T_shelf_H0 << " K (T0)\n"
        << "    T_shelf(H=2000) analytic = " << T_shelf_H2000
        << " K (pressure melting point, linear in H)\n"
        << "    max |shelfbtemp - analytic| over the map = " << max_err_T_shelf << " K\n"
        << "    mass_flux(H=0) analytic = " << mf_H0 << " kg m^-2 s^-1\n"
        << "    mass_flux(H=2000) analytic = " << mf_H2000
        << " kg m^-2 s^-1 = " << mf_H2000 / rho_i * units::convert(sys, 1.0, "year", "second")
        << " m/yr ice-equivalent\n"
        << "    max |shelfbmassflux - analytic| over the map = " << max_err_mf
        << " kg m^-2 s^-1\n"
        << "    water column pressure (H=2000, bed=-3000) analytic = " << P_sample << " Pa\n"
        << "    water column pressure (H=2500, bed=-1509.5, grounded) analytic = " << P_grounded
        << " Pa\n"
        << "    max |awcp - analytic| over the map = " << max_err_P << " Pa\n"
        << "  OCEAN CONSTANT (cross-check):\n"
        << "    shelf_base_temperature = same pressure melting point (max err "
        << max_err_T_shelf_const << " K)\n"
        << "    shelf_base_mass_flux = " << mf_const_value
        << " kg m^-2 s^-1 everywhere (max err " << max_err_mf_const << ")\n"
        << "    awcp same as PIK (max err " << max_err_P_const << " Pa)\n"
        << "  SURFACE PIK:\n"
        << "    max |mass_flux - (accumulation - runoff)| over the map = "
        << max_err_smb_identity << " kg m^-2 s^-1 (must be 0)\n"
        << "    max |surface_temperature - atmosphere air_temperature| over the map = "
        << max_err_surf_T << " K (both martin formula, must be 0)\n";

    std::cout << "Surface/ocean instrumentation complete. Dumps in " << out_dir << "\n";

  } catch (...) {
    handle_fatal_errors(com);
    return 1;
  }

  return 0;
}
