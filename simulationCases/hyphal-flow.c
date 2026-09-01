/**
# hyphal-flow.c

Canonical single-branch hypha simulation. The outer phase is an
incompressible finite-strain neo-Hookean Kelvin--Voigt-type solid, while the
carrier liquid and drop are finite-relaxation Oldroyd-B fluids. Separate VoF
fields prevent the two interfaces from becoming one tracer; the overlap
diagnostic still has to demonstrate that the carrier film remains resolved.

## Author
Vatsal Sanjay

## Version
- August 11, 2024: allow stress relaxation in both drop and
  cytoplasm.
*/

#include "axi.h"
#include "params.h"
#include "navier-stokes/centered.h"
#include "log-conform-viscoelastic.h"
#define FILTERED
#include "three-phase-rheology.h"
#include "tension.h"
#include "reduced-three-phase-nonCoalescing.h"

/**
## Numerical Tolerances
*/
#define fErr (1e-3) // error tolerance in VOF
#define KErr (1e-4) // error tolerance in KAPPA
#define VelErr (1e-2) // error tolerances in velocity
#define AErr (1e-3) // error tolerance in Conformation tensor
#define MINlevel 4 // minimum level
int MAXlevel;
int initial_refine_level;
double tmax;
double snapshot_interval;
double log_interval;
double initial_refine_band;
double initial_volume_tolerance;

/**
## Material Parameters

The Newtonian `Oh_*` contribution is the solvent viscosity. For the two
liquids, `Ec_*` is the elastic modulus and `De_*` the relaxation time, so the
implied Oldroyd-B polymer viscosity is `Ec_* De_*`. The solid uses the same
finite-strain conformation stress in the practical no-relaxation limit; its
Newtonian contribution is the Kelvin--Voigt dashpot.
*/
double Oh_drop, rho_drop, Ec_drop, De_drop;
double Oh_solid, rho_solid, Ec_solid, De_solid;
double Oh_liquid, rho_liquid, Ec_liquid, De_liquid;

/**
## Geometry Helper
*/
#define gap(x,y,radius,width,x0,c0) (y - ((c0+1.) + 0.5*(radius - (c0+1.))*(1. + tanh(sq(x-x0)/width))))

double channel_radius;
double Bond;
double Ldomain;

static void simulation_abort (int status)
{
#if _MPI
  MPI_Abort (MPI_COMM_WORLD, status);
#else
  exit (status);
#endif
}

static int validate_phase (const char * name, double viscosity,
                           double modulus, double relaxation_time)
{
  if (!isfinite (viscosity) || !isfinite (modulus) ||
      !isfinite (relaxation_time) || viscosity <= 0. || modulus < 0. ||
      relaxation_time < 0.) {
    fprintf (stderr, "ERROR: %s viscosity must be positive and Ec/De non-negative\n",
             name);
    return -1;
  }
  if ((modulus > 0.) != (relaxation_time > 0.)) {
    fprintf (stderr,
             "ERROR: %s needs both Ec > 0 and De > 0, or Ec = De = 0\n",
             name);
    return -1;
  }
  return 0;
}

static bool output_due (double current_time, double * next_time,
                        double interval)
{
  if (current_time + 1e-12 < *next_time)
    return false;
  do
    *next_time += interval;
  while (*next_time <= current_time + 1e-12);
  return true;
}

/**
## validate_initial_geometry()

Fail before the first timestep when the freshly constructed VoF geometry is
not the intended resolved axisymmetric state. The direct drop metric volume
has the analytic value

$$
\int_{-1.5}^{1.5}\int_0^{\sqrt{1-(x/1.5)^2}} y\,dy\,dx = 1.
$$

This deliberately uses `f1*dv()` rather than normalised phase weights, so
drop--solid overlap cannot hide a bad initial condition.
*/
static void validate_initial_geometry (void)
{
  const double f_eps = 1e-6;
  double drop_volume = 0., overlap = 0.;
  double drop_xmin = HUGE, drop_xmax = -HUGE, drop_ymax = -HUGE;
  long drop_mixed = 0, solid_mixed = 0;
  int drop_min_level = initial_refine_level + 1;
  int solid_min_level = initial_refine_level + 1;

  foreach (reduction(+:drop_volume) reduction(max:overlap)
           reduction(min:drop_xmin) reduction(max:drop_xmax)
           reduction(max:drop_ymax) reduction(+:drop_mixed)
           reduction(+:solid_mixed) reduction(min:drop_min_level)
           reduction(min:solid_min_level)) {
    drop_volume += f1[]*dv();
    overlap = max (overlap, min (f1[], f2[]));
    if (f1[] > f_eps && f1[] < 1. - f_eps) {
      drop_mixed++;
      drop_min_level = min (drop_min_level, level);
      drop_xmin = min (drop_xmin, x);
      drop_xmax = max (drop_xmax, x);
      drop_ymax = max (drop_ymax, y);
    }
    if (f2[] > f_eps && f2[] < 1. - f_eps) {
      solid_mixed++;
      solid_min_level = min (solid_min_level, level);
    }
  }

  const double geometry_tolerance = 2.*L0/(1 << initial_refine_level);
  const bool valid =
    fabs (drop_volume - 1.) <= initial_volume_tolerance &&
    overlap <= 1e-6 &&
    drop_mixed > 0 && solid_mixed > 0 &&
    drop_min_level == initial_refine_level &&
    solid_min_level == initial_refine_level &&
    fabs (drop_xmin + 1.5) <= geometry_tolerance &&
    fabs (drop_xmax - 1.5) <= geometry_tolerance &&
    fabs (drop_ymax - 1.) <= geometry_tolerance;

  if (pid() == 0)
    fprintf (ferr,
             "INIT_GEOMETRY valid=%d Vd=%.16g overlap=%.16g "
             "drop_mixed=%ld solid_mixed=%ld drop_min_level=%d "
             "solid_min_level=%d xmin=%.16g xmax=%.16g ymax=%.16g\n",
             valid, drop_volume, overlap, drop_mixed, solid_mixed,
             drop_min_level, solid_min_level,
             drop_xmin, drop_xmax, drop_ymax);
  if (!valid) {
    if (pid() == 0)
      fprintf (ferr,
               "ERROR: unresolved or inconsistent initial VoF geometry\n");
    simulation_abort (3);
  }
}

/**
## main()

Initialize properties and forcing, then enter the Basilisk event loop.
*/
int main (int argc, char const * argv[])
{
  // The command is a fixed literal; no parameter text reaches a shell.
  if (system ("mkdir -p intermediate"))
    return 2;

  params_init_from_argv (argc, argv);

  MAXlevel = param_int ("MAXlevel", 12);
  initial_refine_level = param_int ("initial_refine_level", MAXlevel);
  tmax = param_double ("tmax", 200.);
  snapshot_interval = param_double ("snapshot_interval", 0.1);
  log_interval = param_double ("log_interval", 0.01);
  initial_refine_band = param_double ("initial_refine_band", 1.25);
  initial_volume_tolerance = param_double ("initial_volume_tolerance", 0.01);

  Oh_drop = param_double ("Oh_drop", param_double ("Ohd", 1.));
  rho_drop = param_double ("rho_drop", param_double ("RhoR_dc", 1.2));
  Ec_drop = param_double ("Ec_drop", param_double ("Ec_d", 0.1));
  De_drop = param_double ("De_drop", param_double ("De_d", 1.));

  Oh_solid = param_double ("Oh_solid", param_double ("Ohf", 1.));
  rho_solid = param_double ("rho_solid", param_double ("RhoR_hc", 1.));
  Ec_solid = param_double ("Ec_solid", param_double ("Ec_h", 0.1));
  De_solid = param_double ("De_solid", param_double ("De_h", 1e30));

  Oh_liquid = param_double ("Oh_liquid",
                            param_double ("Oh_c", param_double ("Ohc", 0.01)));
  rho_liquid = param_double ("rho_liquid", 1.);
  Ec_liquid = param_double ("Ec_liquid", param_double ("Ec_c", 0.1));
  De_liquid = param_double ("De_liquid", param_double ("De_c", 1.));

  channel_radius = param_double ("channel_radius",
                                 param_double ("hr", param_double ("hf", 0.9)));
  Bond = param_double ("Bond", 1.);
  Ldomain = param_double ("Ldomain", 16.);

  if (MAXlevel < MINlevel || MAXlevel > 20 ||
      initial_refine_level < MINlevel || initial_refine_level > MAXlevel ||
      !isfinite (tmax) || tmax <= 0. ||
      !isfinite (snapshot_interval) || snapshot_interval <= 0. ||
      !isfinite (log_interval) || log_interval <= 0. ||
      !isfinite (initial_refine_band) || initial_refine_band <= 1. ||
      !isfinite (initial_volume_tolerance) ||
      initial_volume_tolerance <= 0. || initial_volume_tolerance >= 0.1 ||
      !isfinite (rho_drop) || !isfinite (rho_solid) ||
      !isfinite (rho_liquid) || rho_drop <= 0. || rho_solid <= 0. ||
      rho_liquid <= 0. || !isfinite (channel_radius) ||
      !isfinite (Ldomain) || !isfinite (Bond) || channel_radius <= 0. ||
      Ldomain < 8. ||
      validate_phase ("drop", Oh_drop, Ec_drop, De_drop) ||
      validate_phase ("outer solid", Oh_solid, Ec_solid, De_solid) ||
      validate_phase ("carrier liquid", Oh_liquid, Ec_liquid, De_liquid))
    return 2;

  if (Ec_solid > 0. && De_solid < 1e6*max (tmax, 1.))
    fprintf (stderr,
             "WARNING: De_solid is finite on the run timescale; the outer "
             "phase is relaxing viscoelastic rather than Kelvin--Voigt-type\n");

  fprintf (ferr,
           "level=%d tmax=%g | drop: Oh=%g Ec=%g De=%g mu_p=%g | "
           "liquid: Oh=%g Ec=%g De=%g mu_p=%g | solid: Oh=%g Ec=%g "
           "De=%g | radius=%g Bond=%g | Ldomain=%g initial_level=%d "
           "initial_band=%g\n",
           MAXlevel, tmax,
           Oh_drop, Ec_drop, De_drop, Ec_drop*De_drop,
           Oh_liquid, Ec_liquid, De_liquid, Ec_liquid*De_liquid,
           Oh_solid, Ec_solid, De_solid, channel_radius, Bond,
           Ldomain, initial_refine_level, initial_refine_band);

  L0 = Ldomain;
  X0 = -4.;
  Y0 = 0.;
  init_grid (1 << MINlevel);
  periodic (right);

  // Phase 1: drop (finite-relaxation Oldroyd-B).
  rho1 = rho_drop; mu1 = Oh_drop; G1 = Ec_drop; lambda1 = De_drop;

  // Phase 2: outer finite-strain Kelvin--Voigt-type solid.
  rho2 = rho_solid; mu2 = Oh_solid; G2 = Ec_solid; lambda2 = De_solid;

  // Phase 3: carrier liquid (finite-relaxation Oldroyd-B).
  rho3 = rho_liquid; mu3 = Oh_liquid; G3 = Ec_liquid; lambda3 = De_liquid;

  Bf1.x = Bond;
  Bf2.x = Bond;

  f1.sigma = 1.; // drop--carrier-liquid interfacial tension
  f2.sigma = 1.; // solid--carrier-liquid interfacial tension

  run ();
  return 0;
}

/**
## init()

Initialize interfaces unless a restart snapshot is available.
*/
event init(t = 0){
  if (!restore (file = "restart")) {
    double width = 2e0; // width of the tanh function
    // Historical wall parameter; the current formula gives a 0.05 minimum gap.
    double clearance = 0.20;
    double x0tanh = 0.0; // midpoint of the tanh function

    // Resolve both interfaces before sampling either implicit geometry.
    refine (y < initial_refine_band && level < initial_refine_level);

    fraction(f1, sq(1e0) - sq(y) - sq(x/1.5));
    fraction(f2, gap(x,y,channel_radius,width,x0tanh,clearance));
    validate_initial_geometry();
  }
}

/**
## adapt()

Adaptive mesh refinement driven by interfaces, curvature, velocity, and
conformation fields.
*/
event adapt(i++){
  scalar KAPPA1[], KAPPA2[], trA[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);
  foreach(){
    trA[] = (conform_p.x.x[] + conform_p.y.y[] + conform_qq[])/3.0;
  }

  adapt_wavelet ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, trA},
  (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, AErr},
  MAXlevel, MINlevel);

  unrefine(x > X0 + L0 - 1.);
  unrefine(y > 1e1);
}
/**
## stop_when_drop_exits()

Stop the run when the leading drop edge approaches the domain outlet.
*/
event stop_when_drop_exits (i++) {

  static double next_check = 0.;
  if (!output_due (t, &next_check, log_interval))
    return 0;

  scalar xpos[];
  coord ex = {1., 0.};
  coord z0 = {0., 0.};
  position (f1, xpos, ex, z0, add = false);

  stats sx = statsf (xpos);
  double xmax = sx.volume > 0. ? sx.max : -HUGE;

  // buffer
  double finest = L0/(1 << MAXlevel);
  double buffer = 2.*finest;

  double x_end = X0 + L0;

  if (pid() == 0)
    fprintf(ferr, "drop front xmax = %.6f, x_end = %.6f\n", xmax, x_end);

  if (xmax > x_end - buffer) {
    if (pid() == 0) {
      fprintf(ferr,
        "\n*** Drop leading edge reached end of domain ***\n"
        "xmax = %.6f, domain end = %.6f\n"
        "Stopping simulation at t = %.6g\n\n",
        xmax, x_end, t);
    }
    dump (file = "final");
    return 1;
  }
}

/**
## writingFiles()

Write periodic restart and snapshot files.
*/
event writingFiles (i++) {
  static double next_snapshot = 0.;
  if (!output_due (t, &next_snapshot, snapshot_interval))
    return 0;
  dump (file = "restart");
  char nameOut[80];
  snprintf (nameOut, sizeof nameOut, "intermediate/snapshot-%5.4f",
            next_snapshot - snapshot_interval);
  dump (file = nameOut);
}

/**
## stop_at_tmax()

Write a terminal snapshot and stop the otherwise open-ended event loop at the
configured physical time. This makes reduced runs genuinely bounded.
*/
event stop_at_tmax (i++; t <= HUGE) {
  if (t + 1e-12 < tmax)
    return 0;
  dump (file = "final");
  return 1;
}

/**
## logWriting()

Log kinetic energy and droplet center-of-mass velocity.
*/
event logWriting (i++) {
  static double next_log = 0.;
  if (!output_due (t, &next_log, log_interval))
    return 0;
  double ke = 0., vcm = 0., wt = 0., overlap = 0.;
  double volume_drop = 0., volume_solid = 0., volume_liquid = 0.;
  double stress_drop = 0., stress_solid = 0., stress_liquid = 0.;
  foreach (reduction(+:ke) reduction(+:vcm) reduction(+:wt)
           reduction(+:volume_drop) reduction(+:volume_solid)
           reduction(+:volume_liquid) reduction(max:overlap)
           reduction(max:stress_drop) reduction(max:stress_solid)
           reduction(max:stress_liquid)) {
    double volume = dv();
    double wd, ws, wl;
    phase_weights (f1[], f2[], &wd, &ws, &wl);
    double stress_norm = sqrt (sq(tau_p.x.x[]) + sq(tau_p.y.y[]) +
                               2.*sq(tau_p.x.y[]) + sq(tau_qq[]));
    ke += 0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[]))*volume;
    vcm += f1[]*u.x[]*volume;
    wt += f1[]*volume;
    volume_drop += wd*volume;
    volume_solid += ws*volume;
    volume_liquid += wl*volume;
    overlap = max (overlap, min (f1[], f2[]));
    stress_drop = max (stress_drop, wd*stress_norm);
    stress_solid = max (stress_solid, ws*stress_norm);
    stress_liquid = max (stress_liquid, wl*stress_norm);
  }
  if (wt > 0.0) vcm /= wt;
  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm overlap Vd Vs Vl Td Ts Tl\n");
      fp = fopen ("log", "w");
      if (!fp) {
        perror ("log");
        simulation_abort (4);
      }
      fprintf (fp, "i dt t ke vcm overlap Vd Vs Vl Td Ts Tl\n");
    } else {
      fp = fopen ("log", "a");
      if (!fp) {
        perror ("log");
        simulation_abort (4);
      }
    }
    fprintf (fp, "%d %g %g %g %5.4e %g %g %g %g %g %g %g\n",
             i, dt, t, ke, vcm, overlap,
             volume_drop, volume_solid, volume_liquid,
             stress_drop, stress_solid, stress_liquid);
    fclose(fp);
    fprintf (ferr, "%d %g %g %g %5.4e %g %g %g %g %g %g %g\n",
             i, dt, t, ke, vcm, overlap,
             volume_drop, volume_solid, volume_liquid,
             stress_drop, stress_solid, stress_liquid);
    if (overlap > 1e-3)
      fprintf (ferr,
               "WARNING: VoF overlap=%g; carrier-film separation is under-resolved\n",
               overlap);
  }
    assert(ke > -1e-10);
  // assert(ke < 1e2);
  // dump(file = "dumpTest");
}


/**
## log_hypha_deformation()

Track maximum hypha interface height using a sub-cell estimate of the
`f2 = 0.5` contour.
*/
event log_hypha_deformation (i++) {

  static double next_deformation_log = 0.;
  if (!output_due (t, &next_deformation_log, log_interval))
    return 0;

  const double f_eps = 1e-6;
  const double dy_eps = 1e-12;
  double y_if_max = -1e9;

  foreach (reduction(max:y_if_max)) {
    double f = f2[];
    if (f <= f_eps || f >= 1.0 - f_eps)
      continue; // not in interfacial band

    const double y_top = y + 0.5*Delta;
    if (y_top <= y_if_max)
      continue; // even clamped estimate cannot beat current local max

    // Default fallback is cell-center estimate.
    double y_if = y;
    double dfdy = (f2[0,1] - f2[0,-1])/(2.*Delta);

    if (fabs(dfdy) > dy_eps) {
      // Solve linearized f2(xc,y) = 0.5 and clamp to current cell bounds.
      const double y_bot = y - 0.5*Delta;
      y_if = y + (0.5 - f)/dfdy;
      if (y_if > y_top) y_if = y_top;
      else if (y_if < y_bot) y_if = y_bot;
    }

    if (y_if > y_if_max)
      y_if_max = y_if;
  }

  if (y_if_max < -1e8)
    y_if_max = NAN;

  if (pid() == 0) {
    static FILE *fh = NULL;
    static int fh_open_failed = 0;
    if (!fh && !fh_open_failed) {
      const bool fresh_run = t <= 1e-12;
      fh = fopen("hypha-def-log", fresh_run ? "w" : "a");
      if (!fh) {
        perror("hypha-def-log");
        fh_open_failed = 1;
      } else if (fresh_run)
        fprintf(fh, "t y_if_max\n");
    }

    if (fh) {
      fprintf(fh, "%g %g\n", t, y_if_max);
      fflush(fh);
    }
  }
}
