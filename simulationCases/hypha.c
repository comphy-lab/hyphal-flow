/**
# hypha.c

Axisymmetric transport of a drop through a single fungal hypha branch.
Three immiscible phases remain non-coalescing: drop (`f1`), hypha film
(`f2`), and cytoplasm (`1 - f1 - f2`). Polymeric stress uses the
MultiRheoFlow log-conformation solver, with per-phase moduli and
relaxation times.

## Driving

Set `DriveMode` in the parameter file:

- `bond` (default): streamwise body force of strength `Bond` on drop
  and film, periodic in $x$. This is the production tube problem.
- `pressure`: open ends with Dirichlet pressures `PL` and `PR`.

A large film Deborah number (`De_h`) recovers a Kelvin--Voigt wall;
a finite value is Maxwell.

## Author
Vatsal Sanjay
Email: vatsal.sanjay@comphy-lab.org
CoMPhy Lab, Durham University
Last updated: 31 Aug 2026
*/

#include "axi.h"
#include "params.h"
#include "navier-stokes/centered.h"
#include "log-conform-viscoelastic-scalar-2D.h"
#define FILTERED
#include "three-phase-nonCoalescing-VE.h"
#include "tension.h"
#include "reduced-three-phase-nonCoalescing.h"

/**
## Numerical Tolerances
*/
#define fErr (1e-3)
#define KErr (1e-4)
#define VelErr (1e-2)
#define AErr (1e-3)
#define MINlevel 4
#define tsnap (1e-1)
#define tsnap2 (1e-2)

int MAXlevel;
double tmax;
int DrivePressure, NanDetector, StopOnDropExit;

/**
## Material Parameters

Drop (`d`), hypha film (`h`), and cytoplasm (`c`).
*/
double Ohd, RhoR_dc, Ec_d, De_d;
double RhoR_hc, Ohf, hf, Ec_h, De_h;
double Ohc, Ec_c, De_c;

/**
## Geometry

`gap()` is the film/cytoplasm interface. Because $\tanh(x^2)\ge 0$,
the film radius is `clearance_gap + 1` at `x = x0tanh_gap` (around the
drop) and `hf` in the far field, so the tube is wider at the drop and
narrower away from it.
*/
#define gap(x,y,hf,width,x0,c0) (y - ((c0+1e0) + 0.5 * (hf - (c0+1e0)) * (1 + tanh(sq(x-x0) / width))))

double Bond, PL, PR, Ldomain, width_gap, clearance_gap, x0tanh_gap;

static inline int bad_number (double a) {
  return isnan(a) || isinf(a);
}

/**
## main()

Load `case.params` (or `argv[1]`), set fluid properties, and enter the
Basilisk event loop.
*/
int main (int argc, char const * argv[]) {
  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system (comm);

  params_init_from_argv (argc, argv);

  MAXlevel = param_int ("MAXlevel", 12);
  tmax = param_double ("tmax", 2e2);

  const char * drive = param_string ("DriveMode", "bond");
  DrivePressure = (drive && (!strcmp (drive, "pressure") ||
                             !strcmp (drive, "Pressure") ||
                             !strcmp (drive, "PRESSURE")));
  NanDetector = param_bool ("NanDetector", false);
  StopOnDropExit = param_bool ("StopOnDropExit", true);

  Ohd = param_double ("Ohd", 1e0);
  RhoR_dc = param_double ("RhoR_dc", 1.2e0);
  Ec_d = param_double ("Ec_d", 0.0);
  De_d = param_double ("De_d", 0.0);

  Ohf = param_double ("Ohf", 1e0);
  hf = param_double ("hr", param_double ("hf", 0.90));
  Ec_h = param_double ("Ec_h", 0.0);
  De_h = param_double ("De_h", 1e30);
  RhoR_hc = param_double ("RhoR_hc", 1e0);

  Ohc = param_double ("Oh_c", param_double ("Ohc", 1e-2));
  Ec_c = param_double ("Ec_c", 0.0);
  De_c = param_double ("De_c", 0.0);

  Bond = param_double ("Bond", 1e0);
  PL = param_double ("PL", 1.0);
  PR = param_double ("PR", 0.0);
  Ldomain = param_double ("Ldomain", 80.0);
  width_gap = param_double ("width_gap", 2.0);
  clearance_gap = param_double ("clearance_gap", 0.20);
  x0tanh_gap = param_double ("x0tanh_gap", 0.0);

  fprintf (ferr,
    "Level %d tmax %g DriveMode %s. "
    "Ohd %3.2f, Ec_d %3.2f, De_d %3.2e, "
    "Ohc %3.2f, Ec_c %3.2f, De_c %3.2e, "
    "Ohf %3.2f, Ec_h %3.2f, De_h %4.3e, hf %3.2f, "
    "Bo %3.2f PL %g PR %g\n",
    MAXlevel, tmax, DrivePressure ? "pressure" : "bond",
    Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf,
    Bond, PL, PR);

  L0 = Ldomain;
  X0 = param_double ("X0", -4.0);
  Y0 = 0.0;
  init_grid (1 << 8);

  if (!DrivePressure)
    periodic (right);

  rho1 = RhoR_dc; mu1 = Ohd; G1 = Ec_d; lambda1 = De_d;
  rho2 = RhoR_hc; mu2 = Ohf; G2 = Ec_h; lambda2 = De_h;
  rho3 = 1e0;     mu3 = Ohc; G3 = Ec_c; lambda3 = De_c;

  if (DrivePressure) {
    u.n[left]  = neumann (0);
    u.t[left]  = neumann (0);
    pf[left]   = dirichlet (PL);
    p[left]    = dirichlet (PL);

    u.n[right] = neumann (0);
    u.t[right] = neumann (0);
    pf[right]  = dirichlet (PR);
    p[right]   = dirichlet (PR);

    Bf1.x = 0.;
    Bf2.x = 0.;
  } else {
    Bf1.x = Bond;
    Bf2.x = Bond;
  }

  f1.sigma = 1e0;
  f2.sigma = 1e0;

  run();
}

/**
## init()

Place the drop and film unless a `restart` snapshot is present.
*/
event init (t = 0) {
  if (!restore (file = "restart")) {
    refine (sq(y) + sq(x/1.5) > 0.81 && sq(y) + sq(x/1.5) < 1.21 &&
            level < MAXlevel);
    fraction (f1, sq(1e0) - sq(y) - sq(x/1.5));
    fraction (f2, gap(x, y, hf, width_gap, x0tanh_gap, clearance_gap));
  }
}

/**
## adapt()

Refine on both interfaces, curvature, velocity, and conformation.
*/
event adapt (i++) {
  scalar KAPPA1[], KAPPA2[], trA[];
  curvature (f1, KAPPA1);
  curvature (f2, KAPPA2);
  foreach() {
#if AXI
    trA[] = (A11[] + A22[] + AThTh[])/3.0;
#else
    trA[] = (A11[] + A22[])/2.0;
#endif
  }

  adapt_wavelet ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, trA,
                            A11, A12, A22},
                 (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, AErr,
                            AErr, AErr, AErr},
                 MAXlevel, MINlevel);

  unrefine (x > X0 + L0 - 1e0);
  unrefine (y > 1e1);
}

/**
## nan_detector()

Optional abort on the first non-finite field value. Enable with
`NanDetector=1`.
*/
event nan_detector (i++) {
  if (NanDetector) {
  if (!(dt > 0.) || bad_number (dt)) {
    fprintf (stderr, "BAD dt: dt=%g t=%g i=%d\n", dt, t, i);
    exit (10);
  }

  foreach() {
    if (bad_number (f1[]) || bad_number (f2[]) ||
        bad_number (u.x[]) || bad_number (u.y[]) || bad_number (p[]) ||
        bad_number (A11[]) || bad_number (A12[]) || bad_number (A22[])) {
      fprintf (stderr,
               "NaN/Inf at x=%g y=%g t=%g i=%d "
               "f1=%g f2=%g u=(%g,%g) p=%g A=(%g,%g,%g)\n",
               x, y, t, i, f1[], f2[], u.x[], u.y[], p[],
               A11[], A12[], A22[]);
      dump (file = "nan");
      exit (11);
    }
  }
  }
}

/**
## stop_when_drop_exits()

Stop when the drop leading edge approaches the outlet.
*/
event stop_when_drop_exits (t += tsnap2) {
  if (StopOnDropExit) {
  scalar xpos[];
  coord ex = {1., 0.};
  coord z0 = {0., 0.};
  position (f1, xpos, ex, z0, add = false);

  stats sx = statsf (xpos);
  double xmax = sx.volume > 0. ? sx.max : -HUGE;
  double finest = L0/(1 << MAXlevel);
  double buffer = 2.*finest;
  double x_end = X0 + L0;

  if (pid() == 0)
    fprintf (ferr, "drop front xmax = %.6f, x_end = %.6f\n", xmax, x_end);

  if (xmax > x_end - buffer) {
    if (pid() == 0) {
      fprintf (ferr,
        "\n*** Drop leading edge reached end of domain ***\n"
        "xmax = %.6f, domain end = %.6f\n"
        "Stopping simulation at t = %.6g\n\n",
        xmax, x_end, t);
    }
    dump (file = "final");
    exit (0);
  }
  }
}

/**
## writingFiles()

Periodic restart and snapshot dumps.
*/
event writingFiles (t = 0, t += tsnap; t <= tmax + tsnap) {
  dump (file = "restart");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

/**
## logWriting()

Kinetic energy and drop centre-of-mass velocity. Volume weights use
`dv()`, so the axisymmetric metric is included.
*/
event logWriting (t = 0, t += tsnap2; t <= tmax + tsnap) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke) reduction(+:vcm) reduction(+:wt)) {
    ke += 0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[]))*dv();
    vcm += f1[]*u.x[]*dv();
    wt += f1[]*dv();
  }
  if (wt > 0.0)
    vcm /= wt;

  if (pid() == 0) {
    static FILE * fp;
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm\n");
      fp = fopen ("log", "w");
      fprintf (fp,
        "Level %d tmax %g DriveMode %s. "
        "Ohd %3.2f, Ec_d %3.2f, De_d %3.2e, "
        "Ohc %3.2f, Ec_c %3.2f, De_c %3.2e, "
        "Ohf %3.2f, Ec_h %3.2f, De_h %4.3e, hf %3.2f, Bo %3.2f\n",
        MAXlevel, tmax, DrivePressure ? "pressure" : "bond",
        Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf, Bond);
      fprintf (fp, "i dt t ke vcm\n");
    } else {
      fp = fopen ("log", "a");
    }
    fprintf (fp, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
    fclose (fp);
    fprintf (ferr, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
  }
  assert (ke > -1e-10);
}

/**
## log_hypha_deformation()

Track the maximum film height from a sub-cell estimate of the
`f2 = 0.5` contour.
*/
event log_hypha_deformation (t = 0, t += tsnap2; t <= tmax + tsnap) {
  const double f_eps = 1e-6;
  const double dy_eps = 1e-12;
  double y_if_max = -1e9;

  foreach (reduction(max:y_if_max)) {
    double f = f2[];
    if (f <= f_eps || f >= 1.0 - f_eps)
      continue;

    const double y_top = y + 0.5*Delta;
    if (y_top <= y_if_max)
      continue;

    double y_if = y;
    double dfdy = (f2[0,1] - f2[0,-1])/(2.*Delta);
    if (fabs(dfdy) > dy_eps) {
      const double y_bot = y - 0.5*Delta;
      y_if = y + (0.5 - f)/dfdy;
      if (y_if > y_top)
        y_if = y_top;
      else if (y_if < y_bot)
        y_if = y_bot;
    }
    if (y_if > y_if_max)
      y_if_max = y_if;
  }

  if (y_if_max < -1e8)
    y_if_max = NAN;

  if (pid() == 0) {
    static FILE * fh = NULL;
    static int fh_open_failed = 0;
    if (!fh && !fh_open_failed) {
      fh = fopen ("hypha-def-log", "w");
      if (!fh) {
        perror ("hypha-def-log");
        fh_open_failed = 1;
      } else {
        fprintf (fh, "t y_if_max\n");
      }
    }
    if (fh) {
      fprintf (fh, "%g %g\n", t, y_if_max);
      fflush (fh);
    }
  }
}
