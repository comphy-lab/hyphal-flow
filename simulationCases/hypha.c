/**
# hypha.c

Baseline single-branch hypha simulation with non-coalescing three-phase
viscoelastic coupling and body-force driven transport.

## Author
Vatsal Sanjay

## Version
- August 11, 2024: allow stress relaxation in both drop and
  cytoplasm.
*/

// #include "axi.h"
#include "params.h"
#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase-nonCoalescing-viscoelastic.h"
#include "log-conform-viscoelastic.h"
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
#define tsnap (1e-1)
#define tsnap2 (1e-2)


int MAXlevel;
double tmax;

/**
## Material Parameters

Drop (`d`), hypha film/wall (`h`), and cytoplasm (`c`) parameter sets.
*/
double Ohd, RhoR_dc, Ec_d, De_d;
double RhoR_hc, Ohf, hf, Ec_h, De_h;
double Ohc, Ec_c, De_c;

/**
## Geometry Helper
*/
#define gap(x,y,hf,width,x0,c0) (y - ((c0+1e0) + 0.5 * (hf - (c0+1e0)) * (1 + tanh(sq(x-x0) / width))))

double Bond; // driving force
double Ldomain;


/**
## main()

Initialize properties and forcing, then enter the Basilisk event loop.
*/
int main(int argc, char const *argv[]) {

  char comm[80];
  sprintf (comm, "mkdir -p intermediate");
  system(comm);

  params_init_from_argv(argc, argv);

  MAXlevel = param_int("MAXlevel", 12);
  tmax = param_double("tmax", 2e2);

  // Drop
  Ohd = param_double("Ohd", 1e0); // <0.000816/sqrt(816*0.017*0.00075) = 0.008>
  RhoR_dc = param_double("RhoR_dc", 1.2e0);
  Ec_d = param_double("Ec_d", 0.0); //atof(argv[3]);
  De_d = param_double("De_d", 0.0); //atof(argv[4]);

  // Hypha
  Ohf = param_double("Ohf", 1e0);
  hf = param_double("hr", param_double("hf", 0.90)); // ratio of the gap thickness to the drop radius, far awaay from the drop.
  Ec_h = param_double("Ec_h", 0.0); //atof(argv[6]); // Elasto-capillary number: 1e-4 (very soft), (0.1 is very similar to 0) to 1e1 (appears to be rigid)
  De_h = param_double("De_h", 1e30); //atof(argv[7]);
  RhoR_hc = 1e0; // density ratio of hypha to cytoplasm

  // cytoplasm
  Ohc = param_double("Oh_c", param_double("Ohc", 1e-2));
  Ec_c = param_double("Ec_c", 0.00);
   //atof(argv[8]);
  De_c = param_double("De_c", 0.0); //atof(argv[9]);

  Bond = param_double("Bond", 1e0); // Bond number: we keep the driving fixed

  Ldomain = param_double("Ldomain", 80.0); // Dimension of the domain: should be large enough to get a steady solution to drop velocity.

  fprintf(ferr, "Level %d tmax %g. Ohd %3.2f, Ec_d %3.2f, De_d %3.2e, Ohc %3.2f, Ec_c %3.2f, De_c %3.2e, Ohf %3.2f, Ec_h %3.2f, De_h %4.3e, hf %3.2f, Bo %3.2f\n", MAXlevel, tmax, Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf, Bond);

  L0=Ldomain;
  X0=-4.0; Y0=0.0;
  init_grid (1 << (8));
  periodic(right);

  // drop
  rho1 = RhoR_dc; mu1 = Ohd; G1 = Ec_d; lambda1 = De_d;

  // Hypha
  rho2 = RhoR_hc; mu2 = Ohf; G2 = Ec_h; lambda2 = De_h;

  // cytoplasm
  rho3 = 1e0; mu3 = Ohc; G3 = Ec_c; lambda3 = De_c;

  Bf1.x = Bond;
  Bf2.x = Bond;

  f1.sigma = 1e0; // drop-cytoplasm interfacial tension
  f2.sigma = 1e0; // hypha-cytoplasm interfacial tension

  run();

}

/**
## init()

Initialize interfaces unless a restart snapshot is available.
*/
event init(t = 0){
  if (!restore (file = "restart")) {
    double width = 2e0; // width of the tanh function
    double clearance = 0.20; // clearance from the drop inside the cytoplasm at t = 0
    double x0tanh = 0.0; // midpoint of the tanh function

    refine (sq(y) + sq(x/1.5) > 0.81 && sq(y) + sq(x/1.5) < 1.21 && level < MAXlevel);

    fraction(f1, sq(1e0) - sq(y) - sq(x/1.5));
    fraction(f2, gap(x,y,hf,width,x0tanh,clearance));

  }
}

/**
## adapt()

Adaptive mesh refinement driven by interfaces, curvature, velocity, and
conformation fields.
*/
event adapt(i++){
  scalar KAPPA1[], KAPPA2[];
  curvature(f1, KAPPA1);
  curvature(f2, KAPPA2);

  adapt_wavelet ((scalar *){f1, f2, KAPPA1, KAPPA2, u.x, u.y, conform_p.x.x, conform_p.x.y, conform_p.y.y},
  (double[]){fErr, fErr, KErr, KErr, VelErr, VelErr, AErr, AErr, AErr},
  MAXlevel, MINlevel);

  unrefine(x > L0-1e0);
  unrefine(y > 1e1);
}
/**
## stop_when_drop_exits()

Stop the run when the leading drop edge approaches the domain outlet.
*/
event stop_when_drop_exits (t += tsnap2) {

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
    exit(0);
  }
}

/**
## writingFiles()

Write periodic restart and snapshot files.
*/
event writingFiles (t = 0, t += tsnap; t <= tmax+tsnap) {
  dump (file = "restart");
  char nameOut[80];
  sprintf (nameOut, "intermediate/snapshot-%5.4f", t);
  dump (file = nameOut);
}

/**
## logWriting()

Log kinetic energy and droplet center-of-mass velocity.
*/
event logWriting (t = 0, t += tsnap2; t <= tmax+tsnap) {
  double ke = 0., vcm = 0., wt = 0.;
  foreach (reduction(+:ke) reduction(+:vcm) reduction(+:wt)){
    ke += (0.5*rho(f1[], f2[])*(sq(u.x[]) + sq(u.y[])))*sq(Delta);
    vcm += (f1[]*u.x[])*sq(Delta);
    wt += f1[]*sq(Delta);
  }
  if (wt > 0.0) vcm /= wt;
  static FILE * fp;

  if (pid() == 0){
    if (i == 0) {
      fprintf (ferr, "i dt t ke vcm\n");
      fp = fopen ("log", "w");
      fprintf(fp, "Level %d tmax %g. Ohd %3.2f, Ec_d %3.2f, De_d %3.2e, Ohc %3.2f, Ec_c %3.2f, De_c %3.2e, Ohf %3.2f, Ec_h %3.2f, De_h %4.3e, hf %3.2f, Bo %3.2f\n", MAXlevel, tmax, Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf, Bond);
      fprintf (fp, "i dt t ke vcm\n");
    } else {
      fp = fopen ("log", "a");
    }
    fprintf (fp, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
    fclose(fp);
    fprintf (ferr, "%d %g %g %g %5.4e\n", i, dt, t, ke, vcm);
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
event log_hypha_deformation (t = 0; t += tsnap2; t <= tmax + tsnap) {

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
      fh = fopen("hypha-def-log","w");
      if (!fh) {
        perror("hypha-def-log");
        fh_open_failed = 1;
      } else
        fprintf(fh, "t y_if_max\n");
    }

    if (fh) {
      fprintf(fh, "%g %g\n", t, y_if_max);
      fflush(fh);
    }
  }
}
