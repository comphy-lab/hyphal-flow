/**
# Hypha-capillary-length.c

Pressure-driven drop transport through a hypha branch with additional
runtime stability diagnostics and deformation tracking.

## Authors
- Vatsal Sanjay
- Peter Croxford
*/

#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase-nonCoalescing-viscoelastic.h"
#include "log-conform-viscoelastic.h"
#include "tension.h"
#include "reduced-three-phase-nonCoalescing.h"

/**
## Numerical Tolerances
*/
#define fErr   (1e-3)
#define KErr   (1e-4)
#define VelErr (1e-2)
#define AErr   (1e-3)
#define MINlevel 4
#define tsnap  (1e-1)
#define tsnap2 (1e-3)

int MAXlevel;
double tmax;
double PL;
double PR;
double Pmax;
double beta_bf;
double mu_ref;
double rho_ref;


/**
## Material Parameters

Drop (`d`), hypha film/wall (`h`), and cytoplasm (`c`) parameter sets.
*/
double Ohd, RhoR_dc, Ec_d, De_d;
double RhoR_hc, Ohf, hf, Ec_h, De_h;
double Ohc, Ec_c, De_c;

double Ldomain;

/**
## Geometry Helper

`gap()` defines the wall profile and `local_gap()` is used by drag and
wall-shear closures.
*/
#define gap(x,y,hf,width,x0,c0) (y - ((c0+1e0) + 0.5 * (hf - (c0+1e0)) * (1 + tanh(sq(x-x0) / width))))
double width_gap  = 2.0;
double clearance_gap = 0.20;
double x0tanh_gap = 0.0;
static inline double local_gap(double x) {
    return gap(x, 0.0, hf, width_gap, x0tanh_gap, clearance_gap);
}
double Ud_global = 0.0;   // droplet velocity
double Uf_global = 0.0;   // carrier fluid velocity

/**
## Diagnostics Utilities
*/
#include <math.h>

static inline int bad (double a) { return isnan(a) || isinf(a); }

/**
## first_nan_detector()

Abort early when invalid values appear in primary state variables,
curvature fields, or conformation tensors.
*/
event first_nan_detector (i++) {

  // Catch bad timestep too
  if (!(dt > 0.) || bad(dt)) {
    fprintf(stderr, "BAD dt: dt=%g t=%g i=%d\n", dt, t, i);
    exit(10);
  }
foreach_face(x) {
    if (isnan(u.x[]) || isinf(u.x[])) {
      fprintf(stderr, "NaN/Inf in FACE u.x at xf=%g y=%g t=%g i=%d ux=%g\n", x, y, t, i, u.x[]);
      exit(30);
    }
  }
  foreach_face(y) {
    if (isnan(u.y[]) || isinf(u.y[])) {
      fprintf(stderr, "NaN/Inf in FACE u.y at x=%g yf=%g t=%g i=%d uy=%g\n", x, y, t, i, u.y[]);
      exit(31);
    }
  }
  /**
  Compute curvatures so they can be included in stability checks.
  */
  scalar K1[], K2[];
  curvature(f1, K1);
  curvature(f2, K2);

  /**
  Scan the domain and terminate on the first invalid value.
  */
  foreach() {

    // Volume fractions
    if (bad(f1[])) {
      fprintf(stderr, "NaN/Inf in f1 at x=%g y=%g t=%g i=%d f1=%g\n", x, y, t, i, f1[]);
      exit(11);
    }
    if (bad(f2[])) {
      fprintf(stderr, "NaN/Inf in f2 at x=%g y=%g t=%g i=%d f2=%g\n", x, y, t, i, f2[]);
      exit(12);
    }

    // Velocity and pressure
    if (bad(u.x[])) {
      fprintf(stderr, "NaN/Inf in u.x at x=%g y=%g t=%g i=%d u.x=%g\n", x, y, t, i, u.x[]);
      // Helpful local context
      fprintf(stderr, "  u.y=%g p=%g f1=%g f2=%g Delta=%g\n", u.y[], p[], f1[], f2[], Delta);
      exit(13);
    }
    if (bad(u.y[])) {
      fprintf(stderr, "NaN/Inf in u.y at x=%g y=%g t=%g i=%d u.y=%g\n", x, y, t, i, u.y[]);
      fprintf(stderr, "  u.x=%g p=%g f1=%g f2=%g Delta=%g\n", u.x[], p[], f1[], f2[], Delta);
      exit(14);
    }
    if (bad(p[])) {
      fprintf(stderr, "NaN/Inf in p at x=%g y=%g t=%g i=%d p=%g\n", x, y, t, i, p[]);
      fprintf(stderr, "  u.x=%g u.y=%g f1=%g f2=%g Delta=%g\n", u.x[], u.y[], f1[], f2[], Delta);
      exit(15);
    }

    // Curvatures (often the first to blow up if interface gets underresolved)
    if (bad(K1[])) {
      fprintf(stderr, "NaN/Inf in curvature K1(f1) at x=%g y=%g t=%g i=%d K1=%g\n", x, y, t, i, K1[]);
      fprintf(stderr, "  f1=%g f2=%g Delta=%g\n", f1[], f2[], Delta);
      exit(16);
    }
    if (bad(K2[])) {
      fprintf(stderr, "NaN/Inf in curvature K2(f2) at x=%g y=%g t=%g i=%d K2=%g\n", x, y, t, i, K2[]);
      fprintf(stderr, "  f1=%g f2=%g Delta=%g\n", f1[], f2[], Delta);
      exit(17);
    }

    // Viscoelastic fields (common NaN source in log-conformation if SPD violated)
    if (bad(conform_p.x.x[])) {
      fprintf(stderr, "NaN/Inf in conform_p.x.x at x=%g y=%g t=%g i=%d Cxx=%g\n", x, y, t, i, conform_p.x.x[]);
      fprintf(stderr, "  Cyy=%g Cxy=%g f2=%g\n", conform_p.y.y[], conform_p.x.y[], f2[]);
      exit(18);
    }
    if (bad(conform_p.y.y[])) {
      fprintf(stderr, "NaN/Inf in conform_p.y.y at x=%g y=%g t=%g i=%d Cyy=%g\n", x, y, t, i, conform_p.y.y[]);
      fprintf(stderr, "  Cxx=%g Cxy=%g f2=%g\n", conform_p.x.x[], conform_p.x.y[], f2[]);
      exit(19);
    }
    if (bad(conform_p.x.y[]) || bad(conform_p.y.x[])) {
      fprintf(stderr, "NaN/Inf in conform_p shear at x=%g y=%g t=%g i=%d Cxy=%g Cyx=%g\n",
              x, y, t, i, conform_p.x.y[], conform_p.y.x[]);
      fprintf(stderr, "  Cxx=%g Cyy=%g f2=%g\n", conform_p.x.x[], conform_p.y.y[], f2[]);
      exit(20);
    }

    // Optional SPD check (only meaningful in phase-2 regions).
    if (f2[] > 0.5) {
      double cxx = conform_p.x.x[];
      double cyy = conform_p.y.y[];
      double cxy = 0.5*(conform_p.x.y[] + conform_p.y.x[]);
      double det = cxx*cyy - cxy*cxy;
      if (!(cxx > 0.) || !(cyy > 0.) || !(det > 0.)) {
        fprintf(stderr, "NON-SPD conformation at x=%g y=%g t=%g i=%d: Cxx=%g Cyy=%g Cxy=%g det=%g\n",
                x, y, t, i, cxx, cyy, cxy, det);
        fprintf(stderr, "  (This often causes log-conformation to NaN soon after)\n");
        exit(21);
      }
    }
  }
}

/**
## main()

Initialize parameters, pressure boundary conditions, and run the
simulation event loop.
*/
int main(int argc, char const *argv[]) {

  system("mkdir -p intermediate");

  MAXlevel = 9;
  tmax = 1e2;
  
  
  // Drop
  Ohd = 1e0;
  RhoR_dc = 1.2;
  Ec_d = 0.0;
  De_d = 0.0;

  // Hypha film + wall
  Ohf = 1e0;
  hf = 0.90;
  Ec_h = 1e0;
  De_h = 0.1;
  RhoR_hc = 1e0;

  // cytoplasm
  Ohc = 1e-2;
  Ec_c = 0.0;
  De_c = 0.0;

  Ldomain = 16.0;

  fprintf(ferr,
    "Level %d tmax %g. Ohd %g Ec_d %g De_d %g Ohc %g Ec_c %g De_c %g Ohf %g Ec_h %g De_h %g hf %g\n",
     MAXlevel, tmax, Ohd, Ec_d, De_d, Ohc, Ec_c, De_c, Ohf, Ec_h, De_h, hf);

  L0 = Ldomain;
  X0 = -4.0;
  Y0 = 0.0;

  init_grid(1 << MINlevel);
  

  // Fluid properties
  rho1 = RhoR_dc; mu1 = Ohd; G1 = Ec_d; lambda1 = De_d;
  rho2 = RhoR_hc; mu2 = Ohf; G2 = Ec_h; lambda2 = De_h;
  rho3 = 1.0;     mu3 = Ohc; G3 = Ec_c; lambda3 = De_c;

rho_ref = rho2;
mu_ref  = mu2;



  /**
  Linear pressure profile: `P(x) = Pmax - (Pmax/L0) x`.
  */
Pmax = 1.0;
beta_bf = 0.01;
  // Re-assert fluid properties used for the boundary closure.
  rho1 = RhoR_dc; mu1 = Ohd; G1 = Ec_d; lambda1 = De_d;
  rho2 = RhoR_hc; mu2 = Ohf; G2 = Ec_h; lambda2 = De_h;
  rho3 = 1.0;     mu3 = Ohc; G3 = Ec_c; lambda3 = De_c;

rho_ref = rho2;
mu_ref  = mu2;

// Helper: local coeff (beta*rho/mu) and backflow indicator (u·n)^-
#define BF_COEFF_REF (beta_bf * rho_ref / (mu_ref + 1e-30))
#define UN_MINUS (max(-u.n[], 0.0))


  // Domain goes from x = X0 to x = X0 + L0.
double xL = X0;
double xR = X0 + L0;

  // Pressure values at the boundaries from the linear law.
PL = Pmax - (Pmax/L0)*xL;
PR = Pmax - (Pmax/L0)*xR;

  // Apply Dirichlet pressure at left/right boundaries.
pf[left]  = dirichlet(PL);
p[left]   = dirichlet(PL);

pf[right] = dirichlet(PR);
p[right]  = dirichlet(PR);

  // Keep velocity as zero-normal-gradient at inlet/outlet.
u.n[left]  = neumann( BF_COEFF_REF * UN_MINUS * u.n[] );
u.n[right] = neumann( BF_COEFF_REF * UN_MINUS * u.n[] );
u.t[left]  = neumann( BF_COEFF_REF * UN_MINUS * u.t[] );
u.t[right] = neumann( BF_COEFF_REF * UN_MINUS * u.t[] );


  // Surface tension
  f1.sigma = 1.0;
  f2.sigma = 1.0;

  run();
}


/**
## init()

Create initial drop and hypha interfaces unless a restart snapshot is
available.
*/
event init(t = 0) {
  if (!restore (file = "restart")) {

    double width     = 2.0;
    double clearance = 0.20;
    double x0tanh    = 0.0;

    refine(y < 1.2 && level < MAXlevel);

    // droplet
    fraction(f1, sq(1.0) - sq(y) - sq(x/1.5));

    // hypha
    fraction(f2, gap(x,y,hf,width,x0tanh,clearance));
  }
}


/**
## adapt()

Adaptive mesh refinement driven by interfaces, curvature, velocity, and
conformation fields.
*/
event adapt(i++) {
  scalar K1[], K2[];
  curvature(f1, K1);
  curvature(f2, K2);

  adapt_wavelet ((scalar *){f1, f2, K1, K2, u.x, u.y,
                            conform_p.x.x, conform_p.x.y, conform_p.y.y},
                 (double[]){fErr,fErr,KErr,KErr,VelErr,VelErr,AErr,AErr,AErr},
                 MAXlevel, MINlevel);

  unrefine(x > L0 - 1.0);
}


/**
## writingFiles()

Write periodic restart and snapshot files.
*/
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump(file="restart");
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}


/**
## logWriting()

Compute kinetic energy and droplet center-of-mass velocity `U_d`.
*/
event logWriting (t = 0; t += tsnap2; t <= tmax+tsnap) {
  double ke = 0., vcm = 0., wt = 0.;

  foreach (reduction(+:ke) reduction(+:vcm) reduction(+:wt)) {
    ke += (0.5*rho(f1[],f2[]) * (sq(u.x[]) + sq(u.y[]))) * sq(Delta);
    vcm += (f1[]*u.x[]) * sq(Delta);
    wt  += f1[] * sq(Delta);
  }
  if (wt > 0.0) vcm /= wt;
  Ud_global = vcm;   // store droplet velocity

  if (pid() == 0) {
    static FILE *fp = NULL;
    if (!fp) {
      fp = fopen("log","w");
      fprintf(fp,"i dt t ke U_d\n");
    } else fp = fopen("log","a");

    fprintf(fp,"%d %g %g %g %g\n", i, dt, t, ke, Ud_global);
    fclose(fp);
  }
}

/**
## log_hypha_deformation()

Track the maximum interface height using a sub-cell linearized estimate
for the `f2 = 0.5` contour.
*/
event log_hypha_deformation (t = 0; t += tsnap2; t <= tmax + tsnap) {

  double y_if_max = -1e9;

  foreach (reduction(max:y_if_max)) {
    // interfacial band
    if (f2[] > 1e-6 && f2[] < 1.0 - 1e-6) {

      // Local gradient of f2 (central differences)
      double dfdx = (f2[1,0] - f2[-1,0])/(2.*Delta);
      double dfdy = (f2[0,1] - f2[0,-1])/(2.*Delta);

      // If gradient is too small, fall back to cell center
      double g = sqrt(dfdx*dfdx + dfdy*dfdy);
      double y_if = y;

      if (g > 1e-12 && fabs(dfdy) > 1e-12) {
        // Linearize: f2(x,y) ~ f2[] + dfdx*(x-xc) + dfdy*(y-yc)
        // Solve for y where f2 = 0.5 along vertical line through cell center:
        y_if = y + (0.5 - f2[])/dfdy;
        // Clamp to cell bounds (avoid crazy jumps)
        if (y_if > y + 0.5*Delta) y_if = y + 0.5*Delta;
        if (y_if < y - 0.5*Delta) y_if = y - 0.5*Delta;
      }

      if (y_if > y_if_max)
        y_if_max = y_if;
    }
  }

  if (y_if_max < -1e8)
    y_if_max = NAN;

  if (pid() == 0) {
    static FILE *fh = NULL;
    if (!fh) {
      fh = fopen("hypha-def-log","w");
      fprintf(fh, "t y_if_max\n");
    } else
      fh = fopen("hypha-def-log","a");

    fprintf(fh, "%g %g\n", t, y_if_max);
    fclose(fh);
  }
}



/**
## measure_Uf()

Compute average carrier-flow velocity `U_f` in the hypha phase.
*/
event measure_Uf (t += tsnap2) {

  double Uf = 0., wt = 0.;
  foreach(reduction(+:Uf) reduction(+:wt)) {
    if (f2[] > 0.5) {         // hypha fluid
      Uf += u.x[] * f2[] * sq(Delta);
      wt += f2[] * sq(Delta);
    }
  }
  if (wt > 0.) Uf /= wt;
  Uf_global = Uf;

  if (pid()==0) {
    static FILE *fU = NULL;
    if (!fU) fU = fopen("Uf-log","w");
    fprintf(fU,"%g %g\n", t, Uf_global);
    fflush(fU);
  }
}


/**
## relative_drag()

Apply a lubrication-inspired body-force correction
$\mu (U_d - U_f)/h_{\mathrm{eff}}(x)$.
*/
#define sign(x) ((x > 0) - (x < 0))
event relative_drag (i++) {

  double Ud = Ud_global;
  double Uf = Uf_global;

  foreach_face(x) {
    // Only in hypha fluid but away from the wall-interface
    if (fm.x[] > 0 && f2[] > 0.1 && f2[] < 0.9) {

      double xloc = x;
      double hloc = fabs(local_gap(xloc));
      if (hloc < 1e-6) hloc = 1e-6;
      double h0 = Delta;
      double h_eff = hloc + h0;

      double mu_loc = mu2;
      double drag = mu_loc * (Ud - Uf) / h_eff;

      // Limit forcing to keep solver stable
      if (fabs(drag) > 10) 
          drag = 10 * sign(drag);

      a.x[] += drag;
    }

  }
}


/**
## wall_shear()

Feed back estimated wall shear into the Kelvin-Voigt stress state via
`conform_p`.
*/
event wall_shear (i++) {
  double Uf = Uf_global;

  foreach() {
    if (f2[] > 0.5) {

      double xloc = x;
      double hloc = fabs(local_gap(xloc));
      if (hloc < 1e-6) hloc = 1e-6;

      double mu_loc = mu2;
      double shear = mu_loc * Uf / hloc;

      conform_p.x.y[] += shear * dt;
      conform_p.y.x[] += shear * dt;
    }
  }
 }
