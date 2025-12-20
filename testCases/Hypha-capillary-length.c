/*
* Title: Flow of a drop through a single hypha branch
   Author: Vatsal Sanjay, Peter Croxford
*/

#include "navier-stokes/centered.h"
#define FILTERED
#include "three-phase-nonCoalescing-viscoelastic.h"
#include "log-conform-viscoelastic.h"
#include "tension.h"
#include "reduced-three-phase-nonCoalescing.h"

// Error tolerances
#define fErr   (1e-3)
#define KErr   (1e-4)
#define VelErr (1e-2)
#define AErr   (1e-3)
#define MINlevel 4
#define tsnap  (1e-1)
#define tsnap2 (1e-3)

int MAXlevel;
double tmax;

// Drop parameters
double Ohd, RhoR_dc, Ec_d, De_d;

// Hypha wall + film parameters
double RhoR_hc, Ohf, hf, Ec_h, De_h;

// Cytoplasm parameters
double Ohc, Ec_c, De_c;

double Ldomain;

// --- existing gap() from original code ---
#define gap(x,y,hf,width,x0,c0) (y - ((c0+1e0) + 0.5 * (hf - (c0+1e0)) * (1 + tanh(sq(x-x0) / width))))

// --- parameters used by local-gap helper ---
double width_gap  = 2.0;
double clearance_gap = 0.20;
double x0tanh_gap = 0.0;

// local film-thickness helper
static inline double local_gap(double x) {
    return gap(x, 0.0, hf, width_gap, x0tanh_gap, clearance_gap);
}

// globals to store velocities
double Ud_global = 0.0;   // droplet velocity
double Uf_global = 0.0;   // carrier fluid velocity


// ===========================================================================
// MAIN
// ===========================================================================
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
  De_h = 1e12;
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
  periodic(right);

  // Fluid properties
  rho1 = RhoR_dc; mu1 = Ohd; G1 = Ec_d; lambda1 = De_d;
  rho2 = RhoR_hc; mu2 = Ohf; G2 = Ec_h; lambda2 = De_h;
  rho3 = 1.0;     mu3 = Ohc; G3 = Ec_c; lambda3 = De_c;

 // -------------------------------
// Linear pressure profile: P(x) = Pmax - (Pmax/L)*x
// -------------------------------
double Pmax = 1.0;

// Domain goes from x = X0 to x = X0 + L0
double xL = X0;
double xR = X0 + L0;

// Pressure values at the boundaries from the desired linear law
double PL = Pmax - (Pmax/L0)*xL;
double PR = Pmax - (Pmax/L0)*xR;

// Apply Dirichlet pressure at left/right so the solver enforces that gradient
pf[left]  = dirichlet(PL);
p[left]   = dirichlet(PL);

pf[right] = dirichlet(PR);
p[right]  = dirichlet(PR);

// Keep velocity as zero-normal-gradient at inlet/outlet (typical for pressure-driven)
u.n[left]  = neumann(0);
u.t[left]  = neumann(0);
u.n[right] = neumann(0);
u.t[right] = neumann(0);


  // Surface tension
  f1.sigma = 1.0;
  f2.sigma = 1.0;

  run();
}


// ===========================================================================
// INITIALIZATION
// ===========================================================================
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


// ===========================================================================
// ADAPTIVITY
// ===========================================================================
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


// ===========================================================================
// OUTPUTS
// ===========================================================================
event writingFiles (t = 0; t += tsnap; t <= tmax+tsnap) {
  dump(file="restart");
  char nameOut[80];
  sprintf(nameOut, "intermediate/snapshot-%5.4f", t);
  dump(file=nameOut);
}


// ===========================================================================
// LOGGING DROPLET VELOCITY (U_d = vcm)
// ===========================================================================
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


// ===========================================================================
// MEASURE FLUID VELOCITY U_f
// ===========================================================================
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


// ===========================================================================
// RELATIVE-LUBRICATION DRAG:  μ (U_d - U_f) / h(x)
// ===========================================================================
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

      double mu_loc = mu2;
      double drag = mu_loc * (Ud - Uf) / hloc;

      // Limit forcing to keep solver stable
      if (fabs(drag) > 10) 
          drag = 10 * sign(drag);

      a.x[] += drag;
    }
  }
}


// ===========================================================================
// WALL-SHEAR FEEDBACK TO KELVIN–VOIGT: τ_xy = μ Uf / h(x)
// ===========================================================================
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
