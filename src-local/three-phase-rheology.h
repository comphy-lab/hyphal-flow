/**
# three-phase-rheology.h

Three-material interfacial model for a viscoelastic drop, a
Kelvin--Voigt outer solid and a viscoelastic carrier liquid. Separate
VoF fields keep the drop and solid interfaces topologically distinct.

## Physical Interpretation
- Drop: `f1 = 1`, `f2 = 0`
- Outer solid: `f2 = 1`, `f1 = 0`
- Carrier liquid: `f1 = f2 = 0`

The phase weights are normalised before any material interpolation. This
keeps density, viscosity and modulus bounded if numerical smearing makes
`f1 + f2` slightly larger than one near contact. Relaxation is mixed through
the modulus-weighted relaxation rate, so the solid's very large relaxation
time does not make an entire interfacial cell artificially non-relaxing.

`Gp` and `lambda` are declared by the log-conform header; this file only
rebuilds `Gpd` and `lambdapd` and rebinds those names.
*/

#include "vof.h"

/**
Current Basilisk invalidates field stencils through `set_prolongation()`. Older
CoMPhy installs expose the same operation through the scalar's `dirty` flag.
The runner detects the installed API and defines `HYPHAL_LEGACY_BASILISK` only
for the latter.
*/
#if TREE
# ifdef HYPHAL_LEGACY_BASILISK
#  define HYPHAL_SET_PROLONGATION(field, method) do { \
     (field).prolongation = (method); \
     (field).dirty = true; \
   } while (0)
# else
#  define HYPHAL_SET_PROLONGATION(field, method) \
     set_prolongation (field, method)
# endif
#endif
/**
Instead of one VoF tracer, we define two, f1 and f2.
*/
scalar f1[], f2[], *interfaces = {f1, f2};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;
double G1 = 0., G2 = 0., G3 = 0.; // elastic moduli
double lambda1 = 0., lambda2 = 0., lambda3 = 0.; // relaxation times
/**
Auxilliary fields are necessary to define the (variable) specific
volume $\alpha=1/\rho$ as well as the cell-centered density. */

face vector alphav[];
scalar rhov[];
scalar Gpd[];
scalar lambdapd[];

event defaults (i = 0) {
  alpha = alphav;
  rho = rhov;
  Gp = Gpd;
  lambda = lambdapd;
  /**
  If the viscosity is non-zero, we need to allocate the face-centered
  viscosity field. */
  mu = new face vector;
}

/**
The density and viscosity are defined using arithmetic averages by
default. The user can overload these definitions to use other types of
averages (i.e. harmonic). The difference comes in how we call these averages.
$$
\hat{A} = (f_1+f_2) + (1-f_1-f_2)\frac{A_g}{A_l}\,\,\,\forall\,\,\,A \in \{\mu,\rho\}
$$
*/

static inline void phase_weights (double fdrop, double fsolid,
                                  double * wdrop, double * wsolid,
                                  double * wliquid)
{
  double wd = clamp (fdrop, 0., 1.);
  double ws = clamp (fsolid, 0., 1.);
  double wl = clamp (1. - wd - ws, 0., 1.);
  double total = wd + ws + wl;

  if (total <= 1e-30)
    wd = ws = 0., wl = 1., total = 1.;

  *wdrop = wd/total;
  *wsolid = ws/total;
  *wliquid = wl/total;
}

static inline double phase_average (double fdrop, double fsolid,
                                    double drop_value,
                                    double solid_value,
                                    double liquid_value)
{
  double wd, ws, wl;
  phase_weights (fdrop, fsolid, &wd, &ws, &wl);
  return wd*drop_value + ws*solid_value + wl*liquid_value;
}

static inline void phase_rheology (double fdrop, double fsolid,
                                   double * modulus,
                                   double * relaxation_time)
{
  double wd, ws, wl;
  phase_weights (fdrop, fsolid, &wd, &ws, &wl);
  double effective_modulus = wd*G1 + ws*G2 + wl*G3;
  double relaxation_rate = 0.;

  if (G1 > 0. && lambda1 > 0.)
    relaxation_rate += wd*G1/lambda1;
  if (G2 > 0. && lambda2 > 0.)
    relaxation_rate += ws*G2/lambda2;
  if (G3 > 0. && lambda3 > 0.)
    relaxation_rate += wl*G3/lambda3;

  *modulus = effective_modulus;
  *relaxation_time = effective_modulus > 0.
    ? (relaxation_rate > 0. ? effective_modulus/relaxation_rate : HUGE)
    : 0.;
}

#ifndef rho
# define rho(fdrop, fsolid) phase_average (fdrop, fsolid, rho1, rho2, rho3)
#endif
#ifndef mu
# define mu(fdrop, fsolid) phase_average (fdrop, fsolid, mu1, mu2, mu3)
#endif

/**
We have the option of using some "smearing" of the density/viscosity
jump. It is modified to take into account that there are two VoF tracers. */

#ifdef FILTERED
scalar sf1[], sf2[], *smearInterfaces = {sf1, sf2};
#else
#define sf1 f1
#define sf2 f2
scalar *smearInterfaces = {sf1, sf2};
#endif

event tracer_advection (i++) {

  /**
  When using smearing of the density jump, we initialise *sf* with the
  vertex-average of *f*. Introduce for loops to ensure that smearing is done properly. */
  #ifdef FILTERED
    int counter1 = 0;
    for (scalar sf in smearInterfaces){
      counter1++;
      int counter2 = 0;
      for (scalar f in interfaces){
        counter2++;
        if (counter1 == counter2){
          // fprintf(ferr, "%s %s\n", sf.name, f.name);
        #if dimension <= 2
            foreach(){
              sf[] = (4.*f[] +
                    2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
                    f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
            }
        #else // dimension == 3
            foreach(){
              sf[] = (8.*f[] +
                    4.*(f[-1] + f[1] + f[0,1] + f[0,-1] + f[0,0,1] + f[0,0,-1]) +
                    2.*(f[-1,1] + f[-1,0,1] + f[-1,0,-1] + f[-1,-1] +
                        f[0,1,1] + f[0,1,-1] + f[0,-1,1] + f[0,-1,-1] +
                        f[1,1] + f[1,0,1] + f[1,-1] + f[1,0,-1]) +
                    f[1,-1,1] + f[-1,1,1] + f[-1,1,-1] + f[1,1,1] +
                    f[1,1,-1] + f[-1,-1,-1] + f[1,-1,-1] + f[-1,-1,1])/64.;
            }
        #endif
        }
      }
    }
    #endif

  #if TREE
    for (scalar sf in smearInterfaces){
      HYPHAL_SET_PROLONGATION (sf, refine_bilinear);
    }
  #endif
}


event properties (i++) {

  foreach_face() {
    double ff1 = (sf1[] + sf1[-1])/2.;
    double ff2 = (sf2[] + sf2[-1])/2.;
    alphav.x[] = fm.x[]/rho(ff1, ff2);
    face vector muv = mu;
    muv.x[] = fm.x[]*mu(ff1, ff2);
  }

  foreach(){
    rhov[] = cm[]*rho(sf1[], sf2[]);

    phase_rheology (sf1[], sf2[], &Gpd[], &lambdapd[]);
  }

#if TREE
  for (scalar sf in smearInterfaces){
    HYPHAL_SET_PROLONGATION (sf, fraction_refine);
  }
#endif
}
