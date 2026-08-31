/**
# Three-Phase Non-Coalescing Viscoelasticity

VoF coupling between MultiRheoFlow log-conformation rheology and three
immiscible phases that do not coalesce.

The two liquid phases occupy distinct volume-fraction fields `f1` and
`f2`. The third phase is the remainder `1 - f1 - f2`. A thin precursor
of phase 3 is assumed to remain between phases 1 and 2, so the two
interfaces never merge.

## Phase Map

- Drop: `f1 = 1`, `f2 = 0`
- Hypha film/wall fluid: `f2 = 1`, `f1 = 0`
- Cytoplasm / ambient: `f1 = f2 = 0`

Densities and solvent viscosities are `rho1`, `mu1`, `rho2`, `mu2`,
`rho3`, `mu3`. Elastic moduli and relaxation times are `G1`, `lambda1`
and likewise for phases 2 and 3. Those scalars are interpolated onto
the log-conform fields `Gp` and `lambda` declared by
[log-conform-viscoelastic-scalar-2D.h](log-conform-viscoelastic-scalar-2D.h).

## Provenance

Adapted from MultiRheoFlow [two-phaseVE.h](https://github.com/comphy-lab/MultiRheoFlow)
(commit 7d9c3df) and the previous hyphal-flow three-phase non-coalescing
header. Include the log-conform header *before* this file.

## Author
Vatsal Sanjay
Email: vatsal.sanjay@comphy-lab.org
CoMPhy Lab, Durham University
Last updated: 31 Aug 2026
*/

#include "vof.h"

scalar f1[], f2[], * interfaces = {f1, f2};

double rho1 = 1., mu1 = 0., rho2 = 1., mu2 = 0., rho3 = 1., mu3 = 0.;
double G1 = 0., G2 = 0., G3 = 0.;
double lambda1 = 0., lambda2 = 0., lambda3 = 0.;
double TOLelastic = 1e-2;

/**
Auxiliary fields define the specific volume $\alpha = 1/\rho$ and the
cell-centred density. `Gpd` and `lambdapd` overwrite the constant
`Gp` and `lambda` placeholders declared by the log-conform header.
*/

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
  If the viscosity is non-zero, allocate the face-centred viscosity
  field. */
  mu = new face vector;
}

/**
Density and solvent viscosity use arithmetic averages of the three
phases. Harmonic means can be supplied by redefining `rho()` / `mu()`.
*/

#ifndef rho
# define rho(f1, f2) (clamp(f1,0.,1.)*rho1 + clamp(f2,0.,1.)*rho2 + clamp(1. - f1 - f2,0.,1.)*rho3)
#endif
#ifndef mu
# define mu(f1, f2)  (clamp(f1,0.,1.)*mu1 + clamp(f2,0.,1.)*mu2 + clamp(1. - f1 - f2,0.,1.)*mu3)
#endif

#ifdef FILTERED
scalar sf1[], sf2[], * smearInterfaces = {sf1, sf2};
#else
# define sf1 f1
# define sf2 f2
scalar * smearInterfaces = {sf1, sf2};
#endif

event tracer_advection (i++) {
#ifdef FILTERED
  int counter1 = 0;
  for (scalar sf in smearInterfaces) {
    counter1++;
    int counter2 = 0;
    for (scalar f in interfaces) {
      counter2++;
      if (counter1 != counter2)
        continue;
#if dimension <= 2
      foreach() {
        sf[] = (4.*f[] +
                2.*(f[0,1] + f[0,-1] + f[1,0] + f[-1,0]) +
                f[-1,-1] + f[1,-1] + f[1,1] + f[-1,1])/16.;
      }
#else
      foreach() {
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
#endif

#if TREE
  for (scalar sf in smearInterfaces)
    set_prolongation (sf, refine_bilinear);
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

  foreach() {
    rhov[] = cm[]*rho(sf1[], sf2[]);

    Gpd[] = 0.;
    lambdapd[] = 0.;

    double c1 = clamp(sf1[], 0., 1.);
    double c2 = clamp(sf2[], 0., 1.);
    double c3 = clamp(1. - c1 - c2, 0., 1.);

    /**
    Arithmetic mixing of $G_p$ and $\lambda$ matches MultiRheoFlow
    `two-phaseVE.h`. A harmonic mix of $1/\lambda$ would be a separate
    constitutive change. */
    if (c1 > TOLelastic) {
      Gpd[] += G1*c1;
      lambdapd[] += lambda1*c1;
    }
    if (c2 > TOLelastic) {
      Gpd[] += G2*c2;
      lambdapd[] += lambda2*c2;
    }
    if (c3 > TOLelastic) {
      Gpd[] += G3*c3;
      lambdapd[] += lambda3*c3;
    }
  }

#if TREE
  for (scalar sf in smearInterfaces)
    set_prolongation (sf, fraction_refine);
#endif
}
