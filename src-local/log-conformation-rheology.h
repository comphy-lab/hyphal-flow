/**
# log-conformation-rheology.h

Compatibility include for the canonical hyphal-flow solver.

The constitutive update is the MultiRheoFlow tensor log-conform
implementation in [log-conform-viscoelastic.h](log-conform-viscoelastic.h),
with a three-phase `Gp` refresh before the stress divergence.
*/

#include "log-conform-viscoelastic.h"
