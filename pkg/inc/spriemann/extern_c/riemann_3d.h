
#pragma once

#include <spriemann/constants/types.h>

// Recall that when Fortran is compiled to
// library code using gfortran, symbols are
// coerced to lower case (My_Name -> my_name)
// and appended with an underscore '_'.

extern "C" {

// this subroutine uses the primitive variables "Vpt" to build the
// extremal wave speeds "S_min, S_max" in the "idirection".
void prim_to_wavspd_(
        int32* idirection,
        // vectors (zone data Fields)
        real64* gamma,
        real64* Vpt,
        // OUTPUT:
        real64* S_min,
        real64* S_max,
        // INPUT:
        real64* B0
);


void prim_to_cons_(
        real64* gamma,
        real64* gamma_euler,
        real64* isospeed,
        real64* barotropic_stiffness,
        real64* molewt,
        real64* gasconst,
        real64* rhofloor,
        real64* prsfloor,
        real64* tempfloor,
        real64* magB_small_cgl,
        real64* eps_bdry_cgl,
        int32* iprs_flag,
        real64* ent,
        real64* B0,
        // OUTPUT:
        real64* Upt,
        // INPUT:
        real64* Vpt
);

void cons_to_prim_(
        real64* gamma,
        real64* gamma_euler,
        real64* isospeed,
        real64* barotropic_stiffness,
        real64* molewt,
        real64* gasconst,
        real64* rhofloor,
        real64* prsfloor,
        real64* tempfloor,
        real64* magB_small_cgl,
        real64* eps_bdry_cgl,
        int32* iprs_flag,
        real64* ent,
        real64* B0,
        // INPUT:
        real64* Upt,
        // OUTPUT:
        real64* Vpt
);

void force_pperp_ppar_to_hyperbolic_domain_(
        real64* pi,
        real64* magB_small_cgl,
        real64* eps_bdry_cgl,
        real64* prsfloor,
        real64* Bsq,
        real64* p_par,
        real64* p_perp,
        real64* updated_p_par,
        real64* updated_p_perp
);




}



