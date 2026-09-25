#pragma once

#include <csriemann/constants/types.h>
#include <csriemann/constants/riemannconstants.h>

#include <csriemann/extern_c/riemann_3d.h>




// Stuff that comes up
// while doing this.

// THE WAVE SPEEDS
struct WaveSpeed {
    real64 S_min;
    real64 S_max;
};





template <
  typename Parameters
>
class FRiemann3D {
public: // types:

    using P = Parameters;

    using RC = RiemannConstants<P>;

    using FluidVector = P::FluidVector;
    using PParPPerp_t = P::PParPPerp_t;


     // "IPDE_SYSTEM == 1" for MHD
     // "IPDE_SYSTEM == 2" for CGL System
     // "IPDE_SYSTEM == 3" for CGL System with p_elec
    static constexpr IPDESystem ipde_system = P::Directives::ipde_system;
    static constexpr int32 numorder = P::Directives::numorder;


public: // fields:



public: // ctor:


    FRiemann3D(


    ) {



    }


public: // API:




    WaveSpeed PrimToWavspd(
            int32 idirection,
            // todo what constructs Vpt???
            FluidVector& Vpt,
            GeoVector& B0
    ) {
        WaveSpeed wavespeed;
        // constructs S_min, S_max using Vpt
        prim_to_wavspd_(
                &idirection,
                &RC::gamma,
                reinterpret_cast<real64*>(&Vpt),
                &wavespeed.S_min,
                &wavespeed.S_max,
                reinterpret_cast<real64*>(&B0)
        );
        return wavespeed;
    }



    FluidVector PrimToCons(
            int32 iprs_flag,
            real64 ent,
            FluidVector& Vpt,
            GeoVector& B0
    ) {
        FluidVector Upt;
        prim_to_cons_(
                &RC::gamma,
                &RC::gamma_euler,
                &RC::isospeed,
                &RC::barotropic_stiffness,
                &RC::molewt,
                &RC::gasconst,
                &RC::rhofloor,
                &RC::prsfloor,
                &RC::tempfloor,
                &RC::magB_small_cgl,
                &RC::eps_bdry_cgl,
                &iprs_flag,
                &ent,
                reinterpret_cast<real64*>(&B0),
                reinterpret_cast<real64*>(&Upt),
                reinterpret_cast<real64*>(&Vpt)
        );
        return Upt;
    }


    FluidVector ConsToPrim(
            int32 iprs_flag,
            real64 ent,
            FluidVector& Upt,
            GeoVector& B0
    ) {
        FluidVector Vpt;
        cons_to_prim_(
                &RC::gamma,
                &RC::gamma_euler,
                &RC::isospeed,
                &RC::barotropic_stiffness,
                &RC::molewt,
                &RC::gasconst,
                &RC::rhofloor,
                &RC::prsfloor,
                &RC::tempfloor,
                &RC::magB_small_cgl,
                &RC::eps_bdry_cgl,
                &iprs_flag,
                &ent,
                reinterpret_cast<real64*>(&B0),
                reinterpret_cast<real64*>(&Upt),
                reinterpret_cast<real64*>(&Vpt)
        );
        return Vpt;
    }

    /*
     * This subroutine forces to p_perp and p_par to remain in the physical domain.
     * Only applicable to CGL system.
     * Inputs: pi, Bsq, p_par, p_perp
     * Outputs: updated_p_par, updated_p_perp
     *
     */
    PParPPerp_t ForceToHyperbolicDomain(
        real64 Bsq,
        GeoVector p_par,
        GeoVector p_perp
    ) {
        PParPPerp_t out;
        force_pperp_ppar_to_hyperbolic_domain_(
            RC::pi,
            RC::magB_small_cgl,
            RC::eps_bdry_cgl,
            RC::prsfloor,
            Bsq,
            reinterpret_cast<real64*>(&out[0]),
            reinterpret_cast<real64*>(&out[1])
        );
        return out;
    }






};


