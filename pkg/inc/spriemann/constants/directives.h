
#pragma once


#include <spriemann/constants/types.h>


struct Empty {};

template <
         IPDESystem IPDE_SYSTEM_ = IPDESystem::CGL,
         int32 NUMORDER_ = 3,
         //------
         // for Riemann_Integrate_CS
         int32 do_loop_for_matmul_ = 0,
         bool debug_geom_ = false,
         bool debug_communication_ = false,
         //------
         // iso2
         int32 NAREA_QUAD_PTS_ = 5,
         int32 NBFLD_QUAD_PTS_ = 5,
         int32 nquick_vol_quad_pts_ = 6,
         //-----
         // radial
         int32 radial_mesh_ = 0, // 0,1,...,9
         //-------
         // etc......
         // Choose whether we want adiabatic or isothermal versions of subroutines. They
         // have slightly different subroutine interfaces. Isothermal case should be run
         // without pressure positivity switch.
         int32 adiabatic_ = 1, // 1


         // $smallnum = 10^{smallnum_power}$
         int32 smallnum_power = -12,
         // todo etc .....
         //------
         // application options
         typename Options_ = Empty
>
class Directives {
public: // constant fields and types:

   // application options (for Helio, CAK, ...)
   using Options = Options_;

   static constexpr IPDESystem IPDE_SYSTEM = IPDE_SYSTEM_;
   static constexpr int32 NUMORDER = NUMORDER_;
   static constexpr int32 nRK_stages = (NUMORDER == 2)? 2: (NUMORDER == 3)? 3: 5;
   static constexpr int32 do_loop_for_matmul = do_loop_for_matmul_;
   static constexpr bool debug_geom = debug_geom_;
   static constexpr bool debug_communication = debug_communication_;

   static constexpr int32 NAREA_QUAD_PTS = NAREA_QUAD_PTS_;
   static constexpr int32 NBFLD_QUAD_PTS = NBFLD_QUAD_PTS_;
   static constexpr int32 nquick_vol_quad_pts = nquick_vol_quad_pts_;
   static constexpr int32 radial_mesh = radial_mesh_;
   static constexpr int32 adiabatic = adiabatic_;

   static constexpr int32 NVAR =
      (IPDE_SYSTEM==IPDESystem::MHD)? 8:
      (IPDE_SYSTEM==IPDESystem::CGL)? 11:
      /*IPDE_SYSTEM==IPDESystem::CGL_p*/ 12;
   static constexpr int32 NFLUID = 0;

   static constexpr bool USE_COOLING_FUNCTION = false;
   static constexpr bool USE_B0B1_FORMULATION = false;


   // For problems with very strong shocks and large magnetic fields, one
   // may have a loss of pressure positivity. The entropy formulation cures
   // that and it is invoked by switching on the pressure positivity switch.
   // It does cause a local loss of conservation in the zone where it is invoked.
   // It is switched on by setting "PRESSURE_POSITIVITY" to 1 or 2 (0, default).
   // When "PRESSURE_POSITIVITY == 1" we only check for positivity of density or
   // pressure. When "PRESSURE_POSITIVITY == 2" we check all the variables.
   static constexpr int32 PRESSURE_POSITIVITY = 0;

   // Whether to keep extra data for silo I/O.
   static constexpr bool keep_silo_extra_data = false;

   // Number of extra variables, written out in the silo file. This helps in
   // monitoring the simulation. "1" write out divergence and other extra variables
   // as needed (see below). "0" write out only the conserved variables.
   static constexpr int32 set_silo_extra_data() {
      int32 out = 0;
      if constexpr (keep_silo_extra_data) {
         if constexpr (PRESSURE_POSITIVITY >= 1) {
            // Write out
            // 1) divergence
            // 2) entropy
            // 3) entropy-use-flag
            // 4) tempfloor-use-flag (*)
            // (*) only if using cooling function.
            out = (USE_COOLING_FUNCTION)? 4: 3;
         }
         else {
            if constexpr (USE_B0B1_FORMULATION) {
               // the silo read/write variable is sized
               // as "n_cc_components + B0 ( 3) + entropy".
               out = 4;
            }
            else {
               // the silo read/write variable is sized
               // as "n_cc_components + entropy".
               out = 1;
            }
         }
      }
      return out;
   }
   static constexpr int32 SILO_EXTRA_DATA = set_silo_extra_data();

   // cell-centered components
   static constexpr int32 n_cc_components = 8;
   // face-centered modes
   static constexpr int32 n_fc_modes = 1;
   // cell-centered modes
   static constexpr int32 n_cc_modes =
      (NUMORDER == 2)? 4:
      (NUMORDER == 3)? 10:
      /* NUMORDER 4 */ 20;

   // Could be a better way
   static constexpr real64 compute_smallnum() {
      real64 out = 1.0;
      for (int i = 0; i < -smallnum_power; ++i) {
         out /= 10.0;
      }
      return out;
   }
   static constexpr real64 smallnum = compute_smallnum();
};


