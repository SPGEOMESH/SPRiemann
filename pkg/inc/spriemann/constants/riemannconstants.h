
#pragma once

#include <csriemann/constants/types.h>





/*
 * TODO these aren't really constants, but
 *   later can allow in more flexibility.
 *
 */
template <
      typename Parameters
>
class RiemannConstants {
public: // types:

   static constexpr real64 pi = M_PI;

   // Polytropic index.
   static constexpr real64 gamma = 5.0/3.0;
   // Polytropic index for the gas.
   static constexpr real64 gamma_euler = 5.0/3.0;
   // Isothermal sound speed; only for isomhd problems.
   static constexpr real64 isospeed = 1.0;
   // Stiffens the isothermal EOS if non-zero.
   static constexpr real64 barotropic_stiffness = 0.0;
   // Set up mole wt. of the gas.
   static constexpr real64 molewt = 1.4;
   // Set up scaled gas constant.
   static constexpr real64 gasconst = 1.0;
   static constexpr real64 gasconst_euler = 1.0;

   // density floor.
   static constexpr real64 rhofloor = 1.0e-12;
   // pressure floor.
   static constexpr real64 prsfloor = 1.0e-9;

   // used only when the cooling fn is called.
   static constexpr real64 tempfloor = 1.0e-5;


   // ------
   // CGL
   // ------

   // CGL relaxation time
   static constexpr real64 source_tau = 1e-10;

   // set p_par = p_perp if magB < magB_small
   static constexpr real64 magB_small_cgl = 1.0e-10;
   // p_par and p_perp should be epsilon distance away
   static constexpr real64 eps_bdry_cgl = 0.05;
   // entropygr ( ixx, iyy, izz, 1)
   // or 0.0
   static constexpr real64 ent = 0.0;



   // ------------------------
   // Problem-specific
   // ------------------------

/*
TODO - sin and cos fields can be initialized as plain static

! Note: The following numbers to be changed for suitable problems.
! If "Mstar" is non-zero then gravity will be switched on in the ADER subroutine.
! If "omega_passive" is non-zero then non-inertial terms will be switched
! on in the ADER subroutine.

      newtonsc = 6.674E-08        ! Newton's constant; in cgs unit
      Mstar = 0.0                 ! Mass of Sun; in grams  =  1.989E33
      omega_passive ( :) = 0.0    ! Rotation axis of star.
      omega_passive ( 3) = 0.0

      omega_active ( :) = 0.0

! ----------
! For multipole simulations.

      iflag_npoles ( 1:4) = 0
      iflag_npoles ( 1)   = 1     ! Enable only dipole by default.

      multipole_Bmag ( 1:4) = 0.0
      multipole_Bmag ( 1)   = 0.0

! ----------
! Define the dipole and qudrapole magnetic field strengths.

      a = pi / 6.0

! monopole.
      iflag_monopole  = 0    ! 0 - no monopole field, 1 - monpole, 2 - split-monopole.
      monopole_u0     = 0.0  ! radial velocity, used in parker spiral solution
      monopole_B ( :) = 0.0

! dipole and quadrupole.
      dipole_B ( 1) = 0.0 * SIN ( a)
      dipole_B ( 2) = 0.0
      dipole_B ( 3) = 0.0 * COS ( a)

      qpole_B ( 1:3, 1:3) = 0.0


! ----------
! Assign the unit vector for the inflow velocity at the outer boundary.
! It is used to correctly set the boundary condition at the outer boundary.
! By default it is assigned as zero.

      inflow_uvec ( 1:3) = 0.0



*/


};


