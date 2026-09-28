
#pragma once

#include <Ibasis.h>

#include <csriemann/csriemann.h>


class IAugust30 {


   using Directives = Directives<>;

   static constexpr size_t ntp_patches = 1;
   static constexpr size_t nr_patches = 1;

   static constexpr size_t ntp_zones = 5;
   static constexpr size_t nr_zones = 5;
   static constexpr size_t NGhost = 2;
   static constexpr size_t ilevel_max = 1;

   // create parameters
   using Parameters = Parameters<
         Directives,
         ntp_zones,
         nr_zones,
         ntp_patches,
         nr_patches,
         NGhost,
         ilevel_max
   >;

   using Simulation = SimulationRK<Parameters>;



public:


   static void start() {

      auto simulation = Simulation();


      int32 ntimesteps = 5;

      // todo these belong here;
      //  decide how they are passed thru
      real64 rmin = 1.5;
      real64 rmax = 2.5;
      // If the ratio = 1.0 then there is
      // no ratioing of the radial mesh,
      // when the geometric mesh is chosen.
      real64 radial_mesh_ratio = 2.0;






//      ! Courant number
//            cfl_coef = 0.2
//
//      ! If the divergence of the velocity at a zone center exceeds
//      ! this coef times the sound speed then we have to flatten the slopes.
//            flatten_coef = 0.4
//
//      ! This is the factor by which the pressure on one or the other side
//      ! of a zone has to be larger. If this condition is met, we look to
//      ! that side to increase the flattener in the zone being considered.
//            prs_flatten_coef = 1.5
//
//      fracdiff = 0.1   ! Set fractional change needed for entropy fix.(deprecated)
//
//      ! If the divergence of the velocity at a zone boundary exceeds this coef times
//      ! the sound speed then we have to blend in some fraction of HLL RS.
//            divratio_einfeldt_rs = 0.4
//
//      ! Set a threshold value that triggers the Lapidus viscosity as well as
//      ! its value. Set "lapidus_coef" between [0.0, 0.2]. Default is 0.0.
//      lapidus_threshold = 0.2
//      lapidus_coef = 0.05
//
//      ! Using 100% multid flux can sometimes mean we take flux contributions from
//      ! a stencil that does not depend on the zones that contribute to that face.
//      ! For that reason, we blend multid RS with 1d RS associated with that face.
//      ! Blending in more multidimensionality also seems to help out with pressure
//      ! positivity without any artificial fixes.
//      ! This value should range between [0.4, 0.6] for good performance.
//
//      ! "blend_multid_flux = 1.0 / 6.0" means we use Simpson rule with 1/6 of the
//      ! flux from each edge.
//      ! "blend_multid_flux = 0.0" means we use full flux from suitable 1d RS.
//      ! The electric field is always fully evaluated by the multid. RS.
//
//            blend_multid_flux = 0.1666
//
//      ! In the fully subsonic case (just) the multi-d HLLC RS can oversteepen
//      ! (almost or fully) stationary contacts. To prevent that we blend in a small
//      ! fraction of the multi-d HLL RS. Ranging between [0.1, 0.2]
//
//      ! While the flattener also plays a role; "blend_mdhll_to_mdhllc = 1.0"
//      ! converts HLLC into a pure HLL. "blend_mdhll_to_mdhllc = 0.0" does the
//      ! reverse in smooth flow.
//
//            blend_mdhll_to_mdhllc = 0.0
//



      real64 cfl_coef = 0.2;

      // todo these seem to me to be separate I files.
//      ! Decide which physical problem you want to initialize with "iproblem".
//      ! iproblem = 0  for testing the MPI messaging or Regular testing
//      ! iproblem = 2  MHD Blast;
//      ! iproblem = 7  Parker wind test;
//      ! iproblem = 9  SOD test;                iproblem = 10 Brio-Wu test
//      ! iproblem = 13 VF slow solar wind;      iproblem = 16 CAK wind
//      ! iproblem = 14 Heliosphere problem;
//
//      iproblem = 14


      simulation.start(ntimesteps);

   }

};

