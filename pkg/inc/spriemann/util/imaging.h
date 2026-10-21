
#pragma once

#include <string>
#include <fstream>

#include <spriemann/constants/types.h>



/*

 ! These parameters control the checkpoint dumps and imaging.

! Start with initial time set to 0.Restart will refresh this from "riemann.in".
      time = 0.0

! Time at which we want to output the next image.Set this yourself for restart.
      timeimage = 0.0

! Number of timesteps you want to take. Code stops after these many steps.
      ntstep = 5000000

! Write out full restart dumps every "ndumpstep" steps. By default, all the
! variables are saved in 64 bit unformatted files called "dump0001" etc.
      ndumpstep = 500

! Write out images at every "nimgstep" timesteps.(not used)
      nimgstep = 10

! Time at which you want the simulation to stop.
      timestop = 0.20

! Choose the number of silo images that are expected to be written-out.
      nsilo_images = 5


 ! Write out images every "dtimage" amount of time. By default, each variable
! is saved in its own file name in 32 bit formatted form, like "rhoa0001" etc.

      dtimage = timestop / REAL ( nsilo_images)


*/




template <
      typename Parameters
>
class Imaging {

   using P = Parameters;

   using r4gr_revt = P::r4gr_revt;


   int32 iproc, nsilo_images;
   // open("riemann_geomesh.out");

   std::string silofile;
   std::string rhofile, vlxfile, vlyfile, vlzfile, prsfile,
   bfxfile, bfyfile, bfzfile, dumpfile, speciefile;

   // variable used for imaging variables.

   r4gr_revt r4gr;

   std::vector<std::string> silo_allvars;


   // todo see riemann_module_3d.
   // this is indexed by the number of patches.
   // the size depends on B0B1_formulation.
   U_silo_amr_revt U_silo_amr;


   // todo see riemann_module_3d.
   // Extra variables for silo output.
   // divergence and entropy, tempfloor flags.
   U_silo_amr_revt divgr_amr, flag_entprs_amr, flag_Tfloor_amr;



public: // API:






};

