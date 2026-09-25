
#pragma once

#include <csriemann/constants/types.h>
#include <csriemann/process.h>
#include <csriemann/connectivity.h>

#include <kiwi/timer.h>
#include <kiwi/mpi/sim.h>

template <
      typename Parameters
>
class SimulationRK {
public: // types:

   using P = Parameters;
   using Process = Process<P>;

   using bcarr_r_revt = P::bcarr_r_revt;

public: // constant fields:

   static constexpr int nRK_stages = P::Directives::nRK_stages;

public: // fields:

   Timer timer;
   real64 timing_info1;
   real64 timing_info2;
   real64 alph_cak;

   // simulation time
   real64 t;

   // Enforcing radial boundary conditions.
   // Used for PAD_BOUNDARY subroutine.
   // todo FINDHOME
   bcarr_r_revt bcarr_r;

   Process process;

public: // ctor:

   SimulationRK(int argc, char** argv):
      timer(Timer()),
      timing_info1(0),
      timing_info2(0),
      alph_cak(0),
      t(0),
      process(Process(argc, argv))
   {

      // At this point, communicators are established
      // but simulation data needs to be initialized.

      // todo how do I know my geometry args???
      real64 rmin = 1.234;
      real64 rmax = 1.234;
      real64 smallnum = 1.234;
      real64 radial_mesh_ratio = 1.234;

      process.initialize_geometry(
            rmin,
            rmax,
            radial_mesh_ratio
      );

      // Possibly read in a restart dump.
      // todo inputs to say whether to do a restart.
//      readin_restart_dump();

      // end of ctor.
   }


   ~SimulationRK() {
      // deinitialize MPI.
      process.~Process();
   }

public: // API:

   void start(int ntimesteps) {

      // At this point geometry is initialized
      // and patches are configured.
      // If there is a restart, that has already
      // occured. Communication is initialized
      // so the timestep is ready to commence.

      for (int ti = 0; ti < ntimesteps; ++ti) {
         timer.reset();
         // TODO:
//         write_image_if_needed();

         // TODO:
//         write_heliodata_if_needed();

         for (int rki = 0; rki < nRK_stages; ++rki) {

            // TODO:
            // load/fetch
            // load/fetch
            // load/fetch

            process.main_loop(rki);

            // TODO
//            process.reduce_to_get_V_min_V_max();

            process.defect_loop1();

//            load_amr_fluxes_at_bdy();
//            load_local();
            // TODO:
            // load/fetch
            // load/fetch
            // load/fetch

//            unpack_defect_fluxes();

            process.defect_loop2();

         }

         timer.done();
         process.print_message(std::format("elapsed time: {}", timer.result()));

//         simulation_wide_cleanup_after_rk_loop();

      }

//      simulation_wide_cleanup_after_time_loop();

   }

};



