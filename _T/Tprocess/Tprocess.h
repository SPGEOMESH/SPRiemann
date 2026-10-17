
#pragma once

#include <Tbasis.h>

#include <csriemann/constants/parameters.h>
#include <csriemann/constants/directives.h>
#include <csriemann/process.h>

#include <patch_global_view.h>

#include <cmath>


/*
 * First-pass tests and sanity checks for Process,
 * that produce only standard out messages or small text files.
 *
 */
class TProcess {
public:



   static void Tdryrun(int argc, char** argv) {
      using Directives = Directives<>;

      static constexpr size_t ntp_patches = 4;
      static constexpr size_t nr_patches = 4;

      static constexpr size_t ntp_zones = 2;
      static constexpr size_t nr_zones = 2;
      static constexpr size_t NGhost = 2;
      static constexpr size_t ilevel_max = 2;

      // create parameters
      using P= Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process = Process<P>;

      auto process = Process();

      cout << process.str() << endl;

   }



   static void Tconnectivityview() {
      using Directives = Directives<>;

      static constexpr size_t ntp_patches = 4;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones = 2;
      static constexpr size_t nr_zones = 2;
      static constexpr size_t NGhost = 2;
      static constexpr size_t ilevel_max = 1;

      // create parameters
      using P= Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process = Process<P>;

      auto process = Process();

      std::map<int32, std::string> symb = {
            {-2,"🚫"},
            {0,"🔵"},
            {1,"🔴"},
            {2, "🟢"},
            {3, "⚪"},
            {4, "🟡"},
            {5, "🟣"},
            {6, "🟤"},
            {7, "🟠"},
            {8, "🔷"},
            {9, "🔶"},
      };

      for (int32 isector = 0; isector < 6; ++isector) {
         for (int32 iz = 0; iz < P::nr_patches; ++iz) {
            std::string filename = std::format("iproc{}conn.isector{}z{}.txt", process.world_rank, isector, iz);
            std::string content;
            for (int32 iy = P::ntp_patches-1; iy >= 0; --iy) {
               for (int32 ix = 0; ix < P::ntp_patches; ++ix) {
                  auto node = process.connectivity.get({isector, 0, ix, iy, iz});
                  auto bot = process.connectivity.get(node.B);
                  auto up = process.connectivity.get(node.U);
                  auto Bi = bot.world_rank;
                  auto Ui = up.world_rank;
                  content += std::format("{}{} {}{} 0⚫     ", (Bi < 0)? 0: Bi, symb[Bi], Ui, symb[Ui]);
               }
               content += "\n";
               for (int32 ix = 0; ix < P::ntp_patches; ++ix) {
                  auto node = process.connectivity.get({isector, 0, ix, iy, iz});
                  auto left = process.connectivity.get(node.L);
                  auto right = process.connectivity.get(node.R);
                  auto Li = left.world_rank;
                  auto Ci = node.world_rank;
                  auto Ri = right.world_rank;
                  content += std::format("{}{} {}{} {}{}     ", Li, symb[Li], Ci, symb[Ci], Ri, symb[Ri]);
                  content += std::format("");
               }
               content += "\n";
               for (int32 ix = 0; ix < P::ntp_patches; ++ix) {
                  auto node = process.connectivity.get({isector, 0, ix, iy, iz});
                  auto dn = process.connectivity.get(node.D);
                  auto top = process.connectivity.get(node.T);
                  auto Di = dn.world_rank;
                  auto Ti = top.world_rank;
                  content += std::format("0⚫ {}{} {}{}     ", Di, symb[Di], (Ti < 0)? 0: Ti, symb[Ti]);
                  content += std::format("");
               }
               content += "\n\n";
            }
            savetxt(filename, content);
         }
      }
   }


   // This is a trivial test of GlobalView,
   // on a very small case. Use it as a low-cost
   // experiment demonstrating functionality of GlobalView::make.
   // Note: GlobalView prints tricorner zones as 3.333333333333333e+33.
   static void Tglobal_view() {
      using Directives = Directives<>;

      static constexpr size_t ntp_patches = 1;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones = 1;
      static constexpr size_t nr_zones = 1;
      static constexpr size_t NGhost = 1;
      static constexpr size_t ilevel_max = 1;

      // create parameters
      using P= Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process = Process<P>;
      using GlobalView = GlobalView<P>;

      auto process = Process();

      for (int32 ii = 0; ii < process.n_patch; ++ii) {
         auto dummy = static_cast<R64>(777.0);
//         auto dummy = static_cast<R64>(process.patches[ii].index.isector);
         for (size_t k = P::NGhost; k < P::NGhost+P::ntp_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost+P::ntp_zones; ++j) {
               for (size_t i = P::NGhost; i < P::NGhost+P::nr_zones; ++i) {
                  process.patches[ii].U[{0, i}][{j, k}] = dummy;
               }
            }
         }
      }

//      GlobalView::make(process.patches.get(), process.n_patch, "", process.world_rank, true);
      GlobalView::make(process.patches.get(), process.n_patch, "");

   }



   // A minimal communication test. Deploy with 6 processes.
   // Sectors numbered 1-6 (Fortran numbering) are printed
   // with data coded by sector, as (sector number)*100.
   static void Texchange_1patch() {
      using Directives = Directives<
            IPDESystem::CGL,
            3,
            0,
            false,
            // debug communication
            true
      >;

      static constexpr size_t ntp_patches = 1;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones = 1;
      static constexpr size_t nr_zones = 1;
      static constexpr size_t NGhost = 1;
      static constexpr size_t ilevel_max = 1;

      // create parameters
      using P= Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process = Process<P>;
      using GlobalView = GlobalView<P>;

      auto process = Process();

      // After initialization, a process has patch memory.

      // > Initialize dummy simulation data.
      for (int32 ii = 0; ii < process.n_patch; ++ii) {
         auto x = static_cast<R64>((process.patches[ii].index.isector+1)*100);
         for (size_t k = P::NGhost; k < P::NGhost+P::nr_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost+P::ntp_zones; ++j) {
               for (size_t i = P::NGhost; i < P::NGhost+P::ntp_zones; ++i) {
                  process.patches[ii].U[{0,k}][{j, i}] = x;
               }
            }
         }
      }

      GlobalView::make(process.patches.get(), process.n_patch, "t1");

      process.halo_exchange();

      GlobalView::make(process.patches.get(), process.n_patch, "t2");

   }



   // A minimal communication test that allows observation of
   // intrasector halo exchange (the exchanges performed by MPI_Put).
   // Deploy with 6 or more processes.
   // Sectors numbered 1-6 (Fortran numbering) are printed
   // with data coded by sector, as (sector number)*100.
   static void Texchange_9patch() {
      using Directives = Directives<
            IPDESystem::CGL,
            3,
            0,
            false,
            // debug communication
            true
      >;

      static constexpr size_t ntp_patches = 3;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones = 1;
      static constexpr size_t nr_zones = 1;
      static constexpr size_t NGhost = 1;
      static constexpr size_t ilevel_max = 1;

      // create parameters
      using P= Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process = Process<P>;
      using GlobalView = GlobalView<P>;

      auto process = Process();

      // After initialization, a process has patch memory.

      // > Initialize dummy simulation data.
      for (int32 ii = 0; ii < process.n_patch; ++ii) {
         auto x = static_cast<R64>((process.patches[ii].index.isector+1)*100);
         for (size_t k = P::NGhost; k < P::NGhost+P::nr_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost+P::ntp_zones; ++j) {
               for (size_t i = P::NGhost; i < P::NGhost+P::ntp_zones; ++i) {
                  process.patches[ii].U[{0,k}][{j, i}] = x;
               }
            }
         }
      }

      GlobalView::make(process.patches.get(), process.n_patch, "t1");

      process.halo_exchange();

      GlobalView::make(process.patches.get(), process.n_patch, "t2");

   }




   // A minimal communication test that checks orientation
   // of ghost zones, intended to serve as the next check
   // after the previous tests establish a foundation.
   // Each sector is divided into 9 patches with four zones in each patch,
   // and 2 ghost zones. This is the minimal case where a depth of 2
   // in the ghost zones is possible.
   //
   // The output is an array of numbers printed like
   //  4.921900900000000e+07
   // in this figure, the first digit is the sector number (4),
   // and two coordinates are printed, in this case, (2,1) and (0,0).
   // These two coordinates are divided by the digit `9`.
   // The rest of the digits are all `0` and are meaningless.
   // This allows the tester to review the orientation in the ghost
   // zones and confirm its orientation/direction is correct,
   // in addition to spotting other kinds of issues.
   // We print this way because the printing routine `GlobalView::make()`
   // is designed to print continuous data when the
   // scale grows larger. It is better to reuse this code
   // than to write, revise, and debug more testing routines.
   static void Texchange_2ghost() {
      using Directives = Directives<
            IPDESystem::CGL,
            3,
            0,
            false,
//             debug communication
            true
      >;

      static constexpr size_t ntp_patches = 3;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones = 2;
      static constexpr size_t nr_zones = 2;
      static constexpr size_t NGhost = 2;
      static constexpr size_t ilevel_max = 1;

      // create parameters
      using P= Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process = Process<P>;
      using GlobalView = GlobalView<P>;

      auto process = Process();

      // After initialization, a process has patch memory.

      // > Initialize dummy simulation data.
      for (int32 ii = 0; ii < process.n_patch; ++ii) {
         auto patch = process.patches[ii];
         auto index = patch.index;
         auto is = index.isector+1;
         for (size_t k = P::NGhost; k < P::NGhost+P::nr_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost+P::ntp_zones; ++j) {
               for (size_t i = P::NGhost; i < P::NGhost+P::ntp_zones; ++i) {
                  process.patches[ii].U[{0,k}][{j, i}] = static_cast<R64>(is*10000000 + 9.0*1000000 + index.ix*100000 + index.iy*10000 + 9.0*1000 + (i-P::NGhost)*100 + (j-P::NGhost)*10 + 9.0*1);
               }
            }
         }
      }

      GlobalView::make(process.patches.get(), process.n_patch, "t1");

      process.halo_exchange();

      GlobalView::make(process.patches.get(), process.n_patch, "t2");

   }


};



