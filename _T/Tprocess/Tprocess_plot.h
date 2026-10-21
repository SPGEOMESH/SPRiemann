
#pragma once

#include <Tbasis.h>

#include <spriemann/constants/parameters.h>
#include <spriemann/constants/directives.h>
#include <spriemann/process.h>

#include <patch_global_view.h>

#include <cmath>

/*
 * First-pass tests and sanity checks for Process,
 * that produce data for post-processing and visualization.
 *
 */
class TProcess_Plot {
public:





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
   static void Texchange_plot() {
      using Directives = Directives<>;

      static constexpr size_t ntp_patches = 3;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones = 40;
      static constexpr size_t nr_zones = 10;
      static constexpr size_t NGhost = 10;
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
         auto& patch = process.patches[ii];
         auto index = patch.index;
         auto si = index.isector+1;
         for (size_t k = P::NGhost; k < P::NGhost+P::nr_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost+P::ntp_zones; ++j) {
               for (size_t i = P::NGhost; i < P::NGhost+P::ntp_zones; ++i) {
                  double754 x1, x2, X, y, umin, umax;
                  if (i < ntp_zones / 2 + P::NGhost) {
                     // X period: 2 vertically
                     X = static_cast<R64>(i - P::NGhost) / R64(ntp_zones) * R64(4 * M_PI);
                  }
                  else {
                     // X period: 8 vertically
                     X = static_cast<R64>(i - P::NGhost) / R64(ntp_zones) * R64(16 * M_PI);
                  }
                  // x1 min 0 max 1
                  x1 = R64(0.5)*(R64(1) - X.cos());
                  // x2 min 0 max 0.25
                  x2 = R64(0.125)*(R64(1) - X.cos());
                  // y min 0 max 1
                  y = static_cast<R64>(j - P::NGhost) / R64(ntp_zones);
                  // set scale in sector-dependent way.
                  // The *total* u variation is from 0.1 to
                  if (si == 1) {
                     umin = R64(1.2);
                     umax = R64(2);
                  }
                  else if (si == 2) {
                     umin = R64(4.0);
                     umax = R64(5.6);
                  }
                  else if (si == 3) {
                     umin = R64(6.2);
                     umax = R64(12.5);
                  }
                  else if (si == 4) {
                     umin = R64(7.8);
                     umax = R64(9);
                  }
                  else if (si == 5) {
                     umin = R64(9.4);
                     umax = R64(11.1);
                  }
                  else if (si == 6) {
                     umin = R64(11);
                     umax = R64(12);
                  }
                  // min umin max umax
                  process.patches[ii].U[{0,k}][{j, i}] = umin + (umax-umin)*y*x1 + (umax-umin)*(R64(1.0) - y)*x2;
               }
            }
         }
      }

      GlobalView::make(process.patches.get(), process.n_patch, "t1");

      process.halo_exchange();

      GlobalView::make(process.patches.get(), process.n_patch, "t2");

   }





   // This illustrates a bug that arose in the Texchange_plot test.
   // Illustrates how easy it is to make a mistake like this from misuse of `auto`,
   // and end up with a difficult bug.
   static void Texchange_plot_clinic() {
      using Directives = Directives<>;

      static constexpr size_t ntp_patches = 3;
      static constexpr size_t nr_patches = 1;

      static constexpr size_t ntp_zones1 = 3;
      static constexpr size_t ntp_zones2 = 30;
      static constexpr size_t nr_zones = 2;
      static constexpr size_t NGhost = 2;
      static constexpr size_t ilevel_max = 1;

      // create parameters
      using P1= Parameters<
            Directives,
            ntp_zones1,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;
      using P2= Parameters<
            Directives,
            ntp_zones2,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
      >;

      using Process1 = Process<P1>;
      using Process2 = Process<P2>;

      cout << "Process does not scale with the number of zones per patch." << endl;
      cout << "Patch does scale with the number of zones per patch." << endl;
      cout << std::format("Process {} Connectivity {}", sizeof(Process1), sizeof(Process1::Connectivity)) << endl;
      cout << std::format("Process {} Connectivity {}", sizeof(Process2), sizeof(Process2::Connectivity)) << endl;
      cout << std::format("Patch {}\nPatch {}", sizeof(Process1::Patch), sizeof(Process2::Patch)) << endl;

      auto process2 = Process2();

      // This seemingly innocuous routine can cause a SIGSEGV (segmentation fault) due to stack overflow.
      // In a recent test the `auto patch` declaration in loop [1] pushes 46 MB on the stack, exceeding the typical limit of 8 MB.
      // An amazing effect of simply creating what appears like a reference (if you've been writing too much Python or Fortran recently),
      // in order to read a single integer (the isector).
      // The bug is made more diabolical because (in a recent test) the stack pointer is pushed during the function's prologue,
      // redirecting away from the loop where the issue is hiding.
      // The revised loop [2] raises the stack by a pointer size, instead of 46 MB.
#if 0
      // [1]
      for (int32 ii = 0; ii < process2.n_patch; ++ii) {
         auto patch = process2.patches[ii];
         auto is = patch.index.isector+1;
         for (size_t k = P2::NGhost; k < P2::NGhost+P2::nr_zones; ++k) {
            for (size_t j = P2::NGhost; j < P2::NGhost+P2::ntp_zones; ++j) {
               for (size_t i = P2::NGhost; i < P2::NGhost+P2::ntp_zones; ++i) {
                  auto x = static_cast<R64>(i-P2::NGhost)/R64(ntp_zones2)*R64(M_PI);
                  auto y = static_cast<R64>(j-P2::NGhost)/R64(ntp_zones2)*R64(2*M_PI);
                  process2.patches[ii].U[{0,k}][{j, i}] = R64(is) + y.cos()*x.cos();
               }
            }
         }
      }
#else
      // [2]
      for (int32 ii = 0; ii < process2.n_patch; ++ii) {
         auto& patch = process2.patches[ii];
         auto is = patch.index.isector+1;
         for (size_t k = P2::NGhost; k < P2::NGhost+P2::nr_zones; ++k) {
            for (size_t j = P2::NGhost; j < P2::NGhost+P2::ntp_zones; ++j) {
               for (size_t i = P2::NGhost; i < P2::NGhost+P2::ntp_zones; ++i) {
                  auto x = static_cast<R64>(i-P2::NGhost)/R64(ntp_zones2)*R64(M_PI);
                  auto y = static_cast<R64>(j-P2::NGhost)/R64(ntp_zones2)*R64(2*M_PI);
                  process2.patches[ii].U[{0,k}][{j, i}] = R64(is) + y.cos()*x.cos();
               }
            }
         }
      }
#endif


   }









};



