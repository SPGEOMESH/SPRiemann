
#pragma once

#include <Tbasis.h>

#include <spriemann/patch/patch.h>

/*
 *
 * Provide a global view of a patch (in storage)
 * via GlobalView::make().
 */
template <typename Parameters>
class GlobalView {

   using P = Parameters;
   using Patch = Patch<P>;

public:

   // Store a single cell-centered variable
   // with index `icc` in a set of files.
   // Only the zones with zone (p,t,r) r index `ir`
   // and patch z index `iz` will be printed.
   // Recall (x,y,z) denotes patch index location
   // while (p,t,r) or (phi, theta, r)
   // denotes zone index location.
   // Should be called by a process. If all
   // processes execute, all patches will be stored.
   static void make(
         Patch* patches,
         int32 n_patch,
         const std::string& prefix = "",
         int32 iproc = 0,
         bool distinguish_proc = false,
         int32 icc = 0,
         int32 iz = 0,
         int32 ir = P::NGhost,
         int32 ilevel = 0
   ) {
      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto& patch = patches[ii];
         auto isector = patch.index.isector;
         auto ix = patch.index.ix;
         auto iy = patch.index.iy;
         std::string filename;
         if (distinguish_proc) {
            filename = std::format("proc{}.{}icc{}r{}z{}l{}isector{}.x{}y{}.dat", iproc, prefix, icc, ir, iz, ilevel, isector, ix, iy);
         }
         else {
            filename = std::format("{}icc{}r{}z{}l{}isector{}.x{}y{}.dat", prefix, icc, ir, iz, ilevel, isector, ix, iy);
         }
         std::string content;
         std::string tricorner = "3.333333333333333e+33 ";
         for (int32 ip = P::Ntp_zones-1; ip >= 0; --ip) {
            for (size_t it = 0; it < P::Ntp_zones; ++it) {
               // check for tricorners.
               if (P::is_tricorner_zone(ix, iy, it, ip)) {
                  content += tricorner;
               }
               else {
                  size_t icc_ = static_cast<size_t>(icc);
                  size_t ir_ = static_cast<size_t>(ir);
                  auto ip_ = static_cast<size_t>(ip);
                  content += std::format("{:.15e} ", patch.U[{icc_, ir_}][{ip_, it}].under());
               }
            }
            content += "\n";
         }
         // for testing: this dumps a raw image of the patch U
//         content += patch.U.str();
         savetxt(filename, content);
      }
   }


   [[nodiscard]]
   static std::string make_3d(
         Patch& patch
   ) {
      auto ix = patch.index.ix;
      auto iy = patch.index.iy;
      std::string content;
      std::string tricorner = "3.333333333333333e+33 ";
      for (int32 ir = 0; ir < P::Nr_zones; ++ir) {
         content += std::format("// ir = {}\n", ir);
         for (int32 ip = P::Ntp_zones-1; ip >= 0; --ip) {
            for (size_t it = 0; it < P::Ntp_zones; ++it) {
               // check for tricorners.
               if (P::is_tricorner_zone(ix, iy, it, ip)) {
                  content += tricorner;
               }
               else {
                  auto ir_ = static_cast<size_t>(ir);
                  auto ip_ = static_cast<size_t>(ip);
                  content += std::format("{:.15e} ", patch.U[{0, ir_}][{ip_, it}].under());
               }
            }
            content += "\n";
         }
      }
      return content;
   }

};



