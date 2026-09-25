
#pragma once

#include <Tbasis.h>

#include <csriemann/patch/patch.h>

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
   static void make(Patch* patches, int32 n_patch, const std::string& prefix = "", int32 iproc = 0, bool distinguish_proc = false, int32 icc = 0, int32 iz = 0, int32 ir = P::NGhost, int32 ilevel = 0) {
      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto patch = patches[ii];
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
         size_t ir_ = static_cast<size_t>(ir);
         size_t icc_ = static_cast<size_t>(icc);
         for (int32 ip = P::Ntp_zones-1; ip >= 0; --ip) {
            for (size_t it = 0; it < P::Ntp_zones; ++it) {
               auto ip_ = static_cast<size_t>(ip);
               // check for tricorners.
               if ((ix == 0 && iy == 0 && it < P::NGhost && ip < P::NGhost) ||
                     (ix == 0 && iy == P::ntp_patches-1 && it < P::NGhost && ip >= P::NGhost+P::ntp_zones) ||
                     (ix == P::ntp_patches-1 && iy == 0 && it >= P::NGhost+P::ntp_zones && ip < P::NGhost) ||
                     (ix == P::ntp_patches-1 && iy == P::ntp_patches-1 && it >= P::NGhost+P::ntp_zones && ip >= P::NGhost+P::ntp_zones)
               ) {
                  content += tricorner;
               }
               else {
                  content += std::format("{:.15e} ", patch.U[{icc_, ir_}][{ip_, it}].under());
               }
            }
            content += "\n";
         }
//         content += patch.U.str();
         savetxt(filename, content);
      }
   }

};



