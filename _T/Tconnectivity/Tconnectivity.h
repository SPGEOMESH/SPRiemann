
#pragma once

#include <Tbasis.h>

#include <spriemann/constants/parameters.h>
#include <spriemann/constants/directives.h>
#include <spriemann/connectivity.h>


class TConnectivity {
public:

   using Directives = Directives<>;

   static constexpr size_t ntp_patches = 4;
   static constexpr size_t nr_patches = ntp_patches;

   static constexpr size_t ntp_zones = 1;
   static constexpr size_t nr_zones = 1;
   static constexpr size_t NGhost = 0;
   static constexpr size_t ilevel_max = 1;

//   static constexpr size_t ntp_global = ntp_patches*ntp_zones;
//   static constexpr size_t nr_global = nr_patches*nr_zones;
//   static constexpr size_t NGhost_global = NGhost;

   using P = Parameters<
         Directives,
         ntp_zones,
         nr_zones,
         ntp_patches,
         nr_patches,
         NGhost,
         ilevel_max
   >;


   static void Tdryrun() {

      auto conn = Connectivity<P>();

      conn.initialize();

      cout << conn.str() << endl;

   }


   static void print_info(int32 ix, int32 iy) {
   }

   static void Tvisualize_xyz() {

      auto conn = Connectivity<P>();

      auto pid = conn.initialize();

      std::vector<int32> isectors{
         0,
         1,
         2,
         3,
         4,
         5,
      };

      // We print the C indices
      for (int32 ilevel = 0; ilevel < ilevel_max; ++ilevel) {
         cout << std::format("=-=-=-=-=-=-\nlevel {}:\n=-=-=-=-=-=-", ilevel) << endl;
         for (auto isector: isectors) {
            cout << std::format("-----------\nsector {}\n-----------", isector) << endl;
            for (int32 ix = 0; ix < ntp_patches; ++ix) {
               for (int32 iy = 0; iy < ntp_patches; ++iy) {
                  for (int32 iz = 0; iz < nr_patches; ++iz) {
                     auto node = conn.get({isector, ilevel, ix, iy, iz});
                     cout << "Using idx " << conn.get_idx({isector, ilevel, ix, iy, iz}) << endl;
                     cout << ix << ":" << iy << ":" << iz << ":\t\t(" << node.U->index.ix << ", " << node.U->index.iy << ", " << node.U->index.iz << ")  [" << node.T->index.ix << ", " << node.T->index.iy << ", " << node.T->index.iz << "]" << endl;
                     cout << "\t(" << node.L->index.ix << ", " << node.L->index.iy << ", " << node.L->index.iz << ")";
                     cout << " (" << node.index.ix << ", " << node.index.iy << ", " << node.index.iz << ")";
                     cout << " (" << node.R->index.ix << ", " << node.R->index.iy << ", " << node.R->index.iz << ")" << endl;
                     cout << "\t[" << node.B->index.ix << ", " << node.B->index.iy << ", " << node.B->index.iz << "]  (" << node.D->index.ix << ", " << node.D->index.iy << ", " << node.D->index.iz << ")" << endl;
                  }
               }
            }
         }
      }

   }

   static void Tvisualize_Fortran() {

      auto conn = Connectivity<P>();

      auto pid = conn.initialize();

      std::vector<int32> isectors{
            0,
            1,
            2,
            3,
            4,
            5,
      };

      // We print the Fortran indices
      for (int32 ilevel = 0; ilevel < ilevel_max; ++ilevel) {
         cout << std::format("=-=-=-=-=-=-\nlevel {}:\n=-=-=-=-=-=-", ilevel+1) << endl;
         for (auto isector: isectors) {
            cout << std::format("-----------\nsector {}\n-----------", isector+1) << endl;
            for (int32 ix = 0; ix < ntp_patches; ++ix) {
               for (int32 iy = 0; iy < ntp_patches; ++iy) {
                  for (int32 iz = 1; iz < 2; ++iz) {
                     auto node = conn.get({isector, ilevel, ix, iy, iz});
                     cout << "\t\t\t(" << node.U->index.ix+1 << ", " << node.U->index.iy+1 << ")" << endl;
                     cout << "\t(" << node.L->index.ix+1 << ", " << node.L->index.iy+1 << ")";
                     cout << " (" << node.index.ix+1 << ", " << node.index.iy+1 << ")";
                     cout << " (" << node.R->index.ix+1 << ", " << node.R->index.iy+1 << ")" << endl;
                     cout << "\t\t\t(" << node.D->index.ix+1 << ", " << node.D->index.iy+1 << ")" << endl;
                     cout << endl;
                  }
               }
            }
         }
      }

   }



};



