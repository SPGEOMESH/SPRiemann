
#pragma once

#include <spriemann/extern_c/com.h>
#include <spriemann/extern_c/mesh.h>
#include <spriemann/constants/types.h>


template <
   typename Parameters
>
class Mesh {
public: // types:

    using P = Parameters;

    using coordinates_gr_global_revt = P::coordinates_gr_global_revt;
    using rb_save_global_revt = P::rb_save_global_revt;
    using dr_save_global_revt = P::dr_save_global_revt;
    using rc_save_global_revt = P::rc_save_global_revt;

    using vert_coords_patch_revt = P::vert_coords_patch_revt;
    using imask_active_patch_revt = P::imask_active_patch_revt;

    using coordinates_gr_revt = P::coordinates_gr_revt;
    using vertices_revt = P::vertices_revt;

    using cc_mode_revt = P::cc_mode_revt;


public: // fields:

    // avoid ODR-use of constexpr
    // static values defined in P.
    int32 ntp_global = P::ntp_global;
    int32 nr_global = P::nr_global;
    int32 NGhost_global = P::NGhost_global;
    int32 ntp_zones = P::ntp_zones;
    int32 nr_zones = P::nr_zones;
    int32 ntp_patches = P::ntp_patches;
    int32 nr_patches = P::nr_patches;
    int32 NGhost = P::NGhost;
    int32 ilevel_max = P::ilevel_max;

    int32 n_cc_modes = P::n_cc_modes;
    // NB: this is nvar+nfluid (Fields size)
    int32 n_cc_components = P::n_cc_components;
    int32 n_fc_modes = P::n_fc_modes;

public: // fields:

    // todo review: these structured types
    //  are byte-equivalent to vertices_revt and coordinates_gr_global
    // todo review Fortran - (4, 3) making 4 coordinates for each of 6 sectors - review
    Vector<Coordinate, 6> vertices;
    Vector<coordinates_gr_revt, 6> coordinates;
    // todo for ilevel_max > 1 - 5D vector with ilevel argument - ntp_global doubles in dimension!!!

   real64 rmin;
   real64 rmax;

    // todo a coordinates_gr_global for each ilevel (CHECK)
    //    ...very inefficient (?)

    // radial mesh
    rb_save_global_revt rb;
    dr_save_global_revt dr;
    rc_save_global_revt rc;


    cc_mode_revt ii_mode;
    cc_mode_revt jj_mode;
    cc_mode_revt kk_mode;


public: // ctor:

   explicit Mesh(
      real64 rmin,
      real64 rmax,
      real64 radial_mesh_ratio
   ) {
      rmin = rmin;
      rmax = rmax;
      // initialize shared memory block on F side.
      init_riemann_();
      // mesh_sector_general has some flexibility
      // that we do not wish to make use of,
      // instead we wish to use it to generate
      // global coordinates for the mesh.
      int32 ntp_patches_tmp = 1;
      int32 ix_patch = 1;
      int32 iy_patch = 1;
      int32 isector;
      real64 smallnum = P::smallnum;
      for (isector = 1; isector <= 6; ++isector) {
         mesh_sector_general_(
            &isector,
            &ntp_patches_tmp,
            &ix_patch,
            &iy_patch,
            &ntp_global,
            &NGhost_global,
            // OUTPUT:
            reinterpret_cast<real64*>(&vertices),
            reinterpret_cast<real64*>(&coordinates[isector-1])
         );
      }
      // todo for general ilevel_max > 1, call mesh_sector_general with
      //  ntp_global *= 2
      //  for each ilevel
      // done filling each sector, now fill ghost zones
      for (isector = 1; isector <= 6; ++isector) {
         fill_ghost2d_global_(
            &isector,
            &ntp_global,
            &NGhost_global,
            // OUTPUT:
            reinterpret_cast<real64*>(&coordinates)
         );
      }
      // radial mesh:
      if constexpr (P::Directives::radial_mesh == 0) {
         generate_global_radial_uniform_mesh_(
            &nr_global,
            &NGhost_global,
            &rmin,
            &rmax,
            &smallnum,
            &radial_mesh_ratio,
            // OUTPUT:
            reinterpret_cast<real64*>(&rb),
            reinterpret_cast<real64*>(&dr),
            reinterpret_cast<real64*>(&rc)
         );
      }
      else if constexpr (P::Directives::radial_mesh == 1) {
         generate_global_radial_log_mesh_(
            &nr_global,
            &NGhost_global,
            &rmin,
            &rmax,
            &smallnum,
            &radial_mesh_ratio,
            // OUTPUT:
            reinterpret_cast<real64*>(&rb),
            reinterpret_cast<real64*>(&dr),
            reinterpret_cast<real64*>(&rc)
         );
      }
      else if constexpr (P::Directives::radial_mesh == 2) {
         generate_global_radial_stretched_mesh_(
            &nr_global,
            &NGhost_global,
            &rmin,
            &rmax,
            &smallnum,
            &radial_mesh_ratio,
            // OUTPUT:
            reinterpret_cast<real64*>(&rb),
            reinterpret_cast<real64*>(&dr),
            reinterpret_cast<real64*>(&rc)
         );
      }
      else if constexpr (P::Directives::radial_mesh == 3) {
         generate_global_radial_hybrid_mesh_(
            &nr_global,
            &NGhost_global,
            &rmin,
            &rmax,
            &smallnum,
            &radial_mesh_ratio,
            // OUTPUT:
            reinterpret_cast<real64*>(&rb),
            reinterpret_cast<real64*>(&dr),
            reinterpret_cast<real64*>(&rc)
         );
      }
      else if constexpr (P::Directives::radial_mesh == 4) {
         generate_global_radial_powerlaw_mesh_(
            &nr_global,
            &NGhost_global,
            &rmin,
            &rmax,
            &smallnum,
            &radial_mesh_ratio,
            // OUTPUT:
            reinterpret_cast<real64*>(&rb),
            reinterpret_cast<real64*>(&dr),
            reinterpret_cast<real64*>(&rc)
         );
      }
      else if constexpr (P::Directives::radial_mesh == 5) {
         generate_global_radial_unif_n_log_mesh_(
            &nr_global,
            &NGhost_global,
            &rmin,
            &rmax,
            &smallnum,
            &radial_mesh_ratio,
            // OUTPUT:
            reinterpret_cast<real64*>(&rb),
            reinterpret_cast<real64*>(&dr),
            reinterpret_cast<real64*>(&rc)
         );
      }
      // Builds the indices that enable us to zip through the
      // zone-centered and face-centered reconstructions.
      // powers of polynomials (useful at higher orders)
      build_mode_indices_(
      &n_cc_modes,
      reinterpret_cast<int32*>(&ii_mode),
      reinterpret_cast<int32*>(&jj_mode),
      reinterpret_cast<int32*>(&kk_mode)
      );
   }

};


