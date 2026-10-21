


#pragma once

#include <Tbasis.h>

#include <spriemann/constants/parameters.h>
#include <spriemann/constants/directives.h>
#include <spriemann/mesh.h>
#include <spriemann/patch/patch.h>


class Tquadpoints {
public:

    using Directives = Directives<>;

    static constexpr size_t ntp_patches = 1;
    static constexpr size_t nr_patches = 1;

    static constexpr size_t ntp_zones = 1;
    static constexpr size_t nr_zones = 1;
    static constexpr size_t NGhost = 0;
    static constexpr size_t ilevel_max = 3;
    static constexpr size_t ilevel = 2;
    static constexpr size_t iscale = DGCALC::ipow2<size_t>(ilevel_max-ilevel);

    static constexpr size_t ntp_global = ntp_patches*ntp_zones*iscale;
    static constexpr size_t nr_global = nr_patches*nr_zones*iscale;
    static constexpr size_t NGhost_global = NGhost*iscale;

    static constexpr size_t n_cc_components = 8;
    static constexpr size_t n_fc_modes = 1;


    using P = Parameters<
            Directives,
            ntp_zones,
            nr_zones,
            ntp_patches,
            nr_patches,
            NGhost,
            ilevel_max
    >;

   // This test mocks the Fortran source code test:
   // TEST_QUADPOINTS
    static void test_quadpoints() {

        real64 rmin = 10.0;
        real64 rmax = 11.0;
        real64 radial_mesh_ratio = 1.0;

        auto mesh = Mesh<P>(
                rmin,
                rmax,
                radial_mesh_ratio
        );

        int32 isector = 1;
        int32 ix = 1;
        int32 iy = 1;
        int32 iz = 1;
        int32 ilevel_ = ilevel;

        auto patch = Patch<P>();

        patch.initialize_geometry(
             mesh,
             isector,
             ix,
             iy,
             iz,
             ilevel_
        );

//        cout << patch.vert_coords.str() << endl;
        cout << patch.imask_active.str() << endl;
//        cout << patch.xf_quadpts.str() << endl;

    }


};



