
#pragma once

#include <csriemann/constants/parameters.h>
#include <csriemann/mesh.h>
#include <csriemann/patch/patch.h>


class TPatch {
public:

    using Directives = Directives<>;

    static constexpr size_t ntp_patches = 1;
    static constexpr size_t nr_patches = 1;

    static constexpr size_t ntp_zones = 5;
    static constexpr size_t nr_zones = 5;
    static constexpr size_t NGhost = 2;
    static constexpr size_t ilevel_max = 1;

    static constexpr size_t ntp_global = ntp_patches*ntp_zones;
    static constexpr size_t nr_global = nr_patches*nr_zones;
    static constexpr size_t NGhost_global = NGhost;

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
        int32 ilevel = 1;

        auto patch = Patch<P>(
                mesh,
                isector,
                ix,
                iy,
                iz,
                ilevel
        );


        // CHECK
        // extract_vert_coords_patch:

//        cout << patch.vert_coords.str() << endl;
//        cout << patch.imask_active.str() << endl;

        // CHECK
        // stencil_indexing_tricorners:

//        cout << patch.itri_LD << endl;
//        cout << patch.itri_RD << endl;
//        cout << patch.itri_RU << endl;
//        cout << patch.itri_LU << endl;
//        cout << patch.istencil.str() << endl;

        // CHECK
        // coords_to_wt_quadpts_faces_iso2:

        // x
//        cout << patch.xf_quadpts.str() << endl;
//        cout << patch.xf_quadwt.str() << endl;
//        cout << patch.xf_quadnorml.str() << endl;
//        cout << patch.xf_area.str() << endl;
//        cout << patch.xf_centroid.str() << endl;
        // y
//        cout << patch.yf_quadpts.str() << endl;
//        cout << patch.yf_quadwt.str() << endl;
//        cout << patch.yf_quadnorml.str() << endl;
//        cout << patch.yf_area.str() << endl;
//        cout << patch.yf_centroid.str() << endl;
//        // z
//        cout << patch.zf_quadpts.str() << endl;
//        cout << patch.zf_quadwt.str() << endl;
//        cout << patch.zf_quadnorml.str() << endl;
//        cout << patch.zf_area.str() << endl;
//        cout << patch.zf_centroid.str() << endl;
//        //
//        cout << patch.v_scale_length.str() << endl;


        // CHECK
        // coords_to_wt_quadpts_bfld_iso2:

//        // x
//        cout << patch.xb_quadpts.str() << endl;
//        cout << patch.xb_quadwt.str() << endl;
//        cout << patch.xb_quadnorml.str() << endl;
//        // y
//        cout << patch.yb_quadpts.str() << endl;
//        cout << patch.yb_quadwt.str() << endl;
//        cout << patch.yb_quadnorml.str() << endl;
//        // z
//        cout << patch.zb_quadpts.str() << endl;
//        cout << patch.zb_quadwt.str() << endl;
//        cout << patch.zb_quadnorml.str() << endl;

        // CHECK
        // coords_to_wt_quadpts_volume_iso2:

//        cout << "v_quadpts\n" << patch.v_quadpts.str() << endl;
//        cout << "v_quadwt\n" << patch.v_quadwt.str() << endl;
//        cout << "v_quadtangent\n" << patch.v_quadtangent.str() << endl;
//        cout << "v_volume\n" << patch.v_volume.str() << endl;
//        cout << "v_centroid\n" << patch.v_centroid.str() << endl;

        // CHECK
        // init_basis_function_coefs:
        cout << "coef_basis_moments\n" << patch.Coef_Basis_Moments.str() << endl;


    }


    static void Txf_quadpts_patch() {

        real64 rmin = 1.5;
        real64 rmax = 2.5;
        real64 radial_mesh_ratio = 2.0;

        auto mesh = Mesh<P>(
                rmin,
                rmax,
                radial_mesh_ratio
        );

        int32 isector = 1;
        int32 ix = 1;
        int32 iy = 1;
        int32 iz = 1;
        int32 ilevel = 1;

        auto patch = Patch<P>(
                mesh,
                isector,
                ix,
                iy,
                iz,
                ilevel
        );


        // Print while spoofing output from a Fortran test.
        // Use to confirm all is well.


        auto vcp = patch.vert_coords;
        for (size_t k = 0; k < P::Ntp_zones+1; ++k) {
            for (size_t j = 0; j < P::Ntp_zones+1; ++j) {
                for (size_t i = 0; i < P::Nr_zones+1; ++i) {
                    cout << "i j k " << P::N(i) << " " << P::N(j) << " " << P::N(k) << endl;
                    cout << "vcp " << "\t\t";
                    for (size_t ii = 0; ii < 3; ++ii) {
                        cout << vcp[{ii,k}][{j,i}].str_dat() << "\t\t";
                    }
                    cout << "\n";
                }
            }
        }


        auto xfq = patch.xf_quadpts;
        for (size_t k = 0; k < P::Ntp_zones+1; ++k) {
            for (size_t j = 0; j < P::Ntp_zones; ++j) {
                for (size_t i = 0; i < P::Nr_zones; ++i) {
                    cout << "i j k " << P::N(i) << " " << P::N(j+1) << " " << P::N(k+1) << " " << P::NAREA_QUAD_PTS << " area quadrature points" << endl;
                    for (size_t quadi = 0; quadi < P::NAREA_QUAD_PTS; ++quadi) {
                        cout << "quadi " << quadi << "\t\t";
                        for (size_t ii = 0; ii < 3; ++ii) {
                            cout << xfq[quadi][{ii,k}][{j,i}].str_dat() << "\t\t";
                        }
                        cout << "\n";
                    }
                }
            }
        }


    }


};



