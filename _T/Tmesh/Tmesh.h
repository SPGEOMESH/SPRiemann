
#include <Tbasis.h>

#include <spriemann/constants/parameters.h>
#include <spriemann/constants/directives.h>
#include <spriemann/mesh.h>


class TMesh {
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

        real64 rmin = 1.0;
        real64 rmax = 2.0;
        real64 radial_mesh_ratio = 2.0;

        auto mesh = Mesh<P>(
                rmin,
                rmax,
                radial_mesh_ratio
        );

        cout << mesh.vertices.str() << endl;
        cout << mesh.coordinates.str() << endl;
//        cout << mesh.rb.str() << endl;
//        cout << mesh.dr.str() << endl;
//        cout << mesh.rc.str() << endl;

        cout << mesh.ii_mode.str() << endl;
        cout << mesh.jj_mode.str() << endl;
        cout << mesh.kk_mode.str() << endl;


    }


};

