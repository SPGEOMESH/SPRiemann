
#pragma once

#include <csriemann/constants/types.h>

// Recall that when Fortran is compiled to
// library code using gfortran, symbols are
// coerced to lower case (My_Name -> my_name)
// and appended with an underscore '_'.

extern "C" {


void mesh_sector_general_(
        int32* isector,
        int32* ntp_patches,
        int32* ix_patch,
        int32* iy_patch,
        int32* ntp_zones,
        int32* NGhost,
        // OUTPUT:
        real64* vertices,
        real64* coordinates_gr
);

void fill_ghost2d_global_(
        int32* isector,
        int32* ntp_global,
        int32* NGhost_global,
        real64* coordinates_gr_global
);



void build_mode_indices_(
        int32* n_cc_modes,
        // OUTPUT:
        int32* ii_mode,
        int32* jj_mode,
        int32* kk_mode
);





void generate_global_radial_uniform_mesh_(
        int32* nr_global,
        int32* NGhost_global,
        real64* rmin,
        real64* rmax,
        real64* smallnum,
        real64* radial_mesh_ratio,
        // OUTPUT:
        real64* rb_save_global,
        real64* dr_save_global,
        real64* rc_save_global
);
void generate_global_radial_log_mesh_(
        int32* nr_global,
        int32* NGhost_global,
        real64* rmin,
        real64* rmax,
        real64* smallnum,
        real64* radial_mesh_ratio,
        // OUTPUT:
        real64* rb_save_global,
        real64* dr_save_global,
        real64* rc_save_global
);
void generate_global_radial_stretched_mesh_(
        int32* nr_global,
        int32* NGhost_global,
        real64* rmin,
        real64* rmax,
        real64* smallnum,
        real64* radial_mesh_ratio,
        // OUTPUT:
        real64* rb_save_global,
        real64* dr_save_global,
        real64* rc_save_global
);
void generate_global_radial_hybrid_mesh_(
        int32* nr_global,
        int32* NGhost_global,
        real64* rmin,
        real64* rmax,
        real64* smallnum,
        real64* radial_mesh_ratio,
        // OUTPUT:
        real64* rb_save_global,
        real64* dr_save_global,
        real64* rc_save_global
);
void generate_global_radial_powerlaw_mesh_(
        int32* nr_global,
        int32* NGhost_global,
        real64* rmin,
        real64* rmax,
        real64* smallnum,
        real64* radial_mesh_ratio,
        // OUTPUT:
        real64* rb_save_global,
        real64* dr_save_global,
        real64* rc_save_global
);
void generate_global_radial_unif_n_log_mesh_(
        int32* nr_global,
        int32* NGhost_global,
        real64* rmin,
        real64* rmax,
        real64* smallnum,
        real64* radial_mesh_ratio,
        // OUTPUT:
        real64* rb_save_global,
        real64* dr_save_global,
        real64* rc_save_global
);



void extract_vert_coords_patch_(
        // constexpr:
        int32* ntp_global,
        int32* nr_global,
        int32* NGhost_global,
        int32* ntp_zones,
        int32* nr_zones,
        int32* ntp_patches,
        int32* nr_patches,
        int32* NGhost,
        // int args:
        int32* isector,
        int32* ix_patch,
        int32* iy_patch,
        int32* iz_patch,
        int32* ilevel,
        // constexpr:
        int32* ilevel_max,
        // arrays:
        real64* coordinates_gr_global,
        real64* rb_save_global,
        // OUTPUT:
        // extracted vertex coordinates for this patch.
        real64* vert_coords_patch,
        // mask which is zone-centered and 1 if zone is active, 0 if inactive.
        int32* imask_active_patch
);




void stencil_indexing_tricorners_(
        int32* ntp_zones,
        int32* nr_zones,
        int32* ntp_patches,
        int32* nr_patches,
        int32* NGhost,
        int32* ix_patch,
        int32* iy_patch,
        // OUTPUT:
        int32* itri_patch_LD,
        int32* itri_patch_RD,
        int32* itri_patch_RU,
        int32* itri_patch_LU,
        // INPUT:
        int32* imask_active_patch,
        // OUTPUT:
        int32* istencil_patch
);



void coords_to_wt_quadpts_faces_iso2_(
        int32* nr_zones,
        int32* ntp_zones,
        int32* NGhost,
        // arrays:
        real64* vert_coords_patch,
        int32* imask_active_patch,
        // OUTPUT:
        // x
        real64* xf_quadpts_patch,
        real64* xf_quadwt_patch,
        real64* xf_quadnorml_patch,
        real64* xf_area_patch,
        real64* xf_centroid_patch,
        // y
        real64* yf_quadpts_patch,
        real64* yf_quadwt_patch,
        real64* yf_quadnorml_patch,
        real64* yf_area_patch,
        real64* yf_centroid_patch,
        // z
        real64* zf_quadpts_patch,
        real64* zf_quadwt_patch,
        real64* zf_quadnorml_patch,
        real64* zf_area_patch,
        real64* zf_centroid_patch,
        // OUTPUT:
        // a scale length for each zone
        real64* v_scale_length
);



void coords_to_wt_quadpts_bfld_iso2_(
        int32* nr_zones,
        int32* ntp_zones,
        int32* NGhost,
        // arrays:
        real64* vert_coords_patch,
        int32* imask_active_patch,
        // OUTPUT:
        // x
        real64* xb_quadpts_patch,
        real64* xb_quadwt_patch,
        real64* xb_quadnorml_patch,
        // y
        real64* yb_quadpts_patch,
        real64* yb_quadwt_patch,
        real64* yb_quadnorml_patch,
        // z
        real64* zb_quadpts_patch,
        real64* zb_quadwt_patch,
        real64* zb_quadnorml_patch
);



void coords_to_wt_quadpts_volume_iso2_(
        int32* nr_zones,
        int32* ntp_zones,
        int32* NGhost,
        // arrays:
        real64* vert_coords_patch,
        int32* imask_active_patch,
        // OUTPUT:
        real64* v_quadpts_patch,
        real64* v_quadwt_patch,
        real64* v_quadtangent_patch,
        real64* v_volume_patch,
        real64* v_centroid_patch
);



void quick_quadpts_volume_iso2_(
        int32* nr_zones,
        int32* ntp_zones,
        int32* NGhost,
        // arrays:
        real64* vert_coords_patch,
        int32* imask_active_patch,
        // OUTPUT:
        real64* vq_quadpts_patch,
        real64* vq_quadwt_patch
);



void init_basis_function_coefs_(
        int32* nr_zones,
        int32* ntp_zones,
        int32* NGhost,
        // the size of the `mode` arrays
        int32* n_cc_modes,
        // `mode` arrays:
        int32* ii_mode,
        int32* jj_mode,
        int32* kk_mode,
        // arrays:
        int32* imask_active_patch,
        // arrays:
        real64* v_quadpts_patch,
        real64* v_quadwt_patch,
        real64* v_volume_patch,
        real64* v_centroid_patch,
        real64* v_scale_length,
        // OUTPUT:
        real64* Coef_Basis_Moments
);




void reconstruction_coefs_(
        int32* iorder,
        int32* n_cc_modes,
        real64* x,
        real64* y,
        real64* z,
        real64* ell_0,
        real64* ell_j,
        // arrays:
        real64* Coef_0,
        real64* Coef_j,
        // OUTPUT:
        real64* Amat_Row
);



void interpolation_coefs_(
        int32* iorder,
        int32* n_cc_modes,
        real64* x,
        real64* y,
        real64* z,
        // arrays:
        real64* Coef_0,
        // OUTPUT:
        real64* Amat_Row10
);



void cmat_row_coefs_p2_(
        real64* x,
        real64* y,
        real64* z,
        real64* nx,
        real64* ny,
        real64* nz,
        // arrays:
        real64* Coef_0,
        // OUTPUT:
        real64* Cmat_Row_B
);


void init_patch_cube_sphere_(
      int32* nr_zones,
      int32* ntp_zones,
      int32* NGhost,
      // nvar+nfluid
      int32* n_cc_components,
      // 4 or 10 or 20 depending on numorder
      int32* n_cc_modes,
      // 1
      int32* n_fc_modes,
      int32* iproblem,
      real64* time,
      real64* dtcur,
      real64* rmin,
      real64* rmax,
      // integer arrays:
      real64* vert_coords_patch,
      int32* imask_active_patch,
      // real volume arrays:
      real64* v_scale_length,
      // arrays used for curl-A calculation.
      // x
      real64* xf_quadpts_patch,
      real64* xf_quadwt_patch,
      real64* xf_quadnorml_patch,
      real64* xf_area_patch,
      real64* xf_centroid_patch,
      // y
      real64* yf_quadpts_patch,
      real64* yf_quadwt_patch,
      real64* yf_quadnorml_patch,
      real64* yf_area_patch,
      real64* yf_centroid_patch,
      // z
      real64* zf_quadpts_patch,
      real64* zf_quadwt_patch,
      real64* zf_quadnorml_patch,
      real64* zf_area_patch,
      real64* zf_centroid_patch,
      // arrays used for volume integration.
      // v
      real64* v_quadpts_patch,
      real64* v_quadwt_patch,
      // -----------------------------
      // Arrays representing the MHD variables
      // in volume and on faces.
      real64* U,
      real64* V,
      real64* entropygr,
      real64* Bxf,
      real64* Byf,
      real64* Bzf,
      real64* B0field,
      real64* gamma_eos_gr,
      real64* V_min_local,
      real64* V_max_local
);



// todo - refactor on fortran side, extract
//  from module and formulate signature
void build_base_patches_(

      // ...

);




}




