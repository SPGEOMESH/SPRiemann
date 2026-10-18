
#pragma once

#include <csriemann/constants/types.h>

// Recall that when Fortran is compiled to
// library code using gfortran, symbols are
// coerced to lower case (My_Name -> my_name)
// and appended with an underscore '_'.

extern "C" {

// numorder == 2
void limit_all_slopes_p1(
      int32* nr_zones,
      int32* ntp_zones,
      int32* NGhost,
      // integers:
      int32* n_cc_components,
      int32* n_cc_modes,
      int32* n_fc_modes,
      int32* iz_patch,
      int32* nr_patches,
      // arrays:
      int32* ii_mode,
      int32* jj_mode,
      int32* kk_mode,
      // structures:
      real64* vert_coords_patch,
      real64* imask_active_patch,
      // xf
      real64* xf_quadpts_patch,
      real64* xf_quadwt_patch,
      real64* xf_quadnorml_patch,
      real64* xf_area_patch,
      real64* xf_centroid_patch,
      // yf
      real64* yf_quadpts_patch,
      real64* yf_quadwt_patch,
      real64* yf_quadnorml_patch,
      real64* yf_area_patch,
      real64* yf_centroid_patch,
      // zf
      real64* zf_quadpts_patch,
      real64* zf_quadwt_patch,
      real64* zf_quadnorml_patch,
      real64* zf_area_patch,
      real64* zf_centroid_patch,
      // xb
      real64* xb_quadpts_patch,
      real64* xb_quadwt_patch,
      real64* xb_quadnorml_patch,
      // yb
      real64* yb_quadpts_patch,
      real64* yb_quadwt_patch,
      real64* yb_quadnorml_patch,
      // zb
      real64* zb_quadpts_patch,
      real64* zb_quadwt_patch,
      real64* zb_quadnorml_patch,
      // v
      real64* v_scale_length,
      real64* v_volume_patch,
      real64* v_centroid_patch,
      // structure:
      real64* Coef_Basis_Moments,
      // integers:
      int32* itri_patch_LD,
      int32* itri_patch_RD,
      int32* itri_patch_RU,
      int32* itri_patch_LU,
      // structure:
      int32* istencil_patch,
      //---------------
      // Fluid Variables:
      real64* U,
      real64* entropygr,
      real64* Bxf,
      real64* Byf,
      real64* Bzf,
      real64* V,
      real64* B1_zc,
      real64* B2_zc,
      real64* B3_zc,
      real64* Ex_zc,
      real64* Ey_zc,
      real64* Ez_zc,
      real64* Wzc,
      real64* gamma_eos_gr,
      real64* gamma_soundspeed_gr,
      real64* mean_mol_wt_gr,
      real64* divvelgr,
      real64* msonicgr,
      real64* flattengr,
      real64* workgr,
      real64* Bxf_quadpts,
      real64* Byf_quadpts,
      real64* Bzf_quadpts,
      real64* V_min_local,
      real64* V_max_local
);



// numorder == 3
// NB: function signature is identical to p1!
void limit_all_slopes_p2_(
      int32* nr_zones,
      int32* ntp_zones,
      int32* NGhost,
      // integers:
      int32* n_cc_components,
      int32* n_cc_modes,
      int32* n_fc_modes,
      int32* iz_patch,
      int32* nr_patches,
      // arrays:
      int32* ii_mode,
      int32* jj_mode,
      int32* kk_mode,
      // structures:
      real64* vert_coords_patch,
      real64* imask_active_patch,
      // xf
      real64* xf_quadpts_patch,
      real64* xf_quadwt_patch,
      real64* xf_quadnorml_patch,
      real64* xf_area_patch,
      real64* xf_centroid_patch,
      // yf
      real64* yf_quadpts_patch,
      real64* yf_quadwt_patch,
      real64* yf_quadnorml_patch,
      real64* yf_area_patch,
      real64* yf_centroid_patch,
      // zf
      real64* zf_quadpts_patch,
      real64* zf_quadwt_patch,
      real64* zf_quadnorml_patch,
      real64* zf_area_patch,
      real64* zf_centroid_patch,
      // xb
      real64* xb_quadpts_patch,
      real64* xb_quadwt_patch,
      real64* xb_quadnorml_patch,
      // yb
      real64* yb_quadpts_patch,
      real64* yb_quadwt_patch,
      real64* yb_quadnorml_patch,
      // zb
      real64* zb_quadpts_patch,
      real64* zb_quadwt_patch,
      real64* zb_quadnorml_patch,
      // v
      real64* v_scale_length,
      real64* v_volume_patch,
      real64* v_centroid_patch,
      // structure:
      real64* Coef_Basis_Moments,
      // integers:
      int32* itri_patch_LD,
      int32* itri_patch_RD,
      int32* itri_patch_RU,
      int32* itri_patch_LU,
      // structure:
      int32* istencil_patch,
      //---------------
      // Fluid Variables:
      real64* U,
      real64* entropygr,
      real64* Bxf,
      real64* Byf,
      real64* Bzf,
      real64* V,
      real64* B1_zc,
      real64* B2_zc,
      real64* B3_zc,
      real64* Ex_zc,
      real64* Ey_zc,
      real64* Ez_zc,
      real64* Wzc,
      real64* gamma_eos_gr,
      real64* gamma_soundspeed_gr,
      real64* mean_mol_wt_gr,
      real64* divvelgr,
      real64* msonicgr,
      real64* flattengr,
      real64* workgr,
      real64* Bxf_quadpts,
      real64* Byf_quadpts,
      real64* Bzf_quadpts,
      real64* V_min_local,
      real64* V_max_local
);



// numorder == 3
// NB: function signature is identical to p1!
void limit_all_slopes_p12_(
      int32* nr_zones,
      int32* ntp_zones,
      int32* NGhost,
      // integers:
      int32* n_cc_components,
      int32* n_cc_modes,
      int32* n_fc_modes,
      int32* iz_patch,
      int32* nr_patches,
      // arrays:
      int32* ii_mode,
      int32* jj_mode,
      int32* kk_mode,
      // structures:
      real64* vert_coords_patch,
      real64* imask_active_patch,
      // xf
      real64* xf_quadpts_patch,
      real64* xf_quadwt_patch,
      real64* xf_quadnorml_patch,
      real64* xf_area_patch,
      real64* xf_centroid_patch,
      // yf
      real64* yf_quadpts_patch,
      real64* yf_quadwt_patch,
      real64* yf_quadnorml_patch,
      real64* yf_area_patch,
      real64* yf_centroid_patch,
      // zf
      real64* zf_quadpts_patch,
      real64* zf_quadwt_patch,
      real64* zf_quadnorml_patch,
      real64* zf_area_patch,
      real64* zf_centroid_patch,
      // xb
      real64* xb_quadpts_patch,
      real64* xb_quadwt_patch,
      real64* xb_quadnorml_patch,
      // yb
      real64* yb_quadpts_patch,
      real64* yb_quadwt_patch,
      real64* yb_quadnorml_patch,
      // zb
      real64* zb_quadpts_patch,
      real64* zb_quadwt_patch,
      real64* zb_quadnorml_patch,
      // v
      real64* v_scale_length,
      real64* v_volume_patch,
      real64* v_centroid_patch,
      // structure:
      real64* Coef_Basis_Moments,
      // integers:
      int32* itri_patch_LD,
      int32* itri_patch_RD,
      int32* itri_patch_RU,
      int32* itri_patch_LU,
      // structure:
      int32* istencil_patch,
      //---------------
      // Fluid Variables:
      real64* U,
      real64* entropygr,
      real64* Bxf,
      real64* Byf,
      real64* Bzf,
      real64* V,
      real64* B1_zc,
      real64* B2_zc,
      real64* B3_zc,
      real64* Ex_zc,
      real64* Ey_zc,
      real64* Ez_zc,
      real64* Wzc,
      real64* gamma_eos_gr,
      real64* gamma_soundspeed_gr,
      real64* mean_mol_wt_gr,
      real64* divvelgr,
      real64* msonicgr,
      real64* flattengr,
      real64* workgr,
      real64* Bxf_quadpts,
      real64* Byf_quadpts,
      real64* Bzf_quadpts,
      real64* V_min_local,
      real64* V_max_local
);



// Enforce positivity of the reconstructed solution.
void enforce_pressure_positivity_prim_(
      int32* nr_zones,
      int32* ntp_zones,
      int32* NGhost,
      int32* n_cc_components,
      int32* n_cc_modes,
      int32* n_fc_modes,
      //
      real64* vert_coords_patch,
      int32* imask_active_patch,
      real64* xf_centroid_patch,
      real64* yf_centroid_patch,
      real64* zf_centroid_patch,
      real64* v_scale_length,
      real64* v_centroid_patch,
      //
      int32* istencil_patch,
      //
      real64* V,
      real64* entropygr,
      real64* B1_zc,
      real64* B2_zc,
      real64* B3_zc,
      real64* flattengr,
      real64* Coef_Basis_Moments
);


void make_facial_fluxes_();

void make_edge_efields_p1();

void build_non_stiff_sources_();

void update_imex_rk_step_p1_();

void cleanup_after_update_();

void update_defect_imex_rk_();

void evaluate_timestep_();

void operator_split_source_();

void cooling_ptwise_();



}



