
#pragma once

#include <csriemann/mesh.h>
#include <csriemann/patch/patch_index.h>

/*
 * A data structure
 * consisting of a logically cubic
 * layout of zones with edge-,
 * face-, and zone-centered quantities,
 * with update routines that are
 * mostly implemented in Fortran
 * and linked.
 *
 */
template <
      typename Parameters
>
class Patch {
public: // types:

   using P = Parameters;
   using Mesh = Mesh<P>;

   static constexpr size_t Ntp_zones = P::Ntp_zones;
   static constexpr size_t Nr_zones = P::Nr_zones;
   static constexpr size_t NAREA_QUAD_PTS = P::NAREA_QUAD_PTS;
   static constexpr size_t NBFLD_QUAD_PTS = P::NBFLD_QUAD_PTS;
   static constexpr size_t SILO_EXTRA_DATA = P::SILO_EXTRA_DATA;
   static constexpr size_t NUMORDER = P::NUMORDER;
   static constexpr size_t NUMORDER3 = P::NUMORDER3;
   static constexpr size_t NGhost = P::NGhost;
   static constexpr size_t n_cc_modes = P::n_cc_modes;
   static constexpr size_t n_cc_components = P::n_cc_components;

   static constexpr int32 nRK_stage = P::nRK_stage;

public: // global mesh buffer fields:

   using coords_3 = QuadMatrix<R64, 3, Nr_zones+1, Ntp_zones+1, Ntp_zones+1>;
   using mask_1 = TriMatrix<N32, Nr_zones, Ntp_zones, Ntp_zones>;
   using index_12_2 = QuadMatrix<N32, 12, 2, Ntp_zones, Ntp_zones>;

   coords_3 vert_coords;
   mask_1 imask_active;
   index_12_2 istencil;

public: // zone-computing helper fields:

// ----------------------------------------------------------------------
// These are types for the facial areas, facial quadrature points, facial centroids, facial normals
// and facial weights. The actual variables that live in the x-, y-, z-faces can be
// specified with these types.

   using xf_3 =  TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using xf_43 = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using xf_4N = QuadMatrix<R64, NAREA_QUAD_PTS, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using xf_5 =  PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;

   xf_3    xf_area;
   xf_43   xf_quadnorml, xf_centroid;
   xf_4N   xf_quadwt;
   xf_5    xf_quadpts;

   using yf_3 =  TriMatrix<R64, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using yf_43 = QuadMatrix<R64, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using yf_4N = QuadMatrix<R64, NAREA_QUAD_PTS, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using yf_5 =  PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;

   yf_3    yf_area;
   yf_43   yf_quadnorml, yf_centroid;
   yf_4N   yf_quadwt;
   yf_5    yf_quadpts;

   using zf_3 =  TriMatrix<R64, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using zf_43 = QuadMatrix<R64, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using zf_4N = QuadMatrix<R64, NAREA_QUAD_PTS, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using zf_5 =  PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;

   zf_3    zf_area;
   zf_43   zf_centroid;
   zf_4N   zf_quadwt;
   zf_5    zf_quadpts, zf_quadnorml;

// ----------------------------------------------------------------------
// Integration points for magnetic field in the multid riemann solver.
// These are also the locations where the magnetic field is directly reconstructed.

   using xb_43 = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using xb_4N = QuadMatrix<R64, NBFLD_QUAD_PTS, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using xb_5 =  PentaMatrix<R64, NBFLD_QUAD_PTS, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;

   xb_43   xb_quadnorml;
   xb_4N   xb_quadwt;
   xb_5    xb_quadpts;

   using yb_43 = QuadMatrix<R64, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using yb_4N = QuadMatrix<R64, NBFLD_QUAD_PTS, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using yb_5  = PentaMatrix<R64, NBFLD_QUAD_PTS, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;

   yb_43   yb_quadnorml;
   yb_4N   yb_quadwt;
   yb_5    yb_quadpts;

   using zb_4N = QuadMatrix<R64, NBFLD_QUAD_PTS, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using zb_5 =  PentaMatrix<R64, NBFLD_QUAD_PTS, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;

   zb_4N   zb_quadwt;
   zb_5    zb_quadpts, zb_quadnorml;

// ----------------------------------------------------------------------
// These are types for the volumetric scale length, volume, volumetric quadrature
// volumetric quadrature weights, centroids etc.

   using v_3 =  TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_43 = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_4N = QuadMatrix<R64, NUMORDER3, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_4n_cc = QuadMatrix<R64, n_cc_modes, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_5 =  PentaMatrix<R64, NUMORDER3, 3, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_6 = HexaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, 3, 3, NUMORDER3>;

   v_3    v_scale_length, v_volume;
   v_43   v_centroid;
   v_4N   v_quadwt;
   v_6 v_quadtangent;
   v_4n_cc Coef_Basis_Moments;
   v_5    v_quadpts;

   // ----------------------------------------------------------------------
// Minimal amounts of MHD information that is stored for each of the
// patches that reside on this processor.

   using mhd_zc_4 = QuadMatrix<R64, n_cc_components, Nr_zones, Ntp_zones, Ntp_zones>;
   using mhd_zc_43 = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones>;

   using mhd_fc_4 = QuadMatrix<R64, n_cc_components, Nr_zones+1, Ntp_zones+1, Ntp_zones+1>;

   using mhd_zc_3 = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;

   using mhd_fc_3 = TriMatrix<R64, Nr_zones+1, Ntp_zones+1, Ntp_zones+1>;

   using mhd_silo_zc_4 = QuadMatrix<R64, n_cc_components + SILO_EXTRA_DATA,
         Nr_zones, Ntp_zones, Ntp_zones>;

   using mhd_fc_rlud_2 = Matrix<R64, Nr_zones+1, Ntp_zones+1>;
   using mhd_fc_tb_2 = Matrix<R64, Ntp_zones+1, Ntp_zones+1>;

   mhd_zc_4   U, V, U0, Source;

   mhd_silo_zc_4  U_silo;

   mhd_zc_43  B0field;

   // These are all the fluctuations and their defects.
   mhd_fc_4   Fluct_x_Pl, Fluct_x_Mi, Fluct_y_Pl, Fluct_y_Mi, Fluct_z_Pl, Fluct_z_Mi,
         Defect_Fluct_x_Pl, Defect_Fluct_x_Mi, Defect_Fluct_y_Pl, Defect_Fluct_y_Mi,
         Defect_Fluct_z_Pl, Defect_Fluct_z_Mi;

   mhd_zc_3   entropygr, entropygr0, divgr, flag_entprs, flag_Tfloor,
         CGL_fbeta_gr_amr, Pelectron_gr_amr,
         gamma_eos_gr;

   // Zone-centered flattener for all patches on this processor.
   mhd_zc_3 flattengr;

   mhd_fc_3   Bxf, Byf, Bzf, Bxf0, Byf0, Bzf0,   // facial magnetic field components.
   Efx, Efy, Efz, Defect_Efx, Defect_Efy, Defect_Efz,   // Electric fields and defects.
   EntFlux_x, EntFlux_y, EntFlux_z,
         Defect_EntFlux_x, Defect_EntFlux_y, Defect_EntFlux_z;  // Entropy fluxes and defects.

   // todo ask if this is window-comm buffer or some other kind of buffer.

   mhd_fc_rlud_2 Bxf_rt_buf, Bxf_lf_buf, Byf_up_buf, Byf_dn_buf,   // Buffer for Bxf, Byf
   Efy_rt_buf, Efz_rt_buf, Efy_lf_buf, Efz_lf_buf,      // Buffer for Efy and Efz at x-faces.
   Efx_up_buf, Efz_up_buf, Efx_dn_buf, Efz_dn_buf;      // Buffer for Efx and Efz at y-faces.

   mhd_fc_tb_2  Bz_tp_buf, Bz_bt_buf,                          // Buffer for Bzf
   Efx_tp_buf, Efy_tp_buf, Efx_bt_buf, Efy_bt_buf;   // Buffer for Efx and Efy at z-baces.


public: // patch-wide quantities and patch attributes:

   using V_1 = Vector<R64, n_cc_components+1>;

   // markers for tricorners (0: no tricorner, 1: tricorner)
   int32 itri_LD;
   int32 itri_RD;
   int32 itri_RU;
   int32 itri_LU;

   // time t for this patch.
   real64 time;
   // dtcur (dt-current, the current timestep distance) for this patch.
   real64 dtcur;

   // A lookup key for the patch inventory.
   PatchIndex index;

   // Is top or bottom patch.
   // i.e., has border at rmin or rmax.
   // Just needed for ADD
   // and PAD boundary conditions.
   int32 is_bottom_or_top;


   // Patch-level predicted timestep ("dt-current plus one")
   real64 dtcurp1;

   // min and max, for each conserved
   // physical quantity, across the entire patch.
   V_1 V_min;
   V_1 V_max;



public: // ctor:

   // constructor for blank initialization,
   // because we have to go through some MPI
   // before we are ready to construct.
   explicit Patch() {
      // do nothing.
   }

   Patch(
         Mesh& mesh,
         int32 isector,
         int32 ix,
         int32 iy,
         int32 iz,
         int32 ilevel
   ) {
      extract_vert_coords_patch_(
         &mesh.ntp_global,
         &mesh.nr_global,
         &mesh.NGhost_global,
         &mesh.ntp_zones,
         &mesh.nr_zones,
         &mesh.ntp_patches,
         &mesh.nr_patches,
         &mesh.NGhost,
         &isector,
         &ix,
         &iy,
         &iz,
         &ilevel,
         &mesh.ilevel_max,
         reinterpret_cast<real64*>(&mesh.coordinates),
         reinterpret_cast<real64*>(&mesh.rb),
         // OUTPUT:
         reinterpret_cast<real64*>(&vert_coords),
         reinterpret_cast<int32*>(&imask_active)
      );
      stencil_indexing_tricorners_(
         &mesh.ntp_zones,
         &mesh.nr_zones,
         &mesh.ntp_patches,
         &mesh.nr_patches,
         &mesh.NGhost,
         &ix,
         &iy,
         // OUTPUT:
         &itri_LD,
         &itri_RD,
         &itri_RU,
         &itri_LU,
         // INPUT:
         reinterpret_cast<int32*>(&imask_active),
         // OUTPUT:
         reinterpret_cast<int32*>(&istencil)
      );
      coords_to_wt_quadpts_faces_iso2_(
          &mesh.nr_global,
          &mesh.ntp_zones,
          &mesh.NGhost,
          reinterpret_cast<real64*>(&vert_coords),
          reinterpret_cast<int32*>(&imask_active),
          // OUTPUT:
          // x
          reinterpret_cast<real64*>(&xf_quadpts),
          reinterpret_cast<real64*>(&xf_quadwt),
          reinterpret_cast<real64*>(&xf_quadnorml),
          reinterpret_cast<real64*>(&xf_area),
          reinterpret_cast<real64*>(&xf_centroid),
          // y
          reinterpret_cast<real64*>(&yf_quadpts),
          reinterpret_cast<real64*>(&yf_quadwt),
          reinterpret_cast<real64*>(&yf_quadnorml),
          reinterpret_cast<real64*>(&yf_area),
          reinterpret_cast<real64*>(&yf_centroid),
          // z
          reinterpret_cast<real64*>(&zf_quadpts),
          reinterpret_cast<real64*>(&zf_quadwt),
          reinterpret_cast<real64*>(&zf_quadnorml),
          reinterpret_cast<real64*>(&zf_area),
          reinterpret_cast<real64*>(&zf_centroid),
          // OUTPUT:
          reinterpret_cast<real64*>(&v_scale_length)
      );
      coords_to_wt_quadpts_bfld_iso2_(
          &mesh.nr_zones,
          &mesh.ntp_zones,
          &mesh.NGhost,
          reinterpret_cast<real64*>(&vert_coords),
          reinterpret_cast<int32*>(&imask_active),
          // OUTPUT:
          // x
          reinterpret_cast<real64*>(&xb_quadpts),
          reinterpret_cast<real64*>(&xb_quadwt),
          reinterpret_cast<real64*>(&xb_quadnorml),
          // y
          reinterpret_cast<real64*>(&yb_quadpts),
          reinterpret_cast<real64*>(&yb_quadwt),
          reinterpret_cast<real64*>(&yb_quadnorml),
          // z
          reinterpret_cast<real64*>(&zb_quadpts),
          reinterpret_cast<real64*>(&zb_quadwt),
          reinterpret_cast<real64*>(&zb_quadnorml)
      );
      coords_to_wt_quadpts_volume_iso2_(
          &mesh.nr_zones,
          &mesh.ntp_zones,
          &mesh.NGhost,
          reinterpret_cast<real64*>(&vert_coords),
          reinterpret_cast<int32*>(&imask_active),
          // OUTPUT:
          reinterpret_cast<real64*>(&v_quadpts),
          reinterpret_cast<real64*>(&v_quadwt),
          reinterpret_cast<real64*>(&v_quadtangent),
          reinterpret_cast<real64*>(&v_volume),
          reinterpret_cast<real64*>(&v_centroid)
      );
      init_basis_function_coefs_(
          &mesh.nr_zones,
          &mesh.ntp_zones,
          &mesh.NGhost,
          &mesh.n_cc_modes,
          reinterpret_cast<int32*>(&mesh.ii_mode),
          reinterpret_cast<int32*>(&mesh.jj_mode),
          reinterpret_cast<int32*>(&mesh.kk_mode),
          reinterpret_cast<int32*>(&imask_active),
          reinterpret_cast<real64*>(&v_quadpts),
          reinterpret_cast<real64*>(&v_quadwt),
          reinterpret_cast<real64*>(&v_volume),
          reinterpret_cast<real64*>(&v_centroid),
          reinterpret_cast<real64*>(&v_scale_length),
          reinterpret_cast<real64*>(&Coef_Basis_Moments)
      );
      // By this point, all geometry is built. Therefore, we are now in a position
      // to initialize the MHD variables for this patch.
      int32 iproblem = P::iproblem;
      init_patch_cube_sphere_(
            &mesh.nr_zones,
            &mesh.ntp_zones,
            &mesh.NGhost,
            &mesh.n_cc_components,
            &mesh.n_cc_modes,
            &mesh.n_fc_modes,
            &iproblem,
            // This method uses time t.
            &time,
            // This method will initialize dtcur for this patch.
            &dtcur,
            //
            &mesh.rmin,
            &mesh.rmax,
            reinterpret_cast<real64*>(&vert_coords),
            reinterpret_cast<int32*>(&imask_active),
            // x
            reinterpret_cast<real64*>(&v_scale_length),
            reinterpret_cast<real64*>(&xf_quadpts),
            reinterpret_cast<real64*>(&xf_quadwt),
            reinterpret_cast<real64*>(&xf_quadnorml),
            reinterpret_cast<real64*>(&xf_area),
            reinterpret_cast<real64*>(&xf_centroid),
            // y
            reinterpret_cast<real64*>(&yf_quadpts),
            reinterpret_cast<real64*>(&yf_quadwt),
            reinterpret_cast<real64*>(&yf_quadnorml),
            reinterpret_cast<real64*>(&yf_area),
            reinterpret_cast<real64*>(&yf_centroid),
            // z
            reinterpret_cast<real64*>(&zf_quadpts),
            reinterpret_cast<real64*>(&zf_quadwt),
            reinterpret_cast<real64*>(&zf_quadnorml),
            reinterpret_cast<real64*>(&zf_area),
            reinterpret_cast<real64*>(&zf_centroid),
            // v
            reinterpret_cast<real64*>(&v_quadpts),
            reinterpret_cast<real64*>(&v_quadwt),
            //
            reinterpret_cast<real64*>(&U),
            reinterpret_cast<real64*>(&V),
            reinterpret_cast<real64*>(&entropygr),
            reinterpret_cast<real64*>(&Bxf),
            reinterpret_cast<real64*>(&Byf),
            reinterpret_cast<real64*>(&Bzf),
            reinterpret_cast<real64*>(&B0field),
            reinterpret_cast<real64*>(&gamma_eos_gr),
            reinterpret_cast<real64*>(&V_min),
            reinterpret_cast<real64*>(&V_max)
      );
   }


private: // methods:

   // This is called during every time step update
   // in order to ()()() ???
   void main_boundary_logic() {
      // todo use ipatch to point Fortran to the right patch
//      bcarr_r = 1.234;
//      pad_boundary_r();
//      ibdy_inner();
//      ibdy_outer();
//      add_boundary();
//      operation_split_source();
   }

   /*
    * Process the slope limiting.
    * This is broken out so we can
    * fumble with the order.
    *
    */
   void main_limit_all_slopes() {
      if constexpr (NUMORDER == 2) {
//         limit_all_slopes_p1_();
      }
      else if constexpr (NUMORDER == 3) {
//         limit_all_slopes_p12_();
      } else {
         // ?
      }
   }

public: // API:


   // This is the primary ordered computational sequence
   // during which the patch data is updated.
   // It will be called for each managed patch
   // by the managing processor.
   void main_logic(int32 rki, int32 ipatch) {

      main_boundary_logic();

      main_limit_all_slopes();

//      make_facial_fluxes();

//      make_edge_efields_p1();

//      build_non_stiff_source();

      if (rki == 0) {
         // todo
         ;
      }
      if (rki == 1) {
         // todo
         ;
      }
//      update_imex_rk_step_p1();
      // are all variables in physical bounds?
//      cleanup_after_update();
      if (rki == nRK_stage) {
//         evaluate_timestep();
//         operator_split_source();
      }
      // todo fill_flux_snd_patch... skinny variables...
//      copy__data__back__from__gpu();
      // todo end
      if (rki == 0) {
//         reset_cgl();
//         reset_pelectron();
      }
      // end of main logic.
   }


   // This is the first of two steps
   // related to defect updates during
   // the cleanup stage of the update.
   void defect_logic1() {
      // zero out a ton of defect fields
//      defect_ent_flux_x = 1.234;
   }

   // This is the second of two steps
   // related to defect updates during
   // the cleanup stage of the update.
   void defect_logic2() {
//      update_defect_imex_rk();
   }


};


