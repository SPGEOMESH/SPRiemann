

#pragma once


#include <csriemann/constants/types.h>




template <
        typename Directives_,
        int32 ntp_zones_,
        int32 nr_zones_,
        int32 ntp_patches_,
        int32 nr_patches_,
        // number of ghost zones
        int32 NGhost_,
        int32 ilevel_max_
>
class Parameters {
public: // fields:

   using Directives = Directives_;

   static constexpr int32 NUMORDER = Directives::NUMORDER;
   static constexpr int32 NUMORDER3 = NUMORDER*NUMORDER*NUMORDER;

   static constexpr int32 NGhost = NGhost_;
   static constexpr int32 NGhost_global = NGhost;

   static constexpr int32 ntp_zones = ntp_zones_;
   static constexpr int32 nr_zones = nr_zones_;
   static constexpr int32 ntp_patches = ntp_patches_;
   static constexpr int32 nr_patches = nr_patches_;
   // In Fortran, the levels are 1, 2, ..., ilevel_max,
   // so ilevel_max is also the number of levels.
   static constexpr int32 ilevel_max = ilevel_max_;

   static constexpr int32 NVAR = Directives::NVAR;
   static constexpr int32 NFLUID = Directives::NFLUID;

   static constexpr real64 smallnum = Directives::smallnum;

   static_assert(ilevel_max >= 1, "The ilevel_max must be at least 1.");
   static_assert(NGhost <= ntp_zones, "It is required that NGhost ≤ ntp_zones.");
   static_assert(NGhost <= nr_zones, "It is required that NGhost ≤ nr_zones.");
   static_assert(ntp_patches >= 1, "It is required that ntp_patches ≥ 1.");
   static_assert(nr_patches >= 1, "It is required that nr_patches ≥ 1.");

public: // derived fields:

   static constexpr int32 ntp_global = ntp_zones*ntp_patches;
   static constexpr int32 nr_global = nr_zones*nr_patches;

   // A general convention in this class (used above)
   // is that lowercase n means "without ghost positions"
   // and upper-case N means "with ghost positions."

   static constexpr int32 Ntp_zones = ntp_zones+2*NGhost;
   static constexpr int32 Nr_zones = nr_zones+2*NGhost;
   static constexpr int32 Ntp_patches = ntp_patches+2*NGhost;
   static constexpr int32 Nr_patches = nr_patches+2*NGhost;

   static constexpr int32 Ntp_global = ntp_global+2*NGhost_global;
   static constexpr int32 Nr_global = nr_global+2*NGhost_global;

   //  scale between base and finest meshes
   static constexpr int32 iscale = DGCALC::ipow2<int32>(ilevel_max-1);

   // elsewhere NGhost_global = NGhost*iscale, but I cannot convince myself of this.

   // Number of zones within a sector on the base mesh
   static constexpr int32 ntp_global_base = ntp_zones * ntp_patches;
   static constexpr int32 nr_global_base = nr_zones * nr_patches;

   // Max theoretical patches within a sector for "ilevel_max" hierarchy.
   // Also, number of patches at base level within a sector and max patches
   // actually allowed on a processor.
   static constexpr int32 Max_base_patches_per_sector = ntp_patches * ntp_patches * nr_patches;

   // sum geometric series with ratio 8
   static constexpr int32 Max_Theoretical_patches_per_sector = (DGCALC::ipow2<int32>(ilevel_max)-1) / 7;

   // The number of RK stages, depends on the order
   static constexpr int32 nRK_stage = (NUMORDER == 2)? 2: (NUMORDER == 3)? 3: 5;


   static constexpr int32 NAREA_QUAD_PTS = Directives::NAREA_QUAD_PTS;
   static constexpr int32 NBFLD_QUAD_PTS = Directives::NBFLD_QUAD_PTS;
   static constexpr int32 SILO_EXTRA_DATA = Directives::SILO_EXTRA_DATA;

   // todo Fields?
   static constexpr int32 n_cc_components = Directives::NVAR + Directives::NFLUID;
   static constexpr int32 n_fc_components = 1;
   static constexpr int32 n_fc_modes = 1;

   static constexpr int32 n_cc_modes = (Directives::NUMORDER == 2)? 4:
                               (Directives::NUMORDER == 3)? 10:
                               /* numorder == 4 */ 20;

   static constexpr int32 iproblem = (Directives::IPDE_SYSTEM == IPDESystem::MHD)? 1:
                                     (Directives::IPDE_SYSTEM == IPDESystem::CGL)? 2:
                                     /* (Directives::IPDE_SYSTEM == IPDESystem::CGL_p_elec) ? */ 3;

   // todo
   //#if( USE_B0B1_FORMULATION == 1 )
   //    CHARACTER*SILO_VARLEN silo_allvars ( NVAR + NFLUID + 3
   //     1                                   + SILO_EXTRA_DATA)
   //#else
   //    CHARACTER*SILO_VARLEN silo_allvars ( NVAR + NFLUID
   //    1                                   + SILO_EXTRA_DATA)
   //#endif




public: // types:



   using mhd_zc_3 = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;








public: // older version of types: (deprecated)


   //-----------------------------------
   // Globally-defined Quantities.
   //-----------------------------------

   // NB, there is no radial coordinate.
   // The storage is for vertices; there is one more vertex than there are centers.

   // "coordinates_gr_global" has ghost zones filled in and has the coordinates on the
   // unit sphere. I.e., it is 2D in theta-phi space. It is the global mesh (in the
   // theta-phi direction) that we would have obtained on the finest level, if it were
   // fully filled in. The finest level is "ilevel_max".


   using coordinates_gr_global_t = QuadMatrix<R64, Ntp_global+1, Ntp_global+1, 3, 6>;
   using coordinates_gr_global_revt = QuadMatrix<R64, 6, 3, Ntp_global+1, Ntp_global+1>;

   // sector-specific vertex coordinates. (no radial coordinate)
   using coordinates_gr_t = TriMatrix<R64, Ntp_zones+1, Ntp_zones+1, 3>;
   using coordinates_gr_revt = TriMatrix<R64, 3, Ntp_zones+1, Ntp_zones+1>;

   using vertices_t = Matrix<R64, 4, 3>;
   using vertices_revt = Matrix<R64, 3, 4>;


   // These are the radial coordinates. (vertex counting.)
   using rb_save_global_t = Vector<R64, Nr_global+1>;
   // interval width of radial cell (center counting)
   using dr_save_global_t = Vector<R64, Nr_global>;
   // vertex of radial cell (center counting)
   using rc_save_global_t = dr_save_global_t;
   // reversed arrays (same for vectors)
   using rb_save_global_revt = rb_save_global_t;
   using dr_save_global_revt = dr_save_global_t;
   using rc_save_global_revt = rc_save_global_t;


   using cc_mode_t = Vector<N32, n_cc_modes>;
   using cc_mode_revt = cc_mode_t;


   //---------------------------------------------------
   // Geometric Quantities at Zones in a Patch.
   //---------------------------------------------------

   using vert_coords_patch_t = QuadMatrix<R64, Ntp_zones+1, Ntp_zones+1, Nr_zones+1, 3>;
   using vert_coords_patch_revt = QuadMatrix<R64, 3, Nr_zones+1, Ntp_zones+1, Ntp_zones+1>;

   using imask_active_patch_t = TriMatrix<N32, Ntp_zones, Ntp_zones, Nr_zones>;
   using imask_active_patch_revt = TriMatrix<N32, Nr_zones, Ntp_zones, Ntp_zones>;

   using istencil_patch_t = QuadMatrix<N32, Ntp_zones, Ntp_zones, 2, 12>;
   using istencil_patch_revt = QuadMatrix<N32, 12, 2, Ntp_zones, Ntp_zones>;

   // NBNBNB: "f" quantities are face-counted (±) in the target direction, and center-counted in the remaining two directions.
   // Let NA = narea_quad_pts, the number of quadpoints in the logical xy plane.

   // 3xNA-MATRIX AT ± FACE (vector for each quadpoint at ± face)
   using xf_quadpts_patch_t = PentaMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones, 3, NAREA_QUAD_PTS>;
   using yf_quadpts_patch_t = PentaMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones, 3, NAREA_QUAD_PTS>;
   using zf_quadpts_patch_t = PentaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1, 3, NAREA_QUAD_PTS>;
   // NA-VECTOR AT ± FACE (scalar for each quadpoint at ± face)
   using xf_quadwt_patch_t = QuadMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones, NAREA_QUAD_PTS>;
   using yf_quadwt_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones, NAREA_QUAD_PTS>;
   using zf_quadwt_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1, NAREA_QUAD_PTS>;
   // GEOVECTOR AT ± FACE
   // the unit normal in the ± face. All are identical if the face is flat.
   // Under our assumptions, this is the case in the x,y directions but not in z direction.
   using xf_quadnorml_patch_t = QuadMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones, 3>;
   using yf_quadnorml_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones, 3>;
   using zf_quadnorml_patch_t = PentaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1, 3, NAREA_QUAD_PTS>;
   // SCALAR AT ± FACE
   using xf_area_patch_t = TriMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones>;
   using yf_area_patch_t = TriMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones>;
   using zf_area_patch_t = TriMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1>;
   // (another geovector at ± face)
   using xf_centroid_patch_t = xf_quadnorml_patch_t;
   using yf_centroid_patch_t = yf_quadnorml_patch_t;
   using zf_centroid_patch_t = zf_quadnorml_patch_t;

   // ...and, the reverses
   using xf_quadpts_patch_revt = PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yf_quadpts_patch_revt = PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zf_quadpts_patch_revt = PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using xf_quadwt_patch_revt = QuadMatrix<R64, NAREA_QUAD_PTS, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yf_quadwt_patch_revt = QuadMatrix<R64, NAREA_QUAD_PTS, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zf_quadwt_patch_revt = QuadMatrix<R64, NAREA_QUAD_PTS, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using xf_quadnorml_patch_revt = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yf_quadnorml_patch_revt = QuadMatrix<R64, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zf_quadnorml_patch_revt = PentaMatrix<R64, NAREA_QUAD_PTS, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using xf_area_patch_revt = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yf_area_patch_revt = TriMatrix<R64, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zf_area_patch_revt = TriMatrix<R64, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using xf_centroid_patch_revt = xf_quadnorml_patch_revt;
   using yf_centroid_patch_revt = yf_quadnorml_patch_revt;
   using zf_centroid_patch_revt = zf_quadnorml_patch_revt;

   // SCALAR AT CENTER
   using v_scale_length_t = TriMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones>;
   using v_scale_length_revt = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;


   // bfld types
   using xb_quadpts_patch_t = PentaMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones, 3, NBFLD_QUAD_PTS>;
   using yb_quadpts_patch_t = PentaMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones, 3, NBFLD_QUAD_PTS>;
   using zb_quadpts_patch_t = PentaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1, 3, NBFLD_QUAD_PTS>;
   using xb_quadwt_patch_t = QuadMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones, NBFLD_QUAD_PTS>;
   using yb_quadwt_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones, NBFLD_QUAD_PTS>;
   using zb_quadwt_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1, NBFLD_QUAD_PTS>;
   using xb_quadnorml_patch_t = QuadMatrix<R64, Ntp_zones+1, Ntp_zones, Nr_zones, 3>;
   using yb_quadnorml_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones+1, Nr_zones, 3>;
   using zb_quadnorml_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones+1, 3>;

   // ...and, the reverses
   using xb_quadpts_patch_revt = PentaMatrix<R64, NBFLD_QUAD_PTS, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yb_quadpts_patch_revt = PentaMatrix<R64, NBFLD_QUAD_PTS, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zb_quadpts_patch_revt = PentaMatrix<R64, NBFLD_QUAD_PTS, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using xb_quadwt_patch_revt = QuadMatrix<R64, NBFLD_QUAD_PTS, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yb_quadwt_patch_revt = QuadMatrix<R64, NBFLD_QUAD_PTS, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zb_quadwt_patch_revt = QuadMatrix<R64, NBFLD_QUAD_PTS, Nr_zones+1, Ntp_zones, Ntp_zones>;
   using xb_quadnorml_patch_revt = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones+1>;
   using yb_quadnorml_patch_revt = QuadMatrix<R64, 3, Nr_zones, Ntp_zones+1, Ntp_zones>;
   using zb_quadnorml_patch_revt = QuadMatrix<R64, 3, Nr_zones+1, Ntp_zones, Ntp_zones>;

   using v_quadpts_patch_t = PentaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, 3, NUMORDER3>;
   using v_quadwt_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, NUMORDER3>;
   using v_quadtangent_patch_t = HexaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, 3, 3, NUMORDER3>;
   using v_volume_patch_t = TriMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones>;
   using v_centroid_patch_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, 3>;

   // ...and, the reverses
   using v_quadpts_patch_revt = PentaMatrix<R64, NUMORDER3, 3, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_quadwt_patch_revt = QuadMatrix<R64, NUMORDER3, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_quadtangent_patch_revt = HexaMatrix<R64, NUMORDER3, 3, 3, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_volume_patch_revt = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;
   using v_centroid_patch_revt = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones>;


   using Coef_Basis_Moments_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, n_cc_modes>;
   using Coef_Basis_Moments_revt = QuadMatrix<R64, n_cc_modes, Nr_zones, Ntp_zones, Ntp_zones>;

   //------------------------------------------------
   // Physical Quantities at Zones in a Patch.
   //------------------------------------------------

   using zonescalar_t = TriMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones>;
   using zonescalar_revt = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;

   using U_t = PentaMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, n_cc_components, n_cc_modes>;
   using V_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, n_cc_components>;
   using entropygr_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, 4>;
   // Bxf, Byf, Bzf
   using Bf_t = QuadMatrix<R64, Ntp_zones+1, Ntp_zones+1, Nr_zones+1, n_fc_modes>;
   using B0field_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, 3>;
   using gamma_eos_gr_t = TriMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones>;

   // ...and, the reverses
   using U_revt = PentaMatrix<R64, n_cc_modes, n_cc_components, Nr_zones, Ntp_zones, Ntp_zones>;
   using V_revt = QuadMatrix<R64, n_cc_components, Nr_zones, Ntp_zones, Ntp_zones>;
   using entropygr_revt = QuadMatrix<R64, 4, Nr_zones, Ntp_zones, Ntp_zones>;
   // Bxf, Byf, Bzf
   using Bf_revt = QuadMatrix<R64, n_fc_modes, Nr_zones+1, Ntp_zones+1, Ntp_zones+1>;
   using B0field_revt = QuadMatrix<R64, 3, Nr_zones, Ntp_zones, Ntp_zones>;
   using gamma_eos_gr_revt = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;


   using iactive_FromSector_t = TriMatrix<N32, 6, Max_Theoretical_patches_per_sector, 3>;

   // Zone-centered flattener for for all patches.
   using flatten_t = zonescalar_t;
   using flatten_revt = zonescalar_revt;

   // Zone-centered Source for for all patches on this processor.
   using Source_t = QuadMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones, n_cc_components>;
   using Source_revt = QuadMatrix<R64, n_cc_components, Nr_zones, Ntp_zones, Ntp_zones>;

   //
   using CGL_fbeta_gr_t = zonescalar_t;
   using CGL_fbeta_gr_revt = zonescalar_revt;
   using Pelectron_gr_t = zonescalar_t;
   using Pelectron_gr_revt = zonescalar_revt;

   // V_min_local, V_max_local
   using V_minmax_t = Vector<R64, n_cc_components+1>;
   // V_min_local, V_max_local
   using V_minmax_revt = V_minmax_t;




   //-----------------------------------------
    // Physical Quantities at one Zone.
    //-----------------------------------------

    // a fluid vector
    // nvar+fluid: sic (from fortran side, I thought it should be multiplied)
    using FluidVector = Vector<R64, NVAR+NFLUID>;

    // a pair of geovectors: PPar and PPerp (CGL)
    using PParPPerp_t = Vector<GeoVector, 2>;



   // Enforcing radial boundary conditions.
   // Used for PAD_BOUNDARY subroutine.
   using bcarr_r_t = TriMatrix<N32, Ntp_zones, Ntp_zones, 2>;
   using bcarr_r_revt = TriMatrix<N32, 2, Ntp_zones, Ntp_zones>;

   //-----------------------------------------
   // I/O Types (SILO)
   //-----------------------------------------

   // variable used for imaging variables.
   // on the Fortran side, it is a REAL*4 (float)
   using r4gr_t = TriMatrix<R64, Ntp_zones, Ntp_zones, Nr_zones>;
   using r4gr_revt = TriMatrix<R64, Nr_zones, Ntp_zones, Ntp_zones>;


public: // API:


    // These are helpers for the Fortran interface.

    // Recall that a general convention in this class (used above)
    // is that lowercase n means "without ghost positions"
    // and upper-case N means "with ghost positions."

   // Cast an index to an integer (signed).
   int32 n(size_t i) {
      return static_cast<int32>(i);
   }

   // Shift an index back by the number of ghost zones.
   static int32 N(size_t i) {
      return static_cast<int32>(i) - static_cast<int32>(NGhost);
   }

   // These are helpers for indexing.
   // The zone address (SECt, SECp, SECr) is independent
   // of ntp_patches and nr_patches,
   // and is therefore used in testing of the grid.

   // sector-wide theta index, as though
   // the entire sector were a single patch with halo regions.
   static int32 SECt(int32 x, int32 theta) {
      return ntp_zones*x+theta;
   }

   // sector-wide phi index, as though
   // the entire sector were a single patch with halo regions.
   static int32 SECp(int32 y, int32 phi) {
      return ntp_zones*y+phi;
   }

   // sector-wide r index, as though
   // the entire sector were a single patch with halo regions.
   static int32 SECr(int32 z, int32 r) {
      return nr_zones*z+r;
   }

   // return true if the zone address represents
   // a pseudozone created by a point where three sectors meet,
   // i.e. a "tricorner". For example, viewed at any 2d slice,
   // consider the zone in the top left corner of a sector.
   // The zone located at relative address "up, left" is a pseudozone.
   // Whatever is stored in that position is junk data.
   static bool is_tricorner_zone(int32 ix, int32 iy, int32 theta, int32 phi) {
      return (ix == 0 && iy == 0 && theta < NGhost && phi < NGhost) ||
      (ix == 0 && iy == ntp_patches-1 && theta < NGhost && phi >= NGhost+ntp_zones) ||
      (ix == ntp_patches-1 && iy == 0 && theta >= NGhost+ntp_zones && phi < NGhost) ||
      (ix == ntp_patches-1 && iy == ntp_patches-1 && theta >= NGhost+ntp_zones && phi >= NGhost+ntp_zones);
   }

   // Same as `is_tricorner_zone` but used when spoofing
   // one pair (ntp_zones, ntp_patches) with another pair, for testing.
   static bool is_tricorner_zone_spoof(int32 ix, int32 iy, int32 theta, int32 phi, int32 ntp_zones_spoof, int32 ntp_patches_spoof) {
      return (ix == 0 && iy == 0 && theta < NGhost && phi < NGhost) ||
             (ix == 0 && iy == ntp_patches_spoof-1 && theta < NGhost && phi >= NGhost+ntp_zones_spoof) ||
             (ix == ntp_patches_spoof-1 && iy == 0 && theta >= NGhost+ntp_zones_spoof && phi < NGhost) ||
             (ix == ntp_patches_spoof-1 && iy == ntp_patches_spoof-1 && theta >= NGhost+ntp_zones_spoof && phi >= NGhost+ntp_zones_spoof);
   }



};
