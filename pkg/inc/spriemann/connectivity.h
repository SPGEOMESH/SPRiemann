
#pragma once

#include <csriemann/patch/patch_node.h>
#include <csriemann/patch/sector_neighbor.h>
#include <csriemann/patch/patch_initialization_data.h>

#include <kiwi/mpi/sim.h>

#include <map>


/*
 * A data management class for each process
 * to receive patches to process, and to use
 * to lookup information about neighboring patches.
 *
 *
 */
template <
      typename Parameters
>
class Connectivity {
   friend class TConnectivity;
public: // types:

   using P = Parameters;
   using PatchInitializationData = PatchInitializationData<P::ilevel_max>;

public: // constant fields:

   static constexpr int32 ilevel_max = P::ilevel_max;
   static constexpr int32 ntp_patches = P::ntp_patches;
   static constexpr int32 nr_patches = P::nr_patches;

   // The number of patches in each sector
   // (theoretical upper limit, if some are disabled by AMR).
   static constexpr size_t compute_n_patch_sector() {
      size_t out = 0;
      int32 Nxy = ntp_patches;
      int32 Nz = nr_patches;
      for (int32 ilevel = 0; ilevel < ilevel_max; ++ilevel) {
         out += Nxy*Nxy*Nz;
         Nxy *= 2;
         Nz *= 2;
      }
      return out;
   }
   static constexpr size_t n_patch_sector = compute_n_patch_sector();
   static constexpr size_t n_patch_global = n_patch_sector*6;

   // A simple lookup table for sector index of a neighboring sector.
   static constexpr std::array<SectorNeighbor, 6> compute_sector() {
      std::array<SectorNeighbor, 6> out{{
            {5,2,4,1},
            {3,0,5,2},
            {4,1,3,0},
            {1,4,2,5},
            {2,5,0,3},
            {0,3,1,4}
      }};
      // Fortran indexing 1-6: subtract one for C array index values.
      // This notation is equivalent:
//      out[0].L = 5; out[0].R = 2; out[0].D = 4; out[0].U = 1;
//      out[1].L = 3; out[1].R = 0; out[1].D = 5; out[1].U = 2;
//      out[2].L = 4; out[2].R = 1; out[2].D = 3; out[2].U = 0;
//      out[3].L = 1; out[3].R = 4; out[3].D = 2; out[3].U = 5;
//      out[4].L = 2; out[4].R = 5; out[4].D = 0; out[4].U = 3;
//      out[5].L = 0; out[5].R = 3; out[5].D = 1; out[5].U = 4;
      return out;
   }
   static constexpr std::array<SectorNeighbor, 6> sector = compute_sector();

private: // fields:

   // The actual data. The structure is private.
   // The index in this array is internal to the Manifest logic.
   // Processors can use `get` to retrieve a PatchNode.
   std::array<PatchNode, n_patch_global+1> nodelist;

   std::map<ProcessIndex, size_t> lookup;

   int32 Max_patches_per_proc_;

private: // methods:

   size_t get_idx(PatchIndex pindex) const {
      int32 pow = DGCALC::ipow2<int32>(pindex.ilevel);
      int32 Nxy = pow*ntp_patches;
      int32 Nz = pow*nr_patches;
      int32 N2 = Nxy*Nxy;
      size_t out = pindex.ix + Nxy*pindex.iy + N2*pindex.iz + N2*Nz*pindex.ilevel
            + N2*Nz*ilevel_max*pindex.isector;
      return out;
   }

   size_t get_idx_sector(PatchIndex pindex) const {
      int32 pow = DGCALC::ipow2<int32>(pindex.ilevel);
      int32 Nxy = pow*ntp_patches;
      int32 Nz = pow*nr_patches;
      int32 N2 = Nxy*Nxy;
      size_t out = pindex.ix + Nxy*pindex.iy + N2*pindex.iz + N2*Nz*pindex.ilevel;
      return out;
   }

public: // API:

   // Initialize the connectivity list
   // and return the number of patches assigned
   // to the current rank, and other information passed to the process.
   //
   // In order to allow the process (i.e., the caller) to
   // loop over all patches *per level* among the
   // assigned patches, there are `ilevel_max+1` integer
   // values returned. These integers are ordered in increasing
   // order, and limit the range of assigned patches by level.
   // For example: if `ilevel_max` is 3, these 4 integers could be:
   //  0 = base pointer for level 1.
   //                 ...level 1 patches are patch 0...3 inclusive
   //  4 = base pointer for level 2.
   //                 ...level 2 patches are patch 4...8 inclusive
   //  9 = base pointer for level 3.
   //                  ...level 3 patches are patch 9...11 inclusive
   //  12 = n_patch, the number of patches assigned to the process.
   // In order to lessen the risk of confusion or error, the value of `n_patch` is also
   // populated separately in the return type struct `PatchInitializationData`.
   //
   // The final piece of data returned is a list of indices with the length `n_patch`.
   // The caller simply unpacks these and uses them to initialize the patches themselves.
   PatchInitializationData initialize(
         int32 world_size = 1,
         int32 world_rank = 0,
         int32 sector_color = 0,
         int32 sector_size = 1,
         int32 sector_rank = 0,
         int32 sector_leader = 0,
         MPI_Comm sector_comm = MPI_COMM_NULL,
         MPI_Comm world_comm = MPI_COMM_NULL,
         UserDefinedPatchList& existinglist = no_static_amr
   ) {

      PatchInitializationData out;

      // The steps are:

      // (1) As the patches are generated (with isector, ix, iy, iz, ilevel)
      // they can be pushed on `nodelist`. All their data can be populated now.
      // (We have only `exists` as data currently.) If static AMR, use the
      // `existinglist`, otherwise set existing only for level 1.

      // (2) Step through the entire set of PatchNode objects that were
      // just generated, and set each pointer to correctly point to another
      // PatchNode - this requires geometric reasoning.

      // (3) Push PatchIndex objects into the `lookup` map and ensure
      // that the values point to *that* patch.

      // (4) Assign patches to processes, and build the return value for
      // this function.

      // Done.

      // Cognitively speaking, these are four "steps", but it is most
      // natural to do all the work in a single pass.

      // set the out-of-domain node. This way the program does not crash
      // if links are checked for diagnostic purposes.
      {
         PatchNode node;
         // Dummy index values - if you see negative index values
         // then you have read this node! The values are all different
         // so you know which position you've read, in case it is useful.
         node.index = {-1, -4, -5, -6, -7};
         node.status = PatchStatus::inactive;
         node.world_rank = -2;
         node.sector_rank = -3;
         node.L = &nodelist[n_patch_global];
         node.R = &nodelist[n_patch_global];
         node.D = &nodelist[n_patch_global];
         node.U = &nodelist[n_patch_global];
         node.B = &nodelist[n_patch_global];
         node.T = &nodelist[n_patch_global];
         nodelist[n_patch_global] = node;
      }

      out.patch_upper[0] = 0;
      int32 upperi = 0;
      int32 proc_ii = 0;
      int32 prev_ilevel = 0;
      int32 idx = 0;
      for (int32 isector = 0; isector < 6; ++isector) {
         int32 Nxy = ntp_patches;
         int32 Nz = nr_patches;
         int32 N2 = Nxy*Nxy;
         for (int32 ilevel = 0; ilevel < ilevel_max; ++ilevel) {
            for (int32 iz = 0; iz < Nz; ++iz) {
               for (int32 iy = 0; iy < Nxy; ++iy) {
                  for (int32 ix = 0; ix < Nxy; ++ix) {
                     PatchNode node;
                     node.index = {isector, ilevel, ix, iy, iz};
                     node.world_rank = 0;
                     node.sector_rank = 0;
                     node.status = (ilevel == 0)? PatchStatus::active: PatchStatus::inactive;
                     if (!existinglist.empty()) {
                        if (ilevel != 0) {
                           if (existinglist.contains(node.index)) {
                              node.status = PatchStatus::active;
                           }
                           else {
                              node.status = PatchStatus::inactive;
                           }
                        }
                     }
                     // set neighbor pointers.
                     if constexpr (ntp_patches == 1) {
                        // both L and R out of sector
                        node.L = &nodelist[get_idx({sector[isector].L, ilevel, 0, Nxy-1-iy, iz})];
                        node.R = &nodelist[get_idx({sector[isector].R, ilevel, iy, Nxy-1, iz})];
                     }
                     else {
                        if (ix == 0) {
                           // L out of sector
                           node.L = &nodelist[get_idx({sector[isector].L, ilevel, 0, Nxy-1-iy, iz})];
                           node.R = &nodelist[idx+1];
                        }
                        else if (ix == Nxy-1) {
                           // R out of sector
                           node.L = &nodelist[idx-1];
                           node.R = &nodelist[get_idx({sector[isector].R, ilevel, iy, Nxy-1, iz})];
                        }
                        else {
                           node.L = &nodelist[idx-1];
                           node.R = &nodelist[idx+1];
                        }
                     }
                     if constexpr (ntp_patches == 1) {
                        // both D and U out of sector
                        node.D = &nodelist[get_idx({sector[isector].D, ilevel, Nxy-1-ix, 0, iz})];
                        node.U = &nodelist[get_idx({sector[isector].U, ilevel, Nxy-1, ix, iz})];
                     }
                     else {
                        if (iy == 0) {
                           // D out of sector
                           node.D = &nodelist[get_idx({sector[isector].D, ilevel, Nxy-1-ix, 0, iz})];
                           node.U = &nodelist[idx + Nxy];
                        }
                        else if (iy == Nxy-1) {
                           // U out of sector
                           node.D = &nodelist[idx - Nxy];
                           node.U = &nodelist[get_idx({sector[isector].U, ilevel, Nxy-1, ix, iz})];
                        }
                        else {
                           node.D = &nodelist[idx - Nxy];
                           node.U = &nodelist[idx + Nxy];
                        }
                     }
                     if constexpr (nr_patches == 1) {
                        // both B and T out of sector
                        node.B = &nodelist[n_patch_global];
                        node.T = &nodelist[n_patch_global];
                     }
                     else {
                        if (iz == 0) {
                           // B out of domain
                           node.B = &nodelist[n_patch_global];
                           node.T = &nodelist[idx + N2];
                        } else if (iz == Nz-1) {
                           // T out of domain
                           node.B = &nodelist[idx - N2];
                           node.T = &nodelist[n_patch_global];
                        }
                        else {
                           node.B = &nodelist[idx - N2];
                           node.T = &nodelist[idx + N2];
                        }
                     }


                     // Globally-defined information:
                     int32 k = get_idx_sector(node.index);
                     int32 m = sector_size;
                     node.sector_rank = k%m;
                     // a sector_rank-specific section of code:
                     if (node.index.isector == sector_color%world_size
                        && node.sector_rank == sector_rank) {
                        // this patch is assigned to me.
                        node.world_rank = world_rank;
                        node.proc_ii = proc_ii;
                        out.indexes.push_back(node.index);
                        if (prev_ilevel != node.index.ilevel) {
                           out.patch_upper[upperi] = proc_ii;
                           upperi++;
                           prev_ilevel = node.index.ilevel;
                        }
                        proc_ii++;
                     }
                     // For our convenience and sanity, we now ensure that
                     // the globally visible connectivity data structure contains
                     // all the information we could possibly need about who owns which patch.
                     int is_initialized = 0;
                     MPI_Initialized(&is_initialized);
                     // Gated so that we can run this method without MPI during testing.
                     if (is_initialized) {
                        // First ensure that the process with world_rank isector has the information.
                        MPI_Bcast(&node.world_rank, 1, MPI_INT32_T, node.sector_rank, sector_comm);
                        MPI_Bcast(&node.proc_ii, 1, MPI_INT32_T, node.sector_rank, sector_comm);
                        // Send the information to everyone.
                        MPI_Bcast(&node.world_rank, 1, MPI_INT32_T, isector, world_comm);
                        MPI_Bcast(&node.proc_ii, 1, MPI_INT32_T, isector, world_comm);
                     }
                     // write node to nodelist
                     nodelist[idx] = node;
                     lookup[{node.world_rank, node.proc_ii}] = idx;
                     // CHECK:
//                     if (idx != get_idx({isector, ilevel, ix, iy, iz})) {
//                        cout << std::format("[Connectivity][Warning] idx {} idx' {}\n", idx, get_idx({isector, ilevel, ix, iy, iz}));
//                     }
                     idx++;
                  }
               }
            }
            // next level: the number of patches is doubled.
            Nxy *= 2;
            Nz *= 2;
            N2 = Nxy*Nxy;
         }
      }
      // Now the number of patches assigned is simple.
      out.n_patch = out.indexes.size();
      out.patch_upper[ilevel_max] = out.n_patch;
      return out;
   }

   // Get a PatchNode data structure
   // from a ProcessIndex data structure.
   PatchNode get(const ProcessIndex& processindex) const {
      auto idx = lookup[processindex];
      return nodelist[idx];
   }

   // Get a PatchNode data structure
   // from a PatchIndex data structure.
   PatchNode get(const PatchIndex& patchindex) const {
      auto idx = get_idx(patchindex);
      return nodelist[idx];
   }

   // Get a PatchNode data structure
   // from a simple linearly ordered patch index.
   PatchNode get(int32 ipatch) const {
      return nodelist[ipatch];
   }

   // Get a PatchNode data structure
   // from a PatchNode pointer.
   PatchNode get(PatchNode* patchnode) const {
      return *patchnode;
   }

   int32 Max_patches_per_proc() const {
      return Max_patches_per_proc_;
   }

   bool inside_domain(PatchNode* node_p) const {
      return (node_p != &nodelist[n_patch_global]);
   }

   std::string str() const {
      std::string out = std::format("n_patch_global {}\nMax_patches_per_proc {}\n", n_patch_global, Max_patches_per_proc());
      for (int32 idx = 0; idx < n_patch_global; ++idx) {
         out += nodelist[idx].str() + "\n";
      }
      return out;
   }

};


