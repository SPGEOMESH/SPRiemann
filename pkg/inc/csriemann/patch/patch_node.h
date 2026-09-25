
#pragma once

#include <csriemann/constants/types.h>
#include <csriemann/patch/patch_index.h>
#include <csriemann/patch/patch_status.h>

/*
 * Data needed about patch,
 * to be kept in the manifest.
 */
struct PatchNode {
   // What is the patch's index?
   // (Needed for lookup of neighbor patches.)
   PatchIndex index;
   // rank in the global communicator of the
   // proces assigned this patch.
   int32 world_rank;
   // the patch-counting index in the patch-list
   // of the process who was assigned this patch.
   int32 proc_ii;
   // rank in the sector communicator of the
   // process assigned this patch.
   int32 sector_rank;
   // Status, or mutable state during simulation
   PatchStatus status;

   // pointers
   PatchNode* L;
   PatchNode* R;
   PatchNode* U;
   PatchNode* D;
   PatchNode* T;
   PatchNode* B;

   std::string str() const {
      std::string out = std::format("index {}\n", index.str());
      out += std::format("world_rank {} proc_ii {} sector_rank {}\n", world_rank, proc_ii, sector_rank);
      out += std::format(" status {}\n", (status == PatchStatus::active)? "active": "inactive");
      return out;
   }

};



/*
 * Alternative identifier of a
 * patch, convenient for a processor
 * as an entry point to the manifest.
 *
 * It is simple, we wrap it up in a struct
 * just so we can use it in a std::map
 * as key for key-value pairs.
 *
 *
 */
struct ProcessIndex {
   // the processor's global rank.
   int32 world_rank;
   // the index in the processor's own list.
   int32 ii;

   // This is needed so that maps can be
   // created with ProcessIndex lookup values.
   bool operator<(const ProcessIndex& rhs) const {
      // dictionary ordering
      if (world_rank < rhs.world_rank) return true;
      else if (ii < rhs.ii) return true;
      else return false;
   }

};


