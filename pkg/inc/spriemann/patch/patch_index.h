
#pragma once

#include <map>
#include <spriemann/constants/types.h>


/*
 * Geometric specification of
 * the identity of a given patch.
 *
 */
struct PatchIndex {
   int32 isector;
   int32 ilevel;
   int32 ix;
   int32 iy;
   int32 iz;

   // This is needed so that maps can be
   // created with PatchIndex lookup values.
   bool operator<(const PatchIndex & rhs) const {
      // dictionary ordering
      if (isector < rhs.isector) return true;
      else if (ilevel < rhs.ilevel) return true;
      else if (ix < rhs.ix) return true;
      else if (iy < rhs.iy) return true;
      else if (iz < rhs.iz) return true;
      else return false;
   }

   std::string str() const {
      std::string out = std::format("s {} l {} x {} y {} z {}", isector, ilevel, ix, iy, iz);
      return out;
   }

};


// A type for the user who specifies a static AMR mesh.
// The semantic interpretation is the same whether
// the patch index is in the map and set to false, or is not in the map.
// In either case, the patch is statically disabled.
using UserDefinedPatchList = std::map<PatchIndex, bool>;

// default value for no static AMR
inline UserDefinedPatchList no_static_amr;


