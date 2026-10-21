
#pragma once


#include <array>
#include <spriemann/constants/types.h>
#include <spriemann/patch/patch_index.h>


template <int32 ilevel_max>
struct PatchInitializationData {

   std::vector<PatchIndex> indexes;

   // Reminder:
   // the first value is always 0 [sic].
   // the final value is always n_patch [sic].
   std::array<int32, ilevel_max+1> patch_upper;

   int32 n_patch;

};


