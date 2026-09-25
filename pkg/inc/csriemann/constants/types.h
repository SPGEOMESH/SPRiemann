
#pragma once


#include <cstdint>

#include <DGCALC.h>

#include <tennsix/linear.h>
#include <scalar/754.h>
#include <types/integers.h>

// These names are more natural for a Fortran programmer.

using int32 = int32_t;
//using int64 = int64_t;

using real32 = float;
using real64 = double;

// todo incorporate into Tennsix types
using R64 = double754;
using R32 = float754;

// todo GeoVector from Fields
// This improves type semantics
using Coordinate = Vector<R64, 3>;
using GeoVector = Vector<R64, 3>;


enum class IPDESystem {
    MHD,
    CGL,
    CGLp_elec
};

