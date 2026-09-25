
#pragma once

#include <csriemann/constants/types.h>

/*
 * Buffers for halo exchanges.
 *
 * WARNING:
 * MPI calls that occur during halo exchange may use `target_disp`
 * values that depend on the structure of these data structures.
 * Do not edit these data structures without carefully reviewing the
 * implementation of halo exchange in `process.h`.
 *
 */






// Zone-centered buffer
template <
      typename Parameters
>
struct ZCBuffer {

   using P = Parameters;

   static constexpr int32 size_ZCBuffer = P::n_cc_components*P::Nr_zones*P::Ntp_zones*P::NGhost;

   static constexpr int32 size_left_right = P::n_cc_components*P::nr_zones*P::ntp_zones*P::NGhost;
   static constexpr int32 size_down_up = P::n_cc_components*P::nr_zones*P::Ntp_zones*P::NGhost;
   static constexpr int32 size_bottom_top = P::n_cc_components*P::Nr_zones*P::Ntp_zones*P::NGhost;

   using zc_U = std::array<real64, size_ZCBuffer>;

   zc_U zc_buffer;

};




// Face-centered buffer
template <
      typename Parameters
>
struct FCBuffer {

   using P = Parameters;

   static constexpr int32 Ntp_zones = P::Ntp_zones;
   static constexpr int32 Nr_zones = P::Nr_zones;
   static constexpr int32 NGhost = P::NGhost;
   static constexpr int32 n_cc_components = P::n_cc_components;

   using fc_B_xy = std::array<real64, (Nr_zones+1)*(Ntp_zones+1)*(NGhost+1)>;
   using fc_B_z = std::array<real64, (Ntp_zones+1)*(Ntp_zones+1)*(NGhost+1)>;

   fc_B_xy Bxf_r, Bxf_l, Byf_r, Byf_l, Bzf_l, Bzf_r,
         Bxf_u, Bxf_d, Byf_u, Byf_d, Bzf_u, Bzf_d;
   fc_B_z Bxf_t, Bxf_b, Byf_t, Byf_b, Bzf_t, Bzf_b;

};






// Edge-centered buffer
template <
      typename Parameters
>
struct ECBuffer {

   using P = Parameters;

   static constexpr int32 Ntp_zones = P::Ntp_zones;
   static constexpr int32 Nr_zones = P::Nr_zones;
   static constexpr int32 NGhost = P::NGhost;
   static constexpr int32 n_cc_components = P::n_cc_components;

   // todo

};



