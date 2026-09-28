
#pragma once

#include <csriemann/patch/patch.h>
#include <csriemann/connectivity.h>
#include <csriemann/mesh.h>
#include <csriemann/patch/buffers.h>

#include <csriemann/extern_c/riemann_integrate.h>

#include <kiwi/mpi/sim.h>

#include <iostream>
#include <unistd.h>


/*
 * Patches that exist in the current global
 * rank's memory address space. Patches
 * that do not live here are accessed via
 * MPI messaging to/from other processes,
 * in a flat, global-communicator communication
 * hierarchy.
 *
 *
 */
template <
      typename Parameters
>
class Process {
private: // types:
   using P = Parameters;
public: // types:
   using Mesh = Mesh<P>;
   using Connectivity = Connectivity<P>;
   using Patch = Patch<P>;

   using ZCBuffer = ZCBuffer<P>;
   using zc_U = ZCBuffer::zc_U;

   using FCBuffer = FCBuffer<P>;
   using fc_B_xy = FCBuffer::fc_B_xy;
   using fc_B_z = FCBuffer::fc_B_z;

   using ECBuffer = ECBuffer<P>;
   // todo

   // deprecated
//   static_assert(sizeof(ZCBuffer) == sizeof(zc_U_xy)+sizeof(zc_U_z), "ZCBuffer is not contiguous on this architecture.");
   static_assert(sizeof(FCBuffer) == sizeof(fc_B_xy)*6*2+sizeof(fc_B_z)*6, "FCBuffer is not contiguous on this architecture.");
   static_assert(sizeof(ECBuffer) == 1, "ECBuffer is not contiguous on this architecture.");

   static constexpr int32 NUMORDER = P::NUMORDER;
   static constexpr int32 ilevel_max = P::ilevel_max;

public: // fields:

   // error codes for MPI calls
   int32 ierr;

   // world comm size == number of images.
   // rank is processor rank in world comm.
   int32 world_size;
   int32 world_rank;

   // size and rank for the sector communicator.
   int32 sector_size;
   int32 sector_rank;

   // The number of patches managed by this processor.
   int32 n_patch;

   // The window to buffer memory within the sector.
   MPI_Win sector_win;

   // global communicator
   MPI_Comm world_comm;
   // sector-wide communicator
   MPI_Comm sector_comm;

   // The global mesh, used to populate patch geometry.
   std::unique_ptr<Mesh> global_mesh;
   // Patches assigned to this processor.
   std::unique_ptr<Patch[]> patches;
   // sector-internal buffer for halo exchange: recv (receive)
   std::unique_ptr<ZCBuffer[]> zc_recv_buf;
   // buffers for halo exchange: send and put (inter-sector and intra-sector, respectively)
   std::unique_ptr<zc_U[]> zc_send_buf;
   std::unique_ptr<zc_U[]> zc_put_buf;
   // Patch connectivity, i.e., lists with indexing
   // and connectivity data. Size scales with global number of patches.
   std::unique_ptr<Connectivity> c_ptr;
   Connectivity& connectivity;


   std::array<int32, ilevel_max+1> patch_upper;


   // sectors assigned to this process
   int32 sector_color;

   // Record of current time of simulation on this processor.
   real64 time;

   // Reduced minimum of next predicted
   // timestep from all patches on this processor.
   real64 dtcurp1;

   // Reduced min and max variables from
   // all patches on this processor.
   real64 V_min;
   real64 V_max;

public: // ctor:

   /*
    * Create and initialize a processor, the worker/agent
    * for all processing in the simulation.
    */
   Process(
         UserDefinedPatchList& existinglist = no_static_amr
   ):
      ierr(0),
      world_size(0),
      world_rank(0),
      sector_size(0),
      sector_rank(0),
      n_patch(0),
      sector_win(MPI_Win()),
      sector_comm(MPI_Comm()),
      c_ptr(std::make_unique<Connectivity>()),
      connectivity(*c_ptr)
   {
      int is_initialized = 0;
      MPI_Initialized(&is_initialized);
      if (!is_initialized) {
         std::cout << "Creating a Process before initializing MPI is not permitted." << endl;
         exit(-1);
      }
      // nimages == comm_world_size TODO
      ierr = MPI_Comm_dup(MPI_COMM_WORLD, &world_comm);
      ierr = MPI_Comm_size(world_comm, &world_size);
      ierr = MPI_Comm_rank(world_comm, &world_rank);

      // There are always 6 sectors.
      if (world_size >= 6) {
         // Realistic case. The sector color simply assigns each process to a sector.
         sector_color = world_rank % 6;
      }
      else {
         // All other cases.
         // In this case, communicators must be obtained based on sector_color%world_size,
         // not sector_color, as in `if (isector == sector_color%world_size) {use_communicator();}`
         sector_color = world_rank;
      }
      MPI_Comm_split(world_comm, sector_color, world_rank, &sector_comm);
      MPI_Comm_size(sector_comm, &sector_size);
      MPI_Comm_rank(sector_comm, &sector_rank);
      // The sector communicators are defined.

      int32 sector_leader;
      if (sector_rank == 0) {
         sector_leader = world_rank;
      }
      MPI_Bcast(&sector_leader, 1, MPI_INT32_T, 0, sector_comm);

      PatchInitializationData pid = connectivity.initialize(
            world_size,
            world_rank,
            sector_color,
            sector_size,
            sector_rank,
            sector_leader,
            sector_comm,
            world_comm,
            existinglist
      );

      // The expected case is that there is more work than available compute.
      // That's why we are parallelizing the code in the first place.
      // If this happens to not be true for any reason,
      // then the parallelization procedure is allowed, but it is inefficient.
      // In this case, warn the user.
      // Note: if world_size < 6, then this warning will never be triggered
      // because in that case sector_size == 1 and n_patch_sector >= 1.
      int32 n_patch_sector = connectivity.n_patch_sector;
      if (n_patch_sector < sector_size) {
         print_warning(std::format("There are {} patches in sector {}, "
                                   "but {} processes have been assigned to that sector. "
                                   "Some processes have been assigned zero patches.", n_patch_sector, sector_color, sector_size));
      }

      // Unpack the Patch Initialization Data.
      // set the n_patches once and for all.
      // This is the actual number of patches this process has been assigned.
      n_patch = pid.n_patch;
      // copy in the beginning and ending indices
      // for level-by-level looping.
      patch_upper = pid.patch_upper;
      // create the patch memory.
      patches = std::make_unique<Patch[]>(n_patch);
      for (int ii = 0; ii < n_patch; ++ii) {
         // pass the indexes to the patches
         // so they can be initialized.
         patches[ii].index = pid.indexes[ii];
      }

      // create the buffer memory for halo exchanges.
      zc_recv_buf = std::make_unique<ZCBuffer[]>(n_patch);
      zc_send_buf = std::make_unique<zc_U[]>(n_patch);
      zc_put_buf = std::make_unique<zc_U[]>(n_patch);
      MPI_Aint size_sector_buffers = sizeof(ZCBuffer)*n_patch;
      MPI_Win_create(zc_recv_buf.get(), size_sector_buffers, 1, MPI_INFO_NULL, sector_comm, &sector_win);
      // The intra-sector communication window is defined (buffers for communication internal to sector).

      // We end up never needing this, but print it as a quiet update.
      int32 Max_patches_per_proc;
      MPI_Allreduce(&n_patch, &Max_patches_per_proc, 1, MPI_INT32_T, MPI_MAX, world_comm);
      print_message(std::format("Max_patches_per_proc {}", Max_patches_per_proc));

      // imaging (I/O)
//      initialize_Imaging_class();
   }

public:

   ~Process() {
      MPI_Win_free(&sector_win);
   }

public: // API:

   void initialize_geometry(
         real64 rmin,
         real64 rmax,
         real64 radial_mesh_ratio
   ) {
      global_mesh = std::make_unique<Mesh>(rmin, rmax, radial_mesh_ratio);
      for (int ii = 0; ii < n_patch; ++ii) {
         PatchIndex idx = patches[ii].index;
         patches[ii] = Patch(
               *global_mesh,
               idx.isector,
               idx.ix,
               idx.iy,
               idx.iz,
               idx.ilevel
         );
      }

      // todo review the function of this.
//      build_base_patches();

      // todo this is performed in a separate proc loop - is that required??
//      is_bottom_or_top();

      // TODO review and test
   }

   void main_loop(int32 rki) {
      for (int ii = 0; ii < n_patch; ++ii) {
         patches[ii].main_logic(rki, ii);
      }
   }

   void defect_loop1() {
      for (int ii = 0; ii < n_patch; ++ii) {
         patches[ii].defect_logic1();
      }
   }

   void defect_loop2() {
      for (int ii = 0; ii < n_patch; ++ii) {
         patches[ii].defect_logic2();
      }
   }


public: // communication:

   // The main halo exchange procedure. Warning: there is no hope of
   // comprehending the logic of the indexing that appears here from the
   // source code alone, without geometric reasoning and knowledge of the
   // choices made in defining the sector relationships.
   // Consult documentation, or seek assistance from project developers/maintainers.
   void halo_exchange() {
      // Open an active synchronization epoch for intra-sector RMA.
      MPI_Win_fence(0, sector_win);
      // messages marked "step_communication" are only printed when
      // the `debug_communication` flag has been set in Directives.
      step_communication("Starting halo exchange.", 0);
      load_local_buffers_UX_step1();
      load_local_buffers_UY_step2_step3();
      load_local_buffers_UZ_step4();
      step_communication("Finished halo exchange.", 0);
   }

   // Step 1 establishes ghost zones in left and right halo regions.
   void load_local_buffers_UX_step1() {
      halo_exchange_send_left(); halo_exchange_writeback_left();
      halo_exchange_send_right(); halo_exchange_writeback_right();
   }

   // Step 2 establishes ghost zones in down and up halo regions.
   // Step 3 establishes ghost zones that are missed during up exchange.
   void load_local_buffers_UY_step2_step3() {
      halo_exchange_send_down(); halo_exchange_writeback_down();
      halo_exchange_send_up(); halo_exchange_writeback_up();
   }

   // Step 4 establishes ghost zones in bottom and top halo regions.
   void load_local_buffers_UZ_step4() {
      halo_exchange_send_bottom(); halo_exchange_writeback_bottom();
      halo_exchange_send_top(); halo_exchange_writeback_top();
   }

   // The halo exchange operation in which data is sent to the left neighbor.
   void halo_exchange_send_left() {
      std::vector<MPI_Request> reqs(2*n_patch);
      int32 reqs_used = 0;

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto up_node = connectivity.get(node.U);
         auto left_node = connectivity.get(node.L);
         if (up_node.index.isector != node.index.isector) {
            size_t idx = 0;
            // Ordinary loop order is zyx. This loop is zxy.
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t i = P::NGhost; i < P::NGhost + P::ntp_zones; ++i) {
                  for (int32 j = P::ntp_zones+P::NGhost-1; j >= P::ntp_zones; --j) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto j_ = static_cast<size_t>(j);
                        auto x = patches[ii].U[{cci, k}][{j_, i}].under();
                        zc_send_buf[ii][idx++] = x;
                     }
                  }
               }
            }
         }
         if (left_node.index.isector == node.index.isector) {
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::NGhost; j < P::NGhost + P::ntp_zones; ++j) {
                  for (size_t i = P::NGhost; i < P::NGhost + P::NGhost; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = patches[ii].U[{cci, k}][{j, i}].under();
                        zc_put_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            step_communication(std::format("[{}] → [win:{}:{}]", world_rank, left_node.sector_rank, left_node.proc_ii));
            // Put zc_put_buf into zc_recv_buf.
            MPI_Put(
                  zc_put_buf[ii].data(),
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  left_node.sector_rank,
                  left_node.proc_ii*sizeof(ZCBuffer),
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  sector_win
            );
         }
      }

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto right_node = connectivity.get(node.R);
         auto up_node = connectivity.get(node.U);
         // multiple patches along the boundary may belong
         // to the same process, therefore tag the payload.
         // Use the receiver's ii integer.
         if (right_node.index.isector != node.index.isector) {
            int32 tag = ii;
            // Receive data into zc_recv_buf.
            step_communication(std::format("{} → [{}] (tag {})", right_node.world_rank, world_rank, tag));
            MPI_Irecv(
                  &zc_recv_buf[ii].zc_buffer[0],
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  right_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
         }
         if (up_node.index.isector != node.index.isector) {
            int32 tag = up_node.proc_ii;
            // Send zc_send_buf.
            step_communication(std::format("[{}] → {} (tag {})", world_rank, up_node.world_rank, tag));
            MPI_Isend(
                  zc_send_buf[ii].data(),
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  up_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
         }
      }
      // Wait until all receive and send tasks are retired from the network.
      MPI_Waitall(reqs_used, reqs.data(), MPI_STATUSES_IGNORE);
      MPI_Win_fence(0, sector_win);

   }

   // The halo exchange operation in which data is sent to the right neighbor.
   void halo_exchange_send_right() {
      std::vector<MPI_Request> reqs(2*n_patch);
      int32 reqs_used = 0;

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto left_node = connectivity.get(node.L);
         auto right_node = connectivity.get(node.R);
         if (left_node.index.isector != node.index.isector) {
            size_t idx = 0;
            // Ordinary loop order is zyx.
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (int32 j = P::NGhost + P::ntp_zones - 1; j >= P::NGhost; --j) {
                  for (int32 i = P::NGhost+P::NGhost-1; i >= P::NGhost; --i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto i_ = static_cast<size_t>(i);
                        auto j_ = static_cast<size_t>(j);
                        auto x = patches[ii].U[{cci, k}][{j_, i_}].under();
                        zc_send_buf[ii][idx++] = x;
                     }
                  }
               }
            }
         }
         if (right_node.index.isector == node.index.isector)  {
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::NGhost; j < P::NGhost + P::ntp_zones; ++j) {
                  for (size_t i = P::ntp_zones; i < P::ntp_zones + P::NGhost; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = patches[ii].U[{cci, k}][{j, i}].under();
                        zc_put_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            step_communication(std::format("[{}] → [win:{}:{}]", world_rank, right_node.sector_rank, right_node.proc_ii));
            // Put zc_put_buf into zc_recv_buf.
            MPI_Put(
                  zc_put_buf[ii].data(),
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  right_node.sector_rank,
                  right_node.proc_ii*sizeof(ZCBuffer),
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  sector_win
            );
         }
      }

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto left_node = connectivity.get(node.L);
         // multiple patches along the boundary may belong
         // to the same process, therefore tag the payload.
         // Use the receiver's ii integer.
         if (left_node.index.isector != node.index.isector) {
            int32 tag = ii;
            // Receive data into zc_recv_buf.
            step_communication(std::format("{} → [{}] (tag {})", left_node.world_rank, world_rank, tag));
            MPI_Irecv(
                  &zc_recv_buf[ii].zc_buffer[0],
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  left_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
         }
         if (left_node.index.isector != node.index.isector) {
            int32 tag = left_node.proc_ii;
            // Send zc_send_buf.
            step_communication(std::format("[{}] → {} (tag {})", world_rank, left_node.world_rank, tag));
            MPI_Isend(
                  zc_send_buf[ii].data(),
                  ZCBuffer::size_left_right,
                  MPI_DOUBLE,
                  left_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
         }
      }
      // Wait until all receive and send tasks are retired from the network.
      MPI_Waitall(reqs_used, reqs.data(), MPI_STATUSES_IGNORE);
      MPI_Win_fence(0, sector_win);

   }

   // The halo exchange operation in which data is sent to the down neighbor.
   void halo_exchange_send_down() {
      std::vector<MPI_Request> reqs(2*3*n_patch);
      int32 reqs_used = 0;

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto right_node = connectivity.get(node.R);
         auto down_node = connectivity.get(node.D);
         if (right_node.index.isector != node.index.isector) {
            // we cannot use the overlapping (OL) communication pattern,
            // because the orientation has not shifted to up-down along the sector boundary.
            size_t idx = 0;
            // Ordinary loop order is zyx. This loop is zxy.
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (int32 i = P::ntp_zones+P::NGhost-1; i >= P::ntp_zones; --i) {
                  for (size_t j = P::NGhost; j < P::NGhost + P::ntp_zones; ++j) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto i_ = static_cast<size_t>(i);
                        auto x = patches[ii].U[{cci, k}][{j, i_}].under();
                        zc_send_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            // up diagonal
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (int32 i = P::ntp_zones+P::NGhost-1; i >= P::ntp_zones; --i) {
                  for (size_t j = P::ntp_zones; j < P::NGhost + P::ntp_zones; ++j) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto i_ = static_cast<size_t>(i);
                        auto x = patches[ii].U[{cci, k}][{j, i_}].under();
                        zc_send_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            // down diagonal
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (int32 i = P::ntp_zones+P::NGhost-1; i >= P::ntp_zones; --i) {
                  for (size_t j = P::NGhost; j < P::NGhost + P::NGhost; ++j) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto i_ = static_cast<size_t>(i);
                        auto x = patches[ii].U[{cci, k}][{j, i_}].under();
                        step_communication(std::format("i {} j {} x {}", i_, j, x), 4);
                        zc_send_buf[ii][idx++] = x;
                     }
                  }
               }
            }
         }
         if (down_node.index.isector == node.index.isector) {
            // we can use the overlapping (OL) communication pattern.
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::NGhost; j < P::NGhost + P::NGhost; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = patches[ii].U[{cci, k}][{j, i}].under();
                        zc_put_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            step_communication(std::format("[{}] → [win:{}:{}]", world_rank, down_node.sector_rank, down_node.proc_ii));
            // Put zc_put_buf into zc_recv_buf.
            MPI_Put(
                  zc_put_buf[ii].data(),
                  ZCBuffer::size_down_up,
                  MPI_DOUBLE,
                  down_node.sector_rank,
                  down_node.proc_ii*sizeof(ZCBuffer),
                  ZCBuffer::size_down_up,
                  MPI_DOUBLE,
                  sector_win
            );
         }
      }


      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto right_node = connectivity.get(node.R);
         auto up_node = connectivity.get(node.U);
         // multiple patches along the boundary may belong
         // to the same process, therefore tag the payload.
         // Use the receiver's ii integer.
         if (up_node.index.isector != node.index.isector) {
            int32 tag = 3*ii;
            // Receive data into zc_recv_buf.
            step_communication(std::format("{} → [{}] (tag {})", up_node.world_rank, world_rank, tag));
            MPI_Irecv(
                  &zc_recv_buf[ii].zc_buffer[0],
                  P::n_cc_components*P::nr_zones*P::ntp_zones*P::NGhost,
                  MPI_DOUBLE,
                  up_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
            // Tricorner-adjacent patches have no partner.
            if (!((node.index.ix == 0) && (node.index.iy == P::ntp_patches-1))) {
               auto up_down_node = connectivity.get(up_node.D);
               int32 up_diagonal_tag = 3*ii+1;
               step_communication(std::format("{} → [{}] (tag {})", up_down_node.world_rank, world_rank, up_diagonal_tag));
               MPI_Irecv(
                     &zc_recv_buf[ii].zc_buffer[P::n_cc_components*P::nr_zones*P::ntp_zones*P::NGhost],
                     P::n_cc_components*P::nr_zones*P::NGhost*P::NGhost,
                     MPI_DOUBLE,
                     up_down_node.world_rank,
                     up_diagonal_tag,
                     world_comm,
                     &reqs[reqs_used]
               );
               ++reqs_used;
            }
            if (!((node.index.ix == P::ntp_patches-1) && (node.index.iy == P::ntp_patches-1))) {
               auto up_up_node = connectivity.get(up_node.U);
               int32 down_diagonal_tag = 3*ii+2;
               step_communication(std::format("{} → [{}] (tag {})", up_up_node.world_rank, world_rank, down_diagonal_tag));
               MPI_Irecv(
                     &zc_recv_buf[ii].zc_buffer[P::n_cc_components*P::nr_zones*(P::ntp_zones+P::NGhost)*P::NGhost],
                     P::n_cc_components*P::nr_zones*P::NGhost*P::NGhost,
                     MPI_DOUBLE,
                     up_up_node.world_rank,
                     down_diagonal_tag,
                     world_comm,
                     &reqs[reqs_used]
               );
               ++reqs_used;
            }
         }
         if (right_node.index.isector != node.index.isector) {
            int32 tag = 3*right_node.proc_ii;
            // Send zc_send_buf.
            step_communication(std::format("[{}] → {} (tag {})", world_rank, right_node.world_rank, tag));
            MPI_Isend(
                  &zc_send_buf[ii][0],
                  P::n_cc_components*P::nr_zones*P::ntp_zones*P::NGhost,
                  MPI_DOUBLE,
                  right_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
            // Tricorner-adjacent patches have no partner.
            if (!((node.index.ix == P::ntp_patches-1) && (node.index.iy == P::ntp_patches-1))) {
               auto right_right_node = connectivity.get(right_node.R);
               int32 up_diagonal_tag = 3*right_right_node.proc_ii + 1;
               step_communication(std::format("[{}] → {} (tag {})", world_rank, right_right_node.world_rank, up_diagonal_tag));
               MPI_Isend(
                     &zc_send_buf[ii][P::n_cc_components * P::nr_zones * P::ntp_zones * P::NGhost],
                     P::n_cc_components * P::nr_zones * P::NGhost * P::NGhost,
                     MPI_DOUBLE,
                     right_right_node.world_rank,
                     up_diagonal_tag,
                     world_comm,
                     &reqs[reqs_used]
               );
               ++reqs_used;
            }
            if (!((node.index.ix == P::ntp_patches-1) && (node.index.iy == 0))) {
               auto right_left_node = connectivity.get(right_node.L);
               int32 down_diagonal_tag = 3*right_left_node.proc_ii + 2;
               step_communication(std::format("[{}] → {} (tag {})", world_rank, right_left_node.world_rank, down_diagonal_tag));
               MPI_Isend(
                     &zc_send_buf[ii][P::n_cc_components * P::nr_zones * P::ntp_zones * P::NGhost+P::n_cc_components * P::nr_zones * P::NGhost * P::NGhost],
                     P::n_cc_components*P::nr_zones*P::NGhost*P::NGhost,
                     MPI_DOUBLE,
                     right_left_node.world_rank,
                     down_diagonal_tag,
                     world_comm,
                     &reqs[reqs_used]
               );
               ++reqs_used;
            }
         }
      }
      // Wait until all receive and send tasks are retired from the network.
      MPI_Waitall(reqs_used, reqs.data(), MPI_STATUSES_IGNORE);
      MPI_Win_fence(0, sector_win);

   }

   // The halo exchange operation in which data is sent to the up neighbor.
   void halo_exchange_send_up() {
      std::vector<MPI_Request> reqs(2*n_patch);
      int32 reqs_used = 0;

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto down_node = connectivity.get(node.D);
         auto up_node = connectivity.get(node.U);
         if (down_node.index.isector != node.index.isector) {
            size_t idx = 0;
            // Ordinary loop order is zyx.
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (int32 j = 2 * P::NGhost-1; j >= P::NGhost; --j) {
                  for (int32 i = P::Ntp_zones-1; i >= 0; --i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto i_ = static_cast<size_t>(i);
                        auto j_ = static_cast<size_t>(j);
                        auto x = patches[ii].U[{cci, k}][{j_, i_}].under();
                        zc_send_buf[ii][idx++] = x;
                     }
                  }
               }
            }
         }
         if (up_node.index.isector == node.index.isector) {
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::ntp_zones; j < P::ntp_zones + P::NGhost; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = patches[ii].U[{cci, k}][{j, i}].under();
                        zc_put_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            step_communication(std::format("[{}] → [win:{}:{}]", world_rank, up_node.sector_rank, up_node.proc_ii));
            // Put zc_put_buf into zc_recv_buf.
            MPI_Put(
                  zc_put_buf[ii].data(),
                  ZCBuffer::size_down_up,
                  MPI_DOUBLE,
                  up_node.sector_rank,
                  up_node.proc_ii*sizeof(ZCBuffer),
                  ZCBuffer::size_down_up,
                  MPI_DOUBLE,
                  sector_win
            );
         }
      }

      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         auto down_node = connectivity.get(node.D);
         // multiple patches along the boundary may belong
         // to the same process, therefore tag the payload.
         // Use the receiver's ii integer.
         if (down_node.index.isector != node.index.isector) {
            int32 tag = ii;
            // Receive data into zc_recv_buf.
            step_communication(std::format("{} → [{}] (tag {})", down_node.world_rank, world_rank, tag));
            MPI_Irecv(
                  &zc_recv_buf[ii].zc_buffer[0],
                  ZCBuffer::size_down_up,
                  MPI_DOUBLE,
                  down_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
         }
         if (down_node.index.isector != node.index.isector) {
            int32 tag = down_node.proc_ii;
            // Send zc_send_buf.
            step_communication(std::format("[{}] → {} (tag {})", world_rank, down_node.world_rank, tag));
            MPI_Isend(
                  zc_send_buf[ii].data(),
                  ZCBuffer::size_down_up,
                  MPI_DOUBLE,
                  down_node.world_rank,
                  tag,
                  world_comm,
                  &reqs[reqs_used]
            );
            ++reqs_used;
         }
      }
      // Wait until all receive and send tasks are retired from the network.
      MPI_Waitall(reqs_used, reqs.data(), MPI_STATUSES_IGNORE);
      MPI_Win_fence(0, sector_win);

   }

   // The halo exchange operation in which data is sent to the bottom neighbor.
   void halo_exchange_send_bottom() {
      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         if (connectivity.inside_domain(node.B)) {
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::NGhost; ++k) {
               for (size_t j = 0; j < P::Ntp_zones; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = patches[ii].U[{cci, k}][{j, i}].under();
                        zc_put_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            auto bottom_node = connectivity.get(node.B);
            step_communication(std::format("[{}] → [win:{}:{}]", world_rank, bottom_node.sector_rank, bottom_node.proc_ii));
            // Put zc_put_buf into zc_recv_buf.
            MPI_Put(
                  zc_put_buf[ii].data(),
                  ZCBuffer::size_bottom_top,
                  MPI_DOUBLE,
                  bottom_node.sector_rank,
                  bottom_node.proc_ii*sizeof(ZCBuffer),
                  ZCBuffer::size_bottom_top,
                  MPI_DOUBLE,
                  sector_win
            );
         }
      }
      MPI_Win_fence(0, sector_win);
   }

   // The halo exchange operation in which data is sent to the top neighbor.
   void halo_exchange_send_top() {
      for (int32 ii = 0; ii < n_patch; ++ii) {
         auto node = connectivity.get(patches[ii].index);
         if (connectivity.inside_domain(node.T)) {
            size_t idx = 0;
            for (size_t k = P::nr_zones; k < P::nr_zones + P::NGhost; ++k) {
               for (size_t j = 0; j < P::Ntp_zones; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = patches[ii].U[{cci, k}][{j, i}].under();
                        zc_put_buf[ii][idx++] = x;
                     }
                  }
               }
            }
            auto top_node = connectivity.get(node.T);
            step_communication(std::format("[{}] → [win:{}:{}]", world_rank, top_node.sector_rank, top_node.proc_ii));
            // Put zc_put_buf into zc_recv_buf.
            MPI_Put(
                  zc_put_buf[ii].data(),
                  ZCBuffer::size_bottom_top,
                  MPI_DOUBLE,
                  top_node.sector_rank,
                  top_node.proc_ii*sizeof(ZCBuffer),
                  ZCBuffer::size_bottom_top,
                  MPI_DOUBLE,
                  sector_win
            );
         }
      }
      MPI_Win_fence(0, sector_win);
   }


   // load the working patch zone structure with
   // the data received from neighbor patches.
   void halo_exchange_writeback_left() {
      for (int ii = 0; ii < n_patch; ++ii) {
         size_t idx = 0;
         for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost + P::ntp_zones; ++j) {
               for (size_t i = P::ntp_zones+P::NGhost; i < P::Ntp_zones; ++i) {
                  for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                     auto x = zc_recv_buf[ii].zc_buffer[idx++];
                     patches[ii].U[{cci, k}][{j, i}] = R64(x);
                  }
               }
            }
         }
      }
   }

   // load the working patch zone structure with
   // the data received from neighbor patches.
   void halo_exchange_writeback_right() {
      for (int ii = 0; ii < n_patch; ++ii) {
         size_t idx = 0;
         for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
            for (size_t j = P::NGhost; j < P::NGhost + P::ntp_zones; ++j) {
               for (size_t i = 0; i < P::NGhost; ++i) {
                  for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                     auto x = zc_recv_buf[ii].zc_buffer[idx++];
                     patches[ii].U[{cci, k}][{j, i}] = R64(x);
                  }
               }
            }
         }
      }
   }

   // load the working patch zone structure with
   // the data received from neighbor patches.
   void halo_exchange_writeback_down() {
      // unload what is received via Put.
      for (int ii = 0; ii < n_patch; ++ii) {
         auto& patch = patches[ii];
         auto node = connectivity.get(patch.index);
         auto up_node = connectivity.get(node.U);
         if (up_node.index.isector != node.index.isector) {
            // unload what is received via Isend/Irecv.
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::ntp_zones+P::NGhost; j < P::Ntp_zones; ++j) {
                  for (size_t i = P::NGhost; i < P::NGhost + P::ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = zc_recv_buf[ii].zc_buffer[idx++];
                        patches[ii].U[{cci, k}][{j, i}] = R64(x);
                     }
                  }
               }
            }
            // up diagonal
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::ntp_zones+P::NGhost; j < P::Ntp_zones; ++j) {
                  for (size_t i = 0; i < P::NGhost; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = zc_recv_buf[ii].zc_buffer[idx++];
                        patches[ii].U[{cci, k}][{j, i}] = R64(x);
                     }
                  }
               }
            }
            // down diagonal
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::ntp_zones+P::NGhost; j < P::Ntp_zones; ++j) {
                  for (size_t i = P::NGhost + P::ntp_zones; i < P::NGhost + P::NGhost + P::ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = zc_recv_buf[ii].zc_buffer[idx++];
                        patches[ii].U[{cci, k}][{j, i}] = R64(x);
                     }
                  }
               }
            }
         }
         else {
            // unload what is received via Put.
            size_t idx = 0;
            for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
               for (size_t j = P::ntp_zones + P::NGhost; j < P::Ntp_zones; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = zc_recv_buf[ii].zc_buffer[idx++];
                        patches[ii].U[{cci, k}][{j, i}] = R64(x);
                     }
                  }
               }
            }
         }
      }
   }

   // load the working patch zone structure with
   // the data received from neighbor patches.
   void halo_exchange_writeback_up() {
      for (int ii = 0; ii < n_patch; ++ii) {
         size_t idx = 0;
         for (size_t k = P::NGhost; k < P::NGhost + P::nr_zones; ++k) {
            for (size_t j = 0; j < P::NGhost; ++j) {
               for (size_t i = 0; i < P::Ntp_zones; ++i) {
                  for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                     auto x = zc_recv_buf[ii].zc_buffer[idx++];
                     patches[ii].U[{cci, k}][{j, i}] = R64(x);
                  }
               }
            }
         }
      }
   }

   // load the working patch zone structure with
   // the data received from neighbor patches.
   void halo_exchange_writeback_bottom() {
      for (int ii = 0; ii < n_patch; ++ii) {
         auto& patch = patches[ii];
         auto node = connectivity.get(patch.index);
         if (connectivity.inside_domain(node.B)) {
            size_t idx = 0;
            for (size_t k = P::nr_zones+P::NGhost; k < P::Nr_zones; ++k) {
               for (size_t j = 0; j < P::Ntp_zones; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = zc_recv_buf[ii].zc_buffer[idx++];
                        patches[ii].U[{cci, k}][{j, i}] = R64(x);
                     }
                  }
               }
            }
         }
      }
   }

   // load the working patch zone structure with
   // the data received from neighbor patches.
   void halo_exchange_writeback_top() {
      for (int ii = 0; ii < n_patch; ++ii) {
         auto& patch = patches[ii];
         auto node = connectivity.get(patch.index);
         if (connectivity.inside_domain(node.T)) {
            size_t idx = 0;
            for (size_t k = 0; k < P::NGhost; ++k) {
               for (size_t j = 0; j < P::Ntp_zones; ++j) {
                  for (size_t i = 0; i < P::Ntp_zones; ++i) {
                     for (size_t cci = 0; cci < P::n_cc_components; ++cci) {
                        auto x = zc_recv_buf[ii].zc_buffer[idx++];
                        patches[ii].U[{cci, k}][{j, i}] = R64(x);
                     }
                  }
               }
            }
         }
      }
   }

public: // messages:

   void print_panic(const std::string& msg) const {
      if (world_rank == 0) {
         std::cout << "[PANIC] " << msg << std::endl;
      }
   }

   void print_warning(const std::string& msg) const {
      if (world_rank == 0) {
         std::cout << "[Warning] " << msg << std::endl;
      }
   }

   void print_message(const std::string& msg) const {
      if (world_rank == 0) {
         std::cout << "[Status] " << msg << std::endl;
      }
   }

   void step_communication(const std::string& msg, int32 select_rank = -600) const {
         if constexpr (P::Directives::debug_communication) {
            if (select_rank == -600)
               std::cout << "[" << world_rank << "] " << msg << std::endl;
            else {
               if (world_rank == select_rank)
                  std::cout << "[" << world_rank << "] " << msg << std::endl;
            }
         }
   }

public: // debug:

   [[nodiscard]]
   std::string str() const {
      std::string out = std::format("Process rank {}/{}\n", world_rank, world_size);
      out += std::format("sector color: {}\nSector rank {}/{}\n", sector_color, sector_rank, sector_size);
      out += std::format("n_patch {}\n", n_patch);
      return out;
   }

};


