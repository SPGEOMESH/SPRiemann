
#pragma once

#include <Nbasis.h>

#include <csriemann/constants/parameters.h>
#include <csriemann/constants/directives.h>
#include <csriemann/process.h>

#include <patch_global_view.h>

#include <cmath>

/*
 * Automated testing of halo exchange in a
 * setting where data is partitioned into patches and sectors.
 * Basic test elements are explained/illustrated in T-tests.
 * This is an automated test for high-impact agile development.
 */
class NPatch {
public:

   // A Parameter set whose only change is ntp_patches.
   // The number of zones per sector remains the same as `iter` changes.
   // Every iteration of the test, the number of patches is divided by `test_multiple`.
   // That is, in each iteration, the patches are smaller and there are more.
   // See the illustration `Npatch_diagram.png`.
   // It is best if the test_multiple is at least 4, i.e. test_factor ≥ 2.
   // Use test_factor 1 only for testing of the test code.
   // It is important to run the test with ilevel_max nontrivial (> 1),
   // to test the halo exchange on a full adaptive mesh,
   // independent of (and prior to) tests of other aspects of AMR capabilities.
   // todo ilevel_max > 1
   template <typename Directives, int32 iter_, int32 ntp_zones_per_sector_, int32 test_factor_, int32 NGhost_, int32 ilevel_max_>
   struct Parameters_Iter: public Parameters<
         Directives,
         // in each iteration, the number of zones per patch is divided by 2^test_factor.
         ntp_zones_per_sector_/DGCALC::ipow2<int32>(iter_*test_factor_),
         ntp_zones_per_sector_/DGCALC::ipow2<int32>(iter_*test_factor_),
         // in each iteration, the number of patches per sector is multiplied by 2^test_factor.
         DGCALC::ipow2<int32>(iter_*test_factor_),
         DGCALC::ipow2<int32>(iter_*test_factor_),
         NGhost_,
         ilevel_max_
   > {
            static constexpr int32 iter = iter_;
            static constexpr int32 test_multiple = DGCALC::ipow2<int32>(test_factor_);
   };


   // The dummy data must depend only on the sector-wide coordinates of zones.
   // That is, it must be independent of ntp_patches and nr_patches.
   template <typename P, typename itype>
   static R64 dummy_data(itype r, itype phi, itype theta, int32 iz, int32 iy, int32 ix, int32 is, int32 ilevel) {
      auto sx = R64(P::SECt(ix, static_cast<int32>(theta)));
      auto sy = R64(P::SECt(iy, static_cast<int32>(phi)));
      auto sz = R64(P::SECt(iz, static_cast<int32>(r)));
//      auto u = sy.cos()*(R64(-1)*sx*sx).exp() + sy.sin()*(R64(-1)*sz*sz).exp() + R64(is)*R64(ilevel);
      auto u = R64(777) + R64(10000)*sz + R64(1000000)*sy + R64(100000000)*sx + R64(10000000000)*R64(is);
      return u;
   }

   static std::string handle(int32 iter, int32 isector, int32 ilevel, int32 ix, int32 iy, int32 iz) {
      std::string out = std::format("iter{}is{}il{}.x{}y{}z{}", iter, isector, ilevel, ix, iy, iz);
      return out;
   }


   template <typename Parameters_Iter, bool detail = false, bool final = false>
   static void test_iteration() {

      using P = Parameters_Iter;

      constexpr int32 iter = P::iter;

      using Process = Process<P>;

      auto process = Process();

      // After initialization, a process has patch memory.

      // > Initialize dummy simulation data.
      for (int32 ii = 0; ii < process.n_patch; ++ii) {
         auto& patch = process.patches[ii];
         auto index = patch.index;
         auto is = index.isector+1;
         auto ilevel = index.ilevel;
         auto ix = index.ix;
         auto iy = index.iy;
         auto iz = index.iz;
         for (size_t r = P::NGhost; r < P::NGhost+P::nr_zones; ++r) {
            for (size_t phi = P::NGhost; phi < P::NGhost+P::ntp_zones; ++phi) {
               for (size_t theta = P::NGhost; theta < P::NGhost+P::ntp_zones; ++theta) {
                  auto u = dummy_data<P>(r, phi, theta, iz, iy, ix, is, ilevel);
                  process.patches[ii].U[{0,r}][{phi, theta}] = u;
               }
            }
         }
      }

      // > perform halo exchange. This is the section of code being tested.
      process.halo_exchange();

//      GlobalView<P>::make(process.patches.get(), process.n_patch, "TEST");

      if constexpr (iter != 0) {
         // not base case. We know that the previous case already spoofed this case.
         // Perform a straightforward reduction in each patch.
         for (int32 ii = 0; ii < process.n_patch; ++ii) {
            auto& patch = process.patches[ii];
            auto index = patch.index;
            auto isector = index.isector+1;
            auto ilevel = index.ilevel;
            auto ix = index.ix;
            auto iy = index.iy;
            auto iz = index.iz;
            R64 umax(-1e30);
            R64 umin(1e30);
            R64 uavg(0.0);
            R64 ustd(0.0);
            R64 delta1, delta2;
            std::string content;
            R64 size3d = R64(P::Ntp_zones)*R64(P::Ntp_zones)*R64(P::Nr_zones);
            // We deliberately calculate over all "padded" patches
            // that include all the ghost zones around the patch.
            for (size_t r = 0; r < P::Nr_zones; ++r) {
               for (size_t phi = 0; phi < P::Ntp_zones; ++phi) {
                  for (size_t theta = 0; theta < P::Ntp_zones; ++theta) {
                     // Do not include in the reduction if it is a tricorner zone.
                     if (!P::is_tricorner_zone(ix, iy, theta, phi)) {
                        auto u = process.patches[ii].U[{0,r}][{phi, theta}];
                        if (u > umax) umax = u;
                        if (u < umin) umin = u;
                        delta1 = u - uavg;
                        uavg += delta1/size3d;
                        delta2 = u - uavg;
                        ustd += delta1*delta2;
                     }
                  }
               }
            }
            // Treat the tricorner zones as negligible error in the avg/stddev.
            ustd = (size3d < R64(2))? R64(0): (ustd/(size3d-R64(1))).sqrt();
            // Have umax, umin, uavg, ustd.
            if (!detail) {
               std::string filename = handle(iter, isector, ilevel, ix, iy, iz) + ".real.dat";
               content += std::format("{:.15e} {:.15e} {:.15e} {:.15e}\n", umin.under(), umax.under(), uavg.under(), ustd.under());
               savetxt(filename, content);
            }
            else {
               // for a `detail` check, only check one sector.
               if (isector == 1) {
                  std::string filename = handle(iter, isector, ilevel, ix, iy, iz) + ".real.dat";
                  auto view = GlobalView<P>::make_3d(patch);
                  content += view;
                  content += std::format("This is a real patch.\nReal: ntp_zones {} ntp_patches {}\n", P::ntp_zones, P::ntp_patches);
                  content += std::format("{:.15e} {:.15e} {:.15e} {:.15e}\n", umin.under(), umax.under(), uavg.under(), ustd.under());
                  savetxt(filename, content);
               }
            }
         }
      }

      if (!final) {
         // spoof the next iteration.
         for (int32 ii = 0; ii < process.n_patch; ++ii) {
            auto& patch = process.patches[ii];
            auto index = patch.index;
            auto isector = index.isector+1;
            auto ilevel = index.ilevel;
            auto ix = index.ix;
            auto iy = index.iy;
            auto iz = index.iz;
            auto TM = P::test_multiple;
            auto ntp_zones_spoof = P::ntp_zones/TM;
            auto nr_zones_spoof = P::nr_zones/TM;
            auto ntp_patches_spoof = P::ntp_patches*TM;
            // There are TF^3 spoofed patches "inside" this patch.
            for (int32 iz_ = 0; iz_ < TM; ++iz_) {
               for (int32 iy_ = 0; iy_ < TM; ++iy_) {
                  for (int32 ix_ = 0; ix_ < TM; ++ix_) {
                     int32 ix_spoof = TM*ix + ix_;
                     int32 iy_spoof = TM*iy + iy_;
                     int32 iz_spoof = TM*iz + iz_;
                     R64 umax(-1e30);
                     R64 umin(1e30);
                     R64 uavg(0.0);
                     R64 ustd(0.0);
                     R64 delta1, delta2;
                     std::string content;
                     std::string detail_content;
                     R64 size3d = R64((ntp_zones_spoof+2*P::NGhost)*(ntp_zones_spoof+2*P::NGhost)*(nr_zones_spoof+2*P::NGhost));
                     for (size_t r = 0; r < 2*P::NGhost+nr_zones_spoof; ++r) {
                        for (size_t phi = 0; phi < 2*P::NGhost+ntp_zones_spoof; ++phi) {
                           for (size_t theta = 0; theta < 2*P::NGhost+ntp_zones_spoof; ++theta) {
                              // Do not include in the reduction if it is a tricorner zone.
                              if (!P::is_tricorner_zone_spoof(ix_spoof, iy_spoof, theta, phi, ntp_zones_spoof, ntp_patches_spoof)) {
                                 auto u = process.patches[ii].U[{0,nr_zones_spoof*iz_ + r}][{ntp_zones_spoof*iy_ + phi, ntp_zones_spoof*ix_ + theta}];
                                 if (u > umax) umax = u;
                                 if (u < umin) umin = u;
                                 if (detail) {
                                    detail_content += std::format("theta {} phi {} r {} u {}\n", theta, phi, r, u.under());
                                 }
                                 delta1 = u - uavg;
                                 uavg += delta1/size3d;
                                 delta2 = u - uavg;
                                 ustd += delta1*delta2;
                              }
                           }
                        }
                     }
                     // Treat the tricorner zones as negligible error in the avg/stddev.
                     ustd = (size3d < R64(2))? R64(0): (ustd/(size3d-R64(1))).sqrt();
                     // Have umax, umin, uavg, ustd.
                     if (!detail) {
                        std::string filename = handle(iter+1, isector, ilevel, ix_spoof, iy_spoof, iz_spoof) + ".spoof.dat";
                        content += std::format("{:.15e} {:.15e} {:.15e} {:.15e}\n", umin.under(), umax.under(), uavg.under(), ustd.under());
                        savetxt(filename, content);
                     }
                     else {
                        // for a `detail` check, only check one sector.
                        if (isector == 1) {
                           std::string filename = std::format("iter{}is{}il{}.x{}y{}z{}.spoof.dat", iter+1, isector, ilevel, ix_spoof, iy_spoof, iz_spoof);
                           auto view = GlobalView<P>::make_3d(patch);
                           content += view;
                           content += std::format("This is a real patch.\nReal: ntp_zones {} ntp_patches {}\n", P::ntp_zones, P::ntp_patches);
                           // We use this as a satisfactory view of the spoof patch.
                           content += detail_content;
                           content += std::format("This is a spoof patch.\nSpoof: ntp_zones {} ntp_patches {}\n", ntp_zones_spoof, ntp_patches_spoof);
                           content += std::format("{:.15e} {:.15e} {:.15e} {:.15e}\n", umin.under(), umax.under(), uavg.under(), ustd.under());
                           savetxt(filename, content);
                        }
                     }
                  }
               }
            }
         }
      }
      process.print_message(std::format("Done computing iter {}", iter));

   }

   template <typename Parameters, int32 detail>
   static bool test_post_process() {
      // pass in the last iter type used (somewhat awk but ok)
      using P = Parameters;
      bool test_evaluation = false;
      if (!detail) {
         // proceed with automated processing, by comparing the saved data from separate runs.
         int32 iter_max = P::iter;
         int32 isector_max = 6;
         int32 ilevel_max = P::ilevel_max;
         int32 count_bad = 0;
         int32 count_total = 0;
         int32 count_error = 0;
         int32 ix_max = P::ntp_patches;
         int32 iy_max = P::ntp_patches;
         int32 iz_max = P::nr_patches;
         auto TM = P::test_multiple;
         for (int32 iter = iter_max; iter >= 1; --iter) {
            for (int32 isector = 1; isector <= isector_max; ++isector) {
               for (int32 ilevel = 0; ilevel < ilevel_max; ++ilevel) {
                  for (int32 iz = 0; iz < iz_max; ++iz) {
                     for (int32 iy = 0; iy < iy_max; ++iy) {
                        for (int32 ix = 0; ix < ix_max; ++ix) {
                           ++count_total;
                           auto H = handle(iter, isector, ilevel, ix, iy, iz);
                           auto Hreal = H+".real.dat";
                           auto Hspoof = H+".spoof.dat";
                           double minreal, maxreal, avgreal, stdreal;
                           double minspf, maxspf, avgspf, stdspf;
                           std::fstream s;
                           s.exceptions(std::ifstream::failbit | std::ifstream::badbit);
                           try {
                              s.open(Hreal.c_str(), std::fstream::in);
                              s >> minreal; s >> maxreal; s >> avgreal; s >> stdreal;
                              s.close();
                           } catch (const std::exception& e) {
                              std::cerr << "Error: " << e.what() << '\n';
                              std::cerr << "Error while loading file. Case " << H << " could not be completed." << endl;
                              ++count_error;
                              continue;
                           }
                           try {
                              s.open(Hspoof.c_str(), std::fstream::in);
                              s >> minspf; s >> maxspf; s >> avgspf; s >> stdspf;
                              s.close();
                           } catch (const std::exception& e) {
                              std::cerr << "Error: " << e.what() << '\n';
                              std::cerr << "Error while loading file. Case " << H << " could not be completed." << endl;
                              ++count_error;
                              continue;
                           }
                           if (!(isclose(minreal, minspf) && isclose(maxreal, maxspf) && isclose(avgreal, avgspf) && isclose(stdreal, stdspf))) {
                              cout << "Problem Found At Case: " << H << endl;
                              if (!isclose(minreal, minspf)) cout << std::format("min: real {} spoof {}", minreal, minspf) << endl;
                              if (!isclose(maxreal, maxspf)) cout << std::format("max: real {} spoof {}", maxreal, maxspf) << endl;
                              if (!isclose(avgreal, avgspf)) cout << std::format("avg: real {} spoof {}", avgreal, avgspf) << endl;
                              if (!isclose(stdreal, stdspf)) cout << std::format("stddev: real {} spoof {}", stdreal, stdspf) << endl;
                              ++count_bad;
                           }
                           else {
                              bool e;
                              e = delete_file(Hreal);
                              if (!e) cerr << "Could not delete " << Hreal << endl;
                              delete_file(Hspoof);
                              if (!e) cerr << "Could not delete " << Hspoof << endl;
                           }
                        }
                     }
                  }
               }
            }
            ix_max /= TM;
            iy_max /= TM;
            iz_max /= TM;
            cout << "[Post Processing] Done with iter " << iter << endl;
         }
         auto bad_pct = static_cast<double>(count_bad)/static_cast<double>(count_total)*100.0;
         auto error_pct = static_cast<double>(count_error)/static_cast<double>(count_total)*100.0;
         cout << std::format("Number of bad cases: {}/{} ({:.5g}%)", count_bad, count_total, bad_pct) << endl;
         if (count_error > 0) {
            cout << std::format("Number of error cases (skipped due to errors): {}/{} ({:.5g}%)", count_error, count_total, error_pct) << endl;
         }
         if (count_bad == 0 && count_error == 0) {
            cout << "Test: PASS, ZERO CASES SKIPPED (PERFECT PASS)." << endl;
            test_evaluation = true;
         }
         else if (count_bad == 0 && error_pct < 5.0) {
            cout << "Test: PASS." << endl;
            test_evaluation = true;
         }
         else if (count_bad == 0) {
            cout << "Test: THIS TEST HAD ISSUES." << endl;
            test_evaluation = false;
         }
         else {
            cout << "Test: FAIL." << endl;
            test_evaluation = false;
         }
         cout << "Done." << endl;
      }
      return test_evaluation;
   }


   // This is a downscaled test, for manual review purposes.
   static void Npatch_test1() {

      using Directives = Directives<>;

      constexpr size_t NGhost = 1;
      constexpr size_t ilevel_max = 1;
      // test_multiple = 2^test_factor
      constexpr int32 test_factor = 1;
      constexpr int32 ntp_zones_per_sector = 2;

      using Param0 = Parameters_Iter<Directives, 0, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;
      using Param1 = Parameters_Iter<Directives, 1, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;

      // set `detail` to produce zone-by-zone analysis for diagnostic/troubleshooting.
      constexpr bool detail = true;
      // set `final` to stop from spoofing the next iteration.
      constexpr bool final = true;

      test_iteration<Param0, detail>();
      test_iteration<Param1, detail, final>();

   }


   // This is a downscaled test, for manual review purposes.
   static bool Npatch_test2() {

      using Directives = Directives<>;

      constexpr size_t NGhost = 2;
      constexpr size_t ilevel_max = 1;
      // test_multiple = 2^test_factor
      constexpr int32 test_factor = 2;
      constexpr int32 ntp_zones_per_sector = 8;

      using Param0 = Parameters_Iter<Directives, 0, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;
      using Param1 = Parameters_Iter<Directives, 1, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;

      // set `detail` to produce zone-by-zone analysis for diagnostic/troubleshooting.
      constexpr bool detail = true;
      // set `final` to stop from spoofing the next iteration.
      constexpr bool final = true;

      test_iteration<Param0, detail>();
      test_iteration<Param1, detail, final>();

      bool test_evaluation = false;
      if (test_serial_processing()) {
         test_evaluation = test_post_process<Param1, detail>();
      }
      test_end_serial_processing();
      return test_evaluation;
   }


   // This is a *fully automated* test.
   static bool Npatch_test3() {

      using Directives = Directives<>;

      constexpr size_t NGhost = 2;
      constexpr size_t ilevel_max = 1;
      // test_multiple = 2^test_factor
      constexpr int32 test_factor = 2;
      constexpr int32 ntp_zones_per_sector = 32;

      using Param0 = Parameters_Iter<Directives, 0, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;
      using Param1 = Parameters_Iter<Directives, 1, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;
      using Param2 = Parameters_Iter<Directives, 2, ntp_zones_per_sector, test_factor, NGhost, ilevel_max>;

      // set `detail` to produce zone-by-zone analysis for diagnostic/troubleshooting.
      constexpr bool detail = false;
      // set `final` to stop from spoofing the next iteration.
      constexpr bool final = true;

      test_iteration<Param0, detail>();
      test_iteration<Param1, detail>();
      test_iteration<Param2, detail, final>();

      bool test_evaluation = false;
      if (test_serial_processing()) {
         test_evaluation = test_post_process<Param2, detail>();
      }
      test_end_serial_processing();
      return test_evaluation;
   }


};



