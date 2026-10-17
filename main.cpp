

#include <T.h>
#include <I.h>
#include <N.h>

int main(int argc, char** argv) {

   MPI_Init(&argc, &argv);

   // Unit Tests:

//   TMesh::Tdryrun();

//    TPatch::Tdryrun();
//    TPatch::Txf_quadpts_patch();
//    Tquadpoints::Tdryrun();

//   TConnectivity::Tdryrun();
//   TConnectivity::Tvisualize_xyz();
//   TConnectivity::Tvisualize_Fortran();

//   TProcess::Tdryrun();

   // warm-up cases for halo exchange testing:
//   TProcess::Tconnectivityview();
//   TProcess::Tglobal_view();
//   TProcess::Texchange_1patch();
//   TProcess::Texchange_9patch();
//   TProcess::Texchange_2ghost();
//   TProcess::Texchange_plot_clinic();

   // plotting routines for halo exchange testing:
//   TProcess_Plot::Texchange_plot();

   // Integration Tests:

//   IAugust30::start();


   // Tests of automated tests

//   NPatch::Npatch_test1();
//   NPatch::Npatch_test2();
   NPatch::Npatch_test3();

   int finalized = 0;
   // Check if MPI was already shut down elsewhere.
   MPI_Finalized(&finalized);
   if (!finalized) {
      MPI_Finalize();
   }

}

