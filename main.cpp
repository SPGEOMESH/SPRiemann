

#include <T.h>
#include <I.h>

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
   TProcess::Texchange_plot();

   // Integration Tests:

//   IAugust30::start();

   int finalized = 0;
   // Check if MPI was already shut down elsewhere.
   MPI_Finalized(&finalized);
   if (!finalized) {
      MPI_Finalize();
   }

}

