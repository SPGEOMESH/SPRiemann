

#include <T.h>
#include <I.h>

int main(int argc, char** argv) {

    // Unit Tests:

//   TMesh::Tdryrun();

//    TPatch::Tdryrun();
//    TPatch::Txf_quadpts_patch();
//    Tquadpoints::Tdryrun();


//   TConnectivity::Tdryrun();
//   TConnectivity::Tvisualize_xyz();
//   TConnectivity::Tvisualize_Fortran();


//   TProcess::Tdryrun(argc, argv);

   // warm-up cases for halo exchange testing:
//   TProcess::Tconnectivityview(argc, argv);
//   TProcess::Tglobal_view(argc, argv);
//   TProcess::Texchange_1patch(argc, argv);
   TProcess::Texchange_9patch(argc, argv);
//   TProcess::Texchange_2ghost(argc, argv);


   // Integration Tests:

//   IAugust30::start(argc, argv);

}

