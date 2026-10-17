
// Edit this list by adding arguments
// like argv for command line control of Catch2.
/*
 * Example:
 * $ foo -?
 * Corresponds to
 * const char * catch_argv[]{"path_of_executable", "-?", nullptr};
 *
 */
const char * catch_argv[]{
        "path_of_executable",
        "-?",
        nullptr
};

//#include <catch2/catch.hpp>


#include <Tmesh/Tmesh.h>
#include <Tmesh/Tquadpoints.h>

#include <Tconnectivity/Tconnectivity.h>
#include <Tpatch/Tpatch.h>


#include <Tprocess/Tprocess.h>
#include <Tprocess/Tprocess_plot.h>




