cmake_minimum_required(VERSION 3.24)
set(PROJECT_LANGUAGES "CXX" "Fortran")
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_FORTRAN_STANDARD 90)
#set(CMAKE_FORTRAN_STANDARD_REQUIRED ON)
# fixed form, C preprocessed, all REAL as double
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffixed-form -cpp -fdefault-real-8")
set(PROJECT_FAMILY_NAME "SPRiemann")
set(PROJECT_VERSION "0.5.1")
set(PROJECT_VERSION_NAME "NashvilleDolly")
set(PROJECT_VERSION_MONTH "September")
set(PROJECT_VERSION_YEAR "2026")
set(PROJECT_WORKING_RDIRS "_T" "_I" "_N")
set(PROJECT_BUILD_TYPE "Debug") # Release, Debug
set(ACCEL "") # "Cuda", "Metal", "" (none)
set(ACCEL_DEPENDENCIES "") # cuBLAS, cuRAND, ...
set(ACCEL_INSTALL_DIR "/usr/local") # /usr/local, $ENV{CUDA_HOME}, $ENV{METAL_HOME}, ...
set(DEPENDENCY_LIB_TYPE "s") # shared (s) or static (a)
set(GROUPS "DGCALC" "SPGEOMESH")
set(DGCALC_DEPENDENCY_TEST "")
set(DGCALC_DEPENDENCY_SRC ${DGCALC_DEPENDENCY_TEST} "Base" "Types" "Scalar" "Kyte" "Tennsix")
set(DGCALC_DEPENDENCY_LIB "")
set(SPGEOMESH_DEPENDENCY_TEST "")
set(SPGEOMESH_DEPENDENCY_SRC ${SPGEOMESH_DEPENDENCY_TEST} "Kiwi" "Fields")
set(SPGEOMESH_DEPENDENCY_LIB "CubeSphere")

#set(NUM_PROCESSES 6)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address -fno-omit-frame-pointer")

# set(LOCAL_ROOT "")
if (LOCAL_ROOT)
    include(${LOCAL_ROOT}/src/dgcalc.1.cmake)
else()
    execute_process(COMMAND dgcalc source 1 OUTPUT_VARIABLE SOURCE1_CODE)
    cmake_language(EVAL CODE "${SOURCE1_CODE}")
endif()
set(PROJECT_NAME "${PROJECT_FAMILY_NAME}${PROJECT_VERSION}${PROJECT_VERSION_NAME}_${PROJECT_VERSION_MONTH}${PROJECT_VERSION_YEAR}")
project(${PROJECT_NAME} VERSION ${PROJECT_VERSION} LANGUAGES ${PROJECT_LANGUAGES})
if (LOCAL_ROOT)
    include(${LOCAL_ROOT}/src/dgcalc.2.cmake)
else()
    execute_process(COMMAND dgcalc source 2 OUTPUT_VARIABLE SOURCE2_CODE)
    cmake_language(EVAL CODE "${SOURCE2_CODE}")
endif()


