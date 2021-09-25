# Tries to find Gperftools.
#
# Usage of this module as follows:
#
#   find_package(Gperftools)
#
# Variables used by this module,
# they can change the default behaviour and need to be set before calling find_package:
#
#   Gperftools_ROOT_DIRS          Root instalation of HDF5 (list)
#
# Variables defined by this module:
#
#   GPERFTOOLS_FOUND              System has Gperftools libs/headers
#   GPERFTOOLS_LIBRARIES          The Gperftools libraries (tcmalloc & profiler)

find_library(GPERFTOOLS_TCMALLOC
        NAMES tcmalloc libtcmalloc.so.4
        HINTS ${Gperftools_ROOT_DIRS}
        PATH_SUFFIXES lib lib64)

find_library(GPERFTOOLS_PROFILER
        NAMES profiler libprofiler.so.4
        HINTS ${Gperftools_ROOT_DIRS}
        PATH_SUFFIXES lib lib64)

find_library(GPERFTOOLS_TCMALLOC_AND_PROFILER
        NAMES tcmalloc_and_profiler libtcmalloc_and_profiler.so.4
        HINTS ${Gperftools_ROOT_DIRS}
        PATH_SUFFIXES lib lib64)

set(GPERFTOOLS_LIBRARIES ${GPERFTOOLS_TCMALLOC})

include(FindPackageHandleStandardArgs)
set(SUCCESS "TRUE")

find_package_handle_standard_args(
        Gperftools
        FOUND_VAR Gperftools_FOUND
        REQUIRED_VARS SUCCESS GPERFTOOLS_LIBRARIES
        FAIL_MESSAGE "Gperftools (GPT) is not found, add GPT root directory to list in FindGperftools.cmake.")

mark_as_advanced(
        Gperftools_ROOT_DIRS
        GPERFTOOLS_TCMALLOC
        GPERFTOOLS_PROFILER
        GPERFTOOLS_TCMALLOC_AND_PROFILER
        GPERFTOOLS_LIBRARIES)
