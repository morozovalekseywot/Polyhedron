# Tries to find HDF5.
#
# Usage of this module as follows:
#
#   find_package(HDF5)
#
# Variables used by this module,
# they can change the default behaviour and need to be set before calling find_package:
#
#   HDF5_ROOT_DIRS          Root instalation of HDF5 (list)
#
# Variables defined by this module:
#
#   HDF5_FOUND              System has HDF5 libs/headers
#   HDF5_LIBRARIES          The HDF5 library
#   HDF5_INCLUDE_DIR        The location of HDF5 headers

find_library(HDF5_LIBRARIES
        NAMES hdf5
        HINTS ${HDF5_ROOT_DIRS}
        PATH_SUFFIXES lib lib64
        NO_DEFAULT_PATH)

find_path(HDF5_INCLUDE_DIR
        NAMES hdf5.h
        PATHS ${HDF5_ROOT_DIRS}
        PATH_SUFFIXES include
        NO_DEFAULT_PATH)

include(FindPackageHandleStandardArgs)

set(SUCCESS "TRUE")

find_package_handle_standard_args(
        HDF5
        FOUND_VAR HDF5_FOUND
        REQUIRED_VARS SUCCESS HDF5_LIBRARIES HDF5_INCLUDE_DIR
        FAIL_MESSAGE "HDF5 package is not found. Set right HDF5_ROOT_DIRS in CMakeLists.txt.\n  ")

mark_as_advanced(
        HDF5_ROOT_DIRS
        HDF5_LIBRARIES
        HDF5_INCLUDE_DIR)