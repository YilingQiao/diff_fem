# - Try to find Eigen3 lib
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(Eigen3 3.1.2)
# to require version 3.1.2 or newer of Eigen3.
#
# Once done this will define
#
#  EIGEN_FOUND - system has eigen lib with correct version
#  EIGEN_INCLUDE_DIR - the eigen include directory
#  EIGEN_VERSION - eigen version

# Copyright (c) 2006, 2007 Montel Laurent, <montel@kde.org>
# Copyright (c) 2008, 2009 Gael Guennebaud, <g.gael@free.fr>
# Copyright (c) 2009 Benoit Jacob <jacob.benoit.1@gmail.com>
# Redistribution and use is allowed according to the terms of the 2-clause BSD license.

if(NOT Eigen_FIND_VERSION)
  if(NOT Eigen_FIND_VERSION_MAJOR)
    set(Eigen_FIND_VERSION_MAJOR 2)
  endif(NOT Eigen_FIND_VERSION_MAJOR)
  if(NOT Eigen_FIND_VERSION_MINOR)
    set(Eigen_FIND_VERSION_MINOR 91)
  endif(NOT Eigen_FIND_VERSION_MINOR)
  if(NOT Eigen_FIND_VERSION_PATCH)
    set(Eigen_FIND_VERSION_PATCH 0)
  endif(NOT Eigen_FIND_VERSION_PATCH)

  set(Eigen_FIND_VERSION "${Eigen_FIND_VERSION_MAJOR}.${Eigen_FIND_VERSION_MINOR}.${Eigen_FIND_VERSION_PATCH}")
endif(NOT Eigen_FIND_VERSION)


# Construct consistent error messages for use below.
set(EIGEN_DIR_DESCRIPTION "directory containing the file 'Eigen/Core', e.g. <EIGEN_INCLUDE_DIR>=PREFIX/include/eigen3.")
set(EIGEN_DIR_MESSAGE "EIGEN not found.  Set the EIGEN_INCLUDE_DIR cmake cache entry to the ${EIGEN_DIR_DESCRIPTION}")

macro(_eigen3_check_version)

  file(READ "${EIGEN_INCLUDE_DIR}/Eigen/src/Core/util/Macros.h" _eigen3_version_header)

  string(REGEX MATCH "define[ \t]+EIGEN_WORLD_VERSION[ \t]+([0-9]+)" _eigen3_world_version_match "${_eigen3_version_header}")
  set(EIGEN_WORLD_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MAJOR_VERSION[ \t]+([0-9]+)" _eigen3_major_version_match "${_eigen3_version_header}")
  set(EIGEN_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+EIGEN_MINOR_VERSION[ \t]+([0-9]+)" _eigen3_minor_version_match "${_eigen3_version_header}")
  set(EIGEN_MINOR_VERSION "${CMAKE_MATCH_1}")

  set(EIGEN_VERSION ${EIGEN_WORLD_VERSION}.${EIGEN_MAJOR_VERSION}.${EIGEN_MINOR_VERSION})
  if(${EIGEN_VERSION} VERSION_LESS ${Eigen_FIND_VERSION})
    set(EIGEN_VERSION_OK FALSE)
  else(${EIGEN_VERSION} VERSION_LESS ${Eigen_FIND_VERSION})
    set(EIGEN_VERSION_OK TRUE)
  endif(${EIGEN_VERSION} VERSION_LESS ${Eigen_FIND_VERSION})

  if(NOT EIGEN_VERSION_OK)

    message(STATUS "Eigen version ${EIGEN_VERSION} found in ${EIGEN_INCLUDE_DIR}, "
                   "but at least version ${Eigen_FIND_VERSION} is required")
  endif(NOT EIGEN_VERSION_OK)
endmacro(_eigen3_check_version)

if(NOT EIGEN_FOUND)

  # Look for signature_of_eigen3_matrix_library in build trees or under <prefix>/include/eigen3.
  find_path(EIGEN_INCLUDE_DIR
    NAMES Eigen/Core
    PATH_SUFFIXES eigen3
    HINTS ENV EIGEN_DIR

    PATHS

    # Help the user find it if we cannot.
    DOC "The ${EIGEN_DIR_DESCRIPTION}"
    )

  if(EIGEN_INCLUDE_DIR)
    if(EXISTS ${EIGEN_INCLUDE_DIR}/Eigen/Core)
      set(EIGEN_FOUND 1)
    else()
      set(EIGEN_INCLUDE_DIR "EIGEN_INCLUDE_DIR-NOTFOUND" CACHE PATH "The ${EIGEN_DIR_DESCRIPTION}" FORCE)
    endif()
  endif()

endif()



if (EIGEN_FOUND)

  _eigen3_check_version()
  set(EIGEN_FOUND ${EIGEN_VERSION_OK})

else ()

  # Eigen not found, explain to the user how to specify its location.
  if(EIGEN_FIND_REQUIRED)
    message(FATAL_ERROR ${EIGEN_DIR_MESSAGE})
  else()
    if(NOT EIGEN_FIND_QUIETLY)
      message(STATUS ${EIGEN_DIR_MESSAGE})
    endif()
  endif()

endif()


