cmake_minimum_required(VERSION 2.8.6)
project(pydifem)

set(CMAKE_CXX_STANDARD 14)

add_compile_options(-fPIC)
add_compile_options(-std=gnu++11)
add_compile_options(-Wdeprecated-declarations)

add_definitions(-DSOFTCON_DIR="${CMAKE_HOME_DIRECTORY}")
set(CMAKE_MODULE_PATH "${CMAKE_HOME_DIRECTORY}/cmake_module")


# CGAL
# set(CGAL_DONT_OVERRIDE_CMAKE_FLAGS true)
# find_package(CGAL COMPONENTS Core)
# include(${CGAL_USE_FILE})


# Tiny Obj Loader
set(TOL_INCLUDE_DIRS
  "${CMAKE_CURRENT_SOURCE_DIR}/extern/tinyObjLoader"
  CACHE PATH "TOL")
include_directories(${TOL_INCLUDE_DIRS})

# Eigen
find_package(Eigen3 3.3.3 REQUIRED)
if (EIGEN3_FOUND)
   include_directories(${EIGEN3_INCLUDE_DIR})
endif()

# Rapidjson
if (NOT RAPIDJSON_INCLUDEDIR)
  set(RAPIDJSON_INCLUDEDIR "${CMAKE_CURRENT_SOURCE_DIR}/extern/rapidjson/include")
endif()
find_package(rapidjson REQUIRED)
include_directories("${RAPIDJSON_INCLUDE_DIRS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${RAPIDJSON_CXX_FLAGS}")


message("Flags : " ${CMAKE_CXX_FLAGS})


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)

link_directories(sim/)
include_directories(sim/)
link_directories(sim/fem/)
include_directories(sim/fem/)
include_directories(render/)
link_directories(render/)

add_subdirectory( sim )
add_subdirectory( render )
add_subdirectory( learn )


add_subdirectory(extern/pybind11)
pybind11_add_module(pydifem python/pydifem.cc ${srcs})

target_compile_definitions(pydifem PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
#add_library(render SHARED ${srcs})
message("boost : " ${Boost_LIBRARIES})
target_link_libraries(pydifem ${OPENGL_LIBRARIES} ${GLUT_LIBRARIES} ${Boost_LIBRARIES} 
    ${PYTHON_LIBRARIES} sim  )
