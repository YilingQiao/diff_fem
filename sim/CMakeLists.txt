cmake_minimum_required(VERSION 2.8.6)
project(sim)

add_compile_options(-std=gnu++11)
add_compile_options(-Wdeprecated-declarations)
SET(CMAKE_BUILD_TYPE Release CACHE STRING
	"Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
#	FORCE
	)

set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR})


find_package(TinyXML REQUIRED)
find_package(OpenMP REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(PythonLibs REQUIRED)
find_package(Boost 1.66.0 COMPONENTS filesystem python3 numpy3)



add_compile_options(-DHAVE_CSTDDEF)
add_definitions(-DOCTOCON_DIR="${CMAKE_HOME_DIRECTORY}")
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${TinyXML_INCLUDE_DIRS})
include_directories(${PYTHON_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIRS})

include_directories(../sim/)
link_directories(./fem/)
link_directories(./utils/)


file(GLOB_RECURSE srcs "*.h" "*.cpp" "*.hpp")
file(GLOB csrcs "*.cpp" )
file(GLOB TDS_HDRS "*tiny_*.h" )
include_directories(${TDS_HDRS})

# message("source files: " ${srcs})
# message("csource files: " ${csrcs})
# message("tiny files: " ${TDS_HDRS})

add_library(sim ${srcs})
target_link_libraries(sim OpenMP::OpenMP_CXX ${Boost_LIBRARIES} ${PYTHON_LIBRARIES} #[[${DART_LIBRARIES}]] ${TinyXML_LIBRARIES} )

add_subdirectory(fem)

# add_executable(test_tiny test_tiny.cpp ${TDS_HDRS})
# target_link_libraries(test_tiny sim )