cmake_minimum_required(VERSION 3.5)
# Project Name
set(PROJECT_NAME_STR generate_network)
set(GEN_NETWORK_VERSION 1.0)
project(${PROJECT_NAME_STR} C CXX)

#
#
set(CMAKE_VERBOSE_MAKEFILE ON)
set(BUILD_SHARED_LIBS ON)

# option to build test
option(BUILD_TESTING "Build tests" ON) # build tests

#
# Requires at least C++17
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Initialize CXXFLAGS and CFLAGS.
set(BASE_OPT_FLAGS "-ffast-math -funroll-loops")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-maybe-uninitialized -fPIC")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} ${BASE_OPT_FLAGS}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "${CMAKE_CXX_FLAGS_RELWITHDEBINFO} ${BASE_OPT_FLAGS}")

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wno-maybe-uninitialized -fPIC")
set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -O0")
set(CMAKE_C_FLAGS_RELEASE "${CMAKE_C_FLAGS_RELEASE} ${BASE_OPT_FLAGS}")
set(CMAKE_C_FLAGS_RELWITHDEBINFO "${CMAKE_C_FLAGS_RELWITHDEBINFO} ${BASE_OPT_FLAGS}")

# Specify default build type
IF(NOT CMAKE_BUILD_TYPE)
  SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
      "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
      FORCE)
ENDIF(NOT CMAKE_BUILD_TYPE)

set(JSON_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/json")
set(FMT_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/fmt")
set(MIO_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/mio")
set(CLI11_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/CLI11")
set(CATCH2_SRC_DIR "${PROJECT_SOURCE_DIR}/ext/Catch2")

# include directories
include_directories("${PROJECT_SOURCE_DIR}/include/")

# Additional cmake files
# list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
# message(STATUS "MODULE PATH   ${CMAKE_MODULE_PATH}")

# Interface to Header files 
add_library(json INTERFACE)
target_include_directories(json INTERFACE ${JSON_SRC_DIR}/include)
add_library(fmt INTERFACE)
target_include_directories(fmt INTERFACE ${FMT_SRC_DIR}/include)
add_library(mio INTERFACE)
target_include_directories(mio INTERFACE ${MIO_SRC_DIR}/include)
add_library(cli11 INTERFACE)
target_include_directories(cli11 INTERFACE ${CLI11_SRC_DIR}/include)
add_library(catch2 INTERFACE)
target_include_directories(catch2 INTERFACE ${CATCH2_SRC_DIR}/include/catch2)

# Check if OpenMP is available
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
  message(STATUS "OpenMP version:  ${OpenMP_CXX_VERSION}")
endif()

# Compile
set(PROJECT_EXEC ${PROJECT_NAME_STR}.x)
set(PROJECT_TEST_EXEC test_${PROJECT_NAME_STR}.x)

# Executable
add_executable(${PROJECT_EXEC} ${FMT_SRC_DIR}/src/format.cc src/generate_network.cpp)
add_compile_definitions(GEN_NETWORK_VERSION=${GEN_NETWORK_VERSION})
# OpenMP
include_directories(${OpenMP_CXX_INCLUDE_DIRS})
target_link_libraries(${PROJECT_EXEC} PUBLIC
    OpenMP::OpenMP_CXX json fmt mio cli11 ${CMAKE_DL_LIBS})

if(${BUILD_TESTING})
    add_executable(${PROJECT_TEST_EXEC} ${FMT_SRC_DIR}/src/format.cc
        ${CATCH2_SRC_DIR}/src/catch_amalgamated.cpp
        tests/test_generate_network.cpp)
    # target_compile_options(${PROJECT_TEST_EXEC}
    #         --coverage -g -O0 -fprofile-arcs -ftest-coverage ${OpenMP_CXX_FLAGS})
    target_link_libraries(${PROJECT_TEST_EXEC} PRIVATE
                          OpenMP::OpenMP_CXX json fmt mio catch2
                          ${CMAKE_DL_LIBS})
endif(${BUILD_TESTING})
