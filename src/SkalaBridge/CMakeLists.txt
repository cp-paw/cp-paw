cmake_minimum_required(VERSION 3.20)

project(CPPAWSkalaBridge VERSION 0.1.0 LANGUAGES CXX Fortran)

include(FetchContent)
include(GNUInstallDirs)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

if(APPLE)
  set(CPPAW_SKALA_LIBRARY_RPATH "@loader_path")
  set(CPPAW_SKALA_EXECUTABLE_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
  set(CPPAW_SKALA_LIBRARY_RPATH "$ORIGIN")
  set(CPPAW_SKALA_EXECUTABLE_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
  list(APPEND CPPAW_SKALA_LIBRARY_RPATH ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
  list(APPEND CPPAW_SKALA_EXECUTABLE_RPATH ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
endif()

set(CPPAW_SKALA_FTORCH_URL
    "https://github.com/Cambridge-ICCS/FTorch/archive/refs/tags/v1.0.0.tar.gz"
    CACHE STRING "Pinned FTorch source archive")
set(CPPAW_SKALA_FTORCH_SHA256
    "c4b6741e582623b7ecaecd59d02f779e8a6f6017f8068c85da8a034f468df375"
    CACHE STRING "SHA-256 of the pinned FTorch source archive")
set(CPPAW_SKALA_DEVICE "CPU" CACHE STRING "Torch device compiled by FTorch (CPU or CUDA)")
set_property(CACHE CPPAW_SKALA_DEVICE PROPERTY STRINGS CPU CUDA)

string(TOUPPER "${CPPAW_SKALA_DEVICE}" CPPAW_SKALA_DEVICE)
if(NOT CPPAW_SKALA_DEVICE MATCHES "^(CPU|CUDA)$")
  message(FATAL_ERROR "CPPAW_SKALA_DEVICE must be CPU or CUDA")
endif()

find_package(Torch REQUIRED)
find_package(FTorch CONFIG QUIET)
if(NOT FTorch_FOUND)
  message(STATUS "FTorch was not found; fetching the pinned v1.0.0 source archive")
  set(CMAKE_BUILD_TESTS OFF CACHE BOOL "Build FTorch tests" FORCE)
  if(CPPAW_SKALA_DEVICE STREQUAL "CUDA")
    set(GPU_DEVICE CUDA CACHE STRING "FTorch target device" FORCE)
  else()
    set(GPU_DEVICE NONE CACHE STRING "FTorch target device" FORCE)
  endif()
  FetchContent_Declare(
    ftorch
    URL "${CPPAW_SKALA_FTORCH_URL}"
    URL_HASH "SHA256=${CPPAW_SKALA_FTORCH_SHA256}"
    DOWNLOAD_EXTRACT_TIMESTAMP ON)
  FetchContent_MakeAvailable(ftorch)
endif()
if(TARGET ftorch AND CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
  set_property(TARGET ftorch APPEND PROPERTY
    INSTALL_RPATH "${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES}")
endif()

add_library(cppaw_skala_bridge SHARED
  paw_skala_bridge.f90
  paw_skala_protocol.cxx)
target_link_libraries(cppaw_skala_bridge
  PUBLIC FTorch::ftorch
  PRIVATE ${TORCH_LIBRARIES})
target_include_directories(cppaw_skala_bridge
  PRIVATE "${TORCH_INCLUDE_DIRS}"
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/modules>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
set_target_properties(cppaw_skala_bridge PROPERTIES
  Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/modules"
  INSTALL_RPATH "${CPPAW_SKALA_LIBRARY_RPATH}")
if(CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
  target_compile_definitions(cppaw_skala_bridge PRIVATE CPPAW_SKALA_NVHPC_FORTRAN=1)
endif()

if(TORCH_CXX_FLAGS)
  separate_arguments(CPPAW_TORCH_CXX_FLAGS NATIVE_COMMAND "${TORCH_CXX_FLAGS}")
  target_compile_options(cppaw_skala_bridge PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:${CPPAW_TORCH_CXX_FLAGS}>)
endif()

add_executable(cppaw_skala_smoke skala_bridge_smoke.f90)
target_link_libraries(cppaw_skala_smoke PRIVATE cppaw_skala_bridge)
set_target_properties(cppaw_skala_smoke PROPERTIES
  INSTALL_RPATH "${CPPAW_SKALA_EXECUTABLE_RPATH}")
if(CMAKE_SYSTEM_NAME STREQUAL "Linux" AND
   CMAKE_Fortran_COMPILER_ID MATCHES "^(NVHPC|PGI)$")
  target_link_options(cppaw_skala_smoke PRIVATE "LINKER:--disable-new-dtags")
endif()

install(TARGETS cppaw_skala_bridge cppaw_skala_smoke
  LIBRARY DESTINATION "${CMAKE_INSTALL_LIBDIR}"
  RUNTIME DESTINATION "${CMAKE_INSTALL_BINDIR}")
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/modules/paw_skala_bridge.mod"
  DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}")
