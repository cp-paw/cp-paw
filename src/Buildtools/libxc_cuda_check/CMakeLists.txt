cmake_minimum_required(VERSION 3.20)
project(CheckLibxcCUDA LANGUAGES NONE)
find_package(Libxc 7.1 CONFIG REQUIRED COMPONENTS C Fortran
             PATHS "${XC_PREFIX}" NO_DEFAULT_PATH)
if(NOT Libxc_ENABLE_CUDA)
  message(FATAL_ERROR "The selected Libxc installation was not built with ENABLE_CUDA=ON")
endif()
