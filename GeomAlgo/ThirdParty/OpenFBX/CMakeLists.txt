cmake_minimum_required(VERSION 3.2 FATAL_ERROR)
cmake_policy(VERSION 3.2)

project(OpenFBX)

include_directories(${CMAKE_SOURCE_DIR}/ThirdParty/OpenFBX/src)

add_library(OpenFBX
    src/ofbx.cpp
    src/ofbx.h
    src/miniz.c
    src/miniz.h
)