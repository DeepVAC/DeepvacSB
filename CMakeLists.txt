# Copyright (c) 2020 DeepVAC
# This file is part of libdeepvac, licensed under the GPLv3 (the "License")
# You may not use this file except in compliance with the License.
cmake_minimum_required(VERSION 3.1.0)

if(${CMAKE_VERSION} VERSION_LESS 3.11)
  cmake_policy(VERSION ${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION})
else()
  cmake_policy(VERSION 3.11)
endif()

project(DeepvacSB CXX)

set(CMAKE_BUILD_TYPE "RELEASE")

file(GLOB SB_SOURCES src/*.cpp)
file(GLOB SB_HEADERS include/*.h)

add_library(DeepvacSB STATIC ${SB_SOURCES})

target_include_directories(DeepvacSB PUBLIC 
    $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

install(TARGETS DeepvacSB EXPORT DeepvacSB DESTINATION "lib")
install(FILES ${SB_HEADERS} DESTINATION "include")