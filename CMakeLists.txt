cmake_minimum_required(VERSION 3.12)

project(fastbc)

add_subdirectory( libfastbc )

#####################################################################
# Fast-BC compile options
#####################################################################

add_executable(fbc main.cpp)

set(CMAKE_BUILD_TYPE Debug)
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg")
#SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pg")
#SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -pg")

target_include_directories(fbc 
  PRIVATE
  fastbc )

target_link_libraries(fbc 
  fastbc )

set_target_properties(fbc PROPERTIES CXX_STANDARD 17 )
