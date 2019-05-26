cmake_minimum_required(VERSION 3.12)

project(fastbc)

add_subdirectory( libfastbc )

#####################################################################
# Fast-BC compile options
#####################################################################

add_executable(fbc main.cpp)

target_include_directories(fbc 
  PUBLIC
  fastbc )

target_link_libraries(fbc 
  fastbc )

set_target_properties(fbc PROPERTIES 
  CXX_STANDARD 17
  CXX_FLAGS "-march=native -fopt-info -ffast-math -ftree-vectorizer-verbose=2" )
