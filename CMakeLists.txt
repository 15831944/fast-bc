cmake_minimum_required(VERSION 3.12)

project(fastbc)

add_subdirectory( libfastbc )

#####################################################################
# Fast-BC compile options
#####################################################################

add_executable(fbc main.cpp)

find_package(OpenMP REQUIRED)

target_include_directories(fbc 
  PRIVATE
  fastbc )

target_link_libraries(fbc 
  fastbc ${OpenMP_CXX_LIBRARIES})

set_target_properties(fbc PROPERTIES CXX_STANDARD 17 )
