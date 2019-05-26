#####################################################################
#   Fast-bc library directory
#####################################################################

file(GLOB_RECURSE LIBFASTBC_SRCS "src/*")

add_library(fastbc STATIC ${LIBFASTBC_SRCS})

set_target_properties(fastbc PROPERTIES 
  CXX_STANDARD 17
  CXX_FLAGS "-march=native -fopt-info -ffast-math" )

target_include_directories(fastbc PUBLIC include/)

#target_link_libraries(fastbc PUBLIC OpenMP::OpenMP_CXX)
#add_subdirectory(test)