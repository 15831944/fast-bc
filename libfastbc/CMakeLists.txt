#####################################################################
#   Fast-bc library directory
#####################################################################

file(GLOB_RECURSE LIBFASTBC_SRCS "src/*")

add_library(fastbc INTERFACE )

target_include_directories(fastbc INTERFACE "include")

