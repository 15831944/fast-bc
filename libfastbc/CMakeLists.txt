#####################################################################
#   Fast-bc library directory
#####################################################################

find_package(spdlog REQUIRED CONFIG)
find_package(OpenMP REQUIRED)

file(GLOB_RECURSE LIBFASTBC_SRCS "src/*")

add_library(fastbc INTERFACE)

target_include_directories(fastbc INTERFACE "include")

target_link_libraries(fastbc 
	INTERFACE 
	spdlog::spdlog 
	OpenMP::OpenMP_CXX )

add_subdirectory(test)