cmake_minimum_required(VERSION 3.12)

project(fastbc)

add_subdirectory( libfastbc )

#####################################################################
# Fast-BC compile options
#####################################################################
find_package(spdlog REQUIRED CONFIG)

include(CheckIPOSupported)

add_executable(fbc main.cpp)

set_target_properties(fbc PROPERTIES CXX_STANDARD 17 )

check_ipo_supported(RESULT result)
if(result)
	set_target_properties(fbc PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
	message( STATUS "Interprocedural optimization enabled")
endif()

if(CMAKE_BUILD_TYPE MATCHES "[Rr]elease")
	target_compile_definitions(fbc PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO")
else()
	target_compile_definitions(fbc PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG")
endif()

target_link_libraries(fbc 
	PRIVATE
	fastbc 
	spdlog::spdlog )