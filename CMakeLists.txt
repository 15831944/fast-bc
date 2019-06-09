cmake_minimum_required(VERSION 3.12)

project(fastbc)

add_subdirectory( libfastbc )

#####################################################################
# Fast-BC compile options
#####################################################################
find_package(spdlog REQUIRED CONFIG)

add_executable(fbc main.cpp)

set_target_properties(fbc PROPERTIES CXX_STANDARD 17 )

if(CMAKE_BUILD_TYPE MATCHES "[Dd]ebug")
	target_compile_definitions(fbc PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG")
else()
	target_compile_definitions(fbc PRIVATE "SPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_WARN")
endif()

target_link_libraries(fbc 
	PRIVATE
	fastbc 
	spdlog::spdlog )