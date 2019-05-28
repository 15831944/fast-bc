cmake_minimum_required(VERSION 3.12)

project(fastbc)

add_subdirectory( libfastbc )

#####################################################################
# Fast-BC compile options
#####################################################################

# Remove after debug
set(CMAKE_BUILD_TYPE Debug)

set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost 1.68.0 COMPONENTS graph) 

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 

    add_executable(fbc main.cpp)

    target_include_directories(fbc PUBLIC fastbc)
    
    target_link_libraries(fbc ${Boost_LIBRARIES} fastbc)

    set_target_properties(fbc PROPERTIES CXX_STANDARD 17 CXX_FLAGS "-march=native -fopt-info -ffast-math -ftree-vectorizer-verbose=2" )
endif()

