find_package(Catch2 CONFIG REQUIRED)
include(CTest)
include(Catch)

#########################################################################################
#	fastbc library tests compile options
#########################################################################################

add_executable(fastbctests 
	test.cpp
	DirectedWeightedGraph.cpp
	SubGraph.cpp )

set_property(TARGET fastbctests PROPERTY CXX_STANDARD 17)

target_include_directories(fastbctests INTERFACE fastbc)

target_link_libraries(fastbctests 
    Catch2::Catch2
    fastbc )

configure_file("${CMAKE_CURRENT_SOURCE_DIR}/DWGtext.txt" 
	"${CMAKE_CURRENT_BINARY_DIR}/DWGtext.txt" COPYONLY)


#########################################################################################
#	fastbc library tests subdirectories
#########################################################################################

add_subdirectory(brandes)

catch_discover_tests(fastbctests)