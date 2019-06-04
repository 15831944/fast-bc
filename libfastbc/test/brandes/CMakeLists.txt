#########################################################################################
#	Brandes tests directory
#########################################################################################

target_sources(fastbctests PRIVATE 
    brandes/DijkstraClusterEvaluator.cpp
	brandes/VertexInfo.cpp
	brandes/VertexInfoPivotSelector.cpp
	brandes/DijkstraSSBrandesBC.cpp )