cmake_minimum_required(VERSION 3.7)

project(tree_similarity)

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2 -lpthread")
endif()
if(UNIX AND NOT APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O2 -pthread")
endif()

add_library(TreeSimilarity INTERFACE)
target_include_directories(TreeSimilarity INTERFACE
        src/node
        src/label
        src/cost_model
        src/ted
        src/data_structures
        src/parser
        src/join
        src/join/naive
        src/join/tjoin
        src/join/tang
        src/join/tminjoin
        src/join/tpartjoin
        src/join/degree_histogram
        src/join/leaf_dist_histogram
        src/join/label_histogram
        src/join/binary_branches
        src/join/tballjoin
        src/join/histogram
        src/join/guha
        src/tree_generator
        src/ted_ub
        src/ted_lb
        )

# For using add_test().
include(CTest)

# Let CMake know about subdirectories.
add_subdirectory(src/)
add_subdirectory(test/)
