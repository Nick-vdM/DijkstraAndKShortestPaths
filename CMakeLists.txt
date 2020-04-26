cmake_minimum_required(VERSION 3.15)
project(3_graphs)

set(CMAKE_CXX_STANDARD 14)

add_executable(3_graphs main.cpp)
add_executable(yan-qi k-shortest-paths-cpp-version/src/DijkstraShortestPathAlg.cpp
        k-shortest-paths-cpp-version/src/DijkstraShortestPathAlg.h
        k-shortest-paths-cpp-version/src/Graph.cpp
        k-shortest-paths-cpp-version/src/Graph.h
        k-shortest-paths-cpp-version/src/GraphElements.h
        k-shortest-paths-cpp-version/src/MainP.cpp
        k-shortest-paths-cpp-version/src/YenTopKShortestPathsAlg.cpp
        k-shortest-paths-cpp-version/src/YenTopKShortestPathsAlg.h)