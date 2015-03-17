#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "Vertex.hpp"
#include "Edge.hpp"
#include "UndirectedGraph.hpp"
#include <limits>

/**
 * Entry point into the netplan program.
 *
 * -Reads a file from the filesystem according to the specification for
 *  PA3, creating an UndirectedGraph.
 * -Finds the total cost & ping time of the graph as presented in the input
 *  file.
 * -Determines the minimum cost graph from the original graph.
 * -Finds the total cost & ping time of the minimum cost graph.
 * -Finds the change of cost & ping time from the original graph to the
 *  minimum cost graph.
 * -Prints the results to stdout.
 *
 * Usage:
 *   ./netplan infile
 *
 */
int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " infile" << std::endl;
        return EXIT_FAILURE;
    }
    
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "Unable to open file for reading." << std::endl;
        return EXIT_FAILURE;
    }
    
    // Implementation here
    
    UndirectedGraph g;
    std::string v1;
    std::string v2;
    double cost;
    double time;
    while(in.good()){
      in >> v1 >> v2 >> cost >> time;
      if(!in.good()) break;
      g.addEdge(v1, v2, cost, time);
    }

    //minimum spanning tree    
    UndirectedGraph minTree;
    //total cost of all edges for graph g
    std::cout << g.totalEdgeCost() << std::endl;
    //minimum spanning tree for g
    minTree = g.minSpanningTree();
    //total cost of all edges for minimum spanning tree minTree
    std::cout << minTree.totalEdgeCost() << std::endl;
    //difference of total edge cost for normal graph and minimum spanning tree
    std::cout << g.totalEdgeCost() - minTree.totalEdgeCost() << std::endl;
    //Dijkstra's algorithm to calculate total distance on original graph
    std::cout << g.totalDistance() << std::endl;
    //Dijkstra's algorithm to calculate total distance on minimum spanning tree
    std::cout << minTree.totalDistance() << std::endl;
    //difference between total time costs.
    std::cout << minTree.totalDistance() - g.totalDistance() << std::endl;
    return EXIT_SUCCESS;
}
