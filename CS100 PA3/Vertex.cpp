#include "Vertex.hpp"
#include "Edge.hpp"
#include <unordered_map>
#include <iostream>
#include <vector>

// Method implementations here

using namespace std;
  
  /**
   * Vertex constructor. Initializes vertex with paramater name
   */ 
  Vertex::Vertex(const std::string &name):
    name(name) {}

  /**
   * adds an edge from this vertex to the parameter "to" vertex
   */ 
  bool Vertex::addEdge(Vertex *to, unsigned int cost, 
    unsigned int length){
    //creates new edge
    Edge e1 = Edge(this, to, cost, length); 
    //inserts the edge into edges map
    edges.insert(std::make_pair(to->getName(), e1)); 
    return true;
  }

  /**
   * returns the name of the vertex
   */
  const std::string & Vertex::getName() const{
    return name;
  }

  //returns the distance of the vertex
  unsigned int Vertex::getDistance() const{
    return distance;
  }

  //sets the distance of the vertex
  void Vertex::setDistance(unsigned int d){
    distance = d;
  }

  //checks to see if vertex was visited.
  //returns true if it was, false if not
  bool Vertex::wasVisited() const{
    return visited;
  }

  //sets the vertex's visited state to true if visited false if not
  void Vertex::setVisited(bool v){
    visited = v;
  }
  
  //clears all the edges of the edge map
  void Vertex::clearEdges(){
    edges.clear();
  }

  //calculates the total cost of all the edges at the vertex
  unsigned int Vertex::totalEdgeCost() const{
    unsigned int t = 0;
    std::unordered_map<std::string, Edge>::iterator it1;
    //loop through all the edges of the vertex
    for(auto it1 = edges.begin(); it1!=edges.end(); ++it1){
      t += it1->second.getCost();
    }
    return t;
  }

  //returns the edge map
  const std::unordered_map<std::string, Edge> & Vertex::getEdges() const{
    return edges;    
  }
