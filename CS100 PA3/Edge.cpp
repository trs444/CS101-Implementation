#include "Edge.hpp"

// Method implementations here

using namespace std;

  /*
   * Constructs an edge between two vertices and sets the cost and length
   */  
  Edge::Edge( Vertex *from, Vertex *to, 
        unsigned int cost,
        unsigned int length)
    : from(from), to(to), cost(cost), length(length) { }

  //returns a pointer to the origin of the edge
  Vertex* Edge::getFrom() const{
    return from;
  }

  //returns a pointer to where the edge ends
  Vertex* Edge::getTo() const{
    return to;
  }

  //sets the cost of the edge
  void Edge::setCost(unsigned int c){
    cost = c;
  }

  // returns the cost of the edge
  unsigned int Edge::getCost() const{
    return cost;
  }

  //sets the length of the edge
  void Edge::setLength(unsigned int l){
    length = l;
  }

  //gets the length of the edge
  unsigned int Edge::getLength() const{
    return length;
  }

  //compares the edges to see which edge is lesser and puts the lesser one
  //as priority
  bool Edge::operator<(const Edge &right) const{
    if(cost != right.cost)
      return right.cost < cost;

    return true;
  } 
  
