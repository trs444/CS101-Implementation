#include "UndirectedGraph.hpp"
#include "Edge.hpp"
#include "Vertex.hpp"
#include <vector>
#include <iostream>
#include <queue>
#include <limits>
// Method implementations here

using namespace std;

  
  /*
   * inserts an edge into the graph with the given parameters 
   * if the edge already exists in the graph, it updates the edge with the new
   * parameters
   */ 
  void UndirectedGraph::addEdge(const std::string &from, const 
    std::string &to, unsigned int cost, unsigned int length){
    Vertex *fr;
    Vertex *t;
      //checks to see if the "from" vertex is in graph
      //if not found, create new vertex and update vertices map
      if(vertices.find(from) == vertices.end()){
        fr = new Vertex(from);
        std::pair<std::string,Vertex*> vert = std::make_pair(from, fr);
        vertices.insert(vert); 
      }
      //if vertex is found, then make pointer point to vertex
      else{
        std::unordered_map<std::string,Vertex*>::iterator i1 = 
                 vertices.find(from);
        fr = i1->second;
      }
   
      //checks to see if the "to" vertex is in the graph
      //if not found, create new vertex and update vertices map
      if(vertices.find(to) == vertices.end()){
        t = new Vertex(to);
        std::pair<std::string,Vertex*> vert1 = std::make_pair(to,t);
        vertices.insert(vert1);  
      }
      //if vertex is found, then make pointer point to vertex
      else{
        std::unordered_map<std::string,Vertex*>::iterator i2 = 
                vertices.find(to);
        t = i2->second;
      }
      //creates an edge between the two vertices going both ways
      fr->addEdge(t, cost, length);
      t->addEdge(fr, cost, length);
  }

  /*
   * calculates the cost of all the edges on the graph
   */  
  unsigned int UndirectedGraph::totalEdgeCost() const{
    unsigned int t = 0;
    //loop through all th edges and adds up the cost
    for(auto it = vertices.begin(); it!=vertices.end(); ++it){
      t += it->second->totalEdgeCost();
    }
    //because it is undirected, divide by 2
    return t/2;
  
  }

  /*
   *  removes all of the edges from the graph except those necessary to form
   *  a minimum spanning tree using Prim's algortihm
   */ 
  UndirectedGraph UndirectedGraph::minSpanningTree(){
    //priority queue to store the edges
    std::priority_queue<Edge,std::vector<Edge>, EdgeComp> pq;
    //graph used to make the tree
    UndirectedGraph t;
    //set all the visited flags to false
    for(auto it = vertices.begin(); it!=vertices.end(); ++it)
      it->second->setVisited(false);
    //first vertex to visit
    Vertex *first = vertices.begin()->second;
    first->setVisited(true);
    //push all the edges into the priority queue
    for(auto it2 = first->getEdges().begin(); 
              it2!=first->getEdges().end(); ++it2){
        pq.push(it2->second);
    }
    //sort through the priority queue until the priority queue is empty
    while(!pq.empty()){
      //gets the top of the queue and takes it out of the queue
      Edge exp = pq.top();
      pq.pop();
      //if it has not been visited, set the visited flag to true, add an edge
      //to the min spanning tree and pushes it to the priority queue
      if(exp.getTo()->wasVisited() == false){
        exp.getTo()->setVisited(true);
        t.addEdge(exp.getFrom()->getName(), exp.getTo()->getName(), 
                    exp.getCost(), exp.getLength());
        for(auto it3 = exp.getTo()->getEdges().begin(); 
                 it3!=exp.getTo()->getEdges().end(); ++it3){
            pq.push(it3->second);
        }
      }
    }
    return t;
  }
  
  /*
   * determines the combined distance from the given vertex to all other 
   * vertices using Dijkstra's algorithm
   */ 
  unsigned int UndirectedGraph::totalDistance(const std::string &from){
    //if the vertex is not in the graph or it is not reachable, return max
    //distance
    if(vertices[from] == NULL || vertices[from]->getEdges().empty()){
      return std::numeric_limits<int>::max();
    }
    //priority queue to store the vertices
    std::priority_queue<std::pair<Vertex*,unsigned int>, 
     std::vector<std::pair<Vertex*,unsigned int>>,DijkstraVertexComparator> pq;
    //loop through all the vertices and sets the visited flags to false and 
    //distance to 0
    for(auto it = vertices.begin(); it!=vertices.end(); ++it){
      it->second->setVisited(false);
      it->second->setDistance(std::numeric_limits<int>::max());
    }
    //finds the wanted vertex and sets the pointer to it
    Vertex *first = vertices.find(from)->second;
    //sets the first vertex to zero
    first->setDistance(0);
    //pair used to place into priority queue
    std::pair<Vertex*, unsigned int> v;
    //distance variable used to calcualte the distance
    unsigned int newDistance = 0;
    //push the first vertex onto pqueue
    pq.push(std::make_pair(first, first->getDistance()));
    //loop until pqueue is empty
    while( !pq.empty()){
      //grabs the top vertex of the pqueue 
      v = pq.top();
      pq.pop();
      //loops through the vertex's edges and finds the distance of it
      //and adds the edge distance to the first one
      for(auto it = v.first->getEdges().begin(); 
                it!=v.first->getEdges().end(); ++it){
        newDistance = v.first->getDistance() + it->second.getLength();
        //update the distance of the vertex if the newDistance is
        //less than the old one then pushes the vertex back onto the pqueue
        if(newDistance < it->second.getTo()->getDistance()){
          it->second.getTo()->setDistance(newDistance);
          std::pair<Vertex*,unsigned int>vert1 = 
                  std::make_pair(it->second.getTo(),newDistance);
          pq.push(vert1);
        }
      }
    }
    //combined distance of all of the vertices
    unsigned int combDist = 0;
    for(auto it2 = vertices.begin(); it2!= vertices.end(); ++it2){
      combDist += it2->second->getDistance();
    }
    return combDist;
  }

  //total distance of all the vertices to all the other vertices in the grpah
  unsigned int UndirectedGraph::totalDistance(){
    unsigned int maxDist = 0;
    for(auto it = vertices.begin(); it != vertices.end(); ++it){
      maxDist += totalDistance(it->first);
    }
    return maxDist;
  }
