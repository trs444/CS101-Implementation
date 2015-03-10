#include "MST.h"
#include "point.h"
#include <stack> //std::stack
#include <algorithm> //std::find

MST::MST(float** input, int size) {
	adjacentMatrix = input;
	key = new int[size];   
  mstSet = new bool[size];  
	parent = new int[size];
  oddDegrees = new int[size];
	N = size;
  isOdd = new bool[N];
}

MST::~MST() {

}

//use Prim's algorithm or Kruskal algorithm. Copied from 'http://www.geeksforgeeks.org/greedy-algorithms-set-5-prims-minimum-spanning-tree-mst-2/'
void MST::makeTree() { 
     // Initialize all keys as INFINITE
     for (int i = 0; i < N; i++)
        key[i] = INT_MAX, mstSet[i] = false;
 
     // Always include first 1st vertex in MST.
     key[0] = 0;     // Make key 0 so that this vertex is picked as first vertex
     parent[0] = -1; // First node is always root of MST 
 
     // The MST will have V vertices
     for (int count = 0; count < N-1; count++)
     {
        // Pick the minimum key vertex from the set of vertices
        // not yet included in MST
        int u = minKey(key, mstSet);
 
        // Add the picked vertex to the MST Set
        mstSet[u] = true;
 
        // Update key value and parent index of the adjacent vertices of
        // the picked vertex. Consider only those vertices which are not yet
        // included in MST
        for (int v = 0; v < N; v++)
           // mstSet[v] is false for vertices not yet included in MST
           // Update the key only if adjacentMatrix[u][v] is smaller than key[v]
          if (adjacentMatrix[u][v] && mstSet[v] == false && adjacentMatrix[u][v] <  key[v])
             parent[v]  = u, key[v] = adjacentMatrix[u][v];
     }
}

// A utility function to find the vertex with minimum key value, from
// the set of vertices not yet included in MST
int MST::minKey(int key[], bool mstSet[])
{
   // Initialize min value
   int min = INT_MAX, min_index;
 
   for (int v = 0; v < N; v++)
     if (mstSet[v] == false && key[v] < min)
         min = key[v], min_index = v;
 
   return min_index;
}

// A utility function to print the constructed MST stored in parent[]
void MST::printMST() {
	cout<<endl;
	cout<<"Minimum spanning tree from the adjacency matrix"<<endl;
	cout<<"Edge   Weight"<<endl;
	for (int i = 1; i < N; i++) {
		cout<<parent[i]<<" - "<<i<<"  "<<adjacentMatrix[i][parent[i]]<<endl;
	}
}

//calculate mean of all edges in the MST
float MST::calMean(int option) {
float mean = 0.0;
float sum = 0;
  if(option == MST_1) {
    //calculate
    for( int i = 0; i < N; i++)
       sum = sum + key[i];
       mean = sum/(N-1);
       cout << "Mean: " << mean << endl;
  }else if(option == TSP2) {
	} else if(option == TSP1_5) {

	}

	return mean;
}

//calculate standard deviation of all edges in the MST
float MST::calStd(int option) {
	float std = 0.0;

  float mean = calMean(option);
    if(option == MST_1) {
      //calculate
      for(int i = 0; i < N; i++)
         std += (key[i]-mean) * (key[i]-mean);
         std = sqrt(std/N-1);
         cout << "std: " << std << endl;
  } else if(option == TSP2) {
	} else if(option == TSP1_5) {

	}

	return std;
}

float MST::makeTSP2(int startV, Point p) {
	//make a Eulerian tour by DFS
  cout << "ENTERING MAKETSP2" << endl;


  int visited [this->N]; //keeps track of what nodes we've seen
  for(int zz = 0; zz < this->N; zz++) {
    visited[zz] = -99;
  }
  size_t visitedSize = sizeof(visited)/sizeof(int);
  int *visitedEnd = visited + visitedSize; //size of array
  std::stack<int> order;

  order.push(startV); //push starting node on stack

  std::vector<pair<int, int>> mstVect;
  std::pair<int, int> mstPair;
  int j = 0;
  for( int i = 0; i < this->N; i++) {
    mstPair = std::make_pair(parent[i], i);
    mstVect.push_back(mstPair);
  }
 
  while(!order.empty()) { //when stack isn't emtpy, pop it
    startV = order.top();
    order.pop();

    if(std::find(visited, visitedEnd, startV) == visitedEnd) {
      for(std::vector<pair<int, int>>::iterator it = mstVect.begin(); it != mstVect.end(); ++it) {
        if(it->first == startV && std::find(visited, visitedEnd, it->second) == visitedEnd) {
           order.push(it->second);
          }
        if(it->second == startV && std::find(visited, visitedEnd, it->first) == visitedEnd) {
           order.push(it->first);
          }
        }
      visited[j] = startV;
      j++;
    }
  }

  //print order visited
  cout << "NODES VISITED IN THIS ORDER:" << endl;
  for(int k = 0; k < (sizeof(visited)/sizeof(int)); k++) {
    cout << visited[k] << endl;
  }

  int coordinates [2*(this->N)]; //Store coordinates for points
  int g = 0;

  float cost = 0; //initialize tsp2 cost

  //add shortcuts if a vertex has no detours
  //calculate heuristic TSP cost

  for(int m = 0; m < this->N-1; m++) {
    cost += adjacentMatrix[visited[m]][visited[m+1]];
    //cout << "edge cost between " << visited[m] << " and " << visited[m+1] << " = " << p.getEuclideanDistance(coordinates[visited[m]*2], coordinates[(visited[m]*2)+1], coordinates[visited[m+1]*2], coordinates[(visited[m+1]*2)+1]) << endl;
    
  }
  cost += (float)adjacentMatrix[visited[0]][visited[this->N-1]];
  //cout << "final edge cost is (correct) " << adjacentMatrix[visited[0]][visited[this->N-1]] << endl;
                                                                                   

  int keySum = 0;
  for(int z = 0; z < this->N; z++) {
    keySum += key[z];
  }



  cout << "MST TOTAL COST IS: " << keySum <<endl;
  cout << "TSP2 COST IS: " << cost << endl;
  return cost;

}



int MST::numOfOdd() {
 
  
  int odd = 0;
  for(int i = 0; i < N; i++) {
    if(isOdd[i] == 1) {
      odd++;
    }
  }
  return odd;

}

int* MST::whichOddVertices() {
  int oddNum = numOfOdd();

  oddDegrees[oddNum];
  int j = 0;
  for( int i = 0; i < N; i++) {
    if(isOdd[i] == true) {
      oddDegrees[j] = i;
      j++;
    }
  }
  return oddDegrees;
  

}
void MST::findOddNode() {
  int * count = new int[N];
  for( int i = 0; i < N; i++) {
    count[i] = 0;
    isOdd[i] = false;
  }
  for( int i = 1; i < N; i++) {
    count[i]++;
  }
  for(int i = 0; i < N; i++) {
    count[parent[i]]++;
  }
  for( int i = 0; i < N; i++) {
    if( count[i] % 2 != 0 ) {
      isOdd[i] = true;
    }
  }
  /*for( int i = 0; i < N; i++) {
    cout << isOdd[i] << endl;
  }*/
}


void MST::makeTSP1_5() {
	
	//construct minimum-weight-matching for the given MST
	minimumMatching();

	//make all edges has even degree by combining mimimum-weight matching and MST
	combine();

	//calculate heuristic TSP cost 
}

void MST::minimumMatching() { //if you choose O(n^2)
	//find minimum-weight matching for the MST. 
	
	//you should carefully choose a matching algorithm to optimize the TSP cost.
}

std::vector<pair<int, int>> MST::combine() {
	//combine minimum-weight matching with the MST to get a multigraph which has vertices with even degree
  std::vector<pair<int, int>> result;

  /*for (int i = 1; i < N; i++) {
    cout<<parent[i]<<" - "<<i<<"  "<<adjacentMatrix[i][parent[i]]<<endl;
  }*/
  std::pair<int, int> bestPair;
  for( int i = 1; i < N; i++) {
    bestPair = std::make_pair(parent[i], i);
    result.push_back(bestPair);
  }



  return result;













}
