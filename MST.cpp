#include "MST.h"
#include "point.h"
#include <stack> //std::stack
#include <algorithm> //std::find

MST::MST(float** input, int size) {
	adjacentMatrix = input;
	key = new int[size];   
  mstSet = new bool[size];  
	parent = new int[size];
	N = size;
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
  //cout << "visited size is " << visitedSize << endl;
  int *visitedEnd = visited + visitedSize; //size of array
  std::stack<int> order;
  set< pair<int,int> > points = p.getPointset();

  order.push(startV); //push starting node on stack

  std::vector<pair<int, int>> mstVect;
  std::pair<int, int> mstPair;
  int j = 0;
  for( int i = 0; i < this->N; i++) {
    mstPair = std::make_pair(parent[i], i);
    mstVect.push_back(mstPair);
  }
 
  //cout << "ENTERING WHILE LOOP startV = " << startV << endl;
  while(!order.empty()) { //when stack isn't emtpy, pop it
    startV = order.top();
    //cout << "stack: " << order.top() << endl;
    order.pop();

    //cout << "startV = " << startV <<endl;

    //cout << "Entering if statement 1" << endl;

    if(std::find(visited, visitedEnd, startV) == visitedEnd) {
      //visited[j] = startV;
      //cout << visited[j] << " IS BEING VISITED**********************************************" << endl;
      //j++;
      //for(int i = 0; i < (this->N); i++) { //wrong?

       // cout << "Entering if statement 2" << endl;

        /*if((this->adjacentMatrix[startV][i] != 0) && (this->parent[i] == startV)) {
          order.push(i);
          //cout << "node " << i << " pushed to stack from node " << startV << endl;
        }*/


        //for(int ab = 0; ab <this->N; ab++) {
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
        //}

      //}
      //iterate through and find which vertices startV is connected to
    }

  }

  //print order visited
  cout << "NODES VISITED IN THIS ORDER:" << endl;
  for(int k = 0; k < (sizeof(visited)/sizeof(int)); k++) {
    //cout << visited[k] << endl;
  }

  int coordinates [2*(this->N)]; //Store coordinates for points
  int g = 0;

  for(set< pair<int,int> >::iterator it = points.begin() ; it != points.end() ; ++it) {
    coordinates[g] = (*it).first; //x-coord first
    coordinates[g+1] = (*it).second; //y-coord
    g += 2;
  }

  float cost = 0; //initialize tsp2 cost

  //add shortcuts if a vertex has no detours
  //calculate heuristic TSP cost
  std::pair<int, int> pair1;
  std::pair<int, int> pair2;

  for(int m = 0; m < this->N-1; m++) {
    //cost += p.getEuclideanDistance(coordinates[visited[m]*2], coordinates[(visited[m]*2)+1], coordinates[visited[m+1]*2], coordinates[(visited[m+1]*2)+1]);
    cost += adjacentMatrix[visited[m]][visited[m+1]];
    //cout << "edge cost between " << visited[m] << " and " << visited[m+1] << " = " << p.getEuclideanDistance(coordinates[visited[m]*2], coordinates[(visited[m]*2)+1], coordinates[visited[m+1]*2], coordinates[(visited[m+1]*2)+1]) << endl;
    
  }
  cost += (float)adjacentMatrix[visited[0]][visited[this->N-1]];
  //cout << "final edge cost is (correct) " << adjacentMatrix[visited[0]][visited[this->N-1]] << endl;
  //cost += p.getEuclideanDistance(coordinates[0], coordinates[1], coordinates[(visited[this->N-1]*2)], coordinates[visited[((this->N-1)*2)+1]]);
  //cout << "cost of final edge is " << p.getEuclideanDistance(coordinates[0], coordinates[1], coordinates[(visited[this->N-1]*2)], coordinates[visited[((this->N-1)*2)+1]]) << endl;
    //cout << "TSP2 Cost: " << cost << endl;
    //cout << "blah" << endl;
    /*for( std::vector<pair<int, int>>::iterator it = mstVect.begin(); it != mstVect.end(); ++it) {
      pair1 = std::make_pair(visited[m], visited[m+1]);
      pair2 = std::make_pair(visited[m+1], visited[m]);
      //cout << "pair1 first: " << pair1.first << endl;
      //cout << "pair1 second: " << pair1.second << endl;
      //cout << "visited[m]" << visited[m] << endl;
      //cout << "visited[m+1] " << visited[m+1] << endl;
      if(mstVect[m] == pair1 || mstVect[m] == pair2) {
        cost += (float)adjacentMatrix[pair1.first][pair1.second];
        //cost += (float)key[visited[m+1]];
        break;
      }
    /*if(parent[visited[m+1]] == visited[m]) {
      cost += (float)key[visited[m+1]];
      cout << "if" << endl;
    //}
    else{
      //cout << "else" << endl;
      //cout << "no edge, creating" << endl;
      cost += p.getEuclideanDistance(coordinates[visited[m]*2], coordinates[(visited[m]*2)+1], coordinates[visited[m+1]*2], coordinates[(visited[m+1]*2)+1]);
      break;
      
      //cout << coordinates[visited[m]*2] << endl;
      //cout << coordinates[(visited[m]*2)+1] << endl;
      //cout << coordinates[visited[m+1]*2] << endl;
      //cout << coordinates[(visited[m+1]*2)+1] << endl;
    }
  }*/
  //}

  //check case if there's an edge form first to last node


  /*for(int n = 0; n < this->N-1; n++) {
    for( std::vector<pair<int, int>>::iterator it = mstVect.begin(); it != mstVect.end(); ++it) {
        pair1 = std::make_pair(visited[0], visited[this->N-1]);
        pair2 = std::make_pair(visited[this->N-1], visited[0]);
        if(mstVect[n] == pair1 || mstVect[n] == pair2) {
          cout << "THERE'S ALREADY AN EXISTING EDGE FROM THE FIRST TO THE LAST NODE" << endl;
        }
        else{                                                                       
          //cost += p.getEuclideanDistance(coordinates[0], coordinates[1], coordinates[(this->N*2)-2], coordinates[(this->N*2)-1]);
          break;
        }     
    }
  }*/  


  /*if(parent[visited[this->N-1]] == visited[0]) {
    cout << "THERE'S ALREADY AN EXISTING EDGE FROM THE FIRST TO THE LAST NODE" << endl;
  }      
  else{      
        cout << "blah" << endl;                                                                 
        //cost += p.getEuclideanDistance(coordinates[0], coordinates[1], coordinates[(this->N*2)-2], coordinates[(this->N*2)-1]);
      }*/                                                                       
            

  int keySum = 0;
  for(int z = 0; z < this->N; z++) {
    //cout << "key at " << z << " is: " << key[z] << endl;
    keySum += key[z];
  }



  cout << "MST TOTAL COST IS: " << keySum <<endl;
  cout << "TSP2 COST IS: " << cost << endl;
  return cost;

}



int MST::findOddVertices() {
 
  int degree [this->N];

  for(int i = 0; i < this->N; i++) {
    degree[i] = 1;
  }

  int currNode;
  for(int j = 2; j < this->N; j++) {
    currNode = parent[j];
    degree[currNode] += 1;
  }

  int numOddDegree = 0;
  for(int k = 0; k < this->N; k++) {
    if(degree[k]%2 == 1) {
      numOddDegree++;
    }
  }
  return numOddDegree;


}

int* MST::whichOddVertices() {
  
  int degree [this->N];
  
  for(int i = 0; i < this->N; i++) {
    degree[i] = 0;
  } 
  int currNode; 
  for (int j = 1; j < this->N; j++) {
    degree[parent[j]] += 1;
    degree[j] += 1; 
  }

  int numOddDegree = 0;
  for(int k = 0; k < this->N; k++) {
     if((degree[k]%2 == 1) || (degree[k] == 0)) {
       numOddDegree++;
     }
  }

 int oddDegrees[numOddDegree];
 int x = 0;
  for(int y = 0; y< this->N; y++) {
    if(degree[y]%2 == 1 || degree[y] == 0) {
      oddDegrees[x] = y;
      x++;
    }
  }

  for( int q = 0; q < numOddDegree; q++) {
    //cout << "odd node: " << oddDegrees[q] << endl;
  }
  return oddDegrees;

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
