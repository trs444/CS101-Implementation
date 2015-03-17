#include "common.h"
#include "point.h"
#include "MST.h"
#include "Minmatching/PerfectMatching.h"
#include <iterator>
#include <vector>
#include <iostream>

/*
This project is a starter code and wrappers for CSE101W15 Implementation project.

point.h - uniform random pointset generator

MST.h - minimum spanning tree

PerfectMatching.h - interface to min cost perfect matching code 

-------------------------------------
PerfectMatching is from the paper:

Vladimir Kolmogorov. "Blossom V: A new implementation of a minimum cost perfect matching algorithm."
In Mathematical Programming Computation (MPC), July 2009, 1(1):43-67.

sourcecode : pub.ist.ac.at/~vnk/software/blossom5-v2.05.src.tar.gz

*/

void LoadInput(int& node_num, int& edge_num, int*& edges, int*& weights, float** adjacentMatrix, int numOddDegrees, int* oddArray) {
	int e = 0;
	node_num = numOddDegrees;
	edge_num = numOddDegrees*(numOddDegrees-1)/2 ; //complete graph

	edges = new int[2*edge_num];
	weights = new int[edge_num];

	for(int i = 0; i < numOddDegrees; ++i) {
		for(int j = i+1 ; j< numOddDegrees ; ++j) {
			edges[2*e] = i;
			edges[2*e+1] = j;
			weights[e] = adjacentMatrix[oddArray[i]][oddArray[j]];
			e++;
		}
	}

    /*for( int i = 0; i < edge_num; i++) {
    	cout << "weight: " << weights[i] << endl;
    }
    for ( int j = 0; j < edge_num*2; j++) {
    	cout << "edges: " << edges[j] << endl;
    }*/
	if (e != edge_num) { 
		cout<<"the number of edge is wrong"<<endl;

		exit(1); 
	}
}

void PrintMatching(int node_num, PerfectMatching* pm, int* oddArray) {
	int i, j;

	for (i = 0; i < node_num; i++) {
		j = pm->GetMatch(i);
		//cout << "i: " << oddArray[i] << endl;
		//cout << "j: " << oddArray[j] << endl;
		//if (i < j) printf("%d %d\n", oddArray[i], oddArray[j]);
	}
	/*for( int i = 0; i < node_num; i++) {
		cout << "oddArray: " << oddArray[i] << endl;
	}*/
}

int main() {
	set< pair<int,int> > generatedPointset;
	float** adjacentMatrix;
	int W, H, N;
	Point pointset;
	int numOddDegrees;
	int* oddArray;

	W = 100;
	H = 100;
	N = 8;

	cout<<"W: "<<W<<" H: "<<H<<" N:"<<N<<endl;

	pointset.generatePoint(W, H, N); //max(W,H,N) should be < 20000 because of memory limitation
	pointset.printPointset();

	generatedPointset = pointset.getPointset();
	adjacentMatrix = pointset.getAdjacentMatrix();

	//Deliverable A: From pointset and adjacentMatrix, you should construct MST with Prim or Kruskal
	MST mst(adjacentMatrix, N);
	mst.makeTree();
	mst.printMST();
	mst.findOddNode();
	numOddDegrees = mst.numOfOdd();

	oddArray = mst.whichOddVertices();
	//cout<<"NUM ODD DEGREE:" << numOddDegrees<<endl;

	//Deliverable B: Find TSP2 path from the constructed MST
	//You won't need any wrappers for B.
    float TSP2Cost = mst.makeTSP2(0, pointset);


	//Deliverable C: Find TSP1.5 path from the constructed MST
	
	//Find the perfect minimum-weight matching 
	struct PerfectMatching::Options options;
	int i, e, node_num = numOddDegrees, edge_num = numOddDegrees*(numOddDegrees-1)/2;
	int* edges;
	int* weights;
	PerfectMatching *pm = new PerfectMatching(node_num, edge_num);

	LoadInput(node_num, edge_num, edges, weights, adjacentMatrix, numOddDegrees, oddArray);

	for (e=0; e<edge_num; e++) {
		pm->AddEdge(edges[2*e], edges[2*e+1], weights[e]);
	}

	pm->options = options;
	pm->Solve();

	double cost = ComputePerfectMatchingCost(node_num, edge_num, edges, weights, pm);
	printf("Total cost of the perfect min-weight matching = %.1f\n", cost);
	
	PrintMatching(node_num, pm, oddArray);

	/*std::vector<pair<int, int>> nGraph;
	nGraph = mst.combine();
	std::pair<int, int> nPair;
	int q, t;

	for (q=0; q<node_num; q++) {
		t = pm->GetMatch(q);
		//if (i < j) printf("%d %d\n", oddArray[i], oddArray[j]);
		if(q < t) {
			nPair = std::make_pair(oddArray[q], oddArray[t]);
			//cout << "oddQ: " << oddArray[q] << endl;
			//cout << "oddT: " << oddArray[t] << endl;
			nGraph.push_back(nPair);
		}
	}
	for(std::vector<pair<int, int>>::iterator it = nGraph.begin(); it != nGraph.end(); ++it) {
		cout << "first value: " << it->first << endl;
		cout << "second value: " << it->second << endl;
	}*/
    //std::vector<pair<int, int>> nGraph;
    //nGraph = mst.combine(pm, oddArray);
    mst.makeTSP1_5(pm, oddArray, pointset);

	delete pm;
	delete [] edges;
	delete [] weights;
	
	return 0;
}
