#include "common.h"
#include "point.h"

#pragma once

class MST {
public:
	float** adjacentMatrix;
	int* parent; //Array to store constructed MST
	int* key; //Key values used to pick minimum weight edge in cut
	bool* mstSet; //To represent set of vertices not yet included in MST
	int N; //the size of pointset

	MST(float** adjacentMatrix, int size);
	~MST();

	//deliverable a
	void makeTree();
	void printMST();

	//deliverable b
	float makeTSP2(int startV, Point p);

	//deliverable c
	void makeTSP1_5();
	
	float calMean(int option);
	float calStd(int option);
    int findOddVertices();
    int* whichOddVertices();

private:
	void minimumMatching();
	void combine();
	int minKey(int key[], bool mstSet[]);

};
