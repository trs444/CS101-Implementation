#include "common.h"
#include "point.h"
#include <stack>
#include <set>
#include <list>
#include "Minmatching/PerfectMatching.h"
#include <vector>
#pragma once

class MST {
public:
	float** adjacentMatrix;
	int* parent; //Array to store constructed MST
	int* key; //Key values used to pick minimum weight edge in cut
	bool* mstSet; //To represent set of vertices not yet included in MST
	int N; //the size of pointset
	int* oddDegrees;
	bool* isOdd;
	std::stack<float> mystack;
	std::vector<pair<int,int>> edge;
	std::vector<pair<int,int>> cyclePath;
	std::vector<int> vertexPath;
	std::vector<int> finalVPath;
	bool pathflag = false;
	bool combineFlag = true;
	bool end = false;
	bool here = false;

	MST(float** adjacentMatrix, int size);
	~MST();

	//deliverable a
	void makeTree();
	void printMST();

	//deliverable b
	float makeTSP2(int startV, Point p);

	//deliverable c
	void makeTSP1_5(PerfectMatching* pm, int* oddNodes,Point p1);
	void findOddNode();
	int numOfOdd();
	int* whichOddVertices();
	void minimumMatching();
	//deliverable c
	//void makeTSP1_5(PerfectMatching* pm, int* oddNodes,Point p1);
	//int* findOddNodes();
	//int oddNodeDegrees();

	//float calMean(int option);
	//float calStd(int option);
	std::vector<pair<int,int>> findCycle(int start, bool visited[], std::vector<pair<int,int>> cycleEdge);
	void findPath(int begin, int target, bool visited[]);
	bool remEdge(int vertex, bool visited[]);
	void combineCycles(int start, std::list<int> negMark);
	std::list<int> sortNegMark(int start, std::list<int> tempMark);
	bool isValidEdge(int vertex1, int vertex2, bool visited[]);
	std::vector<int> shortcut(std::vector<int> vertexPath);
	
private:	
	float calMean(int option);
	float calStd(int option);
    
    
	vector<std::pair<int, int>> combine(PerfectMatching *pm, int* oddNodes);
	int minKey(int key[], bool mstSet[]);

};
