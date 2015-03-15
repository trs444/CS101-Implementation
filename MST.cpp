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
  //cout << "NODES VISITED IN THIS ORDER:" << endl;
  for(int k = 0; k < (sizeof(visited)/sizeof(int)); k++) {
    //cout << visited[k] << endl;
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

bool compare_case (const int first, const int second)
{
  if(first<second) 
    return true;
  else if(first>second)
    return false;
}
std::vector<pair<int,int>> MST::findCycle(int start, bool visited[], std::vector<pair<int,int>> cycleEdge){
  int parent;
  std::vector<pair<int,int>> travelEdge = cycleEdge;
  
  //loops through the edge to find a cycle
  for(int i = 0; i< edge.size(); i++){
    if((edge[i].first == start||edge[i].second == start) && visited[i] == false){
      //cout<<"THERE IS A EDGE THAT'S NOT VISITED"<<endl;
      visited[i] = true;
      parent = edge[i].second;
      cyclePath.push_back(edge[i]);
      vertexPath.push_back(edge[i].first);
      //cout<<"3pushing in: "<<edge[i].first <<" , "<< edge[i].second<<endl;
      pathflag = false;
      findPath(parent, start, visited);
      //cout<<"ay"<<endl;
      break;
      
    }
  }

  
 //cout<<"ALKSFJALDK"<<endl;
  // for(int s = 0; s<vertexPath.size(); s++){
  //  cout<<"VERTEXPATH "<< s <<" is: "<< vertexPath[s]<<endl;
  // }

  // for(int y = 0; y<travelEdge.size(); y++){
  //  cout<<"TRAVELEDGE at "<<y<<" is: "<<travelEdge[y].first <<" , "<< travelEdge[y].second<<endl;
  // }

  // for(int x = 0; x<edge.size(); x++){
  //  cout<<"visited array index "<< x <<" is: "<< visited[x]<<endl;
  // }

  for(int s = 0; s<vertexPath.size(); s++){

    //check if edge is removeable
    if(remEdge(vertexPath[s], visited)){
      //cout<<"removing vertex "<<endl;
      //cout<<"cycle size: "<<cycleEdge.size()<<endl;
      for(int n = 0; n<travelEdge.size(); n++){
        //cout<<"Travel edge at "<< n << " is: "<< travelEdge[n].first<<", "<< travelEdge[n].second<<endl;
        //cout<<"COMPARING VERTEXPATH "<<vertexPath[s]<<endl;
        if(vertexPath[s] == travelEdge[n].first || vertexPath[s] == travelEdge[n].second){
          //cout<<"yellll"<<endl;
          if(isValidEdge(travelEdge[n].first, travelEdge[n].second, visited)){
            //cout<<"removing edge: "<<travelEdge[n].first <<" "<< travelEdge[n].second<<endl;
            travelEdge.erase(travelEdge.begin()+ n);
            break;
            
          }
        }
      }
      //cout<<"what"<<endl;
    }
  }


  //cout<<"hai"<<endl;
  // for(int y = 0; y<travelEdge.size(); y++){
  //  cout<<"AFTER REMOVAL OF EDGES travelEdge at "<<y<<" is: "<<travelEdge[y].first <<" , "<< travelEdge[y].second<<endl;
  // }

  return travelEdge;

}


void MST::findPath(int begin, int target, bool visited[]){
  for(int i = 0; i< edge.size(); i++){
    if(edge[i].first == begin && visited[i] == false){
      visited[i] = true;
      cyclePath.push_back(edge[i]);
      vertexPath.push_back(edge[i].first);
      //cout<<"1pushing in: "<<edge[i].first <<" , "<< edge[i].second<<endl;
      
      //cout<<"TARGET: "<<target<<endl;

      if(edge[i].second == target){
        //cout<<"pushing target1"<<endl;
        vertexPath.push_back(edge[i].second);
        pathflag = true;
        return;
      }
      else{
        findPath(edge[i].second, target, visited);
        if(!pathflag){
          //cout<<"removing1: "<<vertexPath.back()<<endl;
          cyclePath.pop_back();
          vertexPath.pop_back();
          visited[i] = false;
        }

        //cout<<"out of recur1"<<endl;
      }
      
    }
    else if (edge[i].second == begin && visited[i] == false){
      visited[i] = true;
      cyclePath.push_back(edge[i]);
      vertexPath.push_back(edge[i].second);
      //cout<<"2pushing in: "<<edge[i].first <<" , "<< edge[i].second<<endl;
      
      if(edge[i].first == target){
        //cout<<"pushing target2"<<endl;
        vertexPath.push_back(edge[i].first);
        pathflag = true;
        return;
      }
      else{
        findPath(edge[i].first, target, visited);
        if(!pathflag){
          //cout<<"removing2: "<<vertexPath.back()<<endl;
          cyclePath.pop_back();
          vertexPath.pop_back();
          visited[i] = false;
        }
        //cout<<"out of recur2"<<endl;
      }
      
    }
    if(pathflag == true)
      break;
  }
        
  return;

}

bool MST::remEdge(int vertex, bool visited[]){
  bool flag1 = false;

    for(int i = 0; i< edge.size(); i++){
      if(edge[i].first == vertex || edge[i].second == vertex){
        if(visited[i] == true){
          flag1 = true;
          break;
        }
        else
          flag1 = false;
      }
    }

  return flag1;
}

void MST::combineCycles(int start, std::list<int> negMark){

  /*for (std::list<int>::iterator it=negMark.begin(); it!=negMark.end() ; ++it)
      cout << "combineSORTED NEG MARK IS: " << *it<<endl;*/
  
  std::list<int> negMark1 = negMark;
  int firstNeg = negMark1.front();
  
  int secondNeg;  
  // cout<<"firstNeg: "<<firstNeg<<endl;
  // cout<<"secondNeg: "<<secondNeg<<endl;

  // cout<<"COMBINE: "<<combineFlag<<endl;
  if(!negMark1.empty()){
    negMark1.pop_front();

    if(!negMark1.empty()){

      secondNeg = negMark1.front();
      //cout<<"inside secondNeg: "<<secondNeg<<endl;
      
    }
    else
      secondNeg = -1;

    //cout<<"NEW secondNeg: "<<secondNeg<<endl;
    //cout<<"NEGATIVE MARK NOT EMPTY"<<endl;

    for(int b = firstNeg+1; b<vertexPath.size(); b++){
      for(int a = firstNeg+1; a<vertexPath.size(); a++){
        if(vertexPath[a] != -1){
          //cout<<"comboflage set to false"<<endl;
          combineFlag = false;
        }
      }

      //cout<<"COMBINEFLAGE: "<<combineFlag<<endl;

      //cout<<"second neg: "<< secondNeg<< " vertex: "<<vertexPath[secondNeg+1]<<endl;
      //cout<<"vertex b: "<< vertexPath[b]<<endl;
      
      
        
      
      
      if (vertexPath[b] == vertexPath[secondNeg+1]){
        
        //cout<<"RECURRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR"<<endl;
        //combineFlag = true;
        if(negMark1.size()>1){
          
          negMark1.pop_front();

        }

        negMark = sortNegMark(secondNeg+1, negMark1);
        for (std::list<int>::iterator it=negMark.begin(); it!=negMark.end() ; ++it)
            cout << "second SORTED NEG MARK IS: " << *it<<endl;

        combineCycles(secondNeg+1, negMark);

        // negMark.pop_front();
        // cout<<"1ALSO PATH PUSHING: "<<vertexPath[b]<<endl;
        // finalVPath.push_back(vertexPath[b]);

        
      }
      

      else if(vertexPath[b] != -1){
        //cout<<"1VERTEX RECURRRRRRRR PATH PUSHING: "<<vertexPath[b]<<endl;
        finalVPath.push_back(vertexPath[b]);
      }
      else if(vertexPath[b] == -1){
        if(negMark1.size()>1){
          b = negMark1.front();
          negMark1.pop_front();
          //secondNeg = negMark1.front();

        }
        else
          break;
      }
      //cout<<"welp"<<endl;
    }
    //cout<<"end of for loop"<<endl;
  }

  if(negMark.empty() || secondNeg == -1){
    //cout<<"NEGATIVE MARK IS EMPTTTTTTTTTTY"<<endl;
    combineFlag = true;
    for(int b = start; b<vertexPath.size(); b++){
      if (vertexPath[b] == vertexPath[secondNeg+1]){
        //negMark.pop_front();
        cout<<"2ALSO PATH PUSHING: "<<vertexPath[b]<<endl;
        finalVPath.push_back(vertexPath[b]);
        //combineFlag = true;
        //combineCycle(secondNeg+1, negMark);
      }
      if(vertexPath[b] == -1){
        break;
      }

      
      // n
    }
  }

}

std::list<int>MST::sortNegMark(int start, std::list<int> tempMark){
  std::list<int> copy = tempMark;

  //first value = index on path 
  //second value = index of neg 1
  std::vector<pair<int,int>> s1;
  std::list<int> s2;
  std::vector<int> v2;
  std::list<int> finalMark1;
  std::list<int> finalMark;
  int size = tempMark.front();
  bool flag = false;
  int index = 0;
  int vertSize = vertexPath.size();

  // std::cout << "copy contains:";
 //   for (std::list<int>::iterator it=copy.begin(); it!=copy.end() ; ++it)
 //     std::cout << ' ' << *it<<endl;
  
  while(!copy.empty()){

    // cout<<"STARRRRRRRRRRRRRRRRRRRRRRRT:"<<start<<endl;
    // cout<<"STARRRRRRRRRRRRRRRRRRRRRRRT:"<<copy.front()<<endl;
    //cout<<"aaaacopy is: "<<start<<" "<<vertexPath[copy.front()]<<endl;
    for(index = start; index<vertexPath.size(); index++){
      //cout<<"copy is: "<<start<<" "<<vertexPath[copy.front()+1]<<endl;
      // cout<<"vertex"<<index<<": "<<vertexPath[index]<<endl;
      //cout<<"reached"<<endl;

      if(vertexPath[index] == -1)
        break;
      if(copy.empty()){
          //cout<<"in hereeeeeeeee"<<endl;
          flag = true;
          break;
      }
      
      else if(vertexPath[index] == vertexPath[copy.front()+1]){
        s1.push_back(make_pair(index, copy.front()));
        //cout<<"s2 pushing in: "<<index<<" and "<<copy.front()+1<<endl;
        s2.push_back(index);
        //cout<<"what"<<copy.size()<<endl;
        if(!copy.empty()){
          copy.pop_front();
          index = start;
          
        }
        // else
        //  break;
      }
      
    }
    
    //cout<<"noooo"<<endl;
    if(flag){
      //cout<<"noooo"<<endl;
      while(!copy.empty()){
        index++;
        //cout<<"COPY FRONT IS "<< copy.front()+1<<endl;
        s1.push_back(make_pair(index, copy.front()));
        s2.push_back(index);
        copy.pop_front();
        
      }
      break;
    }
    else if(!copy.empty()){
      vertSize++;
      //cout<<"not in cycle pushing in: "<<vertSize<<" and "<<copy.front()+1<<endl;
      s1.push_back(make_pair(vertSize, copy.front()));
      s2.push_back(vertSize);
      copy.pop_front();
    }
  }

  // std::cout << "AFTER copy contains:";
 //   for (std::list<int>::iterator it=copy.begin(); it!=copy.end() ; ++it)
 //     std::cout << ' ' << *it<<endl;

  // std::cout << "s2 contains:";
 //   for (std::list<int>::iterator it=s2.begin(); it!=s2.end() ; ++it)
 //     std::cout << ' ' << *it<<endl;
 

  s2.sort(compare_case);
  std::list<int> sorteds2 = s2;

  while(!sorteds2.empty()){
    //cout<<"pushing v2: "<<sorteds2.front()<<endl;
    v2.push_back(sorteds2.front());
    sorteds2.pop_front();
  }

  // std::cout << "s2 contains:";
 //   for (std::list<int>::iterator it=s2.begin(); it!=s2.end() ; ++it)
 //     std::cout << ' ' << *it<<endl;

  for(int s2start = 0; s2start<tempMark.size(); s2start++){
    for(int markstart = 0; markstart<tempMark.size(); markstart++){
      // cout<<"v2 at "<<s2start<<" : "<<v2[s2start]<<endl;
      // cout<<"s1 at "<<markstart<<" : "<<s1[markstart].first<<endl;
      if(v2[s2start] == s1[markstart].first){
        //cout<<"pushing in final mark: "<<s1[markstart].second<<endl;
        finalMark1.push_back(s1[markstart].second);
        //break;
      }
    }
  }

  //finalMark1.push_back(-1);
  // while(!finalMark1.empty()){
  //  finalMark.push_back(finalMark1.back());
  //  finalMark1.pop_back();
  // }

  return finalMark1;

}

bool MST::isValidEdge(int vertex1, int vertex2, bool visited[]){
  bool flag1;
  bool flag2;

  for(int i = 0; i<edge.size(); i++){
    if(edge[i].first == vertex1 || edge[i].second == vertex1){
      if(visited[i] == false){
        flag1 = true;
        break;
      }
      else
        flag1 = false;
    }
  }


  for(int k = 0; k<edge.size(); k++){
    if(edge[k].first == vertex2 || edge[k].second == vertex2){
      if(visited[k] == false){
        flag2 = true;
        break;
      }
      else
        flag2 = false;
    }
  }

  if(flag1 && flag2)
    return false;
  else
    return true;
  
}

std::vector<int> MST::shortcut(std::vector<int> startPath){
  std::vector<int> revPath;
  std::vector<int> finishedPath;
  std::vector<int> vPath = startPath;
  bool flag = false;
  //cout<<"readched"<<endl;
  //vPath.pop_back();

  while(!startPath.empty()){
    //cout<<"BEING PUSHING IN: "<<revPath.back()<<endl;
    vPath.push_back(startPath.back());
    startPath.pop_back();
  }

  while(!vPath.empty()){
    //cout<<"in while"<<endl;
    //cout<<"size "<< revPath.size() -1 <<endl;
    for(int i = revPath.size() -1; i> 0; i--){
      //cout<<"in for"<<endl;
      //cout<<"rev path is "<< revPath[i]<<endl;
      //cout<<"vpath path is "<< vPath.back() <<endl;
      if(vPath.back() == revPath[i]){
        //cout<<"MARRRRRRRRRRRRRRRRRRRY POOOOOOOOOOOOOPPPINGS"<<endl;
        vPath.pop_back();
        //cout<<"uh oh"<<endl;
        flag = true;
        break;
      }
      //cout<<"bye"<<endl;

    }

    //cout<<"am i done?"<<endl;

    if(flag == true){
      flag = false;
    }
    else if(flag == false){
      //cout<<"pushing in "<<vPath.back()<<endl;
      revPath.push_back(vPath.back());
      // if(!vPath.empty()){
      //  cout<<"POOOOOOOOOOOOOOOOOOOOOOOOOOP"<<endl;
      //  vPath.pop_back();
      //  cout<<"made it"<<endl;
      // }
    }
    //cout<<"what what"<<endl;

  }
  revPath.push_back(0);
  
  while(!revPath.empty()){
    //cout<<"BEING PUSHING IN: "<<revPath.back()<<endl;
    finishedPath.push_back(revPath.back());
    revPath.pop_back();
  }
  


  return finishedPath;

}

void MST::makeTSP1_5(PerfectMatching* pm, int* oddNodes, Point p1) {
  std::vector<pair<int,int>> remainingEdge;

  int start = 0;
  bool canExit= false;

  //store index of -1 to seperate the cycles
  std::list<int> negMark;
  std::list<int> tempMark;


  //make all edges have even degrees by combining mimimum-weight matching and MST edges
  edge = combine(pm, oddNodes);

  bool *visited = new bool[edge.size()];

  //set all the edges to false
  for(int z = 0; z < edge.size(); z++){
    visited[z] = false;
    //cout<<"visited is false"<<endl;
  }

  std::vector<pair<int,int>> cycleEdge = edge;
  
  remainingEdge = findCycle(start, visited, cycleEdge);
  // for(int w = 0; w<remainingEdge.size(); w++){
  //    cout<<"WELLL THIS IS Edge "<< w <<" is: "<< remainingEdge[w].first<<", "<<remainingEdge[w].second<<endl;
  //  }
  //cout<<"here"<<endl;

  while(!remainingEdge.empty()){
    canExit = false;
    
    //cout<<"WHILE LOOP CYCLE"<<endl;
    
    //cout<<"REMAINING EDGE SIZE IS: "<< remainingEdge.size()<<endl;
    //cout<<"VERTEXPATH SIZE IS: "<< vertexPath.size()<<endl;
    //check to see if start vertex is in path
    for(int rev = 0; rev<remainingEdge.size(); rev++){
      for(int d = 0; d<vertexPath.size(); d++){
        // cout<<"FIRST VERTEX OF REMAINING EDGE AT "<<rev<<" is: "<<remainingEdge[rev].first<<endl;
        // cout<<"SECOND VERTEX OF REMAINING EDGE AT "<<rev<<" is: "<<remainingEdge[rev].second<<endl;
        // cout<<"VERTEX PATH AT "<<d<<" is: "<<vertexPath[d]<<endl;
        if(remainingEdge[rev].first == vertexPath[d]){
          start = remainingEdge[rev].first;
          canExit = true;
          //cout<<"using a different start"<<endl;
          break;
        }
        else if(remainingEdge[rev].second == vertexPath[d]){
          start = remainingEdge[rev].second;
          canExit = true;
          //cout<<"using a different start"<<endl;
          break;
        }
      }

      if(canExit)
        break;
    }

    //cout<<"THE START OF IT ALL IS: "<<start<<endl;
    vertexPath.push_back(-1);
    //cout<<"222here"<<endl;

    remainingEdge = findCycle(start, visited, remainingEdge);

    //cout<<"22here"<<endl;
    // for(int w = 0; w<remainingEdge.size(); w++){
    //  cout<<"THE REMAINING EDGE"<< w <<" is: "<< remainingEdge[w].first<<", "<<remainingEdge[w].second<<endl;
    // }

  }

  //cout<<"vertex path size: "<<vertexPath.size()<<endl;

  for(int v = 0; v<vertexPath.size(); v++){
    if(vertexPath[v] == -1){
      //cout<<"UNSORTED NEG MARK: "<< v <<endl;
      tempMark.push_back(v);
    }
  }

  //cout<<"333here"<<endl;
  //sort the negative cycles
  //negMark.sort(compare_case);
  if(!tempMark.empty()){

    negMark = sortNegMark(0,tempMark);
    here = true;
  }
  else 
    negMark = tempMark;

  
    
  for (std::list<int>::iterator it=negMark.begin(); it!=negMark.end() ; ++it)
      //cout << "hereereeeeeeeeeeeeeeeeSORTED NEG MARK IS: " << *it<<endl;



  
//cout<<"222here"<<endl;


  if(!negMark.empty()){
    int nextV = negMark.front();
    int sizeV = negMark.front();
    int tears = 0;
    int startV = 0;

    //cout<<"FIRST NEGATIVE FOUND IS: "<<nextV<<endl;
    //combining the final vertex path
    for(startV = 0; startV<sizeV; startV++){
      tears++;
      //cout<<"startV at "<<startV<<" is: "<< vertexPath[startV]<<endl;
      //cout<<"nextV at "<<nextV<<" is: "<< vertexPath[nextV+1]<<endl;
      if(vertexPath[startV] == vertexPath[nextV+1] && !negMark.empty()){
        //cout<<"REMOVING"<<nextV<<endl;
        if(!here){
          negMark = sortNegMark(nextV+1, negMark);
          for (std::list<int>::iterator it=negMark.begin(); it!=negMark.end() ; ++it)
              cout << "bbafasdlkSORTED NEG MARK IS: " << *it<<endl;
          here = false;
        }
        combineCycles(nextV+1, negMark);

        
        //cout<<"OLDDDDDDD nextV at "<<nextV<<" is: "<< vertexPath[nextV+1]<<endl;
        // if(!negMark.empty()){
        //  negMark.pop_front();
        //  if(!negMark.empty())
        //    nextV = negMark.front();
        //  else
        //    break;
        // }
        //cout<<"NEWWWWWWW nextV at "<<nextV<<" is: "<< vertexPath[nextV+1]<<endl;

      }
      
      else if(vertexPath[startV] != -1){
        //cout<<"VERTEX PATH PUSHING: "<<vertexPath[startV]<<endl;
        finalVPath.push_back(vertexPath[startV]);
      }
      else if(vertexPath[startV] == -1){
        startV = negMark.front();
        negMark.pop_front();
        nextV = negMark.front();
        //cout<<"here startV at "<<startV<<" is: "<< vertexPath[startV]<<endl;
        //cout<<"here nextV at "<<nextV<<" is: "<< vertexPath[nextV+1]<<endl;
        //break;
      }
    }
      // cout<<"tears: "<<tears<<endl;
      // cout<<"sizev: "<<sizeV<<endl;
  }


  // for(int h = 0; h<finalVPath.size(); h++){
  //  cout<<"BEEEEEEFORE final v path is: "<< finalVPath[h]<<endl;
  // }

  for(int h = 0; h<vertexPath.size(); h++){
    //cout<<"vertexPath is at "<<h<<" : "<< vertexPath[h]<<endl;
  }

  // for(int x = 0; x<edge.size(); x++){
  //  cout<<"visited array index "<< x <<" is: "<< visited[x]<<endl;
  // }
  if(finalVPath.empty()){

    finalVPath = vertexPath;
  }
  // cout<<"blah"<<endl;
  for(int h = 0; h<finalVPath.size(); h++){
    //cout<<"final v path is: "<< finalVPath[h]<<endl;
  }

  std::vector<int> shortcutPath = shortcut(finalVPath);

  // for(int h = 0; h<shortcutPath.size(); h++){
  //  cout<<"SHORTCUT path is: "<< shortcutPath[h]<<endl;
  // }


  //cout<<"here"<<endl;
  //ADDING THE TSP COST


  std::vector<pair<int,int>> setvector;
  std::set< pair<int,int> > myset = p1.getPointset();

  for(set< pair<int,int> >::iterator it = myset.begin() ; it != myset.end() ; ++it) {
     //cout<<(*it).first<<" , "<<(*it).second<<endl;
     std::pair <int,int> temppair  = make_pair((*it).first,(*it).second);
     setvector.push_back(temppair);  
  }


  int x1;
    int x2;
    int y1;
    int y2;
    int var1;
    int var2;
    float total;
    float tempDist;

    var1 = shortcutPath.back();
    shortcutPath.pop_back();

    //p1.printPointset();

    while(!shortcutPath.empty()){
      var2 = shortcutPath.back();
      shortcutPath.pop_back();

      //cout<<"edge "<<var1<<" - "<< var2<<endl;
      
                total = total + adjacentMatrix[var1][var2];
                //cout<<"adding existing dist:"<<adjacentMatrix[var1][var2]<<endl;
        

      var1 = var2;

    }



    cout<<"Total cost of heuristic TSP1_5: "<< total<<endl;



  //calculate heuristic TSP cost 
}

void MST::minimumMatching() { //if you choose O(n^2)
	//find minimum-weight matching for the MST. 
	
	//you should carefully choose a matching algorithm to optimize the TSP cost.
}

std::vector<pair<int, int>> MST::combine(PerfectMatching *pm, int* oddNodes) {
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

  int q, t;
  int oddNum = numOfOdd();
  for (q=0; q<oddNum; q++) {
    t = pm->GetMatch(q);
    //if (i < j) printf("%d %d\n", oddArray[i], oddArray[j]);
    if(q < t) {
      //printf("There's an Edge from %d to %d \n", oddNodes[q], oddNodes[t]);
      std::pair<int, int> nPair = std::make_pair(oddNodes[q], oddNodes[t]);
      //cout << "oddQ: " << oddArray[q] << endl;
      //cout << "oddT: " << oddArray[t] << endl;
      result.push_back(nPair);
    }
  }
  for(vector<pair<int,int>>::iterator it = result.begin(); it!=result.end();++it){
    //cout<<"Result pairs are: " << (*it).first << " - " << (*it).second<<endl;
    //cout<<"WEIGHT IS: "<<adjacentMatrix[(*it).first][(*it).second]<<endl;
  }

  return result;

}
