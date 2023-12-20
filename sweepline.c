#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <string>
#include <queue>
#include <vector>
using namespace std;

struct Point{
  double x, y, z; //coords of the point
  int fwdEdge, bwdEdge; //index of the forward and backward edges
  int index; //index of the edge starting at this point
  int type; //classifies how to deal with which lines to add to active set
  //fwdEdge and index should be the same

  //sorting by right end point
  friend bool operator<(const Point& l, const Point& r)
  {
    if(l.x < r.x){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator<=(const Point& l, const Point& r)
  {
    if(l.x <= r.x){
      return true;
    } else {
      return false;
    }
  }
  friend bool operator>(const Point& l, const Point& r)
  {
    if(l.x > r.x){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator==(const Point& l, const Point& r)
  {
    if(l.x == r.x){
      return true;
    } else {
      return false;
    }
  }
  
};

struct Edge{
  int AdjFwd, AdjBwd; //index of the edges that this edge is adjacent to - forward and backward
  int index; // the index of this edge within the edges array - use to get point
  int isOriented; // gives orientation: 1=left->right, 2 = right->left
  double xs, ys, zs; //x,y,z coords of start point
  double xe, ye, ze; //x,y,z coords of end point
};

struct CEdge{
  int thisI, edge, back, forward;
  friend bool operator<(const CEdge& l, const CEdge& r)
  {
    if(l.edge < r.edge){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator<=(const CEdge& l, const CEdge& r)
  {
    if(l.edge <= r.edge){
      return true;
    } else {
      return false;
    }
  }
  friend bool operator>(const CEdge& l, const CEdge& r)
  {
    if(l.edge > r.edge){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator==(const CEdge& l, const CEdge& r)
  {
    if(l.edge == r.edge){
      return true;
    } else {
      return false;
    }
  }
};

struct Crossing{
  int thisE, otherE;
  int crossingIndex;
  double uA, uB;
  int seenBefore; //0 if first time spotting it
  double isAbove; //states if this edge of the crossing is on top (1 is over)
  int isPositive; //0 if not positive, 1 if positive crossing

  friend bool operator<(const Crossing& l, const Crossing& r)
  {
    int lmin, rmin;
    if(l.thisE < l.otherE){
      lmin = l.thisE;
    } else {
      lmin = l.otherE;
    }
    if(r.thisE < r.otherE){
      rmin = r.thisE;
    } else {
      rmin = r.otherE;
    }
    if(lmin < rmin){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator<=(const Crossing& l, const Crossing& r)
  {
    int lmin, rmin;
    if(l.thisE < l.otherE){
      lmin = l.thisE;
    } else {
      lmin = l.otherE;
    }
    if(r.thisE < r.otherE){
      rmin = r.thisE;
    } else {
      rmin = r.otherE;
    }
    if(lmin <= rmin){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator>(const Crossing& l, const Crossing& r)
  {
    int lmin, rmin;
    if(l.thisE < l.otherE){
      lmin = l.thisE;
    } else {
      lmin = l.otherE;
    }
    if(r.thisE < r.otherE){
      rmin = r.thisE;
    } else {
      rmin = r.otherE;
    }
    if(lmin > rmin){
      return true;
    } else {
      return false;
    }
  }

  friend bool operator==(const Crossing& l, const Crossing& r)
  {
    int lmin, rmin;
    if(l.thisE < l.otherE){
      lmin = l.thisE;
    } else {
      lmin = l.otherE;
    }
    if(r.thisE < r.otherE){
      rmin = r.thisE;
    } else {
      rmin = r.otherE;
    }
    if(lmin == rmin){
      return true;
    } else {
      return false;
    }
  }
  
};

int compareCrossing (const void * a, const void * b)
{
  if ( *(Crossing*)a <  *(Crossing*)b ) return -1;
  if ( *(Crossing*)a == *(Crossing*)b ) return 0;
  if ( *(Crossing*)a >  *(Crossing*)b ) return 1;
  else return 0;
}

int comparePoint (const void * a, const void * b)
{
  if ( *(Point*)a <  *(Point*)b ) return -1;
  if ( *(Point*)a == *(Point*)b ) return 0;
  if ( *(Point*)a >  *(Point*)b ) return 1;
  else return 0;
}

int compareCEdge (const void * a, const void * b)
{
  if ( *(CEdge*)a <  *(CEdge*)b ) return -1;
  if ( *(CEdge*)a == *(CEdge*)b ) return 0;
  if ( *(CEdge*)a >  *(CEdge*)b ) return 1;
  else return 0;
}


class TinyQueue
{ public:
  struct Edge* edges;
  int size;

  TinyQueue(struct Edge e[1000]){
    this->edges = e;
    this->size = 0;
  }

  int getSize(){
    return size;
  }

  //return -1 if l1 start left of l2
  int compare(struct Edge e1, struct Edge e2){
    if(e1.xe < e2.xe){
      return 1;
    } else {
      return -1;
    }
  }
  
  bool push(struct Edge e){
    for(int i=0; i<size; i++){
      if(e.index == edges[i].index){
	return false;
      }
    }
			  
    edges[size] = e;
    size+=1;
    
    //bubble process
    //should bubble the smallest - the leftmost to the top
    int child = size-1;
    int parent = (child-1)/2;
    while(parent >= 0 && compare(edges[child], edges[parent]) > 0 && child != parent){
      struct Edge temp = edges[parent];
      edges[parent] = edges[child];
      edges[child] = temp;
      child = parent;
      parent = (child-1)/2;
    }
    return true;
  }
  
  bool isEmpty(){
    return size ==0;
  }
  
  struct Edge pop(){
    struct Edge e;
    if(isEmpty()){
      printf("EMPTY ARRAY - CANNOT POP FURTHER");
      return e;
    }
    struct Edge toReturn = edges[0];
    edges[0] = edges[size-1];
    edges[size-1] = e;
    size-=1;
    
    int parent = 0;
    while(true){
      int leftChild = 2*parent + 1;
      int rightChild = leftChild+1;
      if(leftChild >= size){
	break;
      }
      int maxChild = leftChild;
      if(rightChild < size && compare(edges[rightChild], edges[leftChild]) > 0) {
	maxChild = rightChild;
      }
      if(compare(edges[maxChild], edges[parent]) > 0){
	struct Edge temp = edges[parent];
	edges[parent] = edges[maxChild];
	edges[maxChild] = temp;
	parent = maxChild;
      } else {
	break;
      }
    }
    return toReturn;
  }

  struct Edge peek(){
    if(isEmpty()){
      printf("EMPTY ARRAY - CANNOT PEEK");
      struct Edge e;
      return e;
    }
    return edges[0];
  }

};



void merge(struct Point arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
 
    /* create temp arrays */
    struct Point L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
 
/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(struct Point arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
 
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
 
        merge(arr, l, m, r);
    }
}

void merge(struct CEdge arr[], int l, int m, int r)
{
    int i, j, k;
    int n1 = m - l + 1;
    int n2 = r - m;
 
    /* create temp arrays */
    struct CEdge L[n1], R[n2];
 
    /* Copy data to temp arrays L[] and R[] */
    for (i = 0; i < n1; i++)
        L[i] = arr[l + i];
    for (j = 0; j < n2; j++)
        R[j] = arr[m + 1 + j];
 
    /* Merge the temp arrays back into arr[l..r]*/
    i = 0; // Initial index of first subarray
    j = 0; // Initial index of second subarray
    k = l; // Initial index of merged subarray
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            arr[k] = L[i];
	    
            i++;
        }
        else {
            arr[k] = R[j];
            j++;
        }
        k++;
    }
 
    /* Copy the remaining elements of L[], if there
    are any */
    while (i < n1) {
        arr[k] = L[i];
        i++;
        k++;
    }
 
    /* Copy the remaining elements of R[], if there
    are any */
    while (j < n2) {
        arr[k] = R[j];
        j++;
        k++;
    }
}
 
/* l is for left index and r is right index of the
sub-array of arr to be sorted */
void mergeSort(struct CEdge arr[], int l, int r)
{
    if (l < r) {
        // Same as (l+r)/2, but avoids overflow for
        // large l and h
        int m = l + (r - l) / 2;
 
        // Sort first and second halves
        mergeSort(arr, l, m);
        mergeSort(arr, m + 1, r);
 
        merge(arr, l, m, r);
    }
}


bool testSegmentIntersect (struct Edge seg1,struct Edge seg2, vector<Crossing> &crossings) {
  //if (seg1 == NULL || seg2 == NULL){ return false;}

        double x1 = seg1.xs;
        double y1 = seg1.ys;
        double x2 = seg1.xe;
        double y2 = seg1.ye;
        double x3 = seg2.xs;
        double y3 = seg2.ys;
        double x4 = seg2.xe;
        double y4 = seg2.ye;
	if((x1==x3&&y1==y3) || (x1==x4&&y1==y4) || (x2==x3&&y2==y3) || (x2==x4&&y2==y4)){
	  return false;
	}

        double denom = ((y4 - y3) * (x2 - x1)) - ((x4 - x3) * (y2 - y1));
        double numeA = ((x4 - x3) * (y1 - y3)) - ((y4 - y3) * (x1 - x3));
        double numeB = ((x2 - x1) * (y1 - y3)) - ((y2 - y1) * (x1 - x3));

        if (denom == 0) {
	  if (numeA == 0 && numeB ==0){ return false;}
	  return false;
        }

        double uA = numeA / denom;
        double uB = numeB / denom;
	double positive;

        if (uA >= 0 && uA <= 1 && uB >= 0 && uB <= 1) {
            double x = x1 + (uA * (x2 - x1));
            double y = y1 + (uA * (y2 - y1));
	    //printf("%f, %f\n", x, y);
	    //modified
	    struct Crossing c1;
	    struct Crossing c2;
	    c1.thisE = seg1.index;
	    c1.otherE = seg2.index;
	    c2.thisE = seg2.index;
	    c2.otherE = seg1.index;
	    c1.uA = uA;
	    c1.uB = uB;
	    c2.uA = uB;
	    c2.uB = uA;
	    c1.seenBefore = 0;
	    c2.seenBefore = 0;
	    //calculate which is over and positive
	    double over;
	    over = (uA*seg1.ze+(1.0-uA)*seg1.zs) - (uB*seg2.ze+(1.0-uB)*seg2.zs);
	    if(over>0.0){
	      //seg1 is over
	      c1.isAbove = 1;
	      c2.isAbove = 0;
	      if(seg1.isOriented == 1 && seg2.isOriented == 1){
		//seg1 = l->r, seg2 = l->r
		positive = (seg1.xe-seg1.xs)*(seg2.ye-seg2.ys) - (seg1.ye-seg1.ys)*(seg2.xe-seg2.xs);
	      } else if(seg1.isOriented == 1 && seg2.isOriented == 2){
		//seg1 = l->r, seg2 = r->l
		positive = (seg1.xe-seg1.xs)*(seg2.ys-seg2.ye) - (seg1.ye-seg1.ys)*(seg2.xs-seg2.xe);		
	      } else if(seg1.isOriented == 2 && seg2.isOriented == 1){
		//seg1 = r->l, seg2 = l->r
		positive = (seg1.xs-seg1.xe)*(seg2.ye-seg2.ys) - (seg1.ys-seg1.ye)*(seg2.xe-seg2.xs);
	      } else {
		//seg1 = r->l, seg2 = r->l
		positive = (seg1.xs-seg1.xe)*(seg2.ys-seg2.ye) - (seg1.ys-seg1.ye)*(seg2.xs-seg2.xe);
	      }
	      if(positive<0.01){
		c1.isPositive = 0;
		c2.isPositive = 0;
	      } else {
		c1.isPositive = 1;
		c2.isPositive = 1;
	      }
	    } else {
	      //seg2 is over
	      c1.isAbove = 0;
	      c2.isAbove = 1;

	      if(seg1.isOriented == 1 && seg2.isOriented == 1){
		//seg1 = l->r, seg2 = l->r
		positive = (seg2.xe-seg2.xs)*(seg1.ye-seg1.ys) - (seg2.ye-seg2.ys)*(seg1.xe-seg1.xs);
	      } else if(seg2.isOriented == 1 && seg1.isOriented == 2){
		//seg1 = r->l, seg2 = l->r
		positive = (seg2.xe-seg2.xs)*(seg1.ys-seg1.ye) - (seg2.ye-seg2.ys)*(seg1.xs-seg1.xe);
	      } else if(seg2.isOriented == 2 && seg1.isOriented == 1){
		//seg1 = l->r, seg2 = r->l
		positive = (seg2.xs-seg2.xe)*(seg1.ye-seg1.ys) - (seg2.ys-seg2.ye)*(seg1.xe-seg1.xs);
	      } else {
		//seg1 = r->l, seg2 = r->l
		positive = (seg2.xs-seg2.xe)*(seg1.ys-seg1.ye) - (seg2.ys-seg2.ye)*(seg1.xs-seg1.xe);
	      }
	      if(positive<0.1){
		c1.isPositive = 0;
		c2.isPositive = 0;
	      } else {
		c1.isPositive = 1;
		c2.isPositive = 1;
	      }
	    }
	    
	    crossings.push_back(c1);
	    crossings.push_back(c2);
	    return true;
        }
        return false;
}



void runCheck(TinyQueue active, struct Edge e, vector<Crossing> &crossings){
  for(int i=0; i<active.getSize(); i++){
    if(active.edges[i].ys < active.edges[i].ye){
      int belowLeft = e.ys < active.edges[i].ys && e.ye < active.edges[i].ys;
      int aboveRight = e.ys  > active.edges[i].ye && e.ye > active.edges[i].ye;
      if(belowLeft || aboveRight){
	continue;
      }
    } else {
      int belowRight = e.ys  < active.edges[i].ye && e.ye < active.edges[i].ye;
      int aboveLeft = e.ys > active.edges[i].ys && e.ye > active.edges[i].ys;
      if(belowRight || aboveLeft){
	continue;
      }
      }

    //do the check for intersection - e with all of the active set
    if(e.index < active.edges[i].index){
      bool found = testSegmentIntersect(e, active.edges[i], crossings);
    } else {
      bool found = testSegmentIntersect(active.edges[i], e, crossings);
    }
  }
}


string crossingCode(string code[][4], struct Crossing thisC, struct Crossing otherC){
  //determine the code/letter to go along with the 0/1/2...a
  char thisLetter, otherLetter;
  // printf("Positive: %d, Above: %f\n", c.isPositive, c.isAbove);
  if(thisC.isAbove == 1){
    thisLetter = 'a';
  } else {
    //this crossing is under
    if(thisC.isPositive == 1){
      thisLetter = 'd';
    } else {
      thisLetter = 'b';
    }
  }

  if(otherC.isAbove == 1){
    otherLetter = 'c';
  } else{
    if(otherC.isPositive == 1){
      otherLetter = 'b';
    } else {
      otherLetter = 'd';
    }
  }
  string otherCode = to_string(otherC.crossingIndex) + otherLetter;
  string thisCode = to_string(thisC.crossingIndex) + thisLetter;
  if(thisLetter == 'a'){
    code[thisC.crossingIndex][0] = otherCode.c_str();
  } else if (thisLetter == 'b'){
    code[thisC.crossingIndex][1] = otherCode;
  } else {
    code[thisC.crossingIndex][3] = otherCode;
  }
  
  if(otherLetter == 'b'){
    code[otherC.crossingIndex][1] = thisCode;
  } else if (otherLetter == 'c'){
    code[otherC.crossingIndex][2] = thisCode;
  } else {
    code[otherC.crossingIndex][3] = thisCode;
  }
  return otherCode;
}


int main(){
  double x, y, z;
  //START FILE READING
  char line[200]; // create line var
  cin.getline(line, 500); // gets first line
  if(strcmp(line, "VECT")!=0){ // compares first line to "VECT" to ensure correct file
    printf("Wrong file type - enter a vect file");
    exit(1);
  }
  //gets 2nd line which contains line number
  cin >> line; // get 1st part - # of components
  if(line[0]!='1'){ // compares within second line to make sure that only 1 component
    printf("Not 1 component");
    exit(1);
  }
  int nr_lines;
  cin >> nr_lines; //get nr of lines - 2nd part of 2nd line
  cin.getline(line, 500); // throw rest of line2 away
  cin.getline(line, 500); //get line 3
  if(line[0] == '-'){ //check for closed component + the correct nr of vertices
  }
  cin.getline(line, 500); // skip the colors per component

  struct Point points[nr_lines+1];
  struct Edge edges[nr_lines+1];
  //Deal with 1st point as well here
  cin >> x >> y >> z;
  
  for(int i=0; i < nr_lines; i++){ // dont go to last line - colour info
    if(cin.peek()=='#'){
      i = i-1; //decrement to not count this line
      cin.getline(line, 500);//get line + throw it away
      continue;
    }
    
    //now put all the coordinates into the struct
    points[i].x = x;
    points[i].y = y;
    points[i].z = z;
    points[i].index = i;
    points[i].fwdEdge = i;
    points[i].bwdEdge = i-1;
    if(i==0){
      points[i].bwdEdge = nr_lines-1;
    }

    //ADD TO EDGES ARRAY
    //need to get the next point to add to edges
    cin >> x >> y >> z;
    
    edges[i].index = i;
    edges[i].AdjFwd = i+1;
    edges[i].AdjBwd = i-1;
    points[i].index = i;
    
    if(points[i].x < x){
      // point at i left of i+1 = l->r
      if(edges[i-1].isOriented == 2 || i == 0){
	//if previous edge is r->l 
	points[i].type = 0;
      } else {
	//l->r
	points[i].type = 1;
      }
      edges[i].xs = points[i].x;
      edges[i].ys = points[i].y;
      edges[i].zs = points[i].z;
      edges[i].xe = x;
      edges[i].ye = y;
      edges[i].ze = z;
      edges[i].isOriented = 1;
    } else {
      //r->l
      if(edges[i-1].isOriented == 1){
	//if previous edge is l->r
	points[i].type = 3;
      } else {
	points[i].type = 2;
      }
      edges[i].xs = x;
      edges[i].ys = y;
      edges[i].zs = z;
      edges[i].xe = points[i].x;
      edges[i].ye = points[i].y;
      edges[i].ze = points[i].z;
      edges[i].isOriented = 2;
    }
    cin.getline(line, 500); //get rid of everything at the end of the line
  }

  //put first coord at end as well
  points[nr_lines] = points[0];
  points[nr_lines].index = 0;
  points[nr_lines].bwdEdge = 0;
  edges[nr_lines-1].index = nr_lines-1;

  if(edges[nr_lines-1].xs < edges[nr_lines-1].xe){
    //l->r
    edges[nr_lines-1].xs = points[nr_lines-1].x;
    edges[nr_lines-1].ys = points[nr_lines-1].y;
    edges[nr_lines-1].zs = points[nr_lines-1].z;
    edges[nr_lines-1].xe = points[nr_lines].x;
    edges[nr_lines-1].ye = points[nr_lines].y;
    edges[nr_lines-1].ze = points[nr_lines].z;
    edges[nr_lines-1].isOriented = 1;
  } else {
    edges[nr_lines-1].xs = points[nr_lines].x;
    edges[nr_lines-1].ys = points[nr_lines].y;
    edges[nr_lines-1].zs = points[nr_lines].z;
    edges[nr_lines-1].xe = points[nr_lines-1].x;
    edges[nr_lines-1].ye = points[nr_lines-1].y;
    edges[nr_lines-1].ze = points[nr_lines-1].z;
    edges[nr_lines-1].isOriented = 2;
  }
  for(int i=0; i<nr_lines; i++){
    printf("%d ", points[i].index);
  }
  //                                                               SOMEWHERE HERE PUT A ROTATION IF THE Z VALUES ARE TOO CLOSE TOGETHER
  //mergeSort(points, 0, nr_lines-1);
  qsort(&points, nr_lines+1, sizeof(Point), comparePoint);
  
  // SORTED POINTS ARRAY FROM LEFT TO RIGHT

  
  
  //HERE THE INTERSECTION
  struct Edge a[1000];
  TinyQueue active = TinyQueue(a);
  vector<Crossing> crossings;
  for(int i=0; i<nr_lines+1;i++){
    //look at which type to figure out to add/remove
    int type = points[i].type;
    //active.push(edges[points[i].index]);
    if(type == 0){
      runCheck(active, edges[points[i].index], crossings);
      active.push(edges[points[i].index]);
      runCheck(active, edges[points[i].bwdEdge], crossings);
      active.push(edges[points[i].bwdEdge]);
      for(int j=0; j < active.size; j++){
	printf("0:%d,", active.edges[j].index);
      }
    } else if(type == 1){
      runCheck(active, edges[points[i].index], crossings);
      active.push(edges[points[i].index]);
      active.pop();
      for(int j=0; j < active.size; j++){
	printf("1:%d,", active.edges[j].index);
      }
      
    } else if(type == 2){
      runCheck(active, edges[points[i].bwdEdge], crossings);
      active.push(edges[points[i].bwdEdge]);
      active.pop();
      for(int j=0; j < active.size; j++){
	printf("2:%d,", active.edges[j].index);
      }
    } else {
      
      active.pop();
      active.pop();
      for(int j=0; j < active.size; j++){
	printf("%d,", active.edges[j].index);
      }
    }
   
    printf("\n");
  }
  int counter = 0;
  CEdge es[crossings.size()];
  //printf("sizeof: %lu", sizeof(Crossing));
  //qsort(&crossings[0], crossings.size(), sizeof(Crossing), compareCrossing);

  for(int i=0; i<crossings.size(); i++){
    printf("crossing:%d, this:%d\nOther:%d, Positive:%d\n\n", crossings[i].crossingIndex, crossings[i].thisE, crossings[i].otherE, crossings[i].isPositive);
    //printf("edge: %d, index: %d\n", es[i].edge, es[i].thisI);
    }
  
  for(int i=0; i<crossings.size(); i++){
    if(crossings[i].seenBefore == 0){
      crossings[i].crossingIndex = counter;
      crossings[i].seenBefore = 1;
      crossings[i+1].seenBefore = 1;
      crossings[i+1].crossingIndex = counter;
      counter++;
      es[i].thisI = i;
      es[i].edge = crossings[i].thisE;
      es[i+1].thisI = i+1;
      es[i+1].edge = crossings[i].otherE;
    }
  }

  qsort(&es, crossings.size(), sizeof(CEdge), compareCEdge);
  
  string code[crossings.size()/2][4];
  string thisCode, otherCode;
  int ci;
  for(int i = 0; i <crossings.size(); i++){
    //printf("Loop %d\n",i);
    // go forward one and look at that crossing to determine the code
    ci = crossings[es[i].thisI].crossingIndex; //this crossing index
    thisCode = to_string(ci)+"a";
    
    if(i+1 == crossings.size()){ //at end of array so 'forward' is the first index
      otherCode = crossingCode(code, crossings[es[i].thisI], crossings[es[0].thisI]);
    }
    //need to look at the edge forward of this one to determine code
    else {
      otherCode = crossingCode(code, crossings[es[i].thisI], crossings[es[i+1].thisI]);
    }
  }

  for(int i=0; i<crossings.size()/2; i++){
    printf("%d", i);
    if(crossings[i*2+1].isPositive == 1){
      printf("+");
    } else {
      printf("-");
    }
    for(int j=0; j<4; j++){
      printf("%s", code[i][j].c_str());
    }
    printf("\n");
  }
  
  return 0;
}
