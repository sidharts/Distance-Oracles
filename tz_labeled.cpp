// THORUP AND ZWICK'S DISTANCE ORACLE FOR VERTEX-LABELED GRAPHS
// priority_queue in STL
#include <bits/stdc++.h>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
using namespace std;

#include"tz_labeled.h"

#define MAXINT 99999
//# define INF 0x3f3f3f3f
#define FILE "graphList.txt"
#define LABEL_FILE "labels.txt"
 
// iPair ==>  Integer Pair
typedef pair<int, int> iPair;


// Number of vertices in the graph

void printSolution(int path[]);
int* size;
int* distances;
int * dist;
int ** distforall;
int ** P;
int ** DA;
int * label_size;
int * label;
int * label_count;
int * max_label_number;

int * getlabel_count(){
   return label_count;
}

 
// Allocates memory for adjacency list
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<iPair> [V];
}
 
void Graph::addEdge(int u, int v, int w)
{
    adj[u].push_back(make_pair(v, w));
}
 
// Prints shortest path from src to dest
int Graph::dijkstra_priorityq(int src, int dest){
    int u, v, weight;
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
 
    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    vector<int> dist(V, MAXINT);
 
    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(make_pair(0, src));
    dist[src] = 0;
 
    /* Looping till priority queue becomes empty (or all
      distances are not finalized) */
    while (!pq.empty()){
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        u = pq.top().second;
        pq.pop();

        if(u==dest){
           return dist[u];
        }
 
        // 'i' is used to get all adjacent vertices of a vertex
        list< pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i){
            // Get vertex label and weight of current adjacent
            // of u.
            v = (*i).first;
            weight = (*i).second;
 
            //  If there is shorter path to v through u.
            if (dist[v] > dist[u] + weight){
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
            }
        }
    }
 
    return dist[dest];
}

void Graph::dijkstra_priorityq_all(int src, int level){
    int u, v, weight;
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
 
    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    vector<int> dist(V+1, MAXINT);
 
    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(make_pair(0, src));
    dist[src] = 0;
 
    /* Looping till priority queue becomes empty (or all
      distances are not finalized) */
    while (!pq.empty()){
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        u = pq.top().second;
        pq.pop();
 
        // 'i' is used to get all adjacent vertices of a vertex
        list< pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i){
            // Get vertex label and weight of current adjacent
            // of u.
            v = (*i).first;
            weight = (*i).second;
            //  If there is shorter path to v through u.
            if (dist[v] > dist[u] + weight){
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
                if(u==src){
                   P[level][v] = v;
                }
                else{
                   P[level][v] = P[level][u];
                }
                DA[level][v] = dist[v];
            }
        }
    }
}

unordered_map<int,int> Graph::dijkstra_priorityq_cluster(int src, int level){
    int u, v, weight;
    unordered_map<int, int> singlecluster;
    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
 
    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    vector<int> dist(V+1, MAXINT);
 
    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(make_pair(0, src));
    dist[src] = 0;
 
    /* Looping till priority queue becomes empty (or all
      distances are not finalized) */
    while (!pq.empty()){
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        u = pq.top().second;
        pq.pop();
 
        // 'i' is used to get all adjacent vertices of a vertex
        list< pair<int, int> >::iterator i;
        for (i = adj[u].begin(); i != adj[u].end(); ++i){
            // Get vertex label and weight of current adjacent
            // of u.
            v = (*i).first;
            weight = (*i).second;
 
            //  If there is shorter path to v through u.
            if (dist[v] > dist[u] + weight && DA[level+1][v] > dist[u] + weight){
               // Updating distance of v
               dist[v] = dist[u] + weight;
               pq.push(make_pair(dist[v], v));
               singlecluster[v] = dist[v];
                
            }
        }
    }

    return singlecluster;
}

void readFromFile1(string labelfile){
   
    int l, n, i, j, count, max;
    ifstream file(labelfile);
    file>>n;
    label = new int[n];

    for (j = 0; j < n; j++) {
         file>>l;
         label[j]=l;
         if(l>max) max = l;
    }

    label_size = new int(max+1);
    label_count = new int[max+1];
    max_label_number = new int(max);
    for (j = 0; j <= max; j++) {
         label_count[j] = 0;
    }
    for (j = 0; j < n; j++) {
         label_count[label[j]]++;
    }

    file.close();

}

Graph readFromFile2(int k, string graphfile ,string labelfile) {

     readFromFile1(labelfile);

     int edges, a, b, wt, n, i;
     
     ifstream file2(graphfile);
     file2>>n;
     file2>>edges;
     size = new int(n);
     Graph g(n+k-1 + (*max_label_number+1) );
     
     while(edges-->0) {
        file2>>a;
        file2>>b;
        file2>>wt;
        g.addEdge(a-1, b-1, wt);
     }

     file2.close();
     P = new int*[k];
     for(i = 0; i < k; ++i)
        P[i] = new int[n];
     DA = new int*[k+1];
     for(i = 0; i < k+1; ++i)
        DA[i] = new int[n];

     return g;

}


int closest(Graph g, int * neighbors, int v, int n){
  int closest = -1, closest_dist = MAXINT+1, curr_neighbor, temp_dist;
  for(int i = 0; i < n; i++){
    curr_neighbor = neighbors[i];
    temp_dist = g.dijkstra_priorityq(curr_neighbor, v);
    if(temp_dist < closest_dist){
      closest_dist = temp_dist;
      closest = curr_neighbor;
    }
  }
  
  return closest;
}

void addSuperNodes(vector<int>* A, Graph g, int k){
   int i = 0, j= 0, V = *size, l = *max_label_number+1;
   vector <int> :: iterator it;
   for(i=1;i<k;i++){
      for(it = A[i].begin(); it!= A[i].end(); it++){
         g.addEdge(V+i-1, *it, 0);
      }
   }

   for (i = 0; i < l; ++i){
       for (j = 0; j < V; ++j){
          if (i == label[j]){
             g.addEdge(V+k-1+i, j, 0);
          }   
       }
    }
}

unordered_map<int,unordered_map<int, int> > preprocess_k(Graph g, int k){
  vector <int> A[k+1];
  vector <int> :: iterator it;
  vector <int> :: iterator it2;
    
  int N = *size;
  
  float sample_size = pow(N,1.0/k);
  int s = sample_size * 10000;

  
  for(int i = 0; i < N; i++ ){
    A[0].push_back(i);
  }
  
  while(A[k-1].empty()){
    for(int i = 1; i < k; i++){
      for(it = A[i-1].begin(); it != A[i-1].end(); it++){
          if(rand()%s<10000){
            A[i].push_back(*it);
          }
      }
    }
  }

  /*cout<<"\n=================Ais========================\n";
  for(int i = k-1; i >= 0; i--){
    for (it = A[i].begin(); it != A[i].end(); ++it)
         cout << *it << " ";
         cout<<"\n";
  }
  cout<<"\n=========================================\n";*/

  addSuperNodes(A, g, k);

  unordered_map<int,unordered_map<int, int> > BunchMap;
  unordered_map<int,unordered_map<int, int> > ClusterMap;
  
  bool exists;
  int temp_dist;

  for(int j = 1; j < k; j++){
     g.dijkstra_priorityq_all(*size+j-1, j);
  }

  for(int v = 0; v < N; v++){
    P[0][v] = v;
    DA[0][v] = 0;
    DA[k][v] = MAXINT+10;
  }
    
    
  for(int i = k-1; i >= 0; i --){
      for(it = A[i].begin(); it!= A[i].end(); it++){
          exists=0;
          for(it2 = A[i+1].begin(); it2!= A[i+1].end(); it2++){
             if(*it==*it2){
                exists=1;
                break;
             }
          }
          if(exists==0){
             ClusterMap[*it] = g.dijkstra_priorityq_cluster(*it, i);
          } 
      }
  }

  for (unordered_map<int,unordered_map<int, int> >::iterator it = ClusterMap.begin(); it != ClusterMap.end(); it++){
      //cout << it->first << ":";
      for (unordered_map<int, int>::iterator lit = it->second.begin(); lit  != it->second.end(); lit++){
          BunchMap[lit->first][it->first] = lit->second;
      }
  }
  /*
  cout<<"\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% BUNCHES OF VERTICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n";
  for (unordered_map<int,unordered_map<int, int> >::iterator it = BunchMap.begin(); it != BunchMap.end(); it++){
      cout << it->first << ":";
      for (unordered_map<int, int>::iterator lit = it->second.begin(); lit  != it->second.end(); lit++)
          cout << " " << lit->first;
      cout << "\n";
  }
  cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%";
  */
  return BunchMap;
}

int query(int u , int v, unordered_map<int,unordered_map<int, int> > BunchMap){
  int w = u,temp = 0, i = 0;
  unordered_map<int, int> singlebunch = BunchMap[v];
  
  while(singlebunch.find(w) == singlebunch.end() ){
    i++;
    temp = v;
    v = u;
    u = temp;
    w = P[i][u];
    singlebunch = BunchMap[v];
  }
  
  return DA[i][u] + singlebunch[w];
}

int nearestLabel (Graph g, int v, int l){
   if(label[v] == l){
      return 0;
   }
   int N = *size;
   int dis = MAXINT, temp_dist;
   for (int i = 0; i < N; ++i){
      if(label[i]==l){
         temp_dist = g.dijkstra_priorityq(v, i);
         if (dis > temp_dist){
           dis = temp_dist;
         }
      }
   }

   return dis;
}

int querylabel (int u , int l, unordered_map<int,unordered_map<int, int> > BunchMap){
  int flag = 0;
   if(label[u] == l){
      return 0;
   }
   int N = *size;
   int dis = MAXINT;
   int est = 0;
   for (int i = 0; i < N; ++i){
      if(label[i]==l){
         est = query(u, i, BunchMap);
         if (dis > est){
           dis = est;
         }
      }
   }

   return dis;
}



