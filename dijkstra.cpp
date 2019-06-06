// Program to find Dijkstra's shortest path using
// priority_queue in STL
#include <bits/stdc++.h>
#include <iostream>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */
using namespace std;
# define INF 99999//0x3f3f3f3f
#include"dijkstra.h"
 
// iPair ==>  Integer Pair
typedef pair<int, int> iPair;

int* size;
int* label;
int* label_size;
int * label_count;
int * max_label_number;

int* getlabel_count(){
   return &label_count[0];
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

Graph readFromFile2(string graphfile ,string labelfile) {

     readFromFile1(labelfile);

     int edges, a, b, wt, n, i;
     
     ifstream file2(graphfile);
     file2>>n;
     file2>>edges;
     size = new int(n);
     Graph g(n+(*max_label_number+1) );
     
     while(edges-->0) {
        file2>>a;
        file2>>b;
        file2>>wt;
        g.addEdge(a-1, b-1, wt);
        g.addEdge_toreverse(b-1, a-1, wt);
     }

     // Adding supernodes
    for (int i = 0; i < (*max_label_number+1); ++i){
       for (int j = 0; j < n; ++j){
          if (i == label[j]){
             g.addEdge(j, n+i, 0);
             g.addEdge_toreverse(n+i, j, 0);
          }   
       }
    }

     file2.close();

     return g;

}

 
// Allocates memory for adjacency list
Graph::Graph(int V)
{
    this->V = V;
    adj = new list<iPair> [V];
    adj_reverse = new list<iPair> [V];
}
 
void Graph::addEdge(int u, int v, int w)
{
    adj[u].push_back(make_pair(v, w));
}

void Graph::addEdge_toreverse(int u, int v, int w)
{
    adj_reverse[u].push_back(make_pair(v, w));
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
    vector<int> dist(V, INF);
 
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
 
            //  If there is shorted path to v through u.
            if (dist[v] > dist[u] + weight){
                // Updating distance of v
                dist[v] = dist[u] + weight;
                pq.push(make_pair(dist[v], v));
            }
        }
    }
    return dist[dest];
}

int Graph::bidijkstra_priorityq(int src, int dest)
{
    int u1, u2, v, weight, forward, backward, move, currDist;
    int minFinal = INF;

    bool sptSet1[1050];
    bool sptSet2[1050];

    for (int i = 0; i < V; ++i){
       sptSet1[i] = false;
       sptSet2[i] = false;
    }

    // 'i' is used to get all adjacent vertices of a vertex
    list< pair<int, int> >::iterator i;
    list< pair<int, int> >::iterator j;

    // Create a priority queue to store vertices that
    // are being preprocessed. This is weird syntax in C++.
    // Refer below link for details of this syntax
    // https://www.geeksforgeeks.org/implement-min-heap-using-stl/
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq;
    priority_queue< iPair, vector <iPair> , greater<iPair> > pq_reverse;
 
    // Create a vector for distances and initialize all
    // distances as infinite (INF)
    vector<int> dist1(V, INF);
    vector<int> dist2(V, INF);
 
    // Insert source itself in priority queue and initialize
    // its distance as 0.
    pq.push(make_pair(0, src));
    pq_reverse.push(make_pair(0, dest));
    dist1[src] = 0;
    dist2[dest] = 0;

    /* Looping till priority queue becomes empty (or all
      distances are not finalized) */
    while (!pq.empty() && !pq_reverse.empty()){
        // The first vertex in pair is the minimum distance
        // vertex, extract it from priority queue.
        // vertex label is stored in second of pair (it
        // has to be done this way to keep the vertices
        // sorted distance (distance must be first item
        // in pair)
        
        forward = pq.top().first;
        
        backward = pq_reverse.top().first;
        
        u1 = pq.top().second;

        u2 = pq_reverse.top().second;

        if(forward < backward){
            move = 1;
            //u = pq.top().second;
            if (sptSet2[u1]){
               sptSet1[u1] = true;
               break;
            }
        }
        else if(forward > backward){
            move = 2;
            //u = pq_reverse.top().second;
            if (sptSet1[u2]){
               sptSet2[u2] = true;
               break;
            }
        }
        else{
            if(u1 == u2)
               move = 3;
            else{
               //u = pq.top().second;
               if (sptSet2[u1] == true){
                  sptSet1[u1] = true;
                  move = 1;
                  break;
               }
               //u = pq_reverse.top().second;
               if (sptSet1[u2] == true){
                  sptSet2[u2] = true;
                  move = 2;
                  break;
               }
               move = 1;
            }
        }

        if (move == 1){
            //u = pq.top().second;
            pq.pop();
            if(u1==dest){
               return dist1[u1];
            }
            sptSet1[u1] = true;
        }

        if (move == 2){
            //u = pq_reverse.top().second;
            pq_reverse.pop();
            if(u2==src){
               return dist2[u2];
            }
            sptSet2[u2] = true;
        }

        if(move == 3){
          //u = pq.top().second;
          //In this case, u1 = u2
          if (u1 == dest){
             return dist1[u1];
          }
          if (u2 == src){
             return dist2[u2];
          }
          sptSet1[u1] = true;
          sptSet2[u2] = true;
          break;
        }
        
        if (move == 1){
           for (i = adj[u1].begin(); i != adj[u1].end(); ++i){
              // Get vertex label and weight of current adjacent
              // of u.
              v = (*i).first;
              weight = (*i).second;
 
              //  If there is shorted path to v through u.
              if (dist1[v] > dist1[u1] + weight){
                 // Updating distance of v
                 dist1[v] = dist1[u1] + weight;
                 pq.push(make_pair(dist1[v], v));
              }
           }
        }
        else if (move == 2){
           for (i = adj_reverse[u2].begin(); i != adj_reverse[u2].end(); ++i){
              // Get vertex label and weight of current adjacent
              // of u.
              v = (*i).first;
              weight = (*i).second;
 
              //  If there is shorted path to v through u.
              if (dist2[v] > dist2[u2] + weight){
                 // Updating distance of v
                 dist2[v] = dist2[u2] + weight;
                 pq_reverse.push(make_pair(dist2[v], v));
              }
           }
        }
        
    }
    
    if(move == 2) u1 = u2;
    
    minFinal = dist1[u1]+dist2[u1];
    for (u1 = 0; u1 < V; ++u1){
        if (sptSet1[u1]){
           for (j = adj[u1].begin(); j != adj[u1].end(); ++j){
              v = (*j).first;
              weight = (*j).second;
              if (sptSet2[v]){
                 currDist = dist1[u1]+weight+dist2[v];
                 if (currDist < minFinal){
                    minFinal = currDist;
                 }
              }
           }
        }
    }

     if (minFinal>INF || minFinal < 0){
        minFinal = INF;
     }
    return minFinal;
}
 
int nearestLabel_dijkstra (int v, int l, Graph g){
   int V = *size;
   //return v;
   return g.dijkstra_priorityq(v, V + l);
}

int nearestLabel_bidijkstra (int v, int l, Graph g){
   int V = *size;
   return g.bidijkstra_priorityq(v, V + l);
}