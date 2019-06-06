int* getlabel_count();


// This class represents a directed graph using
// adjacency list representation
class Graph
{
    public:
    int V;    // No. of vertices
 
    // In a weighted graph, we need to store vertex
    // and weight pair for every edge
    list< pair<int, int> > *adj;
    list< pair<int, int> > *adj_reverse;
 
   
    Graph(int V);  // Constructor
 
    // function to add an edge to graph
    void addEdge(int u, int v, int w);
    void addEdge_toreverse(int u, int v, int w);
    int dijkstra_priorityq(int s, int v);
    int bidijkstra_priorityq(int s, int v);

};

void readFromFile1(string a);
     
Graph readFromFile2( string a, string b);

int nearestLabel_dijkstra (int v, int l, Graph g);

int nearestLabel_bidijkstra (int v, int l, Graph g);

void clean(int k, Graph g);
/*
void Graph::dijkstra_priorityq_all(int src, int level);
unordered_map<int,int> Graph::dijkstra_priorityq_cluster(int src, int level);
void Graph::dijkstra_priorityq_label(int src, int k, int * Ak, int nk);
*/
