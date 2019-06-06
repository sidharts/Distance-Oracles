// This class represents a directed graph using
// adjacency list representation
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
 
    int dijkstra_priorityq(int s, int v);
    void dijkstra_priorityq_all(int s, int level);
    unordered_map<int,int> dijkstra_priorityq_cluster(int s, int level);
};

void readFromFile1(string a);
     
Graph readFromFile2(int k, string a, string b);

unordered_map<int,unordered_map<int, int> > preprocess_k(Graph g, int k);

void preprocess_record_k(Graph g, int k);

int querylabel (int u , int l, unordered_map<int,unordered_map<int, int> > BunchMap);

//int getV();

//int getL();

void clean(int k, Graph g);
/*
void Graph::dijkstra_priorityq_all(int src, int level);
unordered_unordered_map<int,int> Graph::dijkstra_priorityq_cluster(int src, int level);
void Graph::dijkstra_priorityq_label(int src, int k, int * Ak, int nk);
*/
