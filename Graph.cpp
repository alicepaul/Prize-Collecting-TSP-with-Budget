//
//  Graph.cpp
//  Prize Collecting TSP
//
//  Created by Alice Paul on 3/14/17.
//  Copyright © 2017 Alice Paul. All rights reserved.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//

#include "Graph.hpp"

/* -------------------------EDGE--------------------------*/

// Create an edge
Edge::Edge(int h, int t, double w)
{
    //std::cout << "Creating Edge (" << h << "," << t << ") \n";
    head = h;
    tail = t;
    weight = w;
}

//Delete an edge
Edge::~Edge()
{
    //std::cout << "Deleting Edge (" << head << "," << tail << ") \n";
}

// Find the other end of an edge
int Edge::getOther(int id)
{
    if (head == id){
        return tail;
    }
    else{
        return head;
    }
}


// Comparison of two edges
bool operator== (const Edge &e1, const Edge &e2)
{
    return (((e1.getHead()==e2.getHead()) && (e1.getTail() == e2.getTail())) && (e1.getWeight() == e2.getWeight()));
}

// How to print an edge
std::ostream& operator<< (std::ostream &out, const Edge &edge)
{
    out << " ( " << edge.head << " , " << edge.tail << " , " << edge.weight << ") ";
    return out;
}

/* -------------------------VERTEX--------------------------*/

// Create a vertex (edge list is empty)
Vertex::Vertex(int i)
{
    //std::cout << "Making vertex " << i << "\n";
    degree = 0;
    id = i;
    prize = 1; // Default prize
}

// Create a vertex (edge list is empty)
Vertex::Vertex(int i, int p)
{
    //std::cout << "Making vertex " << i << "\n";
    degree = 0;
    id = i;
    prize = p;
}

// Delete a vertex
Vertex::~Vertex()
{
    //std::cout << "Destroying vertex " << id << "\n";
}

// Add an edge to the list of neighbors
void Vertex::addEdge(std::shared_ptr<Edge> e)
{
    neighbors.push_back(e);
    if (e->getWeight() < 0){
        degree -= e->getWeight();
    }
    else{
        degree += e->getWeight();
    }
}

// Remove an edge from the list of neighbors
void Vertex::removeEdge(std::shared_ptr<Edge> e)
{
    neighbors.remove(e);
    if (e->getWeight() < 0){
        degree += e->getWeight();
    }
    else{
        degree -= e->getWeight();
    }
}



/* -------------------------GRAPH--------------------------*/

Graph::Graph(){
    W = 0.0; // Nothing else to do
    P = 0.0;
}

Graph::~Graph(){
    // Nothing to do
}

// copy constructor
Graph::Graph(const Graph &G){
    W = G.W;
    P = G.P;
    std::list<int>::const_iterator it;
    for(it = G.vertices.begin(); it != G.vertices.end(); it++){
        vertices.push_back(*it);
        int i = *it;
        int p = G.getVertex(i)->getPrize();
        std::shared_ptr<Vertex> v = std::make_shared<Vertex>(i,p);
        vertexmap[*it] = v;
    }
    
    std::list<std::shared_ptr<Edge>>::const_iterator it2;
    for(it2 = G.edges.begin(); it2 != G.edges.end(); it2++){
        int id1 = (**it2).getHead();
        int id2 = (**it2).getTail();
        double w = (**it2).getWeight();
        std::shared_ptr<Edge> e = std::make_shared<Edge>(id1,id2,w);
        vertexmap[id1]->addEdge(e);
        vertexmap[id2]->addEdge(e);
        edges.push_back(e);
    }
}

// subgraph constructor
Graph::Graph(const Graph &G, const std::list<int> &S){
    W = 0;
    P = 0;
    // add vertices in S
    for (auto x : S){
        int p = G.getVertex(x)->getPrize();
        addVertex(x,p);
    }
    
    // add edges with both endpts in S
    for(auto e : G.getEdges()){
        int id1 = e->getHead(), id2 = e->getTail();
        double w = e->getWeight();
        
        if ((vertexmap.count(id1) > 0) && (vertexmap.count(id2) > 0)){
            addEdge(id1, id2, w);
        }
    }
}


// Add a vertex to a graph with a prize
void Graph::addVertex(int id, int p)
{
    vertices.push_back(id);
    std::shared_ptr<Vertex> v = std::make_shared<Vertex>(id,p);
    vertexmap[id] = v;
    P += p;
}

// Add a vertex to a graph
void Graph::addVertex(int id)
{
    addVertex(id,1);
}


// Add an edge to a graph
void Graph::addEdge(int id1, int id2, double weight)
{
    if ((vertexmap[id1]==NULL) || (vertexmap[id2]==NULL)){
        return;
    }
    std::shared_ptr<Edge> e = std::make_shared<Edge>(id1, id2, weight);
    vertexmap[id1]->addEdge(e);
    vertexmap[id2]->addEdge(e);
    edges.push_back(e);
    if (weight < 0){
        W -= weight;
    }
    else{
        W += weight;
    }
}



// Return ptr to vertex
std::shared_ptr<Vertex> const & Graph::getVertex(int id) const
{
    if (vertexmap.count(id) == 0){
        throw std::invalid_argument("Vertex does not exist");
    }
    return vertexmap.at(id);
}

// Return vertex degree
double Graph::getVertexDegree(int id) const
{
    if (vertexmap.count(id) == 0){
        throw std::invalid_argument("Vertex does not exist");
    }
    return vertexmap.at(id)->getDegree();
}


// Remove vertex from list of vertices and vertex map
void Graph::removeVertexLists(int id){
    vertices.remove(id);
    vertexmap.erase(id);
}

// Delete a vertex from a graph
void Graph::deleteVertex(int id)
{
    std::shared_ptr<Vertex> p_v = vertexmap[id];
    if (p_v == NULL){
        return;
    }
    
    // Update prize
    P -= p_v->getPrize();
    
    // If no edges to go through just remove references and delete vertex
    if (p_v->getIncEdges().size() == 0){
        removeVertexLists(id);
        return;
    }
    
    // Iterate through edges
    std::list<std::shared_ptr<Edge>>::iterator it;
    for (it = edges.begin(); it != edges.end();)
    {
        if (((*it)->getHead() == id) || ((*it)->getTail() == id)){
            // Delete other vertice's pointer to the edge and the lists's pointer
            if ((*it)->getWeight() < 0){
                W += (*it)->getWeight();
            }
            else{
                W -= (*it)->getWeight();
            }
            int j = (*it)->getOther(id);
            vertexmap[j]->removeEdge(*it);
            edges.erase(it++);
        }
        else{
            ++it;
        }
    }
    
    // Delete references in the map/list
    removeVertexLists(id);
}


// Print a graph
std::ostream& operator<< (std::ostream &out, const Graph &graph)
{
    // out << "Vertices: " << "\n";
    std::list<int>::const_iterator it;
    for (it = graph.vertices.begin(); it != graph.vertices.end(); ++it){
        out << *it << ", ";
    }
    out << "\n";
    //out << "Edges: " << "\n";
    //std::list<std::shared_ptr<Edge>>::const_iterator it2;
    //for (it2 = graph.edges.begin(); it2 != graph.edges.end(); ++it2){
    //    out << **it2;
    //}
    //out << "\n";
    return out;
}


// Minimum spanning tree
double Graph::MST(std::list<std::shared_ptr<Edge>> & edges) const
{
    std::map<int,bool> inTree;
    std::map<int,double> key;
    std::map<int,std::shared_ptr<Edge>> edge_keys;
    int numInTree = 0;
    double weightTree = 0;
    
    // initialize all vertices
    for (auto i : vertices){
        key[i] = INT_MAX, inTree[i] = false, edge_keys[i] = NULL;
    }
    
    // add first vertex
    key[vertices.front()] = 0;
    
    // while not all vertices are in the tree
    while (numInTree < vertices.size()){
        // find vertex with min key
        double min = INT_MAX;
        int min_v;
        
        for (auto v : vertices){
            if ((inTree[v] == false) && (key[v] < min)){
                min = key[v], min_v = v;
            }
        }
        
        // add min vertex to tree and edge
        inTree[min_v] = true;
        if (edge_keys[min_v] != NULL){
            edges.push_back(edge_keys[min_v]);
            weightTree += edge_keys[min_v]->getWeight();
        }
        
        // update key values
        for (auto e: getVertex(min_v)->getIncEdges()){
            int u = e->getOther(min_v);
            if (e->getWeight() < key[u]){
                key[u] = e->getWeight();
                edge_keys[u] = e;
            }
        }
        numInTree += 1;
        
    }
    
    return weightTree;
}


// Calculate tour length
double getTourLength(const Graph & G, const std::vector<int> & tour)
{
    int n = tour.size();
    double l = 0;
    for (int i = 0; i < tour.size()-1; i++){
        double w = shortestPath(G, tour[i], tour[i+1]);
        l += w;
        //std::cout << tour[i] << "," << tour[i+1] << "," << w << "\n";
    }
    l += shortestPath(G, tour[n-1], tour[0]);
    return l;
}

// DFS
void DFS(std::list<std::shared_ptr<Edge>> & edges, std::map<int,bool> & visited,
                   std::vector<int> & tour, int v)
{
    // Mark v as visited and add to tour
    visited[v] = true;
    tour.push_back(v);
    
    // Iterate through edges
    for (auto e:edges)
    {
        if ((e->getTail() == v) || (e->getHead() == v)){
            int u = e->getOther(v);
            if (visited[u] == false){
                DFS(edges, visited, tour, u);
            }
        }
    }
    
}

// Find ordered list of tour given a tree of edges
std::vector<int> tourList(std::list<std::shared_ptr<Edge>> & edges)
{
    // Find starting vertex and create bool map
    int v = edges.front()->getHead();
    std::map<int,bool> visited;
    for (auto e:edges){
        visited[e->getHead()] = false, visited[e->getTail()] = false;
    }
    std::vector<int> tour;
    
//    std::cout << "Constructing tour of length " << edges.size() +1 << "\n";
//    
//    for (auto e:edges){
//        std::cout << *e << ",";
//    }
//    std::cout << "\n";
    
    DFS(edges, visited, tour, v);
    return tour;
}

// Find length of shortest path from i to j in G
double shortestPath(const Graph & G, int i, int j)
{
    // Set up distance maps
    std::map<int,double> distances;
    std::map<int,bool> setDist;
    for (auto v:G.getVertices()){
        setDist[v] = false;
        distances[v] = INT_MAX;
    }
    
    // Start at i
    distances[i] = 0;
    
    // While j is not set
    while (setDist[j] == false){
        // Find min unmarked
        int minIndex = INT_MAX, u;
        for (auto v:G.getVertices()){
            if ((setDist[v] == false) && (distances[v] < minIndex)){
                u = v;
            }
        }
        
        // Set distances of u
        setDist[u] = true;
        
        // Iterate through incident edges to u
        if (distances[u] != INT_MAX){
            for (auto e : G.getVertex(u)->getIncEdges()){
                int v = e->getOther(u);
                double w = e->getWeight();
                if ((setDist[v] == false) && (distances[v] > distances[u]+w)){
                    distances[v] = distances[u]+w;
                }
            }
        }
    }
    return distances[j];
    
}
