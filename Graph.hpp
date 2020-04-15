//
//  Graph.hpp
//  Prize Collecting TSP
//
//  Created by Alice Paul on 3/14/17.
//  Copyright © 2017 Alice Paul. All rights reserved.
//

#ifndef Graph_hpp
#define Graph_hpp

#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <set>
#include <climits>
#include <vector>

/* -------------------------EDGE--------------------------*/


class Edge
{
private:
    int head;                               // head vertex id
    int tail;                               // tail vertex id
    double weight;                          // weight of edge
    
public:
    // Constructors and Destructors
    Edge(int h, int t, double w);
    ~Edge();
    
    // Get Functions
    int getHead() const {return head;}
    int getTail() const {return tail;}
    double getWeight() const {return weight;}
    void setWeight(double w) {weight = w;}
    int getOther(int id);
    
    // Print Function
    friend std::ostream& operator<< (std::ostream &out, const Edge &edge);
    
    // Comparison Function
    friend bool operator==(const Edge &e1, const Edge &e2);
};


/* -------------------------VERTEX--------------------------*/


class Vertex
{
private:
    int id;                                     // id
    int prize;                                  // prize value
    std::list<std::shared_ptr<Edge>> neighbors; // list of incident edge pointers
    double degree;                              // sum of weights for incident edges
    
public:
    // Constructors and Destructors
    Vertex(int i);
    Vertex(int i, int p);
    ~Vertex();
    
    // Get Functions
    int getId() const {return id;}
    int getPrize() const {return prize;}
    std::list<std::shared_ptr<Edge>> const & getIncEdges() const {return neighbors;}
    double getDegree() const {return degree;}
    
    // Add and Remove Edges
    void addEdge(std::shared_ptr<Edge> e);
    void removeEdge(std::shared_ptr<Edge> e);
};


/* -------------------------GRAPH--------------------------*/


class Graph
{
private:
    std::map<int, std::shared_ptr<Vertex>> vertexmap;   // map of vertex ids to vertex data structures
    std::list<int> vertices;                            // list of vertex ids
    std::list<std::shared_ptr<Edge>> edges;             // list of edge pointers
    double W;                                           // total weight of edges
    int P;                                              // sum of prizes of vertices
    
    
    // Used internally for removing vertices
    void removeVertexLists(int id);
    
public:
    // Constructors and Destructors
    Graph();
    ~Graph();
    Graph(const Graph &G);                              // copy constructor
    Graph(const Graph &G, const std::list<int> &S);            // subgraph G(S)
    
    // Get Functions
    double getWeight() const {return W;}
    int getPrize() const {return P;}
    std::list<int> const & getVertices() const {return vertices;}
    std::list<std::shared_ptr<Edge>> const & getEdges() const {return edges;}
    std::shared_ptr<Vertex> const & getVertex(int i) const;
    double getVertexDegree(int it) const;
    
    // Add and Remove Functions
    void addVertex(int id);
    void addVertex(int id, int prize);
    void addEdge(int id1, int id2, double weight);
    void deleteVertex(int id);
    
    // Print Function
    friend std::ostream& operator<< (std::ostream &out, const Graph &graph);
    
    // Minimum spanning tree
    // Returns weight and tree is saved to edges
    double MST(std::list<std::shared_ptr<Edge>> & edges) const;
    
};

// Calculate tour length
double getTourLength(const Graph & G, const std::vector<int> & tour);

// DFS
void DFS(std::list<std::shared_ptr<Edge>> & edges, std::map<int,bool> & visited,
         std::list<int> & tour, int v);

// Find ordered list of tour given a tree of edges
std::vector<int> tourList(std::list<std::shared_ptr<Edge>> & edges);

// Find length of shortest path from i to j in G
double shortestPath(const Graph & G, int i, int j);






#endif /* Graph_hpp */
