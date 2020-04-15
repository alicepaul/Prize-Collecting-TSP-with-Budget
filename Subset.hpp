//
//  Subsets.hpp
//  Prize-Collecting TSP
//
//  Created by Alice Paul on 3/13/17.
//  Copyright © 2017 Alice Paul. All rights reserved.
//

#ifndef Subsets_hpp
#define Subsets_hpp

#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <set>
#include <climits>
#include "Graph.hpp"


/* -------------------------SUBSETS--------------------------*/

class Subset
{
private:
    std::list<int> vertices;                        // vertices in set
    std::shared_ptr<Subset> parent1;                // parent subset 1
    std::shared_ptr<Subset> parent2;                // parent subset 2
    std::shared_ptr<Edge> edge;                     // tight edge between parents
    std::shared_ptr<Edge> alt_edge;                     // edge tight at lambda+eps
    
    double potential;                               // potential of S
    double y_val;                                   // y_S dual value
    double edge_total;                              // sum of weights of edges
    int prize;                                      // prize of vertices
    bool active;                                    // whether marked active in PD alg
    bool tied;                                      // mark as tied with an edge going tight
    
public:
    // Constructors and Destructors
    Subset(int v, double pot, int pr = 1);
    Subset(const std::shared_ptr<Subset> p1, const std::shared_ptr<Subset> p2,
           const std::shared_ptr<Edge> e);          // Use for merges
    Subset(const std::shared_ptr<Subset> p1, const std::shared_ptr<Subset> p2,
           const std::shared_ptr<Edge> e, const std::shared_ptr<Edge> alt_e);
                                                    // Use for merges with alt edges
    ~Subset();
    
    // Get Functions
    std::list<int> const & getVertices() const {return vertices;}
    std::shared_ptr<Subset> const & getParent1() const {return parent1;}
    std::shared_ptr<Subset> const & getParent2() const {return parent2;}
    std::shared_ptr<Edge> const & getEdge() const {return edge;}
    std::shared_ptr<Edge> const & getAltEdge() const {return alt_edge;}
    double getPotential() const {return potential;}
    double getY() const {return y_val;}
    double getEdgeTotal() const {return edge_total;}
    int getPrize() const {return prize;}
    bool getActive() const {return active;}
    bool getTied() const {return tied;}
    
    // Set Functions
    void swapEdges();
    void setY(double y) {y_val = y;}
    void setActive(bool b) {active = b;}
    void setTied(bool b) {tied = b;}
    void setPotential(double p) {potential = p;}
    
    // Print Function
    friend std::ostream& operator<< (std::ostream &out, const Subset &S);

};


/* -------------------------FUNCTIONS--------------------------*/

// Find subset with max potential among list of subsets and ancestors of those subsets
std::shared_ptr<Subset> findMaxPotential(const std::list<std::shared_ptr<Subset>> & subsets, int prize = 0);

// Find all maximal laminar sets among list of subsets and ancestors of those subsets that have potential > p
std::list<std::shared_ptr<Subset>> findHighPotential( const std::list<std::shared_ptr<Subset>> & subsets, double p);

// Find the set of maximum potential among list of subsets and ancestors of those subessts
// that contains tree
// Note: the tree is only made up of edges or alt_edges used when merging subsets
std::shared_ptr<Subset> findMaxSuperset(std::shared_ptr<Subset> & s,
                                        const std::list<std::shared_ptr<Edge>> & edges);

// Runs the pick routine which returns contiguous edges in s from vertex v with weight at most limit
// Keep tracks of spanned vertices in visited and edges used in edges, returns weight of edges added
double pick(std::shared_ptr<Subset> & s, std::set<int> & visited, double limit, int v,
            std::list<std::shared_ptr<Edge>> & edges, std::shared_ptr<Edge> & last_e);











#endif /* Subsets_hpp */
