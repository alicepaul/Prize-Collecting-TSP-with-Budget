//
//  Subsets.cpp
//  Prize-Collecting TSP
//
//  Created by Alice Paul on 3/13/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//

#include "Subset.hpp"


/* -------------------------SUBSETS--------------------------*/

// Constructors and Destructors

// Constructor without parents
Subset::Subset(int v, double pot, int pr)
{
    // Set parents to NULL and vertices to v
    vertices.push_back(v);
    parent1 = NULL;
    parent2 = NULL;
    
    potential = pot;
    y_val = 0;
    edge_total = 0;
    prize = pr;
    active = true;
    tied = false;
}

// Constructor with parents but no alt edge
Subset::Subset(const std::shared_ptr<Subset> p1, const std::shared_ptr<Subset> p2,
           const std::shared_ptr<Edge> e)
{
    vertices = p1->getVertices();
    vertices.insert(vertices.end(), p2->getVertices().begin(), p2->getVertices().end());
    parent1 = p1;
    parent2 = p2;
    edge = e;
    alt_edge = NULL;
    
    // Set values
    active = true;
    tied = false;
    y_val = 0;
    potential = p1->getPotential()+p2->getPotential()-2*(p1->getY()+p2->getY());
    edge_total = p1->getEdgeTotal()+p2->getEdgeTotal()+e->getWeight();
    prize = p1->getPrize() + p2->getPrize();
}


// Constructor with parents and an alt edge
Subset::Subset(const std::shared_ptr<Subset> p1, const std::shared_ptr<Subset> p2,
               const std::shared_ptr<Edge> e, const std::shared_ptr<Edge> alt_e)
{
    vertices = p1->getVertices();
    vertices.insert(vertices.end(), p2->getVertices().begin(), p2->getVertices().end());
    parent1 = p1;
    parent2 = p2;
    edge = e;
    alt_edge = alt_e;
    
    // Set values
    active = true;
    tied = false;
    y_val = 0;
    potential = p1->getPotential()+p2->getPotential()-2*(p1->getY()+p2->getY());
    edge_total = p1->getEdgeTotal()+p2->getEdgeTotal()+e->getWeight();
    prize = p1->getPrize() + p2->getPrize();
}

// Destructor
Subset::~Subset()
{
    // Will ensure shared ptrs are safely deleted
}

// Swap function
void Subset::swapEdges()
{
    std::shared_ptr<Edge> temp_e = alt_edge;
    alt_edge = edge;
    edge = temp_e;
    edge_total += edge->getWeight();
    edge_total -= alt_edge->getWeight();
}

// Print Function
std::ostream& operator<< (std::ostream &out, const Subset &S)
{
    out << "VERTICES: ";
    for (auto v: S.vertices){
        out << v << ", ";
    }
    out << "\n";
    //if (S.parent1 != NULL){
    //    out << "EDGE: " << *S.edge << "\n";
    //}
    //out << "Y: " << S.y_val << "\n";
    //out << "POTENTIAL: " << S.potential << "\n";
    //out << "ACTIVE: " << S.active << "\n";
    return out;
}


/* -------------------------FUNCTIONS--------------------------*/

// Find max potential subset among list of subsets and parents of these subsets
std::shared_ptr<Subset> findMaxPotential(const std::list<std::shared_ptr<Subset>> & subsets, int p){
    double max = -INT_MAX;
    std::shared_ptr<Subset> max_s = NULL;
    
    // Iterate through subsets to update max
    for (auto s:subsets){
        // If best potential seen and has prize at least p
        if ((s->getPotential() > max) && (s->getPrize() >= p)){
            max = s->getPotential();
            max_s = s;
        }
        
        // Recurse on parents
        if (s->getParent1() != NULL){
            
            // make list of parents
            std::list<std::shared_ptr<Subset>> parents;
            parents.push_back(s->getParent1());
            parents.push_back(s->getParent2());
            
            std::shared_ptr<Subset> test = findMaxPotential(parents);
            if ((test != NULL) && (test->getPotential() > max)){
                max = test->getPotential();
                max_s = test;
            }
        }
    }
    return max_s;
}

// Find all maximal laminar sets in subsets (inc ancestors) that has potential >= p
std::list<std::shared_ptr<Subset>> findHighPotential( const std::list<std::shared_ptr<Subset>> & subsets, double p)
{
    // List to return
    std::list<std::shared_ptr<Subset>> highS;
    
    // Iterate through subsets to add to list
    for (auto s:subsets){
        
        // If s has high enough potential then add and don't recurse on parents
        if (s->getPotential() > p+0.0001){
            highS.push_back(s);
        }
        
        // Otherwise recurse on parents (if non-null)
        else if (s->getParent1() != NULL){
            
            // make list of parents
            std::list<std::shared_ptr<Subset>> parents;
            parents.push_back(s->getParent1());
            parents.push_back(s->getParent2());
            
            // Recurse and add lists to highS
            std::list<std::shared_ptr<Subset>> test_list = findHighPotential(parents,p);
            highS.insert(highS.end(), test_list.begin(), test_list.end());
        }
    }
    return highS;
}

// Find  set in subset (inc ancestors) that contains edges with highest potential
std::shared_ptr<Subset> findMaxSuperset(std::shared_ptr<Subset> & s,
                                        const std::list<std::shared_ptr<Edge>> & edges)
{
    // Base case: if no parents return NULL
    if (s->getParent1() == NULL){
        return NULL;
    }
    
    // Base case: no edges
    if (edges.size() == 0){
        std::list<std::shared_ptr<Subset>> subsets;
        subsets.push_back(s);
        return findMaxPotential(subsets);
    }
        
    // Otherwise find parents and connecting edge(s)
    std::shared_ptr<Subset> p1 = s->getParent1(), p2 = s->getParent2();
    std::shared_ptr<Edge> e = s->getEdge(), alt_e = s->getAltEdge();
    
    // If e or alt_e in edges then return s
    for (auto test_e : edges){
        if ((test_e == e) || (test_e == alt_e)){
            return s;
        }
    }
    
    // If not then check both parents to see which contains s and has highest potential
    std::shared_ptr<Subset> s1 = findMaxSuperset(p1,edges);
    std::shared_ptr<Subset> s2 = findMaxSuperset(p2,edges);

    if ((s1 != NULL) && (s1->getPotential() > s->getPotential())){
        return s1;
    }
    else if ((s2 != NULL) && (s2->getPotential() > s->getPotential())){
        return s2;
    }
    else if ((s1 == NULL) && (s2 == NULL)){
        return NULL;
    }
    else{
        return s;
    }
}

// Pick routine which returns contiguous edges in s from vertex v with weight at most limit
double pick(std::shared_ptr<Subset> & s, std::set<int> & visited, double limit, int v, std::list<std::shared_ptr<Edge>> & edges, std::shared_ptr<Edge> & last_e)
{
    // If s has no parents then just return because no edges to add
    if (s->getParent1() == NULL){
        return 0;
    }
    
    // Otherwise check which parent has v and which parent has the head of the connecting edge
    int head = s->getEdge()->getHead();
    bool vinP1 = false, headinP1 = false;
    for (auto u: s->getParent1()->getVertices()){
        if (u == v)
            vinP1 = true;
        if (u == head)
            headinP1 = true;
    }
    
    // Set p1, p2 and connecting vertices
    std::shared_ptr<Subset> p1 = s->getParent1(), p2 = s->getParent2();
    int v1 = head, v2 = s->getEdge()->getTail();
    if (vinP1 == false)
        p1 = s->getParent2(), p2 = s->getParent1();
    if (headinP1 == false)
        v1 = s->getEdge()->getTail(), v2=head;
    
    // Get weights
    double w1 = p1->getEdgeTotal(), w2 = p2->getEdgeTotal(), w = s->getEdge()->getWeight();
    
    // If w1 above limit then just recurse on p1
    if (w1 > limit){
        return pick(p1, visited, limit, v, edges,last_e);
    }
    
    // If w1+w above limit then just recurse on p1
    else if (w1+w > limit){
        last_e = s->getEdge();
        std::shared_ptr<Edge> temp_e = NULL;
        return pick(p1, visited, limit, v, edges,temp_e);
    }
    
    // Otherwise recurse on both and add edge between p1 and p2
    else{
        // Add edge
        double total = w;
        edges.push_back(s->getEdge());
        visited.insert(v2);
        visited.insert(v1);
        
        // Recurse on p1 and add all edges
        std::shared_ptr<Edge> temp_e = NULL;
        total += pick(p1, visited, w1+1, v1, edges,temp_e);
        
        // Recurse on p2 and only add as many edges as possible
        total += pick(p2, visited, limit-w1-w, v2, edges,last_e);
        return total;
    }
}



