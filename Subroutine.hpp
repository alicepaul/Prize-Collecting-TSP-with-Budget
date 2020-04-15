//
//  Subroutine.hpp
//  Prize-Collecting TSP
//
//  Created by Alice Paul on 3/14/17.
//  Copyright © 2017 Alice Paul. All rights reserved.
//

#ifndef Subroutine_hpp
#define Subroutine_hpp

#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <set>
#include <vector>
#include <cmath>
#include <climits>
#include "Graph.hpp"
#include "Subset.hpp"




/* ------------------------- HELPER FUNCTIONS --------------------------*/

// Returns whether or not an edge is active
// Endpoints must be in separate components and at least one component must be active
bool activeEdge(std::shared_ptr<Edge> & e, std::map<int,std::shared_ptr<Subset>> & vertexSubs);

// Returns 0.5 if both endpoints are active and 1.0 otherwise
double edgeFactor(std::shared_ptr<Edge> & e, std::map<int,std::shared_ptr<Subset>> & vertexSubs);


// Find min event to occur next
// Stores minimum subset pointer in min_s and minimum edge pointer in min_e
// Also finds event at lambda*(1+eps) (marked with _p)

void findMinEvent(double lambda, const std::list<std::shared_ptr<Subset>> & subsets,
                  std::vector<std::shared_ptr<Edge>> & E, std::map<int,std::shared_ptr<Subset>> & vertexSubs,
                  std::map<std::shared_ptr<Edge>,std::vector<double>> & lin_e,
                  std::map<std::shared_ptr<Subset>,std::vector<double>> & lin_s,
                  std::map<std::shared_ptr<Edge>,std::vector<double>> & lin_e_p,
                  std::map<std::shared_ptr<Subset>,std::vector<double>> & lin_s_p,
                  std::shared_ptr<Subset> & min_s, std::shared_ptr<Edge> & min_e,
                  std::vector<double> & lin_val,
                  std::shared_ptr<Edge> & alt_e,
                  std::vector<double> & lin_val_p1, std::vector<double> & lin_val_p2);


/* ------------------------- MAIN FUNCTIONS --------------------------*/

// Grow Function
// Runs the PD subroutine with lambda_1 = lambda and returns the end subsets

std::list<std::shared_ptr<Subset>> growSubsets(const Graph & G, double lambda);

// Prune Function
// Prune inactive components from s and add all non-pruned edges to edges
// If l_plus = true then use alternate edges and consider tied subsets as inactive

double prune(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
             bool l_plus = false, bool swap = true);


// Modification of the prune function which finds a path of maximal pruned subsets that are
// not pruned when s is marked active. Path is saved to prunedS and prunedE

void maxPrunedS(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
                std::shared_ptr<Subset> & tied_S, std::vector<std::shared_ptr<Subset>> & prunedS,
                std::vector<std::shared_ptr<Edge>> & prunedE);

// Modification of the prune function which finds a path of maximal pruned subsets that are
// not pruned when s uses alt edge. Path is saved to prunedS and prunedE

void maxPrunedE(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
                std::shared_ptr<Subset> & tied_S, std::vector<std::shared_ptr<Subset>> & prunedS,
                std::vector<std::shared_ptr<Edge>> & prunedE);

// Reverse Delete Function
// For each component in subsets, call prune to find remaining tight edges
// Returns the largest sum of tight edge weights and updates edges to contain those edges

double reverseDelete(std::list<std::shared_ptr<Subset>> & subsets, std::list<std::shared_ptr<Edge>> & edges,
                     std::shared_ptr<Subset> & largestSub, bool l_plus = false, bool swap = true);


/* ------------------------- ALTERNATE FUNCTIONS --------------------------*/

// To not keep track of edges
double reverseDelete(std::list<std::shared_ptr<Subset>> & subsets, bool l_plus = false, bool swap = true);
// Same as above but with a single component
double reverseDelete(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
                     bool l_plus = false, bool swap = true);
// Same as above with a single component but don't keep track of edges
double reverseDelete(std::shared_ptr<Subset> & s, bool l_plus = false, bool swap = true);


#endif /* Subroutine_hpp */
