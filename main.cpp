//
//  main.cpp
//  Prize Collecting TSP
//
//  Created by Alice Paul on 3/14/17.

//  MIT License
//  Copyright (c) 2020 alicepaul


#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include "Graph.hpp"
#include "Subset.hpp"
#include "Subroutine.hpp"
#include "PD.hpp"
#include "ReadFile.hpp"




int main(int argc, char* argv[]){
    
    std::string file = "bier127.tsp";
    
    Graph G;
    double meanEdgeWeight;
    int numNodes;
    std::string name;
    // read the graph from the file and get some stats from it
    graphFromFile("Instances/"+file,G,name,meanEdgeWeight,numNodes);
    //std::cout<<"mean edge weight:"<<meanEdgeWeight<<"\n"<<"number of nodes:"<<numNodes<<std::endl;
    //std::cout<< "Graph: "<<G<<std::endl;
    

    double l = 0.0020055;
    double r = 0.0020057;
    std::list<std::shared_ptr<Subset>> subsetsL = growSubsets(G,l);
    std::list<std::shared_ptr<Subset>> subsetsR = growSubsets(G,r);
    double wplus = reverseDelete(subsetsL,false);
    double wminus = reverseDelete(subsetsR,false);
    

    std::list<std::shared_ptr<Edge>> mst;
    double mst_w = G.MST(mst);
    double D = 0.5*mst_w;
    std::cout << "D is " << D << "\n";
    
    std::list<std::shared_ptr<Edge>> edges;
    double upper = 0;
    int recursions = 0;
    double lambda;
    bool found;
    PD(G,D,edges,upper,recursions,lambda,found,true);

}
