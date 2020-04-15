//
//  ReadFile.cpp
//  
//
//  Created by Alice Paul on 4/10/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//
//
//

#include "ReadFile.hpp"

void tokenize(const std::string& str,
              std::vector<std::string>& tokens,
              const std::string& delimiters)
{
    // Skip delimiters at beginning.
    std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
    // Find first "non-delimiter".
    std::string::size_type pos     = str.find_first_of(delimiters, lastPos);
    
    while (std::string::npos != pos || std::string::npos != lastPos)
    {
        // Found a token, add it to the vector.
        tokens.push_back(str.substr(lastPos, pos - lastPos));
        // Skip delimiters.  Note the "not_of"
        lastPos = str.find_first_not_of(delimiters, pos);
        // Find next "non-delimiter"
        pos = str.find_first_of(delimiters, lastPos);
    }
}



bool has_suffix(const std::string &str, const std::string &suffix)
{
    return str.size() >= suffix.size() &&
    str.compare(str.size() - suffix.size(), suffix.size(), suffix) == 0;
}

bool has_prefix(const std::string &str, const std::string &prefix)
{
    auto res = std::mismatch(prefix.begin(), prefix.end(), str.begin());
    return res.first == prefix.end();
}

// bool is_number(const std::string& s)
// {
//     return !s.empty() && std::find_if(s.begin(),
//         s.end(), [](char c) { return !std::isdigit(c); }) == s.end();
// }

bool is_number(const std::string& s)
{
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}


// reads the graphs from a file
int graphFromFile(const std::string &fileName, Graph &G,std::string &graphName, double & meanEdgeWeight, int & numNodes){
    std::string line;
    std::ifstream infile(fileName);
    // get name
    std::getline(infile, line);
    std::istringstream iss(line);
    std::string ignoredValue;
    std::string type;
    iss>>ignoredValue>>graphName;
    infile>>ignoredValue>>type;
    double totalEdgeWeight = 0;
    int nodesAdded = 0;
    int edgesAdded = 0;
    if(type.compare("BikeTSP") == 0){
        // parse bike tsp
        int vertex;
        int head;
        int tail;
        double distance;
        std::set<int> nodesWithEdges;
        std::set<int> vertexSet;
        
        // skip ahead to node section
        do
        {
            infile>>ignoredValue;
        } while (ignoredValue.compare("DISPLAY_DATA_SECTION") != 0);
        
        // parse vertices, skip DISPLAY_DATA_SECTION line
        std::getline(infile, line);
        std::getline(infile, line);
        while (line.compare("EDGE_WEIGHT_SECTION")!=0)
        {
            std::istringstream iss(line);
            iss >> vertex;
            G.addVertex(vertex);
            vertexSet.insert(vertex);
            nodesAdded++;
            std::getline(infile, line);
        }
        
        
        while(std::getline(infile,line))
        {
            std::istringstream iss(line);
            iss >> head >> tail >> distance;
            if ((distance >= 5000) && (nodesWithEdges.find(head) != nodesWithEdges.end()) && (nodesWithEdges.find(tail) != nodesWithEdges.end()))
                continue;
            
            // if one of the two endpoints are not in the vertex set, continue
            if ((vertexSet.find(head) == vertexSet.end()) || (vertexSet.find(tail) == vertexSet.end())){
                std::cout<<"vertex not in graph"<<std::endl;
                continue;
            }
            nodesWithEdges.insert(head);
            nodesWithEdges.insert(tail);
            G.addEdge(head,tail,distance);
            totalEdgeWeight += distance;
            edgesAdded++;
            // std::cout<<head<<", "<<tail<<" "<<std::endl;
        }
    }
    else if(type.compare("BikeVertexWeighted") == 0){
        // parse bike vertex weighted tsp
        int vertex;
        int prize;
        int head;
        int tail;
        double distance;
        std::set<int> nodesWithEdges;
        std::set<int> vertexSet;
        
        // skip ahead to node section
        do
        {
            infile>>ignoredValue;
        } while (ignoredValue.compare("DISPLAY_DATA_SECTION") != 0);
        
        // parse vertices, skip DISPLAY_DATA_SECTION line
        std::getline(infile, line);
        std::getline(infile, line);
        while (line.compare("EDGE_WEIGHT_SECTION")!=0)
        {
            std::istringstream iss(line);
            iss >> vertex >> prize;
            G.addVertex(vertex,prize);
            vertexSet.insert(vertex);
            nodesAdded++;
            std::getline(infile, line);
        }
        
        
        while(std::getline(infile,line))
        {
            std::istringstream iss(line);
            iss >> head >> tail >> distance;
            if ((distance >= 5000) && (nodesWithEdges.find(head) != nodesWithEdges.end()) && (nodesWithEdges.find(tail) != nodesWithEdges.end()))
                continue;
            
            // if one of the two endpoints are not in the vertex set, continue
            if ((vertexSet.find(head) == vertexSet.end()) || (vertexSet.find(tail) == vertexSet.end())){
                std::cout<<"vertex not in graph"<<std::endl;
                continue;
            }
            nodesWithEdges.insert(head);
            nodesWithEdges.insert(tail);
            G.addEdge(head,tail,distance);
            totalEdgeWeight += distance;
            edgesAdded++;
            // std::cout<<head<<", "<<tail<<" "<<std::endl;
        }


    }
    else{
        // parse tsp from python
        int dimension;
        int edgeWeight;
        // std::getline(infile, line);
        // std::getline(infile, line);
        infile>>ignoredValue>>dimension;
        
        for (int i = 0; i < dimension; ++i){
            G.addVertex(i);
            nodesAdded++;
        }
        std::getline(infile,line);
        int head;
        int tail;
        double distance;
        std::string line;
        // for (int i=0; i<dimension; ++i){
        // for (int j=0; j<dimension; ++j){
        while(std::getline(infile,line)){
            std::istringstream iss(line);
            iss >> head >> tail >> distance;
            // std::cout<<"adding edge:"<<head<<","<<tail<<","<<distance<<std::endl;
            G.addEdge(head,tail,distance);
            totalEdgeWeight += distance;
            edgesAdded++;
        }
    }
    meanEdgeWeight = totalEdgeWeight/edgesAdded;
    //std::cout<<meanEdgeWeight<<std::endl;
    
    //std::cout<<"graph:"<<graphName<<std::endl;
    numNodes = nodesAdded;
    std::cout<<"nodes added:"<<nodesAdded<<"\n"
    <<"edges added:"<<edgesAdded<<std::endl;
    return 0;
    
}

// read distances from csv file
std::map<std::pair<int,int>,double> readDistances(const std::string &fileName)
{
    std::map<std::pair<int,int>,double> dists;
    std::ifstream infile;
    infile.open(fileName);
    std::string nextline;
    
    if (infile.is_open())
    {
        while ( getline (infile,nextline) ){
            // Get line and break by spaces
            std::vector<std::string> words;
            tokenize(nextline, words, ",");

            int i = 100*std::stoi(words[0]);
            int j = 100*std::stoi(words[1]);
            double weight = std::stod(words[2]);
            dists[std::make_pair(i,j)] = weight;
            dists[std::make_pair(j,i)] = weight;
        }
        infile.close();
    }
    
    return dists;
}
