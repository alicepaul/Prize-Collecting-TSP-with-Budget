//
//  ReadFile.hpp
//  
//
//  Created by Alice Paul on 4/10/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//
#include <stdio.h>
#include <map>
#include <list>
#include <iostream>
#include <fstream>
#include <memory>
#include <set>
#include <string>
#include <vector>
#include <sstream>
#include <fstream>

#include <sys/types.h>
#include <dirent.h>
#include <unistd.h>


#include "Graph.hpp"
#include "Subset.hpp"
#include "Subroutine.hpp"
#include "PD.hpp"

#ifndef ReadFile_hpp
#define ReadFile_hpp

#include <stdio.h>

void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);

bool has_suffix(const std::string &str, const std::string &suffix);

bool has_prefix(const std::string &str, const std::string &prefix);

bool is_number(const std::string& s);

int graphFromFile(const std::string &fileName, Graph &G,std::string &graphName, double & meanEdgeWeight, int & numNodes);

std::map<std::pair<int,int>,double> readDistances(const std::string &fileName);



#endif /* ReadFile_hpp */
