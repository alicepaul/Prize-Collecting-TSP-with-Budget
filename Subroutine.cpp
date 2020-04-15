//
//  Subroutine.cpp
//  Prize-Collecting TSP
//
//  Created by Alice Paul on 3/14/17.
//
//  MIT License
//  Copyright (c) 2020 alicepaul
//


#include "Subroutine.hpp"
double eps = 1.0e-15; // Level of precision for general ties
double tieeps = 0.001; // Level of precisions to count ties

/* ------------------------- HELPER FUNCTIONS --------------------------*/


// Returns whether or not an edge is active
bool activeEdge(std::shared_ptr<Edge> & e, std::map<int,std::shared_ptr<Subset>> & vertexSubs)
{
    std::shared_ptr<Subset> p1 = vertexSubs[e->getHead()], p2=vertexSubs[e->getTail()];
    if (p1 == p2){
        return false;
    }
    else if ((p1->getActive() == false) && (p2->getActive() == false)){
        return false;
    }
    else{
        return true;
    }
}

// Returns 0.5 if both endpoints are active and 1.0 otherwise
double edgeFactor(std::shared_ptr<Edge> & e, std::map<int,std::shared_ptr<Subset>> & vertexSubs)
{
    std::shared_ptr<Subset> p1 = vertexSubs[e->getHead()], p2=vertexSubs[e->getTail()];
    bool p1active = p1->getActive(), p2active = p2->getActive();

    if (p1active && p2active) {
        return 0.5;
    }
    return 1.0;
}

// Find min event to occur next
void findMinEvent(double lambda, const std::list<std::shared_ptr<Subset>> & subsets,
                                 std::vector<std::shared_ptr<Edge>> & E, std::map<int,std::shared_ptr<Subset>> & vertexSubs,
                                 std::map<std::shared_ptr<Edge>,std::vector<double>> & lin_e,
                                 std::map<std::shared_ptr<Subset>,std::vector<double>> & lin_s,
                                 std::map<std::shared_ptr<Edge>,std::vector<double>> & lin_e_p,
                                 std::map<std::shared_ptr<Subset>,std::vector<double>> & lin_s_p,
                                 std::shared_ptr<Subset> & min_s, std::shared_ptr<Edge> & min_e,
                                 std::vector<double> & lin_val,
                                 std::shared_ptr<Edge> & alt_e,
                                 std::vector<double> & lin_val_p1, std::vector<double> & lin_val_p2)
{
    min_e = NULL;
    min_s = NULL;
    
    // Iterate through active subsets
    double time_s = INT_MAX;
    std::vector<double> lin_val_s = {0.,0.};
    for (auto s:subsets){
        if (s->getActive()){
            // Time to go tight for subset s
            double test = lambda*(1-tieeps)*lin_s[s][0]+lin_s[s][1];
            
            // If new min
            if (test < time_s-eps){ // break ties by slope
                time_s = test;
                lin_val_s = lin_s[s];
                min_s = s;
            }
        }
    }

    
    // Iterate through active edges
    double time_e = INT_MAX;
    std::vector<double> lin_val_e = {0.,0.};
    for (auto e:E){
        if (activeEdge(e,vertexSubs)){
            // Time to go tight for edge e
            double factor = edgeFactor(e,vertexSubs);
            double test = factor*(lambda*(1-tieeps)*lin_e[e][0]+lin_e[e][1]);
        
            // If new min
            if (test < time_e-eps) {
                time_e = test;
                lin_val_e = {factor*lin_e[e][0], factor*lin_e[e][1]};
                min_e = e;
            }
        }
    }
    
    // Find which event is first and which tied
    if (time_e < time_s+eps){
        min_s = NULL;
        lin_val = lin_val_e;
    }
    else{
        min_e = NULL;
        lin_val = lin_val_s;
    }
    
    // If an edge event then we have to check for ties
    lin_val_p1 = {0.,0.}, lin_val_p2 = {lin_val[0], lin_val[1]};
    alt_e = min_e;
    if (min_e != NULL){
        // Find event time for min_e at lambda*(1+tieeps)
        double time_p = edgeFactor(min_e,vertexSubs)*(lambda*(1+tieeps)*lin_e_p[min_e][0]+lin_e_p[min_e][1]);
        
        // Find minimum tied edge between same subsets at lambda*(1+tieeps)
        std::shared_ptr<Subset> p1 = vertexSubs[min_e->getHead()], p2 = vertexSubs[min_e->getTail()];
        for (auto e:E){
            std::shared_ptr<Subset> s1 = vertexSubs[e->getHead()], s2 = vertexSubs[e->getTail()];
            if (((s1==p1)&&(s2==p2)) || ((s1==p2)&&(s2==p1))){
            
                // Time to go tight for edge e
                double factor = edgeFactor(e,vertexSubs);
                double test = factor*(lambda*(1+tieeps)*lin_e_p[e][0]+lin_e_p[e][1]);
                double test_min = factor*(lambda*(1-tieeps)*lin_e[e][0]+lin_e[e][1]);
                
                // If new min at lambda*(1+tieeps)
                if (test < time_p-eps) {
                    lin_val_p2 = {factor*lin_e_p[e][0], factor*lin_e_p[e][1]};
                    alt_e = e;
                    time_p = test;
                }
            }
        }
        
        // Find if parents go neutral before the edge at lambda*(1+eps)
        double testp1 = INT_MAX, testp2 = INT_MAX;
        bool tiedp1 = false, tiedp2 = false;
        if (p1->getActive()){
            testp1 = lambda*(1+tieeps)*lin_s_p[p1][0]+lin_s_p[p1][1];
            tiedp1 = (testp1 < time_p-eps);
            p1->setTied(tiedp1);
        }
        if (p2->getActive()){
            testp2 = lambda*(1+tieeps)*lin_s_p[p2][0]+lin_s_p[p2][1];
            tiedp2 = (testp2 < time_p-eps);
            p2->setTied(tiedp2);
        }
            
//        if (tiedp1){
//            std::cout << "Tied1 " << *p1;
//        }
//        if (tiedp2){
//            std::cout << "Tied2 " << *p2;
//        }

        // P1 ties and P2 is active - edge will go tight next
        if ((tiedp1) && (tiedp2 == false) && (p2->getActive())){
            lin_val_p1 = lin_s_p[p1];
            lin_val_p2 = {lin_e_p[min_e][0]-2*lin_s_p[p1][0], lin_e_p[min_e][1]-2*lin_s_p[p1][1]};
        }
        // P2 ties and P1 is active - edge will go tight next
        else if ((tiedp2) && (tiedp1 == false) && (p1->getActive())){
            lin_val_p1 = lin_s_p[p2]; 
            lin_val_p2 = {lin_e_p[min_e][0]-2*lin_s_p[p2][0], lin_e_p[min_e][1]-lin_s_p[p2][1]};
        }
        // P1 ties and P2 is inactive or goes neutral right after - edge between two inactive components
        else if ((tiedp1) && (testp1 < testp2)){
            alt_e = min_e;
            lin_val_p1 = {lin_val[0], lin_val[1]}, lin_val_p2 = {0.,0.};
        }
        // P2 ties and P1 is inactive or goes neutral right after - edge between two inactive components
        else if ((tiedp2) && (testp2 < testp1)){
            alt_e = min_e;
            lin_val_p1 = {lin_val[0], lin_val[1]}, lin_val_p2 = {0.,0.};
        }
        
    }
}


/* ------------------------- MAIN FUNCTIONS --------------------------*/


// Grow Function
std::list<std::shared_ptr<Subset>> growSubsets(const Graph & G, double lambda)
{
    const std::list<int>& V = G.getVertices();                      // vertices in G
    std::vector<std::shared_ptr<Edge>> E;                           // edges in G
    for (auto e: G.getEdges()){
        E.push_back(e);
    }
    
    std::list<std::shared_ptr<Subset>> subsets;                     // list current subsets
    std::map<int,std::shared_ptr<Subset>> vertexSubs;               // map vertices to subset that contains them
    
    // Time to go tight is alpha*lambda+beta
    std::map<std::shared_ptr<Subset>,std::vector<double>> lin_s;        // a and b values for subsets
    std::map<std::shared_ptr<Edge>,std::vector<double>> lin_e;          // a and b values for edges
    std::map<std::shared_ptr<Subset>,std::vector<double>> lin_s_p;      // a and b values for subsets at lambda+
    std::map<std::shared_ptr<Edge>,std::vector<double>> lin_e_p;        // a and b values for edges at lambda+
    
    // First create an active subset for each vertex and intialize a_s and b_s
    for (auto v: V){
        int prize = G.getVertex(v)->getPrize();
        std::shared_ptr<Subset> Sp = std::make_shared<Subset>(v,prize,prize);
        vertexSubs[v] = Sp;
        subsets.push_back(Sp);
        lin_s[Sp] = {0, 0.5*prize};               // time_s = 0*lambda+0.5
        lin_s_p[Sp] = {0, 0.5*prize};
    }
    // Next initialize a_e and b_e for each edge
    for (auto e:E){
        lin_e[e] = {e->getWeight(), 0};     // time_e = w*lambda + 0
        lin_e_p[e] = {e->getWeight(),0};
    }
    
    while(true)
    {
        //std::cout << "---------Next iteration-----------\n";
        
        // Find next event
        std::shared_ptr<Subset> min_s = NULL;                      // First subset to go tight
        std::shared_ptr<Edge> min_e = NULL, alt_e = NULL;          // First edge to go tight and alt edge
        std::vector<double> lin_val = {0,0};                       // Linear functions for time to go tight
        std::vector<double> lin_val_p1 = {0,0}, lin_val_p2 = {0,0};// Linear functions to go tight at lambda+
        findMinEvent(lambda, subsets, E, vertexSubs, lin_e, lin_s, lin_e_p, lin_s_p, min_s, min_e, lin_val, alt_e, lin_val_p1, lin_val_p2);
        
        
        // If nothing to go tight - then algorithm is done
        if ((min_s == NULL) && (min_e == NULL)){
            return subsets;
        }

        // Update y_vals and times to go tight
        for (auto s : subsets){
            if (s->getActive()){
                s->setY(s->getY() + lin_val[0]*lambda*(1-tieeps)+lin_val[1]);
                lin_s[s][0] -= lin_val[0], lin_s[s][1] -= lin_val[1];
                
                // If tied then only raise first linear amount p1
                if (s->getTied()){
                    lin_s_p[s][0] -= lin_val_p1[0], lin_s_p[s][1] -= lin_val_p1[1];
                }
                else{
                    lin_s_p[s][0] -= (lin_val_p1[0]+lin_val_p2[0]);
                    lin_s_p[s][1] -= (lin_val_p1[1]+lin_val_p2[1]);
                }
            }
        }
        std::vector<std::shared_ptr<Edge>>::iterator it;
        for (it=E.begin(); it != E.end();){
            std::shared_ptr<Edge> e = *it;
            if (activeEdge(e,vertexSubs)){
                double factor = edgeFactor(e,vertexSubs);
                lin_e[e][0] -= lin_val[0]/factor, lin_e[e][1] -= lin_val[1]/factor;
                
                // If parent 1 is tied only raise p1 otherwise p1+p2
                if (vertexSubs[e->getHead()]->getTied()){
                    lin_e_p[e][0] -= lin_val_p1[0], lin_e_p[e][1] -= lin_val_p1[1];
                }
                else if (vertexSubs[e->getHead()]->getActive()){
                    lin_e_p[e][0] -= (lin_val_p1[0]+lin_val_p2[0]);
                    lin_e_p[e][1] -= (lin_val_p1[1]+lin_val_p2[1]);
                }
            
                // If parent 2 is tied only raise p1 otherwise p1+p2
                if (vertexSubs[e->getTail()]->getTied()){
                    lin_e_p[e][0] -= lin_val_p1[0], lin_e_p[e][1] -= lin_val_p1[1];
                }
                else if (vertexSubs[e->getTail()]->getActive()){
                    lin_e_p[e][0] -= (lin_val_p1[0]+lin_val_p2[0]);
                    lin_e_p[e][1] -= (lin_val_p1[1]+lin_val_p2[1]);
                }
            }
            if (vertexSubs[e->getHead()] == vertexSubs[e->getTail()]){
                E.erase(it);
            }
            else{
                ++it;
            }
        }

        //std::cout << "------------ EVENT ------------ \n";

        // If subset goes neutral first - mark inactive
        if (min_s != NULL){
            //std::cout << " Subset Event " << *min_s;
            min_s->setActive(false);
        }
        
        // Else an edge goes tight first - merge subsets
        else{
            //if (min_e != alt_e)
            //    std::cout << " Edge Event " << *min_e  << " with alt edge " << *alt_e << "\n";
            //else
            //    std::cout << "Edge Event " << *min_e << "\n";
            std::shared_ptr<Subset> S1 = vertexSubs[min_e->getHead()];
            std::shared_ptr<Subset> S2 = vertexSubs[min_e->getTail()];
            std::shared_ptr<Subset> S = std::make_shared<Subset>(S1,S2,min_e,alt_e);
            lin_s[S] = {lin_s[S1][0]+lin_s[S2][0], lin_s[S1][1]+lin_s[S2][1]};
            lin_s_p[S] = {lin_s_p[S1][0]+lin_s_p[S2][0], lin_s_p[S1][1]+lin_s_p[S2][1]};
            
            subsets.remove(S1);
            subsets.remove(S2);
            subsets.push_back(S);
            
            // Update vertexSubs
            const std::list<int>& vinS = S->getVertices();
            for (auto v:vinS){
                vertexSubs[v] = S;
            }
        }
        
    }
}


// Prune Function
double prune(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges, bool l_plus, bool swap)
{
    // Base case is that s does not have a parent
    if (s->getParent1() == NULL){
        return 0;
    }
    
    // Otherwise find parents and number of endpoints
    std::shared_ptr<Subset> p1 = s->getParent1(), p2 = s->getParent2();
    std::shared_ptr<Edge> e = s->getEdge(), alt_e = s->getAltEdge();
    if ((l_plus==true) && (alt_e != NULL) && (swap == true)){
        e = alt_e;
    }
    
    // Find degree of p1 given current edges
    int connected1 = 0;
    for (auto v : p1->getVertices()){
        for (auto e:edges){
            if ((e->getHead() == v) || (e->getTail() == v))
                connected1 += 1;
        }
    }
    
    // Find degree of p2 given current edges
    int connected2 = 0;
    for (auto v : p2->getVertices()){
        for (auto e:edges){
            if ((e->getHead() == v) || (e->getTail() == v))
                connected2 += 1;
        }
    }
    
    // Find if active. If l_plus = true don't include those that tied
    bool p1active = p1->getActive();
    bool p2active = p2->getActive();
    if (l_plus){
        if (p1->getTied()){
            p1active = false;
        }
        if (p2->getTied()){
            p2active = false;
        }
    }
    
    // If p1 and p2 are active or neither needs to be pruned, add edge between and recurse
    if (((p1active && p2active) || (p1active && connected2))
        || (p2active && connected1) || ((connected1 > 0) && (connected2 > 0)))
    {
        edges.push_back(e);
        double w1 = prune(p1,edges,l_plus,swap);
        double w2 = prune(p2,edges,l_plus,swap);
        return w1+w2+e->getWeight();
    }
    
    // Else if p1 is active or does need to be pruned, recurse on p1
    else if ( (p1active) || (connected1 > 1))
    {
        //std::cout << "Prune " << *p2;
        return prune(p1,edges,l_plus,swap);
    }
    
    // Else if p2 is active or does need to be pruned, recurse on p2
    else if ( (p2active) || (connected2 > 1))
    {
        //std::cout << "Prune " << *p1;
        return prune(p2,edges,l_plus,swap);
    }
    
    // Else if no edges are added yet - find largest and return
    else if (edges.size() == 0){
        std::list<std::shared_ptr<Edge>> edges1, edges2;
        double w1 = prune(p1,edges1,l_plus,swap), w2 = prune(p2,edges2,l_plus,swap);
        if (w1 > w2){
            edges = edges1;
            return w1;
        }
        else{
            edges = edges2;
            return w2;
        }
    }
    
    return 0;
}

// Find path of maximal pruned subsets that are not pruned when s is marked active
void maxPrunedS(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
                std::shared_ptr<Subset> & tied_S, std::vector<std::shared_ptr<Subset>> & prunedS,
                std::vector<std::shared_ptr<Edge>> & prunedE)
{
    // Base case is that s does not have a parent
    if (s->getParent1() == NULL){
        return;
    }
    
    // Otherwise find parents and number of endpoints
    std::shared_ptr<Subset> p1 = s->getParent1(), p2 = s->getParent2();
    int u = tied_S->getVertices().back();
    
    // Find if p1 has another incident edge already included and searched pruned edges
    int connected1 = 0, connectedTied1 = 0;
    std::shared_ptr<Edge> connectedE1 = NULL;
    bool inp1 = false;
    for (auto v : p1->getVertices()){
        for (auto e:edges){
            if ((e->getHead() == v) || (e->getTail() == v)){
                connectedE1 = e;
                connected1 += 1;
            }
        }
        if (v == u)
            inp1 = true;
        for (auto e: prunedE){
            if ((e->getHead() == v) || (e->getTail() == v))
                connectedTied1 += 1;
        }
    }
    
    // Find if p2 has another incident edge already included and search pruned edges
    int connected2 = 0, connectedTied2 = 0;
    std::shared_ptr<Edge> connectedE2 = NULL;
    bool inp2 = false;
    for (auto v : p2->getVertices()){
        for (auto e:edges){
            if ((e->getHead() == v) || (e->getTail() == v)){
                connectedE2 = e;
                connected2 += 1;
            }
        }
        if (v == u)
            inp2 = true;
        for (auto e: prunedE){
            if ((e->getHead() == v) || (e->getTail() == v))
                connectedTied2 += 1;
        }
    }
    
    // Find if active. If l_plus = true don't include those that tied
    bool p1active = p1->getActive();
    bool p2active = p2->getActive();
    
    // If p1 and p2 are active or neither needs to be pruned, add edge between and recurse
    if (((p1active && p2active) || (p1active && connected2))
        || (p2active && connected1) || ((connected1 > 0) && (connected2 > 0)))
    {
        //std::cout << "Case 1 \n";
        //std::cout << "Edge added " << *(s->getEdge()) << "\n";
        edges.push_back(s->getEdge());
        maxPrunedS(p1,edges,tied_S,prunedS,prunedE);
        maxPrunedS(p2,edges,tied_S,prunedS,prunedE);
    }
    
    // Else if p1 is active, recurse on p1 - check whether p2 would have been pruned in other case
    else if ( (p1active) || (connected1 > 1))
    {
        //std::cout << "Case 2 \n";
        if ((connectedTied2 > 0) || (tied_S == p2)){
            prunedE.insert(prunedE.begin(),s->getEdge());
            prunedS.insert(prunedS.begin(),p2);
        }
        maxPrunedS(p1,edges,tied_S,prunedS,prunedE);
    }
    
    // Else if p2 is active, recurse on p2 - check whether p1 would have been pruned in other case
    else if ( (p2active) || (connected2 > 1))
    {
        //std::cout << "Case 3 \n";
        if ((connectedTied1 > 0) || (tied_S == p1)){
            prunedE.insert(prunedE.begin(),s->getEdge());
            prunedS.insert(prunedS.begin(),p1);
        }
        maxPrunedS(p2,edges,tied_S,prunedS,prunedE);
    }
    
    // Else if p1 == tied_S and connected by a single edge - would have kept otherwise
    else if ((p1 == tied_S) && (connected1 == 1)){
        //std::cout << "Case 4 \n";
        prunedE.insert(prunedE.begin(),connectedE1);
        prunedS.insert(prunedS.begin(),p1);
    }
    
    // Else if p1 == tied_S and p2 is inactive and connected by a single edge
    else if ((p1 == tied_S) && (connected2 == 1)){
        //std::cout << "Case 5 \n";
        prunedE.push_back(connectedE2);
        prunedS.push_back(p2);
        prunedE.push_back(s->getEdge());
        prunedS.push_back(p1);
    }
    
    // Else if p1 == tied_S and connected by a single edge - would have kept otherwise
    else if ((p2 == tied_S) && (connected2 == 1)){
        //std::cout << "Case 6 \n";
        prunedE.insert(prunedE.begin(),connectedE2);
        prunedS.insert(prunedS.begin(),p2);
    }
    
    // Else if p2 == tied_S and p1 is inactive and connected by a single edge
    else if ((p2 == tied_S) && (connected1 == 1)){
        //std::cout << "Case 7 \n";
        prunedE.push_back(connectedE1);
        prunedS.push_back(p1);
        prunedE.push_back(s->getEdge());
        prunedS.push_back(p2);
    }
    
    // Else if both p1 and p2 are inactive then we need to find which one contains test_s
    else if (inp1){
        maxPrunedS(p1,edges,tied_S,prunedS,prunedE);
    }
    else if (inp2){
        maxPrunedS(p2,edges,tied_S,prunedS,prunedE);
    }
}

// Find path of maximal pruned subsets that are not pruned when s uses alt edge
void maxPrunedE(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
                std::shared_ptr<Subset> & tied_S, std::vector<std::shared_ptr<Subset>> & prunedS,
                std::vector<std::shared_ptr<Edge>> & prunedE)
{
    // Base case is that s does not have a parent
    if (s->getParent1() == NULL){
        return;
    }
    
    // Otherwise find parents and number of endpoints
    std::shared_ptr<Subset> p1 = s->getParent1(), p2 = s->getParent2();
    std::shared_ptr<Edge> alt_e = tied_S->getAltEdge();
    int u = tied_S->getVertices().back();
    
    // Find if p1 has another incident edge already included and searched pruned edges
    int connected1 = 0, connectedTied1 = 0;
    bool inp1 = false;
    for (auto v : p1->getVertices()){
        for (auto e:edges){
            if ((e->getHead() == v) || (e->getTail() == v))
                connected1 += 1;
        }
        if (v == u)
            inp1 = true;
        for (auto e: prunedE){
            if ((e->getHead() == v) || (e->getTail() == v))
                connectedTied1 += 1;
        }
        if ((alt_e->getHead() == v) || (alt_e->getTail() == v))
            connectedTied1 += 1;
    }
    
    // Find if p2 has another incident edge already included and search pruned edges plus the alt edge
    int connected2 = 0, connectedTied2 = 0;
    bool inp2 = false;
    for (auto v : p2->getVertices()){
        for (auto e:edges){
            if ((e->getHead() == v) || (e->getTail() == v))
                connected2 += 1;
        }
        if (v == u)
            inp2 = true;
        for (auto e: prunedE){
            if ((e->getHead() == v) || (e->getTail() == v))
                connectedTied2 += 1;
        }
        if ((alt_e->getHead() == v) || (alt_e->getTail() == v))
            connectedTied2 += 1;
    }
    
    // Find if active. If l_plus = true don't include those that tied
    bool p1active = p1->getActive();
    bool p2active = p2->getActive();
    
    // If s == tied_S then keep e and recurse on both parents
    if (s == tied_S){
        edges.push_back(s->getEdge());
        maxPrunedE(p1,edges,tied_S,prunedS,prunedE);
        prunedE.push_back(alt_e);                       // Push alt_e onto back of path then continue
        maxPrunedE(p2,edges,tied_S,prunedS,prunedE);
    }
    
    // If p1 and p2 are active or neither needs to be pruned, add edge between and recurse
    else if (((p1active && p2active) || (p1active && connected2))
        || (p2active && connected1) || ((connected1 > 0) && (connected2 > 0)))
    {
        edges.push_back(s->getEdge());
        maxPrunedE(p1,edges,tied_S,prunedS,prunedE);
        maxPrunedE(p2,edges,tied_S,prunedS,prunedE);
    }
    
    // Else if p1 is active, recurse on p1 - check whether p2 would have been pruned in other case
    else if ( (p1active) || (connected1 > 1)){
        if (connectedTied2 > 0){
            prunedE.push_back(s->getEdge());
            prunedS.push_back(p2);
        }
        maxPrunedE(p1,edges,tied_S,prunedS,prunedE);
    }
    
    // Else if p2 is active, recurse on p2 - check whether p1 would have been pruned in other case
    else if ( (p2active) || (connected2 > 1)){
        if (connectedTied1 > 0){
            prunedE.push_back(s->getEdge());
            prunedS.push_back(p1);
        }
        maxPrunedE(p2,edges,tied_S,prunedS,prunedE);
    }
    
    // Else if both p1 and p2 are inactive then we need to find which one contains test_s
    else if (inp1){
        maxPrunedE(p1,edges,tied_S,prunedS,prunedE);
    }
    else if (inp2){
        maxPrunedE(p2,edges,tied_S,prunedS,prunedE);
    }
}

// Reverse Delete Function
double reverseDelete(std::list<std::shared_ptr<Subset>> & subsets, std::list<std::shared_ptr<Edge>> & edges,
                     std::shared_ptr<Subset> & largestSub, bool l_plus, bool swap_edges)
{
    // Iterate through components to prune and return the largest
    double largest = -INT_MAX;
    largestSub = NULL;
    for (auto s:subsets){
        std::list<std::shared_ptr<Edge>> edgesS;
        double sizeComp = prune(s,edgesS,l_plus,swap_edges);
        if (sizeComp > largest){
            largestSub = s;
            largest = sizeComp;
            edges = edgesS;
        }
        
    }
    return largest;
}


/* ------------------------- ALTERNATE FUNCTIONS --------------------------*/


double reverseDelete(std::list<std::shared_ptr<Subset>> & subsets, bool l_plus, bool swap_edges)
{
    std::list<std::shared_ptr<Edge>> edges;
    std::shared_ptr<Subset> largestSub = NULL;
    return reverseDelete(subsets, edges, largestSub, l_plus, swap_edges);
}

double reverseDelete(std::shared_ptr<Subset> & s, std::list<std::shared_ptr<Edge>> & edges,
                     bool l_plus, bool swap_edges)
{
    double w = prune(s,edges,l_plus,swap_edges);
    return w;
}

double reverseDelete(std::shared_ptr<Subset> & s, bool l_plus, bool swap_edges)
{
    std::list<std::shared_ptr<Edge>> edges;
    return reverseDelete(s, edges, l_plus, swap_edges);
}






