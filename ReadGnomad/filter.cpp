//
//  filter.cpp
//  ReadGnomad
//
//  Created by Gang Peng on 4/26/20.
//  Copyright © 2020 Gang Peng. All rights reserved.
//

#include "filter.hpp"
#include "normal.h"
#include <iostream>

using namespace std;


int checkFilter(std::map<std::string, std::vector<std::string> > cmLine,
                 std::map<std::string, double> & filterLarger,
                 std::map<std::string, double> & filterSmaller,
                 std::map<std::string, std::vector<std::string> > & filterNotEqual) {
    int rlt = 0;
    filterLarger.clear();
    filterSmaller.clear();
    filterNotEqual.clear();
    
    map<string, vector<string> >::iterator it;
    it = cmLine.find("-filterLarger");
    if(it != cmLine.end()){
        if(it->second.size() > 0) {
            if(it->second.size() % 2 != 0){
                cout << "The parameter set in -filterLarger is not paired. It should be -filterLarger id1 cutoff1 id2 cutoff2 ..." << endl;
                return -1;
            }
            
            for(int i=0; i<it->second.size()/2; i++){
                filterLarger[it->second[i*2]] = stod(it->second[i*2+1]);
            }
            
            rlt++;
        }
    }
    
    
    it = cmLine.find("-filterSmaller");
    if(it != cmLine.end()){
        if(it->second.size() > 0) {
            if(it->second.size() % 2 != 0){
                cout << "The parameter set in -filterSmaller is not paired. It should be -filterSmaller id1 cutoff1 id2 cutoff2 ..." << endl;
                return -1;
            }
            
            for(int i=0; i<it->second.size()/2; i++){
                filterSmaller[it->second[i*2]] = stod(it->second[i*2+1]);
            }
            
            rlt++;
        }
    }
    
    it = cmLine.find("-filterNotEqual");
    if(it != cmLine.end()){
        if(it->second.size() > 0) {
            if(it->second.size() % 2 != 0){
                cout << "The parameter set in -filterNotEqual is not paired. It should be -filterNotEqual id1 var1 id2 var2 ..." << endl;
                return -1;
            }
            
            for(int i=0; i<it->second.size()/2; i++){
                filterNotEqual[it->second[i*2]] = split(it->second[i*2+1], "&", true);
            }
            
            rlt++;
        }
    }
    
    return rlt;
}


int checkVEPFilter(std::map<std::string, std::vector<std::string> > cmLine,
                   std::map<std::string, std::vector<std::string> > & filterVEPNotEqual){
    int rlt = 0;
    
    map<string, vector<string> >::iterator it;
    it = cmLine.find("-filterVEPNotEqual");
    if(it != cmLine.end()){
        if(it->second.size() > 0) {
            if(it->second.size() % 2 != 0){
                cout << "The parameter set in -filterVEPNotEqual is not paired. It should be -filterVEPNotEqual id1 var1 id2 var2 ..." << endl;
                return -1;
            }
            
            for(int i=0; i<it->second.size()/2; i++){
                filterVEPNotEqual[it->second[i*2]] = split(it->second[i*2+1], "&", true);
            }
            
            rlt++;
        }
    }
    
    return rlt;
}
