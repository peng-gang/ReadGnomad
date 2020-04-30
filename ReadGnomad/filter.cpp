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



/* check for range
 tbx_itr_next(htsfp, tbx, itr, r)   hts_itr_next(hts_get_bgzfp(htsfp), (itr), (r), (tbx))

 static inline int bcf_itr_next(htsFile *htsfp, hts_itr_t *itr, void *r) {
         if (htsfp->is_bgzf)
             return hts_itr_next(htsfp->fp.bgzf, itr, r, 0);

         hts_log_error("Only bgzf compressed files can be used with iterators");
         errno = EINVAL;
         return -2;
     }


  #define tbx_itr_querys(tbx, s) hts_itr_querys((tbx)->idx, (s), (hts_name2id_f)(tbx_name2id), (tbx), hts_itr_query, tbx_readrec)

 #define bcf_itr_querys(idx, hdr, s) hts_itr_querys((idx), (s), (hts_name2id_f)(bcf_hdr_name2id), (hdr), hts_itr_query, bcf_readrec)


 int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end);
 int tbx_readrec(BGZF *fp, void *tbxv, void *sv, int *tid, hts_pos_t *beg, hts_pos_t *end)
 {
     tbx_t *tbx = (tbx_t *) tbxv;
     kstring_t *s = (kstring_t *) sv;
     int ret;
     if ((ret = bgzf_getline(fp, '\n', s)) >= 0) {
         tbx_intv_t intv;
         if (get_intv(tbx, s, &intv, 0) < 0)
             return -2;
         *tid = intv.tid; *beg = intv.beg; *end = intv.end;
     }
     return ret;
 }


 int bcf_readrec(BGZF *fp, void *null, void *v, int *tid, hts_pos_t *beg, hts_pos_t *end);
 int bcf_readrec(BGZF *fp, void *null, void *vv, int *tid, int *beg, int *end)
 {
     bcf1_t *v = (bcf1_t *) vv;
     int ret;
     if ((ret = bcf_read1_core(fp, v)) >= 0)
         *tid = v->rid, *beg = v->pos, *end = v->pos + v->rlen;
     return ret;
 }
 */
