//#include <config.h>

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

#include "normal.h"
#include "filter.hpp"

using namespace std;

void error(const char *format, ...)
{
    va_list ap;
    va_start(ap, format);
    vfprintf(stderr, format, ap);
    va_end(ap);
    if (strrchr(format, '\n') == NULL) fputc('\n', stderr);
    exit(-1);
}

// how to use index file
// bcf_itr_querys
// bcf_itr_queryi
// hts_itr_query
// hts_itr_next
// check later

int main(int argc, char **argv){
/*
 * -i input.vcf.bgz
 * -o output.txt (tab delimited)
 * -info AF AC ...  // id in information section to include in output file, please check VCF header file for details
 * -vep SYMBOL Gene ... // id in VEP section to include in output file, please check VCF header for details
 * -filterLarger AF 0.2 ... // filter out all variants with AF > 0.2
 * -filterSmaller AF 0.01 ...   // filter out all variants with AF < 0.01
 * -filterNotEqual variant_type snv ...  // include snv only
 * -filterVEPNotEqual SYMBOL ACE&ACE2 ...  // include gene ACE and ACE2 only
 * -range chr1:1000-20000  // only consider the variants in chromosome 1 from 1000 to 20000
 */
    vector<string> mustOptions = {"-i", "-o"};
    vector<string> allOptions = {"-i", "-o", "-info", "-vep", "-filterLarger", "-filterSmaller", "-filterNotEqual", "-filterVEPNotEqual", "-range"};
    
    
    // process filter parameter
    bool filter = false;
    bool filterVEP = false;
    
    
    map<string, vector<string> > cmLine = parseCMLine(argc, argv, allOptions, mustOptions);
    for(map<string, vector<string> >::iterator it = cmLine.begin(); it != cmLine.end(); it++){
        cout << it->first << endl;
        for(vector<string>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
            cout << *it2 << "\t";
        }
        cout << endl << endl;
    }
    
    
    // check filter
    std::map<std::string, double> filterLarger;
    std::map<std::string, double> filterSmaller;
    std::map<std::string, std::vector<std::string> > filterNotEqual;
    
    int rt = checkFilter(cmLine, filterLarger, filterSmaller, filterNotEqual);
    if(rt > 0){
        filter = true;
    } else if(rt < 0){
        cout << "Error during parsing input filter parameters" << endl;
        return -1;
    }
    
    // check VEP filter
    std::map<std::string, std::vector<std::string> > filterVEPNotEqual;
    rt = checkVEPFilter(cmLine, filterVEPNotEqual);
    if(rt > 0){
        filterVEP = true;
    } else if(rt < 0){
        cout << "Error during parsing input VEP filter parameters" << endl;
        return -1;
    }
    
    
    string fname = cmLine["-i"][0];
    htsFile *fp    = hts_open(fname.c_str(),"rb");
    if (!fp) error("Failed to open \"%s\" : %s", fname.c_str(), strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) error("bcf_hdr_read : %s", strerror(errno));
    
    
    // check header information
    vector<int> idxInfo;
    if(cmLine.find("-info") != cmLine.end()){
        for(size_t i=0; i<cmLine["-info"].size(); i++){
            int idxTmp = bcf_hdr_id2int(hdr, BCF_DT_ID, cmLine["-info"][i].c_str());
            if(idxTmp < 0){
                cout << "ID " << cmLine["-info"][i] << " is not included in the file. Please recheck header of vcf file." << endl;
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                return -1;
            }
            
            idxInfo.push_back(idxTmp);
        }
    }
    
    // get filter idx in header
    vector<int> idxFilterLarger;
    vector<int> idxFilterSmaller;
    vector<int> idxFilterNotEqual;
    if(filter){
        if(filterLarger.size() > 0){
            for(std::map<std::string, double>::iterator it = filterLarger.begin(); it!=filterLarger.end(); it++){
                int idxTmp = bcf_hdr_id2int(hdr, BCF_DT_ID, (it->first).c_str());
                if(idxTmp < 0){
                    cout << "Filter ID " << it->first << " is not included in the file. Please recheck header of vcf file.";
                    bcf_hdr_destroy(hdr);
                    hts_close(fp);
                    return -1;
                }
                idxFilterLarger.push_back(idxTmp);
            }
        }
        
        if(filterSmaller.size() > 0){
            for(std::map<std::string, double>::iterator it = filterSmaller.begin(); it!=filterSmaller.end(); it++){
                int idxTmp = bcf_hdr_id2int(hdr, BCF_DT_ID, (it->first).c_str());
                if(idxTmp < 0){
                    cout << "Filter ID " << it->first << " is not included in the file. Please recheck header of vcf file.";
                    bcf_hdr_destroy(hdr);
                    hts_close(fp);
                    return -1;
                }
                idxFilterSmaller.push_back(idxTmp);
            }
        }
        
        if(filterNotEqual.size() > 0){
            for(std::map<std::string, std::vector<std::string> >::iterator it = filterNotEqual.begin(); it!=filterNotEqual.end(); it++){
                int idxTmp = bcf_hdr_id2int(hdr, BCF_DT_ID, (it->first).c_str());
                if(idxTmp < 0){
                    cout << "Filter ID " << it->first << " is not included in the file. Please recheck header of vcf file.";
                    bcf_hdr_destroy(hdr);
                    hts_close(fp);
                    return -1;
                }
                idxFilterNotEqual.push_back(idxTmp);
            }
        }
        cout << "filter" << endl;
    }
    
    vector<size_t> idxVepTag;
    int idxVep = -1;
    if(cmLine.find("-vep") != cmLine.end()){
        if(cmLine["-vep"].size() > 0){
            idxVep = bcf_hdr_id2int(hdr, BCF_DT_ID, "vep");
            if(idxVep < 0){
                cout << "vep is not included in the file. Please recheck header of vcf file." << endl;
                bcf_hdr_destroy(hdr);
                hts_close(fp);
                return -1;
            }
            
            bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "vep", NULL);
            for(int i=0; i<hrec->nkeys; i++){
                if(string(hrec->keys[i]) == "Description"){
                    string vepDesp(hrec->vals[i]);
                    size_t idxFormat = vepDesp.find("Format:");
                    if(idxFormat == string::npos){
                        cout << "Cannot find the format of vep section in header of vcf file." << endl;
                        bcf_hdr_destroy(hdr);
                        hts_close(fp);
                        return -1;
                    }
                    vector<string> vsVep = split(vepDesp.substr(idxFormat + 8), "|");
                    for(size_t j=0; j<cmLine["-vep"].size(); j++){
                        vector<string>::iterator it = find(vsVep.begin(), vsVep.end(), cmLine["-vep"][j]);
                        if(it == vsVep.end()){
                            cout << cmLine["-vep"][j] << " is not included in the file. Please recheck header of vcf file." << endl;
                            bcf_hdr_destroy(hdr);
                            hts_close(fp);
                            return -1;
                        } else {
                            idxVepTag.push_back(std::distance(vsVep.begin(), it));
                        }
                    }
                    break;
                }
                //cout << hrec->keys[i] << " = " << hrec->vals [i] << endl;
            }
            //hrec is a pointer to somewhere in hdr, no need to free now
            //bcf_hrec_destroy(hrec);
        }
    }
    
    // vep filter index
    vector<size_t> idxFilterVEPNotEqual;
    if(filterVEP){
        if(idxVep < 0){
            cout << "vep information is not included in the file. Please recheck header of vcf file." << endl;
            bcf_hdr_destroy(hdr);
            hts_close(fp);
            return -1;
        }
        
        bcf_hrec_t *hrec = bcf_hdr_get_hrec(hdr, BCF_HL_INFO, "ID", "vep", NULL);
        for(int i=0; i<hrec->nkeys; i++){
            if(string(hrec->keys[i]) == "Description"){
                string vepDesp(hrec->vals[i]);
                size_t idxFormat = vepDesp.find("Format:");
                if(idxFormat == string::npos){
                    cout << "Cannot find the format of vep section in header of vcf file." << endl;
                    bcf_hdr_destroy(hdr);
                    hts_close(fp);
                    return -1;
                }
                vector<string> vsVep = split(vepDesp.substr(idxFormat + 8), "|");
                
                if(filterVEPNotEqual.size() > 0){
                    for(std::map<std::string, std::vector<std::string> >::iterator it = filterVEPNotEqual.begin(); it!=filterVEPNotEqual.end(); it++){
                        vector<string>::iterator itTmp = find(vsVep.begin(), vsVep.end(), it->first);
                        if(itTmp == vsVep.end()){
                            cout << "Filter ID " << it->first << " is not included in the VEP part of the file. Please recheck header of vcf file.";
                            bcf_hdr_destroy(hdr);
                            hts_close(fp);
                            return -1;
                        }
                        idxFilterVEPNotEqual.push_back(std::distance(vsVep.begin(), itTmp));
                    }
                }
                break;
            }
            //cout << hrec->keys[i] << " = " << hrec->vals [i] << endl;
        }
        cout << "filterVEP" << endl;
    }
    
    
    bcf1_t *rec = bcf_init1();
    if (!rec) error("bcf_init1 : %s", strerror(errno));
    
    ofstream fout(cmLine["-o"][0].c_str());
    if(!fout.is_open()){
        cout << "Cannot open file " << cmLine["-o"][0] << " for output" << endl;
        return -1;
    }
    
    //output header
    fout << "Chr\tPos\tID\tref\talt";
    for(size_t i=0; i<cmLine["-info"].size(); i++){
        fout << "\t" << cmLine["-info"][i];
    }
    for(size_t i=0; i<cmLine["-vep"].size(); i++){
        fout << "\t" << cmLine["-vep"][i];
    }
    fout << endl;
    
    
    
    int r;
    while((r = bcf_read1(fp, hdr, rec)) >=0){
        /*
        float *values = 0;
        char *dst = NULL;
        int ct = 0;
        int rt = 0;
         */
        
        bcf_info_t *info;
        bcf_unpack(rec, BCF_UN_STR);
        
        //filter
        if(filter){
            //cout << "start filter" << endl;
            bool include = true;
            int idxTmp = 0;
            if(filterLarger.size() > 0){
                for(map<string, double>::iterator it = filterLarger.begin(); it!=filterLarger.end(); it++){
                    info = bcf_get_info_id(rec, idxFilterLarger[idxTmp]);
                    idxTmp++;
                    
                    if(!info){
                        continue;
                    }
                    
                    double valTmp;
                    switch (info->type) {
                        case BCF_BT_INT8:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_INT16:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_INT32:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_INT64:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_FLOAT:
                            valTmp = info->v1.f;
                            break;
                        case BCF_BT_CHAR:
                        {
                            cout << it->first << " is a string instead of a number. Please recheck vcf file." << endl;
                            return  -1;
                        }
                        default:
                            cout << "Unrecognized info type!" << endl;
                            return  -1;
                    }
                    
                    if(valTmp > it->second){
                        include = false;
                        break;
                    }
                }
            }
            
            if(!include){
                continue;
            }
            
            
            
            
            idxTmp = 0;
            if(filterSmaller.size() > 0){
                for(map<string, double>::iterator it = filterSmaller.begin(); it!=filterSmaller.end(); it++){
                    info = bcf_get_info_id(rec, idxFilterSmaller[idxTmp]);
                    idxTmp++;
                    
                    if(!info){
                        continue;
                    }
                    
                    double valTmp;
                    switch (info->type) {
                        case BCF_BT_INT8:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_INT16:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_INT32:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_INT64:
                            valTmp = info->v1.i;
                            break;
                        case BCF_BT_FLOAT:
                            valTmp = info->v1.f;
                            break;
                        case BCF_BT_CHAR:
                        {
                            cout << it->first << " is a string instead of a number. Please recheck vcf file." << endl;
                            return  -1;
                        }
                        default:
                            cout << "Unrecognized info type!" << endl;
                            return  -1;
                    }
                    
                    if(valTmp < it->second){
                        include = false;
                        break;
                    }
                }
            }
            
            if(!include){
                continue;
            }
            
            
            
            idxTmp = 0;
            if(filterNotEqual.size() > 0){
                for(map<string, vector<string> > ::iterator it = filterNotEqual.begin(); it!=filterNotEqual.end(); it++){
                    info = bcf_get_info_id(rec, idxFilterNotEqual[idxTmp]);
                    idxTmp++;
                    
                    if(!info){
                        continue;
                    }
                    
                    string valTmp;
                    switch (info->type) {
                        case BCF_BT_INT8:
                            valTmp = to_string(info->v1.i);
                            break;
                        case BCF_BT_INT16:
                            valTmp = to_string(info->v1.i);
                            break;
                        case BCF_BT_INT32:
                            valTmp = to_string(info->v1.i);
                            break;
                        case BCF_BT_INT64:
                            valTmp = to_string(info->v1.i);
                            break;
                        case BCF_BT_FLOAT:
                            valTmp = to_string(info->v1.f);
                            break;
                        case BCF_BT_CHAR:
                        {
                            valTmp = string((char*) info->vptr, info->len);
                            break;
                        }
                        default:
                            cout << "Unrecognized info type!" << endl;
                            return  -1;
                    }
                    
                    if(find(it->second.begin(), it->second.end(), valTmp) == it->second.end()){
                        include = false;
                        break;
                    }
                }
            }
            
            if(!include){
                continue;
            }
        }
        
        //VEP filter
        if(filterVEP){
            //cout << "start vep filter" << endl;
            
            info = bcf_get_info_id(rec, idxVep);
            
            string infoTmp((char*) info->vptr, info->len);
            vector<string> vsTmp = split(infoTmp, "|");
            
            bool include = true;
            int idxTmp = 0;
            if(filterVEPNotEqual.size() > 0){
                for(map<string, vector<string> > ::iterator it = filterVEPNotEqual.begin(); it!=filterVEPNotEqual.end(); it++){
                    
                    string valTmp = vsTmp[idxFilterVEPNotEqual[idxTmp]];
                    idxTmp++;
                    
                    if(find(it->second.begin(), it->second.end(), valTmp) == it->second.end()){
                        include = false;
                        break;
                    }
                }
            }
            if(!include){
                continue;
            }
        }
        
        // basic variation information
        string chr = bcf_hdr_id2name(hdr, rec->rid);
        long long position =  rec->pos + 1;
        string rsID = rec->d.id;
        string ref = rec->d.allele[0];
        string alt = ".";
        if(rec->d.m_allele > 1){
            alt = rec->d.allele[1];
        }
        
        
        fout << chr << "\t" << position << "\t" << rsID << "\t" << ref << "\t" << alt;
        
        for(size_t i=0; i<idxInfo.size(); i++){
            // a pointer to some place in rec, no need to free
            info = bcf_get_info_id(rec, idxInfo[i]);
            
            if(!info){
                fout << "\t" << "NA";
                continue;
            }
            
            switch (info->type) {
                case BCF_BT_INT8:
                    fout << "\t" << info->v1.i;
                    break;
                case BCF_BT_INT16:
                    fout << "\t" << info->v1.i;
                    break;
                case BCF_BT_INT32:
                    fout << "\t" << info->v1.i;
                    break;
                case BCF_BT_INT64:
                    fout << "\t" << info->v1.i;
                    break;
                case BCF_BT_FLOAT:
                    fout << "\t" << info->v1.f;
                    break;
                case BCF_BT_CHAR:
                {
                    string infoTmp((char*) info->vptr, info->len);
                    fout << "\t" << infoTmp;
                    break;
                }
                default:
                    cout << "Unrecognized info type!" << endl;
                    return  -1;
            }
        }

        if(idxVep >=0){
            //BCF_BT_CHAR
            info = bcf_get_info_id(rec, idxVep);
            
            string infoTmp((char*) info->vptr, info->len);
            vector<string> vsTmp = split(infoTmp, "|");
            for(size_t i=0; i<idxVepTag.size(); i++){
                fout << "\t" << vsTmp[idxVepTag[i]];
            }
        }
        fout << endl;
        
        
        
        /*
         * another method to get infomation from vcf.bgz file
        rt = bcf_get_info_float(hdr, rec, "controls_AF", &values, &ct);
        if(rt >= 0){
            cout << rt << endl;
            cout << ct << endl;
            cout << values[0] << endl;
        }
        
        
        rt = bcf_get_info_string(hdr, rec, "vep", &dst, &ct);
        if ( rt >= 0 ) {
            vector<string> vsVep = split(dst, "|");
        } else {
            
        }
        
        free(dst);
        free(values);
         */
        
        //kstring_t str = {0,0,0};
        //vcf_format(hdr, rec, &str);
        //bcf_unpack(rec, BCF_UN_ALL);
        //string tmp = str.s;
        //cout << tmp << endl;
    }
    
    if (r < -1) error("bcf_read1");
    
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname.c_str(),ret);
        exit(ret);
    }
    
    cout << "Finished processing!" << endl;
    
    return 0;
}
