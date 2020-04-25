//#include <config.h>

#include <stdio.h>
#include <iostream>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/kstring.h>
#include <htslib/kseq.h>

#include "normal.h"

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

int main(int argc, char **argv){
/*
 * -i input.vcf.bgz
 * -o output.txt (tab delimited)
 * -info AF AC ...  // id in information section, please check VCF header file for details
 * -vep SYMBOL Gene ... // id in VEP section, please check VCF header for details
 */
    map<string, vector<string> > cmLine = parseCMLine(argc, argv);
    for(map<string, vector<string> >::iterator it = cmLine.begin(); it != cmLine.end(); it++){
        cout << it->first << endl;
        for(vector<string>::iterator it2 = it->second.begin(); it2 != it->second.end(); it2++){
            cout << *it2 << "\t";
        }
        cout << endl << endl;
    }
    
    /*
    char* fname = (char*) "/Users/gpeng/Work/MetaRace/src/CPP/ReadGnomad/data/gnomad.genomes.r3.0.sites.chrY.vcf.bgz";
    htsFile *fp    = hts_open(fname,"rb");
    if (!fp) error("Failed to open \"%s\" : %s", fname, strerror(errno));
    bcf_hdr_t *hdr = bcf_hdr_read(fp);
    if (!hdr) error("bcf_hdr_read : %s", strerror(errno));
    bcf1_t *rec    = bcf_init1();
    if (!rec) error("bcf_init1 : %s", strerror(errno));
     
    
    int r;
    float *values = 0;
    char *dst = NULL;
    int ct = 0;
    int rt = 0;
    while((r = bcf_read1(fp, hdr, rec)) >=0){
        bcf_unpack(rec, BCF_UN_STR);
        
        string chr = bcf_hdr_id2name(hdr, rec->rid);
        long long position =  rec->pos + 1;
        string rsID = rec->d.id;
        string ref = rec->d.allele[0];
        string alt = ".";
        if(rec->d.m_allele > 1){
            alt = rec->d.allele[1];
        }
        
        
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
           
        
        //kstring_t str = {0,0,0};
        //vcf_format(hdr, rec, &str);
        //bcf_unpack(rec, BCF_UN_ALL);
        //string tmp = str.s;
        //cout << tmp << endl;
    }
    
    if (r < -1) error("bcf_read1");
    
    
    free(dst);
    free(values);
    bcf_destroy1(rec);
    bcf_hdr_destroy(hdr);
    int ret;
    if ( (ret=hts_close(fp)) )
    {
        fprintf(stderr,"hts_close(%s): non-zero status %d\n",fname,ret);
        exit(ret);
    }
     */
    
    return 0;
}
