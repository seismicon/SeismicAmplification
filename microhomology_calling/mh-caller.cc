#include <stdio.h>
#include <stdint.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <unordered_map>
#include <sstream>

#include "nibtools.h"
#include "installdir.h"

using namespace std;



string mh_caller_help = "\n\t --- mh-caller: Calling Micro-Homologies from Rearrangment Data ---\n \
      \t\t    Version 0.1  10/December/2019\n \
      \t     University of Cologne \n \
      \t     Author: Peifer Lab\n\n\
      SYNOPSIS \n \t mh-caller  <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -pl       \t rearrangement data file from peiflyne\n \
     \t -icgc     \t rearrangement data file from ICGC PCAWG project\n \
     \t -o        \t sample name\n\
     \t -std      \t write output to standard out (flag)\n \
     \t -fs       \t flanksize [10]\n\
     \t -md       \t maximal breakpoint distance [2]\n\
     \t -ms       \t minimal microhomology size [3]\n\
     \t -build    \t genome build (possible: hg19) [hg19]\n\n";


typedef struct
{
  string chr1;
  int chr1_file;
  long pos1;
  string chr2;
  int chr2_file;
  long pos2;
  long nreads;
  char strand1;
  char strand2;
  char type;
  int hom_length;
  string hom_seq;
  int dist_from_bkp1;
  int dist_from_bkp2;
} rearrangement_data;

void init_ref_genome(nib *ref,string build,unordered_map<string,int> chrmap,int nchr);
void read_icgc(string icgc_file,unordered_map<string,int> chrmap,vector<rearrangement_data> &rearr);
void output_stdout(string sample,vector<rearrangement_data> rearr);
void output_file(string sample,vector<rearrangement_data> rearr);
void call_mh(nib *ref,int flanksize,int max_bkp_dist,int min_mh_size,rearrangement_data &rearr);
void read_peiflyne(string peiflyne_file,unordered_map<string,int> chrmap,vector<rearrangement_data> &rearr);

int main(int argc, char *argv[])
{
  int longindex,opt;
  // option definition

  static struct option longopts[]={
    {"help"     , 0, 0,  'h'},
    {"help"     , 0, 0,  '?'},
    {"pl"       , 1, 0,    1},
    {"icgc"     , 1, 0,    2},
    {"o"        , 1, 0,    3},
    {"build"    , 1, 0,    4},
    {"fs"       , 1, 0,    5},
    {"md"       , 1, 0,    6},
    {"ms"       , 1, 0,    7},
    {"std"      , 0, 0,    8},
    {0, 0, 0, 0}
  }; 

  string pl_file="";
  string icgc_file="";
  string build="hg19";
  string sample="";
  int flanksize=10;
  int max_bkp_dist=2;
  int min_mh_size=3;
  bool std=0;

  optind=0;
  //parse command line arguments
  while((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
	{
        case 'h':
	  cerr << mh_caller_help;
	  exit(1);
	  break;
      	case '?':
	  cerr << mh_caller_help;
	  exit(1);  
	  break;
        case 1:
	  pl_file=(string)optarg;
          break;
        case 2:
	  icgc_file=(string)optarg;
          break;
        case 3:
	  sample=(string)optarg;
          break;
	case 4:
	  build = (string)optarg;
	  if(build != "hg19")
	    {
	      cerr << "Error: please use hg19. Falling back to default value: hg19.\n"; 
	      build = "hg19";
	    }
	  break;
	case 5:
	  flanksize=atoi(optarg);
	  break;
	case 6:
	  max_bkp_dist=atoi(optarg);
	  break;
	case 7:
	  min_mh_size=atoi(optarg);
	  break;
	case 8:
	  std=1;
	  break;
	default:
          cerr << "Error: cannot parse arguments.\n";
          exit(1);
	  break;
	}
    }

  if((pl_file=="" && icgc_file=="") || (pl_file!="" && icgc_file!=""))
    {
      cerr << "Error: please specifiy either and icgc pcawg or a peiflyne rearrangement file.\n";
      cerr << mh_caller_help;
      exit(1);
    }

  if(sample=="")
    {
      cerr << "Error: please provide sample name.\n";
      cerr << mh_caller_help;
      exit(1);
    }

  unordered_map<string,int> chrmap;
  int nchr;
  nib *ref;
  vector<rearrangement_data> rearr;

  if(build=="hg19")
    {
      chrmap["chr1"]=0;   chrmap["chr2"]=1;   chrmap["chr3"]=2;   chrmap["chr4"]=3;   chrmap["chr5"]=4;
      chrmap["chr6"]=5;   chrmap["chr7"]=6;   chrmap["chr8"]=7;   chrmap["chr9"]=8;   chrmap["chr10"]=9;
      chrmap["chr11"]=10; chrmap["chr12"]=11; chrmap["chr13"]=12; chrmap["chr14"]=13; chrmap["chr15"]=14;
      chrmap["chr16"]=15; chrmap["chr17"]=16; chrmap["chr18"]=17; chrmap["chr19"]=18; chrmap["chr20"]=19;
      chrmap["chr21"]=20; chrmap["chr22"]=21; chrmap["chrX"]=22;  chrmap["chrY"]=23;
      nchr=24;
    }
  
  ref = new nib[nchr];
  
  init_ref_genome(ref,build,chrmap,nchr);
  
  if(icgc_file!="")
    read_icgc(icgc_file,chrmap,rearr);
  if(pl_file!="")
    read_peiflyne(pl_file,chrmap,rearr);
 
  for(long i=0;i<rearr.size();++i)
    call_mh(ref,flanksize,max_bkp_dist,min_mh_size,rearr[i]);

  delete[] ref;
  
  if(std==1)
    output_stdout(sample,rearr);
  else
    output_file(sample,rearr);
}

char complement(char base)
{
  if(base=='A')
    return('T');
  else if(base=='C')
    return('G');
  else if(base=='G')
    return('C');
  else if(base=='T')
    return('A');
  else
    return('N');
}

string rev_compl(string seq)
{
  string seq_rc="";
  

  for(long i=seq.size()-1;i>=0;--i)
    seq_rc+=complement(seq[i]);
  
  return(seq_rc);
    
}

bool comp_base(char base1,char base2)
{
  if(base1=='N' || base2=='N')
    return(0);
  
  if(base1==base2)
    return(1);
  else
    return(0);
}

void get_contig_matrix(int flanksize,int max_bkp_dist,int min_mh_size,string seq1,string seq2,vector< vector<int> > &contig)
{
  vector< vector<bool> > match(2*max_bkp_dist+1,vector<bool>(flanksize+1,0));
  vector< vector<int> > lm(2*max_bkp_dist+1,vector<int>(flanksize+1,0));
  vector< vector<int> > rm(2*max_bkp_dist+1,vector<int>(flanksize+1,0));

  for(int h=-max_bkp_dist;h<=max_bkp_dist;++h)
    {
      for(int i=0;i<=flanksize;++i)
	{
	  match[max_bkp_dist+h][i]=comp_base(seq1[i],seq2[i+max_bkp_dist+h]);
	}

      lm[max_bkp_dist+h][0]=(int)match[max_bkp_dist+h][0];
      rm[max_bkp_dist+h][flanksize]=(int)match[max_bkp_dist+h][flanksize];
      for(int i=1;i<=flanksize;++i)
	{
	  if(match[max_bkp_dist+h][i]!=0)
	    lm[max_bkp_dist+h][i]=lm[max_bkp_dist+h][i-1]+1;
	  else if(i==flanksize+1)
	    lm[max_bkp_dist+h][i]=(int)match[max_bkp_dist+h][i];
	  else
	    lm[max_bkp_dist+h][i]=0;

	  if(match[max_bkp_dist+h][flanksize-i]!=0)
	    rm[max_bkp_dist+h][flanksize-i]=rm[max_bkp_dist+h][flanksize-i+1]+1;
	  else if(i==flanksize)
	    rm[max_bkp_dist+h][flanksize-i]=(int)match[max_bkp_dist+h][flanksize-i];
	  else
	    rm[max_bkp_dist+h][flanksize-i]=0;
	}

      for(int i=0;i<=flanksize;++i)
	{
	  if(lm[max_bkp_dist+h][i]+rm[max_bkp_dist+h][i]-1 >= min_mh_size)
	    contig[max_bkp_dist+h][i]=lm[max_bkp_dist+h][i]+rm[max_bkp_dist+h][i]-1;
	  else
	    contig[max_bkp_dist+h][i]=0;
	}
    }

}

typedef struct{
  int x;
  int y;
  int length;
} contig_data;

void call_mh(nib *ref,int flanksize,int max_bkp_dist,int min_mh_size,rearrangement_data &rearr)
{

  string seq1_l,seq2_l,seq1_r,seq2_r,tmp_l="",tmp_r="";
  vector< vector<int> > contig_l(2*max_bkp_dist+1,vector<int>(flanksize+1,0));
  vector< vector<int> > contig_r(2*max_bkp_dist+1,vector<int>(flanksize+1,0));
  contig_data max_contig_l,max_contig_r;

  max_contig_l.x=-1; max_contig_l.y=-1; max_contig_l.length=0;
  max_contig_r.x=-1; max_contig_r.y=-1; max_contig_r.length=0;

  ref[rearr.chr1_file].getSeq(seq1_l,rearr.pos1-flanksize,rearr.pos1+1);
  ref[rearr.chr2_file].getSeq(seq2_l,rearr.pos2-flanksize-max_bkp_dist,rearr.pos2+max_bkp_dist+1);

  ref[rearr.chr1_file].getSeq(seq1_r,rearr.pos1,rearr.pos1+flanksize+1);
  ref[rearr.chr2_file].getSeq(seq2_r,rearr.pos2-max_bkp_dist,rearr.pos2+flanksize+max_bkp_dist+1);

  if(rearr.type=='h' || rearr.type=='t')
    {
      tmp_r=rev_compl(seq2_l);
      tmp_l=rev_compl(seq2_r);
      seq2_l=tmp_l;
      seq2_r=tmp_r;
    }

  get_contig_matrix(flanksize,max_bkp_dist,min_mh_size,seq1_l,seq2_l,contig_l);
  get_contig_matrix(flanksize,max_bkp_dist,min_mh_size,seq1_r,seq2_r,contig_r);

  //find max contig in bounds
  //left matrix
  for(int y=flanksize;y>=flanksize-max_bkp_dist;--y)
    {
      for(int x=flanksize-y;x<=max_bkp_dist;++x)
	{
	  
	  if(contig_l[x][y] > max_contig_l.length)
	    {
	      max_contig_l.x=x;
	      max_contig_l.y=y;
	      max_contig_l.length=contig_l[x][y];
	    }
	}
    }
  //right matrix
  for(int y=0;y<=max_bkp_dist;++y)
    {
      for(int x=max_bkp_dist;x<=2*max_bkp_dist-y;++x)
	{
	  if(contig_r[x][y] > max_contig_r.length)
	    {
	      max_contig_r.x=x;
	      max_contig_r.y=y;
	      max_contig_r.length=contig_r[x][y];
	    }
	}
    }

  if(max_contig_l.length==0 && max_contig_r.length==0)
    {
      rearr.hom_length=0;
      rearr.hom_seq=".";
      rearr.dist_from_bkp1=-1;
      rearr.dist_from_bkp2=-1;
    }
  else
    {
      if(max_contig_l.length >= max_contig_r.length)
	{
	  rearr.dist_from_bkp1=(flanksize-max_contig_l.y);
	  rearr.dist_from_bkp2=abs(max_bkp_dist-max_contig_l.x);
	  rearr.hom_length=max_contig_l.length;
	  rearr.hom_seq=seq1_l.substr(max_contig_l.y-max_contig_l.length+1,max_contig_l.length);
	}
      else
	{
	  rearr.dist_from_bkp1=(max_contig_r.y);
	  rearr.dist_from_bkp2=abs(max_bkp_dist-max_contig_r.x);
	  rearr.hom_length=max_contig_r.length;
	  rearr.hom_seq=seq1_r.substr(max_contig_r.y,max_contig_r.length);
	}
    }
}

void read_icgc(string icgc_file,unordered_map<string,int> chrmap,vector<rearrangement_data> &rearr)
{
  ifstream in;
  string tmp;
  stringstream line;
  rearrangement_data rearr_tmp;
  unordered_map<string,int>::iterator it;

  rearr.clear();
  in.open(icgc_file.c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open: " << icgc_file << ".\n";
      exit(1);
    }

  getline(in,tmp,'\n');
  while(getline(in,tmp,'\n'))
    {
      rearr_tmp.hom_length=-1;
      rearr_tmp.hom_seq=".";
      rearr_tmp.chr1_file=-1;
      rearr_tmp.chr2_file=-1;
      line.str(""); line.clear();
      line << tmp;
      getline(line,tmp,'\t');
      rearr_tmp.chr1="chr"+tmp;
      getline(line,tmp,'\t');
      rearr_tmp.pos1=atol(tmp.c_str());
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.chr2="chr"+tmp;
      getline(line,tmp,'\t');
      rearr_tmp.pos2=atol(tmp.c_str());
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.nreads=atol(tmp.c_str());
      getline(line,tmp,'\t');
      rearr_tmp.strand1=tmp[0];
      getline(line,tmp,'\t');
      rearr_tmp.strand2=tmp[0];
      
      if((rearr_tmp.strand1=='+' && rearr_tmp.strand2=='-') || (rearr_tmp.strand1=='-' && rearr_tmp.strand2=='+'))
	rearr_tmp.type='n';
      else if(rearr_tmp.strand1=='-' && rearr_tmp.strand2=='-')
	rearr_tmp.type='t';
      else if(rearr_tmp.strand1=='+' && rearr_tmp.strand2=='+')
	rearr_tmp.type='h';
      else
	rearr_tmp.type='f';

      it=chrmap.find(rearr_tmp.chr1);
      if(it!=chrmap.end())
	rearr_tmp.chr1_file=it->second;
      it=chrmap.find(rearr_tmp.chr2);
      if(it!=chrmap.end())
	rearr_tmp.chr2_file=it->second;
      
      if(rearr_tmp.type!='f' && rearr_tmp.chr1_file!=-1 && rearr_tmp.chr2_file!=-1)
	rearr.push_back(rearr_tmp);
    }
  in.close();
}

void read_peiflyne(string peiflyne_file,unordered_map<string,int> chrmap,vector<rearrangement_data> &rearr)
{
  ifstream in;
  string tmp;
  stringstream line;
  rearrangement_data rearr_tmp;
  unordered_map<string,int>::iterator it;

  rearr.clear();
  in.open(peiflyne_file.c_str());
  if(!in.is_open())
    {
      cerr << "Error: cannot open: " << peiflyne_file << ".\n";
      exit(1);
    }

  getline(in,tmp,'\n');
  while(getline(in,tmp,'\n'))
    {
      rearr_tmp.hom_length=-1;
      rearr_tmp.hom_seq=".";
      rearr_tmp.chr1_file=-1;
      rearr_tmp.chr2_file=-1;
      line.str(""); line.clear();
      line << tmp;
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.chr1=tmp;
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.pos1=atol(tmp.c_str());
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.chr2=tmp;
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.pos2=atol(tmp.c_str());
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      rearr_tmp.nreads=atol(tmp.c_str());
      getline(line,tmp,'\t');
      rearr_tmp.type=tmp[0];
 
      if(rearr_tmp.type=='n')
	{
	  rearr_tmp.strand1='+'; rearr_tmp.strand2='-';
	}
      else if(rearr_tmp.type=='t')
	{
	  rearr_tmp.strand1='-'; rearr_tmp.strand2='-';
	}
      else if(rearr_tmp.type=='h')
	{
	  rearr_tmp.strand1='+'; rearr_tmp.strand2='+';
	}
      else
	rearr_tmp.type='f';	

      it=chrmap.find(rearr_tmp.chr1);
      if(it!=chrmap.end())
	rearr_tmp.chr1_file=it->second;
      it=chrmap.find(rearr_tmp.chr2);
      if(it!=chrmap.end())
	rearr_tmp.chr2_file=it->second;
      
      if(rearr_tmp.type!='f' && rearr_tmp.chr1_file!=-1 && rearr_tmp.chr2_file!=-1)
	rearr.push_back(rearr_tmp);
    }
  in.close();
}

void output_stdout(string sample,vector<rearrangement_data> rearr)
{
  cout << "sample\tchr_1\tpos_1\tchr_2\tpos_2\tnreads\tstrand_1\tstrand_2\ttype\tmh_length\tmh_seq\tdist_from_bkp1\tdist_from_bkp2\n";
  for(int i=0;i<rearr.size();++i)
    {
      cout << sample << "\t";
      cout << rearr[i].chr1 << "\t" << rearr[i].pos1 << "\t" << rearr[i].chr2 << "\t" << rearr[i].pos2 << "\t" ;
      cout << rearr[i].nreads << "\t" << rearr[i].strand1 << "\t" << rearr[i].strand2 << "\t" << rearr[i].type << "\t" ;
      cout << rearr[i].hom_length << "\t" << rearr[i].hom_seq << "\t";
      if(rearr[i].dist_from_bkp1==-1 || rearr[i].dist_from_bkp2==-1)
	cout << ".\t.\n";
      else
	cout << rearr[i].dist_from_bkp1 << "\t" << rearr[i].dist_from_bkp2 << endl;
    }
}

void output_file(string sample,vector<rearrangement_data> rearr)
{

  ofstream out;

  out.open((sample+"_rearrangemnents_mh.txt").c_str());

  out << "sample\tchr_1\tpos_1\tchr_2\tpos_2\tnreads\tstrand_1\tstrand_2\ttype\tmh_length\tmh_seq\tdist_from_bkp1\tdist_from_bkp2\n";
  for(int i=0;i<rearr.size();++i)
    {
      out << sample << "\t";
      out << rearr[i].chr1 << "\t" << rearr[i].pos1 << "\t" << rearr[i].chr2 << "\t" << rearr[i].pos2 << "\t" ;
      out << rearr[i].nreads << "\t" << rearr[i].strand1 << "\t" << rearr[i].strand2 << "\t" << rearr[i].type << "\t" ;
      out << rearr[i].hom_length << "\t" << rearr[i].hom_seq << "\t";
      if(rearr[i].dist_from_bkp1==-1 || rearr[i].dist_from_bkp2==-1)
	out << ".\t.\n";
      else
	out << rearr[i].dist_from_bkp1 << "\t" << rearr[i].dist_from_bkp2 << endl;
    }

  out.close();
}

void init_ref_genome(nib *ref,string build,unordered_map<string,int> chrmap,int nchr)
{
  stringstream line;
  
  for(unordered_map<string,int>::iterator it=chrmap.begin();it!=chrmap.end();++it)
    {
      line.str(""); line.clear();
      line << INSTALLDIR << "/" << build << "/" << build << "_" << it->first << ".nib";
      
      ref[it->second].open(line.str());
    }
}
