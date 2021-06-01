#include <stdio.h>
#include <stdint.h>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <getopt.h>
#include <math.h>
#include <time.h>
#include <list>
#include <sstream>

using namespace std;

typedef struct {
  long start;
  long end;
} interval;

void read_data(string filename,vector<interval> &data,vector<string> &names,vector<int> &cn);

int main(int argc, char *argv[])
{
  if(argc!=2)
    {
      cerr << "Usage: <seg-file>\n";
      exit(1);
    }

  string filename=argv[1];
  vector<interval> seg,seg_avg;
  interval tmp_seg;
  vector<double> cn_avg;
  vector<string> ring_names;
  vector<int> cn;
  long glength=133851895;
  vector<double> avg_profile(glength+1,0);
  long n_rings=1;
  
  read_data(filename,seg,ring_names,cn);
  for(long i=2;i<ring_names.size();++i)
    {
      if(ring_names[i]!=ring_names[i-1])
	++n_rings;
    }
  for(long i=0;i<seg.size();++i)
    {
      for(long h=seg[i].start;h<=seg[i].end;++h)
	avg_profile[h]+=cn[i];
    }
  
  for(long i=0;i<glength+1;++i)
    avg_profile[i]*=1.0/(double)n_rings;

  tmp_seg.start=1;
  for(long i=1;i<glength;++i)
    {
      if(avg_profile[i]!=avg_profile[i+1])
	{
	  tmp_seg.end=i;
	  seg_avg.push_back(tmp_seg);
	  cn_avg.push_back(avg_profile[i]);
	  tmp_seg.start=i+1;
	}
    }
  tmp_seg.end=glength;
  seg_avg.push_back(tmp_seg);
  cn_avg.push_back(avg_profile[glength]);

  cout << "Sample\tChromosome\tStart\tEnd\tnMarker\tCN\n";
  for(long i=0;i<seg_avg.size();++i)
    {
      cout << filename.substr(0,filename.find(".seg")) << "\tchr12\t" << seg_avg[i].start << "\t" << seg_avg[i].end << "\t" << seg_avg[i].end-seg_avg[i].start+1 << "\t" << cn_avg[i] << endl;
    }
  
}

void read_data(string filename,vector<interval> &data,vector<string> &names,vector<int> &cn)
{
  ifstream in;
  string tmp;
  stringstream line;
  interval tmp_data;

  in.open(filename.c_str());
  if(!in.is_open())
    {
      cerr << "Cannot open: " << filename << endl;
      exit(1);
    }

  data.clear(); names.clear(); cn.clear();

  getline(in,tmp,'\n');
  while(getline(in,tmp,'\n'))
    {
      line.str(""); line.clear();
      line << tmp;

      getline(line,tmp,'\t');
      names.push_back(tmp);
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      tmp_data.start=atol(tmp.c_str());
      getline(line,tmp,'\t');
      tmp_data.end=atol(tmp.c_str());
      data.push_back(tmp_data);
      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      cn.push_back(atoi(tmp.c_str()));
    }
  in.close();
}
