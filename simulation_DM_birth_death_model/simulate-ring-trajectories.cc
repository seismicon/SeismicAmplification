#include <stdio.h>
#include <stdint.h>
#include <string>
#include <string.h>
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

//definitions for ran1 --> start
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

double ran1(long *idum);
//end <--

using namespace std;

typedef struct {
  long start;
  long end;
} interval;

typedef struct {
  double fitness;
  long size;
  long cycle_tot;
  long cycle_even;
  long cycle_odd;
  vector<interval> c_ring;
  vector<interval> seg;
  vector<int> cn_seg;
  vector<int> cn_regions;
} ring_data;

typedef struct {
  long glength;
  long max_size;
  double br;
  double br0;
  double dr;
  double N_max;
  double K;
  long N;
  long N_burnin;
  double dt;
  long N_fitness;
  double R_fitness;
  double p_cross_even;
  double dt_out;
  vector<double> fitness_function;
  vector<interval> sel_regions;
  vector<string> sel_regions_names;
  vector<interval> regions_amp;
  string outname;
  string seg_folder;
} sim_par;

int sample(long *idum,sim_par par,double fitness);
void read_data(string filename,vector<interval> &data,vector<string> &names,bool extract_names);
void calc_cn(sim_par par,ring_data &ring,vector<int> &cn);
void ring_size(ring_data &ring);
void crossover_odd(long *idum,sim_par par,ring_data ring,ring_data &ring_1);
void crossover_even(long *idum,sim_par par,ring_data ring,ring_data &ring_1,ring_data &ring_2);
double calc_fitness(sim_par par,ring_data ring);
void output_stats(double t,sim_par &par, vector<ring_data> ring); 

string help = "\n      SYNOPSIS \n \t simulate-ring-trajectories  <options>\n\n \
     DESCRIPTION\n \
     \t -h -? -help \t help\n \
     \t -o        \t output name (required)\n \
     \t -c        \t chromothripsis file (required)\n \
     \t -l        \t region list (1,2,3) [1]\n \
     \t -br       \t birth rate (maximal) [1]\n \
     \t -br0      \t neutral birth rate [0.5]\n \
     \t -dr       \t death rate [0.1]\n \
     \t -pe       \t probablilty of even crossover [0.9]\n \
     \t -maxS     \t maximal ring size in Mbp [25]\n \
     \t -maxN     \t maximal average ring number [1000]\n \
     \t -tmax     \t maximal time of simulation [200]\n \
   \n";

int main(int argc, char *argv[])
{
  int longindex,opt;
  // option definition

  static struct option longopts[]={
    {"help"     , 0, 0,  'h'},
    {"help"     , 0, 0,  '?'},
    {"o"        , 1, 0,    1},
    {"c"        , 1, 0,    2},
    {"l"        , 1, 0,    3},
    {"br"       , 1, 0,    4},
    {"br0"      , 1, 0,    5},
    {"dr"       , 1, 0,    6},
    {"pe"       , 1, 0,    7},
    {"maxS"     , 1, 0,    8},
    {"maxN"     , 1, 0,    9},
    {"tmax"     , 1, 0,    10},
    {0, 0, 0, 0}
  }; 

  sim_par par;

  par.br=1; //bith rate
  par.dr=0.1; //death rate
  par.br0=0.5; //zero birth rate 
  par.N_max=1000;
  par.K=par.N_max/(1.0-par.dr/par.br);
  par.N_burnin=8;
  par.dt=0.01;
  par.glength=133851895;
  par.N_fitness=16;
  par.R_fitness=4;
  par.max_size=25000000;
  par.p_cross_even=0.9;
  par.outname="";
  par.dt_out=1;

  double t_max=200;
  double t=0,t_out=1;
  string file_chromothripsis="";
  string list="1";

  optind=0;
  //parse command line arguments
  while((opt=getopt_long_only(argc,argv,"h?",longopts,&longindex)) != -1)
    {
      switch(opt)
	{
        case 'h':
	  cerr << help;
	  cerr << endl;
	  exit(1);
	  break;
      	case '?':
	  cerr << help;
	  cerr << endl;
	  exit(1);
	  break;
        case 1:
	  par.outname=(string)optarg;
          break;
	case 2:
	  file_chromothripsis=(string)optarg;
          break;
	case 3:
	  if(strncmp(optarg,"1",1)==0 || strncmp(optarg,"2",1)==0 || strncmp(optarg,"3",1)==0)
	    {
	      list=(string)optarg;
	    }
	  else
	    {
	      cerr << "Error: invalid list number\n";
	      exit(1);
	    }
	  break;
	case 4:
	  par.br=atof(optarg);
	  break;
	case 5:
	  par.br0=atof(optarg);
	  break;
	case 6:
	  par.dr=atof(optarg);
	  break;
	case 7:
	  par.p_cross_even=atof(optarg);
	  break;
	case 8:
	  par.max_size=atoi(optarg)*1000000;
	  break;
	case 9:
	  par.N_max=atol(optarg);
	  break;
	case 10:
	  t_max=atof(optarg);
	  break;
	default:
          cerr << "Error: cannot parse arguments.\n";
          exit(1);
	  break;
	}
    }

  if((argc-optind)!=0)
    {
      cerr << help;
      cerr << endl;
      exit(1);
    }

  if(par.outname=="")
    {
      cerr << "Error: please specify output name.\n";
      cerr << help;
      cerr << endl;
      exit(1);
    }

  if(file_chromothripsis=="")
    {
      cerr << "Error: please specify chromothripsis file.\n";
      cerr << help;
      cerr << endl;
      exit(1);
    }

  par.K=par.N_max/(1.0-par.dr/par.br);
  par.seg_folder=par.outname+"_copy_number_seg_files";
  string file_sel_regions="list-"+list+"-selected_regions.txt";
  string file_regions_amp="list-"+list+"-regions_amp.txt";

  long idum=-time(NULL);
  ofstream out;
  ring_data tmp_ring;
  ring_data tmp_ring_1;
  ring_data tmp_ring_2;
  vector<ring_data> ring;
  vector<ring_data> ring_new;
  vector<string> tmp_names;
  vector<interval> tmp_data;
  vector<int> cn;
  int action;
  
  system(("rm -rf "+ par.seg_folder).c_str());
  system(("mkdir "+ par.seg_folder).c_str());

  out.open((par.outname+"_statistics.txt").c_str());
  out << "time\tN\tavg_fitness\tavg_fitness_not_zero\tN_not_zero\tavg_ring_size[Mb]\tavg_cycles_tot\tavg_cycles_odd\tavg_cycles_even\n";
  out.close();

  //compute fitness function
  par.fitness_function=vector<double> (par.N_fitness+1,0);
  par.fitness_function[1]=(par.R_fitness*(double)par.N_fitness+1.0)/((double)par.N_fitness-1.0);
  for(long j=2;j<=par.N_fitness;++j)
    par.fitness_function[j]=par.fitness_function[j-1]*((par.R_fitness*(double)par.N_fitness+1.0)-(par.R_fitness-1)*j)/((double)par.N_fitness-1.0);
  for(long j=0;j<=par.N_fitness;++j)
    {
      par.fitness_function[j]*=1.0/par.fitness_function[par.N_fitness];
    }
  //end compute fitness function
  
  read_data(file_regions_amp,tmp_data,tmp_names,0);
  par.regions_amp=tmp_data;
  read_data(file_sel_regions,tmp_data,tmp_names,1);
  par.sel_regions=tmp_data; par.sel_regions_names=tmp_names;
  read_data(file_chromothripsis,tmp_data,tmp_names,0);
  tmp_ring.c_ring=tmp_data;
  tmp_ring.fitness=par.br0/par.br;
  tmp_ring.cycle_tot=0;
  tmp_ring.cycle_even=0;
  tmp_ring.cycle_odd=0;
  calc_cn(par,tmp_ring,cn);
  ring_size(tmp_ring);
  ring.push_back(tmp_ring);

  for(t=0;t<=t_max;t=t+par.dt)
    {
      par.N=ring.size();
      ring_new.clear();
      for(long i=0;i<par.N;++i)
	{
	  action=sample(&idum,par,ring[i].fitness);
	  if(action==0)
	    ring_new.push_back(ring[i]);
	  if(action==1)
	    {
	      if(ran1(&idum) < par.p_cross_even)
		{
		  crossover_even(&idum,par,ring[i],tmp_ring_1,tmp_ring_2);
		  ring_new.push_back(tmp_ring_1);
		  ring_new.push_back(tmp_ring_2);
		}
	      else
		{
		  crossover_odd(&idum,par,ring[i],tmp_ring_1);
		  ring_new.push_back(tmp_ring_1);
		}
	    }
	}
      ring=ring_new;
      par.N=ring.size();
      //output statistics
      if(t_out<=t+0.001)
	{
	  output_stats(t,par,ring);
	  t_out+=par.dt_out;
	}
    }
  

  return(0);
}

int sample(long *idum,sim_par par,double fitness)
{
  //output: -1 death, 1 proliferate, 0 do nothing
  int out;
  double sample;
  double pd=par.dr*par.dt;
  double pb=par.br*fitness*par.dt*(1.0-(double)par.N/par.K);

  if(pb<0)
    pb=0;

  if(par.N <= par.N_burnin)
    pd=0;
  
  sample=ran1(idum);

  if(sample < pb)
    out=1;
  else if(sample >= pb && sample < pb+pd)
    out=-1;
  else
    out=0;

  return(out);
}

double ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  double temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--) {
      k=(*idum)/IQ;
      *idum=IA*(*idum-k*IQ)-IR*k;
      if (*idum < 0) *idum += IM;
      if (j < NTAB) iv[j] = *idum;
    }
    iy=iv[0];
  }
  k=(*idum)/IQ;
  *idum=IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0) *idum += IM;
  j=iy/NDIV;
  iy=iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}

void read_data(string filename,vector<interval> &data,vector<string> &names,bool extract_names)
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

  data.clear(); names.clear();

  getline(in,tmp,'\n');
  while(getline(in,tmp,'\n'))
    {
      line.str(""); line.clear();
      line << tmp;

      getline(line,tmp,'\t');
      getline(line,tmp,'\t');
      tmp_data.start=atol(tmp.c_str());
      getline(line,tmp,'\t');
      tmp_data.end=atol(tmp.c_str());
      data.push_back(tmp_data);
      if(extract_names==1)
	{
	   getline(line,tmp,'\t');
	   getline(line,tmp,'\t');
	   getline(line,tmp,'\t');
	   names.push_back(tmp);
	}
    }
  in.close();
}

void calc_cn(sim_par par,ring_data &ring,vector<int> &cn)
{
  cn=vector<int>(par.glength+1,1);
  interval tmp_seg;
  ring.cn_regions=vector<int>(par.sel_regions.size(),1);

  for(long i=0;i<par.regions_amp.size();++i)
    {
      for(long h=par.regions_amp[i].start;h<=par.regions_amp[i].end;++h)
	cn[h]=0;
    }

  for(long i=0;i<ring.c_ring.size();++i)
    {
      for(long h=ring.c_ring[i].start;h<=ring.c_ring[i].end;++h)
	{
	  ++cn[h];
	}
    }
  
  for(long i=0;i<par.sel_regions.size();++i)
    {
      for(long h=par.sel_regions[i].start;h<=par.sel_regions[i].end;++h)
	{
	  if(ring.cn_regions[i] < cn[h])
	    ring.cn_regions[i]=cn[h];
	}
    }
  //calculate segments
  ring.seg.clear();
  ring.cn_seg.clear();
  tmp_seg.start=1;
  for(long i=1;i<par.glength;++i)
    {
      if(cn[i]!=cn[i+1])
	{
	  tmp_seg.end=i;
	  ring.seg.push_back(tmp_seg);
	  ring.cn_seg.push_back(cn[i]);
	  tmp_seg.start=i+1;
	}
    }
  tmp_seg.end=par.glength;
  ring.seg.push_back(tmp_seg);
  ring.cn_seg.push_back(cn[par.glength]);

}

void ring_size(ring_data &ring)
{
  ring.size=0;
  for(long i=0;i<ring.c_ring.size();++i)
    ring.size+=(ring.c_ring[i].end-ring.c_ring[i].start+1);
}

double calc_fitness(sim_par par,ring_data ring)
{
  double F,f=1;
   if(ring.size > par.max_size)
     return(0);
   else
     {
       for(long i=0;i<ring.cn_regions.size();++i)
	 {
	   if(ring.cn_regions[i] < par.N_fitness)
	     f=f*par.fitness_function[ring.cn_regions[i]];
	 }
       F=(1-par.br0/par.br)*f+par.br0/par.br;
       return(F);
     }

}

void get_bkp(long *idum,vector<interval> c_ring,long &idx,long &bp)
{
  vector<long> usable_idx;
  vector<double> p;
  double rn,seg_size,sum;

  seg_size=0;
  for(long i=0;i<c_ring.size();++i)
    {
      if(c_ring[i].end-c_ring[i].start+1 > 3)//only segments with more than 3bp are usable
	{
	  usable_idx.push_back(i);
	  p.push_back(c_ring[i].end-c_ring[i].start+1);
	  seg_size+=c_ring[i].end-c_ring[i].start+1;
	}
    }

  sum=0;
  rn=ran1(idum);
  for(long i=0;i<p.size();++i)
    {
      p[i]=p[i]/seg_size;
      sum+=p[i];
      if(sum>rn)
	{
	  idx=i;
	  break;
	}
    }
  
  idx=usable_idx[idx];
  rn=ran1(idum);
  seg_size=c_ring[idx].end-c_ring[idx].start-1;
  bp=c_ring[idx].start+(long)(rn*seg_size);
}

void add_bkp_recenter(long idx,long bp,vector<interval> &c_ring)
{
  interval tmp_seg;
  vector<interval> tmp_c_ring;
  long pos;

  tmp_seg.end=c_ring[idx].end;
  c_ring[idx].end=bp;
  tmp_seg.start=bp+1;
  c_ring.insert(c_ring.begin()+idx+1,tmp_seg);
  //recenter
  pos=idx+1;
  for(long i=0;i<c_ring.size();++i)
    {
      tmp_c_ring.push_back(c_ring[pos]);
      ++pos;
      if(pos>=c_ring.size())
	pos=0;
    }
  c_ring=tmp_c_ring;
}

void crossover_odd(long *idum,sim_par par,ring_data ring,ring_data &ring_1)
{
  vector<interval> tmp_pos_1;
  vector<int> cn;
  long idx,bp,n;
  
  ++ring.cycle_tot;
  ++ring.cycle_odd;

  ring_1=ring;
  get_bkp(idum,ring_1.c_ring,idx,bp);
  add_bkp_recenter(idx,bp,ring_1.c_ring);

  //duplicate ring
  n=ring_1.c_ring.size();
  for(long i=0;i<n;++i)
    {
      ring_1.c_ring.push_back(ring_1.c_ring[i]);
    }
  
  ring_size(ring_1);
  calc_cn(par,ring_1,cn);
  ring_1.fitness=calc_fitness(par,ring_1);
}

void get_insert(long *idum,long r_size,vector<interval> &c_ring,vector<interval> &insert)
{
  long max_insert_size=(long)((double)r_size/10.0);
  long insert_size=(long)(ran1(idum)*(double)max_insert_size);
  long sum,h,pos;
  interval tmp_seg;
  vector<interval> tmp_c_ring;
  insert.clear();

  sum=0;
  for(h=0;h<c_ring.size();h++)
    {
      if(sum+c_ring[h].end-c_ring[h].start+1 > insert_size)
	break;
      sum+=c_ring[h].end-c_ring[h].start+1;
    }
  pos=c_ring[h].start+insert_size-sum;

  tmp_seg.end=c_ring[h].end;
  c_ring[h].end=pos;
  tmp_seg.start=pos+1;
  c_ring.insert(c_ring.begin()+h+1,tmp_seg);

  for(long i=0;i<=h;++i)
    insert.push_back(c_ring[i]);

  for(long i=h+1;i<c_ring.size();++i)
    tmp_c_ring.push_back(c_ring[i]);
  
  c_ring=tmp_c_ring;

}

void crossover_even(long *idum,sim_par par,ring_data ring,ring_data &ring_1,ring_data &ring_2)
{
  vector<interval> insert_1,insert_2;
  vector<int> cn;
  long idx,bp;

  ++ring.cycle_tot;
  ++ring.cycle_even;

  ring_1=ring;
  get_bkp(idum,ring_1.c_ring,idx,bp);
  add_bkp_recenter(idx,bp,ring_1.c_ring);
  get_insert(idum,ring_1.size,ring_1.c_ring,insert_1);

  ring_2=ring;
  get_bkp(idum,ring_2.c_ring,idx,bp);
  add_bkp_recenter(idx,bp,ring_2.c_ring);
  get_insert(idum,ring_2.size,ring_2.c_ring,insert_2);

  //assemble rings R1+I2 and R2+I1
  for(long i=0;i<insert_2.size();++i)
    ring_1.c_ring.push_back(insert_2[i]);

  for(long i=0;i<insert_1.size();++i)
    ring_2.c_ring.push_back(insert_1[i]);
  
  ring_size(ring_1);
  calc_cn(par,ring_1,cn);
  ring_1.fitness=calc_fitness(par,ring_1);

  ring_size(ring_2);
  calc_cn(par,ring_2,cn);
  ring_2.fitness=calc_fitness(par,ring_2);
}


void output_stats(double t,sim_par &par, vector<ring_data> ring)
{
  ofstream out;
  vector<int> cn;
  vector<double> avg_cn(par.glength+1,1);
  vector<double> std_cn(par.glength+1,0);
  double avg_cycles_tot=0;
  double avg_cycles_even=0;
  double avg_cycles_odd=0;
  double avg_fitness=0;
  double avg_fitness_nozero=0;
  double avg_ring_size=0;
  long n_nozero=0;
  
  out.open((par.outname+"_statistics.txt").c_str(),ios::out | ios::app);

  long n=ring.size();

  for(long i=0;i<n;++i)
    {
      avg_cycles_tot+=(double)ring[i].cycle_tot;
      avg_cycles_even+=(double)ring[i].cycle_even;
      avg_cycles_odd+=(double)ring[i].cycle_odd;
      avg_ring_size+=(double)ring[i].size;
      avg_fitness+=ring[i].fitness;
      if(ring[i].fitness>0)
	{
	  ++n_nozero;
	  avg_fitness_nozero+=ring[i].fitness;
	}
    }
  avg_cycles_tot*=1.0/(double)n;
  avg_cycles_even*=1.0/(double)n;
  avg_cycles_odd*=1.0/(double)n;
  avg_ring_size*=1.0/(double)n;
  avg_fitness*=1.0/(double)n;
  avg_fitness_nozero*=1.0/(double)n_nozero;


  cout << "time=" << t << "\tN=" << n << "\tavg_fitness (all="  << avg_fitness << "; not_zero=" << avg_fitness_nozero << " (N=" << n_nozero << ")) ";
  cout << "\tavg_ring_size[Mb]=" << avg_ring_size/1e6;
  cout << "\tavg_cycles (tot=" << avg_cycles_tot << "; odd=" << avg_cycles_odd << "; even=" << avg_cycles_even << ")\n";
  out << t << "\t" << n << "\t" << avg_fitness << "\t" << avg_fitness_nozero << "\t" << n_nozero << "\t" <<  avg_ring_size/1e6 << "\t";
  out << avg_cycles_tot << "\t" << avg_cycles_odd << "\t" << avg_cycles_even << endl;
  out.close();

  // write seg files
  stringstream line;
  line.str(""); line.clear();
  line << par.seg_folder << "/" << par.outname << "-time_" << t << ".seg";
  out.open((line.str()).c_str());
  if(!out.is_open())
    {
      cerr << "Error: cannot open seg file.\n";
      exit(1);
    }

  out << "Sample\tChromosome\tStart\tEnd\tnMarker\tCN\n";
  
  for(long i=0;i<n;++i)
    {
      for(long h=0;h<ring[i].seg.size();++h)
	{
	  out << "R" << i+1 << "\tchr12\t" << ring[i].seg[h].start << "\t";
	  out << ring[i].seg[h].end << "\t" << ring[i].seg[h].end-ring[i].seg[h].start+1 << "\t" << ring[i].cn_seg[h]+1 << endl;
	}
    }

  out.close();
}
