// C++ Includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <cmath>

// ROOT Includes
#include "TROOT.h"
#include "ROOT/TThreadedObject.hxx"
#include <ROOT/RDataFrame.hxx>
#include "TTree.h"
#include "TFile.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TH1.h"
#include "TH2.h"
#include "TTreeReader.h"
#include "TTreeReaderArray.h"
#include "TVirtualFFT.h"
#include "TFile.h"
#include "TChain.h"
#include "TSystemDirectory.h"

using namespace std;

// Function Declaration
void AddFiles(TChain *ch, const char *fileList);
int get_tpc_ew(float a, float b);
int get_tpc(int cryo, int tpc_ew);
ROOT::RVec<float> get_track_hits(ROOT::RVec<bool> htraj, ROOT::RVec<float> hits, int tpc_ew, ROOT::RVec<unsigned short> tpcs);
vector<float> get_min_time_xyz(ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z);
vector<float> get_max_time_xyz(ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z);
ROOT::RVec<ROOT::RVec<float>> get_offsets(ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, ROOT::RVec<float> ts, float min_time_x, float min_time_z, float slope);

//////////////////////////////////////////////////
/// Main Analysis
//////////////////////////////////////////////////
int main(int argc, char** argv)
{

  // Output

  // Load data
  Char_t *inputfilename = (Char_t*)"";
  char *cryo;
  inputfilename = argv[1];
  cryo = (char*) argv[2];

  char* treename;
  if (!strcmp(cryo, "E"))
    {
      treename = (char*) "caloskimE/TrackCaloSkim";
    }
  else if (!strcmp(cryo, "W"))
    {
      treename = (char*) "caloskimW/TrackCaloSkim";
    }
      
  TChain* inputfiles = new TChain(treename);
  AddFiles(inputfiles, inputfilename);

  // Create rdataframe
  //ROOT::EnableImplicitMT();
  ROOT::RDataFrame rdf(*inputfiles);  

  // TPC info
  auto rdf_tpc = rdf.Define("tpc_ew", [&] (float hit_max_time_p2_tpcE, float hit_min_time_p2_tpcE, float hit_max_time_p2_tpcW, float hit_min_time_p2_tpcW)->int{return get_tpc_ew(hit_max_time_p2_tpcE - hit_min_time_p2_tpcE, hit_max_time_p2_tpcW - hit_min_time_p2_tpcW);}, {"hit_max_time_p2_tpcE", "hit_min_time_p2_tpcE", "hit_max_time_p2_tpcW", "hit_min_time_p2_tpcW"}).Define("tpc", [&] (int cryostat, int tpc_ew)->int{return get_tpc(cryostat, tpc_ew);}, {"cryostat", "tpc_ew"});
  
  // Anode-cathode crossers
  auto rdf_ac = rdf_tpc.Filter("selected == 1")
                       .Filter("tpc != -1");


  // Hit filtering
  auto rdf_hits = rdf_ac.Define("hit_xs", [&] (ROOT::RVec<bool> htraj, ROOT::RVec<float> hits, int tpc_ew, ROOT::RVec<unsigned short> tpcs)->ROOT::RVec<float>{return get_track_hits(htraj, hits, tpc_ew, tpcs);}, {"hits2.ontraj", "hits2.h.sp.x", "tpc_ew", "hits2.h.tpc"})
                        .Define("hit_ys", [&] (ROOT::RVec<bool> htraj, ROOT::RVec<float> hits, int tpc_ew, ROOT::RVec<unsigned short> tpcs)->ROOT::RVec<float>{return get_track_hits(htraj, hits, tpc_ew, tpcs);}, {"hits2.ontraj", "hits2.h.sp.y", "tpc_ew", "hits2.h.tpc"})
                        .Define("hit_zs", [&] (ROOT::RVec<bool> htraj, ROOT::RVec<float> hits, int tpc_ew, ROOT::RVec<unsigned short> tpcs)->ROOT::RVec<float>{return get_track_hits(htraj, hits, tpc_ew, tpcs);}, {"hits2.ontraj", "hits2.h.sp.z", "tpc_ew", "hits2.h.tpc"})
                        .Define("hit_ts", [&] (ROOT::RVec<bool> htraj, ROOT::RVec<float> hits, int tpc_ew, ROOT::RVec<unsigned short> tpcs)->ROOT::RVec<float>{return get_track_hits(htraj, hits, tpc_ew, tpcs);}, {"hits2.ontraj", "hits2.h.time", "tpc_ew", "hits2.h.tpc"});

  // Get track endpoints
  auto rdf_endpoints = rdf_hits.Define("min_time_t", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_min_time_xyz(t,x,y,z)[0];}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("min_time_x", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_min_time_xyz(t,x,y,z)[1]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("min_time_y", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_min_time_xyz(t,x,y,z)[2]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("min_time_z", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_min_time_xyz(t,x,y,z)[3]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("max_time_t", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_max_time_xyz(t,x,y,z)[0]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("max_time_x", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_max_time_xyz(t,x,y,z)[1]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("max_time_y", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_max_time_xyz(t,x,y,z)[2]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"})
    .Define("max_time_z", [&] (ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)->float{return get_max_time_xyz(t,x,y,z)[3]\
	  ;}, {"hit_ts", "hit_xs", "hit_ys", "hit_zs"});
  
  // Calculate offsets
  auto rdf_dx = rdf_endpoints.Define("slope", [&] (float max_time_x, float min_time_x, float max_time_z, float min_time_z)->float{return (max_time_x - min_time_x)/(max_time_z - min_time_z);}, {"max_time_x", "min_time_x", "max_time_z", "min_time_z"})
    .Define("sel_xs", [&] (ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, ROOT::RVec<float> ts, float min_time_x, float min_time_z, float slope)->ROOT::RVec<float>{return get_offsets(xs, ys, zs, ts, min_time_x, min_time_z, slope)[0];}, {"hit_xs", "hit_ys", "hit_zs", "hit_ts", "min_time_x", "min_time_z", "slope"})
    .Define("sel_ys", [&] (ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, ROOT::RVec<float> ts, float min_time_x, float min_time_z, float slope)->ROOT::RVec<float>{return get_offsets(xs, ys, zs, ts, min_time_x, min_time_z, slope)[1];}, {"hit_xs", "hit_ys", "hit_zs", "hit_ts", "min_time_x", "min_time_z", "slope"})
    .Define("sel_zs", [&] (ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, ROOT::RVec<float> ts, float min_time_x, float min_time_z, float slope)->ROOT::RVec<float>{return get_offsets(xs, ys, zs, ts, min_time_x, min_time_z, slope)[2];}, {"hit_xs", "hit_ys", "hit_zs", "hit_ts", "min_time_x", "min_time_z", "slope"})
    .Define("offsets", [&] (ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, ROOT::RVec<float> ts, float min_time_x, float min_time_z, float slope)->ROOT::RVec<float>{return get_offsets(xs, ys, zs, ts, min_time_x, min_time_z, slope)[3];}, {"hit_xs", "hit_ys", "hit_zs", "hit_ts", "min_time_x", "min_time_z", "slope"});

  // Save new output tree
  string output_file_text = string("offsets_") + cryo + string(".root");
  rdf_dx.Snapshot("offsettree", output_file_text.c_str(), {"tpc", "min_time_y", "min_time_z", "max_time_y", "max_time_z", "sel_xs", "sel_ys", "sel_zs", "offsets"});

  return 0;

} // end main

/////////////////////////////////////////////////
/// Functions
/////////////////////////////////////////////////

// Add intput files to TChain
void AddFiles(TChain *ch, const char *fileList)
{
  ifstream InFile(fileList);
  vector<string> fileVec;
  string line;
  int n_files = 0;
  while(getline(InFile, line))
    {
      fileVec.push_back(line);
      n_files++;
    }

  for(unsigned int iFile=0; iFile<fileVec.size(); ++iFile)
    {
      //cout << fileVec[iFile].c_str() << endl;
      ch->AddFile(fileVec[iFile].c_str());
    }
  return;
}

int get_tpc_ew(float a, float b)
{
  if(a > b)
    {
      return 0;
    }
  else
    {
      return 1;
    }
}

int get_tpc(int cryo, int tpc_ew)
{
  if(cryo == 0 && tpc_ew == 0)
    {
      return 0;
    }
  else if(cryo == 0 && tpc_ew == 1)
    {
      return 1;
    }
  else if(cryo == 1 && tpc_ew == 0)
    {
      return 2;
    }
  else if (cryo == 1 && tpc_ew == 1)
    {
      return 3;
    }
  else
    {
      return -1;
    }
}

ROOT::RVec<float> get_track_hits(ROOT::RVec<bool> htraj, ROOT::RVec<float> hits, int tpc_ew, ROOT::RVec<unsigned short> tpcs)
{
  ROOT::RVec<float> sel_hits;
  for(int i=0; i < hits.size(); i++)
    {
      if((htraj[i] == true) && (((tpc_ew == 0) && (tpcs[i] <2)) || ((tpc_ew == 1) && (tpcs[i] >= 2))))
	{
          sel_hits.push_back(hits[i]);
        }
    }
  return sel_hits;
}

vector<float> get_min_time_xyz(ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)
{
  vector<float> txyz;
  float min_time = 99999;
  float min_time_x;
  float min_time_y;
  float min_time_z;
  for(int i = 0; i < t.size(); i++)
    {
      if(t[i] < min_time)
        {
          min_time = t[i];
          min_time_x = x[i];
          min_time_y = y[i];
          min_time_z = z[i];
        }
    }

  txyz.push_back(min_time);
  txyz.push_back(min_time_x);
  txyz.push_back(min_time_y);
  txyz.push_back(min_time_z);
  return txyz;
}

vector<float> get_max_time_xyz(ROOT::RVec<float> t, ROOT::RVec<float> x, ROOT::RVec<float> y, ROOT::RVec<float> z)
{
  vector<float> txyz;
  float max_time = -99999;
  float max_time_x;
  float max_time_y;
  float max_time_z;
  for(int i = 0; i < t.size(); i++)
    {
      if(t[i] > max_time)
        {
          max_time = t[i];
          max_time_x = x[i];
          max_time_y = y[i];
          max_time_z = z[i];
        }
    }

  txyz.push_back(max_time);
  txyz.push_back(max_time_x);
  txyz.push_back(max_time_y);
  txyz.push_back(max_time_z);
  return txyz;
}

ROOT::RVec<ROOT::RVec<float>> get_offsets(ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, ROOT::RVec<float> ts, float min_time_x, float min_time_z, float slope)
{
  ROOT::RVec<ROOT::RVec<float>> xyz_dx;

  ROOT::RVec<float> sel_xs;
  ROOT::RVec<float> sel_ys;
  ROOT::RVec<float> sel_zs;
  ROOT::RVec<float> offsets;
  for(int i = 0; i < xs.size(); i++)
    {
      float ideal_x = slope*(zs[i] - min_time_z) + min_time_x;
      float dx = xs[i] - ideal_x;
      
      if(dx == 0.0) continue;

      sel_xs.push_back(xs[i]);
      sel_ys.push_back(ys[i]);
      sel_zs.push_back(zs[i]);
      offsets.push_back(dx);
    }
  xyz_dx.push_back(sel_xs);
  xyz_dx.push_back(sel_ys);
  xyz_dx.push_back(sel_zs);
  xyz_dx.push_back(offsets);
  return xyz_dx;
}



ROOT::RVec<ROOT::RVec<float>> get_offsets(ROOT::RVec<float> ts, ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, float slope, float min_time_t, float max_time_t, float min_time_x, float min_time_z)
{
  ROOT::RVec<float> sel_xs;
  ROOT::RVec<float> sel_ys;
  ROOT::RVec<float> sel_zs;
  ROOT::RVec<float> offsets;
  ROOT::RVec<ROOT::RVec<float>> xyz_dx;
  for(int i = 0; i < xs.size(); i++)
    {
      float ideal_x = slope*(zs[i] - min_time_z) + min_time_x;
      float dx = xs[i] - ideal_x;
      if(dx == 0.0)
        {
          continue;
        }
      sel_xs.push_back(xs[i]);
      sel_ys.push_back(ys[i]);
      sel_zs.push_back(zs[i]);
      offsets.push_back(dx);
    }

  xyz_dx.push_back(sel_xs);
  xyz_dx.push_back(sel_ys);
  xyz_dx.push_back(sel_zs);
  xyz_dx.push_back(offsets);
  return xyz_dx;
}
