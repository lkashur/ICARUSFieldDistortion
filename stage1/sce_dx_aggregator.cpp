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
ROOT::RVec<int> get_3D_bins(ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, TH3F *DummyHist3D);
void create_bin_map(ROOT::RVec<int> bin_ids, ROOT::RVec<float> offsets, unordered_map< int, vector<float> > bin_vals);
float Median(vector<float> v, int n);

//////////////////////////////////////////////////
/// Main Analysis
//////////////////////////////////////////////////
int main(int argc, char** argv)
{

  // Load data
  Char_t *inputfilename = (Char_t*)"";
  inputfilename = argv[1];

  TChain* inputfiles = new TChain("offsettree");
  AddFiles(inputfiles, inputfilename);

  // Create rdataframe
  ROOT::RDataFrame rdf(*inputfiles);

  // Output
  char *tpc;
  char *sim_or_reco;
  float xrange0;
  float xrange1;
  tpc = (char*) argv[2];
  sim_or_reco = (char*) argv[3];
  int sel_tpc;
  if (!strcmp(tpc, "EE"))
    {
      xrange0=-358.49;
      xrange1=-210.29;
      sel_tpc = 0;
    }
  else if (!strcmp(tpc, "EW"))
    {
      xrange0=-210.14;
      xrange1=-61.94;
      sel_tpc = 1;
    }
  else if (!strcmp(tpc, "WE"))
    {
      xrange0=61.94;
      xrange1=210.14;
      sel_tpc = 2;
    }
  else if(!strcmp(tpc, "WW"))
    {
      xrange0=210.29;
      xrange1=358.49;
      sel_tpc = 3;
    }
  int Nxbins = 15;
  int Nybins = 10; //10
  int Nzbins = 8; //15
  TH3F *DummyHist3D = new TH3F("DummyHist3D","",Nxbins, xrange0, xrange1, Nybins, -181.86, 134.96, Nzbins, -894.9515, 894.9515);

  ofstream bindetails;
  ofstream alloffsets;
  string bindetails_text = string("bindetails_") + tpc + string(".txt");
  bindetails.open(bindetails_text);
  string alloffsets_text = string("alloffsets_") + tpc + string(".txt");
  alloffsets.open(alloffsets_text);
  string aggoffsets_text = string("aggoffsets_") + tpc + string(".root");
  TFile outf(aggoffsets_text.c_str(), "recreate");
  unordered_map< int, vector<float> > bin_vals;

  // Cuts on anode-cathode side of track
  auto rdf_cuts = rdf.Filter([&](int tpc){ return tpc == sel_tpc; }, {"tpc"} )
                     .Filter("min_time_y < 134.96 - 15")
                     .Filter("min_time_y > -181.86 + 15")
                     .Filter("min_time_z > -894.9515 + 100")
                     .Filter("min_time_z < 894.9515 - 100")
                     .Filter("max_time_y < 134.96 - 85")
                     .Filter("max_time_y > -181.86 + 85")
                     .Filter("max_time_z > -894.9515 + 100")
                     .Filter("max_time_z < 894.9515 - 100");
  
  // Bin in 3D
  bindetails << "gbn,yzbin,bin_low_x,bin_high_x,bin_low_y,bin_high_y,bin_low_z,bin_high_z" << endl;
  auto rdf_bins = rdf_cuts.Define("bin_ids", [&] (ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs)->ROOT::RVec<int>{return get_3D_bins(xs, ys, zs, DummyHist3D);}, {"sel_xs", "sel_ys", "sel_zs"});

  for (int i = 0; i <= DummyHist3D->GetNbinsX()+1; ++i)
    {
      float bin_low_x = DummyHist3D->GetXaxis()->GetBinLowEdge(i);
      float bin_high_x = DummyHist3D->GetXaxis()->GetBinUpEdge(i);

      for (int j = 0; j <= DummyHist3D->GetNbinsY()+1; ++j)
        {
          float bin_low_y = DummyHist3D->GetYaxis()->GetBinLowEdge(j);
          float bin_high_y = DummyHist3D->GetYaxis()->GetBinUpEdge(j);

          for (int k = 0; k <= DummyHist3D->GetNbinsZ()+1; ++k)
            {
              float bin_low_z = DummyHist3D->GetZaxis()->GetBinLowEdge(k);
              float bin_high_z = DummyHist3D->GetZaxis()->GetBinUpEdge(k);

              int gbn = DummyHist3D->GetBin(i, j, k);
              string yzstring = to_string(j) + "_" + to_string(k);
              bindetails << gbn << "," << yzstring << "," << bin_low_x << "," << bin_high_x << "," << bin_low_y << "," << bin_high_y << "," << bin_low_z << "," << bin_high_z << endl;

              // Create emtpy vector element for every bin
              bin_vals.insert(pair<int,vector<float> >(gbn, vector<float>()));

            } // end z loop
        } // end y loop
    } // end x loop 
  

  // Fill offset vectors for aggregation
  // Reco
  if(!strcmp(sim_or_reco, "reco"))
    {
      rdf_bins.Foreach([&] (ROOT::RVec<int> bin_ids, ROOT::RVec<float> offsets) {
	  for(int i=0; i < offsets.size(); i++)
	    {
	      //cout << bin_ids[i] << " " << offsets[i] << endl;
	      bin_vals[bin_ids[i]].push_back(offsets[i]);
	    }
	}, {"bin_ids", "offsets"}); 
    }
  // Sim
  else if(!strcmp(sim_or_reco, "sim"))
    {
      unique_ptr<TFile> simFile(new TFile("/cvmfs/icarus.opensciencegrid.org/products/icarus/icarus_data/v09_84_00/icarus_data/SCEoffsets/SCEoffsets_ICARUS_E500_voxelTH3_2DSigSimHack.root", "READ"));
      TH3F* hTrueBkwdX = (TH3F*) simFile->Get("TrueBkwd_Displacement_X");
      rdf_bins.Foreach([&] (ROOT::RVec<int> bin_ids, ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs) {
	  for(int i=0; i < xs.size(); i++)
	    {
	      bin_vals[bin_ids[i]].push_back(-hTrueBkwdX->Interpolate(xs[i],ys[i],zs[i])); 
	    }
	}, {"bin_ids", "sel_xs", "sel_ys", "sel_zs"});
      simFile->Close();
    }
  // Aggregation
  alloffsets << "gbn,offsets" <<endl;
  for (auto const& b : bin_vals)
    {
      float bin_agg = -99999.0; 
      if(!b.second.empty())
        { 
          // Calculate aggregated offset
          bin_agg = Median(b.second, b.second.size());
	  
          // Store every bin's in output text file
          alloffsets << b.first;
          for(auto val : b.second)
	    {
	      alloffsets << "," << val;
	    }
          alloffsets << endl;
        }

      // Store aggregagted offsets in text file/root TH3
      //aggoffsets << b.first << "," << bin_agg << endl;
      DummyHist3D->SetBinContent(b.first, bin_agg);
    }
  
  bindetails.close();
  alloffsets.close();
  outf.cd();
  DummyHist3D->Write();
  outf.Close();

  //rdf_bins.Display({"tpc", "sel_xs", "sel_ys", "sel_zs", "offsets", "bin_ids"})->Print();
  cout << "No. tracks: " << rdf_bins.Count().GetValue() << endl;

  
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

ROOT::RVec<int> get_3D_bins(ROOT::RVec<float> xs, ROOT::RVec<float> ys, ROOT::RVec<float> zs, TH3F *DummyHist3D)
{
  ROOT::RVec<int> bin_ids;
  for(int i=0; i < xs.size(); i++)
    {
      bin_ids.push_back(DummyHist3D->FindBin(xs[i], ys[i], zs[i]));
    }
  
  return bin_ids;
}

void create_bin_map(ROOT::RVec<int> bin_ids, ROOT::RVec<float> offsets, unordered_map< int, vector<float> > bin_vals)
{
  for(int i=0; i < offsets.size(); i++)
    {
      //cout << bin_ids[i] << " " << offsets[i] << endl;
      bin_vals[bin_ids[i]].push_back(offsets[i]);
    }
  return;
} 

float Median(vector<float> v, int n)
{
  // Sort the vector
  sort(v.begin(), v.end());

  // Check if the number of elements is odd
  if (n % 2 != 0)
    {
      return (float)v[n / 2];
    }
  else
    {
      return (float)(v[(n - 1) / 2] + v[n / 2]) / 2.0;
    }
}
