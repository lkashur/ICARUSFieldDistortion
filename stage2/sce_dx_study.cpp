// C++ Includes
#include <iostream>
#include <iomanip>
#include <fstream>
#include <unordered_map>
#include <cmath>

// ROOT Includes
#include "TROOT.h"
#include <TStyle.h>
#include "ROOT/TThreadedObject.hxx"
#include <ROOT/RDataFrame.hxx>
#include "TTree.h"
#include "TFile.h"
#include "TSystemFile.h"
#include "TString.h"
#include "TF1.h"
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
vector<double> linspace(double start, double end, int n);

//////////////////////////////////////////////////
// Main Analysis
//////////////////////////////////////////////////
int main(int argc, char** argv)
{
  // Drift Velocity vs. E-Field fit
  TF1 driftVelFit("driftVelFit","(1.0 - 0.0184 * (87.5 - 89.0)) * pol5",-0.05,1.1);
  driftVelFit.SetParameter(0,0.0);
  driftVelFit.SetParameter(1,5.53416);
  driftVelFit.SetParameter(2,-6.53093);
  driftVelFit.SetParameter(3,3.20752);
  driftVelFit.SetParameter(4,0.389696);
  driftVelFit.SetParameter(5,-0.556184);

  // Load data
  Char_t *inputfilename = (Char_t*)"";
  inputfilename = argv[1];
  unique_ptr<TFile> infile( TFile::Open(inputfilename) );
  unique_ptr<TH3> hAgg(infile->Get<TH3>("DummyHist3D"));
  
  char *tpc;
  float xrange0;
  float xrange1;
  tpc = (char*) argv[2];
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

  // Output file
  string scedxstudy_text = string("scedxstudy_") + tpc + string(".root");
  TFile outf(scedxstudy_text.c_str(), "recreate");
  outf.cd();

  // Make output TH3, starting from a copy of input TH3
  TH3D *hRho3D = new TH3D("hRho3D","",hAgg->GetNbinsX(), xrange0, xrange1, hAgg->GetNbinsY(), -181.86, 134.96, hAgg->GetNbinsZ(), -894.9515, 894.9515);

  // Loop over {y,z} bins
  for (int j = 0; j <= hAgg->GetNbinsY()+1; j++)
    {
      float bin_low_y = hAgg->GetYaxis()->GetBinLowEdge(j);
      float bin_high_y = hAgg->GetYaxis()->GetBinUpEdge(j);

      for (int k = 0; k <= hAgg->GetNbinsZ()+1; k++)
        {

          float bin_low_z = hAgg->GetZaxis()->GetBinLowEdge(k);
          float bin_high_z = hAgg->GetZaxis()->GetBinUpEdge(k);

	  // Create unique (output) directory for every {y,z} bin
	  string dirstring = to_string(j) + "_" + to_string(k);
          TDirectory* yzbindir = outf.mkdir(dirstring.c_str());
          yzbindir->cd();

	  ////////////////////////////////////////////////////////////
	  /// Spatial Offsets
	  ////////////////////////////////////////////////////////////

	  // "Project" a single {y,,z} bin along x to get dx vs. x histo
	  TH1D *hdx = hAgg->ProjectionX("hdx", j, j, k, k, "o");

	  // Fit dx vs. x with polynomial
	  TF1 *dxfit = new TF1("dxfit","pol3",xrange0,xrange1);
          hdx->Fit("dxfit", "MRQF");
	  hdx->SetTitle("Spatial Offsets vs. Drift Length");
	  hdx->GetXaxis()->SetTitle("Drift Length [cm]");
	  hdx->GetYaxis()->SetTitle("Spatial Offset dX [cm]");
	  hdx->SetMarkerStyle(21);
	  hdx->Write();
	  
	  /////////////////////////////////////////////////////////
	  /// Local Drift Velocity
	  /////////////////////////////////////////////////////////
	  
	  // Sample dx vs. x fit on a linspace
	  vector<double> tempx = linspace(xrange0, xrange1, 100);
          vector<double> tempy;
          for(const auto& x : tempx)
            {
              tempy.push_back(dxfit->Eval(x));
            }
	  
	  // Calculate difference between adjacent dxs
	  vector<double> ddxs;
          adjacent_difference(tempy.begin(), tempy.end(), back_inserter(ddxs));
          ddxs.erase(ddxs.begin());
	  
	  // Velocity calculation
	  vector<double> loc_v;
	  double spacing = tempx[1] - tempx[0];
          for(const auto& ddx : ddxs)
            {
              float vel = 0.157565 * (1 - ddx / spacing);
              loc_v.push_back(vel);
            }
	  
	  vector<double> xadj;
          for(const auto& x : tempx)
            {
              xadj.push_back(x + spacing/2);
            }
          xadj.pop_back();

	  TGraph *gVel = new TGraph(xadj.size(), xadj.data(), loc_v.data());
          gVel->SetName("gVel");
          gVel->SetTitle("Drift Velocity vs. Drift Length");
	  gVel->GetXaxis()->SetTitle("Drift Length [cm]");
	  gVel->GetYaxis()->SetTitle("Drift Velocity [cm/#mus]");
	  gVel->Write();

	  /////////////////////////////////////////////////////////////
	  /// Local Electric Field
	  /////////////////////////////////////////////////////////////
	  
	  // Get E-field from known relationship with drift velocity
	  vector<double> loc_E;
          for(const auto& v : loc_v)
            {
              loc_E.push_back(driftVelFit.GetX(10*v));
            }

	  // E-field vs. drift length
	  TGraph *gEfield = new TGraph(xadj.size(), xadj.data(), loc_E.data());
          gEfield->SetName("gEfield");
	  gEfield->SetTitle("Electric Field vs. Drift Length");
	  gEfield->GetXaxis()->SetTitle("Drift Length [cm]");
	  gEfield->GetYaxis()->SetTitle("Electric Field [kV/cm]");
          TF1 *efit = new TF1("efit","pol3",xrange0,xrange1);
          gEfield->Fit("efit", "MRQF");
	  gEfield->Write();

	  ///////////////////////////////////////////////////////////////////
	  /// Local Space Charge Density
	  ///////////////////////////////////////////////////////////////////

	  // Sample dx vs. x fit on a linspace
	  vector<double> rhovals;
          for(const auto& x : tempx)
            {
              rhovals.push_back((pow(10.0,9))*0.00013316416*(efit->Derivative(x)));
            }	  

	  // Space charge density vs. drift length
          TGraph *gRho = new TGraph(tempx.size(), tempx.data(), rhovals.data());
          gRho->SetName("gRho");
	  gRho->SetTitle("Space Charge Density vs. Drift Length");
	  gRho->GetXaxis()->SetTitle("Drift Length [cm]");
	  gRho->GetYaxis()->SetTitle("Space Charge Density [nC/m^{3}]");
	  gRho->Write();

	  // Average space charge density over bins of x
	  for (int i = 0; i <= hAgg->GetNbinsX()+1; i++)
            {
              float bin_low_x = hAgg->GetXaxis()->GetBinLowEdge(i);
              float bin_high_x = hAgg->GetXaxis()->GetBinUpEdge(i);

              int gbn = hAgg->GetBin(i, j, k);
	      
	      // Sum rho entries for specific x-bin
              vector<double> xbin_rhovals;
              for(int s=0; s<rhovals.size(); s++)
                {
                  if(tempx[s] > bin_low_x && tempx[s] < bin_high_x)
                    {
                      xbin_rhovals.push_back(rhovals[s]);
                    }
                }

	      // Average caculation
              double sum = 0;
              for(const auto& val : xbin_rhovals)
		{
                  sum += val;
                }
	      double xbin_rhomean = sum / xbin_rhovals.size();
	      hRho3D->SetBinContent(gbn, xbin_rhomean);
	      
	    } // end x loop
	} // end z loop
    } // end y loop
  

  /////////////////////////////////////////////////////////////
  /// Clean up ... Using hRho3D, fill ONLY good bins of new TH3
  //////////////////////////////////////////////////////////////
  /*
  for (int i = 0; i <= hRho3D->GetNbinsX()+1; i++)
    {
      float bin_cen_x = hRho3D->GetXaxis()->GetBinCenter(i);
      for (int j = 0; j <= hRho3D->GetNbinsY()+1; j++)
	{
	  float bin_cen_y = hRho3D->GetYaxis()->GetBinCenter(j);
	  for (int k = 0; k <= hRho3D->GetNbinsY()+1; k++)
	    {
	      float bin_cen_z = hRho3D->GetZaxis()->GetBinLowEdge(k);
	      int gbn = hRho3D->GetBin(i, j, k);

	      // exclude nasty bins
	      if(bin_cen_y > -130 && bin_cen_y < 85 && bin_cen_z > -795 && bin_cen_z < 795)
		{
		  hRho3D_clean->SetBinContent(gbn, hRho3D->GetBinContent(gbn));
		}
	      
	    }
	}
    }
  */

  outf.cd();
  hRho3D->GetXaxis()->SetTitle("x [cm]");
  hRho3D->GetYaxis()->SetTitle("y [cm]");
  hRho3D->GetZaxis()->SetTitle("z [cm]");
  hRho3D->Write();
  TH2D *hRhoXY = (TH2D*) hRho3D->Project3D("yx");
  hRhoXY->SetTitle("Average Space Charge Density [nC/m^{3}]");
  hRhoXY->GetXaxis()->SetTitle("x [cm]");
  hRhoXY->GetYaxis()->SetTitle("y [cm]");
  hRhoXY->Write();

  /*
  hRho3D_clean->GetXaxis()->SetTitle("x [cm]");
  hRho3D_clean->GetYaxis()->SetTitle("y [cm]");
  hRho3D_clean->GetZaxis()->SetTitle("z [cm]");
  hRho3D_clean->Write();
  TH2D *hRhoXY_clean = (TH2D*) hRho3D_clean->Project3D("yx");
  hRhoXY_clean->SetTitle("Average Space Charge Density [nC/m^{3}]");
  hRhoXY_clean->GetXaxis()->SetTitle("x [cm]");
  hRhoXY_clean->GetYaxis()->SetTitle("y [cm]");
  hRhoXY_clean->Write();
  outf.Close();
  */
  return 0;

} // end main

/////////////////////////////////////////////////
/// Functions
/////////////////////////////////////////////////

// Create linspace, assume n > 1
vector<double> linspace(double start, double end, int n)
{
  // output
  vector<double> xs;
  
  double delta = (end - start) / (n - 1);
  for(int i=0; i < n-1; ++i)
    {
      xs.push_back(start + delta * i);
    }

  xs.push_back(end);
  return xs;
}


