# imports
import ROOT
import numpy as np
import pandas as pd
import argparse

def main(args):
    
    # Load/prepare histograms
    outfile = ROOT.TFile.Open('teststudy.root', 'RECREATE')
    histfile = ROOT.TFile.Open(args.in_file_hist, 'READ')
    dummyhist = histfile.Get('DummyHist3D')
    nbins_x = dummyhist.GetNbinsX()
    nbins_y = dummyhist.GetNbinsY()
    nbins_z = dummyhist.GetNbinsZ()
    hreco = dummyhist.Clone('hreco')
    hsim = dummyhist.Clone('hsim')

    # Load reco/sim aggregated offsets into pandas dataframes
    dfs = []
    for arg in vars(args):
        if (arg == 'in_file_reco' or arg == 'in_file_sim') and (getattr(args, arg)):
            df = pd.read_csv(getattr(args,arg))
            df['gbn'] = df['gbn'].astype(int)
            dfs.append(df)

    # Fill TH3s with aggregated offsets
    for i in range(dummyhist.GetNbinsX() + 2):
        for j in range(dummyhist.GetNbinsY() + 2):
            for k in range(dummyhist.GetNbinsZ() + 2):
                gbn = dummyhist.GetBin(i, j, k)
                # Reco
                if gbn in dfs[0]['gbn'].values:
                    hreco.SetBinContent(gbn, dfs[0][dfs[0]['gbn'] == gbn]['med'].values[0])
                # Sim
                if (len(dfs) > 1) and (gbn in dfs[1]['gbn'].values):
                    hsim.SetBinContent(gbn, dfs[1][dfs[1]['gbn'] == gbn]['med'].values[0])

    # Analysis
    outfile.cd()
    print(hreco.GetNbinsZ())
    for j in range(dummyhist.GetNbinsY() + 2):
        for k in range(dummyhist.GetNbinsZ() + 2):
            
            print(j,k)

            # We are now working in a distinct y-z bin
            # Create a directory for this bin
            dirstring = str(j) + "_" + str(k)
            yzbindir = outfile.mkdir(dirstring)
            yzbindir.cd()

            ##################################################
            ### dx vs x
            ##################################################
            # reco
            hdx = hreco.ProjectionX('hdx', j, j, k, k, 'o')
            gx, gx_low, gx_high = [], [], []
            gy, gy_low, gy_high = [], [], []
            for i in range(hreco.GetNbinsX() + 2):
                gbn = hreco.FindBin(hdx.GetXaxis().GetBinCenter(i), hreco.GetYaxis().GetBinCenter(j), hreco.GetZaxis().GetBinCenter(k))
                gx.append(hdx.GetXaxis().GetBinCenter(i))
                gx_low.append(hdx.GetXaxis().GetBinWidth(i)/2)
                gx_high.append(hdx.GetXaxis().GetBinWidth(i)/2)
                gy.append(hdx.GetBinContent(i))
                if gbn in dfs[0].gbn.values:
                    gy_low.append(hdx.GetBinContent(i) - dfs[0][dfs[0]['gbn'] == gbn]['med_low'].values[0])
                    gy_high.append(dfs[0][dfs[0]['gbn'] == gbn]['med_high'].values[0] - hdx.GetBinContent(i))
                else:
                    gy_low.append(hdx.GetBinContent(i))
                    gy_high.append(hdx.GetBinContent(i))

            gdx = ROOT.TGraphAsymmErrors(len(gx), np.array(gx, dtype='float'), np.array(gy, dtype='float'), np.array(gx_low, dtype='float'), np.array(gx_high, dtype='float'), np.array(gy_low, dtype='float'), np.array(gy_high, dtype='float'))
            gdx.SetName('gdx')
            dxfit = ROOT.TF1('dxfit','pol3',hreco.GetXaxis().GetBinLowEdge(1), hreco.GetXaxis().GetBinLowEdge(hreco.GetNbinsX() + 1)); # fit function
            gdx.Fit('dxfit', 'MRQFEX0')
            gdx.Write()

            # sim
            if len(dfs) > 1:
                hdx_sim = hsim.ProjectionX('hdx_sim', j, j, k, k, 'o')
                gx_sim, gy_sim = [], []
                for i in range(hsim.GetNbinsX() + 2):
                    gbn = hsim.FindBin(hdx_sim.GetXaxis().GetBinCenter(i), hsim.GetYaxis().GetBinCenter(j), hsim.GetZaxis().GetBinCenter(k))
                    gx_sim.append(hdx_sim.GetXaxis().GetBinCenter(i))
                    gy_sim.append(hdx_sim.GetBinContent(i))
                gdx_sim = ROOT.TGraph(len(gx_sim), np.array(gx_sim, dtype='float'), np.array(gy_sim, dtype='float'))
                gdx_sim.SetName('gdx_sim')
                gdx_sim.Write()
               

            ##################################################
            ### Local drift velocity
            ##################################################
            tempx = np.linspace(hreco.GetXaxis().GetBinLowEdge(1), hreco.GetXaxis().GetBinLowEdge(hreco.GetNbinsX() + 1), 100)
            tempy = [dxfit.Eval(x) for x in tempx]
            
            # Calculate difference between adjacent dxs
            ddxs = np.diff(tempy)
            spacing = np.diff(tempx)[0]
            loc_v = [0.157565 * (1 - ddx / spacing) for ddx in ddxs]
            xadj = [x + spacing/2 for x in tempx]
            del xadj[-1]
            
            gvel = ROOT.TGraph(len(xadj), np.array(xadj, dtype='float'), np.array(loc_v, dtype='float'))
            gvel.SetName('gvel')
            gvel.Write()

            ##################################################
            ### Local efield
            ##################################################
            # Fit for efield using know relationship with drift velocity
            driftvelfit = ROOT.TF1('driftvelfit','(1.0 - 0.0184 * (87.5 - 89.0)) * pol5', -0.05, 1.1)
            driftvelfit.SetParameter(0, 0.0)
            driftvelfit.SetParameter(1, 5.53416)
            driftvelfit.SetParameter(2, -6.53093)
            driftvelfit.SetParameter(3, 3.20752)
            driftvelfit.SetParameter(4, 0.389696)
            driftvelfit.SetParameter(5, -0.556184)
            loc_e = [driftvelfit.GetX(10*v) for v in loc_v]
            
            # Efield vs. x
            gefield = ROOT.TGraph(len(xadj), np.array(xadj, dtype='float'), np.array(loc_e, dtype='float'))
            # Fit curve to efield, so that we can take its derivative to find rho
            efit = ROOT.TF1('efit', 'pol3', hreco.GetXaxis().GetBinLowEdge(1), hreco.GetXaxis().GetBinLowEdge(hreco.GetNbinsX() + 1))
            gefield.Fit("efit", "MRQF");
            gefield.SetName('gefield')
            gefield.Write()

            ##################################################
            ### Local space charge density
            ##################################################
            rho = [10**9 * 0.00013316416 * efit.Derivative(x) for x in xadj]
            
            # rho vs. x
            grho = ROOT.TGraph(len(xadj), np.array(xadj, dtype='float'), np.array(rho, dtype='float'))
            grho.SetName('grho')
            grho.Write()

            
    
    # Plot
    yz_bin = '6_5'
    gdxs = [outfile.Get(yz_bin + '/gdx'), outfile.Get(yz_bin + '/gdx_sim')]
    plot_yz_bin_dx_vs_x(gdxs)

    histfile.Close()
    outfile.Close()




def plot_yz_bin_dx_vs_x(gdxs):
    c = ROOT.TCanvas('c', '', 800, 800)
    c.cd()

    c.SetGrid();
    #mg = ROOT.TMultiGraph('mg', '')
    
    # reco
    g0 = gdxs[0]
    g0.GetFunction('dxfit').Delete()
    g0.SetTitle('Reco')
    g0.SetMarkerStyle(21)
    g0.SetMarkerColor(39)
    g0.GetHistogram().SetMinimum(-0.1)
    g0.GetHistogram().SetMaximum(1.0)
    g0.Draw('AP')
    #mg.Add(g0)
    
    legend = ROOT.TLegend(0.55,0.65,0.76,0.82);
    legend.SetTextSize(0.03)
    legend.AddEntry(g0, 'Reco', 'lp')

    # sim
    if len(gdxs) > 1:
        g1 = gdxs[1]
        g1.SetTitle('Simulation')
        g1.SetDrawOption('L')
        g1.Draw('SAME')
        legend.AddEntry(g1, 'Sim', 'lp')

    legend.Draw()

    c.SaveAs('test.png')            


  
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--in_file_hist', type=str, required=True)
    parser.add_argument('--in_file_reco', type=str, required=True)
    parser.add_argument('--in_file_sim', type=str, required=False)
    parser.add_argument('--yz_bin', type=str, required=False)
    args = parser.parse_args()
    main(args)
