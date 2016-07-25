#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooDataHist.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "TRandom.h"
#include <stdio.h>
#include "TMath.h"

using namespace RooFit ;

// Script uses a previously made .root file containing a histogram and
// fits an exponentially decaying function to it. Script then plots (and
// saves) the data and fit in .pfd format.


void lifetimeFit()
{
  gROOT->SetBatch(kTRUE) ;
  // Import file.
  //TFile f("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/histoTAU_Xi_cplus_SigOnly_cut05_20bins_snglgaus.root") ;
  //TFile f("~/Documents/uni/LHCb_CharmSummerProj/learning_root/histoTAU_Xi_cplus_SigOnly_20bins_IPCHI3_PM50_snglgaus.root") ;
TFile f("~/Documents/uni/LHCb_CharmSummerProj/learning_root/histoTAU_Xi_cplus_SigOnly_20bins_IPCHI3_PM70_snglgaus.root") ;
  
  TH1F *h_signal = (TH1F*)f.Get("h_Signal");  // gets the histogram of signal yields from the .root file


  // Read histogram into RooFit dataset
  RooRealVar tau("tau", "tau (ns)", 0.00025, 0.002) ;
  RooDataHist binnedData("binnedData", "binnedData", RooArgSet(tau), h_signal) ;



  // Build exponential PDF.
  RooRealVar expoParam("expoParam", "expoParam", -5000., -5000., 0.) ;
  RooRealVar blindConst("blindConst", "blindConst", -5000., -5000., 0.) ;

  RooUnblindOffset decayConstUnblind("decayConstUnblind", "Unblind decay rate", "xiblindingstring", 500, blindConst) ;

  //RooExponential lifetimePDF("lifetimePDF", "lifetimePDF", tau, expoParam) ; 
  RooExponential lifetimePDF("lifetimePDF", "lifetimePDF", tau, decayConstUnblind) ;


  // Fit to data
  lifetimePDF.fitTo(binnedData, RooFit::SumW2Error(false)) ;


  // print the blinded decay rate:
  cout << "Decay rate: " << blindConst.getVal() << " +/- " << blindConst.getError() << endl ;


  // Plot
  TCanvas* c1 = new TCanvas("taufit","taufit",900,600) ;
  //TCanvas c1 ;
  RooPlot* lifetimePlot = tau.frame(Title(" ")) ;
  binnedData.plotOn(lifetimePlot) ;
  lifetimePDF.plotOn(lifetimePlot);
  lifetimePlot->GetYaxis()->SetTitleOffset(1.4);
  lifetimePlot->SetXTitle("time #font[12]{t} (#font[12]{ns})") ;
  lifetimePlot->SetYTitle("Entries per bin (bin size: 8.75e-5 ns)") ;
  lifetimePlot->Draw();  
  cout<<"!!!!!!!!!!  chi2 = "<< lifetimePlot->chiSquare() << endl ;

  TLatex Tl ;
  Tl.SetNDC() ;
  Tl.SetTextSize(0.04) ;
  Tl.DrawLatex(0.55, 0.75, "#scale[1.1]{#tau_{#Xi_{c}^{+}}} = 0.525362  #pm 0.027682 #font[12]{ps}") ;
  Tl.DrawLatex(0.55, 0.67, "#chi^{2} = 0.7872") ;
  // Save file as .pdf
  c1->SaveAs("histoTAU_Xi_cplus_SigOnly_IPCHI3_PM70_20bins_snglgausFITTED.pdf") ;
 
  //______________________________________________________
  
  RooHist* hresid = lifetimePlot->residHist() ;
  RooHist* hpull = lifetimePlot->pullHist() ;
  
  RooPlot* frame2 = tau.frame(Title("Residual Distribution")) ;
  frame2->addPlotable(hresid,"P") ;
  
  RooPlot* frame3 = tau.frame(Title(" ")) ;
  frame3->addPlotable(hpull,"P") ;
  
  TCanvas* c2 = new TCanvas("chi2residuals", "chi2residuals", 900,600) ;
  //c2->Divide(2) ;
  //c2->cd(1) ; 
  gPad->SetGrid(0,1) ;
  //gPad->SetLeftMargin(0.15) ; 
  frame3->SetXTitle("time #it{t} (#it{ns})") ;
  frame3->SetYTitle("Pull (#tau_{#Xi_{c}^{+}})") ;
  frame3->GetYaxis()->SetTitleOffset(1.4) ; 
  frame3->Draw() ;
  //c2->cd(2) ; gPad->SetLeftMargin(0.15) ; frame3->GetYaxis()->SetTitleOffset(1.6) ; frame3->Draw() ;

  c2->SaveAs("res_pull_xiTAU_IPCHI3_PM70_20bins_snglgaus.pdf") ;
  //______________________________________________________
  
}



