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

  // Import file.  
  //TFile f("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/histoTAU_Lambda_cplus_SigOnly_cut04_20bins.root") ;
  TFile f("~/Documents/uni/LHCb_CharmSummerProj/learning_root/histoTAU_Lambda_cplus_SigOnly_IPCHI3_PM70_20bins.root") ;  

  TH1F *h_signal = (TH1F*)f.Get("h_Signal"); // gets the histogram of signal yields from the .root file



  // Read histogram into RooFit dataset
  RooRealVar tau("tau", "tau (ns)", 0.00025, 0.002) ;
  RooDataHist binnedData("binnedData", "binnedData", RooArgSet(tau), h_signal) ;



  // Build exponential PDF.
  //RooRealVar expoParam("expoParam", "expoParam", -5000., -5000., 0.) ; //for unblinded analysis
  RooRealVar blindConst("blindConst", "blindConst", -5000., -5000., 0.) ; //for blinded analysis

  RooUnblindOffset decayConstUnblind("decayConstUnblind", "Unblind decay rate", "someblindingstring", 500, blindConst) ; // to make comparisons with previous work.

  //RooExponential lifetimePDF("lifetimePDF", "lifetimePDF", tau, expoParam) ; // non-blinded PDF 
  RooExponential lifetimePDF("lifetimePDF", "lifetimePDF", tau, decayConstUnblind) ; // blinded PDF



  // Fit to data
  lifetimePDF.fitTo(binnedData, RooFit::SumW2Error(false)) ;


  // Print the blinded decay rate to screen
  cout << "Decay rate: " << blindConst.getVal() << " +/- " << blindConst.getError() << endl ;


  // Plot
  TCanvas c1;
  RooPlot* lifetimePlot = tau.frame() ;
  binnedData.plotOn(lifetimePlot) ;
  lifetimePDF.plotOn(lifetimePlot);
  lifetimePlot->GetYaxis()->SetTitleOffset(1.4);
  lifetimePlot->Draw();  


  // Save file as .pdf
  c1.SaveAs("histoTAU_Lambda_cplus_SigOnly_IPCHI3_PM70_20binsFITTED.pdf") ;

}



