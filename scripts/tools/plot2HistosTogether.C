#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TFile.h"
using namespace RooFit ;


void plot_together()
{

  // Create observables
  RooRealVar nSignal_sw("nSignal_sw","nSignal_sw",-2.,2.) ;
  RooRealVar nBkg_sw("nBkg_sw", "nBkg_sw",-2.,2.) ;

  TFile *datafile = TFile::Open("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/sWeightsTree_M_TAU_cut01.root");
 
  TTree* mytree = (TTree*) datafile->Get("sTree") ;
  TCanvas *c1 = new TCanvas("twoHistos","twoHistos",900,600) ;
  sTree->Draw("nSignal_sw>>h1") ;
  sTree->Draw("nBkg_sw>>+h1") ;

  h1->SetTitle(" ") ;
  h1->SetXTitle("sWeight") ;
  h1->SetYTitle("Entries per bin") ;

  gStyle->SetOptStat("") ;


  c1->SaveAs("sWeightsPlot.pdf") ;


}
