#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
using namespace RooFit ;


void basic_plotter()
{

  // S e t u p   m o d e l 
  // ---------------------

  // Create observables
  RooRealVar Lambda_cplus_M("Lambda_cplus_M","Lambda_cplus_M",2220,2540) ;
  RooRealVar Lambda_cplus_TAU("Lambda_cplus_TAU", "Lambda_cplus_TAU",0.00025,0.002) ;
  // Create Gaussian
  RooRealVar sigma("sigma","sigma",1,0.1,10) ;
  RooRealVar mean("mean","mean",-3,-10,10) ;
  RooGaussian gauss("gauss","gauss",Lambda_cplus_M,mean,sigma) ;

  // Generate a sample of 1000 events with sigma=3
  //RooDataSet* data = gauss.generate(Lambda_cplus_M,1000) ;


  TFile *datafile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/learning_root/turbo_2015_data.root");
 
  TTree* mytree = (TTree*) datafile->Get("Lambda_cToKppiTuple/DecayTree;1") ;
  // Define dataset
  //RooDataSet* ds = (RooDataSet*) mytree->Get("ds") ;
  //RooDataSet data("data","data",RooArgSet(Lambda_cplus_M),Import(*mytree), Cut("(0.00025<Lambda_cplus_TAU)&&(Lambda_cplus_TAU<0.002)&&(Lambda_cplus_M<2540)&&(2220<Lambda_cplus_M)&&(Lambda_cplus_IPCHI2_OWNPV<3)")) ;
  RooDataSet data("data","data",RooArgSet(Lambda_cplus_M),Import(*mytree), Cut("(Lambda_cplus_M<2540)&&(2220<Lambda_cplus_M)")) ;


  // Fit pdf to data
  //gauss.fitTo(*data) ;


  // P l o t   p . d . f   a n d   d a t a 
  // -------------------------------------

  // Overlay projection of gauss on data
  RooPlot* frame = Lambda_cplus_M.frame(Name("Lambda_cplus_Mframe"),Title(" "),Bins(100)) ;
  data.plotOn(frame) ;
  //data.statOn(frame, What("N"), Layout(0.55, 0.8,0.8)) ; 
  data.statOn(frame, What("N"), Layout(0.7, 0.94,0.9)) ; 
 //gauss.plotOn(frame) ;
  frame->SetXTitle("Mass (MeV)") ;
  frame->SetYTitle("Entries per bin") ;

  // A d d   b o x   w i t h   p d f   p a r a m e t e r s 
  // -----------------------------------------------------

  // Left edge of box starts at 55% of Xaxis)
  //gauss.paramOn(frame,Layout(0.55)) ;


  // A d d   b o x   w i t h   d a t a   s t a t i s t i c s
  // -------------------------------------------------------  

  // X size of box is from 55% to 99% of Xaxis range, top of box is at 80% of Yaxis range)
  //data.statOn(frame,Layout(0.55,0.99,0.8)) ;


  // A d d   t e x t   a n d   a r r o w 
  // -----------------------------------

  // Add text to frame
  //TText* txt = new TText(2400,45000,"") ;
  //txt->SetTextSize(0.04) ;
  //txt->SetTextColor(kRed) ;
  //frame->addObject(txt) ;

  // Add Xic arrow to frame
  TArrow* Xic_arrow = new TArrow(2435,20000,2465,14000,0.02,"|>") ;
  Xic_arrow->SetLineColor(kRed) ;
  Xic_arrow->SetFillColor(kRed) ;
  Xic_arrow->SetLineWidth(3) ;
  frame->addObject(Xic_arrow) ;

  
  // Add Lambda_cplus arrow to frame
  TArrow* Lambda_arrow = new TArrow(2350,26000,2300,22500,0.02,"|>") ;
  Lambda_arrow->SetLineColor(kRed) ;
  Lambda_arrow->SetFillColor(kRed) ;
  Lambda_arrow->SetLineWidth(3) ;
  frame->addObject(Lambda_arrow) ;
  
  //TFile f("rf106_plotdecoration.root","RECREATE") ;
  //frame->Write() ;
  //f.Close() ;

  // To read back and plot frame with all decorations in clean root session do
  // root> TFile f("rf106_plotdecoration.root") ;
  // root>  xframe->Draw() ;

TCanvas *c1 = new TCanvas("plotwithdecoration","plotwithdecoration",900,600) ;
  gPad->SetLeftMargin(0.15) ; frame->GetYaxis()->SetTitleOffset(1.6) ; frame->Draw() ;
  gPad->SetGrid(1,0) ;
  c1->GetXaxis()->SetAxisColor(17);
  c1->GetYaxis()->SetAxisColor(17);



  TLatex Tl ;
  Tl.SetNDC() ;
  Tl.SetTextSize(0.05) ;
  Tl.DrawLatex(0.6, 0.6, "#scale[1.1]{#Xi^{+}_{c}} resonance") ;
  Tl.DrawLatex(0.45, 0.75, "#scale[1.1]{#Lambda^{+}_{c}} resonance") ;
  //Tl.DrawLatex(0.69, 0.8, "#scale[1.1]{#Lambda^{+}_{c}} (udc)") ;
  //Tl.DrawLatex(0.69, 0.6, "mass resonance") ;
  cout<<"saving..."<<endl ;  
c1->SaveAs("FullData_wArrows01.pdf") ;
}
