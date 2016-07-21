#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "RooPlot.h"
#include "TTree.h"
#include "TH1D.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooPolynomial.h"

using namespace std ;
using namespace RooFit ;


// Script opens a previously made ROOT file which contains a dataset type object (which is a reduced version 
// of the full data, reduced but applying specific cuts in specfic variables) and calculates 'sWeights'.
// sWeights are then saved to a 'dataset' .root file AND/OR to a TTree type ROOT file, along with their
// original data.


void sWeigher() {



  //TFile *datafile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/learning_root/turbo_M_TAU_cut01.root") ; // opens 'reduced' data.
  TFile *datafile = TFile::Open("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/turbo_M_TAU_cut01.root") ;  


  // Define dataset
  RooDataSet* ds = (RooDataSet*)datafile->Get("ds") ;


  // Define TAU variable, get limits.
  RooRealVar Lambda_cplus_TAU("Lambda_cplus_TAU","Lambda_cplus_TAU",0.00025 ,0.002 ,"ns") ; 
  double highestTAU;
  double lowestTAU;
  ds->RooAbsData::getRange(Lambda_cplus_TAU, lowestTAU, highestTAU);

  // Define Mass variable, get limits.
  RooRealVar Lambda_cplus_M("Lambda_cplus_M","Lambda_cplus_M",2216 ,2356, "GeV") ; 
  double highestM;
  double lowestM;
  ds->RooAbsData::getRange(Lambda_cplus_M, lowestM, highestM) ;
  
  // Define IPCHI2 variable
  //RooRealVar Lambda_cplus_IPCHI2_OWNPV("Lambda_cplus_IPCHI2_OWNPV","Lambda_cplus_IPCHI2_OWNPV",-100 ,100) ; 


  // Build combined double Gaussian PDF; "gaussComb". 
  double mass_peak = 2286 ;  // initialize gaussian mean value for fit for Lambda_cplus

  RooRealVar gausMean1("gausMean1", "gausMean1",mass_peak, lowestM, highestM, "GeV") ;
  RooRealVar sigma1("sigma1","sigma1", 6, 0, 50) ;
  RooGaussian gauss1("gauss1","gauss1",Lambda_cplus_M, gausMean1, sigma1) ;

  RooRealVar gausMean2("gausMean2", "gausMean2",mass_peak, lowestM, highestM, "GeV") ;
  RooRealVar sigma2("sigma2","sigma2", 4, 0, 100) ;
  RooGaussian gauss2("gauss2","gauss2",Lambda_cplus_M, gausMean2, sigma2) ;

  RooRealVar nFrac("nFrac", "nFrac", 0.5, 0.0, 1.) ;
  RooAddPdf gaussComb("gaussComb","gaussComb", RooArgList(gauss1, gauss2), RooArgList(nFrac)) ;


   
  // Build exponential PDF; "expo_bkg"
  RooRealVar expoPar("expoPar","expoPar", -0.0001, -1., 0.);
  RooExponential expo_bkg("expo_bkg", "expo_bkg", Lambda_cplus_M, expoPar);


  // Build model PDF
  RooRealVar nSignal("nSignal","nSignal", 200000, 0, 409570); // For CHI2 < 3 cut, the max yield is 147000
  RooRealVar nBkg("nBkg","nBkg", 200000, 0, 409570);          // For CHI2 < 3 cut, the max yield is 147000
  
  //RooAddPdf model("model","model",RooArgList(gauss, pol0),RooArgList(nSignal, pol0_yield));
  //RooAddPdf model("model","model",RooArgList(gauss, expo_bkg),RooArgList(nSignal, nBkg));
  RooAddPdf model("model","model",RooArgList(gaussComb, expo_bkg),RooArgList(nSignal, nBkg));



  // Fit model
  model.fitTo(*ds, Extended()) ;
  params = RooArgSet(gausMean1, gausMean2, sigma1, sigma2, nFrac) ;


  // Plot
  RooPlot *fullDataFit = Lambda_cplus_M.frame(Title(" ")) ;
  //ds->plotOn(fullDataFit, Binning(50)) ; //default is 100 bins
  ds->plotOn(fullDataFit, Name("data"), MarkerColor(kBlack)) ;
  ds->statOn(fullDataFit, Layout(0.65,0.88,0.2), What("N")) ; //NB Layout(xmin,xmax,ymax)
  model.plotOn(fullDataFit, Name("Model"), DrawOption("L")) ;
  model.plotOn(fullDataFit, Components(expo_bkg), LineStyle(kDashed)) ;
  model.plotOn(fullDataFit, Components(gaussComb), LineStyle(kDashed)) ;
  //model.paramOn(fullDataFit,Layout(0.19, 0.45, 0.88)) ;
  //model.paramOn(fullDataFit,Layout(0.19, 0.4, 0.88)) ;
  //model.paramOn(fullDataFit,Layout(0.18, 0.41, 0.88)) ;
  model.paramOn(fullDataFit,Layout(0.15, 0.41, 0.88), Parameters(params)) ;
  

  fullDataFit->getAttText()->SetTextSize(0.028) ; //was 0.025, then 0.03, then 0.027
  fullDataFit->SetXTitle("Mass (GeV)") ;
  fullDataFit->SetYTitle("Number of Events") ;
  fullDataFit->SetTitleOffset(1.3,"Y") ; //was 0.01, then 0.08, then 1.1
  fullDataFit->SetTitleOffset(1.3,"X") ; //was 1.01, then 1.1 



  // Set model fit variables to constants (NOT COEFFs!) ; this is important for the sPlot function below.
  gausMean1.setConstant() ;
  gausMean2.setConstant() ;
  sigma1.setConstant() ;
  sigma2.setConstant() ;
  nFrac.setConstant() ;
  expoPar.setConstant() ;




  // Find signal and background yields in a sub range of the Lambda_cplus_M total range //
  //************************************************************************************//

  double Ntot = ds->numEntries() ;
  double sigma = sigma1.getVal() ;
  double pm = 3.0*sigma ;
  double mean = gausMean1.getVal() ;
  double R1Low = mean - pm ;
  double R1Hi = mean + pm ;

  //cout<<"!!!!!!!!!!!!!  R1Low:"<<R1Low<<" R1Hi:"<<R1Hi<<endl ;
  //cout<<"!!!!!!!!!!!!!  Ntot:"<<Ntot<<" sigma:"<<sigma<<" pm:"<<pm<<" mean:"<<mean<<endl;

  Lambda_cplus_M.setRange("R1", R1Low, R1Hi) ;
  //Lambda_cplus_M.setRange("Rtotal", 2216, 2356) ;

  massArg = RooArgSet(Lambda_cplus_M) ;
  RooAbsReal* bkg_integral = expo_bkg.createIntegral(massArg, NormSet(massArg), Range("R1"));  
  //RooAbsReal* bkg_integral = expo_bkg.createIntegral(Lambda_cplus_M, Range("R1"));  
  double bkg_integral_value = bkg_integral->getVal() ;
  cout<<"bkg_integral_value:"<<bkg_integral_value<<endl ;

  RooAbsReal* sig_integral = gaussComb.createIntegral(massArg, NormSet(massArg), Range("R1"));  
  double sig_integral_value = sig_integral->getVal() ;
  cout<<"sig_integral_value:"<<sig_integral_value<<endl ;

  //RooAbsReal* totalbkg_integral = expo_bkg.createIntegral(massArg, NormSet(massArg), Range("Rtotal"));  
  //double totalbkg_integral_value = totalbkg_integral->getVal() ;
  //cout<<"totalbkg_integral_value:"<<totalbkg_integral_value<<endl ;
  
  //RooAbsReal* totalsig_integral = gaussComb.createIntegral(massArg, NormSet(massArg), Range("Rtotal"));  
  //double totalsig_integral_value = totalsig_integral->getVal() ;
  //cout<<"totalsig_integral_value:"<<totalsig_integral_value<<endl ;

 
  double R1sig = nSignal.getVal()*sig_integral_value ;
  double R1bkg = nBkg.getVal()*bkg_integral_value ;
  double R1tot = R1sig + R1bkg ; 

  //cout<<"!!!!!!!!!!!!! number of entries in 3sigma range of mean1 is: "<<R1tot<<endl ;
  //cout<<"!!!!!!!!!!!!! ...of which "<<R1sig<<" are signal and "<<R1bkg<<" are background."<<endl ;
  //************************************************************************************//




  
  //RooMsgService::instance().setSilentMode(true);

  // Create and record sWeights
  RooStats::SPlot* sDataM = new RooStats::SPlot("sData", "An SPlot", *ds, &model, RooArgList(nSignal, nBkg)) ;

  /*
  // Save dataset with sWeights
   TFile *dsWithWeights = TFile::Open("DataSet_Lambda_TAUmin200fs_max2200fs_Mmin2216_max2356_CutIPCHI2lt3_withWeights.root","RECREATE") ;
  ds->Write() ;
  dsWithWeights->Close() ;
  cout<<endl<<"  Dataset type file created in current directory..."<<endl ;
  */
  /*
  // Save OUTPUT TTree with sWeights
  TFile treefile("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/sWeightsTree_M_TAU_cut01.root", "recreate") ;
  RooTreeDataStore sTree("sTree", "sTree", *ds->get(0), *ds->store()) ;
  TTree& mystree = sTree.tree() ;
  mystree.SetBranchStatus("Lambda_cplus_M", 0) ; // don't want to save the mass variables again.
  mystree.Write() ;
  treefile.Close() ;
  */

  // Print info to screen
  cout<<endl<<endl<<"========== ROOFIT INFO ==========="<<endl<<endl ;
  cout<<":::Double Gaussian fit parameters to full data (signal + background):"<<endl;
  cout<<gausMean1<<endl ;
  cout<<gausMean2<<endl ;
  cout<<sigma1<<endl ;
  cout<<sigma2<<endl ;
  cout<<nFrac<<endl<<endl ;
  cout<<":::Exponential fit parameters to full data (signal + background):"<<endl;
  cout<<expoPar<<endl<<endl<<endl ;
  cout<<"========== Signal and Background yields INFO ========"<<endl ;
  // Check the values of the yields from the fit and the yields from the sWeight are the same
  cout<<"1) Yield of signal is: "<<nSignal.getVal()<<"; from sWeights it is: "
      <<sDataM->GetYieldFromSWeight("nSignal")<<" ...these should match."<<endl ;
  cout<<"2) Yield of background is: "<<nBkg.getVal()<<"; from sWeights it is: "
      <<sDataM->GetYieldFromSWeight("nBkg")<<" ...these should match."<<endl ;
  cout<<"3) Number of entries in 3*sigma range of mean1 is: "<<R1tot<<endl ;
  cout<<"   ...of which "<<R1sig<<" are signal and "<<R1bkg<<" are background."<<endl ;
  // Check to see if total events is the same..
  //cout<<"nSignal + nBkg = "<<nSignal.getVal()<<" + "<<nBkg.getVal()<<" = "<<nSignal.getVal()+nBkg.getVal()<<endl;
  //cout<<"^ This should match the number of events (below): "<<endl ;
  //ds->Print() ;
  cout<<"========== INFO: END ========"<<endl<<endl<<endl ;

  // Save canvas
  TCanvas *c2 = new TCanvas("c2","",900,600) ;
  //gStyle->SetOptStat("") ; // on in 02
  fullDataFit->Draw() ;
  c2->Print("SignalAndBkg06.pdf");

}

