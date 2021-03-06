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
#include "TH1F.h"
#include "TH1.h"
#include "TRandom.h"
#include "RooFitResult.h"
#include "RooAddPdf.h"
#include "RooExponential.h"
#include "RooPolynomial.h"
#include "RooAbsData.h"
#include "TLatex.h"
#include "TFile.h"
#include "TLegend.h"

using namespace RooFit ;


// Script opens a previously made .root file which contains a dataset type object which 
// holds all events in preselected and appropriate ranges for variables of the 
// Lambda_cplus particle.


void binFit() {
  gROOT->SetBatch(kTRUE); // stops "visial confirmation plots" being shown onscreen.
  // Open appropriate .root file.
  //TFile *datafile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/Gedcode/baryon-lifetimes-2015/data/run-II-data/datafileLambda_TAUmin200fs_max2200fs_Mmin2216_max2356.root");
  //TFile *datafile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/Gedcode/baryon-lifetimes-2015/data/run-II-data/datafileLambda_TAUmin200fs_max2200fs_Mmin2216_max2356_CutIPCHI2lt3.root"); 
  //TFile *datafile = TFile::Open("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/TMVA_cut01_firsttry.root") ;
TFile *datafile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/learning_root/TMVA_cut01_firsttry.root") ;


  // Define and fill Tree.
  //TTree* mytree = (TTree*) datafile->Get("Lambda_cToKppiTuple/DecayTree;1");
  TTree* mytree = (TTree*) datafile->Get("TrainTree");


  // Variables and their limits.
  RooRealVar Lambda_cplus_TAU("Lambda_cplus_TAU","Lambda_cplus_TAU",0.00025 ,0.002 ,"ns") ;
  double highestTAU;
  double lowestTAU;

  RooRealVar Lambda_cplus_M("Lambda_cplus_M","Lambda_cplus_M",2216 ,2356, "MeV") ; // +/-70MeV  //CHANGE
  //RooRealVar Lambda_cplus_M("Lambda_cplus_M","Lambda_cplus_M",2236 ,2336, "MeV") ; // +/-50MeV 
  double highestM;
  double lowestM;

  RooRealVar Lambda_cplus_IPCHI2_OWNPV("Lambda_cplus_IPCHI2_OWNPV","Lambda_cplus_IPCHI2_OWNPV",-100 ,100) ; 

  RooRealVar BDTG("BDTG","BDTG",-1. ,1.) ; 


  // Define dataset and make cuts.
  //RooDataSet* ds = (RooDataSet*)datafile->Get("ds") ;
  RooDataSet ds("ds","ds",RooArgSet(Lambda_cplus_TAU, Lambda_cplus_M, Lambda_cplus_IPCHI2_OWNPV, BDTG),Import(*mytree),Cut("(0.00025<Lambda_cplus_TAU)&&(Lambda_cplus_TAU<0.002)&&(Lambda_cplus_M<2356)&&(2216<Lambda_cplus_M)&&(0.0107<BDTG)&&(Lambda_cplus_IPCHI2_OWNPV<3)")) ; // +/-70MeV  //CHANGE
  //RooDataSet ds("ds","ds",RooArgSet(Lambda_cplus_TAU, Lambda_cplus_M, Lambda_cplus_IPCHI2_OWNPV, BDTG),Import(*mytree),Cut("(0.00025<Lambda_cplus_TAU)&&(Lambda_cplus_TAU<0.002)&&(Lambda_cplus_M<2336)&&(2236<Lambda_cplus_M)&&(0.0107<BDTG)&&(Lambda_cplus_IPCHI2_OWNPV<3)")) ; // +/-50MeV  //CHANGE


  // Build probability density functions (PDFs).
  // Combined double Gaussian PDF.

  ds.RooAbsData::getRange(Lambda_cplus_TAU, lowestTAU, highestTAU);
  ds.RooAbsData::getRange(Lambda_cplus_M, lowestM, highestM) ;
  double mass_peak = 2286 ;  // initial mean value for fit (from theory)

  RooRealVar gaus_mean1("gaus_mean1", "gaus_mean1",mass_peak, lowestM, highestM, "MeV") ;
  RooRealVar sigma1("sigma1","sigma1", 6, 0, 50) ;
  RooGaussian gauss1("gauss1","gauss1",Lambda_cplus_M, gaus_mean1, sigma1) ;

  RooRealVar gaus_mean2("gaus_mean2", "gaus_mean2",mass_peak, lowestM, highestM, "MeV") ;
  RooRealVar sigma2("sigma2","sigma2", 4, 0, 100) ;
  RooGaussian gauss2("gauss2","gauss2",Lambda_cplus_M, gaus_mean2, sigma2) ;

  RooRealVar gaus1_fraction("gaus1_fraction", "gaus1_fraction", 0.5, 0.0, 1.) ;
  RooAddPdf gaussComb("gaussComb","gaussComb", RooArgList(gauss1, gauss2), RooArgList(gaus1_fraction)) ;

  //Lambda_cplus_TAU.setRange("R1",0.00018, 0.0012);
   
  // Exponential PDF.
  RooRealVar exp_const("exp_const","exp_const", -0.0001, -1., 0.);
  RooExponential expo_bkg("expo_bkg", "expo_bkg", Lambda_cplus_M, exp_const);

  // Model PDF.
  RooRealVar nSignal("nSignal","nSignal", 23000, 0, 409570);
  RooRealVar nBkg("nBkg","nBkg", 13000, 0, 409570);
  
  //RooAddPdf model("model","model",RooArgList(gauss, pol0),RooArgList(nSignal, pol0_yield));
  //RooAddPdf model("model","model",RooArgList(gauss, expo_bkg),RooArgList(nSignal, nBkg));
  RooAddPdf model("model","model",RooArgList(gaussComb, expo_bkg),RooArgList(nSignal, nBkg));



  // Fit model PDF to data.
  //model.fitTo(*ds, Range("R1"));
  //expo1.fitTo(*ds, Range("R1"));
  model.fitTo(ds, Extended()) ;
  //RooFitResult* rmodel = model.fitTo(*ds,Save()) ;
  //RooFitResult* rexpo_bkg = expo_bkg.fitTo(*ds,Save());




  // Plot the data and the fit (not necesary, but provides a vidual confirmation that
  // things 'look good'.)

  RooPlot *fullDataFit = Lambda_cplus_M.frame(Title(" "));
  params = RooArgSet(gaus1_mean, gaus2_mean, sigma1, sigma2, exp_const) ;
  //ds.plotOn(frame,Binning(25)); //default is 100 bins
  ds.plotOn(fullDataFit, Name("data"), MarkerColor(kBlack)) ;
  ds.statOn(fullDataFit, Layout(0.65,0.88,0.9), What("N")) ; //NB Layout(xmin,xmax,ymax)
  model.plotOn(fullDataFit, Name("Model"), DrawOption("L")) ;
  fullDataFit->GetYaxis()->SetTitleOffset(1.4);
  fullDataFit->SetXTitle("mass (MeV)") ;
  //  fullDataFit->SetYTitle("Entries per bin (bin size: 1 MeV)") ; // +/-50MeV //CHANGE
  fullDataFit->SetYTitle("Entries per bin (bin size: 1.4 MeV)") ;   // +/-70MeV

  //=======Jump to 'here'=============



  // Plot the residuals and Pulls...
  //______________________________________________________

  RooHist* hresid = fullDataFit->residHist() ;
  RooHist* hpull = fullDataFit->pullHist() ;

  //  RooPlot* frame2 = Lambda_cplus_M.frame(Title("Residual Distribution")) ;
  //frame2->addPlotable(hresid,"P") ;

  RooPlot* frame3 = Lambda_cplus_M.frame(Title(" ")) ;
  frame3->addPlotable(hpull,"P") ;



  TCanvas* c2 = new TCanvas("chi2residuals", "chi2residuals", 900,600) ;
  gPad->SetGrid(0,1) ;
  //c2->Divide(2) ;
  //c2->cd(1) ; gPad->SetLeftMargin(0.15) ; frame2->GetYaxis()->SetTitleOffset(1.6) ; frame2->Draw() ;
  //c2->cd(2) ;
  //gPad->SetLeftMargin(0.15) ;
  frame3->SetXTitle("mass (MeV)") ;
  frame3->SetYTitle("Mass Pull (#Lambda_{c}^{+})") ;
  frame3->GetYaxis()->SetTitleOffset(1.4) ;
  frame3->Draw() ;
  //c2->SaveAs("pull_LambdaM_IPCHI3_PM50_20bins_dblgaus.pdf") ;  //CHANGE
  c2->SaveAs("pull_LambdaM_IPCHI3_PM70_20bins_dblgaus.pdf") ;
  //______________________________________________________


  //========='here'==========


  model.plotOn(fullDataFit, Components(expo_bkg), LineStyle(kDashed)) ;
  model.paramOn(fullDataFit,Layout(0.2, 0.5, 0.88)) ; //was 0.4
  fullDataFit->getAttText()->SetTextSize(0.022) ;

  /*
  RooDataHist hist4Chi2("hist4Chi2","hist4Chi2", RooArgSet(Lambda_cplus_M), ds) ;
  Double_t chi2 = fullDataFit->chiSquare("Model","data",7) ;
  */

  // Set model fit variables to constants (NOT COEFFs!)
  gaus_mean1.setConstant() ;
  gaus_mean2.setConstant() ;
  sigma1.setConstant() ;
  sigma2.setConstant() ;
  gaus1_fraction.setConstant() ;
  exp_const.setConstant() ;
 


  // Create array of signal candidates in TAU variable.
  int nBins=20 ;      // If you change this, rememeber to change (*) array sizes below AND (**) h_Signal histogram!
  


  // Create arrays to store signal yield and error for each bin. 
  double signalYield[20] ;  // (*)
  double signalError[20] ;
  double binCent[20] ;      // to record centre of bin.

  // Split data into 'TAU' bins, make mass fits to individual bins and store signal 
  // yields in array.


  int i ;
  for (i=0; i<nBins; i++){
 
    double binBoundaryLo = 0.00025 + (i * ((0.002-0.00025)/nBins)) ;
    double binBoundaryHi = 0.00025 + ((i+1) * ((0.002-0.00025)/nBins)) ;

    std::ostringstream binLowEdgeStr;
    binLowEdgeStr << binBoundaryLo;
    std::ostringstream binHighEdgeStr;
    binHighEdgeStr << binBoundaryHi;

    // Creat dataset for bin i only
    RooDataSet* bindata = ds.reduce(RooFit::Cut(TString(binLowEdgeStr.str()) + " < Lambda_cplus_TAU && Lambda_cplus_TAU <= " + TString(binHighEdgeStr.str()))) ;

    //model.fitTo(*bindata, gErrorIgnoreLevel = kInfo) ;
    model.fitTo(*bindata, "l") ;
    double binSignalYield = nSignal.getVal() ;
    double binError = nSignal.getError() ;
    double binCentre = (binBoundaryLo + binBoundaryHi)/2. ;

    printf("binSignalYield = %f", binSignalYield) ; 
    puts("") ;
    printf("binError = %f", binError) ; 
    puts("") ;
    printf("binCentre = %f", binCentre) ; 
    puts("") ;
    cout<<i<<endl ;

    // Fill arrays
    signalYield[i] = binSignalYield ;
    signalError[i] = binError ;
    binCent[i] = binCentre ;

}

  TH1D *h_Signal = new TH1D("h_Signal","h_Signal",20 ,0.00025 ,0.002) ; // (**)
  int j ;
    
  // Fill histogram
  for (j=0; j<nBins; j++){
  h_Signal->SetBinContent(j+1, signalYield[j]) ;
  h_Signal->SetBinError(j+1, signalError[j]) ;
  cout<<"bin"<<j+1<<endl;
  cout<<"yield:"<<h_Signal->GetBinContent(j+1)<<" error:"<<h_Signal->GetBinError(j+1)<<endl;
}

  // Create and fill histogram with signal yields and errors calculated 
  // from the fits made in the 'for' loop, above.

  //  TFile hf("/afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/histoTAU_Lambda_cplus_SigOnly_cut04_20bins.root", "RECREATE") ;
  TFile hf("~/Documents/uni/LHCb_CharmSummerProj/learning_root/histoTAU_Lambda_cplus_SigOnly_IPCHI3_PM70_20bins_dblgaus.root", "RECREATE") ; //CHANGE
  //TFile hf("~/Documents/uni/LHCb_CharmSummerProj/learning_root/histoTAU_Lambda_cplus_SigOnly_IPCHI3_PM50_20bins_dblgaus.root", "RECREATE") ;
  h_Signal->Write() ;


  // Plot the full sig & bkg fit (after cuts)
  TCanvas *c4 = new TCanvas(" ", " ",900, 600) ;
  fullDataFit->Draw() ;
  //c4->SaveAs("FullFit_Xi_IPCHI3_PM50_snglgaus.pdf") ; //CHANGE
  c4->SaveAs("FullFit_Xi_IPCHI3_PM70_snglgaus.pdf") ;


  //=================NON-ESSENTIAL TO ANALYSIS=================
  // This section draws the overall fit (to sig and background (after cuts) for all data) and draws the signal-only 
  // TAU variable hitogram with a exponential fit. This is simple for a visual confirmation that 
  // initial fit is sensible, and that output histo is exponential in form. 

  RooDataHist vis_binnedData("vis_binnedData", "vis_binnedData", RooArgSet(Lambda_cplus_TAU),h_Signal) ;  
  RooRealVar vis_expo("vis_expo","exponential parameter for visual confirmation", -5000.,-5000., 0.) ;
  RooExponential vis_lifetimePDF("vis_lifetimePDF","vis_lifetimePDF", Lambda_cplus_TAU, vis_expo) ;
  vis_lifetimePDF.fitTo(vis_binnedData) ;
  RooPlot* vis_lifetimePlot = Lambda_cplus_TAU.frame() ;
  vis_binnedData.plotOn(vis_lifetimePlot) ;
  vis_lifetimePDF.plotOn(vis_lifetimePlot) ;


  TCanvas *c102 = new TCanvas("c102","",600,900) ;
  formatCanvas4(c102) ;

  c102_1->cd() ;
  fullDataFit->Draw() ;

  c102_2->cd() ;
  h_Signal->Draw() ;

  c102_3->cd() ;
  vis_lifetimePlot->Draw() ;

  c102->Update() ;
  c102->SaveAs("visConf_Lambda_20bins_IPCHI3_PM70_dblgaus.pdf") ; //CHANGE
  //c102->SaveAs("visConf_Lambda_20bins_IPCHI3_PM50_dblgaus.pdf") ;

  // Print useful info to screen

  cout<<endl<<"   ************info************"<<endl<<endl;
  cout<<"  Fit Quality>>> chi2 = "<< fullDataFit->chiSquare() << endl ;
  cout<<"Gaussian fit parameters to full data (signal + background):"<<endl;
  cout<<"  "<<gaus_mean1<<endl ;
  cout<<"  "<<gaus_mean2<<endl ;
  cout<<"  "<<sigma1<<endl ;
  cout<<"  "<<sigma2<<endl ;
  cout<<"  "<<gaus1_fraction<<endl ;
  cout<<"Exponential fit parameters to full data (signal + background):"<<endl;
  cout<<"  "<<exp_const<<endl<<endl ;



  int k ;
  for (k=0; k<nBins; k++){
    cout<<"Bin number "<<k+1<<" (bin centre = "<<binCent[k]*1000<<"ps):"<<endl;
    cout<<"   "<<"signal yield = "<<signalYield[k]<<endl ;
    cout<<"   "<<"signal error = "<<signalError[k]<<endl ;
  }
 

    /*
        // Chi^2___________________
	Lambda_cplus_M.setBins(100);
	RooDataHist hist("hist","hist", RooArgSet(Lambda_cplus_M), *ds) ;
	RooChi2Var chi2("chi2","chi2", model, hist, true) ;
	RooMinuit minuit(chi2);
	minuit.migrad();
	minuit.hesse();
    */
}


// Canvas formatting (4 plots on landscape A4)
void formatCanvas4(TCanvas* canvas){
  canvas->SetBorderSize(2);
  canvas->SetBorderMode(-1);
  canvas->Clear();
  canvas->Divide(2,2);
  canvas->cd(1);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  canvas->cd(2);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  canvas->cd(3);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
  canvas->cd(4);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.2);
}



void formatHisto(TH1F* hist, const char* xtitle, const char* ytitle, int ndivisionsX=506, int ndivisionsY=506){

  //hist->Sumw2();
  //hist->SetMinimum(0.8);
  //hist->SetMaximum(1.0);
  //hist->SetTitle("");
  //hist->SetStats(0010);
  hist->SetNdivisions(ndivisionsX,"X");
  hist->SetNdivisions(ndivisionsY,"Y");
  hist->GetXaxis()->SetTitle(xtitle);
  hist->GetYaxis()->SetTitle(ytitle);
  hist->GetXaxis()->SetTitleColor(kBlack);
  hist->GetYaxis()->SetTitleColor(kBlack);
  hist->GetXaxis()->SetTitleOffset(1.3);//0.9
  hist->GetYaxis()->SetTitleOffset(1.7);//0.9
  hist->GetXaxis()->SetLabelOffset(0.02);//0.065
  hist->GetYaxis()->SetLabelOffset(0.02);//0.065
  hist->GetXaxis()->SetTitleSize(0.07);//0.065
  hist->GetYaxis()->SetTitleSize(0.06);//0.065
  hist->GetXaxis()->SetLabelSize(0.055);//0.065
  hist->GetYaxis()->SetLabelSize(0.055);//0.065
  hist->SetFillStyle(3005);
  hist->SetFillColor(kBlack);
  hist->SetLineColor(kBlack);
  hist->SetLineWidth(1.1);
}
