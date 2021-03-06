#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "TTree.h"


using namespace std ;

// Script will prepare a TTree for merging if the two TTree to be merged share the same 
// events - but where one has some missing entries (due to some cut having been made). 
// The larger file is reduced in this script to contain the same entries by making the 
// same cuts.

 
void stripper() {

  // Open file to be stripped
  //TFile *originalfile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/learning_root/turbo_2015_data.root") ;
  TFile *originalfile = TFile::Open("~/Documents/uni/LHCb_CharmSummerProj/learning_root/turbo_2015_data_mergedWithAngularVars.root") ; 
 
  // Get TTree to be stripped
  //TTree* originaltree = (TTree*) originalfile->Get("Lambda_cToKppiTuple/DecayTree;1") ;
  TTree* originaltree = (TTree*) originalfile->Get("DecayTree;35") ;


  // Create new 'stripped' file
  //TFile newfile("~/Documents/uni/LHCb_CharmSummerProj/learning_root/strippedToFitData.root", "recreate") ;
  TFile newfile("~/Documents/uni/LHCb_CharmSummerProj/learning_root/strippedToFitData_M_TAU_cut01_wAngVars.root", "recreate") ;
  

  // Strip the TTree acording to the appropriate cut
  TTree* fitdatatree = originaltree->CopyTree("0.00025 < Lambda_cplus_TAU && Lambda_cplus_TAU < 0.002 && 2216 < Lambda_cplus_M && Lambda_cplus_M < 2356") ;
 

  // Write file
  fitdatatree->Write() ;
  newfile.Close() ;

  cout<<"========================================================="<<endl;
  cout<<"To merge this reduced/stripped file with another file containing the same entries (but presumably with different variables) do:  python mergetrees.py <output file name> <input file name 1> <input TTreename 1> <input file name 2> <input TTree name 2>"<<endl<<endl ;

  cout<<"e.g."<<endl;
  cout<<"  python mergetrees.py /afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/turbo_2015_data_wAngVars.root /afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/run-II-data/turbo_2015_data.root Lambda_cToKppiTuple/DecayTree /afs/phas.gla.ac.uk/user/n/nwarrack/public_ppe/myLHCb/Gedcode/LHCb_CharmedHadrons/data/AngularTreefile.root newtree"<<endl<<endl;
cout<<"or..."<<endl;
  cout<<"  python mergetrees.py ~/Documents/uni/LHCb_CharmSummerProj/learning_root/strippedToFitData_M_TAU_cut01_wAngVars_wSWeights.root ~/Documents/uni/LHCb_CharmSummerProj/learning_root/strippedToFitData_wAngVars.root DecayTree ~/Documents/uni/LHCb_CharmSummerProj/learning_root/sWeightsTree_M_TAU_cut01.root sTree"<<endl;
  cout<<"========================================================="<<endl<<endl<<endl ;

}
