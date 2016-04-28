#define B2GTTree_cxx
#include "B2GTTree.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <stdlib.h>
#include <algorithm>

void B2GTTree::Loop()
{
//   In a ROOT session, you can do:
//      root> .L B2GTTree.C
//      root> B2GTTree t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   bookNtuple(); 

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   //   std::cout << met_Phi[0] << std::endl;      //
      T_met_Phi = met_Phi[0];
//Leptons:
      smalltree_nlepton = mu_size + el_size;
      for(int i = 0; i < mu_size; i++){
          if(i == 0){
              if(mu_Charge[i] > 0 ){
                smalltree_lept_flav[i] =13;
              } else{
                smalltree_lept_flav[i] = -13;
              }
              smalltree_lept_pt[i] = mu_Pt[i];
              smalltree_lept_eta[i] = mu_Eta[i];
              smalltree_lept_phi[i] = mu_Phi[i];
              smalltree_lept_iso[i] = mu_Iso04[i];
          }
//jets:
    smalltree_njets = jetAK4_size;
    smalltree_jesup_njets = jetAK4_size;
    smalltree_jesdown_njets = jetAK4_size;
    smalltree_jerup_njets = jetAK4_size;
    smalltree_jerdown_njets = jetAK4_size;
    for(int i = 0; i < std::max((int) jetAK4_size, 7); i++){
        smalltree_jet_pt[i] = jetAK4_Pt[i];
        smalltree_jet_eta[i] = jetAK4_Eta[i];   
        smalltree_jet_phi[i] = jetAK4_Phi[i];   
        smalltree_jet_btagdiscri[i] = jetAK4_CSVv2[i];   
        smalltree_jet_btagdiscri_up[i] = jetAK4_CSVv2[i];   
        smalltree_jet_btagdiscri_down[i] = jetAK4_CSVv2[i];   
        smalltree_jet_flav[i] = jetAK4_HadronFlavour[i];   
        smalltree_jet_jesup_pt[i] = jetAK4_Pt[i];
        smalltree_jet_jesup_eta[i] = jetAK4_Eta[i];   
        smalltree_jet_jesup_phi[i] = jetAK4_Phi[i];   
        smalltree_jet_jesup_btagdiscri[i] = jetAK4_CSVv2[i];   
        smalltree_jet_jesup_flav[i] = jetAK4_HadronFlavour[i];   
        smalltree_jet_jesdown_pt[i] = jetAK4_Pt[i];
        smalltree_jet_jesdown_eta[i] = jetAK4_Eta[i];   
        smalltree_jet_jesdown_phi[i] = jetAK4_Phi[i];   
        smalltree_jet_jesdown_btagdiscri[i] = jetAK4_CSVv2[i];   
        smalltree_jet_jesdown_flav[i] = jetAK4_HadronFlavour[i];   
        smalltree_jet_jerup_pt[i] = jetAK4_Pt[i];
        smalltree_jet_jerup_eta[i] = jetAK4_Eta[i];   
        smalltree_jet_jerup_phi[i] = jetAK4_Phi[i];   
        smalltree_jet_jerup_btagdiscri[i] = jetAK4_CSVv2[i];   
        smalltree_jet_jerup_flav[i] = jetAK4_HadronFlavour[i];   
        smalltree_jet_jerdown_pt[i] = jetAK4_Pt[i];
        smalltree_jet_jerdown_eta[i] = jetAK4_Eta[i];   
        smalltree_jet_jerdown_phi[i] = jetAK4_Phi[i];   
        smalltree_jet_jerdown_btagdiscri[i] = jetAK4_CSVv2[i];   
        smalltree_jet_jerdown_flav[i] = jetAK4_HadronFlavour[i];   
    }
//met
    for (int i = 0; i < met_size; i++){
        smalltree_metnosmear_pt = met_uncorPt[i];
        smalltree_metnosmear_phi = met_uncorPhi[i];

        smalltree_met_pt = met_Pt[i];
        smalltree_met_phi = met_Phi[i];

        smalltree_met_jesup_pt = met_Pt[i];
        smalltree_met_jesup_phi = met_Phi[i];
        smalltree_met_jesdown_pt = met_Pt[i];
        smalltree_met_jesdown_phi = met_Phi[i];

        smalltree_met_jerup_pt = met_Pt[i];
        smalltree_met_jerup_phi = met_Phi[i];
        smalltree_met_jerdown_pt = met_Pt[i];
        smalltree_met_jerdown_phi = met_Phi[i];

        smalltree_met_unclsup_pt = met_Pt[i];
        smalltree_met_unclsup_phi = met_Phi[i];
        smalltree_met_unclsdown_pt = met_Pt[i];
        smalltree_met_unclsdown_phi = met_Phi[i];
    }
//random stuff
         smalltree_weight_trigup = 1;
         smalltree_weight_trigdown = 1;
         smalltree_weight_leptup = 1;
         smalltree_weight_leptdown = 1;
         smalltree_weight_PDFup = 1;
         smalltree_weight_PDFdown = 1;
         smalltree_weight_PUup = 1;
         smalltree_weight_PUdown = 1;
         smalltree_weight_toppt = 1;

         smalltree_evtweight = 1;
         smalltree_tmeme = 1;
         smalltree_nvertex = 1;



      }

      //
      //
      if( 1==1/*smalltree_nlepton ==1*/){
        outputTree->Fill();
      }      //
   }
    outputTree->Write();
   //outputTree->Write(); 
   //outputTree->Print();
}
