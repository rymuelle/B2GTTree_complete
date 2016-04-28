//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Apr 19 21:07:55 2016 by ROOT version 6.02/13
// from TTree B2GTree/B2G TTree Ntuple
// found on file: B2GTTreeNtupleExtra_MC_25ns_76X.root
//////////////////////////////////////////////////////////

#ifndef B2GTTree_h
#define B2GTTree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
// Header file for the classes stored in the TTree if any.
#include "vector"
class B2GTTree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain


// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           evt_npv;
   Double_t        evt_vx;
   Double_t        evt_vy;
   Double_t        evt_vz;
   UInt_t          evt_RunNumber;
   UInt_t          evt_LumiBlock;
   ULong64_t       evt_EventNumber;
   Double_t        evt_rho;
   Int_t           pu_NtrueInt;
   UInt_t          pu_size;
   Int_t           pu_BX[16];   //[pu_size]
   Int_t           pu_NInt[16];   //[pu_size]
   UInt_t          vtx_size;
   Int_t           vtx_ndof[37];   //[vtx_size]
   Float_t         vtx_z[37];   //[vtx_size]
   Float_t         vtx_rho[37];   //[vtx_size]
   Float_t         vtx_chi[37];   //[vtx_size]
   UInt_t          met_size;
   Float_t         met_Pt[1];   //[met_size]
   Float_t         met_Px[1];   //[met_size]
   Float_t         met_Py[1];   //[met_size]
   Float_t         met_Phi[1];   //[met_size]
   Float_t         met_uncorPt[1];   //[met_size]
   Float_t         met_uncorPhi[1];   //[met_size]
   Float_t         met_uncorSumEt[1];   //[met_size]
   UInt_t          el_size;
   Float_t         el_Pt[5];   //[el_size]
   Float_t         el_Eta[5];   //[el_size]
   Float_t         el_Phi[5];   //[el_size]
   Float_t         el_E[5];   //[el_size]
   Float_t         el_Charge[5];   //[el_size]
   Float_t         el_Key[5];   //[el_size]
   Float_t         el_Iso03[5];   //[el_size]
   Float_t         el_Iso03db[5];   //[el_size]
   Float_t         el_MiniIso[5];   //[el_size]
   Float_t         el_rho[5];   //[el_size]
   Float_t         el_EA[5];   //[el_size]
   Float_t         el_sumChargedHadronPt[5];   //[el_size]
   Float_t         el_sumNeutralHadronEt[5];   //[el_size]
   Float_t         el_sumPhotonEt[5];   //[el_size]
   Float_t         el_sumPUPt[5];   //[el_size]
   Float_t         el_D0[5];   //[el_size]
   Float_t         el_Dz[5];   //[el_size]
   Float_t         el_dEtaIn[5];   //[el_size]
   Float_t         el_dPhiIn[5];   //[el_size]
   Float_t         el_HoE[5];   //[el_size]
   Float_t         el_full5x5siee[5];   //[el_size]
   Float_t         el_ooEmooP[5];   //[el_size]
   Float_t         el_missHits[5];   //[el_size]
   Float_t         el_hasMatchedConVeto[5];   //[el_size]
   Float_t         el_vidVeto[5];   //[el_size]
   Float_t         el_vidLoose[5];   //[el_size]
   Float_t         el_vidTight[5];   //[el_size]
   Float_t         el_vidMedium[5];   //[el_size]
   Float_t         el_vidHEEP[5];   //[el_size]
   Float_t         el_vidHEEPnoiso[5];   //[el_size]
   Float_t         el_SCEta[5];   //[el_size]
   Float_t         el_SCPhi[5];   //[el_size]
   UInt_t          mu_size;
   Float_t         mu_Pt[13];   //[mu_size]
   Float_t         mu_Eta[13];   //[mu_size]
   Float_t         mu_Phi[13];   //[mu_size]
   Float_t         mu_E[13];   //[mu_size]
   Float_t         mu_Charge[13];   //[mu_size]
   Float_t         mu_Key[13];   //[mu_size]
   Float_t         mu_Iso04[13];   //[mu_size]
   Float_t         mu_MiniIso[13];   //[mu_size]
   Float_t         mu_D0[13];   //[mu_size]
   Float_t         mu_D0err[13];   //[mu_size]
   Float_t         mu_Dxy[13];   //[mu_size]
   Float_t         mu_Dxyerr[13];   //[mu_size]
   Float_t         mu_Dz[13];   //[mu_size]
   Float_t         mu_Dzerr[13];   //[mu_size]
   Float_t         mu_IsSoftMuon[13];   //[mu_size]
   Float_t         mu_IsLooseMuon[13];   //[mu_size]
   Float_t         mu_IsMediumMuon[13];   //[mu_size]
   Float_t         mu_IsTightMuon[13];   //[mu_size]
   Float_t         mu_IsHighPtMuon[13];   //[mu_size]
   Float_t         mu_IsPFMuon[13];   //[mu_size]
   Float_t         mu_IsGlobalMuon[13];   //[mu_size]
   Float_t         mu_IsTrackerMuon[13];   //[mu_size]
   Float_t         mu_GlbTrkNormChi2[13];   //[mu_size]
   Float_t         mu_NumberValidMuonHits[13];   //[mu_size]
   Float_t         mu_NumberMatchedStations[13];   //[mu_size]
   Float_t         mu_NumberValidPixelHits[13];   //[mu_size]
   Float_t         mu_NumberTrackerLayers[13];   //[mu_size]
   Float_t         mu_NumberOfValidTrackerHits[13];   //[mu_size]
   Float_t         mu_NumberOfPixelLayers[13];   //[mu_size]
   Float_t         mu_InTrkNormChi2[13];   //[mu_size]
   Float_t         mu_SumChargedHadronPt[13];   //[mu_size]
   Float_t         mu_SumNeutralHadronPt[13];   //[mu_size]
   Float_t         mu_SumPhotonPt[13];   //[mu_size]
   Float_t         mu_SumPUPt[13];   //[mu_size]
   Float_t         mu_GenMuonY[13];   //[mu_size]
   Float_t         mu_GenMuonEta[13];   //[mu_size]
   Float_t         mu_GenMuonPhi[13];   //[mu_size]
   Float_t         mu_GenMuonPt[13];   //[mu_size]
   Float_t         mu_GenMuonE[13];   //[mu_size]
   Float_t         mu_GenMuonCharge[13];   //[mu_size]
   UInt_t          jetAK4_size;
   Float_t         jetAK4_Pt[13];   //[jetAK4_size]
   Float_t         jetAK4_Eta[13];   //[jetAK4_size]
   Float_t         jetAK4_Phi[13];   //[jetAK4_size]
   Float_t         jetAK4_E[13];   //[jetAK4_size]
   Float_t         jetAK4_Charge[13];   //[jetAK4_size]
   Float_t         jetAK4_CSVv2[13];   //[jetAK4_size]
   Float_t         jetAK4_DoubleBAK8[13];   //[jetAK4_size]
   Float_t         jetAK4_DoubleBCA15[13];   //[jetAK4_size]
   Float_t         jetAK4_CMVAv2[13];   //[jetAK4_size]
   Float_t         jetAK4_CvsL[13];   //[jetAK4_size]
   Float_t         jetAK4_CvsB[13];   //[jetAK4_size]
   Float_t         jetAK4_CMVA[13];   //[jetAK4_size]
   Float_t         jetAK4_GenPartonY[13];   //[jetAK4_size]
   Float_t         jetAK4_GenPartonEta[13];   //[jetAK4_size]
   Float_t         jetAK4_GenPartonPhi[13];   //[jetAK4_size]
   Float_t         jetAK4_GenPartonPt[13];   //[jetAK4_size]
   Float_t         jetAK4_GenPartonE[13];   //[jetAK4_size]
   Float_t         jetAK4_GenPartonCharge[13];   //[jetAK4_size]
   Float_t         jetAK4_PartonFlavour[13];   //[jetAK4_size]
   Float_t         jetAK4_HadronFlavour[13];   //[jetAK4_size]
   Float_t         jetAK4_GenJetY[13];   //[jetAK4_size]
   Float_t         jetAK4_GenJetEta[13];   //[jetAK4_size]
   Float_t         jetAK4_GenJetPhi[13];   //[jetAK4_size]
   Float_t         jetAK4_GenJetPt[13];   //[jetAK4_size]
   Float_t         jetAK4_GenJetE[13];   //[jetAK4_size]
   Float_t         jetAK4_GenJetCharge[13];   //[jetAK4_size]
   Float_t         jetAK4_muonMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_PhotonEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_ElectronEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_MuonEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_HFHadronEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_HFEMEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_ChargedHadronMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_numberOfDaughters[13];   //[jetAK4_size]
   Float_t         jetAK4_neutralHadronMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_neutralHadronEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_neutralEmEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_chargedEmEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_chargedHadronEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_photonMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_electronMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_HFHadronMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_HFEMMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_ChargeMuEnergy[13];   //[jetAK4_size]
   Float_t         jetAK4_neutralMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_neutralHadronEnergyFrac[13];   //[jetAK4_size]
   Float_t         jetAK4_neutralEmEnergyFrac[13];   //[jetAK4_size]
   Float_t         jetAK4_chargedHadronEnergyFrac[13];   //[jetAK4_size]
   Float_t         jetAK4_muonEnergyFrac[13];   //[jetAK4_size]
   Float_t         jetAK4_chargedEmEnergyFrac[13];   //[jetAK4_size]
   Float_t         jetAK4_chargedMultiplicity[13];   //[jetAK4_size]
   Float_t         jetAK4_NumConstituents[13];   //[jetAK4_size]
   Float_t         jetAK4_jecFactor0[13];   //[jetAK4_size]
   Float_t         jetAK4_jecFactorL3Absolute[13];   //[jetAK4_size]
   Float_t         jetAK4_jetArea[13];   //[jetAK4_size]
   Float_t         jetAK4_nSV[13];   //[jetAK4_size]
   Float_t         jetAK4_SV0mass[13];   //[jetAK4_size]
   Float_t         jetAK4_SV1mass[13];   //[jetAK4_size]
   Float_t         jetAK4_QGL[13];   //[jetAK4_size]
   Float_t         jetAK4_jecUncertainty[13];   //[jetAK4_size]
   Float_t         jetAK4_PtResolution[13];   //[jetAK4_size]
   Float_t         jetAK4_JERSF[13];   //[jetAK4_size]
   Float_t         jetAK4_JERSFUp[13];   //[jetAK4_size]
   Float_t         jetAK4_JERSFDown[13];   //[jetAK4_size]
   Float_t         jetAK4_SmearedPt[13];   //[jetAK4_size]
   Float_t         jetAK4_SmearedPEta[13];   //[jetAK4_size]
   Float_t         jetAK4_SmearedPhi[13];   //[jetAK4_size]
   Float_t         jetAK4_SmearedE[13];   //[jetAK4_size]
   vector<vector<int> > *jetAK4_Keys;
   UInt_t          jetAK4Puppi_size;
   Float_t         jetAK4Puppi_Pt[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_Eta[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_Phi[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_E[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_Charge[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_CSVv2[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_DoubleBAK8[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_DoubleBCA15[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_CMVAv2[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_CvsL[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_CvsB[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_CMVA[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenPartonY[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenPartonEta[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenPartonPhi[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenPartonPt[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenPartonE[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenPartonCharge[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_PartonFlavour[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_HadronFlavour[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenJetY[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenJetEta[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenJetPhi[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenJetPt[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenJetE[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_GenJetCharge[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_muonMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_PhotonEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_ElectronEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_MuonEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_HFHadronEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_HFEMEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_ChargedHadronMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_numberOfDaughters[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_neutralHadronMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_neutralHadronEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_neutralEmEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_chargedEmEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_chargedHadronEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_photonMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_electronMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_HFHadronMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_HFEMMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_ChargeMuEnergy[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_neutralMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_neutralHadronEnergyFrac[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_neutralEmEnergyFrac[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_chargedHadronEnergyFrac[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_muonEnergyFrac[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_chargedEmEnergyFrac[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_chargedMultiplicity[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_NumConstituents[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_jecFactor0[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_jecFactorL3Absolute[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_jetArea[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_nSV[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_SV0mass[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_SV1mass[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_jecUncertainty[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_PtResolution[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_JERSF[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_JERSFUp[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_JERSFDown[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_SmearedPt[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_SmearedPEta[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_SmearedPhi[12];   //[jetAK4Puppi_size]
   Float_t         jetAK4Puppi_SmearedE[12];   //[jetAK4Puppi_size]
   vector<vector<int> > *jetAK4Puppi_Keys;
   UInt_t          jetAK8_size;
   Float_t         jetAK8_Pt[5];   //[jetAK8_size]
   Float_t         jetAK8_Eta[5];   //[jetAK8_size]
   Float_t         jetAK8_Phi[5];   //[jetAK8_size]
   Float_t         jetAK8_E[5];   //[jetAK8_size]
   Float_t         jetAK8_Charge[5];   //[jetAK8_size]
   Float_t         jetAK8_CSVv2[5];   //[jetAK8_size]
   Float_t         jetAK8_DoubleBAK8[5];   //[jetAK8_size]
   Float_t         jetAK8_DoubleBCA15[5];   //[jetAK8_size]
   Float_t         jetAK8_CMVAv2[5];   //[jetAK8_size]
   Float_t         jetAK8_CvsL[5];   //[jetAK8_size]
   Float_t         jetAK8_CvsB[5];   //[jetAK8_size]
   Float_t         jetAK8_CMVA[5];   //[jetAK8_size]
   Float_t         jetAK8_GenPartonY[5];   //[jetAK8_size]
   Float_t         jetAK8_GenPartonEta[5];   //[jetAK8_size]
   Float_t         jetAK8_GenPartonPhi[5];   //[jetAK8_size]
   Float_t         jetAK8_GenPartonPt[5];   //[jetAK8_size]
   Float_t         jetAK8_GenPartonE[5];   //[jetAK8_size]
   Float_t         jetAK8_GenPartonCharge[5];   //[jetAK8_size]
   Float_t         jetAK8_PartonFlavour[5];   //[jetAK8_size]
   Float_t         jetAK8_HadronFlavour[5];   //[jetAK8_size]
   Float_t         jetAK8_GenJetY[5];   //[jetAK8_size]
   Float_t         jetAK8_GenJetEta[5];   //[jetAK8_size]
   Float_t         jetAK8_GenJetPhi[5];   //[jetAK8_size]
   Float_t         jetAK8_GenJetPt[5];   //[jetAK8_size]
   Float_t         jetAK8_GenJetE[5];   //[jetAK8_size]
   Float_t         jetAK8_GenJetCharge[5];   //[jetAK8_size]
   Float_t         jetAK8_muonMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_PhotonEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_ElectronEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_MuonEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_HFHadronEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_HFEMEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_ChargedHadronMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_numberOfDaughters[5];   //[jetAK8_size]
   Float_t         jetAK8_neutralHadronMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_neutralHadronEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_neutralEmEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_chargedEmEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_chargedHadronEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_photonMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_electronMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_HFHadronMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_HFEMMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_ChargeMuEnergy[5];   //[jetAK8_size]
   Float_t         jetAK8_neutralMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_neutralHadronEnergyFrac[5];   //[jetAK8_size]
   Float_t         jetAK8_neutralEmEnergyFrac[5];   //[jetAK8_size]
   Float_t         jetAK8_chargedHadronEnergyFrac[5];   //[jetAK8_size]
   Float_t         jetAK8_muonEnergyFrac[5];   //[jetAK8_size]
   Float_t         jetAK8_chargedEmEnergyFrac[5];   //[jetAK8_size]
   Float_t         jetAK8_chargedMultiplicity[5];   //[jetAK8_size]
   Float_t         jetAK8_NumConstituents[5];   //[jetAK8_size]
   Float_t         jetAK8_jecFactor0[5];   //[jetAK8_size]
   Float_t         jetAK8_jecFactorL3Absolute[5];   //[jetAK8_size]
   Float_t         jetAK8_jetArea[5];   //[jetAK8_size]
   Float_t         jetAK8_nSV[5];   //[jetAK8_size]
   Float_t         jetAK8_SV0mass[5];   //[jetAK8_size]
   Float_t         jetAK8_SV1mass[5];   //[jetAK8_size]
   Float_t         jetAK8_jecUncertainty[5];   //[jetAK8_size]
   Float_t         jetAK8_PtResolution[5];   //[jetAK8_size]
   Float_t         jetAK8_JERSF[5];   //[jetAK8_size]
   Float_t         jetAK8_JERSFUp[5];   //[jetAK8_size]
   Float_t         jetAK8_JERSFDown[5];   //[jetAK8_size]
   Float_t         jetAK8_SmearedPt[5];   //[jetAK8_size]
   Float_t         jetAK8_SmearedPEta[5];   //[jetAK8_size]
   Float_t         jetAK8_SmearedPhi[5];   //[jetAK8_size]
   Float_t         jetAK8_SmearedE[5];   //[jetAK8_size]
   Float_t         jetAK8_vSubjetIndex0[5];   //[jetAK8_size]
   Float_t         jetAK8_vSubjetIndex1[5];   //[jetAK8_size]
   Float_t         jetAK8_tau1[5];   //[jetAK8_size]
   Float_t         jetAK8_tau2[5];   //[jetAK8_size]
   Float_t         jetAK8_tau3[5];   //[jetAK8_size]
   Float_t         jetAK8_softDropMass[5];   //[jetAK8_size]
   Float_t         jetAK8_trimmedMass[5];   //[jetAK8_size]
   Float_t         jetAK8_prunedMass[5];   //[jetAK8_size]
   Float_t         jetAK8_filteredMass[5];   //[jetAK8_size]
   vector<vector<int> > *jetAK8_Keys;
   UInt_t          jetAK8Puppi_size;
   Float_t         jetAK8Puppi_Pt[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_Eta[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_Phi[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_E[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_Charge[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_CSVv2[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_DoubleBAK8[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_DoubleBCA15[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_CMVAv2[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_CvsL[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_CvsB[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_CMVA[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenPartonY[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenPartonEta[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenPartonPhi[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenPartonPt[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenPartonE[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenPartonCharge[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_PartonFlavour[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_HadronFlavour[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenJetY[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenJetEta[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenJetPhi[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenJetPt[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenJetE[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_GenJetCharge[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_muonMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_PhotonEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_ElectronEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_MuonEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_HFHadronEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_HFEMEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_ChargedHadronMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_numberOfDaughters[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_neutralHadronMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_neutralHadronEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_neutralEmEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_chargedEmEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_chargedHadronEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_photonMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_electronMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_HFHadronMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_HFEMMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_ChargeMuEnergy[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_neutralMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_neutralHadronEnergyFrac[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_neutralEmEnergyFrac[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_chargedHadronEnergyFrac[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_muonEnergyFrac[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_chargedEmEnergyFrac[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_chargedMultiplicity[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_NumConstituents[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_jecFactor0[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_jecFactorL3Absolute[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_jetArea[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_nSV[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_SV0mass[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_SV1mass[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_jecUncertainty[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_PtResolution[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_JERSF[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_JERSFUp[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_JERSFDown[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_SmearedPt[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_SmearedPEta[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_SmearedPhi[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_SmearedE[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_vSubjetIndex0[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_vSubjetIndex1[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_tau1[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_tau2[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_tau3[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_softDropMass[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_trimmedMass[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_prunedMass[5];   //[jetAK8Puppi_size]
   Float_t         jetAK8Puppi_filteredMass[5];   //[jetAK8Puppi_size]
   vector<vector<int> > *jetAK8Puppi_Keys;
   UInt_t          subjetAK8_size;
   Float_t         subjetAK8_Pt[6];   //[subjetAK8_size]
   Float_t         subjetAK8_Eta[6];   //[subjetAK8_size]
   Float_t         subjetAK8_Phi[6];   //[subjetAK8_size]
   Float_t         subjetAK8_E[6];   //[subjetAK8_size]
   Float_t         subjetAK8_Charge[6];   //[subjetAK8_size]
   Float_t         subjetAK8_CSVv2[6];   //[subjetAK8_size]
   Float_t         subjetAK8_DoubleBAK8[6];   //[subjetAK8_size]
   Float_t         subjetAK8_DoubleBCA15[6];   //[subjetAK8_size]
   Float_t         subjetAK8_CMVAv2[6];   //[subjetAK8_size]
   Float_t         subjetAK8_CvsL[6];   //[subjetAK8_size]
   Float_t         subjetAK8_CvsB[6];   //[subjetAK8_size]
   Float_t         subjetAK8_CMVA[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenPartonY[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenPartonEta[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenPartonPhi[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenPartonPt[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenPartonE[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenPartonCharge[6];   //[subjetAK8_size]
   Float_t         subjetAK8_PartonFlavour[6];   //[subjetAK8_size]
   Float_t         subjetAK8_HadronFlavour[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenJetY[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenJetEta[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenJetPhi[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenJetPt[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenJetE[6];   //[subjetAK8_size]
   Float_t         subjetAK8_GenJetCharge[6];   //[subjetAK8_size]
   Float_t         subjetAK8_muonMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_PhotonEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_ElectronEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_MuonEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_HFHadronEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_HFEMEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_ChargedHadronMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_numberOfDaughters[6];   //[subjetAK8_size]
   Float_t         subjetAK8_neutralHadronMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_neutralHadronEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_neutralEmEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_chargedEmEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_chargedHadronEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_photonMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_electronMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_HFHadronMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_HFEMMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_ChargeMuEnergy[6];   //[subjetAK8_size]
   Float_t         subjetAK8_neutralMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_neutralHadronEnergyFrac[6];   //[subjetAK8_size]
   Float_t         subjetAK8_neutralEmEnergyFrac[6];   //[subjetAK8_size]
   Float_t         subjetAK8_chargedHadronEnergyFrac[6];   //[subjetAK8_size]
   Float_t         subjetAK8_muonEnergyFrac[6];   //[subjetAK8_size]
   Float_t         subjetAK8_chargedEmEnergyFrac[6];   //[subjetAK8_size]
   Float_t         subjetAK8_chargedMultiplicity[6];   //[subjetAK8_size]
   Float_t         subjetAK8_NumConstituents[6];   //[subjetAK8_size]
   Float_t         subjetAK8_jecFactor0[6];   //[subjetAK8_size]
   Float_t         subjetAK8_jecFactorL3Absolute[6];   //[subjetAK8_size]
   Float_t         subjetAK8_jetArea[6];   //[subjetAK8_size]
   Float_t         subjetAK8_nSV[6];   //[subjetAK8_size]
   Float_t         subjetAK8_SV0mass[6];   //[subjetAK8_size]
   Float_t         subjetAK8_SV1mass[6];   //[subjetAK8_size]
   vector<vector<int> > *subjetAK8_Keys;
   UInt_t          subjetAK8Puppi_size;
   Float_t         subjetAK8Puppi_Pt[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_Eta[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_Phi[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_E[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_Charge[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_CSVv2[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_DoubleBAK8[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_DoubleBCA15[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_CMVAv2[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_CvsL[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_CvsB[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_CMVA[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenPartonY[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenPartonEta[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenPartonPhi[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenPartonPt[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenPartonE[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenPartonCharge[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_PartonFlavour[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_HadronFlavour[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenJetY[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenJetEta[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenJetPhi[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenJetPt[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenJetE[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_GenJetCharge[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_muonMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_PhotonEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_ElectronEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_MuonEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_HFHadronEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_HFEMEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_ChargedHadronMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_numberOfDaughters[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_neutralHadronMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_neutralHadronEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_neutralEmEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_chargedEmEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_chargedHadronEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_photonMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_electronMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_HFHadronMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_HFEMMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_ChargeMuEnergy[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_neutralMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_neutralHadronEnergyFrac[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_neutralEmEnergyFrac[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_chargedHadronEnergyFrac[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_muonEnergyFrac[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_chargedEmEnergyFrac[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_chargedMultiplicity[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_NumConstituents[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_jecFactor0[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_jecFactorL3Absolute[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_jetArea[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_nSV[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_SV0mass[6];   //[subjetAK8Puppi_size]
   Float_t         subjetAK8Puppi_SV1mass[6];   //[subjetAK8Puppi_size]
   vector<vector<int> > *subjetAK8Puppi_Keys;
   UInt_t          genjetAK8SD_size;
   Float_t         genjetAK8SD_Pt[39];   //[genjetAK8SD_size]
   Float_t         genjetAK8SD_Eta[39];   //[genjetAK8SD_size]
   Float_t         genjetAK8SD_Phi[39];   //[genjetAK8SD_size]
   Float_t         genjetAK8SD_E[39];   //[genjetAK8SD_size]
   Float_t         genjetAK8SD_Charge[39];   //[genjetAK8SD_size]
   Int_t           Flag_HBHEIsoNoiseFilterResult;
   Int_t           Flag_HBHENoiseFilterResult;
   Int_t           Flag_HBHENoiseFilterResultRun1;
   Int_t           Flag_HBHENoiseFilterResultRun2Loose;
   Int_t           Flag_HBHENoiseFilterResultRun2Tight;
   Int_t           Flag_goodVertices;
   Int_t           Flag_eeBadScFilter;
   Int_t           Flag_ecalLaserCorrFilter;
   Int_t           Flag_EcalDeadCellTriggerPrimitiveFilter;
   Int_t           Flag_EcalDeadCellBoundaryEnergyFilter;
   Int_t           Flag_HcalStripHaloFilter;
   Int_t           Flag_hcalLaserEventFilter;
   Int_t           Flag_HBHENoiseFilter;
   Int_t           Flag_HBHENoiseIsoFilter;
   Int_t           Flag_CSCTightHaloFilter;
   Int_t           Flag_CSCTightHaloTrkMuUnvetoFilter;
   Int_t           Flag_CSCTightHalo2015Filter;
   Int_t           Flag_muonBadTrackFilter;
   Int_t           Flag_chargedHadronTrackResolutionFilter;
   Int_t           Flag_trkPOGFilters;
   Int_t           Flag_trkPOG_manystripclus53X;
   Int_t           Flag_trkPOG_toomanystripclus53X;
   Int_t           Flag_trkPOG_logErrorTooManyClusters;
   Int_t           Flag_METFiltersHLT_AK8PFJet360_TrimMass30;
   Int_t           HLT_PFJet40;
   Int_t           HLT_PFJet40_prescale;
   Int_t           HLT_PFJet60;
   Int_t           HLT_PFJet60_prescale;
   Int_t           HLT_PFJet80;
   Int_t           HLT_PFJet80_prescale;
   Int_t           HLT_PFJet140;
   Int_t           HLT_PFJet140_prescale;
   Int_t           HLT_PFJet200;
   Int_t           HLT_PFJet200_prescale;
   Int_t           HLT_PFJet260;
   Int_t           HLT_PFJet260_prescale;
   Int_t           HLT_PFJet320;
   Int_t           HLT_PFJet320_prescale;
   Int_t           HLT_PFJet400;
   Int_t           HLT_PFJet400_prescale;
   Int_t           HLT_PFJet450;
   Int_t           HLT_PFJet450_prescale;
   Int_t           HLT_PFJet500;
   Int_t           HLT_PFJet500_prescale;
   Int_t           HLT_DiPFJetAve40;
   Int_t           HLT_DiPFJetAve40_prescale;
   Int_t           HLT_DiPFJetAve60;
   Int_t           HLT_DiPFJetAve60_prescale;
   Int_t           HLT_DiPFJetAve80;
   Int_t           HLT_DiPFJetAve80_prescale;
   Int_t           HLT_DiPFJetAve140;
   Int_t           HLT_DiPFJetAve140_prescale;
   Int_t           HLT_DiPFJetAve200;
   Int_t           HLT_DiPFJetAve200_prescale;
   Int_t           HLT_DiPFJetAve260;
   Int_t           HLT_DiPFJetAve260_prescale;
   Int_t           HLT_DiPFJetAve320;
   Int_t           HLT_DiPFJetAve320_prescale;
   Int_t           HLT_DiPFJetAve400;
   Int_t           HLT_DiPFJetAve400_prescale;
   Int_t           HLT_DiPFJetAve500;
   Int_t           HLT_DiPFJetAve500_prescale;
   Int_t           HLT_DiPFJetAve60_HFJEC;
   Int_t           HLT_DiPFJetAve60_HFJEC_prescale;
   Int_t           HLT_DiPFJetAve80_HFJEC;
   Int_t           HLT_DiPFJetAve80_HFJEC_prescale;
   Int_t           HLT_DiPFJetAve100_HFJEC;
   Int_t           HLT_DiPFJetAve100_HFJEC_prescale;
   Int_t           HLT_DiPFJetAve160_HFJEC;
   Int_t           HLT_DiPFJetAve160_HFJEC_prescale;
   Int_t           HLT_DiPFJetAve220_HFJEC;
   Int_t           HLT_DiPFJetAve220_HFJEC_prescale;
   Int_t           HLT_DiPFJetAve300_HFJEC;
   Int_t           HLT_DiPFJetAve300_HFJEC_prescale;
   Int_t           HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45;
   Int_t           HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale;
   Int_t           HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45;
   Int_t           HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale;
   Int_t           HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45;
   Int_t           HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale;
   Int_t           HLT_AK8PFHT650_TrimR0p1PT0p03Mass50;
   Int_t           HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_prescale;
   Int_t           HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;
   Int_t           HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_prescale;
   Int_t           HLT_PFHT550_4JetPt50;
   Int_t           HLT_PFHT550_4JetPt50_prescale;
   Int_t           HLT_PFHT650_4JetPt50;
   Int_t           HLT_PFHT650_4JetPt50_prescale;
   Int_t           HLT_PFHT750_4JetPt50;
   Int_t           HLT_PFHT750_4JetPt50_prescale;
   Int_t           HLT_ECALHT800;
   Int_t           HLT_ECALHT800_prescale;
   Int_t           HLT_PFHT600;
   Int_t           HLT_PFHT600_prescale;
   Int_t           HLT_PFHT650;
   Int_t           HLT_PFHT650_prescale;
   Int_t           HLT_PFHT800;
   Int_t           HLT_PFHT800_prescale;
   Int_t           HLT_PFHT200;
   Int_t           HLT_PFHT200_prescale;
   Int_t           HLT_PFHT250;
   Int_t           HLT_PFHT250_prescale;
   Int_t           HLT_PFHT300;
   Int_t           HLT_PFHT300_prescale;
   Int_t           HLT_PFHT350;
   Int_t           HLT_PFHT350_prescale;
   Int_t           HLT_PFHT400;
   Int_t           HLT_PFHT400_prescale;
   Int_t           HLT_PFHT475;
   Int_t           HLT_PFHT475_prescale;
   Int_t           HLT_Rsq0p25;
   Int_t           HLT_Rsq0p25_prescale;
   Int_t           HLT_Rsq0p30;
   Int_t           HLT_Rsq0p30_prescale;
   Int_t           HLT_RsqMR240_Rsq0p09_MR200;
   Int_t           HLT_RsqMR240_Rsq0p09_MR200_prescale;
   Int_t           HLT_RsqMR240_Rsq0p09_MR200_4jet;
   Int_t           HLT_RsqMR240_Rsq0p09_MR200_4jet_prescale;
   Int_t           HLT_RsqMR270_Rsq0p09_MR200;
   Int_t           HLT_RsqMR270_Rsq0p09_MR200_prescale;
   Int_t           HLT_RsqMR270_Rsq0p09_MR200_4jet;
   Int_t           HLT_RsqMR270_Rsq0p09_MR200_4jet_prescale;
   Int_t           HLT_Mu30_eta2p1_PFJet150_PFJet50;
   Int_t           HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale;
   Int_t           HLT_Mu40_eta2p1_PFJet200_PFJet50;
   Int_t           HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale;
   Int_t           HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;
   Int_t           HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale;
   Int_t           HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;
   Int_t           HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale;
   Int_t           HLT_Ele8_CaloIdM_TrackIdM_PFJet30;
   Int_t           HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale;
   Int_t           HLT_Ele12_CaloIdM_TrackIdM_PFJet30;
   Int_t           HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale;
   Int_t           HLT_Ele18_CaloIdM_TrackIdM_PFJet30;
   Int_t           HLT_Ele18_CaloIdM_TrackIdM_PFJet30_prescale;
   Int_t           HLT_Ele23_CaloIdM_TrackIdM_PFJet30;
   Int_t           HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale;
   Int_t           HLT_Ele33_CaloIdM_TrackIdM_PFJet30;
   Int_t           HLT_Ele33_CaloIdM_TrackIdM_PFJet30_prescale;
   Int_t           HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Int_t           HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;
   Int_t           HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Int_t           HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;
   Int_t           HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30;
   Int_t           HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;
   Int_t           HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07;
   Int_t           HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;
   Int_t           HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07;
   Int_t           HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;
   Int_t           HLT_Ele27_eta2p1_WPLoose_Gsf_HT200;
   Int_t           HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_prescale;
   Int_t           HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF;
   Int_t           HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale;
   Int_t           HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400;
   Int_t           HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;
   Int_t           HLT_Ele15_IsoVVVL_PFHT350_PFMET50;
   Int_t           HLT_Ele15_IsoVVVL_PFHT350_PFMET50_prescale;
   Int_t           HLT_Ele15_IsoVVVL_PFHT600;
   Int_t           HLT_Ele15_IsoVVVL_PFHT600_prescale;
   Int_t           HLT_Ele15_IsoVVVL_PFHT350;
   Int_t           HLT_Ele15_IsoVVVL_PFHT350_prescale;
   Int_t           HLT_Mu3er_PFHT140_PFMET125;
   Int_t           HLT_Mu3er_PFHT140_PFMET125_prescale;
   Int_t           HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72;
   Int_t           HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale;
   Int_t           HLT_Mu6_PFHT200_PFMET100;
   Int_t           HLT_Mu6_PFHT200_PFMET100_prescale;
   Int_t           HLT_Mu10_CentralPFJet30_BTagCSV0p54PF;
   Int_t           HLT_Mu10_CentralPFJet30_BTagCSV0p54PF_prescale;
   Int_t           HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400;
   Int_t           HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;
   Int_t           HLT_Mu15_IsoVVVL_PFHT350_PFMET50;
   Int_t           HLT_Mu15_IsoVVVL_PFHT350_PFMET50_prescale;
   Int_t           HLT_Mu15_IsoVVVL_PFHT600;
   Int_t           HLT_Mu15_IsoVVVL_PFHT600_prescale;
   Int_t           HLT_Mu15_IsoVVVL_PFHT350;
   Int_t           HLT_Mu15_IsoVVVL_PFHT350_prescale;
   Int_t           HLT_Mu8;
   Int_t           HLT_Mu8_prescale;
   Int_t           HLT_Mu17;
   Int_t           HLT_Mu17_prescale;
   Int_t           HLT_Mu20;
   Int_t           HLT_Mu20_prescale;
   Int_t           HLT_TkMu20;
   Int_t           HLT_TkMu20_prescale;
   Int_t           HLT_Mu24_eta2p1;
   Int_t           HLT_Mu24_eta2p1_prescale;
   Int_t           HLT_TkMu24_eta2p1;
   Int_t           HLT_TkMu24_eta2p1_prescale;
   Int_t           HLT_Mu27;
   Int_t           HLT_Mu27_prescale;
   Int_t           HLT_TkMu27;
   Int_t           HLT_TkMu27_prescale;
   Int_t           HLT_Mu45_eta2p1;
   Int_t           HLT_Mu45_eta2p1_prescale;
   Int_t           HLT_Mu50;
   Int_t           HLT_Mu50_prescale;
   Int_t           HLT_Mu55;
   Int_t           HLT_Mu55_prescale;
   Int_t           HLT_Mu300;
   Int_t           HLT_Mu300_prescale;
   Int_t           HLT_Mu350;
   Int_t           HLT_Mu350_prescale;
   Int_t           HLT_Ele105_CaloIdVT_GsfTrkIdT;
   Int_t           HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;
   Int_t           HLT_Ele115_CaloIdVT_GsfTrkIdT;
   Int_t           HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;
   Int_t           HLT_IsoMu17_eta2p1;
   Int_t           HLT_IsoMu17_eta2p1_prescale;
   Int_t           HLT_IsoMu18;
   Int_t           HLT_IsoMu18_prescale;
   Int_t           HLT_OldIsoMu18;
   Int_t           HLT_OldIsoMu18_prescale;
   Int_t           HLT_IsoTkMu18;
   Int_t           HLT_IsoTkMu18_prescale;
   Int_t           HLT_OldIsoTkMu18;
   Int_t           HLT_OldIsoTkMu18_prescale;
   Int_t           HLT_IsoMu20;
   Int_t           HLT_IsoMu20_prescale;
   Int_t           HLT_IsoMu22;
   Int_t           HLT_IsoMu22_prescale;
   Int_t           HLT_IsoMu27;
   Int_t           HLT_IsoMu27_prescale;
   Int_t           HLT_IsoTkMu24_eta2p1;
   Int_t           HLT_IsoTkMu24_eta2p1_prescale;
   Int_t           HLT_IsoTkMu27;
   Int_t           HLT_IsoTkMu27_prescale;
   Int_t           HLT_Mu8_TrkIsoVVL;
   Int_t           HLT_Mu8_TrkIsoVVL_prescale;
   Int_t           HLT_Mu17_TrkIsoVVL;
   Int_t           HLT_Mu17_TrkIsoVVL_prescale;
   Int_t           HLT_Ele22_eta2p1_WPLoose_Gsf;
   Int_t           HLT_Ele22_eta2p1_WPLoose_Gsf_prescale;
   Int_t           HLT_Ele22_eta2p1_WPTight_Gsf;
   Int_t           HLT_Ele22_eta2p1_WPTight_Gsf_prescale;
   Int_t           HLT_Ele23_WPLoose_Gsf;
   Int_t           HLT_Ele23_WPLoose_Gsf_prescale;
   Int_t           HLT_Ele27_WPLoose_Gsf;
   Int_t           HLT_Ele27_WPLoose_Gsf_prescale;
   Int_t           HLT_Ele27_eta2p1_WPLoose_Gsf;
   Int_t           HLT_Ele27_eta2p1_WPLoose_Gsf_prescale;
   Int_t           HLT_Ele27_eta2p1_WPTight_Gsf;
   Int_t           HLT_Ele27_eta2p1_WPTight_Gsf_prescale;
   Int_t           HLT_Ele32_eta2p1_WPTight_Gsf;
   Int_t           HLT_Ele32_eta2p1_WPTight_Gsf_prescale;
   Int_t           HLT_Ele12_CaloIdL_TrackIdL_IsoVL;
   Int_t           HLT_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;
   Int_t           HLT_Ele17_CaloIdL_TrackIdL_IsoVL;
   Int_t           HLT_Ele17_CaloIdL_TrackIdL_IsoVL_prescale;
   Int_t           HLT_Ele23_CaloIdL_TrackIdL_IsoVL;
   Int_t           HLT_Ele23_CaloIdL_TrackIdL_IsoVL_prescale;
   Int_t           evt_NGoodVtx;
   Int_t           evt_NLep;
   Int_t           evt_NTopHad;
   Int_t           evt_NTopHadPreTag;
   Int_t           evt_NTopLep;
   Int_t           evt_NTop;
   Int_t           evt_LHA_PDF_ID;
   Float_t         evt_HtLep;
   Float_t         evt_HtTop;
   Float_t         evt_Ht;
   Float_t         evt_HtAll;
   Float_t         evt_HtEx;
   Float_t         evt_HtExFr;
   Float_t         evt_HtTopFr;
   Float_t         evt_TTHadDR;
   Float_t         evt_TTHadDPhi;
   Float_t         evt_TTHadDEta;
   Float_t         evt_TTHadMass;
   Float_t         evt_TTHadSumPt;
   Float_t         evt_TTHadPz;
   Float_t         evt_TTHadHz;
   Float_t         evt_TTHadDPz;
   Float_t         evt_TTHadMR;
   Float_t         evt_TTHadMTR;
   Float_t         evt_TTHadR;
   Float_t         evt_TTHadR2;
   Float_t         evt_MR;
   Float_t         evt_MTR;
   Float_t         evt_R;
   Float_t         evt_R2;
   Float_t         evt_AK8Puppi_MR;
   Float_t         evt_AK8Puppi_MTR;
   Float_t         evt_AK8Puppi_R;
   Float_t         evt_AK8Puppi_R2;
   Float_t         evt_AK4_MR;
   Float_t         evt_AK4_MTR;
   Float_t         evt_AK4_R;
   Float_t         evt_AK4_R2;
   Float_t         evt_XSec;
   Float_t         evt_NEvent_Corr;
   Float_t         evt_Lumi_Weight;
   Float_t         evt_Gen_Weight;
   Float_t         evt_Gen_Ht;
   Float_t         SUSY_Gluino_Mass;
   Float_t         SUSY_LSP_Mass;
   UInt_t          gen_size;
   Int_t           gen_ID[85];   //[gen_size]
   Int_t           gen_Status[85];   //[gen_size]
   Int_t           gen_Mom0ID[85];   //[gen_size]
   Int_t           gen_Mom0Status[85];   //[gen_size]
   Int_t           gen_Mom1ID[85];   //[gen_size]
   Int_t           gen_Mom1Status[85];   //[gen_size]
   Int_t           gen_Dau0ID[85];   //[gen_size]
   Int_t           gen_Dau0Status[85];   //[gen_size]
   Int_t           gen_Dau1ID[85];   //[gen_size]
   Int_t           gen_Dau1Status[85];   //[gen_size]
   Int_t           jetAK8_HasNearGenTop[5];   //[jetAK8_size]
   Int_t           jetAK8_NearGenTopIsHadronic[5];   //[jetAK8_size]
   Int_t           jetAK8_NearGenWIsHadronic[5];   //[jetAK8_size]
   Int_t           jetAK8_NearGenWToENu[5];   //[jetAK8_size]
   Int_t           jetAK8_NearGenWToMuNu[5];   //[jetAK8_size]
   Int_t           jetAK8_NearGenWToTauNu[5];   //[jetAK8_size]
   Int_t           jetAK8_PassTopTag[5];   //[jetAK8_size]
   Int_t           el_IsPartOfNearAK4Jet[5];   //[el_size]
   Int_t           el_IsPartOfNearAK8Jet[5];   //[el_size]
   Int_t           el_IsPartOfNearSubjet[5];   //[el_size]
   Int_t           mu_IsPartOfNearAK4Jet[13];   //[mu_size]
   Int_t           mu_IsPartOfNearAK8Jet[13];   //[mu_size]
   Int_t           mu_IsPartOfNearSubjet[13];   //[mu_size]
   UInt_t          scale_size;
   Float_t         scale_Weights[6];   //[scale_size]
   UInt_t          pdf_size;
   Float_t         pdf_Weights[100];   //[pdf_size]
   UInt_t          alphas_size;
   Float_t         alphas_Weights[2];   //[alphas_size]
   Float_t         gen_Pt[85];   //[gen_size]
   Float_t         gen_Eta[85];   //[gen_size]
   Float_t         gen_Phi[85];   //[gen_size]
   Float_t         gen_E[85];   //[gen_size]
   Float_t         gen_Mass[85];   //[gen_size]
   Float_t         gen_Charge[85];   //[gen_size]
   Float_t         jetAK8_DRNearGenTop[5];   //[jetAK8_size]
   Float_t         jetAK8_DRNearGenWFromTop[5];   //[jetAK8_size]
   Float_t         jetAK8_DRNearGenBFromTop[5];   //[jetAK8_size]
   Float_t         jetAK8_DRNearGenLepFromSLTop[5];   //[jetAK8_size]
   Float_t         jetAK8_DRNearGenNuFromSLTop[5];   //[jetAK8_size]
   Float_t         jetAK8_PtNearGenTop[5];   //[jetAK8_size]
   Float_t         jetAK8_PtNearGenBFromTop[5];   //[jetAK8_size]
   Float_t         jetAK8_PtNearGenWFromTop[5];   //[jetAK8_size]
   Float_t         jetAK8_PtNearGenLepFromSLTop[5];   //[jetAK8_size]
   Float_t         jetAK8_PtNearGenNuFromSLTop[5];   //[jetAK8_size]
   Float_t         el_DRNearGenEleFromSLTop[5];   //[el_size]
   Float_t         el_PtNearGenEleFromSLTop[5];   //[el_size]
   Float_t         el_PtNearGenTop[5];   //[el_size]
   Float_t         el_LepAK4JetFrac[5];   //[el_size]
   Float_t         el_LepAK8JetFrac[5];   //[el_size]
   Float_t         el_LepSubjetFrac[5];   //[el_size]
   Float_t         el_LepAK4JetMassDrop[5];   //[el_size]
   Float_t         el_LepAK8JetMassDrop[5];   //[el_size]
   Float_t         el_LepSubjetMassDrop[5];   //[el_size]
   Float_t         el_AK4JetV1DR[5];   //[el_size]
   Float_t         el_AK4JetV2DR[5];   //[el_size]
   Float_t         el_AK4JetV3DR[5];   //[el_size]
   Float_t         el_AK8JetV1DR[5];   //[el_size]
   Float_t         el_AK8JetV2DR[5];   //[el_size]
   Float_t         el_AK8JetV3DR[5];   //[el_size]
   Float_t         el_SubjetV1DR[5];   //[el_size]
   Float_t         el_SubjetV2DR[5];   //[el_size]
   Float_t         el_SubjetV3DR[5];   //[el_size]
   Float_t         el_AK4JetV1PtRel[5];   //[el_size]
   Float_t         el_AK4JetV2PtRel[5];   //[el_size]
   Float_t         el_AK4JetV3PtRel[5];   //[el_size]
   Float_t         el_AK8JetV1PtRel[5];   //[el_size]
   Float_t         el_AK8JetV2PtRel[5];   //[el_size]
   Float_t         el_AK8JetV3PtRel[5];   //[el_size]
   Float_t         el_SubjetV1PtRel[5];   //[el_size]
   Float_t         el_SubjetV2PtRel[5];   //[el_size]
   Float_t         el_SubjetV3PtRel[5];   //[el_size]
   Float_t         mu_DRNearGenMuFromSLTop[13];   //[mu_size]
   Float_t         mu_PtNearGenMuFromSLTop[13];   //[mu_size]
   Float_t         mu_PtNearGenTop[13];   //[mu_size]
   Float_t         mu_LepAK4JetFrac[13];   //[mu_size]
   Float_t         mu_LepAK8JetFrac[13];   //[mu_size]
   Float_t         mu_LepSubjetFrac[13];   //[mu_size]
   Float_t         mu_LepAK4JetMassDrop[13];   //[mu_size]
   Float_t         mu_LepAK8JetMassDrop[13];   //[mu_size]
   Float_t         mu_LepSubjetMassDrop[13];   //[mu_size]
   Float_t         mu_AK4JetV1DR[13];   //[mu_size]
   Float_t         mu_AK4JetV2DR[13];   //[mu_size]
   Float_t         mu_AK4JetV3DR[13];   //[mu_size]
   Float_t         mu_AK8JetV1DR[13];   //[mu_size]
   Float_t         mu_AK8JetV2DR[13];   //[mu_size]
   Float_t         mu_AK8JetV3DR[13];   //[mu_size]
   Float_t         mu_SubjetV1DR[13];   //[mu_size]
   Float_t         mu_SubjetV2DR[13];   //[mu_size]
   Float_t         mu_SubjetV3DR[13];   //[mu_size]
   Float_t         mu_AK4JetV1PtRel[13];   //[mu_size]
   Float_t         mu_AK4JetV2PtRel[13];   //[mu_size]
   Float_t         mu_AK4JetV3PtRel[13];   //[mu_size]
   Float_t         mu_AK8JetV1PtRel[13];   //[mu_size]
   Float_t         mu_AK8JetV2PtRel[13];   //[mu_size]
   Float_t         mu_AK8JetV3PtRel[13];   //[mu_size]
   Float_t         mu_SubjetV1PtRel[13];   //[mu_size]
   Float_t         mu_SubjetV2PtRel[13];   //[mu_size]
   Float_t         mu_SubjetV3PtRel[13];   //[mu_size]

   // List of branches
   TBranch        *b_evt_npv;   //!
   TBranch        *b_evt_vx;   //!
   TBranch        *b_evt_vy;   //!
   TBranch        *b_evt_vz;   //!
   TBranch        *b_evt_RunNumber;   //!
   TBranch        *b_evt_LumiBlock;   //!
   TBranch        *b_evt_EventNumber;   //!
   TBranch        *b_evt_rho;   //!
   TBranch        *b_pu_NtrueInt;   //!
   TBranch        *b_pu_size;   //!
   TBranch        *b_pu_BX;   //!
   TBranch        *b_pu_NInt;   //!
   TBranch        *b_vtx_size;   //!
   TBranch        *b_vtx_ndof;   //!
   TBranch        *b_vtx_z;   //!
   TBranch        *b_vtx_rho;   //!
   TBranch        *b_vtx_chi;   //!
   TBranch        *b_met_size;   //!
   TBranch        *b_met_Pt;   //!
   TBranch        *b_met_Px;   //!
   TBranch        *b_met_Py;   //!
   TBranch        *b_met_Phi;   //!
   TBranch        *b_met_uncorPt;   //!
   TBranch        *b_met_uncorPhi;   //!
   TBranch        *b_met_uncorSumEt;   //!
   TBranch        *b_el_size;   //!
   TBranch        *b_el_Pt;   //!
   TBranch        *b_el_Eta;   //!
   TBranch        *b_el_Phi;   //!
   TBranch        *b_el_E;   //!
   TBranch        *b_el_Charge;   //!
   TBranch        *b_el_Key;   //!
   TBranch        *b_el_Iso03;   //!
   TBranch        *b_el_Iso03db;   //!
   TBranch        *b_el_MiniIso;   //!
   TBranch        *b_el_rho;   //!
   TBranch        *b_el_EA;   //!
   TBranch        *b_el_sumChargedHadronPt;   //!
   TBranch        *b_el_sumNeutralHadronEt;   //!
   TBranch        *b_el_sumPhotonEt;   //!
   TBranch        *b_el_sumPUPt;   //!
   TBranch        *b_el_D0;   //!
   TBranch        *b_el_Dz;   //!
   TBranch        *b_el_dEtaIn;   //!
   TBranch        *b_el_dPhiIn;   //!
   TBranch        *b_el_HoE;   //!
   TBranch        *b_el_full5x5siee;   //!
   TBranch        *b_el_ooEmooP;   //!
   TBranch        *b_el_missHits;   //!
   TBranch        *b_el_hasMatchedConVeto;   //!
   TBranch        *b_el_vidVeto;   //!
   TBranch        *b_el_vidLoose;   //!
   TBranch        *b_el_vidTight;   //!
   TBranch        *b_el_vidMedium;   //!
   TBranch        *b_el_vidHEEP;   //!
   TBranch        *b_el_vidHEEPnoiso;   //!
   TBranch        *b_el_SCEta;   //!
   TBranch        *b_el_SCPhi;   //!
   TBranch        *b_mu_size;   //!
   TBranch        *b_mu_Pt;   //!
   TBranch        *b_mu_Eta;   //!
   TBranch        *b_mu_Phi;   //!
   TBranch        *b_mu_E;   //!
   TBranch        *b_mu_Charge;   //!
   TBranch        *b_mu_Key;   //!
   TBranch        *b_mu_Iso04;   //!
   TBranch        *b_mu_MiniIso;   //!
   TBranch        *b_mu_D0;   //!
   TBranch        *b_mu_D0err;   //!
   TBranch        *b_mu_Dxy;   //!
   TBranch        *b_mu_Dxyerr;   //!
   TBranch        *b_mu_Dz;   //!
   TBranch        *b_mu_Dzerr;   //!
   TBranch        *b_mu_IsSoftMuon;   //!
   TBranch        *b_mu_IsLooseMuon;   //!
   TBranch        *b_mu_IsMediumMuon;   //!
   TBranch        *b_mu_IsTightMuon;   //!
   TBranch        *b_mu_IsHighPtMuon;   //!
   TBranch        *b_mu_IsPFMuon;   //!
   TBranch        *b_mu_IsGlobalMuon;   //!
   TBranch        *b_mu_IsTrackerMuon;   //!
   TBranch        *b_mu_GlbTrkNormChi2;   //!
   TBranch        *b_mu_NumberValidMuonHits;   //!
   TBranch        *b_mu_NumberMatchedStations;   //!
   TBranch        *b_mu_NumberValidPixelHits;   //!
   TBranch        *b_mu_NumberTrackerLayers;   //!
   TBranch        *b_mu_NumberOfValidTrackerHits;   //!
   TBranch        *b_mu_NumberOfPixelLayers;   //!
   TBranch        *b_mu_InTrkNormChi2;   //!
   TBranch        *b_mu_SumChargedHadronPt;   //!
   TBranch        *b_mu_SumNeutralHadronPt;   //!
   TBranch        *b_mu_SumPhotonPt;   //!
   TBranch        *b_mu_SumPUPt;   //!
   TBranch        *b_mu_GenMuonY;   //!
   TBranch        *b_mu_GenMuonEta;   //!
   TBranch        *b_mu_GenMuonPhi;   //!
   TBranch        *b_mu_GenMuonPt;   //!
   TBranch        *b_mu_GenMuonE;   //!
   TBranch        *b_mu_GenMuonCharge;   //!
   TBranch        *b_jetAK4_size;   //!
   TBranch        *b_jetAK4_Pt;   //!
   TBranch        *b_jetAK4_Eta;   //!
   TBranch        *b_jetAK4_Phi;   //!
   TBranch        *b_jetAK4_E;   //!
   TBranch        *b_jetAK4_Charge;   //!
   TBranch        *b_jetAK4_CSVv2;   //!
   TBranch        *b_jetAK4_DoubleBAK8;   //!
   TBranch        *b_jetAK4_DoubleBCA15;   //!
   TBranch        *b_jetAK4_CMVAv2;   //!
   TBranch        *b_jetAK4_CvsL;   //!
   TBranch        *b_jetAK4_CvsB;   //!
   TBranch        *b_jetAK4_CMVA;   //!
   TBranch        *b_jetAK4_GenPartonY;   //!
   TBranch        *b_jetAK4_GenPartonEta;   //!
   TBranch        *b_jetAK4_GenPartonPhi;   //!
   TBranch        *b_jetAK4_GenPartonPt;   //!
   TBranch        *b_jetAK4_GenPartonE;   //!
   TBranch        *b_jetAK4_GenPartonCharge;   //!
   TBranch        *b_jetAK4_PartonFlavour;   //!
   TBranch        *b_jetAK4_HadronFlavour;   //!
   TBranch        *b_jetAK4_GenJetY;   //!
   TBranch        *b_jetAK4_GenJetEta;   //!
   TBranch        *b_jetAK4_GenJetPhi;   //!
   TBranch        *b_jetAK4_GenJetPt;   //!
   TBranch        *b_jetAK4_GenJetE;   //!
   TBranch        *b_jetAK4_GenJetCharge;   //!
   TBranch        *b_jetAK4_muonMultiplicity;   //!
   TBranch        *b_jetAK4_PhotonEnergy;   //!
   TBranch        *b_jetAK4_ElectronEnergy;   //!
   TBranch        *b_jetAK4_MuonEnergy;   //!
   TBranch        *b_jetAK4_HFHadronEnergy;   //!
   TBranch        *b_jetAK4_HFEMEnergy;   //!
   TBranch        *b_jetAK4_ChargedHadronMultiplicity;   //!
   TBranch        *b_jetAK4_numberOfDaughters;   //!
   TBranch        *b_jetAK4_neutralHadronMultiplicity;   //!
   TBranch        *b_jetAK4_neutralHadronEnergy;   //!
   TBranch        *b_jetAK4_neutralEmEnergy;   //!
   TBranch        *b_jetAK4_chargedEmEnergy;   //!
   TBranch        *b_jetAK4_chargedHadronEnergy;   //!
   TBranch        *b_jetAK4_photonMultiplicity;   //!
   TBranch        *b_jetAK4_electronMultiplicity;   //!
   TBranch        *b_jetAK4_HFHadronMultiplicity;   //!
   TBranch        *b_jetAK4_HFEMMultiplicity;   //!
   TBranch        *b_jetAK4_ChargeMuEnergy;   //!
   TBranch        *b_jetAK4_neutralMultiplicity;   //!
   TBranch        *b_jetAK4_neutralHadronEnergyFrac;   //!
   TBranch        *b_jetAK4_neutralEmEnergyFrac;   //!
   TBranch        *b_jetAK4_chargedHadronEnergyFrac;   //!
   TBranch        *b_jetAK4_muonEnergyFrac;   //!
   TBranch        *b_jetAK4_chargedEmEnergyFrac;   //!
   TBranch        *b_jetAK4_chargedMultiplicity;   //!
   TBranch        *b_jetAK4_NumConstituents;   //!
   TBranch        *b_jetAK4_jecFactor0;   //!
   TBranch        *b_jetAK4_jecFactorL3Absolute;   //!
   TBranch        *b_jetAK4_jetArea;   //!
   TBranch        *b_jetAK4_nSV;   //!
   TBranch        *b_jetAK4_SV0mass;   //!
   TBranch        *b_jetAK4_SV1mass;   //!
   TBranch        *b_jetAK4_QGL;   //!
   TBranch        *b_jetAK4_jecUncertainty;   //!
   TBranch        *b_jetAK4_PtResolution;   //!
   TBranch        *b_jetAK4_JERSF;   //!
   TBranch        *b_jetAK4_JERSFUp;   //!
   TBranch        *b_jetAK4_JERSFDown;   //!
   TBranch        *b_jetAK4_SmearedPt;   //!
   TBranch        *b_jetAK4_SmearedPEta;   //!
   TBranch        *b_jetAK4_SmearedPhi;   //!
   TBranch        *b_jetAK4_SmearedE;   //!
   TBranch        *b_jetAK4_Keys;   //!
   TBranch        *b_jetAK4Puppi_size;   //!
   TBranch        *b_jetAK4Puppi_Pt;   //!
   TBranch        *b_jetAK4Puppi_Eta;   //!
   TBranch        *b_jetAK4Puppi_Phi;   //!
   TBranch        *b_jetAK4Puppi_E;   //!
   TBranch        *b_jetAK4Puppi_Charge;   //!
   TBranch        *b_jetAK4Puppi_CSVv2;   //!
   TBranch        *b_jetAK4Puppi_DoubleBAK8;   //!
   TBranch        *b_jetAK4Puppi_DoubleBCA15;   //!
   TBranch        *b_jetAK4Puppi_CMVAv2;   //!
   TBranch        *b_jetAK4Puppi_CvsL;   //!
   TBranch        *b_jetAK4Puppi_CvsB;   //!
   TBranch        *b_jetAK4Puppi_CMVA;   //!
   TBranch        *b_jetAK4Puppi_GenPartonY;   //!
   TBranch        *b_jetAK4Puppi_GenPartonEta;   //!
   TBranch        *b_jetAK4Puppi_GenPartonPhi;   //!
   TBranch        *b_jetAK4Puppi_GenPartonPt;   //!
   TBranch        *b_jetAK4Puppi_GenPartonE;   //!
   TBranch        *b_jetAK4Puppi_GenPartonCharge;   //!
   TBranch        *b_jetAK4Puppi_PartonFlavour;   //!
   TBranch        *b_jetAK4Puppi_HadronFlavour;   //!
   TBranch        *b_jetAK4Puppi_GenJetY;   //!
   TBranch        *b_jetAK4Puppi_GenJetEta;   //!
   TBranch        *b_jetAK4Puppi_GenJetPhi;   //!
   TBranch        *b_jetAK4Puppi_GenJetPt;   //!
   TBranch        *b_jetAK4Puppi_GenJetE;   //!
   TBranch        *b_jetAK4Puppi_GenJetCharge;   //!
   TBranch        *b_jetAK4Puppi_muonMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_PhotonEnergy;   //!
   TBranch        *b_jetAK4Puppi_ElectronEnergy;   //!
   TBranch        *b_jetAK4Puppi_MuonEnergy;   //!
   TBranch        *b_jetAK4Puppi_HFHadronEnergy;   //!
   TBranch        *b_jetAK4Puppi_HFEMEnergy;   //!
   TBranch        *b_jetAK4Puppi_ChargedHadronMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_numberOfDaughters;   //!
   TBranch        *b_jetAK4Puppi_neutralHadronMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_neutralHadronEnergy;   //!
   TBranch        *b_jetAK4Puppi_neutralEmEnergy;   //!
   TBranch        *b_jetAK4Puppi_chargedEmEnergy;   //!
   TBranch        *b_jetAK4Puppi_chargedHadronEnergy;   //!
   TBranch        *b_jetAK4Puppi_photonMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_electronMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_HFHadronMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_HFEMMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_ChargeMuEnergy;   //!
   TBranch        *b_jetAK4Puppi_neutralMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_neutralHadronEnergyFrac;   //!
   TBranch        *b_jetAK4Puppi_neutralEmEnergyFrac;   //!
   TBranch        *b_jetAK4Puppi_chargedHadronEnergyFrac;   //!
   TBranch        *b_jetAK4Puppi_muonEnergyFrac;   //!
   TBranch        *b_jetAK4Puppi_chargedEmEnergyFrac;   //!
   TBranch        *b_jetAK4Puppi_chargedMultiplicity;   //!
   TBranch        *b_jetAK4Puppi_NumConstituents;   //!
   TBranch        *b_jetAK4Puppi_jecFactor0;   //!
   TBranch        *b_jetAK4Puppi_jecFactorL3Absolute;   //!
   TBranch        *b_jetAK4Puppi_jetArea;   //!
   TBranch        *b_jetAK4Puppi_nSV;   //!
   TBranch        *b_jetAK4Puppi_SV0mass;   //!
   TBranch        *b_jetAK4Puppi_SV1mass;   //!
   TBranch        *b_jetAK4Puppi_jecUncertainty;   //!
   TBranch        *b_jetAK4Puppi_PtResolution;   //!
   TBranch        *b_jetAK4Puppi_JERSF;   //!
   TBranch        *b_jetAK4Puppi_JERSFUp;   //!
   TBranch        *b_jetAK4Puppi_JERSFDown;   //!
   TBranch        *b_jetAK4Puppi_SmearedPt;   //!
   TBranch        *b_jetAK4Puppi_SmearedPEta;   //!
   TBranch        *b_jetAK4Puppi_SmearedPhi;   //!
   TBranch        *b_jetAK4Puppi_SmearedE;   //!
   TBranch        *b_jetAK4Puppi_Keys;   //!
   TBranch        *b_jetAK8_size;   //!
   TBranch        *b_jetAK8_Pt;   //!
   TBranch        *b_jetAK8_Eta;   //!
   TBranch        *b_jetAK8_Phi;   //!
   TBranch        *b_jetAK8_E;   //!
   TBranch        *b_jetAK8_Charge;   //!
   TBranch        *b_jetAK8_CSVv2;   //!
   TBranch        *b_jetAK8_DoubleBAK8;   //!
   TBranch        *b_jetAK8_DoubleBCA15;   //!
   TBranch        *b_jetAK8_CMVAv2;   //!
   TBranch        *b_jetAK8_CvsL;   //!
   TBranch        *b_jetAK8_CvsB;   //!
   TBranch        *b_jetAK8_CMVA;   //!
   TBranch        *b_jetAK8_GenPartonY;   //!
   TBranch        *b_jetAK8_GenPartonEta;   //!
   TBranch        *b_jetAK8_GenPartonPhi;   //!
   TBranch        *b_jetAK8_GenPartonPt;   //!
   TBranch        *b_jetAK8_GenPartonE;   //!
   TBranch        *b_jetAK8_GenPartonCharge;   //!
   TBranch        *b_jetAK8_PartonFlavour;   //!
   TBranch        *b_jetAK8_HadronFlavour;   //!
   TBranch        *b_jetAK8_GenJetY;   //!
   TBranch        *b_jetAK8_GenJetEta;   //!
   TBranch        *b_jetAK8_GenJetPhi;   //!
   TBranch        *b_jetAK8_GenJetPt;   //!
   TBranch        *b_jetAK8_GenJetE;   //!
   TBranch        *b_jetAK8_GenJetCharge;   //!
   TBranch        *b_jetAK8_muonMultiplicity;   //!
   TBranch        *b_jetAK8_PhotonEnergy;   //!
   TBranch        *b_jetAK8_ElectronEnergy;   //!
   TBranch        *b_jetAK8_MuonEnergy;   //!
   TBranch        *b_jetAK8_HFHadronEnergy;   //!
   TBranch        *b_jetAK8_HFEMEnergy;   //!
   TBranch        *b_jetAK8_ChargedHadronMultiplicity;   //!
   TBranch        *b_jetAK8_numberOfDaughters;   //!
   TBranch        *b_jetAK8_neutralHadronMultiplicity;   //!
   TBranch        *b_jetAK8_neutralHadronEnergy;   //!
   TBranch        *b_jetAK8_neutralEmEnergy;   //!
   TBranch        *b_jetAK8_chargedEmEnergy;   //!
   TBranch        *b_jetAK8_chargedHadronEnergy;   //!
   TBranch        *b_jetAK8_photonMultiplicity;   //!
   TBranch        *b_jetAK8_electronMultiplicity;   //!
   TBranch        *b_jetAK8_HFHadronMultiplicity;   //!
   TBranch        *b_jetAK8_HFEMMultiplicity;   //!
   TBranch        *b_jetAK8_ChargeMuEnergy;   //!
   TBranch        *b_jetAK8_neutralMultiplicity;   //!
   TBranch        *b_jetAK8_neutralHadronEnergyFrac;   //!
   TBranch        *b_jetAK8_neutralEmEnergyFrac;   //!
   TBranch        *b_jetAK8_chargedHadronEnergyFrac;   //!
   TBranch        *b_jetAK8_muonEnergyFrac;   //!
   TBranch        *b_jetAK8_chargedEmEnergyFrac;   //!
   TBranch        *b_jetAK8_chargedMultiplicity;   //!
   TBranch        *b_jetAK8_NumConstituents;   //!
   TBranch        *b_jetAK8_jecFactor0;   //!
   TBranch        *b_jetAK8_jecFactorL3Absolute;   //!
   TBranch        *b_jetAK8_jetArea;   //!
   TBranch        *b_jetAK8_nSV;   //!
   TBranch        *b_jetAK8_SV0mass;   //!
   TBranch        *b_jetAK8_SV1mass;   //!
   TBranch        *b_jetAK8_jecUncertainty;   //!
   TBranch        *b_jetAK8_PtResolution;   //!
   TBranch        *b_jetAK8_JERSF;   //!
   TBranch        *b_jetAK8_JERSFUp;   //!
   TBranch        *b_jetAK8_JERSFDown;   //!
   TBranch        *b_jetAK8_SmearedPt;   //!
   TBranch        *b_jetAK8_SmearedPEta;   //!
   TBranch        *b_jetAK8_SmearedPhi;   //!
   TBranch        *b_jetAK8_SmearedE;   //!
   TBranch        *b_jetAK8_vSubjetIndex0;   //!
   TBranch        *b_jetAK8_vSubjetIndex1;   //!
   TBranch        *b_jetAK8_tau1;   //!
   TBranch        *b_jetAK8_tau2;   //!
   TBranch        *b_jetAK8_tau3;   //!
   TBranch        *b_jetAK8_softDropMass;   //!
   TBranch        *b_jetAK8_trimmedMass;   //!
   TBranch        *b_jetAK8_prunedMass;   //!
   TBranch        *b_jetAK8_filteredMass;   //!
   TBranch        *b_jetAK8_Keys;   //!
   TBranch        *b_jetAK8Puppi_size;   //!
   TBranch        *b_jetAK8Puppi_Pt;   //!
   TBranch        *b_jetAK8Puppi_Eta;   //!
   TBranch        *b_jetAK8Puppi_Phi;   //!
   TBranch        *b_jetAK8Puppi_E;   //!
   TBranch        *b_jetAK8Puppi_Charge;   //!
   TBranch        *b_jetAK8Puppi_CSVv2;   //!
   TBranch        *b_jetAK8Puppi_DoubleBAK8;   //!
   TBranch        *b_jetAK8Puppi_DoubleBCA15;   //!
   TBranch        *b_jetAK8Puppi_CMVAv2;   //!
   TBranch        *b_jetAK8Puppi_CvsL;   //!
   TBranch        *b_jetAK8Puppi_CvsB;   //!
   TBranch        *b_jetAK8Puppi_CMVA;   //!
   TBranch        *b_jetAK8Puppi_GenPartonY;   //!
   TBranch        *b_jetAK8Puppi_GenPartonEta;   //!
   TBranch        *b_jetAK8Puppi_GenPartonPhi;   //!
   TBranch        *b_jetAK8Puppi_GenPartonPt;   //!
   TBranch        *b_jetAK8Puppi_GenPartonE;   //!
   TBranch        *b_jetAK8Puppi_GenPartonCharge;   //!
   TBranch        *b_jetAK8Puppi_PartonFlavour;   //!
   TBranch        *b_jetAK8Puppi_HadronFlavour;   //!
   TBranch        *b_jetAK8Puppi_GenJetY;   //!
   TBranch        *b_jetAK8Puppi_GenJetEta;   //!
   TBranch        *b_jetAK8Puppi_GenJetPhi;   //!
   TBranch        *b_jetAK8Puppi_GenJetPt;   //!
   TBranch        *b_jetAK8Puppi_GenJetE;   //!
   TBranch        *b_jetAK8Puppi_GenJetCharge;   //!
   TBranch        *b_jetAK8Puppi_muonMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_PhotonEnergy;   //!
   TBranch        *b_jetAK8Puppi_ElectronEnergy;   //!
   TBranch        *b_jetAK8Puppi_MuonEnergy;   //!
   TBranch        *b_jetAK8Puppi_HFHadronEnergy;   //!
   TBranch        *b_jetAK8Puppi_HFEMEnergy;   //!
   TBranch        *b_jetAK8Puppi_ChargedHadronMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_numberOfDaughters;   //!
   TBranch        *b_jetAK8Puppi_neutralHadronMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_neutralHadronEnergy;   //!
   TBranch        *b_jetAK8Puppi_neutralEmEnergy;   //!
   TBranch        *b_jetAK8Puppi_chargedEmEnergy;   //!
   TBranch        *b_jetAK8Puppi_chargedHadronEnergy;   //!
   TBranch        *b_jetAK8Puppi_photonMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_electronMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_HFHadronMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_HFEMMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_ChargeMuEnergy;   //!
   TBranch        *b_jetAK8Puppi_neutralMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_neutralHadronEnergyFrac;   //!
   TBranch        *b_jetAK8Puppi_neutralEmEnergyFrac;   //!
   TBranch        *b_jetAK8Puppi_chargedHadronEnergyFrac;   //!
   TBranch        *b_jetAK8Puppi_muonEnergyFrac;   //!
   TBranch        *b_jetAK8Puppi_chargedEmEnergyFrac;   //!
   TBranch        *b_jetAK8Puppi_chargedMultiplicity;   //!
   TBranch        *b_jetAK8Puppi_NumConstituents;   //!
   TBranch        *b_jetAK8Puppi_jecFactor0;   //!
   TBranch        *b_jetAK8Puppi_jecFactorL3Absolute;   //!
   TBranch        *b_jetAK8Puppi_jetArea;   //!
   TBranch        *b_jetAK8Puppi_nSV;   //!
   TBranch        *b_jetAK8Puppi_SV0mass;   //!
   TBranch        *b_jetAK8Puppi_SV1mass;   //!
   TBranch        *b_jetAK8Puppi_jecUncertainty;   //!
   TBranch        *b_jetAK8Puppi_PtResolution;   //!
   TBranch        *b_jetAK8Puppi_JERSF;   //!
   TBranch        *b_jetAK8Puppi_JERSFUp;   //!
   TBranch        *b_jetAK8Puppi_JERSFDown;   //!
   TBranch        *b_jetAK8Puppi_SmearedPt;   //!
   TBranch        *b_jetAK8Puppi_SmearedPEta;   //!
   TBranch        *b_jetAK8Puppi_SmearedPhi;   //!
   TBranch        *b_jetAK8Puppi_SmearedE;   //!
   TBranch        *b_jetAK8Puppi_vSubjetIndex0;   //!
   TBranch        *b_jetAK8Puppi_vSubjetIndex1;   //!
   TBranch        *b_jetAK8Puppi_tau1;   //!
   TBranch        *b_jetAK8Puppi_tau2;   //!
   TBranch        *b_jetAK8Puppi_tau3;   //!
   TBranch        *b_jetAK8Puppi_softDropMass;   //!
   TBranch        *b_jetAK8Puppi_trimmedMass;   //!
   TBranch        *b_jetAK8Puppi_prunedMass;   //!
   TBranch        *b_jetAK8Puppi_filteredMass;   //!
   TBranch        *b_jetAK8Puppi_Keys;   //!
   TBranch        *b_subjetAK8_size;   //!
   TBranch        *b_subjetAK8_Pt;   //!
   TBranch        *b_subjetAK8_Eta;   //!
   TBranch        *b_subjetAK8_Phi;   //!
   TBranch        *b_subjetAK8_E;   //!
   TBranch        *b_subjetAK8_Charge;   //!
   TBranch        *b_subjetAK8_CSVv2;   //!
   TBranch        *b_subjetAK8_DoubleBAK8;   //!
   TBranch        *b_subjetAK8_DoubleBCA15;   //!
   TBranch        *b_subjetAK8_CMVAv2;   //!
   TBranch        *b_subjetAK8_CvsL;   //!
   TBranch        *b_subjetAK8_CvsB;   //!
   TBranch        *b_subjetAK8_CMVA;   //!
   TBranch        *b_subjetAK8_GenPartonY;   //!
   TBranch        *b_subjetAK8_GenPartonEta;   //!
   TBranch        *b_subjetAK8_GenPartonPhi;   //!
   TBranch        *b_subjetAK8_GenPartonPt;   //!
   TBranch        *b_subjetAK8_GenPartonE;   //!
   TBranch        *b_subjetAK8_GenPartonCharge;   //!
   TBranch        *b_subjetAK8_PartonFlavour;   //!
   TBranch        *b_subjetAK8_HadronFlavour;   //!
   TBranch        *b_subjetAK8_GenJetY;   //!
   TBranch        *b_subjetAK8_GenJetEta;   //!
   TBranch        *b_subjetAK8_GenJetPhi;   //!
   TBranch        *b_subjetAK8_GenJetPt;   //!
   TBranch        *b_subjetAK8_GenJetE;   //!
   TBranch        *b_subjetAK8_GenJetCharge;   //!
   TBranch        *b_subjetAK8_muonMultiplicity;   //!
   TBranch        *b_subjetAK8_PhotonEnergy;   //!
   TBranch        *b_subjetAK8_ElectronEnergy;   //!
   TBranch        *b_subjetAK8_MuonEnergy;   //!
   TBranch        *b_subjetAK8_HFHadronEnergy;   //!
   TBranch        *b_subjetAK8_HFEMEnergy;   //!
   TBranch        *b_subjetAK8_ChargedHadronMultiplicity;   //!
   TBranch        *b_subjetAK8_numberOfDaughters;   //!
   TBranch        *b_subjetAK8_neutralHadronMultiplicity;   //!
   TBranch        *b_subjetAK8_neutralHadronEnergy;   //!
   TBranch        *b_subjetAK8_neutralEmEnergy;   //!
   TBranch        *b_subjetAK8_chargedEmEnergy;   //!
   TBranch        *b_subjetAK8_chargedHadronEnergy;   //!
   TBranch        *b_subjetAK8_photonMultiplicity;   //!
   TBranch        *b_subjetAK8_electronMultiplicity;   //!
   TBranch        *b_subjetAK8_HFHadronMultiplicity;   //!
   TBranch        *b_subjetAK8_HFEMMultiplicity;   //!
   TBranch        *b_subjetAK8_ChargeMuEnergy;   //!
   TBranch        *b_subjetAK8_neutralMultiplicity;   //!
   TBranch        *b_subjetAK8_neutralHadronEnergyFrac;   //!
   TBranch        *b_subjetAK8_neutralEmEnergyFrac;   //!
   TBranch        *b_subjetAK8_chargedHadronEnergyFrac;   //!
   TBranch        *b_subjetAK8_muonEnergyFrac;   //!
   TBranch        *b_subjetAK8_chargedEmEnergyFrac;   //!
   TBranch        *b_subjetAK8_chargedMultiplicity;   //!
   TBranch        *b_subjetAK8_NumConstituents;   //!
   TBranch        *b_subjetAK8_jecFactor0;   //!
   TBranch        *b_subjetAK8_jecFactorL3Absolute;   //!
   TBranch        *b_subjetAK8_jetArea;   //!
   TBranch        *b_subjetAK8_nSV;   //!
   TBranch        *b_subjetAK8_SV0mass;   //!
   TBranch        *b_subjetAK8_SV1mass;   //!
   TBranch        *b_subjetAK8_Keys;   //!
   TBranch        *b_subjetAK8Puppi_size;   //!
   TBranch        *b_subjetAK8Puppi_Pt;   //!
   TBranch        *b_subjetAK8Puppi_Eta;   //!
   TBranch        *b_subjetAK8Puppi_Phi;   //!
   TBranch        *b_subjetAK8Puppi_E;   //!
   TBranch        *b_subjetAK8Puppi_Charge;   //!
   TBranch        *b_subjetAK8Puppi_CSVv2;   //!
   TBranch        *b_subjetAK8Puppi_DoubleBAK8;   //!
   TBranch        *b_subjetAK8Puppi_DoubleBCA15;   //!
   TBranch        *b_subjetAK8Puppi_CMVAv2;   //!
   TBranch        *b_subjetAK8Puppi_CvsL;   //!
   TBranch        *b_subjetAK8Puppi_CvsB;   //!
   TBranch        *b_subjetAK8Puppi_CMVA;   //!
   TBranch        *b_subjetAK8Puppi_GenPartonY;   //!
   TBranch        *b_subjetAK8Puppi_GenPartonEta;   //!
   TBranch        *b_subjetAK8Puppi_GenPartonPhi;   //!
   TBranch        *b_subjetAK8Puppi_GenPartonPt;   //!
   TBranch        *b_subjetAK8Puppi_GenPartonE;   //!
   TBranch        *b_subjetAK8Puppi_GenPartonCharge;   //!
   TBranch        *b_subjetAK8Puppi_PartonFlavour;   //!
   TBranch        *b_subjetAK8Puppi_HadronFlavour;   //!
   TBranch        *b_subjetAK8Puppi_GenJetY;   //!
   TBranch        *b_subjetAK8Puppi_GenJetEta;   //!
   TBranch        *b_subjetAK8Puppi_GenJetPhi;   //!
   TBranch        *b_subjetAK8Puppi_GenJetPt;   //!
   TBranch        *b_subjetAK8Puppi_GenJetE;   //!
   TBranch        *b_subjetAK8Puppi_GenJetCharge;   //!
   TBranch        *b_subjetAK8Puppi_muonMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_PhotonEnergy;   //!
   TBranch        *b_subjetAK8Puppi_ElectronEnergy;   //!
   TBranch        *b_subjetAK8Puppi_MuonEnergy;   //!
   TBranch        *b_subjetAK8Puppi_HFHadronEnergy;   //!
   TBranch        *b_subjetAK8Puppi_HFEMEnergy;   //!
   TBranch        *b_subjetAK8Puppi_ChargedHadronMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_numberOfDaughters;   //!
   TBranch        *b_subjetAK8Puppi_neutralHadronMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_neutralHadronEnergy;   //!
   TBranch        *b_subjetAK8Puppi_neutralEmEnergy;   //!
   TBranch        *b_subjetAK8Puppi_chargedEmEnergy;   //!
   TBranch        *b_subjetAK8Puppi_chargedHadronEnergy;   //!
   TBranch        *b_subjetAK8Puppi_photonMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_electronMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_HFHadronMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_HFEMMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_ChargeMuEnergy;   //!
   TBranch        *b_subjetAK8Puppi_neutralMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_neutralHadronEnergyFrac;   //!
   TBranch        *b_subjetAK8Puppi_neutralEmEnergyFrac;   //!
   TBranch        *b_subjetAK8Puppi_chargedHadronEnergyFrac;   //!
   TBranch        *b_subjetAK8Puppi_muonEnergyFrac;   //!
   TBranch        *b_subjetAK8Puppi_chargedEmEnergyFrac;   //!
   TBranch        *b_subjetAK8Puppi_chargedMultiplicity;   //!
   TBranch        *b_subjetAK8Puppi_NumConstituents;   //!
   TBranch        *b_subjetAK8Puppi_jecFactor0;   //!
   TBranch        *b_subjetAK8Puppi_jecFactorL3Absolute;   //!
   TBranch        *b_subjetAK8Puppi_jetArea;   //!
   TBranch        *b_subjetAK8Puppi_nSV;   //!
   TBranch        *b_subjetAK8Puppi_SV0mass;   //!
   TBranch        *b_subjetAK8Puppi_SV1mass;   //!
   TBranch        *b_subjetAK8Puppi_Keys;   //!
   TBranch        *b_genjetAK8SD_size;   //!
   TBranch        *b_genjetAK8SD_Pt;   //!
   TBranch        *b_genjetAK8SD_Eta;   //!
   TBranch        *b_genjetAK8SD_Phi;   //!
   TBranch        *b_genjetAK8SD_E;   //!
   TBranch        *b_genjetAK8SD_Charge;   //!
   TBranch        *b_Flag_HBHEIsoNoiseFilterResult;   //!
   TBranch        *b_Flag_HBHENoiseFilterResult;   //!
   TBranch        *b_Flag_HBHENoiseFilterResultRun1;   //!
   TBranch        *b_Flag_HBHENoiseFilterResultRun2Loose;   //!
   TBranch        *b_Flag_HBHENoiseFilterResultRun2Tight;   //!
   TBranch        *b_Flag_goodVertices;   //!
   TBranch        *b_Flag_eeBadScFilter;   //!
   TBranch        *b_Flag_ecalLaserCorrFilter;   //!
   TBranch        *b_Flag_EcalDeadCellTriggerPrimitiveFilter;   //!
   TBranch        *b_Flag_EcalDeadCellBoundaryEnergyFilter;   //!
   TBranch        *b_Flag_HcalStripHaloFilter;   //!
   TBranch        *b_Flag_hcalLaserEventFilter;   //!
   TBranch        *b_Flag_HBHENoiseFilter;   //!
   TBranch        *b_Flag_HBHENoiseIsoFilter;   //!
   TBranch        *b_Flag_CSCTightHaloFilter;   //!
   TBranch        *b_Flag_CSCTightHaloTrkMuUnvetoFilter;   //!
   TBranch        *b_Flag_CSCTightHalo2015Filter;   //!
   TBranch        *b_Flag_muonBadTrackFilter;   //!
   TBranch        *b_Flag_chargedHadronTrackResolutionFilter;   //!
   TBranch        *b_Flag_trkPOGFilters;   //!
   TBranch        *b_Flag_trkPOG_manystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_toomanystripclus53X;   //!
   TBranch        *b_Flag_trkPOG_logErrorTooManyClusters;   //!
   TBranch        *b_Flag_METFiltersHLT_AK8PFJet360_TrimMass30;   //!
   TBranch        *b_HLT_PFJet40;   //!
   TBranch        *b_HLT_PFJet40_prescale;   //!
   TBranch        *b_HLT_PFJet60;   //!
   TBranch        *b_HLT_PFJet60_prescale;   //!
   TBranch        *b_HLT_PFJet80;   //!
   TBranch        *b_HLT_PFJet80_prescale;   //!
   TBranch        *b_HLT_PFJet140;   //!
   TBranch        *b_HLT_PFJet140_prescale;   //!
   TBranch        *b_HLT_PFJet200;   //!
   TBranch        *b_HLT_PFJet200_prescale;   //!
   TBranch        *b_HLT_PFJet260;   //!
   TBranch        *b_HLT_PFJet260_prescale;   //!
   TBranch        *b_HLT_PFJet320;   //!
   TBranch        *b_HLT_PFJet320_prescale;   //!
   TBranch        *b_HLT_PFJet400;   //!
   TBranch        *b_HLT_PFJet400_prescale;   //!
   TBranch        *b_HLT_PFJet450;   //!
   TBranch        *b_HLT_PFJet450_prescale;   //!
   TBranch        *b_HLT_PFJet500;   //!
   TBranch        *b_HLT_PFJet500_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve40;   //!
   TBranch        *b_HLT_DiPFJetAve40_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve60;   //!
   TBranch        *b_HLT_DiPFJetAve60_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve80;   //!
   TBranch        *b_HLT_DiPFJetAve80_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve140;   //!
   TBranch        *b_HLT_DiPFJetAve140_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve200;   //!
   TBranch        *b_HLT_DiPFJetAve200_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve260;   //!
   TBranch        *b_HLT_DiPFJetAve260_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve320;   //!
   TBranch        *b_HLT_DiPFJetAve320_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve400;   //!
   TBranch        *b_HLT_DiPFJetAve400_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve500;   //!
   TBranch        *b_HLT_DiPFJetAve500_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve60_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve60_HFJEC_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve80_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve80_HFJEC_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve100_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve100_HFJEC_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve160_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve160_HFJEC_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve220_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve220_HFJEC_prescale;   //!
   TBranch        *b_HLT_DiPFJetAve300_HFJEC;   //!
   TBranch        *b_HLT_DiPFJetAve300_HFJEC_prescale;   //!
   TBranch        *b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45;   //!
   TBranch        *b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45;   //!
   TBranch        *b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale;   //!
   TBranch        *b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45;   //!
   TBranch        *b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale;   //!
   TBranch        *b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50;   //!
   TBranch        *b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_prescale;   //!
   TBranch        *b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50;   //!
   TBranch        *b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_prescale;   //!
   TBranch        *b_HLT_PFHT550_4JetPt50;   //!
   TBranch        *b_HLT_PFHT550_4JetPt50_prescale;   //!
   TBranch        *b_HLT_PFHT650_4JetPt50;   //!
   TBranch        *b_HLT_PFHT650_4JetPt50_prescale;   //!
   TBranch        *b_HLT_PFHT750_4JetPt50;   //!
   TBranch        *b_HLT_PFHT750_4JetPt50_prescale;   //!
   TBranch        *b_HLT_ECALHT800;   //!
   TBranch        *b_HLT_ECALHT800_prescale;   //!
   TBranch        *b_HLT_PFHT600;   //!
   TBranch        *b_HLT_PFHT600_prescale;   //!
   TBranch        *b_HLT_PFHT650;   //!
   TBranch        *b_HLT_PFHT650_prescale;   //!
   TBranch        *b_HLT_PFHT800;   //!
   TBranch        *b_HLT_PFHT800_prescale;   //!
   TBranch        *b_HLT_PFHT200;   //!
   TBranch        *b_HLT_PFHT200_prescale;   //!
   TBranch        *b_HLT_PFHT250;   //!
   TBranch        *b_HLT_PFHT250_prescale;   //!
   TBranch        *b_HLT_PFHT300;   //!
   TBranch        *b_HLT_PFHT300_prescale;   //!
   TBranch        *b_HLT_PFHT350;   //!
   TBranch        *b_HLT_PFHT350_prescale;   //!
   TBranch        *b_HLT_PFHT400;   //!
   TBranch        *b_HLT_PFHT400_prescale;   //!
   TBranch        *b_HLT_PFHT475;   //!
   TBranch        *b_HLT_PFHT475_prescale;   //!
   TBranch        *b_HLT_Rsq0p25;   //!
   TBranch        *b_HLT_Rsq0p25_prescale;   //!
   TBranch        *b_HLT_Rsq0p30;   //!
   TBranch        *b_HLT_Rsq0p30_prescale;   //!
   TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200;   //!
   TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200_prescale;   //!
   TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_HLT_RsqMR240_Rsq0p09_MR200_4jet_prescale;   //!
   TBranch        *b_HLT_RsqMR270_Rsq0p09_MR200;   //!
   TBranch        *b_HLT_RsqMR270_Rsq0p09_MR200_prescale;   //!
   TBranch        *b_HLT_RsqMR270_Rsq0p09_MR200_4jet;   //!
   TBranch        *b_HLT_RsqMR270_Rsq0p09_MR200_4jet_prescale;   //!
   TBranch        *b_HLT_Mu30_eta2p1_PFJet150_PFJet50;   //!
   TBranch        *b_HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale;   //!
   TBranch        *b_HLT_Mu40_eta2p1_PFJet200_PFJet50;   //!
   TBranch        *b_HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale;   //!
   TBranch        *b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50;   //!
   TBranch        *b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale;   //!
   TBranch        *b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50;   //!
   TBranch        *b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale;   //!
   TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele18_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele18_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele33_CaloIdM_TrackIdM_PFJet30;   //!
   TBranch        *b_HLT_Ele33_CaloIdM_TrackIdM_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30;   //!
   TBranch        *b_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale;   //!
   TBranch        *b_HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07;   //!
   TBranch        *b_HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;   //!
   TBranch        *b_HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07;   //!
   TBranch        *b_HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale;   //!
   TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200;   //!
   TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_prescale;   //!
   TBranch        *b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF;   //!
   TBranch        *b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50_prescale;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT600_prescale;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT350;   //!
   TBranch        *b_HLT_Ele15_IsoVVVL_PFHT350_prescale;   //!
   TBranch        *b_HLT_Mu3er_PFHT140_PFMET125;   //!
   TBranch        *b_HLT_Mu3er_PFHT140_PFMET125_prescale;   //!
   TBranch        *b_HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72;   //!
   TBranch        *b_HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale;   //!
   TBranch        *b_HLT_Mu6_PFHT200_PFMET100;   //!
   TBranch        *b_HLT_Mu6_PFHT200_PFMET100_prescale;   //!
   TBranch        *b_HLT_Mu10_CentralPFJet30_BTagCSV0p54PF;   //!
   TBranch        *b_HLT_Mu10_CentralPFJet30_BTagCSV0p54PF_prescale;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50_prescale;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT600;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT600_prescale;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT350;   //!
   TBranch        *b_HLT_Mu15_IsoVVVL_PFHT350_prescale;   //!
   TBranch        *b_HLT_Mu8;   //!
   TBranch        *b_HLT_Mu8_prescale;   //!
   TBranch        *b_HLT_Mu17;   //!
   TBranch        *b_HLT_Mu17_prescale;   //!
   TBranch        *b_HLT_Mu20;   //!
   TBranch        *b_HLT_Mu20_prescale;   //!
   TBranch        *b_HLT_TkMu20;   //!
   TBranch        *b_HLT_TkMu20_prescale;   //!
   TBranch        *b_HLT_Mu24_eta2p1;   //!
   TBranch        *b_HLT_Mu24_eta2p1_prescale;   //!
   TBranch        *b_HLT_TkMu24_eta2p1;   //!
   TBranch        *b_HLT_TkMu24_eta2p1_prescale;   //!
   TBranch        *b_HLT_Mu27;   //!
   TBranch        *b_HLT_Mu27_prescale;   //!
   TBranch        *b_HLT_TkMu27;   //!
   TBranch        *b_HLT_TkMu27_prescale;   //!
   TBranch        *b_HLT_Mu45_eta2p1;   //!
   TBranch        *b_HLT_Mu45_eta2p1_prescale;   //!
   TBranch        *b_HLT_Mu50;   //!
   TBranch        *b_HLT_Mu50_prescale;   //!
   TBranch        *b_HLT_Mu55;   //!
   TBranch        *b_HLT_Mu55_prescale;   //!
   TBranch        *b_HLT_Mu300;   //!
   TBranch        *b_HLT_Mu300_prescale;   //!
   TBranch        *b_HLT_Mu350;   //!
   TBranch        *b_HLT_Mu350_prescale;   //!
   TBranch        *b_HLT_Ele105_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale;   //!
   TBranch        *b_HLT_Ele115_CaloIdVT_GsfTrkIdT;   //!
   TBranch        *b_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale;   //!
   TBranch        *b_HLT_IsoMu17_eta2p1;   //!
   TBranch        *b_HLT_IsoMu17_eta2p1_prescale;   //!
   TBranch        *b_HLT_IsoMu18;   //!
   TBranch        *b_HLT_IsoMu18_prescale;   //!
   TBranch        *b_HLT_OldIsoMu18;   //!
   TBranch        *b_HLT_OldIsoMu18_prescale;   //!
   TBranch        *b_HLT_IsoTkMu18;   //!
   TBranch        *b_HLT_IsoTkMu18_prescale;   //!
   TBranch        *b_HLT_OldIsoTkMu18;   //!
   TBranch        *b_HLT_OldIsoTkMu18_prescale;   //!
   TBranch        *b_HLT_IsoMu20;   //!
   TBranch        *b_HLT_IsoMu20_prescale;   //!
   TBranch        *b_HLT_IsoMu22;   //!
   TBranch        *b_HLT_IsoMu22_prescale;   //!
   TBranch        *b_HLT_IsoMu27;   //!
   TBranch        *b_HLT_IsoMu27_prescale;   //!
   TBranch        *b_HLT_IsoTkMu24_eta2p1;   //!
   TBranch        *b_HLT_IsoTkMu24_eta2p1_prescale;   //!
   TBranch        *b_HLT_IsoTkMu27;   //!
   TBranch        *b_HLT_IsoTkMu27_prescale;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu8_TrkIsoVVL_prescale;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL;   //!
   TBranch        *b_HLT_Mu17_TrkIsoVVL_prescale;   //!
   TBranch        *b_HLT_Ele22_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele22_eta2p1_WPLoose_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele22_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele22_eta2p1_WPTight_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele23_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele23_WPLoose_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele27_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele27_WPLoose_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf;   //!
   TBranch        *b_HLT_Ele27_eta2p1_WPLoose_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele27_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele27_eta2p1_WPTight_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele32_eta2p1_WPTight_Gsf;   //!
   TBranch        *b_HLT_Ele32_eta2p1_WPTight_Gsf_prescale;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_prescale;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_prescale;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL;   //!
   TBranch        *b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_prescale;   //!
   TBranch        *b_evt_NGoodVtx;   //!
   TBranch        *b_evt_NLep;   //!
   TBranch        *b_evt_NTopHad;   //!
   TBranch        *b_evt_NTopHadPreTag;   //!
   TBranch        *b_evt_NTopLep;   //!
   TBranch        *b_evt_NTop;   //!
   TBranch        *b_evt_LHA_PDF_ID;   //!
   TBranch        *b_evt_HtLep;   //!
   TBranch        *b_evt_HtTop;   //!
   TBranch        *b_evt_Ht;   //!
   TBranch        *b_evt_HtAll;   //!
   TBranch        *b_evt_HtEx;   //!
   TBranch        *b_evt_HtExFr;   //!
   TBranch        *b_evt_HtTopFr;   //!
   TBranch        *b_evt_TTHadDR;   //!
   TBranch        *b_evt_TTHadDPhi;   //!
   TBranch        *b_evt_TTHadDEta;   //!
   TBranch        *b_evt_TTHadMass;   //!
   TBranch        *b_evt_TTHadSumPt;   //!
   TBranch        *b_evt_TTHadPz;   //!
   TBranch        *b_evt_TTHadHz;   //!
   TBranch        *b_evt_TTHadDPz;   //!
   TBranch        *b_evt_TTHadMR;   //!
   TBranch        *b_evt_TTHadMTR;   //!
   TBranch        *b_evt_TTHadR;   //!
   TBranch        *b_evt_TTHadR2;   //!
   TBranch        *b_evt_MR;   //!
   TBranch        *b_evt_MTR;   //!
   TBranch        *b_evt_R;   //!
   TBranch        *b_evt_R2;   //!
   TBranch        *b_evt_AK8Puppi_MR;   //!
   TBranch        *b_evt_AK8Puppi_MTR;   //!
   TBranch        *b_evt_AK8Puppi_R;   //!
   TBranch        *b_evt_AK8Puppi_R2;   //!
   TBranch        *b_evt_AK4_MR;   //!
   TBranch        *b_evt_AK4_MTR;   //!
   TBranch        *b_evt_AK4_R;   //!
   TBranch        *b_evt_AK4_R2;   //!
   TBranch        *b_evt_XSec;   //!
   TBranch        *b_evt_NEvent_Corr;   //!
   TBranch        *b_evt_Lumi_Weight;   //!
   TBranch        *b_evt_Gen_Weight;   //!
   TBranch        *b_evt_Gen_Ht;   //!
   TBranch        *b_SUSY_Gluino_Mass;   //!
   TBranch        *b_SUSY_LSP_Mass;   //!
   TBranch        *b_gen_size;   //!
   TBranch        *b_gen_ID;   //!
   TBranch        *b_gen_Status;   //!
   TBranch        *b_gen_Mom0ID;   //!
   TBranch        *b_gen_Mom0Status;   //!
   TBranch        *b_gen_Mom1ID;   //!
   TBranch        *b_gen_Mom1Status;   //!
   TBranch        *b_gen_Dau0ID;   //!
   TBranch        *b_gen_Dau0Status;   //!
   TBranch        *b_gen_Dau1ID;   //!
   TBranch        *b_gen_Dau1Status;   //!
   TBranch        *b_jetAK8_HasNearGenTop;   //!
   TBranch        *b_jetAK8_NearGenTopIsHadronic;   //!
   TBranch        *b_jetAK8_NearGenWIsHadronic;   //!
   TBranch        *b_jetAK8_NearGenWToENu;   //!
   TBranch        *b_jetAK8_NearGenWToMuNu;   //!
   TBranch        *b_jetAK8_NearGenWToTauNu;   //!
   TBranch        *b_jetAK8_PassTopTag;   //!
   TBranch        *b_el_IsPartOfNearAK4Jet;   //!
   TBranch        *b_el_IsPartOfNearAK8Jet;   //!
   TBranch        *b_el_IsPartOfNearSubjet;   //!
   TBranch        *b_mu_IsPartOfNearAK4Jet;   //!
   TBranch        *b_mu_IsPartOfNearAK8Jet;   //!
   TBranch        *b_mu_IsPartOfNearSubjet;   //!
   TBranch        *b_scale_size;   //!
   TBranch        *b_scale_Weights;   //!
   TBranch        *b_pdf_size;   //!
   TBranch        *b_pdf_Weights;   //!
   TBranch        *b_alphas_size;   //!
   TBranch        *b_alphas_Weights;   //!
   TBranch        *b_gen_Pt;   //!
   TBranch        *b_gen_Eta;   //!
   TBranch        *b_gen_Phi;   //!
   TBranch        *b_gen_E;   //!
   TBranch        *b_gen_Mass;   //!
   TBranch        *b_gen_Charge;   //!
   TBranch        *b_jetAK8_DRNearGenTop;   //!
   TBranch        *b_jetAK8_DRNearGenWFromTop;   //!
   TBranch        *b_jetAK8_DRNearGenBFromTop;   //!
   TBranch        *b_jetAK8_DRNearGenLepFromSLTop;   //!
   TBranch        *b_jetAK8_DRNearGenNuFromSLTop;   //!
   TBranch        *b_jetAK8_PtNearGenTop;   //!
   TBranch        *b_jetAK8_PtNearGenBFromTop;   //!
   TBranch        *b_jetAK8_PtNearGenWFromTop;   //!
   TBranch        *b_jetAK8_PtNearGenLepFromSLTop;   //!
   TBranch        *b_jetAK8_PtNearGenNuFromSLTop;   //!
   TBranch        *b_el_DRNearGenEleFromSLTop;   //!
   TBranch        *b_el_PtNearGenEleFromSLTop;   //!
   TBranch        *b_el_PtNearGenTop;   //!
   TBranch        *b_el_LepAK4JetFrac;   //!
   TBranch        *b_el_LepAK8JetFrac;   //!
   TBranch        *b_el_LepSubjetFrac;   //!
   TBranch        *b_el_LepAK4JetMassDrop;   //!
   TBranch        *b_el_LepAK8JetMassDrop;   //!
   TBranch        *b_el_LepSubjetMassDrop;   //!
   TBranch        *b_el_AK4JetV1DR;   //!
   TBranch        *b_el_AK4JetV2DR;   //!
   TBranch        *b_el_AK4JetV3DR;   //!
   TBranch        *b_el_AK8JetV1DR;   //!
   TBranch        *b_el_AK8JetV2DR;   //!
   TBranch        *b_el_AK8JetV3DR;   //!
   TBranch        *b_el_SubjetV1DR;   //!
   TBranch        *b_el_SubjetV2DR;   //!
   TBranch        *b_el_SubjetV3DR;   //!
   TBranch        *b_el_AK4JetV1PtRel;   //!
   TBranch        *b_el_AK4JetV2PtRel;   //!
   TBranch        *b_el_AK4JetV3PtRel;   //!
   TBranch        *b_el_AK8JetV1PtRel;   //!
   TBranch        *b_el_AK8JetV2PtRel;   //!
   TBranch        *b_el_AK8JetV3PtRel;   //!
   TBranch        *b_el_SubjetV1PtRel;   //!
   TBranch        *b_el_SubjetV2PtRel;   //!
   TBranch        *b_el_SubjetV3PtRel;   //!
   TBranch        *b_mu_DRNearGenMuFromSLTop;   //!
   TBranch        *b_mu_PtNearGenMuFromSLTop;   //!
   TBranch        *b_mu_PtNearGenTop;   //!
   TBranch        *b_mu_LepAK4JetFrac;   //!
   TBranch        *b_mu_LepAK8JetFrac;   //!
   TBranch        *b_mu_LepSubjetFrac;   //!
   TBranch        *b_mu_LepAK4JetMassDrop;   //!
   TBranch        *b_mu_LepAK8JetMassDrop;   //!
   TBranch        *b_mu_LepSubjetMassDrop;   //!
   TBranch        *b_mu_AK4JetV1DR;   //!
   TBranch        *b_mu_AK4JetV2DR;   //!
   TBranch        *b_mu_AK4JetV3DR;   //!
   TBranch        *b_mu_AK8JetV1DR;   //!
   TBranch        *b_mu_AK8JetV2DR;   //!
   TBranch        *b_mu_AK8JetV3DR;   //!
   TBranch        *b_mu_SubjetV1DR;   //!
   TBranch        *b_mu_SubjetV2DR;   //!
   TBranch        *b_mu_SubjetV3DR;   //!
   TBranch        *b_mu_AK4JetV1PtRel;   //!
   TBranch        *b_mu_AK4JetV2PtRel;   //!
   TBranch        *b_mu_AK4JetV3PtRel;   //!
   TBranch        *b_mu_AK8JetV1PtRel;   //!
   TBranch        *b_mu_AK8JetV2PtRel;   //!
   TBranch        *b_mu_AK8JetV3PtRel;   //!
   TBranch        *b_mu_SubjetV1PtRel;   //!
   TBranch        *b_mu_SubjetV2PtRel;   //!
   TBranch        *b_mu_SubjetV3PtRel;   //!

   B2GTTree(TTree *tree=0);
   virtual ~B2GTTree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   TFile            *outFile;
   TTree            *outputTree;
   void             bookNtuple();
   Float_t              T_met_Phi;

   Int_t            smalltree_nlepton;
   Float_t         smalltree_lept_pt[1];   //[smalltree_nlepton]
   Float_t         smalltree_lept_eta[1];   //[smalltree_nlepton]
   Float_t         smalltree_lept_phi[1];   //[smalltree_nlepton]
   Float_t         smalltree_lept_iso[1];   //[smalltree_nlepton]
   Int_t           smalltree_lept_flav[1]; 

   Int_t           smalltree_njets;
   Float_t         smalltree_jet_pt[7];   //[smalltree_njets]
   Float_t         smalltree_jet_eta[7];   //[smalltree_njets]
   Float_t         smalltree_jet_phi[7];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri[7];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri_up[7];   //[smalltree_njets]
   Float_t         smalltree_jet_btagdiscri_down[7];   //[smalltree_njets]
   Int_t           smalltree_jet_flav[7];   //[smalltree_njets]
   Int_t           smalltree_jesup_njets;
   Float_t         smalltree_jet_jesup_pt[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_eta[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_phi[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jesup_btagdiscri[7];   //[smalltree_njets]
   Int_t           smalltree_jet_jesup_flav[7];   //[smalltree_njets]
   Int_t           smalltree_jesdown_njets;
   Float_t         smalltree_jet_jesdown_pt[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_eta[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_phi[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jesdown_btagdiscri[7];   //[smalltree_njets]
   Int_t           smalltree_jet_jesdown_flav[7];   //[smalltree_njets]
   Int_t           smalltree_jerup_njets;
   Float_t         smalltree_jet_jerup_pt[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_eta[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_phi[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jerup_btagdiscri[7];   //[smalltree_njets]
   Int_t           smalltree_jet_jerup_flav[7];   //[smalltree_njets]
   Int_t           smalltree_jerdown_njets;
   Float_t         smalltree_jet_jerdown_pt[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_eta[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_phi[7];   //[smalltree_njets]
   Float_t         smalltree_jet_jerdown_btagdiscri[7];   //[smalltree_njets]
   Int_t           smalltree_jet_jerdown_flav[7];   //[smalltree_njets]

   Float_t         smalltree_metnosmear_pt;
   Float_t         smalltree_metnosmear_phi;
   Float_t         smalltree_met_pt;
   Float_t         smalltree_met_phi;
   Float_t         smalltree_met_jesup_pt;
   Float_t         smalltree_met_jesup_phi;
   Float_t         smalltree_met_jesdown_pt;
   Float_t         smalltree_met_jesdown_phi;
   Float_t         smalltree_met_jerup_pt;
   Float_t         smalltree_met_jerup_phi;
   Float_t         smalltree_met_jerdown_pt;
   Float_t         smalltree_met_jerdown_phi;
   Float_t         smalltree_met_unclsup_pt;
   Float_t         smalltree_met_unclsup_phi;
   Float_t         smalltree_met_unclsdown_pt;
   Float_t         smalltree_met_unclsdown_phi;

   Float_t         smalltree_weight_trigup;
   Float_t         smalltree_weight_trigdown;
   Float_t         smalltree_weight_leptup;
   Float_t         smalltree_weight_leptdown;
   Float_t         smalltree_weight_PDFup;
   Float_t         smalltree_weight_PDFdown;
   Float_t         smalltree_weight_PUup;
   Float_t         smalltree_weight_PUdown;
   Float_t         smalltree_weight_toppt;

   Float_t         smalltree_evtweight;
   Int_t           smalltree_tmeme;
   Int_t           smalltree_nvertex;
};

#endif

void B2GTTree::bookNtuple(){
    outFile = new TFile("SmallTree.root","RECREATE","test");
    outputTree = new TTree("T", "test tree");
    outputTree->Branch("T_met_Phi", &T_met_Phi, "T_met_Phi/F");
    outputTree->Branch("smalltree_nlepton", &smalltree_nlepton, "smalltree_nlepton/I");
    outputTree->Branch("smalltree_lept_pt", &smalltree_lept_pt, "smalltree_lept_pt[1]/F");
    outputTree->Branch("smalltree_lept_eta", &smalltree_lept_eta, "smalltree_lept_eta[1]/F");
    outputTree->Branch("smalltree_lept_phi", &smalltree_lept_phi, "smalltree_lept_phi[1]/F");
    outputTree->Branch("smalltree_lept_iso", &smalltree_lept_iso, "smalltree_lept_iso[1]/F");
    outputTree->Branch("smalltree_lept_flav", &smalltree_lept_flav, "smalltree_lept_flav[1]/I");

    outputTree->Branch("smalltree_njets", &smalltree_njets, "smalltree_njets/I");
    outputTree->Branch("smalltree_jet_pt", &smalltree_jet_pt, "smalltree_jet_pt[7]/F");
    outputTree->Branch("smalltree_jet_eta", &smalltree_jet_eta, "smalltree_jet_eta[7]/F");
    outputTree->Branch("smalltree_jet_phi", &smalltree_jet_phi, "smalltree_jet_phi[7]/F");
    outputTree->Branch("smalltree_jet_btagdiscri", &smalltree_jet_btagdiscri, "smalltree_jet_btagdiscri[7]/F");
    outputTree->Branch("smalltree_jet_btagdiscri_up", &smalltree_jet_btagdiscri_up, "smalltree_jet_btagdiscri_up[7]/F");
    outputTree->Branch("smalltree_jet_btagdiscri_down", &smalltree_jet_btagdiscri_down, "smalltree_jet_btagdiscri_down[7]/F");
    outputTree->Branch("smalltree_jet_flav", &smalltree_jet_flav, "smalltree_jet_flav[7]/I");
    outputTree->Branch("smalltree_jesup_njets", &smalltree_jesup_njets, "smalltree_jesup_njets/I");
    outputTree->Branch("smalltree_jet_jesup_pt", &smalltree_jet_jesup_pt, "smalltree_jet_jesup_pt[7]/F");
    outputTree->Branch("smalltree_jet_jesup_eta", &smalltree_jet_jesup_eta, "smalltree_jet_jesup_eta[7]/F");
    outputTree->Branch("smalltree_jet_jesup_phi", &smalltree_jet_jesup_phi, "smalltree_jet_jesup_phi[7]/F");
    outputTree->Branch("smalltree_jet_jesup_btagdiscri", &smalltree_jet_jesup_btagdiscri, "smalltree_jet_jesup_btagdiscri[7]/F");
    outputTree->Branch("smalltree_jet_jesup_flav", &smalltree_jet_jesup_flav, "smalltree_jet_jesup_flav[7]/I");
    outputTree->Branch("smalltree_jesdown_njets", &smalltree_jesdown_njets, "smalltree_jesdown_njets/I");
    outputTree->Branch("smalltree_jet_jesdown_pt", &smalltree_jet_jesdown_pt, "smalltree_jet_jesdown_pt[7]/F");
    outputTree->Branch("smalltree_jet_jesdown_eta", &smalltree_jet_jesdown_eta, "smalltree_jet_jesdown_eta[7]/F");
    outputTree->Branch("smalltree_jet_jesdown_phi", &smalltree_jet_jesdown_phi, "smalltree_jet_jesdown_phi[7]/F");
    outputTree->Branch("smalltree_jet_jesdown_btagdiscri", &smalltree_jet_jesdown_btagdiscri, "smalltree_jet_jesdown_btagdiscri[7]/F");
    outputTree->Branch("smalltree_jet_jesdown_flav", &smalltree_jet_jesdown_flav, "smalltree_jet_jesdown_flav[7]/I");
    outputTree->Branch("smalltree_jerup_njets", &smalltree_jerup_njets, "smalltree_jerup_njets/I");
    outputTree->Branch("smalltree_jet_jerup_pt", &smalltree_jet_jerup_pt, "smalltree_jet_jerup_pt[7]/F");
    outputTree->Branch("smalltree_jet_jerup_eta", &smalltree_jet_jerup_eta, "smalltree_jet_jerup_eta[7]/F");
    outputTree->Branch("smalltree_jet_jerup_phi", &smalltree_jet_jerup_phi, "smalltree_jet_jerup_phi[7]/F");
    outputTree->Branch("smalltree_jet_jerup_btagdiscri", &smalltree_jet_jerup_btagdiscri, "smalltree_jet_jerup_btagdiscri[7]/F");
    outputTree->Branch("smalltree_jet_jerup_flav", &smalltree_jet_jerup_flav, "smalltree_jet_jerup_flav[7]/I");
    outputTree->Branch("smalltree_jerdown_njets", &smalltree_jerdown_njets, "smalltree_jerdown_njets/I");
    outputTree->Branch("smalltree_jet_jerdown_pt", &smalltree_jet_jerdown_pt, "smalltree_jet_jerdown_pt[7]/F");
    outputTree->Branch("smalltree_jet_jerdown_eta", &smalltree_jet_jerdown_eta, "smalltree_jet_jerdown_eta[7]/F");
    outputTree->Branch("smalltree_jet_jerdown_phi", &smalltree_jet_jerdown_phi, "smalltree_jet_jerdown_phi[7]/F");
    outputTree->Branch("smalltree_jet_jerdown_btagdiscri", &smalltree_jet_jerdown_btagdiscri, "smalltree_jet_jerdown_btagdiscri[7]/F");
    outputTree->Branch("smalltree_jet_jerdown_flav", &smalltree_jet_jerdown_flav, "smalltree_jet_jerdown_flav[7]/I");

    outputTree->Branch("smalltree_metnosmear_pt", &smalltree_metnosmear_pt, "smalltree_metnosmear_pt/F");
    outputTree->Branch("smalltree_metnosmear_phi", &smalltree_metnosmear_phi, "smalltree_metnosmear_phi/F");
    outputTree->Branch("smalltree_met_pt", &smalltree_met_pt, "smalltree_met_pt/F");
    outputTree->Branch("smalltree_met_phi", &smalltree_met_phi, "smalltree_met_phi/F");
    outputTree->Branch("smalltree_met_jesup_pt", &smalltree_met_jesup_pt, "smalltree_met_jesup_pt/F");
    outputTree->Branch("smalltree_met_jesup_phi", &smalltree_met_jesup_phi, "smalltree_met_jesup_phi/F");
    outputTree->Branch("smalltree_met_jesdown_pt", &smalltree_met_jesdown_pt, "smalltree_met_jesdown_pt/F");
    outputTree->Branch("smalltree_met_jesdown_phi", &smalltree_met_jesdown_phi, "smalltree_met_jesdown_phi/F");
    outputTree->Branch("smalltree_met_jerup_pt", &smalltree_met_jerup_pt, "smalltree_met_jerup_pt/F");
    outputTree->Branch("smalltree_met_jerup_phi", &smalltree_met_jerup_phi, "smalltree_met_jerup_phi/F");
    outputTree->Branch("smalltree_met_jerdown_pt", &smalltree_met_jerdown_pt, "smalltree_met_jerdown_pt/F");
    outputTree->Branch("smalltree_met_jerdown_phi", &smalltree_met_jerdown_phi, "smalltree_met_jerdown_phi/F");
    outputTree->Branch("smalltree_met_unclsup_pt", &smalltree_met_unclsup_pt, "smalltree_met_unclsup_pt/F");
    outputTree->Branch("smalltree_met_unclsup_phi", &smalltree_met_unclsup_phi, "smalltree_met_unclsup_phi/F");
    outputTree->Branch("smalltree_met_unclsdown_pt", &smalltree_met_unclsdown_pt, "smalltree_met_unclsdown_pt/F");
    outputTree->Branch("smalltree_met_unclsdown_phi", &smalltree_met_unclsdown_phi, "smalltree_met_unclsdown_phi/F");

    outputTree->Branch("smalltree_weight_trigup", &smalltree_weight_trigup, "smalltree_weight_trigup/F");
    outputTree->Branch("smalltree_weight_trigdown", &smalltree_weight_trigdown, "smalltree_weight_trigdown/F");
    outputTree->Branch("smalltree_weight_leptup", &smalltree_weight_leptup, "smalltree_weight_leptup/F");
    outputTree->Branch("smalltree_weight_leptdown", &smalltree_weight_leptdown, "smalltree_weight_leptdown/F");
    outputTree->Branch("smalltree_weight_PDFup", &smalltree_weight_PDFup, "smalltree_weight_PDFup/F");
    outputTree->Branch("smalltree_weight_PDFdown", &smalltree_weight_PDFdown, "smalltree_weight_PDFdown/F");
    outputTree->Branch("smalltree_weight_PUup", &smalltree_weight_PUup, "smalltree_weight_PUup/F");
    outputTree->Branch("smalltree_weight_PUdown", &smalltree_weight_PUdown, "smalltree_weight_PUdown/F");
    outputTree->Branch("smalltree_weight_toppt", &smalltree_weight_toppt, "smalltree_weight_toppt/F");

    outputTree->Branch("smalltree_evtweight", &smalltree_evtweight, "smalltree_evtweight/F");
    outputTree->Branch("smalltree_tmeme", &smalltree_tmeme, "smalltree_tmeme/I");
    outputTree->Branch("smalltree_nvertex", &smalltree_nvertex, "smalltree_nvertex/I");
    
}


#ifdef B2GTTree_cxx
B2GTTree::B2GTTree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("B2GTTreeNtupleExtra_MC_25ns_76X.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("B2GTTreeNtupleExtra_MC_25ns_76X.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("B2GTTreeNtupleExtra_MC_25ns_76X.root:/B2GTTreeMaker");
      dir->GetObject("B2GTree",tree);

   }
   Init(tree);
}

B2GTTree::~B2GTTree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t B2GTTree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t B2GTTree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void B2GTTree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   jetAK4_Keys = 0;
   jetAK4Puppi_Keys = 0;
   jetAK8_Keys = 0;
   jetAK8Puppi_Keys = 0;
   subjetAK8_Keys = 0;
   subjetAK8Puppi_Keys = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("evt_npv", &evt_npv, &b_evt_npv);
   fChain->SetBranchAddress("evt_vx", &evt_vx, &b_evt_vx);
   fChain->SetBranchAddress("evt_vy", &evt_vy, &b_evt_vy);
   fChain->SetBranchAddress("evt_vz", &evt_vz, &b_evt_vz);
   fChain->SetBranchAddress("evt_RunNumber", &evt_RunNumber, &b_evt_RunNumber);
   fChain->SetBranchAddress("evt_LumiBlock", &evt_LumiBlock, &b_evt_LumiBlock);
   fChain->SetBranchAddress("evt_EventNumber", &evt_EventNumber, &b_evt_EventNumber);
   fChain->SetBranchAddress("evt_rho", &evt_rho, &b_evt_rho);
   fChain->SetBranchAddress("pu_NtrueInt", &pu_NtrueInt, &b_pu_NtrueInt);
   fChain->SetBranchAddress("pu_size", &pu_size, &b_pu_size);
   fChain->SetBranchAddress("pu_BX", pu_BX, &b_pu_BX);
   fChain->SetBranchAddress("pu_NInt", pu_NInt, &b_pu_NInt);
   fChain->SetBranchAddress("vtx_size", &vtx_size, &b_vtx_size);
   fChain->SetBranchAddress("vtx_ndof", vtx_ndof, &b_vtx_ndof);
   fChain->SetBranchAddress("vtx_z", vtx_z, &b_vtx_z);
   fChain->SetBranchAddress("vtx_rho", vtx_rho, &b_vtx_rho);
   fChain->SetBranchAddress("vtx_chi", vtx_chi, &b_vtx_chi);
   fChain->SetBranchAddress("met_size", &met_size, &b_met_size);
   fChain->SetBranchAddress("met_Pt", met_Pt, &b_met_Pt);
   fChain->SetBranchAddress("met_Px", met_Px, &b_met_Px);
   fChain->SetBranchAddress("met_Py", met_Py, &b_met_Py);
   fChain->SetBranchAddress("met_Phi", met_Phi, &b_met_Phi);
   fChain->SetBranchAddress("met_uncorPt", met_uncorPt, &b_met_uncorPt);
   fChain->SetBranchAddress("met_uncorPhi", met_uncorPhi, &b_met_uncorPhi);
   fChain->SetBranchAddress("met_uncorSumEt", met_uncorSumEt, &b_met_uncorSumEt);
   fChain->SetBranchAddress("el_size", &el_size, &b_el_size);
   fChain->SetBranchAddress("el_Pt", el_Pt, &b_el_Pt);
   fChain->SetBranchAddress("el_Eta", el_Eta, &b_el_Eta);
   fChain->SetBranchAddress("el_Phi", el_Phi, &b_el_Phi);
   fChain->SetBranchAddress("el_E", el_E, &b_el_E);
   fChain->SetBranchAddress("el_Charge", el_Charge, &b_el_Charge);
   fChain->SetBranchAddress("el_Key", el_Key, &b_el_Key);
   fChain->SetBranchAddress("el_Iso03", el_Iso03, &b_el_Iso03);
   fChain->SetBranchAddress("el_Iso03db", el_Iso03db, &b_el_Iso03db);
   fChain->SetBranchAddress("el_MiniIso", el_MiniIso, &b_el_MiniIso);
   fChain->SetBranchAddress("el_rho", el_rho, &b_el_rho);
   fChain->SetBranchAddress("el_EA", el_EA, &b_el_EA);
   fChain->SetBranchAddress("el_sumChargedHadronPt", el_sumChargedHadronPt, &b_el_sumChargedHadronPt);
   fChain->SetBranchAddress("el_sumNeutralHadronEt", el_sumNeutralHadronEt, &b_el_sumNeutralHadronEt);
   fChain->SetBranchAddress("el_sumPhotonEt", el_sumPhotonEt, &b_el_sumPhotonEt);
   fChain->SetBranchAddress("el_sumPUPt", el_sumPUPt, &b_el_sumPUPt);
   fChain->SetBranchAddress("el_D0", el_D0, &b_el_D0);
   fChain->SetBranchAddress("el_Dz", el_Dz, &b_el_Dz);
   fChain->SetBranchAddress("el_dEtaIn", el_dEtaIn, &b_el_dEtaIn);
   fChain->SetBranchAddress("el_dPhiIn", el_dPhiIn, &b_el_dPhiIn);
   fChain->SetBranchAddress("el_HoE", el_HoE, &b_el_HoE);
   fChain->SetBranchAddress("el_full5x5siee", el_full5x5siee, &b_el_full5x5siee);
   fChain->SetBranchAddress("el_ooEmooP", el_ooEmooP, &b_el_ooEmooP);
   fChain->SetBranchAddress("el_missHits", el_missHits, &b_el_missHits);
   fChain->SetBranchAddress("el_hasMatchedConVeto", el_hasMatchedConVeto, &b_el_hasMatchedConVeto);
   fChain->SetBranchAddress("el_vidVeto", el_vidVeto, &b_el_vidVeto);
   fChain->SetBranchAddress("el_vidLoose", el_vidLoose, &b_el_vidLoose);
   fChain->SetBranchAddress("el_vidTight", el_vidTight, &b_el_vidTight);
   fChain->SetBranchAddress("el_vidMedium", el_vidMedium, &b_el_vidMedium);
   fChain->SetBranchAddress("el_vidHEEP", el_vidHEEP, &b_el_vidHEEP);
   fChain->SetBranchAddress("el_vidHEEPnoiso", el_vidHEEPnoiso, &b_el_vidHEEPnoiso);
   fChain->SetBranchAddress("el_SCEta", el_SCEta, &b_el_SCEta);
   fChain->SetBranchAddress("el_SCPhi", el_SCPhi, &b_el_SCPhi);
   fChain->SetBranchAddress("mu_size", &mu_size, &b_mu_size);
   fChain->SetBranchAddress("mu_Pt", mu_Pt, &b_mu_Pt);
   fChain->SetBranchAddress("mu_Eta", mu_Eta, &b_mu_Eta);
   fChain->SetBranchAddress("mu_Phi", mu_Phi, &b_mu_Phi);
   fChain->SetBranchAddress("mu_E", mu_E, &b_mu_E);
   fChain->SetBranchAddress("mu_Charge", mu_Charge, &b_mu_Charge);
   fChain->SetBranchAddress("mu_Key", mu_Key, &b_mu_Key);
   fChain->SetBranchAddress("mu_Iso04", mu_Iso04, &b_mu_Iso04);
   fChain->SetBranchAddress("mu_MiniIso", mu_MiniIso, &b_mu_MiniIso);
   fChain->SetBranchAddress("mu_D0", mu_D0, &b_mu_D0);
   fChain->SetBranchAddress("mu_D0err", mu_D0err, &b_mu_D0err);
   fChain->SetBranchAddress("mu_Dxy", mu_Dxy, &b_mu_Dxy);
   fChain->SetBranchAddress("mu_Dxyerr", mu_Dxyerr, &b_mu_Dxyerr);
   fChain->SetBranchAddress("mu_Dz", mu_Dz, &b_mu_Dz);
   fChain->SetBranchAddress("mu_Dzerr", mu_Dzerr, &b_mu_Dzerr);
   fChain->SetBranchAddress("mu_IsSoftMuon", mu_IsSoftMuon, &b_mu_IsSoftMuon);
   fChain->SetBranchAddress("mu_IsLooseMuon", mu_IsLooseMuon, &b_mu_IsLooseMuon);
   fChain->SetBranchAddress("mu_IsMediumMuon", mu_IsMediumMuon, &b_mu_IsMediumMuon);
   fChain->SetBranchAddress("mu_IsTightMuon", mu_IsTightMuon, &b_mu_IsTightMuon);
   fChain->SetBranchAddress("mu_IsHighPtMuon", mu_IsHighPtMuon, &b_mu_IsHighPtMuon);
   fChain->SetBranchAddress("mu_IsPFMuon", mu_IsPFMuon, &b_mu_IsPFMuon);
   fChain->SetBranchAddress("mu_IsGlobalMuon", mu_IsGlobalMuon, &b_mu_IsGlobalMuon);
   fChain->SetBranchAddress("mu_IsTrackerMuon", mu_IsTrackerMuon, &b_mu_IsTrackerMuon);
   fChain->SetBranchAddress("mu_GlbTrkNormChi2", mu_GlbTrkNormChi2, &b_mu_GlbTrkNormChi2);
   fChain->SetBranchAddress("mu_NumberValidMuonHits", mu_NumberValidMuonHits, &b_mu_NumberValidMuonHits);
   fChain->SetBranchAddress("mu_NumberMatchedStations", mu_NumberMatchedStations, &b_mu_NumberMatchedStations);
   fChain->SetBranchAddress("mu_NumberValidPixelHits", mu_NumberValidPixelHits, &b_mu_NumberValidPixelHits);
   fChain->SetBranchAddress("mu_NumberTrackerLayers", mu_NumberTrackerLayers, &b_mu_NumberTrackerLayers);
   fChain->SetBranchAddress("mu_NumberOfValidTrackerHits", mu_NumberOfValidTrackerHits, &b_mu_NumberOfValidTrackerHits);
   fChain->SetBranchAddress("mu_NumberOfPixelLayers", mu_NumberOfPixelLayers, &b_mu_NumberOfPixelLayers);
   fChain->SetBranchAddress("mu_InTrkNormChi2", mu_InTrkNormChi2, &b_mu_InTrkNormChi2);
   fChain->SetBranchAddress("mu_SumChargedHadronPt", mu_SumChargedHadronPt, &b_mu_SumChargedHadronPt);
   fChain->SetBranchAddress("mu_SumNeutralHadronPt", mu_SumNeutralHadronPt, &b_mu_SumNeutralHadronPt);
   fChain->SetBranchAddress("mu_SumPhotonPt", mu_SumPhotonPt, &b_mu_SumPhotonPt);
   fChain->SetBranchAddress("mu_SumPUPt", mu_SumPUPt, &b_mu_SumPUPt);
   fChain->SetBranchAddress("mu_GenMuonY", mu_GenMuonY, &b_mu_GenMuonY);
   fChain->SetBranchAddress("mu_GenMuonEta", mu_GenMuonEta, &b_mu_GenMuonEta);
   fChain->SetBranchAddress("mu_GenMuonPhi", mu_GenMuonPhi, &b_mu_GenMuonPhi);
   fChain->SetBranchAddress("mu_GenMuonPt", mu_GenMuonPt, &b_mu_GenMuonPt);
   fChain->SetBranchAddress("mu_GenMuonE", mu_GenMuonE, &b_mu_GenMuonE);
   fChain->SetBranchAddress("mu_GenMuonCharge", mu_GenMuonCharge, &b_mu_GenMuonCharge);
   fChain->SetBranchAddress("jetAK4_size", &jetAK4_size, &b_jetAK4_size);
   fChain->SetBranchAddress("jetAK4_Pt", jetAK4_Pt, &b_jetAK4_Pt);
   fChain->SetBranchAddress("jetAK4_Eta", jetAK4_Eta, &b_jetAK4_Eta);
   fChain->SetBranchAddress("jetAK4_Phi", jetAK4_Phi, &b_jetAK4_Phi);
   fChain->SetBranchAddress("jetAK4_E", jetAK4_E, &b_jetAK4_E);
   fChain->SetBranchAddress("jetAK4_Charge", jetAK4_Charge, &b_jetAK4_Charge);
   fChain->SetBranchAddress("jetAK4_CSVv2", jetAK4_CSVv2, &b_jetAK4_CSVv2);
   fChain->SetBranchAddress("jetAK4_DoubleBAK8", jetAK4_DoubleBAK8, &b_jetAK4_DoubleBAK8);
   fChain->SetBranchAddress("jetAK4_DoubleBCA15", jetAK4_DoubleBCA15, &b_jetAK4_DoubleBCA15);
   fChain->SetBranchAddress("jetAK4_CMVAv2", jetAK4_CMVAv2, &b_jetAK4_CMVAv2);
   fChain->SetBranchAddress("jetAK4_CvsL", jetAK4_CvsL, &b_jetAK4_CvsL);
   fChain->SetBranchAddress("jetAK4_CvsB", jetAK4_CvsB, &b_jetAK4_CvsB);
   fChain->SetBranchAddress("jetAK4_CMVA", jetAK4_CMVA, &b_jetAK4_CMVA);
   fChain->SetBranchAddress("jetAK4_GenPartonY", jetAK4_GenPartonY, &b_jetAK4_GenPartonY);
   fChain->SetBranchAddress("jetAK4_GenPartonEta", jetAK4_GenPartonEta, &b_jetAK4_GenPartonEta);
   fChain->SetBranchAddress("jetAK4_GenPartonPhi", jetAK4_GenPartonPhi, &b_jetAK4_GenPartonPhi);
   fChain->SetBranchAddress("jetAK4_GenPartonPt", jetAK4_GenPartonPt, &b_jetAK4_GenPartonPt);
   fChain->SetBranchAddress("jetAK4_GenPartonE", jetAK4_GenPartonE, &b_jetAK4_GenPartonE);
   fChain->SetBranchAddress("jetAK4_GenPartonCharge", jetAK4_GenPartonCharge, &b_jetAK4_GenPartonCharge);
   fChain->SetBranchAddress("jetAK4_PartonFlavour", jetAK4_PartonFlavour, &b_jetAK4_PartonFlavour);
   fChain->SetBranchAddress("jetAK4_HadronFlavour", jetAK4_HadronFlavour, &b_jetAK4_HadronFlavour);
   fChain->SetBranchAddress("jetAK4_GenJetY", jetAK4_GenJetY, &b_jetAK4_GenJetY);
   fChain->SetBranchAddress("jetAK4_GenJetEta", jetAK4_GenJetEta, &b_jetAK4_GenJetEta);
   fChain->SetBranchAddress("jetAK4_GenJetPhi", jetAK4_GenJetPhi, &b_jetAK4_GenJetPhi);
   fChain->SetBranchAddress("jetAK4_GenJetPt", jetAK4_GenJetPt, &b_jetAK4_GenJetPt);
   fChain->SetBranchAddress("jetAK4_GenJetE", jetAK4_GenJetE, &b_jetAK4_GenJetE);
   fChain->SetBranchAddress("jetAK4_GenJetCharge", jetAK4_GenJetCharge, &b_jetAK4_GenJetCharge);
   fChain->SetBranchAddress("jetAK4_muonMultiplicity", jetAK4_muonMultiplicity, &b_jetAK4_muonMultiplicity);
   fChain->SetBranchAddress("jetAK4_PhotonEnergy", jetAK4_PhotonEnergy, &b_jetAK4_PhotonEnergy);
   fChain->SetBranchAddress("jetAK4_ElectronEnergy", jetAK4_ElectronEnergy, &b_jetAK4_ElectronEnergy);
   fChain->SetBranchAddress("jetAK4_MuonEnergy", jetAK4_MuonEnergy, &b_jetAK4_MuonEnergy);
   fChain->SetBranchAddress("jetAK4_HFHadronEnergy", jetAK4_HFHadronEnergy, &b_jetAK4_HFHadronEnergy);
   fChain->SetBranchAddress("jetAK4_HFEMEnergy", jetAK4_HFEMEnergy, &b_jetAK4_HFEMEnergy);
   fChain->SetBranchAddress("jetAK4_ChargedHadronMultiplicity", jetAK4_ChargedHadronMultiplicity, &b_jetAK4_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("jetAK4_numberOfDaughters", jetAK4_numberOfDaughters, &b_jetAK4_numberOfDaughters);
   fChain->SetBranchAddress("jetAK4_neutralHadronMultiplicity", jetAK4_neutralHadronMultiplicity, &b_jetAK4_neutralHadronMultiplicity);
   fChain->SetBranchAddress("jetAK4_neutralHadronEnergy", jetAK4_neutralHadronEnergy, &b_jetAK4_neutralHadronEnergy);
   fChain->SetBranchAddress("jetAK4_neutralEmEnergy", jetAK4_neutralEmEnergy, &b_jetAK4_neutralEmEnergy);
   fChain->SetBranchAddress("jetAK4_chargedEmEnergy", jetAK4_chargedEmEnergy, &b_jetAK4_chargedEmEnergy);
   fChain->SetBranchAddress("jetAK4_chargedHadronEnergy", jetAK4_chargedHadronEnergy, &b_jetAK4_chargedHadronEnergy);
   fChain->SetBranchAddress("jetAK4_photonMultiplicity", jetAK4_photonMultiplicity, &b_jetAK4_photonMultiplicity);
   fChain->SetBranchAddress("jetAK4_electronMultiplicity", jetAK4_electronMultiplicity, &b_jetAK4_electronMultiplicity);
   fChain->SetBranchAddress("jetAK4_HFHadronMultiplicity", jetAK4_HFHadronMultiplicity, &b_jetAK4_HFHadronMultiplicity);
   fChain->SetBranchAddress("jetAK4_HFEMMultiplicity", jetAK4_HFEMMultiplicity, &b_jetAK4_HFEMMultiplicity);
   fChain->SetBranchAddress("jetAK4_ChargeMuEnergy", jetAK4_ChargeMuEnergy, &b_jetAK4_ChargeMuEnergy);
   fChain->SetBranchAddress("jetAK4_neutralMultiplicity", jetAK4_neutralMultiplicity, &b_jetAK4_neutralMultiplicity);
   fChain->SetBranchAddress("jetAK4_neutralHadronEnergyFrac", jetAK4_neutralHadronEnergyFrac, &b_jetAK4_neutralHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK4_neutralEmEnergyFrac", jetAK4_neutralEmEnergyFrac, &b_jetAK4_neutralEmEnergyFrac);
   fChain->SetBranchAddress("jetAK4_chargedHadronEnergyFrac", jetAK4_chargedHadronEnergyFrac, &b_jetAK4_chargedHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK4_muonEnergyFrac", jetAK4_muonEnergyFrac, &b_jetAK4_muonEnergyFrac);
   fChain->SetBranchAddress("jetAK4_chargedEmEnergyFrac", jetAK4_chargedEmEnergyFrac, &b_jetAK4_chargedEmEnergyFrac);
   fChain->SetBranchAddress("jetAK4_chargedMultiplicity", jetAK4_chargedMultiplicity, &b_jetAK4_chargedMultiplicity);
   fChain->SetBranchAddress("jetAK4_NumConstituents", jetAK4_NumConstituents, &b_jetAK4_NumConstituents);
   fChain->SetBranchAddress("jetAK4_jecFactor0", jetAK4_jecFactor0, &b_jetAK4_jecFactor0);
   fChain->SetBranchAddress("jetAK4_jecFactorL3Absolute", jetAK4_jecFactorL3Absolute, &b_jetAK4_jecFactorL3Absolute);
   fChain->SetBranchAddress("jetAK4_jetArea", jetAK4_jetArea, &b_jetAK4_jetArea);
   fChain->SetBranchAddress("jetAK4_nSV", jetAK4_nSV, &b_jetAK4_nSV);
   fChain->SetBranchAddress("jetAK4_SV0mass", jetAK4_SV0mass, &b_jetAK4_SV0mass);
   fChain->SetBranchAddress("jetAK4_SV1mass", jetAK4_SV1mass, &b_jetAK4_SV1mass);
   fChain->SetBranchAddress("jetAK4_QGL", jetAK4_QGL, &b_jetAK4_QGL);
   fChain->SetBranchAddress("jetAK4_jecUncertainty", jetAK4_jecUncertainty, &b_jetAK4_jecUncertainty);
   fChain->SetBranchAddress("jetAK4_PtResolution", jetAK4_PtResolution, &b_jetAK4_PtResolution);
   fChain->SetBranchAddress("jetAK4_JERSF", jetAK4_JERSF, &b_jetAK4_JERSF);
   fChain->SetBranchAddress("jetAK4_JERSFUp", jetAK4_JERSFUp, &b_jetAK4_JERSFUp);
   fChain->SetBranchAddress("jetAK4_JERSFDown", jetAK4_JERSFDown, &b_jetAK4_JERSFDown);
   fChain->SetBranchAddress("jetAK4_SmearedPt", jetAK4_SmearedPt, &b_jetAK4_SmearedPt);
   fChain->SetBranchAddress("jetAK4_SmearedPEta", jetAK4_SmearedPEta, &b_jetAK4_SmearedPEta);
   fChain->SetBranchAddress("jetAK4_SmearedPhi", jetAK4_SmearedPhi, &b_jetAK4_SmearedPhi);
   fChain->SetBranchAddress("jetAK4_SmearedE", jetAK4_SmearedE, &b_jetAK4_SmearedE);
   fChain->SetBranchAddress("jetAK4_Keys", &jetAK4_Keys, &b_jetAK4_Keys);
   fChain->SetBranchAddress("jetAK4Puppi_size", &jetAK4Puppi_size, &b_jetAK4Puppi_size);
   fChain->SetBranchAddress("jetAK4Puppi_Pt", jetAK4Puppi_Pt, &b_jetAK4Puppi_Pt);
   fChain->SetBranchAddress("jetAK4Puppi_Eta", jetAK4Puppi_Eta, &b_jetAK4Puppi_Eta);
   fChain->SetBranchAddress("jetAK4Puppi_Phi", jetAK4Puppi_Phi, &b_jetAK4Puppi_Phi);
   fChain->SetBranchAddress("jetAK4Puppi_E", jetAK4Puppi_E, &b_jetAK4Puppi_E);
   fChain->SetBranchAddress("jetAK4Puppi_Charge", jetAK4Puppi_Charge, &b_jetAK4Puppi_Charge);
   fChain->SetBranchAddress("jetAK4Puppi_CSVv2", jetAK4Puppi_CSVv2, &b_jetAK4Puppi_CSVv2);
   fChain->SetBranchAddress("jetAK4Puppi_DoubleBAK8", jetAK4Puppi_DoubleBAK8, &b_jetAK4Puppi_DoubleBAK8);
   fChain->SetBranchAddress("jetAK4Puppi_DoubleBCA15", jetAK4Puppi_DoubleBCA15, &b_jetAK4Puppi_DoubleBCA15);
   fChain->SetBranchAddress("jetAK4Puppi_CMVAv2", jetAK4Puppi_CMVAv2, &b_jetAK4Puppi_CMVAv2);
   fChain->SetBranchAddress("jetAK4Puppi_CvsL", jetAK4Puppi_CvsL, &b_jetAK4Puppi_CvsL);
   fChain->SetBranchAddress("jetAK4Puppi_CvsB", jetAK4Puppi_CvsB, &b_jetAK4Puppi_CvsB);
   fChain->SetBranchAddress("jetAK4Puppi_CMVA", jetAK4Puppi_CMVA, &b_jetAK4Puppi_CMVA);
   fChain->SetBranchAddress("jetAK4Puppi_GenPartonY", jetAK4Puppi_GenPartonY, &b_jetAK4Puppi_GenPartonY);
   fChain->SetBranchAddress("jetAK4Puppi_GenPartonEta", jetAK4Puppi_GenPartonEta, &b_jetAK4Puppi_GenPartonEta);
   fChain->SetBranchAddress("jetAK4Puppi_GenPartonPhi", jetAK4Puppi_GenPartonPhi, &b_jetAK4Puppi_GenPartonPhi);
   fChain->SetBranchAddress("jetAK4Puppi_GenPartonPt", jetAK4Puppi_GenPartonPt, &b_jetAK4Puppi_GenPartonPt);
   fChain->SetBranchAddress("jetAK4Puppi_GenPartonE", jetAK4Puppi_GenPartonE, &b_jetAK4Puppi_GenPartonE);
   fChain->SetBranchAddress("jetAK4Puppi_GenPartonCharge", jetAK4Puppi_GenPartonCharge, &b_jetAK4Puppi_GenPartonCharge);
   fChain->SetBranchAddress("jetAK4Puppi_PartonFlavour", jetAK4Puppi_PartonFlavour, &b_jetAK4Puppi_PartonFlavour);
   fChain->SetBranchAddress("jetAK4Puppi_HadronFlavour", jetAK4Puppi_HadronFlavour, &b_jetAK4Puppi_HadronFlavour);
   fChain->SetBranchAddress("jetAK4Puppi_GenJetY", jetAK4Puppi_GenJetY, &b_jetAK4Puppi_GenJetY);
   fChain->SetBranchAddress("jetAK4Puppi_GenJetEta", jetAK4Puppi_GenJetEta, &b_jetAK4Puppi_GenJetEta);
   fChain->SetBranchAddress("jetAK4Puppi_GenJetPhi", jetAK4Puppi_GenJetPhi, &b_jetAK4Puppi_GenJetPhi);
   fChain->SetBranchAddress("jetAK4Puppi_GenJetPt", jetAK4Puppi_GenJetPt, &b_jetAK4Puppi_GenJetPt);
   fChain->SetBranchAddress("jetAK4Puppi_GenJetE", jetAK4Puppi_GenJetE, &b_jetAK4Puppi_GenJetE);
   fChain->SetBranchAddress("jetAK4Puppi_GenJetCharge", jetAK4Puppi_GenJetCharge, &b_jetAK4Puppi_GenJetCharge);
   fChain->SetBranchAddress("jetAK4Puppi_muonMultiplicity", jetAK4Puppi_muonMultiplicity, &b_jetAK4Puppi_muonMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_PhotonEnergy", jetAK4Puppi_PhotonEnergy, &b_jetAK4Puppi_PhotonEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_ElectronEnergy", jetAK4Puppi_ElectronEnergy, &b_jetAK4Puppi_ElectronEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_MuonEnergy", jetAK4Puppi_MuonEnergy, &b_jetAK4Puppi_MuonEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_HFHadronEnergy", jetAK4Puppi_HFHadronEnergy, &b_jetAK4Puppi_HFHadronEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_HFEMEnergy", jetAK4Puppi_HFEMEnergy, &b_jetAK4Puppi_HFEMEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_ChargedHadronMultiplicity", jetAK4Puppi_ChargedHadronMultiplicity, &b_jetAK4Puppi_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_numberOfDaughters", jetAK4Puppi_numberOfDaughters, &b_jetAK4Puppi_numberOfDaughters);
   fChain->SetBranchAddress("jetAK4Puppi_neutralHadronMultiplicity", jetAK4Puppi_neutralHadronMultiplicity, &b_jetAK4Puppi_neutralHadronMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_neutralHadronEnergy", jetAK4Puppi_neutralHadronEnergy, &b_jetAK4Puppi_neutralHadronEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_neutralEmEnergy", jetAK4Puppi_neutralEmEnergy, &b_jetAK4Puppi_neutralEmEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_chargedEmEnergy", jetAK4Puppi_chargedEmEnergy, &b_jetAK4Puppi_chargedEmEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_chargedHadronEnergy", jetAK4Puppi_chargedHadronEnergy, &b_jetAK4Puppi_chargedHadronEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_photonMultiplicity", jetAK4Puppi_photonMultiplicity, &b_jetAK4Puppi_photonMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_electronMultiplicity", jetAK4Puppi_electronMultiplicity, &b_jetAK4Puppi_electronMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_HFHadronMultiplicity", jetAK4Puppi_HFHadronMultiplicity, &b_jetAK4Puppi_HFHadronMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_HFEMMultiplicity", jetAK4Puppi_HFEMMultiplicity, &b_jetAK4Puppi_HFEMMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_ChargeMuEnergy", jetAK4Puppi_ChargeMuEnergy, &b_jetAK4Puppi_ChargeMuEnergy);
   fChain->SetBranchAddress("jetAK4Puppi_neutralMultiplicity", jetAK4Puppi_neutralMultiplicity, &b_jetAK4Puppi_neutralMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_neutralHadronEnergyFrac", jetAK4Puppi_neutralHadronEnergyFrac, &b_jetAK4Puppi_neutralHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK4Puppi_neutralEmEnergyFrac", jetAK4Puppi_neutralEmEnergyFrac, &b_jetAK4Puppi_neutralEmEnergyFrac);
   fChain->SetBranchAddress("jetAK4Puppi_chargedHadronEnergyFrac", jetAK4Puppi_chargedHadronEnergyFrac, &b_jetAK4Puppi_chargedHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK4Puppi_muonEnergyFrac", jetAK4Puppi_muonEnergyFrac, &b_jetAK4Puppi_muonEnergyFrac);
   fChain->SetBranchAddress("jetAK4Puppi_chargedEmEnergyFrac", jetAK4Puppi_chargedEmEnergyFrac, &b_jetAK4Puppi_chargedEmEnergyFrac);
   fChain->SetBranchAddress("jetAK4Puppi_chargedMultiplicity", jetAK4Puppi_chargedMultiplicity, &b_jetAK4Puppi_chargedMultiplicity);
   fChain->SetBranchAddress("jetAK4Puppi_NumConstituents", jetAK4Puppi_NumConstituents, &b_jetAK4Puppi_NumConstituents);
   fChain->SetBranchAddress("jetAK4Puppi_jecFactor0", jetAK4Puppi_jecFactor0, &b_jetAK4Puppi_jecFactor0);
   fChain->SetBranchAddress("jetAK4Puppi_jecFactorL3Absolute", jetAK4Puppi_jecFactorL3Absolute, &b_jetAK4Puppi_jecFactorL3Absolute);
   fChain->SetBranchAddress("jetAK4Puppi_jetArea", jetAK4Puppi_jetArea, &b_jetAK4Puppi_jetArea);
   fChain->SetBranchAddress("jetAK4Puppi_nSV", jetAK4Puppi_nSV, &b_jetAK4Puppi_nSV);
   fChain->SetBranchAddress("jetAK4Puppi_SV0mass", jetAK4Puppi_SV0mass, &b_jetAK4Puppi_SV0mass);
   fChain->SetBranchAddress("jetAK4Puppi_SV1mass", jetAK4Puppi_SV1mass, &b_jetAK4Puppi_SV1mass);
   fChain->SetBranchAddress("jetAK4Puppi_jecUncertainty", jetAK4Puppi_jecUncertainty, &b_jetAK4Puppi_jecUncertainty);
   fChain->SetBranchAddress("jetAK4Puppi_PtResolution", jetAK4Puppi_PtResolution, &b_jetAK4Puppi_PtResolution);
   fChain->SetBranchAddress("jetAK4Puppi_JERSF", jetAK4Puppi_JERSF, &b_jetAK4Puppi_JERSF);
   fChain->SetBranchAddress("jetAK4Puppi_JERSFUp", jetAK4Puppi_JERSFUp, &b_jetAK4Puppi_JERSFUp);
   fChain->SetBranchAddress("jetAK4Puppi_JERSFDown", jetAK4Puppi_JERSFDown, &b_jetAK4Puppi_JERSFDown);
   fChain->SetBranchAddress("jetAK4Puppi_SmearedPt", jetAK4Puppi_SmearedPt, &b_jetAK4Puppi_SmearedPt);
   fChain->SetBranchAddress("jetAK4Puppi_SmearedPEta", jetAK4Puppi_SmearedPEta, &b_jetAK4Puppi_SmearedPEta);
   fChain->SetBranchAddress("jetAK4Puppi_SmearedPhi", jetAK4Puppi_SmearedPhi, &b_jetAK4Puppi_SmearedPhi);
   fChain->SetBranchAddress("jetAK4Puppi_SmearedE", jetAK4Puppi_SmearedE, &b_jetAK4Puppi_SmearedE);
   fChain->SetBranchAddress("jetAK4Puppi_Keys", &jetAK4Puppi_Keys, &b_jetAK4Puppi_Keys);
   fChain->SetBranchAddress("jetAK8_size", &jetAK8_size, &b_jetAK8_size);
   fChain->SetBranchAddress("jetAK8_Pt", jetAK8_Pt, &b_jetAK8_Pt);
   fChain->SetBranchAddress("jetAK8_Eta", jetAK8_Eta, &b_jetAK8_Eta);
   fChain->SetBranchAddress("jetAK8_Phi", jetAK8_Phi, &b_jetAK8_Phi);
   fChain->SetBranchAddress("jetAK8_E", jetAK8_E, &b_jetAK8_E);
   fChain->SetBranchAddress("jetAK8_Charge", jetAK8_Charge, &b_jetAK8_Charge);
   fChain->SetBranchAddress("jetAK8_CSVv2", jetAK8_CSVv2, &b_jetAK8_CSVv2);
   fChain->SetBranchAddress("jetAK8_DoubleBAK8", jetAK8_DoubleBAK8, &b_jetAK8_DoubleBAK8);
   fChain->SetBranchAddress("jetAK8_DoubleBCA15", jetAK8_DoubleBCA15, &b_jetAK8_DoubleBCA15);
   fChain->SetBranchAddress("jetAK8_CMVAv2", jetAK8_CMVAv2, &b_jetAK8_CMVAv2);
   fChain->SetBranchAddress("jetAK8_CvsL", jetAK8_CvsL, &b_jetAK8_CvsL);
   fChain->SetBranchAddress("jetAK8_CvsB", jetAK8_CvsB, &b_jetAK8_CvsB);
   fChain->SetBranchAddress("jetAK8_CMVA", jetAK8_CMVA, &b_jetAK8_CMVA);
   fChain->SetBranchAddress("jetAK8_GenPartonY", jetAK8_GenPartonY, &b_jetAK8_GenPartonY);
   fChain->SetBranchAddress("jetAK8_GenPartonEta", jetAK8_GenPartonEta, &b_jetAK8_GenPartonEta);
   fChain->SetBranchAddress("jetAK8_GenPartonPhi", jetAK8_GenPartonPhi, &b_jetAK8_GenPartonPhi);
   fChain->SetBranchAddress("jetAK8_GenPartonPt", jetAK8_GenPartonPt, &b_jetAK8_GenPartonPt);
   fChain->SetBranchAddress("jetAK8_GenPartonE", jetAK8_GenPartonE, &b_jetAK8_GenPartonE);
   fChain->SetBranchAddress("jetAK8_GenPartonCharge", jetAK8_GenPartonCharge, &b_jetAK8_GenPartonCharge);
   fChain->SetBranchAddress("jetAK8_PartonFlavour", jetAK8_PartonFlavour, &b_jetAK8_PartonFlavour);
   fChain->SetBranchAddress("jetAK8_HadronFlavour", jetAK8_HadronFlavour, &b_jetAK8_HadronFlavour);
   fChain->SetBranchAddress("jetAK8_GenJetY", jetAK8_GenJetY, &b_jetAK8_GenJetY);
   fChain->SetBranchAddress("jetAK8_GenJetEta", jetAK8_GenJetEta, &b_jetAK8_GenJetEta);
   fChain->SetBranchAddress("jetAK8_GenJetPhi", jetAK8_GenJetPhi, &b_jetAK8_GenJetPhi);
   fChain->SetBranchAddress("jetAK8_GenJetPt", jetAK8_GenJetPt, &b_jetAK8_GenJetPt);
   fChain->SetBranchAddress("jetAK8_GenJetE", jetAK8_GenJetE, &b_jetAK8_GenJetE);
   fChain->SetBranchAddress("jetAK8_GenJetCharge", jetAK8_GenJetCharge, &b_jetAK8_GenJetCharge);
   fChain->SetBranchAddress("jetAK8_muonMultiplicity", jetAK8_muonMultiplicity, &b_jetAK8_muonMultiplicity);
   fChain->SetBranchAddress("jetAK8_PhotonEnergy", jetAK8_PhotonEnergy, &b_jetAK8_PhotonEnergy);
   fChain->SetBranchAddress("jetAK8_ElectronEnergy", jetAK8_ElectronEnergy, &b_jetAK8_ElectronEnergy);
   fChain->SetBranchAddress("jetAK8_MuonEnergy", jetAK8_MuonEnergy, &b_jetAK8_MuonEnergy);
   fChain->SetBranchAddress("jetAK8_HFHadronEnergy", jetAK8_HFHadronEnergy, &b_jetAK8_HFHadronEnergy);
   fChain->SetBranchAddress("jetAK8_HFEMEnergy", jetAK8_HFEMEnergy, &b_jetAK8_HFEMEnergy);
   fChain->SetBranchAddress("jetAK8_ChargedHadronMultiplicity", jetAK8_ChargedHadronMultiplicity, &b_jetAK8_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("jetAK8_numberOfDaughters", jetAK8_numberOfDaughters, &b_jetAK8_numberOfDaughters);
   fChain->SetBranchAddress("jetAK8_neutralHadronMultiplicity", jetAK8_neutralHadronMultiplicity, &b_jetAK8_neutralHadronMultiplicity);
   fChain->SetBranchAddress("jetAK8_neutralHadronEnergy", jetAK8_neutralHadronEnergy, &b_jetAK8_neutralHadronEnergy);
   fChain->SetBranchAddress("jetAK8_neutralEmEnergy", jetAK8_neutralEmEnergy, &b_jetAK8_neutralEmEnergy);
   fChain->SetBranchAddress("jetAK8_chargedEmEnergy", jetAK8_chargedEmEnergy, &b_jetAK8_chargedEmEnergy);
   fChain->SetBranchAddress("jetAK8_chargedHadronEnergy", jetAK8_chargedHadronEnergy, &b_jetAK8_chargedHadronEnergy);
   fChain->SetBranchAddress("jetAK8_photonMultiplicity", jetAK8_photonMultiplicity, &b_jetAK8_photonMultiplicity);
   fChain->SetBranchAddress("jetAK8_electronMultiplicity", jetAK8_electronMultiplicity, &b_jetAK8_electronMultiplicity);
   fChain->SetBranchAddress("jetAK8_HFHadronMultiplicity", jetAK8_HFHadronMultiplicity, &b_jetAK8_HFHadronMultiplicity);
   fChain->SetBranchAddress("jetAK8_HFEMMultiplicity", jetAK8_HFEMMultiplicity, &b_jetAK8_HFEMMultiplicity);
   fChain->SetBranchAddress("jetAK8_ChargeMuEnergy", jetAK8_ChargeMuEnergy, &b_jetAK8_ChargeMuEnergy);
   fChain->SetBranchAddress("jetAK8_neutralMultiplicity", jetAK8_neutralMultiplicity, &b_jetAK8_neutralMultiplicity);
   fChain->SetBranchAddress("jetAK8_neutralHadronEnergyFrac", jetAK8_neutralHadronEnergyFrac, &b_jetAK8_neutralHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK8_neutralEmEnergyFrac", jetAK8_neutralEmEnergyFrac, &b_jetAK8_neutralEmEnergyFrac);
   fChain->SetBranchAddress("jetAK8_chargedHadronEnergyFrac", jetAK8_chargedHadronEnergyFrac, &b_jetAK8_chargedHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK8_muonEnergyFrac", jetAK8_muonEnergyFrac, &b_jetAK8_muonEnergyFrac);
   fChain->SetBranchAddress("jetAK8_chargedEmEnergyFrac", jetAK8_chargedEmEnergyFrac, &b_jetAK8_chargedEmEnergyFrac);
   fChain->SetBranchAddress("jetAK8_chargedMultiplicity", jetAK8_chargedMultiplicity, &b_jetAK8_chargedMultiplicity);
   fChain->SetBranchAddress("jetAK8_NumConstituents", jetAK8_NumConstituents, &b_jetAK8_NumConstituents);
   fChain->SetBranchAddress("jetAK8_jecFactor0", jetAK8_jecFactor0, &b_jetAK8_jecFactor0);
   fChain->SetBranchAddress("jetAK8_jecFactorL3Absolute", jetAK8_jecFactorL3Absolute, &b_jetAK8_jecFactorL3Absolute);
   fChain->SetBranchAddress("jetAK8_jetArea", jetAK8_jetArea, &b_jetAK8_jetArea);
   fChain->SetBranchAddress("jetAK8_nSV", jetAK8_nSV, &b_jetAK8_nSV);
   fChain->SetBranchAddress("jetAK8_SV0mass", jetAK8_SV0mass, &b_jetAK8_SV0mass);
   fChain->SetBranchAddress("jetAK8_SV1mass", jetAK8_SV1mass, &b_jetAK8_SV1mass);
   fChain->SetBranchAddress("jetAK8_jecUncertainty", jetAK8_jecUncertainty, &b_jetAK8_jecUncertainty);
   fChain->SetBranchAddress("jetAK8_PtResolution", jetAK8_PtResolution, &b_jetAK8_PtResolution);
   fChain->SetBranchAddress("jetAK8_JERSF", jetAK8_JERSF, &b_jetAK8_JERSF);
   fChain->SetBranchAddress("jetAK8_JERSFUp", jetAK8_JERSFUp, &b_jetAK8_JERSFUp);
   fChain->SetBranchAddress("jetAK8_JERSFDown", jetAK8_JERSFDown, &b_jetAK8_JERSFDown);
   fChain->SetBranchAddress("jetAK8_SmearedPt", jetAK8_SmearedPt, &b_jetAK8_SmearedPt);
   fChain->SetBranchAddress("jetAK8_SmearedPEta", jetAK8_SmearedPEta, &b_jetAK8_SmearedPEta);
   fChain->SetBranchAddress("jetAK8_SmearedPhi", jetAK8_SmearedPhi, &b_jetAK8_SmearedPhi);
   fChain->SetBranchAddress("jetAK8_SmearedE", jetAK8_SmearedE, &b_jetAK8_SmearedE);
   fChain->SetBranchAddress("jetAK8_vSubjetIndex0", jetAK8_vSubjetIndex0, &b_jetAK8_vSubjetIndex0);
   fChain->SetBranchAddress("jetAK8_vSubjetIndex1", jetAK8_vSubjetIndex1, &b_jetAK8_vSubjetIndex1);
   fChain->SetBranchAddress("jetAK8_tau1", jetAK8_tau1, &b_jetAK8_tau1);
   fChain->SetBranchAddress("jetAK8_tau2", jetAK8_tau2, &b_jetAK8_tau2);
   fChain->SetBranchAddress("jetAK8_tau3", jetAK8_tau3, &b_jetAK8_tau3);
   fChain->SetBranchAddress("jetAK8_softDropMass", jetAK8_softDropMass, &b_jetAK8_softDropMass);
   fChain->SetBranchAddress("jetAK8_trimmedMass", jetAK8_trimmedMass, &b_jetAK8_trimmedMass);
   fChain->SetBranchAddress("jetAK8_prunedMass", jetAK8_prunedMass, &b_jetAK8_prunedMass);
   fChain->SetBranchAddress("jetAK8_filteredMass", jetAK8_filteredMass, &b_jetAK8_filteredMass);
   fChain->SetBranchAddress("jetAK8_Keys", &jetAK8_Keys, &b_jetAK8_Keys);
   fChain->SetBranchAddress("jetAK8Puppi_size", &jetAK8Puppi_size, &b_jetAK8Puppi_size);
   fChain->SetBranchAddress("jetAK8Puppi_Pt", jetAK8Puppi_Pt, &b_jetAK8Puppi_Pt);
   fChain->SetBranchAddress("jetAK8Puppi_Eta", jetAK8Puppi_Eta, &b_jetAK8Puppi_Eta);
   fChain->SetBranchAddress("jetAK8Puppi_Phi", jetAK8Puppi_Phi, &b_jetAK8Puppi_Phi);
   fChain->SetBranchAddress("jetAK8Puppi_E", jetAK8Puppi_E, &b_jetAK8Puppi_E);
   fChain->SetBranchAddress("jetAK8Puppi_Charge", jetAK8Puppi_Charge, &b_jetAK8Puppi_Charge);
   fChain->SetBranchAddress("jetAK8Puppi_CSVv2", jetAK8Puppi_CSVv2, &b_jetAK8Puppi_CSVv2);
   fChain->SetBranchAddress("jetAK8Puppi_DoubleBAK8", jetAK8Puppi_DoubleBAK8, &b_jetAK8Puppi_DoubleBAK8);
   fChain->SetBranchAddress("jetAK8Puppi_DoubleBCA15", jetAK8Puppi_DoubleBCA15, &b_jetAK8Puppi_DoubleBCA15);
   fChain->SetBranchAddress("jetAK8Puppi_CMVAv2", jetAK8Puppi_CMVAv2, &b_jetAK8Puppi_CMVAv2);
   fChain->SetBranchAddress("jetAK8Puppi_CvsL", jetAK8Puppi_CvsL, &b_jetAK8Puppi_CvsL);
   fChain->SetBranchAddress("jetAK8Puppi_CvsB", jetAK8Puppi_CvsB, &b_jetAK8Puppi_CvsB);
   fChain->SetBranchAddress("jetAK8Puppi_CMVA", jetAK8Puppi_CMVA, &b_jetAK8Puppi_CMVA);
   fChain->SetBranchAddress("jetAK8Puppi_GenPartonY", jetAK8Puppi_GenPartonY, &b_jetAK8Puppi_GenPartonY);
   fChain->SetBranchAddress("jetAK8Puppi_GenPartonEta", jetAK8Puppi_GenPartonEta, &b_jetAK8Puppi_GenPartonEta);
   fChain->SetBranchAddress("jetAK8Puppi_GenPartonPhi", jetAK8Puppi_GenPartonPhi, &b_jetAK8Puppi_GenPartonPhi);
   fChain->SetBranchAddress("jetAK8Puppi_GenPartonPt", jetAK8Puppi_GenPartonPt, &b_jetAK8Puppi_GenPartonPt);
   fChain->SetBranchAddress("jetAK8Puppi_GenPartonE", jetAK8Puppi_GenPartonE, &b_jetAK8Puppi_GenPartonE);
   fChain->SetBranchAddress("jetAK8Puppi_GenPartonCharge", jetAK8Puppi_GenPartonCharge, &b_jetAK8Puppi_GenPartonCharge);
   fChain->SetBranchAddress("jetAK8Puppi_PartonFlavour", jetAK8Puppi_PartonFlavour, &b_jetAK8Puppi_PartonFlavour);
   fChain->SetBranchAddress("jetAK8Puppi_HadronFlavour", jetAK8Puppi_HadronFlavour, &b_jetAK8Puppi_HadronFlavour);
   fChain->SetBranchAddress("jetAK8Puppi_GenJetY", jetAK8Puppi_GenJetY, &b_jetAK8Puppi_GenJetY);
   fChain->SetBranchAddress("jetAK8Puppi_GenJetEta", jetAK8Puppi_GenJetEta, &b_jetAK8Puppi_GenJetEta);
   fChain->SetBranchAddress("jetAK8Puppi_GenJetPhi", jetAK8Puppi_GenJetPhi, &b_jetAK8Puppi_GenJetPhi);
   fChain->SetBranchAddress("jetAK8Puppi_GenJetPt", jetAK8Puppi_GenJetPt, &b_jetAK8Puppi_GenJetPt);
   fChain->SetBranchAddress("jetAK8Puppi_GenJetE", jetAK8Puppi_GenJetE, &b_jetAK8Puppi_GenJetE);
   fChain->SetBranchAddress("jetAK8Puppi_GenJetCharge", jetAK8Puppi_GenJetCharge, &b_jetAK8Puppi_GenJetCharge);
   fChain->SetBranchAddress("jetAK8Puppi_muonMultiplicity", jetAK8Puppi_muonMultiplicity, &b_jetAK8Puppi_muonMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_PhotonEnergy", jetAK8Puppi_PhotonEnergy, &b_jetAK8Puppi_PhotonEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_ElectronEnergy", jetAK8Puppi_ElectronEnergy, &b_jetAK8Puppi_ElectronEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_MuonEnergy", jetAK8Puppi_MuonEnergy, &b_jetAK8Puppi_MuonEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_HFHadronEnergy", jetAK8Puppi_HFHadronEnergy, &b_jetAK8Puppi_HFHadronEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_HFEMEnergy", jetAK8Puppi_HFEMEnergy, &b_jetAK8Puppi_HFEMEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_ChargedHadronMultiplicity", jetAK8Puppi_ChargedHadronMultiplicity, &b_jetAK8Puppi_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_numberOfDaughters", jetAK8Puppi_numberOfDaughters, &b_jetAK8Puppi_numberOfDaughters);
   fChain->SetBranchAddress("jetAK8Puppi_neutralHadronMultiplicity", jetAK8Puppi_neutralHadronMultiplicity, &b_jetAK8Puppi_neutralHadronMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_neutralHadronEnergy", jetAK8Puppi_neutralHadronEnergy, &b_jetAK8Puppi_neutralHadronEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_neutralEmEnergy", jetAK8Puppi_neutralEmEnergy, &b_jetAK8Puppi_neutralEmEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_chargedEmEnergy", jetAK8Puppi_chargedEmEnergy, &b_jetAK8Puppi_chargedEmEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_chargedHadronEnergy", jetAK8Puppi_chargedHadronEnergy, &b_jetAK8Puppi_chargedHadronEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_photonMultiplicity", jetAK8Puppi_photonMultiplicity, &b_jetAK8Puppi_photonMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_electronMultiplicity", jetAK8Puppi_electronMultiplicity, &b_jetAK8Puppi_electronMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_HFHadronMultiplicity", jetAK8Puppi_HFHadronMultiplicity, &b_jetAK8Puppi_HFHadronMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_HFEMMultiplicity", jetAK8Puppi_HFEMMultiplicity, &b_jetAK8Puppi_HFEMMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_ChargeMuEnergy", jetAK8Puppi_ChargeMuEnergy, &b_jetAK8Puppi_ChargeMuEnergy);
   fChain->SetBranchAddress("jetAK8Puppi_neutralMultiplicity", jetAK8Puppi_neutralMultiplicity, &b_jetAK8Puppi_neutralMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_neutralHadronEnergyFrac", jetAK8Puppi_neutralHadronEnergyFrac, &b_jetAK8Puppi_neutralHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK8Puppi_neutralEmEnergyFrac", jetAK8Puppi_neutralEmEnergyFrac, &b_jetAK8Puppi_neutralEmEnergyFrac);
   fChain->SetBranchAddress("jetAK8Puppi_chargedHadronEnergyFrac", jetAK8Puppi_chargedHadronEnergyFrac, &b_jetAK8Puppi_chargedHadronEnergyFrac);
   fChain->SetBranchAddress("jetAK8Puppi_muonEnergyFrac", jetAK8Puppi_muonEnergyFrac, &b_jetAK8Puppi_muonEnergyFrac);
   fChain->SetBranchAddress("jetAK8Puppi_chargedEmEnergyFrac", jetAK8Puppi_chargedEmEnergyFrac, &b_jetAK8Puppi_chargedEmEnergyFrac);
   fChain->SetBranchAddress("jetAK8Puppi_chargedMultiplicity", jetAK8Puppi_chargedMultiplicity, &b_jetAK8Puppi_chargedMultiplicity);
   fChain->SetBranchAddress("jetAK8Puppi_NumConstituents", jetAK8Puppi_NumConstituents, &b_jetAK8Puppi_NumConstituents);
   fChain->SetBranchAddress("jetAK8Puppi_jecFactor0", jetAK8Puppi_jecFactor0, &b_jetAK8Puppi_jecFactor0);
   fChain->SetBranchAddress("jetAK8Puppi_jecFactorL3Absolute", jetAK8Puppi_jecFactorL3Absolute, &b_jetAK8Puppi_jecFactorL3Absolute);
   fChain->SetBranchAddress("jetAK8Puppi_jetArea", jetAK8Puppi_jetArea, &b_jetAK8Puppi_jetArea);
   fChain->SetBranchAddress("jetAK8Puppi_nSV", jetAK8Puppi_nSV, &b_jetAK8Puppi_nSV);
   fChain->SetBranchAddress("jetAK8Puppi_SV0mass", jetAK8Puppi_SV0mass, &b_jetAK8Puppi_SV0mass);
   fChain->SetBranchAddress("jetAK8Puppi_SV1mass", jetAK8Puppi_SV1mass, &b_jetAK8Puppi_SV1mass);
   fChain->SetBranchAddress("jetAK8Puppi_jecUncertainty", jetAK8Puppi_jecUncertainty, &b_jetAK8Puppi_jecUncertainty);
   fChain->SetBranchAddress("jetAK8Puppi_PtResolution", jetAK8Puppi_PtResolution, &b_jetAK8Puppi_PtResolution);
   fChain->SetBranchAddress("jetAK8Puppi_JERSF", jetAK8Puppi_JERSF, &b_jetAK8Puppi_JERSF);
   fChain->SetBranchAddress("jetAK8Puppi_JERSFUp", jetAK8Puppi_JERSFUp, &b_jetAK8Puppi_JERSFUp);
   fChain->SetBranchAddress("jetAK8Puppi_JERSFDown", jetAK8Puppi_JERSFDown, &b_jetAK8Puppi_JERSFDown);
   fChain->SetBranchAddress("jetAK8Puppi_SmearedPt", jetAK8Puppi_SmearedPt, &b_jetAK8Puppi_SmearedPt);
   fChain->SetBranchAddress("jetAK8Puppi_SmearedPEta", jetAK8Puppi_SmearedPEta, &b_jetAK8Puppi_SmearedPEta);
   fChain->SetBranchAddress("jetAK8Puppi_SmearedPhi", jetAK8Puppi_SmearedPhi, &b_jetAK8Puppi_SmearedPhi);
   fChain->SetBranchAddress("jetAK8Puppi_SmearedE", jetAK8Puppi_SmearedE, &b_jetAK8Puppi_SmearedE);
   fChain->SetBranchAddress("jetAK8Puppi_vSubjetIndex0", jetAK8Puppi_vSubjetIndex0, &b_jetAK8Puppi_vSubjetIndex0);
   fChain->SetBranchAddress("jetAK8Puppi_vSubjetIndex1", jetAK8Puppi_vSubjetIndex1, &b_jetAK8Puppi_vSubjetIndex1);
   fChain->SetBranchAddress("jetAK8Puppi_tau1", jetAK8Puppi_tau1, &b_jetAK8Puppi_tau1);
   fChain->SetBranchAddress("jetAK8Puppi_tau2", jetAK8Puppi_tau2, &b_jetAK8Puppi_tau2);
   fChain->SetBranchAddress("jetAK8Puppi_tau3", jetAK8Puppi_tau3, &b_jetAK8Puppi_tau3);
   fChain->SetBranchAddress("jetAK8Puppi_softDropMass", jetAK8Puppi_softDropMass, &b_jetAK8Puppi_softDropMass);
   fChain->SetBranchAddress("jetAK8Puppi_trimmedMass", jetAK8Puppi_trimmedMass, &b_jetAK8Puppi_trimmedMass);
   fChain->SetBranchAddress("jetAK8Puppi_prunedMass", jetAK8Puppi_prunedMass, &b_jetAK8Puppi_prunedMass);
   fChain->SetBranchAddress("jetAK8Puppi_filteredMass", jetAK8Puppi_filteredMass, &b_jetAK8Puppi_filteredMass);
   fChain->SetBranchAddress("jetAK8Puppi_Keys", &jetAK8Puppi_Keys, &b_jetAK8Puppi_Keys);
   fChain->SetBranchAddress("subjetAK8_size", &subjetAK8_size, &b_subjetAK8_size);
   fChain->SetBranchAddress("subjetAK8_Pt", subjetAK8_Pt, &b_subjetAK8_Pt);
   fChain->SetBranchAddress("subjetAK8_Eta", subjetAK8_Eta, &b_subjetAK8_Eta);
   fChain->SetBranchAddress("subjetAK8_Phi", subjetAK8_Phi, &b_subjetAK8_Phi);
   fChain->SetBranchAddress("subjetAK8_E", subjetAK8_E, &b_subjetAK8_E);
   fChain->SetBranchAddress("subjetAK8_Charge", subjetAK8_Charge, &b_subjetAK8_Charge);
   fChain->SetBranchAddress("subjetAK8_CSVv2", subjetAK8_CSVv2, &b_subjetAK8_CSVv2);
   fChain->SetBranchAddress("subjetAK8_DoubleBAK8", subjetAK8_DoubleBAK8, &b_subjetAK8_DoubleBAK8);
   fChain->SetBranchAddress("subjetAK8_DoubleBCA15", subjetAK8_DoubleBCA15, &b_subjetAK8_DoubleBCA15);
   fChain->SetBranchAddress("subjetAK8_CMVAv2", subjetAK8_CMVAv2, &b_subjetAK8_CMVAv2);
   fChain->SetBranchAddress("subjetAK8_CvsL", subjetAK8_CvsL, &b_subjetAK8_CvsL);
   fChain->SetBranchAddress("subjetAK8_CvsB", subjetAK8_CvsB, &b_subjetAK8_CvsB);
   fChain->SetBranchAddress("subjetAK8_CMVA", subjetAK8_CMVA, &b_subjetAK8_CMVA);
   fChain->SetBranchAddress("subjetAK8_GenPartonY", subjetAK8_GenPartonY, &b_subjetAK8_GenPartonY);
   fChain->SetBranchAddress("subjetAK8_GenPartonEta", subjetAK8_GenPartonEta, &b_subjetAK8_GenPartonEta);
   fChain->SetBranchAddress("subjetAK8_GenPartonPhi", subjetAK8_GenPartonPhi, &b_subjetAK8_GenPartonPhi);
   fChain->SetBranchAddress("subjetAK8_GenPartonPt", subjetAK8_GenPartonPt, &b_subjetAK8_GenPartonPt);
   fChain->SetBranchAddress("subjetAK8_GenPartonE", subjetAK8_GenPartonE, &b_subjetAK8_GenPartonE);
   fChain->SetBranchAddress("subjetAK8_GenPartonCharge", subjetAK8_GenPartonCharge, &b_subjetAK8_GenPartonCharge);
   fChain->SetBranchAddress("subjetAK8_PartonFlavour", subjetAK8_PartonFlavour, &b_subjetAK8_PartonFlavour);
   fChain->SetBranchAddress("subjetAK8_HadronFlavour", subjetAK8_HadronFlavour, &b_subjetAK8_HadronFlavour);
   fChain->SetBranchAddress("subjetAK8_GenJetY", subjetAK8_GenJetY, &b_subjetAK8_GenJetY);
   fChain->SetBranchAddress("subjetAK8_GenJetEta", subjetAK8_GenJetEta, &b_subjetAK8_GenJetEta);
   fChain->SetBranchAddress("subjetAK8_GenJetPhi", subjetAK8_GenJetPhi, &b_subjetAK8_GenJetPhi);
   fChain->SetBranchAddress("subjetAK8_GenJetPt", subjetAK8_GenJetPt, &b_subjetAK8_GenJetPt);
   fChain->SetBranchAddress("subjetAK8_GenJetE", subjetAK8_GenJetE, &b_subjetAK8_GenJetE);
   fChain->SetBranchAddress("subjetAK8_GenJetCharge", subjetAK8_GenJetCharge, &b_subjetAK8_GenJetCharge);
   fChain->SetBranchAddress("subjetAK8_muonMultiplicity", subjetAK8_muonMultiplicity, &b_subjetAK8_muonMultiplicity);
   fChain->SetBranchAddress("subjetAK8_PhotonEnergy", subjetAK8_PhotonEnergy, &b_subjetAK8_PhotonEnergy);
   fChain->SetBranchAddress("subjetAK8_ElectronEnergy", subjetAK8_ElectronEnergy, &b_subjetAK8_ElectronEnergy);
   fChain->SetBranchAddress("subjetAK8_MuonEnergy", subjetAK8_MuonEnergy, &b_subjetAK8_MuonEnergy);
   fChain->SetBranchAddress("subjetAK8_HFHadronEnergy", subjetAK8_HFHadronEnergy, &b_subjetAK8_HFHadronEnergy);
   fChain->SetBranchAddress("subjetAK8_HFEMEnergy", subjetAK8_HFEMEnergy, &b_subjetAK8_HFEMEnergy);
   fChain->SetBranchAddress("subjetAK8_ChargedHadronMultiplicity", subjetAK8_ChargedHadronMultiplicity, &b_subjetAK8_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("subjetAK8_numberOfDaughters", subjetAK8_numberOfDaughters, &b_subjetAK8_numberOfDaughters);
   fChain->SetBranchAddress("subjetAK8_neutralHadronMultiplicity", subjetAK8_neutralHadronMultiplicity, &b_subjetAK8_neutralHadronMultiplicity);
   fChain->SetBranchAddress("subjetAK8_neutralHadronEnergy", subjetAK8_neutralHadronEnergy, &b_subjetAK8_neutralHadronEnergy);
   fChain->SetBranchAddress("subjetAK8_neutralEmEnergy", subjetAK8_neutralEmEnergy, &b_subjetAK8_neutralEmEnergy);
   fChain->SetBranchAddress("subjetAK8_chargedEmEnergy", subjetAK8_chargedEmEnergy, &b_subjetAK8_chargedEmEnergy);
   fChain->SetBranchAddress("subjetAK8_chargedHadronEnergy", subjetAK8_chargedHadronEnergy, &b_subjetAK8_chargedHadronEnergy);
   fChain->SetBranchAddress("subjetAK8_photonMultiplicity", subjetAK8_photonMultiplicity, &b_subjetAK8_photonMultiplicity);
   fChain->SetBranchAddress("subjetAK8_electronMultiplicity", subjetAK8_electronMultiplicity, &b_subjetAK8_electronMultiplicity);
   fChain->SetBranchAddress("subjetAK8_HFHadronMultiplicity", subjetAK8_HFHadronMultiplicity, &b_subjetAK8_HFHadronMultiplicity);
   fChain->SetBranchAddress("subjetAK8_HFEMMultiplicity", subjetAK8_HFEMMultiplicity, &b_subjetAK8_HFEMMultiplicity);
   fChain->SetBranchAddress("subjetAK8_ChargeMuEnergy", subjetAK8_ChargeMuEnergy, &b_subjetAK8_ChargeMuEnergy);
   fChain->SetBranchAddress("subjetAK8_neutralMultiplicity", subjetAK8_neutralMultiplicity, &b_subjetAK8_neutralMultiplicity);
   fChain->SetBranchAddress("subjetAK8_neutralHadronEnergyFrac", subjetAK8_neutralHadronEnergyFrac, &b_subjetAK8_neutralHadronEnergyFrac);
   fChain->SetBranchAddress("subjetAK8_neutralEmEnergyFrac", subjetAK8_neutralEmEnergyFrac, &b_subjetAK8_neutralEmEnergyFrac);
   fChain->SetBranchAddress("subjetAK8_chargedHadronEnergyFrac", subjetAK8_chargedHadronEnergyFrac, &b_subjetAK8_chargedHadronEnergyFrac);
   fChain->SetBranchAddress("subjetAK8_muonEnergyFrac", subjetAK8_muonEnergyFrac, &b_subjetAK8_muonEnergyFrac);
   fChain->SetBranchAddress("subjetAK8_chargedEmEnergyFrac", subjetAK8_chargedEmEnergyFrac, &b_subjetAK8_chargedEmEnergyFrac);
   fChain->SetBranchAddress("subjetAK8_chargedMultiplicity", subjetAK8_chargedMultiplicity, &b_subjetAK8_chargedMultiplicity);
   fChain->SetBranchAddress("subjetAK8_NumConstituents", subjetAK8_NumConstituents, &b_subjetAK8_NumConstituents);
   fChain->SetBranchAddress("subjetAK8_jecFactor0", subjetAK8_jecFactor0, &b_subjetAK8_jecFactor0);
   fChain->SetBranchAddress("subjetAK8_jecFactorL3Absolute", subjetAK8_jecFactorL3Absolute, &b_subjetAK8_jecFactorL3Absolute);
   fChain->SetBranchAddress("subjetAK8_jetArea", subjetAK8_jetArea, &b_subjetAK8_jetArea);
   fChain->SetBranchAddress("subjetAK8_nSV", subjetAK8_nSV, &b_subjetAK8_nSV);
   fChain->SetBranchAddress("subjetAK8_SV0mass", subjetAK8_SV0mass, &b_subjetAK8_SV0mass);
   fChain->SetBranchAddress("subjetAK8_SV1mass", subjetAK8_SV1mass, &b_subjetAK8_SV1mass);
   fChain->SetBranchAddress("subjetAK8_Keys", &subjetAK8_Keys, &b_subjetAK8_Keys);
   fChain->SetBranchAddress("subjetAK8Puppi_size", &subjetAK8Puppi_size, &b_subjetAK8Puppi_size);
   fChain->SetBranchAddress("subjetAK8Puppi_Pt", subjetAK8Puppi_Pt, &b_subjetAK8Puppi_Pt);
   fChain->SetBranchAddress("subjetAK8Puppi_Eta", subjetAK8Puppi_Eta, &b_subjetAK8Puppi_Eta);
   fChain->SetBranchAddress("subjetAK8Puppi_Phi", subjetAK8Puppi_Phi, &b_subjetAK8Puppi_Phi);
   fChain->SetBranchAddress("subjetAK8Puppi_E", subjetAK8Puppi_E, &b_subjetAK8Puppi_E);
   fChain->SetBranchAddress("subjetAK8Puppi_Charge", subjetAK8Puppi_Charge, &b_subjetAK8Puppi_Charge);
   fChain->SetBranchAddress("subjetAK8Puppi_CSVv2", subjetAK8Puppi_CSVv2, &b_subjetAK8Puppi_CSVv2);
   fChain->SetBranchAddress("subjetAK8Puppi_DoubleBAK8", subjetAK8Puppi_DoubleBAK8, &b_subjetAK8Puppi_DoubleBAK8);
   fChain->SetBranchAddress("subjetAK8Puppi_DoubleBCA15", subjetAK8Puppi_DoubleBCA15, &b_subjetAK8Puppi_DoubleBCA15);
   fChain->SetBranchAddress("subjetAK8Puppi_CMVAv2", subjetAK8Puppi_CMVAv2, &b_subjetAK8Puppi_CMVAv2);
   fChain->SetBranchAddress("subjetAK8Puppi_CvsL", subjetAK8Puppi_CvsL, &b_subjetAK8Puppi_CvsL);
   fChain->SetBranchAddress("subjetAK8Puppi_CvsB", subjetAK8Puppi_CvsB, &b_subjetAK8Puppi_CvsB);
   fChain->SetBranchAddress("subjetAK8Puppi_CMVA", subjetAK8Puppi_CMVA, &b_subjetAK8Puppi_CMVA);
   fChain->SetBranchAddress("subjetAK8Puppi_GenPartonY", subjetAK8Puppi_GenPartonY, &b_subjetAK8Puppi_GenPartonY);
   fChain->SetBranchAddress("subjetAK8Puppi_GenPartonEta", subjetAK8Puppi_GenPartonEta, &b_subjetAK8Puppi_GenPartonEta);
   fChain->SetBranchAddress("subjetAK8Puppi_GenPartonPhi", subjetAK8Puppi_GenPartonPhi, &b_subjetAK8Puppi_GenPartonPhi);
   fChain->SetBranchAddress("subjetAK8Puppi_GenPartonPt", subjetAK8Puppi_GenPartonPt, &b_subjetAK8Puppi_GenPartonPt);
   fChain->SetBranchAddress("subjetAK8Puppi_GenPartonE", subjetAK8Puppi_GenPartonE, &b_subjetAK8Puppi_GenPartonE);
   fChain->SetBranchAddress("subjetAK8Puppi_GenPartonCharge", subjetAK8Puppi_GenPartonCharge, &b_subjetAK8Puppi_GenPartonCharge);
   fChain->SetBranchAddress("subjetAK8Puppi_PartonFlavour", subjetAK8Puppi_PartonFlavour, &b_subjetAK8Puppi_PartonFlavour);
   fChain->SetBranchAddress("subjetAK8Puppi_HadronFlavour", subjetAK8Puppi_HadronFlavour, &b_subjetAK8Puppi_HadronFlavour);
   fChain->SetBranchAddress("subjetAK8Puppi_GenJetY", subjetAK8Puppi_GenJetY, &b_subjetAK8Puppi_GenJetY);
   fChain->SetBranchAddress("subjetAK8Puppi_GenJetEta", subjetAK8Puppi_GenJetEta, &b_subjetAK8Puppi_GenJetEta);
   fChain->SetBranchAddress("subjetAK8Puppi_GenJetPhi", subjetAK8Puppi_GenJetPhi, &b_subjetAK8Puppi_GenJetPhi);
   fChain->SetBranchAddress("subjetAK8Puppi_GenJetPt", subjetAK8Puppi_GenJetPt, &b_subjetAK8Puppi_GenJetPt);
   fChain->SetBranchAddress("subjetAK8Puppi_GenJetE", subjetAK8Puppi_GenJetE, &b_subjetAK8Puppi_GenJetE);
   fChain->SetBranchAddress("subjetAK8Puppi_GenJetCharge", subjetAK8Puppi_GenJetCharge, &b_subjetAK8Puppi_GenJetCharge);
   fChain->SetBranchAddress("subjetAK8Puppi_muonMultiplicity", subjetAK8Puppi_muonMultiplicity, &b_subjetAK8Puppi_muonMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_PhotonEnergy", subjetAK8Puppi_PhotonEnergy, &b_subjetAK8Puppi_PhotonEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_ElectronEnergy", subjetAK8Puppi_ElectronEnergy, &b_subjetAK8Puppi_ElectronEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_MuonEnergy", subjetAK8Puppi_MuonEnergy, &b_subjetAK8Puppi_MuonEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_HFHadronEnergy", subjetAK8Puppi_HFHadronEnergy, &b_subjetAK8Puppi_HFHadronEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_HFEMEnergy", subjetAK8Puppi_HFEMEnergy, &b_subjetAK8Puppi_HFEMEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_ChargedHadronMultiplicity", subjetAK8Puppi_ChargedHadronMultiplicity, &b_subjetAK8Puppi_ChargedHadronMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_numberOfDaughters", subjetAK8Puppi_numberOfDaughters, &b_subjetAK8Puppi_numberOfDaughters);
   fChain->SetBranchAddress("subjetAK8Puppi_neutralHadronMultiplicity", subjetAK8Puppi_neutralHadronMultiplicity, &b_subjetAK8Puppi_neutralHadronMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_neutralHadronEnergy", subjetAK8Puppi_neutralHadronEnergy, &b_subjetAK8Puppi_neutralHadronEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_neutralEmEnergy", subjetAK8Puppi_neutralEmEnergy, &b_subjetAK8Puppi_neutralEmEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_chargedEmEnergy", subjetAK8Puppi_chargedEmEnergy, &b_subjetAK8Puppi_chargedEmEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_chargedHadronEnergy", subjetAK8Puppi_chargedHadronEnergy, &b_subjetAK8Puppi_chargedHadronEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_photonMultiplicity", subjetAK8Puppi_photonMultiplicity, &b_subjetAK8Puppi_photonMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_electronMultiplicity", subjetAK8Puppi_electronMultiplicity, &b_subjetAK8Puppi_electronMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_HFHadronMultiplicity", subjetAK8Puppi_HFHadronMultiplicity, &b_subjetAK8Puppi_HFHadronMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_HFEMMultiplicity", subjetAK8Puppi_HFEMMultiplicity, &b_subjetAK8Puppi_HFEMMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_ChargeMuEnergy", subjetAK8Puppi_ChargeMuEnergy, &b_subjetAK8Puppi_ChargeMuEnergy);
   fChain->SetBranchAddress("subjetAK8Puppi_neutralMultiplicity", subjetAK8Puppi_neutralMultiplicity, &b_subjetAK8Puppi_neutralMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_neutralHadronEnergyFrac", subjetAK8Puppi_neutralHadronEnergyFrac, &b_subjetAK8Puppi_neutralHadronEnergyFrac);
   fChain->SetBranchAddress("subjetAK8Puppi_neutralEmEnergyFrac", subjetAK8Puppi_neutralEmEnergyFrac, &b_subjetAK8Puppi_neutralEmEnergyFrac);
   fChain->SetBranchAddress("subjetAK8Puppi_chargedHadronEnergyFrac", subjetAK8Puppi_chargedHadronEnergyFrac, &b_subjetAK8Puppi_chargedHadronEnergyFrac);
   fChain->SetBranchAddress("subjetAK8Puppi_muonEnergyFrac", subjetAK8Puppi_muonEnergyFrac, &b_subjetAK8Puppi_muonEnergyFrac);
   fChain->SetBranchAddress("subjetAK8Puppi_chargedEmEnergyFrac", subjetAK8Puppi_chargedEmEnergyFrac, &b_subjetAK8Puppi_chargedEmEnergyFrac);
   fChain->SetBranchAddress("subjetAK8Puppi_chargedMultiplicity", subjetAK8Puppi_chargedMultiplicity, &b_subjetAK8Puppi_chargedMultiplicity);
   fChain->SetBranchAddress("subjetAK8Puppi_NumConstituents", subjetAK8Puppi_NumConstituents, &b_subjetAK8Puppi_NumConstituents);
   fChain->SetBranchAddress("subjetAK8Puppi_jecFactor0", subjetAK8Puppi_jecFactor0, &b_subjetAK8Puppi_jecFactor0);
   fChain->SetBranchAddress("subjetAK8Puppi_jecFactorL3Absolute", subjetAK8Puppi_jecFactorL3Absolute, &b_subjetAK8Puppi_jecFactorL3Absolute);
   fChain->SetBranchAddress("subjetAK8Puppi_jetArea", subjetAK8Puppi_jetArea, &b_subjetAK8Puppi_jetArea);
   fChain->SetBranchAddress("subjetAK8Puppi_nSV", subjetAK8Puppi_nSV, &b_subjetAK8Puppi_nSV);
   fChain->SetBranchAddress("subjetAK8Puppi_SV0mass", subjetAK8Puppi_SV0mass, &b_subjetAK8Puppi_SV0mass);
   fChain->SetBranchAddress("subjetAK8Puppi_SV1mass", subjetAK8Puppi_SV1mass, &b_subjetAK8Puppi_SV1mass);
   fChain->SetBranchAddress("subjetAK8Puppi_Keys", &subjetAK8Puppi_Keys, &b_subjetAK8Puppi_Keys);
   fChain->SetBranchAddress("genjetAK8SD_size", &genjetAK8SD_size, &b_genjetAK8SD_size);
   fChain->SetBranchAddress("genjetAK8SD_Pt", genjetAK8SD_Pt, &b_genjetAK8SD_Pt);
   fChain->SetBranchAddress("genjetAK8SD_Eta", genjetAK8SD_Eta, &b_genjetAK8SD_Eta);
   fChain->SetBranchAddress("genjetAK8SD_Phi", genjetAK8SD_Phi, &b_genjetAK8SD_Phi);
   fChain->SetBranchAddress("genjetAK8SD_E", genjetAK8SD_E, &b_genjetAK8SD_E);
   fChain->SetBranchAddress("genjetAK8SD_Charge", genjetAK8SD_Charge, &b_genjetAK8SD_Charge);
   fChain->SetBranchAddress("Flag_HBHEIsoNoiseFilterResult", &Flag_HBHEIsoNoiseFilterResult, &b_Flag_HBHEIsoNoiseFilterResult);
   fChain->SetBranchAddress("Flag_HBHENoiseFilterResult", &Flag_HBHENoiseFilterResult, &b_Flag_HBHENoiseFilterResult);
   fChain->SetBranchAddress("Flag_HBHENoiseFilterResultRun1", &Flag_HBHENoiseFilterResultRun1, &b_Flag_HBHENoiseFilterResultRun1);
   fChain->SetBranchAddress("Flag_HBHENoiseFilterResultRun2Loose", &Flag_HBHENoiseFilterResultRun2Loose, &b_Flag_HBHENoiseFilterResultRun2Loose);
   fChain->SetBranchAddress("Flag_HBHENoiseFilterResultRun2Tight", &Flag_HBHENoiseFilterResultRun2Tight, &b_Flag_HBHENoiseFilterResultRun2Tight);
   fChain->SetBranchAddress("Flag_goodVertices", &Flag_goodVertices, &b_Flag_goodVertices);
   fChain->SetBranchAddress("Flag_eeBadScFilter", &Flag_eeBadScFilter, &b_Flag_eeBadScFilter);
   fChain->SetBranchAddress("Flag_ecalLaserCorrFilter", &Flag_ecalLaserCorrFilter, &b_Flag_ecalLaserCorrFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellTriggerPrimitiveFilter", &Flag_EcalDeadCellTriggerPrimitiveFilter, &b_Flag_EcalDeadCellTriggerPrimitiveFilter);
   fChain->SetBranchAddress("Flag_EcalDeadCellBoundaryEnergyFilter", &Flag_EcalDeadCellBoundaryEnergyFilter, &b_Flag_EcalDeadCellBoundaryEnergyFilter);
   fChain->SetBranchAddress("Flag_HcalStripHaloFilter", &Flag_HcalStripHaloFilter, &b_Flag_HcalStripHaloFilter);
   fChain->SetBranchAddress("Flag_hcalLaserEventFilter", &Flag_hcalLaserEventFilter, &b_Flag_hcalLaserEventFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseFilter", &Flag_HBHENoiseFilter, &b_Flag_HBHENoiseFilter);
   fChain->SetBranchAddress("Flag_HBHENoiseIsoFilter", &Flag_HBHENoiseIsoFilter, &b_Flag_HBHENoiseIsoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloFilter", &Flag_CSCTightHaloFilter, &b_Flag_CSCTightHaloFilter);
   fChain->SetBranchAddress("Flag_CSCTightHaloTrkMuUnvetoFilter", &Flag_CSCTightHaloTrkMuUnvetoFilter, &b_Flag_CSCTightHaloTrkMuUnvetoFilter);
   fChain->SetBranchAddress("Flag_CSCTightHalo2015Filter", &Flag_CSCTightHalo2015Filter, &b_Flag_CSCTightHalo2015Filter);
   fChain->SetBranchAddress("Flag_muonBadTrackFilter", &Flag_muonBadTrackFilter, &b_Flag_muonBadTrackFilter);
   fChain->SetBranchAddress("Flag_chargedHadronTrackResolutionFilter", &Flag_chargedHadronTrackResolutionFilter, &b_Flag_chargedHadronTrackResolutionFilter);
   fChain->SetBranchAddress("Flag_trkPOGFilters", &Flag_trkPOGFilters, &b_Flag_trkPOGFilters);
   fChain->SetBranchAddress("Flag_trkPOG_manystripclus53X", &Flag_trkPOG_manystripclus53X, &b_Flag_trkPOG_manystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_toomanystripclus53X", &Flag_trkPOG_toomanystripclus53X, &b_Flag_trkPOG_toomanystripclus53X);
   fChain->SetBranchAddress("Flag_trkPOG_logErrorTooManyClusters", &Flag_trkPOG_logErrorTooManyClusters, &b_Flag_trkPOG_logErrorTooManyClusters);
   fChain->SetBranchAddress("Flag_METFiltersHLT_AK8PFJet360_TrimMass30", &Flag_METFiltersHLT_AK8PFJet360_TrimMass30, &b_Flag_METFiltersHLT_AK8PFJet360_TrimMass30);
   fChain->SetBranchAddress("HLT_PFJet40", &HLT_PFJet40, &b_HLT_PFJet40);
   fChain->SetBranchAddress("HLT_PFJet40_prescale", &HLT_PFJet40_prescale, &b_HLT_PFJet40_prescale);
   fChain->SetBranchAddress("HLT_PFJet60", &HLT_PFJet60, &b_HLT_PFJet60);
   fChain->SetBranchAddress("HLT_PFJet60_prescale", &HLT_PFJet60_prescale, &b_HLT_PFJet60_prescale);
   fChain->SetBranchAddress("HLT_PFJet80", &HLT_PFJet80, &b_HLT_PFJet80);
   fChain->SetBranchAddress("HLT_PFJet80_prescale", &HLT_PFJet80_prescale, &b_HLT_PFJet80_prescale);
   fChain->SetBranchAddress("HLT_PFJet140", &HLT_PFJet140, &b_HLT_PFJet140);
   fChain->SetBranchAddress("HLT_PFJet140_prescale", &HLT_PFJet140_prescale, &b_HLT_PFJet140_prescale);
   fChain->SetBranchAddress("HLT_PFJet200", &HLT_PFJet200, &b_HLT_PFJet200);
   fChain->SetBranchAddress("HLT_PFJet200_prescale", &HLT_PFJet200_prescale, &b_HLT_PFJet200_prescale);
   fChain->SetBranchAddress("HLT_PFJet260", &HLT_PFJet260, &b_HLT_PFJet260);
   fChain->SetBranchAddress("HLT_PFJet260_prescale", &HLT_PFJet260_prescale, &b_HLT_PFJet260_prescale);
   fChain->SetBranchAddress("HLT_PFJet320", &HLT_PFJet320, &b_HLT_PFJet320);
   fChain->SetBranchAddress("HLT_PFJet320_prescale", &HLT_PFJet320_prescale, &b_HLT_PFJet320_prescale);
   fChain->SetBranchAddress("HLT_PFJet400", &HLT_PFJet400, &b_HLT_PFJet400);
   fChain->SetBranchAddress("HLT_PFJet400_prescale", &HLT_PFJet400_prescale, &b_HLT_PFJet400_prescale);
   fChain->SetBranchAddress("HLT_PFJet450", &HLT_PFJet450, &b_HLT_PFJet450);
   fChain->SetBranchAddress("HLT_PFJet450_prescale", &HLT_PFJet450_prescale, &b_HLT_PFJet450_prescale);
   fChain->SetBranchAddress("HLT_PFJet500", &HLT_PFJet500, &b_HLT_PFJet500);
   fChain->SetBranchAddress("HLT_PFJet500_prescale", &HLT_PFJet500_prescale, &b_HLT_PFJet500_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve40", &HLT_DiPFJetAve40, &b_HLT_DiPFJetAve40);
   fChain->SetBranchAddress("HLT_DiPFJetAve40_prescale", &HLT_DiPFJetAve40_prescale, &b_HLT_DiPFJetAve40_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve60", &HLT_DiPFJetAve60, &b_HLT_DiPFJetAve60);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_prescale", &HLT_DiPFJetAve60_prescale, &b_HLT_DiPFJetAve60_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve80", &HLT_DiPFJetAve80, &b_HLT_DiPFJetAve80);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_prescale", &HLT_DiPFJetAve80_prescale, &b_HLT_DiPFJetAve80_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve140", &HLT_DiPFJetAve140, &b_HLT_DiPFJetAve140);
   fChain->SetBranchAddress("HLT_DiPFJetAve140_prescale", &HLT_DiPFJetAve140_prescale, &b_HLT_DiPFJetAve140_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve200", &HLT_DiPFJetAve200, &b_HLT_DiPFJetAve200);
   fChain->SetBranchAddress("HLT_DiPFJetAve200_prescale", &HLT_DiPFJetAve200_prescale, &b_HLT_DiPFJetAve200_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve260", &HLT_DiPFJetAve260, &b_HLT_DiPFJetAve260);
   fChain->SetBranchAddress("HLT_DiPFJetAve260_prescale", &HLT_DiPFJetAve260_prescale, &b_HLT_DiPFJetAve260_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve320", &HLT_DiPFJetAve320, &b_HLT_DiPFJetAve320);
   fChain->SetBranchAddress("HLT_DiPFJetAve320_prescale", &HLT_DiPFJetAve320_prescale, &b_HLT_DiPFJetAve320_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve400", &HLT_DiPFJetAve400, &b_HLT_DiPFJetAve400);
   fChain->SetBranchAddress("HLT_DiPFJetAve400_prescale", &HLT_DiPFJetAve400_prescale, &b_HLT_DiPFJetAve400_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve500", &HLT_DiPFJetAve500, &b_HLT_DiPFJetAve500);
   fChain->SetBranchAddress("HLT_DiPFJetAve500_prescale", &HLT_DiPFJetAve500_prescale, &b_HLT_DiPFJetAve500_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC", &HLT_DiPFJetAve60_HFJEC, &b_HLT_DiPFJetAve60_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve60_HFJEC_prescale", &HLT_DiPFJetAve60_HFJEC_prescale, &b_HLT_DiPFJetAve60_HFJEC_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC", &HLT_DiPFJetAve80_HFJEC, &b_HLT_DiPFJetAve80_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve80_HFJEC_prescale", &HLT_DiPFJetAve80_HFJEC_prescale, &b_HLT_DiPFJetAve80_HFJEC_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC", &HLT_DiPFJetAve100_HFJEC, &b_HLT_DiPFJetAve100_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve100_HFJEC_prescale", &HLT_DiPFJetAve100_HFJEC_prescale, &b_HLT_DiPFJetAve100_HFJEC_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC", &HLT_DiPFJetAve160_HFJEC, &b_HLT_DiPFJetAve160_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve160_HFJEC_prescale", &HLT_DiPFJetAve160_HFJEC_prescale, &b_HLT_DiPFJetAve160_HFJEC_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC", &HLT_DiPFJetAve220_HFJEC, &b_HLT_DiPFJetAve220_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve220_HFJEC_prescale", &HLT_DiPFJetAve220_HFJEC_prescale, &b_HLT_DiPFJetAve220_HFJEC_prescale);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC", &HLT_DiPFJetAve300_HFJEC, &b_HLT_DiPFJetAve300_HFJEC);
   fChain->SetBranchAddress("HLT_DiPFJetAve300_HFJEC_prescale", &HLT_DiPFJetAve300_HFJEC_prescale, &b_HLT_DiPFJetAve300_HFJEC_prescale);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45, &b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45);
   fChain->SetBranchAddress("HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale", &HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale, &b_HLT_AK8DiPFJet250_200_TrimMass30_BTagCSV0p45_prescale);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45);
   fChain->SetBranchAddress("HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale", &HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale, &b_HLT_AK8DiPFJet280_200_TrimMass30_BTagCSV0p45_prescale);
   fChain->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45, &b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45);
   fChain->SetBranchAddress("HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale", &HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale, &b_HLT_AK8PFHT600_TrimR0p1PT0p03Mass50_BTagCSV0p45_prescale);
   fChain->SetBranchAddress("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT650_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50);
   fChain->SetBranchAddress("HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_prescale", &HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_prescale, &b_HLT_AK8PFHT650_TrimR0p1PT0p03Mass50_prescale);
   fChain->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50, &b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50);
   fChain->SetBranchAddress("HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_prescale", &HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_prescale, &b_HLT_AK8PFHT700_TrimR0p1PT0p03Mass50_prescale);
   fChain->SetBranchAddress("HLT_PFHT550_4JetPt50", &HLT_PFHT550_4JetPt50, &b_HLT_PFHT550_4JetPt50);
   fChain->SetBranchAddress("HLT_PFHT550_4JetPt50_prescale", &HLT_PFHT550_4JetPt50_prescale, &b_HLT_PFHT550_4JetPt50_prescale);
   fChain->SetBranchAddress("HLT_PFHT650_4JetPt50", &HLT_PFHT650_4JetPt50, &b_HLT_PFHT650_4JetPt50);
   fChain->SetBranchAddress("HLT_PFHT650_4JetPt50_prescale", &HLT_PFHT650_4JetPt50_prescale, &b_HLT_PFHT650_4JetPt50_prescale);
   fChain->SetBranchAddress("HLT_PFHT750_4JetPt50", &HLT_PFHT750_4JetPt50, &b_HLT_PFHT750_4JetPt50);
   fChain->SetBranchAddress("HLT_PFHT750_4JetPt50_prescale", &HLT_PFHT750_4JetPt50_prescale, &b_HLT_PFHT750_4JetPt50_prescale);
   fChain->SetBranchAddress("HLT_ECALHT800", &HLT_ECALHT800, &b_HLT_ECALHT800);
   fChain->SetBranchAddress("HLT_ECALHT800_prescale", &HLT_ECALHT800_prescale, &b_HLT_ECALHT800_prescale);
   fChain->SetBranchAddress("HLT_PFHT600", &HLT_PFHT600, &b_HLT_PFHT600);
   fChain->SetBranchAddress("HLT_PFHT600_prescale", &HLT_PFHT600_prescale, &b_HLT_PFHT600_prescale);
   fChain->SetBranchAddress("HLT_PFHT650", &HLT_PFHT650, &b_HLT_PFHT650);
   fChain->SetBranchAddress("HLT_PFHT650_prescale", &HLT_PFHT650_prescale, &b_HLT_PFHT650_prescale);
   fChain->SetBranchAddress("HLT_PFHT800", &HLT_PFHT800, &b_HLT_PFHT800);
   fChain->SetBranchAddress("HLT_PFHT800_prescale", &HLT_PFHT800_prescale, &b_HLT_PFHT800_prescale);
   fChain->SetBranchAddress("HLT_PFHT200", &HLT_PFHT200, &b_HLT_PFHT200);
   fChain->SetBranchAddress("HLT_PFHT200_prescale", &HLT_PFHT200_prescale, &b_HLT_PFHT200_prescale);
   fChain->SetBranchAddress("HLT_PFHT250", &HLT_PFHT250, &b_HLT_PFHT250);
   fChain->SetBranchAddress("HLT_PFHT250_prescale", &HLT_PFHT250_prescale, &b_HLT_PFHT250_prescale);
   fChain->SetBranchAddress("HLT_PFHT300", &HLT_PFHT300, &b_HLT_PFHT300);
   fChain->SetBranchAddress("HLT_PFHT300_prescale", &HLT_PFHT300_prescale, &b_HLT_PFHT300_prescale);
   fChain->SetBranchAddress("HLT_PFHT350", &HLT_PFHT350, &b_HLT_PFHT350);
   fChain->SetBranchAddress("HLT_PFHT350_prescale", &HLT_PFHT350_prescale, &b_HLT_PFHT350_prescale);
   fChain->SetBranchAddress("HLT_PFHT400", &HLT_PFHT400, &b_HLT_PFHT400);
   fChain->SetBranchAddress("HLT_PFHT400_prescale", &HLT_PFHT400_prescale, &b_HLT_PFHT400_prescale);
   fChain->SetBranchAddress("HLT_PFHT475", &HLT_PFHT475, &b_HLT_PFHT475);
   fChain->SetBranchAddress("HLT_PFHT475_prescale", &HLT_PFHT475_prescale, &b_HLT_PFHT475_prescale);
   fChain->SetBranchAddress("HLT_Rsq0p25", &HLT_Rsq0p25, &b_HLT_Rsq0p25);
   fChain->SetBranchAddress("HLT_Rsq0p25_prescale", &HLT_Rsq0p25_prescale, &b_HLT_Rsq0p25_prescale);
   fChain->SetBranchAddress("HLT_Rsq0p30", &HLT_Rsq0p30, &b_HLT_Rsq0p30);
   fChain->SetBranchAddress("HLT_Rsq0p30_prescale", &HLT_Rsq0p30_prescale, &b_HLT_Rsq0p30_prescale);
   fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200", &HLT_RsqMR240_Rsq0p09_MR200, &b_HLT_RsqMR240_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_prescale", &HLT_RsqMR240_Rsq0p09_MR200_prescale, &b_HLT_RsqMR240_Rsq0p09_MR200_prescale);
   fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet", &HLT_RsqMR240_Rsq0p09_MR200_4jet, &b_HLT_RsqMR240_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_RsqMR240_Rsq0p09_MR200_4jet_prescale", &HLT_RsqMR240_Rsq0p09_MR200_4jet_prescale, &b_HLT_RsqMR240_Rsq0p09_MR200_4jet_prescale);
   fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200", &HLT_RsqMR270_Rsq0p09_MR200, &b_HLT_RsqMR270_Rsq0p09_MR200);
   fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_prescale", &HLT_RsqMR270_Rsq0p09_MR200_prescale, &b_HLT_RsqMR270_Rsq0p09_MR200_prescale);
   fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_4jet", &HLT_RsqMR270_Rsq0p09_MR200_4jet, &b_HLT_RsqMR270_Rsq0p09_MR200_4jet);
   fChain->SetBranchAddress("HLT_RsqMR270_Rsq0p09_MR200_4jet_prescale", &HLT_RsqMR270_Rsq0p09_MR200_4jet_prescale, &b_HLT_RsqMR270_Rsq0p09_MR200_4jet_prescale);
   fChain->SetBranchAddress("HLT_Mu30_eta2p1_PFJet150_PFJet50", &HLT_Mu30_eta2p1_PFJet150_PFJet50, &b_HLT_Mu30_eta2p1_PFJet150_PFJet50);
   fChain->SetBranchAddress("HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale", &HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale, &b_HLT_Mu30_eta2p1_PFJet150_PFJet50_prescale);
   fChain->SetBranchAddress("HLT_Mu40_eta2p1_PFJet200_PFJet50", &HLT_Mu40_eta2p1_PFJet200_PFJet50, &b_HLT_Mu40_eta2p1_PFJet200_PFJet50);
   fChain->SetBranchAddress("HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale", &HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale, &b_HLT_Mu40_eta2p1_PFJet200_PFJet50_prescale);
   fChain->SetBranchAddress("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50", &HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50, &b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50);
   fChain->SetBranchAddress("HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale", &HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale, &b_HLT_Ele35_CaloIdVT_GsfTrkIdT_PFJet150_PFJet50_prescale);
   fChain->SetBranchAddress("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50", &HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50, &b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50);
   fChain->SetBranchAddress("HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale", &HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale, &b_HLT_Ele45_CaloIdVT_GsfTrkIdT_PFJet200_PFJet50_prescale);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale", &HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale, &b_HLT_Ele8_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdM_TrackIdM_PFJet30", &HLT_Ele12_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale", &HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale, &b_HLT_Ele12_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele18_CaloIdM_TrackIdM_PFJet30", &HLT_Ele18_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele18_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele18_CaloIdM_TrackIdM_PFJet30_prescale", &HLT_Ele18_CaloIdM_TrackIdM_PFJet30_prescale, &b_HLT_Ele18_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale", &HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale, &b_HLT_Ele23_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele33_CaloIdM_TrackIdM_PFJet30", &HLT_Ele33_CaloIdM_TrackIdM_PFJet30, &b_HLT_Ele33_CaloIdM_TrackIdM_PFJet30);
   fChain->SetBranchAddress("HLT_Ele33_CaloIdM_TrackIdM_PFJet30_prescale", &HLT_Ele33_CaloIdM_TrackIdM_PFJet30_prescale, &b_HLT_Ele33_CaloIdM_TrackIdM_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30", &HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30, &b_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30);
   fChain->SetBranchAddress("HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale", &HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale, &b_HLT_Ele33_CaloIdL_TrackIdL_IsoVL_PFJet30_prescale);
   fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07", &HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07, &b_HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07);
   fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale", &HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale, &b_HLT_Ele23_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale);
   fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07", &HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07, &b_HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07);
   fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale", &HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale, &b_HLT_Ele27_WPLoose_Gsf_CentralPFJet30_BTagCSV07_prescale);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200", &HLT_Ele27_eta2p1_WPLoose_Gsf_HT200, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_prescale", &HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_prescale, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_HT200_prescale);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF", &HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF, &b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF);
   fChain->SetBranchAddress("HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale", &HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale, &b_HLT_Ele10_CaloIdM_TrackIdM_CentralPFJet30_BTagCSV0p54PF_prescale);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400", &HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400, &b_HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale", &HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale, &b_HLT_Ele15_IsoVVVL_BTagCSV0p72_PFHT400_prescale);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_PFMET50", &HLT_Ele15_IsoVVVL_PFHT350_PFMET50, &b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_PFMET50_prescale", &HLT_Ele15_IsoVVVL_PFHT350_PFMET50_prescale, &b_HLT_Ele15_IsoVVVL_PFHT350_PFMET50_prescale);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600", &HLT_Ele15_IsoVVVL_PFHT600, &b_HLT_Ele15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT600_prescale", &HLT_Ele15_IsoVVVL_PFHT600_prescale, &b_HLT_Ele15_IsoVVVL_PFHT600_prescale);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350", &HLT_Ele15_IsoVVVL_PFHT350, &b_HLT_Ele15_IsoVVVL_PFHT350);
   fChain->SetBranchAddress("HLT_Ele15_IsoVVVL_PFHT350_prescale", &HLT_Ele15_IsoVVVL_PFHT350_prescale, &b_HLT_Ele15_IsoVVVL_PFHT350_prescale);
   fChain->SetBranchAddress("HLT_Mu3er_PFHT140_PFMET125", &HLT_Mu3er_PFHT140_PFMET125, &b_HLT_Mu3er_PFHT140_PFMET125);
   fChain->SetBranchAddress("HLT_Mu3er_PFHT140_PFMET125_prescale", &HLT_Mu3er_PFHT140_PFMET125_prescale, &b_HLT_Mu3er_PFHT140_PFMET125_prescale);
   fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72", &HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72, &b_HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72);
   fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale", &HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale, &b_HLT_Mu6_PFHT200_PFMET80_BTagCSV0p72_prescale);
   fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET100", &HLT_Mu6_PFHT200_PFMET100, &b_HLT_Mu6_PFHT200_PFMET100);
   fChain->SetBranchAddress("HLT_Mu6_PFHT200_PFMET100_prescale", &HLT_Mu6_PFHT200_PFMET100_prescale, &b_HLT_Mu6_PFHT200_PFMET100_prescale);
   fChain->SetBranchAddress("HLT_Mu10_CentralPFJet30_BTagCSV0p54PF", &HLT_Mu10_CentralPFJet30_BTagCSV0p54PF, &b_HLT_Mu10_CentralPFJet30_BTagCSV0p54PF);
   fChain->SetBranchAddress("HLT_Mu10_CentralPFJet30_BTagCSV0p54PF_prescale", &HLT_Mu10_CentralPFJet30_BTagCSV0p54PF_prescale, &b_HLT_Mu10_CentralPFJet30_BTagCSV0p54PF_prescale);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400", &HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400, &b_HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale", &HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale, &b_HLT_Mu15_IsoVVVL_BTagCSV0p72_PFHT400_prescale);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_PFMET50", &HLT_Mu15_IsoVVVL_PFHT350_PFMET50, &b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_PFMET50_prescale", &HLT_Mu15_IsoVVVL_PFHT350_PFMET50_prescale, &b_HLT_Mu15_IsoVVVL_PFHT350_PFMET50_prescale);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600", &HLT_Mu15_IsoVVVL_PFHT600, &b_HLT_Mu15_IsoVVVL_PFHT600);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT600_prescale", &HLT_Mu15_IsoVVVL_PFHT600_prescale, &b_HLT_Mu15_IsoVVVL_PFHT600_prescale);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350", &HLT_Mu15_IsoVVVL_PFHT350, &b_HLT_Mu15_IsoVVVL_PFHT350);
   fChain->SetBranchAddress("HLT_Mu15_IsoVVVL_PFHT350_prescale", &HLT_Mu15_IsoVVVL_PFHT350_prescale, &b_HLT_Mu15_IsoVVVL_PFHT350_prescale);
   fChain->SetBranchAddress("HLT_Mu8", &HLT_Mu8, &b_HLT_Mu8);
   fChain->SetBranchAddress("HLT_Mu8_prescale", &HLT_Mu8_prescale, &b_HLT_Mu8_prescale);
   fChain->SetBranchAddress("HLT_Mu17", &HLT_Mu17, &b_HLT_Mu17);
   fChain->SetBranchAddress("HLT_Mu17_prescale", &HLT_Mu17_prescale, &b_HLT_Mu17_prescale);
   fChain->SetBranchAddress("HLT_Mu20", &HLT_Mu20, &b_HLT_Mu20);
   fChain->SetBranchAddress("HLT_Mu20_prescale", &HLT_Mu20_prescale, &b_HLT_Mu20_prescale);
   fChain->SetBranchAddress("HLT_TkMu20", &HLT_TkMu20, &b_HLT_TkMu20);
   fChain->SetBranchAddress("HLT_TkMu20_prescale", &HLT_TkMu20_prescale, &b_HLT_TkMu20_prescale);
   fChain->SetBranchAddress("HLT_Mu24_eta2p1", &HLT_Mu24_eta2p1, &b_HLT_Mu24_eta2p1);
   fChain->SetBranchAddress("HLT_Mu24_eta2p1_prescale", &HLT_Mu24_eta2p1_prescale, &b_HLT_Mu24_eta2p1_prescale);
   fChain->SetBranchAddress("HLT_TkMu24_eta2p1", &HLT_TkMu24_eta2p1, &b_HLT_TkMu24_eta2p1);
   fChain->SetBranchAddress("HLT_TkMu24_eta2p1_prescale", &HLT_TkMu24_eta2p1_prescale, &b_HLT_TkMu24_eta2p1_prescale);
   fChain->SetBranchAddress("HLT_Mu27", &HLT_Mu27, &b_HLT_Mu27);
   fChain->SetBranchAddress("HLT_Mu27_prescale", &HLT_Mu27_prescale, &b_HLT_Mu27_prescale);
   fChain->SetBranchAddress("HLT_TkMu27", &HLT_TkMu27, &b_HLT_TkMu27);
   fChain->SetBranchAddress("HLT_TkMu27_prescale", &HLT_TkMu27_prescale, &b_HLT_TkMu27_prescale);
   fChain->SetBranchAddress("HLT_Mu45_eta2p1", &HLT_Mu45_eta2p1, &b_HLT_Mu45_eta2p1);
   fChain->SetBranchAddress("HLT_Mu45_eta2p1_prescale", &HLT_Mu45_eta2p1_prescale, &b_HLT_Mu45_eta2p1_prescale);
   fChain->SetBranchAddress("HLT_Mu50", &HLT_Mu50, &b_HLT_Mu50);
   fChain->SetBranchAddress("HLT_Mu50_prescale", &HLT_Mu50_prescale, &b_HLT_Mu50_prescale);
   fChain->SetBranchAddress("HLT_Mu55", &HLT_Mu55, &b_HLT_Mu55);
   fChain->SetBranchAddress("HLT_Mu55_prescale", &HLT_Mu55_prescale, &b_HLT_Mu55_prescale);
   fChain->SetBranchAddress("HLT_Mu300", &HLT_Mu300, &b_HLT_Mu300);
   fChain->SetBranchAddress("HLT_Mu300_prescale", &HLT_Mu300_prescale, &b_HLT_Mu300_prescale);
   fChain->SetBranchAddress("HLT_Mu350", &HLT_Mu350, &b_HLT_Mu350);
   fChain->SetBranchAddress("HLT_Mu350_prescale", &HLT_Mu350_prescale, &b_HLT_Mu350_prescale);
   fChain->SetBranchAddress("HLT_Ele105_CaloIdVT_GsfTrkIdT", &HLT_Ele105_CaloIdVT_GsfTrkIdT, &b_HLT_Ele105_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale", &HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale, &b_HLT_Ele105_CaloIdVT_GsfTrkIdT_prescale);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT", &HLT_Ele115_CaloIdVT_GsfTrkIdT, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT);
   fChain->SetBranchAddress("HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale", &HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale, &b_HLT_Ele115_CaloIdVT_GsfTrkIdT_prescale);
   fChain->SetBranchAddress("HLT_IsoMu17_eta2p1", &HLT_IsoMu17_eta2p1, &b_HLT_IsoMu17_eta2p1);
   fChain->SetBranchAddress("HLT_IsoMu17_eta2p1_prescale", &HLT_IsoMu17_eta2p1_prescale, &b_HLT_IsoMu17_eta2p1_prescale);
   fChain->SetBranchAddress("HLT_IsoMu18", &HLT_IsoMu18, &b_HLT_IsoMu18);
   fChain->SetBranchAddress("HLT_IsoMu18_prescale", &HLT_IsoMu18_prescale, &b_HLT_IsoMu18_prescale);
   fChain->SetBranchAddress("HLT_OldIsoMu18", &HLT_OldIsoMu18, &b_HLT_OldIsoMu18);
   fChain->SetBranchAddress("HLT_OldIsoMu18_prescale", &HLT_OldIsoMu18_prescale, &b_HLT_OldIsoMu18_prescale);
   fChain->SetBranchAddress("HLT_IsoTkMu18", &HLT_IsoTkMu18, &b_HLT_IsoTkMu18);
   fChain->SetBranchAddress("HLT_IsoTkMu18_prescale", &HLT_IsoTkMu18_prescale, &b_HLT_IsoTkMu18_prescale);
   fChain->SetBranchAddress("HLT_OldIsoTkMu18", &HLT_OldIsoTkMu18, &b_HLT_OldIsoTkMu18);
   fChain->SetBranchAddress("HLT_OldIsoTkMu18_prescale", &HLT_OldIsoTkMu18_prescale, &b_HLT_OldIsoTkMu18_prescale);
   fChain->SetBranchAddress("HLT_IsoMu20", &HLT_IsoMu20, &b_HLT_IsoMu20);
   fChain->SetBranchAddress("HLT_IsoMu20_prescale", &HLT_IsoMu20_prescale, &b_HLT_IsoMu20_prescale);
   fChain->SetBranchAddress("HLT_IsoMu22", &HLT_IsoMu22, &b_HLT_IsoMu22);
   fChain->SetBranchAddress("HLT_IsoMu22_prescale", &HLT_IsoMu22_prescale, &b_HLT_IsoMu22_prescale);
   fChain->SetBranchAddress("HLT_IsoMu27", &HLT_IsoMu27, &b_HLT_IsoMu27);
   fChain->SetBranchAddress("HLT_IsoMu27_prescale", &HLT_IsoMu27_prescale, &b_HLT_IsoMu27_prescale);
   fChain->SetBranchAddress("HLT_IsoTkMu24_eta2p1", &HLT_IsoTkMu24_eta2p1, &b_HLT_IsoTkMu24_eta2p1);
   fChain->SetBranchAddress("HLT_IsoTkMu24_eta2p1_prescale", &HLT_IsoTkMu24_eta2p1_prescale, &b_HLT_IsoTkMu24_eta2p1_prescale);
   fChain->SetBranchAddress("HLT_IsoTkMu27", &HLT_IsoTkMu27, &b_HLT_IsoTkMu27);
   fChain->SetBranchAddress("HLT_IsoTkMu27_prescale", &HLT_IsoTkMu27_prescale, &b_HLT_IsoTkMu27_prescale);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL", &HLT_Mu8_TrkIsoVVL, &b_HLT_Mu8_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu8_TrkIsoVVL_prescale", &HLT_Mu8_TrkIsoVVL_prescale, &b_HLT_Mu8_TrkIsoVVL_prescale);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL", &HLT_Mu17_TrkIsoVVL, &b_HLT_Mu17_TrkIsoVVL);
   fChain->SetBranchAddress("HLT_Mu17_TrkIsoVVL_prescale", &HLT_Mu17_TrkIsoVVL_prescale, &b_HLT_Mu17_TrkIsoVVL_prescale);
   fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf", &HLT_Ele22_eta2p1_WPLoose_Gsf, &b_HLT_Ele22_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPLoose_Gsf_prescale", &HLT_Ele22_eta2p1_WPLoose_Gsf_prescale, &b_HLT_Ele22_eta2p1_WPLoose_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPTight_Gsf", &HLT_Ele22_eta2p1_WPTight_Gsf, &b_HLT_Ele22_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele22_eta2p1_WPTight_Gsf_prescale", &HLT_Ele22_eta2p1_WPTight_Gsf_prescale, &b_HLT_Ele22_eta2p1_WPTight_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf", &HLT_Ele23_WPLoose_Gsf, &b_HLT_Ele23_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele23_WPLoose_Gsf_prescale", &HLT_Ele23_WPLoose_Gsf_prescale, &b_HLT_Ele23_WPLoose_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf", &HLT_Ele27_WPLoose_Gsf, &b_HLT_Ele27_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_WPLoose_Gsf_prescale", &HLT_Ele27_WPLoose_Gsf_prescale, &b_HLT_Ele27_WPLoose_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf", &HLT_Ele27_eta2p1_WPLoose_Gsf, &b_HLT_Ele27_eta2p1_WPLoose_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPLoose_Gsf_prescale", &HLT_Ele27_eta2p1_WPLoose_Gsf_prescale, &b_HLT_Ele27_eta2p1_WPLoose_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPTight_Gsf", &HLT_Ele27_eta2p1_WPTight_Gsf, &b_HLT_Ele27_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele27_eta2p1_WPTight_Gsf_prescale", &HLT_Ele27_eta2p1_WPTight_Gsf_prescale, &b_HLT_Ele27_eta2p1_WPTight_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf", &HLT_Ele32_eta2p1_WPTight_Gsf, &b_HLT_Ele32_eta2p1_WPTight_Gsf);
   fChain->SetBranchAddress("HLT_Ele32_eta2p1_WPTight_Gsf_prescale", &HLT_Ele32_eta2p1_WPTight_Gsf_prescale, &b_HLT_Ele32_eta2p1_WPTight_Gsf_prescale);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele12_CaloIdL_TrackIdL_IsoVL_prescale", &HLT_Ele12_CaloIdL_TrackIdL_IsoVL_prescale, &b_HLT_Ele12_CaloIdL_TrackIdL_IsoVL_prescale);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele17_CaloIdL_TrackIdL_IsoVL_prescale", &HLT_Ele17_CaloIdL_TrackIdL_IsoVL_prescale, &b_HLT_Ele17_CaloIdL_TrackIdL_IsoVL_prescale);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL);
   fChain->SetBranchAddress("HLT_Ele23_CaloIdL_TrackIdL_IsoVL_prescale", &HLT_Ele23_CaloIdL_TrackIdL_IsoVL_prescale, &b_HLT_Ele23_CaloIdL_TrackIdL_IsoVL_prescale);
   fChain->SetBranchAddress("evt_NGoodVtx", &evt_NGoodVtx, &b_evt_NGoodVtx);
   fChain->SetBranchAddress("evt_NLep", &evt_NLep, &b_evt_NLep);
   fChain->SetBranchAddress("evt_NTopHad", &evt_NTopHad, &b_evt_NTopHad);
   fChain->SetBranchAddress("evt_NTopHadPreTag", &evt_NTopHadPreTag, &b_evt_NTopHadPreTag);
   fChain->SetBranchAddress("evt_NTopLep", &evt_NTopLep, &b_evt_NTopLep);
   fChain->SetBranchAddress("evt_NTop", &evt_NTop, &b_evt_NTop);
   fChain->SetBranchAddress("evt_LHA_PDF_ID", &evt_LHA_PDF_ID, &b_evt_LHA_PDF_ID);
   fChain->SetBranchAddress("evt_HtLep", &evt_HtLep, &b_evt_HtLep);
   fChain->SetBranchAddress("evt_HtTop", &evt_HtTop, &b_evt_HtTop);
   fChain->SetBranchAddress("evt_Ht", &evt_Ht, &b_evt_Ht);
   fChain->SetBranchAddress("evt_HtAll", &evt_HtAll, &b_evt_HtAll);
   fChain->SetBranchAddress("evt_HtEx", &evt_HtEx, &b_evt_HtEx);
   fChain->SetBranchAddress("evt_HtExFr", &evt_HtExFr, &b_evt_HtExFr);
   fChain->SetBranchAddress("evt_HtTopFr", &evt_HtTopFr, &b_evt_HtTopFr);
   fChain->SetBranchAddress("evt_TTHadDR", &evt_TTHadDR, &b_evt_TTHadDR);
   fChain->SetBranchAddress("evt_TTHadDPhi", &evt_TTHadDPhi, &b_evt_TTHadDPhi);
   fChain->SetBranchAddress("evt_TTHadDEta", &evt_TTHadDEta, &b_evt_TTHadDEta);
   fChain->SetBranchAddress("evt_TTHadMass", &evt_TTHadMass, &b_evt_TTHadMass);
   fChain->SetBranchAddress("evt_TTHadSumPt", &evt_TTHadSumPt, &b_evt_TTHadSumPt);
   fChain->SetBranchAddress("evt_TTHadPz", &evt_TTHadPz, &b_evt_TTHadPz);
   fChain->SetBranchAddress("evt_TTHadHz", &evt_TTHadHz, &b_evt_TTHadHz);
   fChain->SetBranchAddress("evt_TTHadDPz", &evt_TTHadDPz, &b_evt_TTHadDPz);
   fChain->SetBranchAddress("evt_TTHadMR", &evt_TTHadMR, &b_evt_TTHadMR);
   fChain->SetBranchAddress("evt_TTHadMTR", &evt_TTHadMTR, &b_evt_TTHadMTR);
   fChain->SetBranchAddress("evt_TTHadR", &evt_TTHadR, &b_evt_TTHadR);
   fChain->SetBranchAddress("evt_TTHadR2", &evt_TTHadR2, &b_evt_TTHadR2);
   fChain->SetBranchAddress("evt_MR", &evt_MR, &b_evt_MR);
   fChain->SetBranchAddress("evt_MTR", &evt_MTR, &b_evt_MTR);
   fChain->SetBranchAddress("evt_R", &evt_R, &b_evt_R);
   fChain->SetBranchAddress("evt_R2", &evt_R2, &b_evt_R2);
   fChain->SetBranchAddress("evt_AK8Puppi_MR", &evt_AK8Puppi_MR, &b_evt_AK8Puppi_MR);
   fChain->SetBranchAddress("evt_AK8Puppi_MTR", &evt_AK8Puppi_MTR, &b_evt_AK8Puppi_MTR);
   fChain->SetBranchAddress("evt_AK8Puppi_R", &evt_AK8Puppi_R, &b_evt_AK8Puppi_R);
   fChain->SetBranchAddress("evt_AK8Puppi_R2", &evt_AK8Puppi_R2, &b_evt_AK8Puppi_R2);
   fChain->SetBranchAddress("evt_AK4_MR", &evt_AK4_MR, &b_evt_AK4_MR);
   fChain->SetBranchAddress("evt_AK4_MTR", &evt_AK4_MTR, &b_evt_AK4_MTR);
   fChain->SetBranchAddress("evt_AK4_R", &evt_AK4_R, &b_evt_AK4_R);
   fChain->SetBranchAddress("evt_AK4_R2", &evt_AK4_R2, &b_evt_AK4_R2);
   fChain->SetBranchAddress("evt_XSec", &evt_XSec, &b_evt_XSec);
   fChain->SetBranchAddress("evt_NEvent_Corr", &evt_NEvent_Corr, &b_evt_NEvent_Corr);
   fChain->SetBranchAddress("evt_Lumi_Weight", &evt_Lumi_Weight, &b_evt_Lumi_Weight);
   fChain->SetBranchAddress("evt_Gen_Weight", &evt_Gen_Weight, &b_evt_Gen_Weight);
   fChain->SetBranchAddress("evt_Gen_Ht", &evt_Gen_Ht, &b_evt_Gen_Ht);
   fChain->SetBranchAddress("SUSY_Gluino_Mass", &SUSY_Gluino_Mass, &b_SUSY_Gluino_Mass);
   fChain->SetBranchAddress("SUSY_LSP_Mass", &SUSY_LSP_Mass, &b_SUSY_LSP_Mass);
   fChain->SetBranchAddress("gen_size", &gen_size, &b_gen_size);
   fChain->SetBranchAddress("gen_ID", gen_ID, &b_gen_ID);
   fChain->SetBranchAddress("gen_Status", gen_Status, &b_gen_Status);
   fChain->SetBranchAddress("gen_Mom0ID", gen_Mom0ID, &b_gen_Mom0ID);
   fChain->SetBranchAddress("gen_Mom0Status", gen_Mom0Status, &b_gen_Mom0Status);
   fChain->SetBranchAddress("gen_Mom1ID", gen_Mom1ID, &b_gen_Mom1ID);
   fChain->SetBranchAddress("gen_Mom1Status", gen_Mom1Status, &b_gen_Mom1Status);
   fChain->SetBranchAddress("gen_Dau0ID", gen_Dau0ID, &b_gen_Dau0ID);
   fChain->SetBranchAddress("gen_Dau0Status", gen_Dau0Status, &b_gen_Dau0Status);
   fChain->SetBranchAddress("gen_Dau1ID", gen_Dau1ID, &b_gen_Dau1ID);
   fChain->SetBranchAddress("gen_Dau1Status", gen_Dau1Status, &b_gen_Dau1Status);
   fChain->SetBranchAddress("jetAK8_HasNearGenTop", jetAK8_HasNearGenTop, &b_jetAK8_HasNearGenTop);
   fChain->SetBranchAddress("jetAK8_NearGenTopIsHadronic", jetAK8_NearGenTopIsHadronic, &b_jetAK8_NearGenTopIsHadronic);
   fChain->SetBranchAddress("jetAK8_NearGenWIsHadronic", jetAK8_NearGenWIsHadronic, &b_jetAK8_NearGenWIsHadronic);
   fChain->SetBranchAddress("jetAK8_NearGenWToENu", jetAK8_NearGenWToENu, &b_jetAK8_NearGenWToENu);
   fChain->SetBranchAddress("jetAK8_NearGenWToMuNu", jetAK8_NearGenWToMuNu, &b_jetAK8_NearGenWToMuNu);
   fChain->SetBranchAddress("jetAK8_NearGenWToTauNu", jetAK8_NearGenWToTauNu, &b_jetAK8_NearGenWToTauNu);
   fChain->SetBranchAddress("jetAK8_PassTopTag", jetAK8_PassTopTag, &b_jetAK8_PassTopTag);
   fChain->SetBranchAddress("el_IsPartOfNearAK4Jet", el_IsPartOfNearAK4Jet, &b_el_IsPartOfNearAK4Jet);
   fChain->SetBranchAddress("el_IsPartOfNearAK8Jet", el_IsPartOfNearAK8Jet, &b_el_IsPartOfNearAK8Jet);
   fChain->SetBranchAddress("el_IsPartOfNearSubjet", el_IsPartOfNearSubjet, &b_el_IsPartOfNearSubjet);
   fChain->SetBranchAddress("mu_IsPartOfNearAK4Jet", mu_IsPartOfNearAK4Jet, &b_mu_IsPartOfNearAK4Jet);
   fChain->SetBranchAddress("mu_IsPartOfNearAK8Jet", mu_IsPartOfNearAK8Jet, &b_mu_IsPartOfNearAK8Jet);
   fChain->SetBranchAddress("mu_IsPartOfNearSubjet", mu_IsPartOfNearSubjet, &b_mu_IsPartOfNearSubjet);
   fChain->SetBranchAddress("scale_size", &scale_size, &b_scale_size);
   fChain->SetBranchAddress("scale_Weights", scale_Weights, &b_scale_Weights);
   fChain->SetBranchAddress("pdf_size", &pdf_size, &b_pdf_size);
   fChain->SetBranchAddress("pdf_Weights", pdf_Weights, &b_pdf_Weights);
   fChain->SetBranchAddress("alphas_size", &alphas_size, &b_alphas_size);
   fChain->SetBranchAddress("alphas_Weights", alphas_Weights, &b_alphas_Weights);
   fChain->SetBranchAddress("gen_Pt", gen_Pt, &b_gen_Pt);
   fChain->SetBranchAddress("gen_Eta", gen_Eta, &b_gen_Eta);
   fChain->SetBranchAddress("gen_Phi", gen_Phi, &b_gen_Phi);
   fChain->SetBranchAddress("gen_E", gen_E, &b_gen_E);
   fChain->SetBranchAddress("gen_Mass", gen_Mass, &b_gen_Mass);
   fChain->SetBranchAddress("gen_Charge", gen_Charge, &b_gen_Charge);
   fChain->SetBranchAddress("jetAK8_DRNearGenTop", jetAK8_DRNearGenTop, &b_jetAK8_DRNearGenTop);
   fChain->SetBranchAddress("jetAK8_DRNearGenWFromTop", jetAK8_DRNearGenWFromTop, &b_jetAK8_DRNearGenWFromTop);
   fChain->SetBranchAddress("jetAK8_DRNearGenBFromTop", jetAK8_DRNearGenBFromTop, &b_jetAK8_DRNearGenBFromTop);
   fChain->SetBranchAddress("jetAK8_DRNearGenLepFromSLTop", jetAK8_DRNearGenLepFromSLTop, &b_jetAK8_DRNearGenLepFromSLTop);
   fChain->SetBranchAddress("jetAK8_DRNearGenNuFromSLTop", jetAK8_DRNearGenNuFromSLTop, &b_jetAK8_DRNearGenNuFromSLTop);
   fChain->SetBranchAddress("jetAK8_PtNearGenTop", jetAK8_PtNearGenTop, &b_jetAK8_PtNearGenTop);
   fChain->SetBranchAddress("jetAK8_PtNearGenBFromTop", jetAK8_PtNearGenBFromTop, &b_jetAK8_PtNearGenBFromTop);
   fChain->SetBranchAddress("jetAK8_PtNearGenWFromTop", jetAK8_PtNearGenWFromTop, &b_jetAK8_PtNearGenWFromTop);
   fChain->SetBranchAddress("jetAK8_PtNearGenLepFromSLTop", jetAK8_PtNearGenLepFromSLTop, &b_jetAK8_PtNearGenLepFromSLTop);
   fChain->SetBranchAddress("jetAK8_PtNearGenNuFromSLTop", jetAK8_PtNearGenNuFromSLTop, &b_jetAK8_PtNearGenNuFromSLTop);
   fChain->SetBranchAddress("el_DRNearGenEleFromSLTop", el_DRNearGenEleFromSLTop, &b_el_DRNearGenEleFromSLTop);
   fChain->SetBranchAddress("el_PtNearGenEleFromSLTop", el_PtNearGenEleFromSLTop, &b_el_PtNearGenEleFromSLTop);
   fChain->SetBranchAddress("el_PtNearGenTop", el_PtNearGenTop, &b_el_PtNearGenTop);
   fChain->SetBranchAddress("el_LepAK4JetFrac", el_LepAK4JetFrac, &b_el_LepAK4JetFrac);
   fChain->SetBranchAddress("el_LepAK8JetFrac", el_LepAK8JetFrac, &b_el_LepAK8JetFrac);
   fChain->SetBranchAddress("el_LepSubjetFrac", el_LepSubjetFrac, &b_el_LepSubjetFrac);
   fChain->SetBranchAddress("el_LepAK4JetMassDrop", el_LepAK4JetMassDrop, &b_el_LepAK4JetMassDrop);
   fChain->SetBranchAddress("el_LepAK8JetMassDrop", el_LepAK8JetMassDrop, &b_el_LepAK8JetMassDrop);
   fChain->SetBranchAddress("el_LepSubjetMassDrop", el_LepSubjetMassDrop, &b_el_LepSubjetMassDrop);
   fChain->SetBranchAddress("el_AK4JetV1DR", el_AK4JetV1DR, &b_el_AK4JetV1DR);
   fChain->SetBranchAddress("el_AK4JetV2DR", el_AK4JetV2DR, &b_el_AK4JetV2DR);
   fChain->SetBranchAddress("el_AK4JetV3DR", el_AK4JetV3DR, &b_el_AK4JetV3DR);
   fChain->SetBranchAddress("el_AK8JetV1DR", el_AK8JetV1DR, &b_el_AK8JetV1DR);
   fChain->SetBranchAddress("el_AK8JetV2DR", el_AK8JetV2DR, &b_el_AK8JetV2DR);
   fChain->SetBranchAddress("el_AK8JetV3DR", el_AK8JetV3DR, &b_el_AK8JetV3DR);
   fChain->SetBranchAddress("el_SubjetV1DR", el_SubjetV1DR, &b_el_SubjetV1DR);
   fChain->SetBranchAddress("el_SubjetV2DR", el_SubjetV2DR, &b_el_SubjetV2DR);
   fChain->SetBranchAddress("el_SubjetV3DR", el_SubjetV3DR, &b_el_SubjetV3DR);
   fChain->SetBranchAddress("el_AK4JetV1PtRel", el_AK4JetV1PtRel, &b_el_AK4JetV1PtRel);
   fChain->SetBranchAddress("el_AK4JetV2PtRel", el_AK4JetV2PtRel, &b_el_AK4JetV2PtRel);
   fChain->SetBranchAddress("el_AK4JetV3PtRel", el_AK4JetV3PtRel, &b_el_AK4JetV3PtRel);
   fChain->SetBranchAddress("el_AK8JetV1PtRel", el_AK8JetV1PtRel, &b_el_AK8JetV1PtRel);
   fChain->SetBranchAddress("el_AK8JetV2PtRel", el_AK8JetV2PtRel, &b_el_AK8JetV2PtRel);
   fChain->SetBranchAddress("el_AK8JetV3PtRel", el_AK8JetV3PtRel, &b_el_AK8JetV3PtRel);
   fChain->SetBranchAddress("el_SubjetV1PtRel", el_SubjetV1PtRel, &b_el_SubjetV1PtRel);
   fChain->SetBranchAddress("el_SubjetV2PtRel", el_SubjetV2PtRel, &b_el_SubjetV2PtRel);
   fChain->SetBranchAddress("el_SubjetV3PtRel", el_SubjetV3PtRel, &b_el_SubjetV3PtRel);
   fChain->SetBranchAddress("mu_DRNearGenMuFromSLTop", mu_DRNearGenMuFromSLTop, &b_mu_DRNearGenMuFromSLTop);
   fChain->SetBranchAddress("mu_PtNearGenMuFromSLTop", mu_PtNearGenMuFromSLTop, &b_mu_PtNearGenMuFromSLTop);
   fChain->SetBranchAddress("mu_PtNearGenTop", mu_PtNearGenTop, &b_mu_PtNearGenTop);
   fChain->SetBranchAddress("mu_LepAK4JetFrac", mu_LepAK4JetFrac, &b_mu_LepAK4JetFrac);
   fChain->SetBranchAddress("mu_LepAK8JetFrac", mu_LepAK8JetFrac, &b_mu_LepAK8JetFrac);
   fChain->SetBranchAddress("mu_LepSubjetFrac", mu_LepSubjetFrac, &b_mu_LepSubjetFrac);
   fChain->SetBranchAddress("mu_LepAK4JetMassDrop", mu_LepAK4JetMassDrop, &b_mu_LepAK4JetMassDrop);
   fChain->SetBranchAddress("mu_LepAK8JetMassDrop", mu_LepAK8JetMassDrop, &b_mu_LepAK8JetMassDrop);
   fChain->SetBranchAddress("mu_LepSubjetMassDrop", mu_LepSubjetMassDrop, &b_mu_LepSubjetMassDrop);
   fChain->SetBranchAddress("mu_AK4JetV1DR", mu_AK4JetV1DR, &b_mu_AK4JetV1DR);
   fChain->SetBranchAddress("mu_AK4JetV2DR", mu_AK4JetV2DR, &b_mu_AK4JetV2DR);
   fChain->SetBranchAddress("mu_AK4JetV3DR", mu_AK4JetV3DR, &b_mu_AK4JetV3DR);
   fChain->SetBranchAddress("mu_AK8JetV1DR", mu_AK8JetV1DR, &b_mu_AK8JetV1DR);
   fChain->SetBranchAddress("mu_AK8JetV2DR", mu_AK8JetV2DR, &b_mu_AK8JetV2DR);
   fChain->SetBranchAddress("mu_AK8JetV3DR", mu_AK8JetV3DR, &b_mu_AK8JetV3DR);
   fChain->SetBranchAddress("mu_SubjetV1DR", mu_SubjetV1DR, &b_mu_SubjetV1DR);
   fChain->SetBranchAddress("mu_SubjetV2DR", mu_SubjetV2DR, &b_mu_SubjetV2DR);
   fChain->SetBranchAddress("mu_SubjetV3DR", mu_SubjetV3DR, &b_mu_SubjetV3DR);
   fChain->SetBranchAddress("mu_AK4JetV1PtRel", mu_AK4JetV1PtRel, &b_mu_AK4JetV1PtRel);
   fChain->SetBranchAddress("mu_AK4JetV2PtRel", mu_AK4JetV2PtRel, &b_mu_AK4JetV2PtRel);
   fChain->SetBranchAddress("mu_AK4JetV3PtRel", mu_AK4JetV3PtRel, &b_mu_AK4JetV3PtRel);
   fChain->SetBranchAddress("mu_AK8JetV1PtRel", mu_AK8JetV1PtRel, &b_mu_AK8JetV1PtRel);
   fChain->SetBranchAddress("mu_AK8JetV2PtRel", mu_AK8JetV2PtRel, &b_mu_AK8JetV2PtRel);
   fChain->SetBranchAddress("mu_AK8JetV3PtRel", mu_AK8JetV3PtRel, &b_mu_AK8JetV3PtRel);
   fChain->SetBranchAddress("mu_SubjetV1PtRel", mu_SubjetV1PtRel, &b_mu_SubjetV1PtRel);
   fChain->SetBranchAddress("mu_SubjetV2PtRel", mu_SubjetV2PtRel, &b_mu_SubjetV2PtRel);
   fChain->SetBranchAddress("mu_SubjetV3PtRel", mu_SubjetV3PtRel, &b_mu_SubjetV3PtRel);
   Notify();
}

Bool_t B2GTTree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void B2GTTree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t B2GTTree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef B2GTTree_cxx
