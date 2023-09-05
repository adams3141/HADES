//includes
//kinfit
#include "Includes.h" // A lot of header files from Root and Hydra
#include "KinFitter.h"
#include "KFitVertexFinder.h"
#include "KFitNeutralCandFinder.h"
//kinfit


#include "ecal_tests.h"
#include "emcdef.h"

#include <hades.h>
#include <hdst.h>
#include <hparasciifileio.h>
#include <hspectrometer.h>
#include <hparticlecand.h>
#include "hemcneutralcand.h"
#include "hemcclustersim.h"
#include "hemccluster.h"

#include <hgeomvolume.h>
#include <hgeomcompositevolume.h>
#include <hgeomvector.h>
#include "hparticleevtinfo.h"
#include "hparticletracksorter.h"
#include "hphysicsconstants.h"
#include "htool.h"
#include "htime.h"
#include "heventmixer.h"
#include "TMultiEventMixer.h" 


#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TH2I.h>
#include <TH3I.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TMatrixD.h>
#include <TVector3.h>
#include <TMath.h>

#include <algorithm>
#include <cstdlib>
#include <vector>

#define PI 3.14159265

long long int total=0;
long long int refited=0;
double trackDistance(HParticleCand* track1, HParticleCand*  track2)
    {
      double dist;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      dist=p_tool.calculateMinimumDistance(base_1,dir_1,base_2,dir_2);
      return dist;
    }

    HGeomVector trackVertex(HParticleCand* track1, HParticleCand*  track2)
    {
      HGeomVector ver;
      HGeomVector base_1, base_2, dir_1, dir_2;
      HParticleTool p_tool;

      p_tool.calcSegVector(track1->getZ(),track1->getR(),TMath::DegToRad()*track1->getPhi(),TMath::DegToRad()*track1->getTheta(),base_1,dir_1);
      p_tool.calcSegVector(track2->getZ(),track2->getR(),TMath::DegToRad()*track2->getPhi(),TMath::DegToRad()*track2->getTheta(),base_2,dir_2);
      ver=p_tool.calcVertexAnalytical(base_1,dir_1,base_2,dir_2);
      return ver;
    }
#define PR(x)                                                                                      \
    std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229;

TLorentzVector TLV(float th, float ph, float ek){


    float px, py, pz, mom;
    float m=0.; 
    TLorentzVector lv;
    mom=ek;
    
    px = mom*TMath::Sin(th*TMath::DegToRad())*TMath::Cos(ph*TMath::DegToRad());
    py = mom*TMath::Sin(th*TMath::DegToRad())*TMath::Sin(ph*TMath::DegToRad());
    pz = mom*TMath::Cos(th*TMath::DegToRad());
    lv.SetPxPyPzE(px,py,pz,ek);

    return lv;
}

TLorentzVector TLV_p(float th, float ph, float mom, float mass){

    float px, py, pz, E;
    
    TLorentzVector lv;
    
    px = mom*TMath::Sin(th*TMath::DegToRad())*TMath::Cos(ph*TMath::DegToRad());
    py = mom*TMath::Sin(th*TMath::DegToRad())*TMath::Sin(ph*TMath::DegToRad());
    pz = mom*TMath::Cos(th*TMath::DegToRad());
    E=sqrt(mom*mom + mass*mass);
    lv.SetPxPyPzE(px,py,pz,E);

    return lv;
}

void FillData(HParticleCand *cand, KFitParticle *outcand, double arr[],
              double mass)//function for filling atributes of KFitParticle
{
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand->SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                         std::cos(cand->getPhi() * deg2rad),
                     cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                         std::sin(cand->getPhi() * deg2rad),
                     cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                     mass);
    outcand->setMomentum(cand->getMomentum());
    outcand->setTheta(cand->getTheta());
    outcand->setPhi(cand->getPhi());
    outcand->setR(cand->getR());
    outcand->setZ(cand->getZ());
    outcand->setCovariance(cov);
    
}

void FillData(HForwardCand *cand, KFitParticle *outcand, double arr[],
              double mass)//function for filling atributes of KFitParticle
{
 
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand->SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                         std::cos(cand->getPhi() * deg2rad),
                     cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                         std::sin(cand->getPhi() * deg2rad),
                     cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                     mass);
    outcand->setMomentum(cand->getMomentum());
    outcand->setTheta(cand->getTheta());
    outcand->setPhi(cand->getPhi());
    outcand->setR(cand->getR());
    outcand->setZ(cand->getZ());
    outcand->setCovariance(cov);
}
    double deg2rad = TMath::DegToRad();
void FillData(HEmcNeutralCand *cand, KFitParticle *outcand, double arr[],double mass)//function for filling atributes of KFitParticle
{
    double deg2rad = TMath::DegToRad();

    TMatrixD cov(5, 5);
    cov(0, 0) = std::pow(arr[0], 2);
    cov(1, 1) = std::pow(arr[1], 2);
    cov(2, 2) = std::pow(arr[2], 2);
    cov(3, 3) = std::pow(arr[3], 2);
    cov(4, 4) = std::pow(arr[4], 2);

    outcand->SetXYZM(cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                         std::cos(cand->getPhi() * deg2rad),
                     cand->getMomentum() * std::sin(cand->getTheta() * deg2rad) *
                         std::sin(cand->getPhi() * deg2rad),
                     cand->getMomentum() * std::cos(cand->getTheta() * deg2rad),
                     mass);
    outcand->setMomentum(cand->getMomentum());
    outcand->setTheta(cand->getTheta());
    outcand->setPhi(cand->getPhi());
    outcand->setR(cand->getR());
    outcand->setZ(cand->getZ());
    outcand->setCovariance(cov);
}
Int_t ecal_tests(HLoop* loop, const AnaParameters& anapars)
{
    // Check if loop was properly initialized
    if (!loop->setInput(""))
    { // reading file structure
        std::cerr << "READBACK: ERROR : cannot read input !" << std::endl;
        std::exit(EXIT_FAILURE);
    }

    // Timer for checking analysis time
    TStopwatch timer;
    timer.Reset(); // Reset timer
    timer.Start(); // Start timer (T0)

    //----------------------------------------------------------------------------------------------
    // Parameters file opening
    //----------------------------------------------------------------------------------------------
     //ofstream out1, out2;
     //----------------------------------------------------------------------------------------------
    // Accessing categories inside input file
    //----------------------------------------------------------------------------------------------

    // Hades Particle Candidates
    HCategory* fParticleCand = HCategoryManager::getCategory(catParticleCand, kTRUE, "catParticleCand");
    if (!fParticleCand) { cout << "No catParticleCand!" << endl; }

    HCategory* fEmcNeutralCand = HCategoryManager::getCategory(catEmcNeutralCand, kTRUE, "catEmcNeutralCand");
    if(!fEmcNeutralCand){cout << "No catEmcNeutralCand!" << endl;}
    
    HCategory* fEmcCluster = HCategoryManager::getCategory(catEmcCluster, 0, "catEmcCluster");
    if(!fEmcCluster){cout << "No catEmcCluster!" << endl;}

    HCategory * fStart2Hit = HCategoryManager::getCategory(catStart2Hit, kTRUE, "catStart2Hit");
    if (!fStart2Hit) { cout << "No catStart2Hit!" << endl; }

    HEnergyLossCorrPar dEdxCorr;
    dEdxCorr.setDefaultPar("feb22");    

    
    //----------------------------------------------------------------------------------------------
    // Setting parameters for loop over events
    //----------------------------------------------------------------------------------------------

    Int_t entries = loop->getEntries(); // Number of entries in loop
    int limit_sta = anapars.start;      // Limit START - Where to start the loop
    int limit_sto = 0;                  // Limit STOP - Where to stop the loop

    if (anapars.events >= 0)
        limit_sto = limit_sta + anapars.events;
    else
        limit_sto = entries;

    if (limit_sto > entries) limit_sto = entries;

    //----------------------------------------------------------------------------------------------
    // Specifying output file
    //----------------------------------------------------------------------------------------------

    TFile* output_file = TFile::Open(anapars.outfile, "RECREATE");
    //TFile* output_file = new TFile("file.root","RECREATE");
    output_file->cd();
    cout << "NEW ROOT TREE " << endl;

    //----------------------------------------------------------------------------------------------
  
  TH1F * Eta_m=  new TH1F("Eta_m","Eta_m",200000,0,2000);
  TH1F * Eta_theta=  new TH1F("Eta_theta","Eta_theta",200000,0,90);
  TH1F * Eta_phi=  new TH1F("Eta_phi","Eta_phi",200000,-360,360);
  
  TH1F * Eta_m_before=  new TH1F("Eta_m_before","Eta_m_before",200000,0,2000);
  TH1F * Eta_theta_before=  new TH1F("Eta_theta_before","Eta_theta_before",200000,0,90);
  TH1F * Eta_phi_before=  new TH1F("Eta_phi_before","Eta_phi_before",200000,-360,360);

  TH1F * Eta_m_selection=  new TH1F("Eta_m_selection","Eta_m_selection",200000,0,2000);
  TH1F * Eta_theta_selection=  new TH1F("Eta_theta-selection","Eta_theta_selection",200000,0,90);
  TH1F * Eta_phi_selection=  new TH1F("Eta_phi_selection","Eta_phi_selection",200000,-360,360);
//

  TH1F * EtaP_m=  new TH1F("EtaP_m","EtaP_m",200000,0,2000);
  TH1F * EtaP_theta=  new TH1F("EtaP_theta","EtaP_theta",200000,0,90);
  TH1F * EtaP_phi=  new TH1F("EtaP_phi","EtaP_phi",200000,-360,360);
  
  
  //
    
    TH2F *hbeta_mom= new TH2F("hbeta_mom","hbeta_mom",1000,-4000.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_p= new TH2F("hbeta_mom_p","hbeta_mom_p",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pip= new TH2F("hbeta_mom_pip","hbeta_mom_pip",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pim= new TH2F("hbeta_mom_pim","hbeta_mom_pim",1000,-4000.,0.,700,0.,1.4);

    TH2F* hmass_mom = new TH2F("hmass_mom","hmass_mom; p*q [MeV/c]; mass",1000,-2000,4000,1000,0,2000);
    TH1F *hmass=new TH1F("hmass","hmass",1000,-2000,2000);
    
    TH1F *hbeta=new TH1F("hbeta","hbeta",100,0,2);
    TH1F *hbeta100=new TH1F("hbeta100","hbeta100",100,0,2);
    TH1F *hbeta150=new TH1F("hbeta150","hbeta150",100,0,2);
    TH1F *hbeta500=new TH1F("hbeta500","hbeta500",100,0,2);
    TH1F *hbeta100_1=new TH1F("hbeta100_1","hbeta100_1",100,0,2);
    TH1F *hbeta500_1=new TH1F("hbeta500_1","hbeta500_1",100,0,2);
    TH1F *htime=new TH1F("htime","htime",2200,-200,2000);
    TH1F *htime_PT2=new TH1F("htime_PT2","htime_PT2",2200,-200,2000);
    TH1F *htime_PT3=new TH1F("htime_PT3","htime_PT3",2200,-200,2000);
    
    TH1F *htracklength=new TH1F("htracklength","htracklength",2000,2000,4000);
    TH1F *hcl_size=new TH1F("hcl_size","hcl_size",20,0,20);

    TH1F *hnumOfNeutral=new TH1F("hnumOfNeutral","hnumOfNeutral",100,0,100);
    TH1F *hnumOfNeutral_g100=new TH1F("hnumOfNeutral_g100","hnumOfNeutral_g100",100,0,100);
    TH1F *hg_energy=new TH1F("hg_energy","hg_energy",2000,0,2000);
    TH1F *hg_energy_cl1=new TH1F("hg_energy_cl1","hg_energy_cl1",2000,0,2000);

    TH1F *hdphi_pp=new TH1F("hdphi_pp","hdphi_pp",360,0,360);
    TH2F *htantan_dphi_pp=new TH2F("htantan_dphi_pp","htantan_dphi_pp",360,0,360,100,0,1);

    TH1F *hMM_pp=new TH1F("hMM_pp","hMM_pp",1400,0,1400);
    TH1F *hMM_pp_pmomCut1=new TH1F("hMM_pp_pmomCut1","hMM_pp_pmomCut1",1400,0,1400);
    TH1F *hMM_pp_pmomCut2=new TH1F("hMM_pp_pmomCut2","hMM_pp_pmomCut2",1400,0,1400);
    TH1F *hMM_pp_pmomCut3=new TH1F("hMM_pp_pmomCut3","hMM_pp_pmomCut3",1400,0,1400);
    TH1F *hMM_pp_pmomCut4=new TH1F("hMM_pp_pmomCut4","hMM_pp_pmomCut4",1400,0,1400);
    TH1F *hMM2_pp=new TH1F("hMM2_pp","hMM2_pp",900,-1,2);
    TH1F *hMM2_pp_PT2=new TH1F("hMM2_pp_PT2","hMM2_pp_PT2",900,-1,2);
    TH1F *hMM2_pp_PT3=new TH1F("hMM2_pp_PT3","hMM2_pp_PT3",900,-1,2);
    TH1F *hMM2_pp_pmomCut1=new TH1F("hMM2_pp_pmomCut1","hMM2_pp_pmomCut1",900,-1,2);
    TH1F *hMM2_pp_pmomCut2=new TH1F("hMM2_pp_pmomCut2","hMM2_pp_pmomCut2",900,-1,2);
    TH1F *hMM2_pp_pmomCut3=new TH1F("hMM2_pp_pmomCut3","hMM2_pp_pmomCut3",900,-1,2);
    TH1F *hMM2_pp_pmomCut4=new TH1F("hMM2_pp_pmomCut4","hMM2_pp_pmomCut4",900,-1,2);


    TH1F *hMM2_pp_noEl=new TH1F("hMM2_pp_noEl","hMM2_pp_noEl",1500,-1,2);
    TH2F *htheta_mom= new TH2F("htheta_mom","htheta_mom",400,0,8,200,0,90);
    TH2F *htheta_mom_1= new TH2F("htheta_mom_1","htheta_mom_1",400,0,8,200,0,90);
    TH2F *htheta_mom_p_elast= new TH2F("htheta_mom_p_elast","htheta_mom_p_elast",200,0,5,200,0.0,90);

    TH1F *hmult_p=new TH1F("hmult_p","hmult_p",10,0,10);
    TH1F *hOA_gg=new TH1F("hOA_gg","hOA_gg",180,0,180);

    TH2F *henergyvscell = new TH2F("henergyvscell", "henergyvscell; sec*255+cell; energy [MeV]", 1530,0,1530,500,0,1500);
    TH2F *htimevscell = new TH2F("htimevscell", "htimevscell; sec*255+cell; time [ns]", 1530,0,1530,110,-10,100);
    TH2F *henergyvscell1 = new TH2F("henergyvscell1", "henergyvscell1; sec*255+cell; energy [MeV]", 1530,0,1530,500,0,1500);
    TH2F *hbetavscell1 = new TH2F("hbetavscell1", "hbetavscell1; sec*255+cell; beta", 1530,0,1530,200,0,2);


    TH1F *hinvM_gg_ETA_PT3=new TH1F("hinvM_gg_ETA_PT3","hinvM_gg_ETA_PT3",700,0,700);
    TH1F *hinvM_gg_omega_PT3=new TH1F("hinvM_gg_omega_PT3","hinvM_gg_omega_PT3",400,300,1500);
    TH1F *hinvM_gg_omega_ETA_PT3=new TH1F("hinvM_gg_omega_ETA_PT3","hinvM_gg_omega_ETA_PT3",700,0,700);

    TH1F *hinvM_gg_PT2=new TH1F("hinvM_gg_PT2","hinvM_gg_PT2",700,0,700);
    TH1F *hinvM_gg_PT3=new TH1F("hinvM_gg_PT3","hinvM_gg_PT3",700,0,700);
    //TH1F *hinvM_gg_oa2=new TH1F("hinvM_gg_oa2","hinvM_gg_oa2",700,0,700);
    TH1F *hinvM_ggMix=new TH1F("hinvM_ggMix","hinvM_ggMix",700,0,700);
    TH1F *hinvM_ggMixS=new TH1F("hinvM_ggMixS","hinvM_ggMixS",700,0,700);
    TH1F *hM3g=new TH1F("hM3g","hM3g; M_{3#gamma} [MeV/c^{2}]",1000,0,2000);
    TH1F* hM2g_3g=new TH1F("hM2g_3g","hM2g_3g; M_{#gamma#gamma} [MeV/c^{2}]",700,0,700);
    TH1F* hM3g_ETA=new TH1F("hM3g_ETA","hM3g_ETA; M_{3#gamma} [MeV/c^{2}]",1000,0,2000);
    TH1F* hM3g_ETAMixS=new TH1F("hM3g_ETAMixS","hM3g_ETAMixS; M_{3#gamma} [MeV/c^{2}]",1000,0,2000);

    TH1F* hMETAg_Mix=new TH1F("hMETAg_Mix","hMETAg_Mix; M_{#pi^{2}#gamma} [MeV/c^{2}]",1000,0,2000);
    
    TH1F *hMpETA_PT2=new TH1F("hMpETA_PT2","hMpETA_PT2; M_{p#pi^{0}} [MeV/c^{2}]",2000,1000,3000);
    TH1F *hMpETA_PT3=new TH1F("hMpETA_PT3","hMpETA_PT3; M_{p#pi^{0}} [MeV/c^{2}]",2000,1000,3000);
    TH1F *hMpETA_Mix=new TH1F("hMpETA_Mix","hMpETA_Mix; M_{p#pi^{0}} [MeV/c^{2}]",2000,1000,3000);
    TH1F *hMMpETA=new TH1F("hMMpETA","hMMpETA; MM_{p#pi^{0}} [MeV/c^{2}]",1000,0,3000);
    TH1F *hMM2pETA=new TH1F("hMM2pETA","hMM2pETA; MM^{2}_{p#pi^{0}} GeV^{2}/c^{4}]",1200,-6.,6);
    TH1F *hMM2pETA_Mix=new TH1F("hMM2pETA_Mix","hMM2pETA_Mix; MM^{2}_{p#pi^{0}} GeV^{2}/c^{4}]",1200,-6.,6);

    TH1F* hMpippimETA_PT2=new TH1F("hMpippimETA_PT2","hMpippimETA_PT2; M_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]",2000,0,2000);
    TH1F* hMpippimETA_PT3=new TH1F("hMpippimETA_PT3","hMpippimETA_PT3; M_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]",2000,0,2000);
    TH1F* hMggpippim_PT2=new TH1F("hMggpippim_PT2","hMggpippim_PT2; M_{#gamma#gamma#pi^{+}#pi^{-}} [MeV/c^{2}]",2000,0,2000);
    TH1F* hMggpippim_PT3=new TH1F("hMggpippim_PT3","hMggpippim_PT3; M_{#gamma#gamma#pi^{+}#pi^{-}} [MeV/c^{2}]",2000,0,2000);
    TH1F* hMpippimETA_Mix =new TH1F("hMpippimETA_Mix","hMpippimETA_Mix; M_{#pi^{+}#pi^{-}#pi^{0}} [MeV/c^{2}]",2000,0,2000);    TH1F* hMggpippim_Mix =new TH1F("hMggpippim_Mix","hMggpippim_Mix; M_{#gamma#gamma#pi^{+}#pi^{-}} [MeV/c^{2}]",2000,0,2000);   
    TH1F *hinvM_gg_pippim_PT3=new TH1F("hinvM_gg_pippim_PT3","hinvM_gg_pippim_PT3",700,0,700);

    TH1F * htheta_reco= new TH1F("htheta_reco","htheta_reco",1000,-180,180);
    TH1F * htheta_reco150= new TH1F("htheta_reco150","htheta_reco150",1000,-180,180);



//////////////////////////////////////////////////moje


    TH1F *hinvM_gg=new TH1F("hinvM_gg","hinvM_gg",700,0,700);
    char name[200];

/////////////////////////////////////////////////////////////HISTROGRAMS
//EMC gen
TH1F * hEPhtoton1EMC_gen= new TH1F("hEPhtoton1EMC_gen","E Phtoton1 EMC gen",100,0,3000);
TH1F * hEPhtoton2EMC_gen= new TH1F("hEPhtoton2EMC_gen","E Phtoton2 EMC gen",100,0,3000);

//PULLS
 TH1F * hPullEInvPhoton1Conv_mass_mix= new TH1F("Pull_E_photon_1_mix","Photon pull 1/E",5000,-100,100);
    TH1F * hPullThetaPhoton1Conv_mass_mix= new TH1F("Pull_theta_photon_1_mix","Photon pull #theta",5000,-100,100);
    TH1F * hPullPhiPhoton1Conv_mass_mix= new TH1F("Pull_phi_photon_1_mix","Photon pull #phi",5000,-100,100);
 


    TH1F * hPullEInvPhoton1Conv_mass= new TH1F("Pull_E_photon_1","Pull E photon 1",5000,-100,100);
    TH1F * hPullThetaPhoton1Conv_mass= new TH1F("Pull_theta_photon_1","Pull theta photon 1",5000,-100,100);
    TH1F * hPullPhiPhoton1Conv_mass= new TH1F("Pull_phi_photon_1","Pull phiphoton 1",5000,-100,100);
    TH1F * hPullRPhoton1Conv_mass= new TH1F("Pull_R_photon_1","Pull R photon 1",10000,-1,1);
    TH1F * hPullZPhoton1Conv_mass= new TH1F("Pull_Z_photon_1","Pull Z photon 1",10000,-1,1);

    TH1F * hPullEInvPhoton2Conv_mass= new TH1F("Pull_E_photon_2","Pull E photon 2",5000,-100,100);
    TH1F * hPullThetaPhoton2Conv_mass= new TH1F("Pull_theta_photon_2","Pull theta photon 2",5000,-100,100);
    TH1F * hPullPhiPhoton2Conv_mass= new TH1F("Pull_phi_photon_2","Pull phi photon 2",5000,-100,100);
    TH1F * hPullRPhoton2Conv_mass= new TH1F("Pull_R_photon_2","Pull R photon 2",10000,-1,1);
    TH1F * hPullZPhoton2Conv_mass= new TH1F("Pull_Z_photon_2","Pull Z photon 2",10000,-1,1);

    //photons before refit
    //photons before refit
  
    TH1F * hEPhtotonbeforeRefitcomb= new TH1F("Before Energy Photon comb","E Phtoton comb Before Refit",100,0,3000);
    TH1F * hEPhtoton1beforeRefit= new TH1F("Before Energy Phtoton1 ","E Phtoton1 Before Refit",100,0,3000);
    TH1F * hEPhtoton2beforeRefit= new TH1F("Before Energy Photon2","E Phtoton2 Before Refit",100,0,3000);

    //photons after refit
    TH1F * hEPhtotonAfterRefit_mix= new TH1F("hEPhtotonAfterRefit_mix","Refit Photon E",100,0,3000);
    TH1F * hPPhtotonAfterRefit_mix= new TH1F("hPPhtotonAfterRefit_mix","Refit Photon P",100,0,3000);
    TH1F * hmassPhtotonAfterRefit_mix= new TH1F("hmassPhtotonAfterRefit_mix","mass Phtoton Refit_mix",100,-100,100);
    TH1F * hthetaPhtotonAfterRefit_mix= new TH1F("hthetaphtotonAfterRefit_mix","Refit Photon #theta",1000,-180,180);
    TH1F * hphiPhtotonAfterRefit_mix= new TH1F("hphiPhtotonAfterRefit_mix","Refit Photon #phi",300,-180,180);


    TH1F * hEPhtoton1AfterRefit= new TH1F("hEPhtoton1AfterRefit","E Phtoton1 Refit",100,0,3000);
    TH1F * hPPhtoton1AfterRefit= new TH1F("hPPhtoton1AfterRefit","P Phtoton1 Refit",100,0,3000);
    TH1F * hmassPhtoton1AfterRefit= new TH1F("hmassPhtoton1AfterRefit","mass Phtoton1 Refit",100,-100,100);
    TH1F * hthetaPhtoton1AfterRefit= new TH1F("hthetahtoton1AfterRefit","theta Phtoton1 Refit",1000,-180,180);
    TH1F * hphiPhtoton1AfterRefit= new TH1F("hphiPhtoton1AfterRefit"," phi Phtoton1 Refit",300,-180,180);
    //REFIT

    TH1F * hEPhtoton2AfterRefit= new TH1F("hEPhtoton2AfterRefit","E Phtoton2 Refit",100,0,3000);
    TH1F * hPPhtoton2AfterRefit= new TH1F("hPPhtoton2AfterRefit","P Phtoton2 Refit",100,0,3000);
    TH1F * hmassPhtoton2AfterRefit= new TH1F("hmassPhtoton2AfterRefit","mass Phtoton2 Refit",100,-100,100);
    TH1F * hthetaPhtoton2AfterRefit= new TH1F("hthetaphtoton2AfterRefit"," theta Phtoton2 Refit",1000,-180,180);
    TH1F * hphiPhtoton2AfterRefit= new TH1F("hphiPhtoton2AfterRefit","phi Phtoton2 Refit",300,-180,180);

    TH1I * hIterations = new TH1I("hIterations","hIterations",10,0,10);
    TH1F * hChi2 = new TH1F("hChi2","hChi2",100,0,1);
    TH1F * hProb = new TH1F("hProb","hProb",100,0,1);


    TH1F * hEETAAfterRefit= new TH1F("hEETAAfterRefit","hEETAAfterRefit",1000,0,4000);
    TH1F * hPETAAfterRefit= new TH1F("hPETAAfterRefit","hPETAAfterRefit",1000,0,4000);
    TH1F * hmassETAAfterRefit= new TH1F("hmassETAAfterRefit","hmassETAAfterRefit",2000000,100,2000);
    TH1F * hthetaETAAfterRefit= new TH1F("hthetahtoton2AfterRefit","hthetaETAAfterRefit",1000,-180,180);
    TH1F * hphiETAAfterRefit= new TH1F("hphiETAAfterRefit","hEphiETAAfterRefit",1000,-180,180);


    //ETA simple- from g+g without refit
     
    TH1F * hEETABeforeRefit_selection= new TH1F("hEETAABeforeRefit","hEETABeforeRefit Selection",1000,0,4000);
    TH1F * hEETABeforeRefit= new TH1F("hEETABeforeRefit","hEETABeforeRefit",1000,0,4000);

    TH1F * hPETABeforeRefit_selection= new TH1F("hPETABeforeRefit_selection","hPETABeforeRefit_selection",100,0,4000);
    TH1F * hPETABeforeRefit= new TH1F("hPETABeforeRefit","hPETABeforeRefit",100,0,4000);
    TH1F * hmassETABeforeRefit_selection= new TH1F("hmassETABeforeRefit_selection","hmassETABeforeRefit_selection",2000000,100,2000);
    TH1F * hmassETABeforeRefit= new TH1F("hmassETABeforeRefit","hmassETABeforeRefit",2000000,0,2000);
    TH1F * hthetahtoton2BeforeRefit_selection= new TH1F("hthetahtoton2BeforeRefit_selection","hthetaETABeforeRefit_selection",1000,-180,180);
    TH1F * hthetahtoton2BeforeRefit= new TH1F("hthetahtoton2BeforeRefit","hthetaETABeforeRefit",300,-180,180);


    TH1F * hthetaETABeforeRefit_selection= new TH1F("hthetaETABeforeRefit_selection","hEphiETABeforeRefit_selection",1000,-180,180);
    TH1F * hthetaETABeforeRefit= new TH1F("hthetaETABeforeReft","hEphiETABeforeRefit",1000,-180,180);

    TH1F * hphiETABeforeRefit_selection= new TH1F("hphiETABeforeRefit_selection","hEphiETABeforeRefit_selection",1000,-180,180);
    TH1F * hphiETABeforeRefit= new TH1F("hphiETABeforeRefit","hEphiETABeforeRefit",1000,-180,180);
    //Residuals
    TH1F * hResidua_Refit_Reco_photon1_E= new TH1F("Residua Refit-Reco 1/E photon1","Residua Refit-Reco E 1/E",100,-1000,1000);
    TH1F * hResidua_Refit_Reco_photon1_theta= new TH1F("Residua Refit-Reco theta photon1","Residua Refit-Reco theta photon1",100,-1000,1000);
    TH1F * hResidua_Refit_Reco_photon1_phi= new TH1F("Residua Refit-Reco phi photon1","Residua Refit-Reco phi photon1",100,-1000,1000);
    TH1F * hResidua_Refit_Reco_photon1_R= new TH1F("Residua Refit-Reco R photon1","Residua Refit-Reco R photon1",100,-1000,1000);
    TH1F * hResidua_Refit_Reco_photon1_Z= new TH1F("Residua Refit-Reco Z photon1","Residua Refit-Reco Z photon1",100,-1000,1000);

     TH1F * hResidua_True_Reco_photon1_E= new TH1F("Residua True-Reco 1/E photon1","Residua True-Reco E 1/E",100,-1000,1000);
    TH1F * hResidua_True_Reco_photon1_theta= new TH1F("Residua True-Reco theta photon1","Residua True-Reco theta photon1",100,-1000,1000);
    TH1F * hResidua_True_Reco_photon1_phi= new TH1F("Residua True-Reco phi photon1","Residua True-Reco phi photon1",100,-1000,1000);
    TH1F * hResidua_True_Reco_photon1_R= new TH1F("Residua True-Reco R photon1","Residua True-Reco R photon1",100,-1000,1000);
    TH1F * hResidua_True_Reco_photon1_Z= new TH1F("Residua True-Reco Z photon1","Residua True-Reco Z photon1",100,-1000,1000);


    TH1F * hResidua_True_Refit_photon1_E= new TH1F("Residua True-Refit 1/E photon1","Residua True-Refit E 1/E",100,-1000,1000);
    TH1F * hResidua_True_Refit_photon1_theta= new TH1F("Residua True-Refit theta photon1","Residua True-Refit theta photon1",100,-1000,1000);
    TH1F * hResidua_True_Refit_photon1_phi= new TH1F("Residua True-Refit phi photon1","Residua True-Refit phi photon1",100,-1000,1000);
    TH1F * hResidua_True_Refit_photon1_R= new TH1F("Residua True-Refit R photon1","Residua True-Refit R photon1",100,-1000,1000);
    TH1F * hResidua_True_Refit_photon1_Z= new TH1F("Residua True-Refit Z photon1","Residua True-Refit Z photon1",100,-1000,1000);

  

    //Refit GEANT- performing refit on geant tracks
    //PULLS
    TH1F * hPullEInvPhoton1Conv_mass_geant= new TH1F("Pull_E_photon_1_GEANT","Pull E photon 1_GEANT",5000,-100,100);
    TH1F * hPullThetaPhoton1Conv_mass_geant= new TH1F("Pull_theta_photon_1_GEANT","Pull theta photon 1_GEANT",5000,-100,100);
    TH1F * hPullPhiPhoton1Conv_mass_geant= new TH1F("Pull_phi_photon_1_GEANT","Pull phiphoton 1_GEANT",5000,-100,100);
    TH1F * hPullRPhoton1Conv_mass_geant= new TH1F("Pull_R_photon_1_GEANT","Pull R photon 1_GEANT",10000,-1,1);
    TH1F * hPullZPhoton1Conv_mass_geant= new TH1F("Pull_Z_photon_1_GEANT","Pull Z photon 1_GEANT",10000,-1,1);

    TH1F * hPullEInvPhoton2Conv_mass_geant= new TH1F("Pull_E_photon_2_GEANT","Pull E photon 2_GEANT",5000,-100,100);
    TH1F * hPullThetaPhoton2Conv_mass_geant= new TH1F("Pull_theta_photon_2_GEANT","Pull theta photon 2_GEANT",5000,-100,100);
    TH1F * hPullPhiPhoton2Conv_mass_geant= new TH1F("Pull_phi_photon_2_GEANT","Pull phi photon 2_GEANT",5000,-100,100);
    TH1F * hPullRPhoton2Conv_mass_geant= new TH1F("Pull_R_photon_2_GEANT","Pull R photon 2_GEANT",10000,-1,1);
    TH1F * hPullZPhoton2Conv_mass_geant= new TH1F("Pull_Z_photon_2_GEANT","Pull Z photon 2_GEANT",10000,-1,1);

    //photons before Refit_GEANT
    TH1F * hEPhtoton1beforeRefit_GEANT= new TH1F("Before Energy Phtoton1 geant ","E Phtoton1 Before Refit_GEANT",100,0,3000);
    TH1F * hEPhtoton2beforeRefit_GEANT= new TH1F("Before Energy Photon2 geant","E Phtoton2 Before Refit_GEANT",100,0,3000);

    //photons after Refit_GEANT
    TH1F * hEPhtoton1AfterRefit_GEANT= new TH1F("hEPhtoton1AfterRefit_GEANT","E Phtoton1 Refit_GEANT",100,0,3000);
    TH1F * hPPhtoton1AfterRefit_GEANT= new TH1F("hPPhtoton1AfterRefit_GEANT","P Phtoton1 Refit_GEANT",100,0,3000);
    TH1F * hmassPhtoton1AfterRefit_GEANT= new TH1F("hmassPhtoton1AfterRefit_GEANT","mass Phtoton1 Refit_GEANT",100,-100,100);
    TH1F * hthetaPhtoton1AfterRefit_GEANT= new TH1F("hthetahtoton1AfterRefit_GEANT","theta Phtoton1 Refit_GEANT",1000,-180,180);
    TH1F * hphiPhtoton1AfterRefit_GEANT= new TH1F("hphiPhtoton1AfterRefit_GEANT"," phi Phtoton1 Refit_GEANT",300,-180,180);
    //Refit_GEANT

    TH1F * hEPhtoton2AfterRefit_GEANT= new TH1F("hEPhtoton2AfterRefit_GEANT","E Phtoton2 Refit_GEANT",100,0,3000);
    TH1F * hPPhtoton2AfterRefit_GEANT= new TH1F("hPPhtoton2AfterRefit_GEANT","P Phtoton2 Refit_GEANT",100,0,3000);
    TH1F * hmassPhtoton2AfterRefit_GEANT= new TH1F("hmassPhtoton2AfterRefit_GEANT","mass Phtoton2 Refit_GEANT",100,-100,100);
    TH1F * hthetaPhtoton2AfterRefit_GEANT= new TH1F("hthetaphtoton2AfterRefit_GEANT"," theta Phtoton2 Refit_GEANT",1000,-180,180);
    TH1F * hphiPhtoton2AfterRefit_GEANT= new TH1F("hphiPhtoton2AfterRefit_GEANT","phi Phtoton2 Refit_GEANT",300,-180,180);

    TH1I * hIterations_GEANT = new TH1I("hIterations_geant","hIterations_geant",10,0,10);
    TH1F * hChi2_GEANT = new TH1F("hChi2_GEANT","hChi2_GEANT",100,0,1);
    TH1F * hProb_GEANT = new TH1F("hProb_GEANT","hProb_GEANT",100,0,1);


    TH1F * hEETAAfterRefit_GEANT= new TH1F("hEETAAfterRefit_GEANT","hEETAAfterRefit_GEANT",100,0,4000);
    TH1F * hPETAAfterRefit_GEANT= new TH1F("hPETAAfterRefit_GEANT","hPETAAfterRefit_GEANT",100,0,4000);
    TH1F * hmassETAAfterRefit_GEANT= new TH1F("hmassETAAfterRefit_GEANT","hmassETAAfterRefit_GEANT",2000000,100,2000);
    TH1F * hthetaETAAfterRefit_GEANT= new TH1F("hthetahtoton2AfterRefit_GEANT","hthetaETAAfterRefit_GEANT",1000,-180,180);
    TH1F * hphiETAAfterRefit_GEANT= new TH1F("hphiETAAfterRefit_GEANT","hEphiETAAfterRefit_GEANT",300,-180,180);

//////////////////////////////////////////////////moje
   // char name[200];
    TH1F *hbeta_sec[6];
    
    for (int j=0;j<6;j++){

      sprintf(name,"hbeta_sec%d",j);  
      hbeta_sec[j]= new TH1F(name,name,100,0,2);

    }
    
    //*********************************************************
    //elastic cut - tight
    TCutG *cutg_th_mom_p = new TCutG("cutg_th_mom_p",10);
    cutg_th_mom_p->SetPoint(0,980.4627,52.15846);
    cutg_th_mom_p->SetPoint(1,2083.29,34.22742);
    cutg_th_mom_p->SetPoint(2,3332.648,19.71182);
    cutg_th_mom_p->SetPoint(3,4011.311,12.13383);
    cutg_th_mom_p->SetPoint(4,4959.897,13.30788);
    cutg_th_mom_p->SetPoint(5,4057.584,22.5936);
    cutg_th_mom_p->SetPoint(6,2854.499,34.22742);
    cutg_th_mom_p->SetPoint(7,1828.792,51.09113);
    cutg_th_mom_p->SetPoint(8,965.0385,52.05172);
    cutg_th_mom_p->SetPoint(9,980.4627,52.15846);

    //elast cut loose
    TCutG *cutg = new TCutG("cutg",6);
    cutg->SetPoint(0,1.232012,45.9507);
    cutg->SetPoint(1,3.955947,11.72535);
    cutg->SetPoint(2,6.474303,13.15141);
    cutg->SetPoint(3,3.544787,51.02113);
    cutg->SetPoint(4,1.254038,45.9507);
    cutg->SetPoint(5,1.232012,45.9507);
    //****************************************
    //PROT PID - mass vs mom

    TCutG *cutP = new TCutG("cutP",11);
    cutP->SetPoint(0,188.0792,1187.5);
    cutP->SetPoint(1,1559.183,1378.472);
    cutP->SetPoint(2,2203.665,1517.361);
    cutP->SetPoint(3,2677.548,1579.861);
    cutP->SetPoint(4,2886.057,1392.361);
    cutP->SetPoint(5,2544.861,1177.083);
    cutP->SetPoint(6,2052.022,888.8889);
    cutP->SetPoint(7,188.0792,593.75);
    cutP->SetPoint(8,86.98397,642.3611);
    cutP->SetPoint(9,200.7161,1190.972);
    cutP->SetPoint(10,188.0792,1187.5);

    //****************************************
    //PI+ PID - mass vs mom
    TCutG *cutPIP = new TCutG("cutPIP",10);
     cutPIP->SetPoint(0,42.75482,218.75);
     cutPIP->SetPoint(1,207.0345,520.8333);
     cutPIP->SetPoint(2,1489.68,791.6667);
     cutPIP->SetPoint(3,1913.016,760.4167);
     cutPIP->SetPoint(4,1495.998,364.5833);
     cutPIP->SetPoint(5,794.6504,27.77775);
     cutPIP->SetPoint(6,598.7784,-20.83336);
     cutPIP->SetPoint(7,30.11792,10.41664);
     cutPIP->SetPoint(8,49.07327,225.6944);
     cutPIP->SetPoint(9,42.75482,218.75);
    //****************************************
    //PI- PID - mass vs mom
     TCutG *cutPIM = new TCutG("cutPIM",11);
     cutPIM->SetPoint(0,-52.02194,208.3333);
     cutPIM->SetPoint(1,-816.5544,684.0278);
     cutPIM->SetPoint(2,-1543.176,836.8056);
     cutPIM->SetPoint(3,-1783.277,701.3889);
     cutPIM->SetPoint(4,-1505.265,395.8333);
     cutPIM->SetPoint(5,-1189.343,170.1389);
     cutPIM->SetPoint(6,-645.9562,3.472193);
     cutPIM->SetPoint(7,-64.65884,-6.944474);
     cutPIM->SetPoint(8,-52.02194,229.1666);
     cutPIM->SetPoint(9,-70.97729,215.2778);
     cutPIM->SetPoint(10,-52.02194,208.3333);
    //*******************************************************
    //Float_t Chi2RkCut      = 500.;
    //Float_t Chi2InCut      = 500.;
    //Float_t Chi2OutCut     = 500.;

    const float oAngleCut=6.;

    
    vector<TLorentzVector> lv_neutr, lv_neutr1, lv_prot, lv_pip, lv_pim ;
    vector<HEmcNeutralCand*> gamma_vector;//vector of pointers to selected gammas
    vector<TLorentzVector> lv_prot1, lv_prot2, lv_prot3, lv_prot4 ;
    vector<TLorentzVector> lv_gMix, lv_ETA;
    vector<HParticleCand*> part_neg, part_pos;
    vector<void*> vgTracks, vpipTracks, vpimTracks;

    //lv_gMix.clear();
    /*
    HGenericEventMixerObj < TLorentzVector > eventmixer;
    eventmixer.setPIDs(1, 1, 7);
    eventmixer.setBuffSize(30);


    HGenericEventMixerObj < TLorentzVector > eventmixer1;
    eventmixer1.setPIDs(1, 7, 24);
    eventmixer1.setBuffSize(30);
    */
  

    
    const double mp=938.27231;
    const double beam_energy = 4536. + mp;
    const double beam_momentum = sqrt(beam_energy*beam_energy-mp*mp);

    TLorentzVector proj(0,0,beam_momentum, beam_energy);
    TLorentzVector targ(0,0,0, mp);
    TLorentzVector beam(0,0,0,0);
    beam = proj + targ;
    
    int tbit_PT2 = 0;
    int tbit_PT3 = 0;
    double ult_z=0,ult_r=0;

std::vector<KFitParticle *> kFit_gamma;
std::vector<KFitParticle> candsFit;
    for (Int_t ev = limit_sta; ev < limit_sto; ev++) // event loop
    {

        if (ev % 10000 == 0)
        {
            printf("Event nr.: %d, progress: %.2f%%\n", ev, (double)(ev - limit_sta) / (limit_sto - limit_sta) * 100.);
        }

        /*Int_t nbytes =*/loop->nextEvent(ev); // get next event. categories will be cleared before

	lv_neutr.clear();
	lv_neutr1.clear();
	part_neg.clear();
	part_pos.clear();
	lv_prot1.clear();
	lv_prot2.clear();
	lv_prot3.clear();
	lv_prot4.clear();
	lv_prot.clear();
	lv_pip.clear();
	lv_pim.clear();
	lv_ETA.clear();
	
	lv_gMix.clear();
	vgTracks.clear();
	vpipTracks.clear();
	vpimTracks.clear();

  kFit_gamma.clear();
  gamma_vector.clear();
  //clearing
	
	HStart2Hit * fstart = nullptr;
	fstart = (HStart2Hit *) fStart2Hit->getObject(0);
	if (!fstart || fstart->getCorrFlag() == -1) continue;
	//fCorrFlag=-1 iTOF
	//fCorrFlag=0 only LGAD
	//fCorrFlag>=0
	
	
        HEventHeader* event_header = NULL;
        if (!(event_header = gHades->getCurrentEvent()->getHeader())) continue;

        Int_t TBit = (Int_t)event_header->getTBit();
        Double_t VertexX = event_header->getVertexReco().getX();
        Double_t VertexY = event_header->getVertexReco().getY();
        Double_t VertexZ = event_header->getVertexReco().getZ();
	if(VertexZ<-200 || VertexZ>-0) continue;


	//         Int_t DF     = (Int_t) event_header->getDownscalingFlag();
        //         Int_t SeqNum = (Int_t) event_header->getEventSeqNumber();
        //         Int_t TDec   = (Int_t) event_header->getTriggerDecision();


	int trigbit=-1;

	if((TBit&4096)==4096){

	  trigbit=0;
	  tbit_PT2++;
	}
	if((TBit&8192)==8192){

	  trigbit=1;
	  tbit_PT3++;


	}


	
	//****************************************************
        
	Int_t nNeutral_ev = fEmcNeutralCand->getEntries();

	hnumOfNeutral->Fill(nNeutral_ev);
	
	for (int j = 0; j < nNeutral_ev; ++j)
	  {
	  
	    HEmcNeutralCand* neutr_cand = HCategoryManager::getObject(neutr_cand, fEmcNeutralCand, j);

		Float_t dist  = neutr_cand->getDistanceToEmc();
		Int_t ind=neutr_cand->getEmcClusterIndex();

		//cout<<j<<" "<<ind<<endl;
		//if(ind<0)continue;

		HEmcCluster *cl=nullptr;
		cl=HCategoryManager::getObject(cl, fEmcCluster, ind);
		Int_t cl_size = cl->getNCells();
		Int_t pid=neutr_cand->getPID();

		Int_t sec = cl->getSector();

		Int_t cel = cl->getCell();		
		if(cel<33) continue;
		

		Double_t energy  = cl->getEnergy();
		HGeomVector trackVec(cl->getXLab()- VertexX, cl->getYLab()- VertexY, cl->getZLab() - VertexZ);
		
		Double_t trackLength = trackVec.length();
		trackVec  *= (energy/trackLength);

		Double_t tof  =cl->getTime();
		//Double_t beta = (trackLength/1000.) / (tof * 1.e-9 * TMath::C());
		Float_t theta = cl->getTheta();
		Float_t phi = cl->getPhi();
		Double_t beta = neutr_cand->getBeta();

		  
		hcl_size->Fill(cl_size);
		hbeta->Fill(beta);
		htime->Fill(tof);

		if(trigbit==0)htime_PT2->Fill(tof);
		if(trigbit==1)htime_PT3->Fill(tof);
		htracklength->Fill(trackLength);
		//cout<<tof<<" "<<trackLength<<endl;
		henergyvscell->Fill(sec*200+cel,energy);
		htimevscell->Fill(sec*200+cel,tof);
		
		if(energy>100 && cl_size==1)hbeta100->Fill(beta);
		if(energy>500 && cl_size==1)hbeta500->Fill(beta);
		
		if(energy>100 && cl_size==1 && tof>0 && tof<50){

		  hbeta100_1->Fill(beta);
		  henergyvscell1->Fill(sec*200+cel,energy);
		  hbetavscell1->Fill(sec*200+cel,beta);
		}
		if(energy>500 && cl_size==1 && tof>0 && tof<50){
		  hbeta500_1->Fill(beta);
		  hbeta_sec[sec]->Fill(beta);
		}
		hg_energy->Fill(energy);
		if(cl_size==1)hg_energy_cl1->Fill(energy);


		

		TLorentzVector lvg1, lvg;
		lvg1.SetXYZM(trackVec.getX(),trackVec.getY(),trackVec.getZ(),0);

		TLorentzVector *lvga=new TLorentzVector();
		lvga->SetXYZM(trackVec.getX(),trackVec.getY(),trackVec.getZ(),0);		

		//if (energy>100 && cl_size==1 && beta>0.7 && beta<1.3){
		//if (beta>0.7 && beta<1.3){
		//if (energy>100 && beta>0.8 && beta<1.2){
      KFitParticle* candidate = new KFitParticle(&lvg1,neutr_cand->getR(),neutr_cand->getZ());  //puste? 
      htheta_reco->Fill(candidate->getTheta()*TMath::RadToDeg());
 
		//if (energy>150 && beta>0.8 && beta<1.2){
      if(energy>100&&pid==1){
       hbeta150->Fill(beta);
                double errors[] = {0.055/(lvg1.E()*std::sqrt(lvg1.E())), 2.5*deg2rad,
                                    2.5*deg2rad,0.00001,0.00001};// need to be adjusted for gammmas /std::pow(lvg.E(),2)
                    FillData(neutr_cand, candidate, errors,0);//filling atributes of kfitparticles
                 //   cout<<lvg1.E()<<endl;
                kFit_gamma.push_back(candidate);
                gamma_vector.push_back(neutr_cand);
                hEPhtotonbeforeRefitcomb->Fill(candidate->Energy());
                htheta_reco150->Fill(candidate->getTheta()*TMath::RadToDeg());

		  lv_neutr1.push_back(lvg1);	    
		  lv_gMix.push_back(lvg1);
		  vgTracks.push_back(new TLorentzVector(*lvga));
		}



	  }
	//*********************************************************

	if (fParticleCand)
	      {
		Int_t nPart_ev  = fParticleCand->getEntries();
		//cout<<"---------->>> "<<endl;
		for (int j = 0; j < nPart_ev; ++j)
		  {
		    HParticleCand* fparticlecand = HCategoryManager::getObject(fparticlecand, fParticleCand, j);
		   
		    Float_t theta = fparticlecand->getTheta();
		    Float_t phi = fparticlecand->getPhi();
		    Float_t mom = fparticlecand->getMomentum();
		    Float_t beta = fparticlecand->getBeta();
		    Float_t tof = fparticlecand->getTof();
		    Int_t sec = fparticlecand->getSector();
		    //Float_t d2meta=fparticlecand->getDistanceToMetaHit();
		    Float_t charge=fparticlecand->getCharge(); 
		    Float_t mass= fparticlecand->getMass();

		    if(!fparticlecand->isFlagBit(kIsUsed))continue;


			hbeta_mom->Fill(charge*mom,beta);
			hmass_mom->Fill(charge*mom,mass);
			hmass->Fill(charge*mass);


					
			//if(charge>0 && mass > 650 && mass<1600){

			if(charge>0 &&  cutP->IsInside(charge*mom,mass) && beta < 1.){//PID protony
			  hbeta_mom_p->Fill(charge*mom,beta);
			  Double_t deltaMom = dEdxCorr.getDeltaMom(14, mom, theta); 
			  fparticlecand->setMomentum(mom + deltaMom);                                     
			  fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(14));

			  
			  part_pos.push_back(fparticlecand);
			  TLorentzVector lv1;
			  lv1=TLV_p(theta, phi, mom, 938.);

			  lv_prot.push_back(lv1);
			}			

			//if(charge>0 && mass < 650){//PID pi+
			if(charge>0 &&  cutPIP->IsInside(charge*mom,mass) && beta < 1.0){//PID pi+
			  
			  TLorentzVector lv1a;
			  lv1a=TLV_p(theta, phi, mom, 140.);
			  Double_t deltaMom = dEdxCorr.getDeltaMom(8, mom, theta); 
			  fparticlecand->setMomentum(mom + deltaMom);                                     
			  fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(8));

			  
			  lv_pip.push_back(lv1a);
			  hbeta_mom_pip->Fill(charge*mom,beta);
			  vpipTracks.push_back(new HParticleCand(*fparticlecand)); // Copy constructor required!
	  
			}			

			//if(charge<0 && mass < 550){//PID pi-
			if(charge<0 && cutPIM->IsInside(charge*mom,mass) && beta < 1.0){//PID pi-
			  
			  TLorentzVector lv1b;
			  lv1b=TLV_p(theta, phi, mom, 140.);

			  Double_t deltaMom = dEdxCorr.getDeltaMom(9, mom, theta); 
			  fparticlecand->setMomentum(mom + deltaMom);                                     
			  fparticlecand->calc4vectorProperties(HPhysicsConstants::mass(9));


			  
			  lv_pim.push_back(lv1b);
			  hbeta_mom_pim->Fill(charge*mom,beta);
			  vpimTracks.push_back(new HParticleCand(*fparticlecand));
			}			
			
			//******************************************	
			//******************************************	
		      
		  } 
    
	      }//fParticleCand

	//*****************************************************
	//**************   gg  *****************

	hnumOfNeutral_g100->Fill(lv_neutr1.size());
	
	    if(lv_neutr1.size()>=2){

	    
	      for (int ii=0;ii<lv_neutr1.size();ii++){
		for (int jj=0;jj<lv_neutr1.size();jj++){
		  if(ii==jj || ii>jj)continue;

		  float oAngle1 = lv_neutr1[ii].Angle(lv_neutr1[jj].Vect())*TMath::RadToDeg();
		  hOA_gg->Fill(oAngle1);

		  TLorentzVector lvg2;
		  lvg2=lv_neutr1[ii]+lv_neutr1[jj];
		  float mass_gg=lvg2.M();
		  if (oAngle1<oAngleCut) continue;
		  
		  //************** ETA ************************
		   if (trigbit==0) hinvM_gg_PT2->Fill(mass_gg);
		   if (trigbit==1) hinvM_gg_PT3->Fill(mass_gg);
		   
		      lv_ETA.push_back(lvg2);

		     //************* p-ETA  ***********************

		     if(lv_prot.size()==2){
		       TLorentzVector lv3, lv3a;
		       for (int i=0;i<lv_prot.size();i++){
			 lv3=lvg2+lv_prot[i];
			 lv3a=beam-lvg2-lv_prot[i];
			 //cout<<i<<" "<<lv3a.M()<<" "<<lvg2.M()<<" "<<lv_prot[i].P()<<endl;
			 
			 if (trigbit==0) hMpETA_PT2->Fill(lv3.M());
			 if (trigbit==1){
			   hMpETA_PT3->Fill(lv3.M());
			   hMMpETA->Fill(lv3a.M());
			   hMM2pETA->Fill(lv3a.M2()*1e-6);

			 }
		       }
		     }
		   //**************** pip-pim-ETA  ********************
		     if(lv_pip.size() && lv_pim.size()){
		       for (int i=0;i<lv_pip.size();i++){
			 for (int j=0;j<lv_pim.size();j++){
			   TLorentzVector lv3b=lvg2+lv_pip[i]+lv_pim[j];
			   float Mgg=lvg2.M();
			   
			   if (trigbit==1) {
			     hinvM_gg_pippim_PT3->Fill(Mgg);
			     hMggpippim_PT3->Fill(lv3b.M());
			   }
			   if (trigbit==0) hMpippimETA_PT2->Fill(lv3b.M());
			   if (trigbit==1 && Mgg>100 && Mgg<200) hMpippimETA_PT3->Fill(lv3b.M());
		      
			 }
		       }
		     }
		     //**********************************************************



		}
	      }
	    }

	    //**********************************************************
	    //*******************************
	    //********************3 gammas ************************************
	    	
	    if(lv_neutr1.size()>=3){
	    
	      TLorentzVector lvg3;
		 lvg3=lv_neutr1[0]+lv_neutr1[1]+lv_neutr1[2];
		 hM3g->Fill(lvg3.M());
		  
		  for (int ii=0;ii<lv_neutr1.size();ii++){
		    for (int jj=0;jj<lv_neutr1.size();jj++){
		      if(ii==jj || ii>jj)continue;

		      TLorentzVector lvg2;
		      lvg2=lv_neutr1[ii]+lv_neutr1[jj];

		      float mass_gg=lvg2.M();

		      hM2g_3g->Fill(mass_gg);

		      if(mass_gg>100 && mass_gg<200){
    
			for (int k=0;k<lv_neutr1.size();k++){
			  if( (k!=ii) && (k!=jj)) {
			    TLorentzVector lv4, lv4a;
			    
			    lv4=lvg2+lv_neutr1[k];
			    lv4a=beam-lvg2-lv_neutr1[k];
			    hM3g_ETA->Fill(lv4.M());

			  }
			}

		      }
		    }
		  }
	    }
	    
	 

	    

	    
 
	    hnumOfNeutral_g100->Fill(lv_neutr1.size());
        int gamma_num=kFit_gamma.size();
        if(gamma_num>=2)//at least 2 photons  were selected
        {
             //vector for photon candidates to refit
            for (int ii=0;ii<gamma_num;ii++)
            {
                for (int jj=ii+1;jj<gamma_num;jj++)//loop over second photon
                {
                    
                    //IZA analysis
                    float oAngle1 = lv_neutr1[ii].Angle(lv_neutr1[jj].Vect())*TMath::RadToDeg();
                    hOA_gg->Fill(oAngle1);

                    TLorentzVector lvg2;
                    lvg2=lv_neutr1[ii]+lv_neutr1[jj];
                    float mass_gg=lvg2.M();
                    hEETABeforeRefit->Fill(lvg2.E());
                    if (oAngle1<oAngleCut) continue;
                    
                    //************** ETA ************************
                    hinvM_gg->Fill(mass_gg);
                   // KFitParticle cand2Dec = *kFit_gamma[ii];
                    
                    KFitParticle *cand1DecPtr = kFit_gamma[ii];
                    KFitParticle cand1Dec = *cand1DecPtr;
                    candsFit.clear();
                    candsFit.push_back(cand1Dec);
                    KFitParticle *cand11DecPtr = kFit_gamma[jj];
                    KFitParticle  cand11Dec = *cand11DecPtr;
                    hEPhtoton1beforeRefit->Fill(cand1Dec.Energy());//energy before refit//gamma_vector[ii]->E());//
                    hEPhtoton2beforeRefit->Fill(cand11Dec.Energy());//gamma_vector[jj]->E());//
                    

                    //// emc gen
                    //double E1=emc_gen_kine[ii]->getTotalMomentum();
                   //double E2=emc_gen_kine[jj]->getTotalMomentum();
                   // hEPhtoton1EMC_gen->Fill(E1);
                    //hEPhtoton2EMC_gen->Fill(E2);

                    ///
                    candsFit.push_back(cand11Dec);//filling vector with second photon cadidate
                    KinFitter FitterMass(candsFit);// fitter declaration
                    FitterMass.setVerbosity(0);//fitter does not print any information
                    FitterMass.addMassConstraint(547.862);//mass constraint
                    FitterMass.setNumberOfIterations(30);
                    FitterMass.setConvergenceCriteria(1,9999.9,9999.9);//1 is for chi2 further ones are constraints for parameters better not to touch
                    //FitterMass.SetVerbosity(2);//fitter prints a lot of infomation
                    FitterMass.fit();
                    TLorentzVector ETA_simple=*kFit_gamma[ii]+*kFit_gamma[jj];
                    hphiETABeforeRefit->Fill(ETA_simple.Phi()*TMath::RadToDeg());
                    hthetaETABeforeRefit->Fill(ETA_simple.Theta()*TMath::RadToDeg());
                     total++;


                      Eta_m_before->Fill(ETA_simple.M());
                      Eta_theta_before->Fill(ETA_simple.Theta()*TMath::RadToDeg());
                      Eta_phi_before->Fill(ETA_simple.Phi()*TMath::RadToDeg());
                     /*  for(int ii_pim=0;ii_pim<lv_pim.size();ii_pim++)
                      {
                        for(int ii_pip=0;ii_pip<lv_pip.size();ii_pip++)
                        {
                          TLorentzVector Eta_before=ETA_simple+lv_pim[ii_pim]+lv_pip[ii_pip];
                          //cout<<EtaP.M()<<"  "<<EtaP.Theta()*TMath::RadToDeg()<<"  "<<EtaP.Phi()*TMath::RadToDeg()<<endl;
                          Eta_m_before->Fill(Eta_before.M());
                          Eta_theta_before->Fill(Eta_before.Theta()*TMath::RadToDeg());
                          Eta_phi_before->Fill(Eta_before.Phi()*TMath::RadToDeg());
                        }
                      }
                      */
                    if(FitterMass.isConverged())
                    {
                        hIterations->Fill(FitterMass.getIteration());
                        hChi2->Fill(FitterMass.getChi2());
                        hProb->Fill(FitterMass.getProb());
                            
                        // Get the pull distributions

                        hPullEInvPhoton1Conv_mass_mix->Fill(FitterMass.getPull(0));
                        hPullThetaPhoton1Conv_mass_mix->Fill(FitterMass.getPull(1));
                        hPullPhiPhoton1Conv_mass_mix->Fill(FitterMass.getPull(2));
                        hPullEInvPhoton1Conv_mass_mix->Fill(FitterMass.getPull(5));
                        hPullThetaPhoton1Conv_mass_mix->Fill(FitterMass.getPull(6));
                        hPullPhiPhoton1Conv_mass_mix->Fill(FitterMass.getPull(7));

                        hPullEInvPhoton1Conv_mass->Fill(FitterMass.getPull(0));
                        hPullThetaPhoton1Conv_mass->Fill(FitterMass.getPull(1));
                        hPullPhiPhoton1Conv_mass->Fill(FitterMass.getPull(2));
                        hPullRPhoton1Conv_mass->Fill(FitterMass.getPull(3));
                        hPullZPhoton1Conv_mass->Fill(FitterMass.getPull(4));

                        hPullEInvPhoton2Conv_mass->Fill(FitterMass.getPull(5));
                        hPullThetaPhoton2Conv_mass->Fill(FitterMass.getPull(6));
                        hPullPhiPhoton2Conv_mass->Fill(FitterMass.getPull(7));
                        hPullRPhoton2Conv_mass->Fill(FitterMass.getPull(8));
                        hPullZPhoton2Conv_mass->Fill(FitterMass.getPull(9));
                        



                        if (FitterMass.getProb() > 0.25)
                        {
                          refited++;


                            KFitParticle cand1mass = FitterMass.getDaughter(0); 
                            KFitParticle cand2mass = FitterMass.getDaughter(1);
                        TLorentzVector ETA_simple_new=cand1mass+cand2mass;
                          hthetaETABeforeRefit_selection->Fill(ETA_simple.Theta()*TMath::RadToDeg());
                                

                            //ETA
                           // cout<<lv_pim.size()<<"  "<<lv_pim.size()<<endl;
                    
                                for(int ii_pim=0;ii_pim<lv_pim.size();ii_pim++)
                                {
                                  for(int ii_pip=0;ii_pip<lv_pip.size();ii_pip++)
                                  {
                                    TLorentzVector EtaP=cand1mass+cand2mass+lv_pim[ii_pim]+lv_pip[ii_pip];
                                    //cout<<EtaP.M()<<"  "<<EtaP.Theta()*TMath::RadToDeg()<<"  "<<EtaP.Phi()*TMath::RadToDeg()<<endl;
                                    EtaP_m->Fill(EtaP.M());
                                    EtaP_theta->Fill(EtaP.Theta()*TMath::RadToDeg());
                                    EtaP_phi->Fill(EtaP.Phi()*TMath::RadToDeg());
                                  }
                                }

                          

                            hEPhtotonAfterRefit_mix->Fill(cand1mass.E());
                            hPPhtotonAfterRefit_mix->Fill(cand1mass.getMomentum());
                            hmassPhtotonAfterRefit_mix->Fill(std::pow(cand1mass.E(),2)-std::pow(cand1mass.getMomentum(),2));
                            hthetaPhtotonAfterRefit_mix->Fill(cand1mass.getTheta()*180.0/M_PI);
                            hphiPhtotonAfterRefit_mix->Fill(cand1mass.Phi()*180.0/M_PI);
                            hEPhtotonAfterRefit_mix->Fill(cand2mass.E());
                            hPPhtotonAfterRefit_mix->Fill(cand2mass.getMomentum());
                            hmassPhtotonAfterRefit_mix->Fill(std::pow(cand2mass.E(),2)-std::pow(cand2mass.getMomentum(),2));
                            hthetaPhtotonAfterRefit_mix->Fill(cand2mass.getTheta()*180.0/M_PI);
                            hphiPhtotonAfterRefit_mix->Fill(cand2mass.Phi()*180.0/M_PI);


                            hEPhtoton1AfterRefit->Fill(cand1mass.E());
                            hPPhtoton1AfterRefit->Fill(cand1mass.getMomentum());
                            hmassPhtoton1AfterRefit->Fill(std::pow(cand1mass.E(),2)-std::pow(cand1mass.getMomentum(),2));
                            hthetaPhtoton1AfterRefit->Fill(cand1mass.getTheta()*180.0/M_PI);
                            hphiPhtoton1AfterRefit->Fill(cand1mass.Phi()*180.0/M_PI);

                            hEPhtoton2AfterRefit->Fill(cand2mass.E());
                            hPPhtoton2AfterRefit->Fill(cand2mass.getMomentum());
                            hmassPhtoton2AfterRefit->Fill(std::pow(cand2mass.E(),2)-std::pow(cand2mass.getMomentum(),2));
                            hthetaPhtoton2AfterRefit->Fill(cand2mass.getTheta()*180.0/M_PI);
                            hphiPhtoton2AfterRefit->Fill(cand2mass.Phi()*180.0/M_PI);


                            //photon residuals
                            hResidua_Refit_Reco_photon1_E->Fill(cand1mass.Energy()-cand1Dec.Energy());
                            hResidua_Refit_Reco_photon1_theta->Fill(cand1mass.getTheta()-cand1Dec.getTheta());
                            hResidua_Refit_Reco_photon1_phi->Fill(cand1mass.getPhi()-cand1Dec.getPhi());
                            hResidua_Refit_Reco_photon1_R->Fill(cand1mass.getR()-cand1Dec.getR());
                            hResidua_Refit_Reco_photon1_Z->Fill(cand1mass.getZ()-cand1Dec.getZ());

                            //hResidua_True_Refit_photon1_E->Fill(geantkine_vec[ii]->getTotalMomentum()-cand1mass.Energy());
                            //cout<<geantkine_vec[ii]->getTotalMomentum()<<"    "<<cand1mass.Energy()<<endl;
                            //hResidua_True_Refit_photon1_theta->Fill(geantkine_vec[ii]->getTheta()-cand1mass.Theta());
                            //hResidua_True_Refit_photon1_phi
                            //hResidua_True_Refit_photon1_R
                            //hResidua_True_Refit_photon1_Z



                            KFitParticle ETA_refit = FitterMass.getMother(); //obtaining ETA mother partivle from refit
                            
                            
                            
                            hEETABeforeRefit_selection->Fill(ETA_simple.E());
                            
                            hPETABeforeRefit_selection->Fill(ETA_simple.P());
                            //hPETABeforeRefit
                            hmassETABeforeRefit_selection->Fill(ETA_simple.Mag());
                            //hmassETABeforeRefit
                            hthetahtoton2BeforeRefit_selection->Fill(ETA_simple.Theta()*TMath::RadToDeg());
                           // cout<<ETA_simple.Theta()*TMath::RadToDeg()<<endl;
                            //hthetahtoton2BeforeRefit
                           // cout<<ETA_simple.Phi()*TMath::RadToDeg()<<endl<<endl;;
                            hphiETABeforeRefit_selection->Fill(ETA_simple.Phi()*TMath::RadToDeg());
                           // hphiETABeforeRefit
                            
                            
                            
                            
                            
                            
                            
                            
                            
                          
                            
                            
                            
                            
                            
                            
                            hEETAAfterRefit->Fill(ETA_simple_new.E());
                            hPETAAfterRefit->Fill(std::sqrt(std::pow(ETA_simple_new.Px(),2)+std::pow(ETA_simple_new.Py(),2)+std::pow(ETA_simple_new.Pz(),2)));
                            double ETA_mass_refit=std::sqrt(std::pow(ETA_simple_new.E(),2)-std::pow(ETA_simple_new.Px(),2)-std::pow(ETA_simple_new.Py(),2)-std::pow(ETA_simple_new.Pz(),2));
                            hmassETAAfterRefit->Fill(ETA_mass_refit);
                            hthetaETAAfterRefit->Fill(ETA_simple_new.Theta()*180.0/M_PI);
                            hphiETAAfterRefit->Fill(ETA_simple_new.Phi()*180.0/M_PI);
                        }
                    }
                    candsFit.pop_back();
                    candsFit.pop_back();
                }
                    
            }
	    }
  
	    //**********************************************	    
    }//end event loop


    
	    
    output_file->cd();

    hbeta_mom->Write();
    hmass_mom->Write();
    hmass->Write();
    hbeta_mom_p->Write();
    hbeta_mom_pip->Write();
    hbeta_mom_pim->Write();

    hMM_pp->Write();
    hMM_pp_pmomCut1->Write();
    hMM_pp_pmomCut2->Write();
    hMM_pp_pmomCut3->Write();
    hMM_pp_pmomCut4->Write();
    
    hMM2_pp_noEl->Write();
    hMM2_pp_PT2->Write();
    hMM2_pp_PT3->Write();

    hMM2_pp->Write();
    hMM2_pp_pmomCut1->Write();
    hMM2_pp_pmomCut2->Write();
    hMM2_pp_pmomCut3->Write();
    hMM2_pp_pmomCut4->Write();
    
    
    hcl_size->Write();
    hnumOfNeutral->Write();  
    hnumOfNeutral_g100->Write();
    hbeta->Write();
    hbeta100->Write();
     hbeta150->Write();
    hbeta100_1->Write();
    htime->Write();
    htime_PT2->Write();
    htime_PT3->Write();
    htracklength->Write();




    htheta_mom->Write();
    htheta_mom_1->Write();
    htheta_mom_p_elast->Write();
    hdphi_pp->Write();
    htantan_dphi_pp->Write();
    
    hg_energy->Write();
    hg_energy_cl1->Write();
    hOA_gg->Write();
    

    for (int j=0;j<6;j++){
      hbeta_sec[j]->Write();
    }

    hmult_p->Write();
    henergyvscell->Write();
    htimevscell->Write();
    henergyvscell1->Write();
    hbetavscell1->Write();
    //*************************
    hinvM_gg->Write();
    hinvM_gg_ETA_PT3->Write();    
    hinvM_gg_omega_PT3->Write();    
    hinvM_gg_omega_ETA_PT3->Write();

    
    hMpETA_PT2->Write();
    hMpETA_PT3->Write();
    hMpETA_Mix->Write();

    hMMpETA->Write();
    hMM2pETA->Write();
    hMM2pETA_Mix->Write();

    
    hM3g->Write();
    hM2g_3g->Write();
    hM3g_ETA->Write();
    hMETAg_Mix->Write();
    hinvM_ggMixS->Write();
    hM3g_ETAMixS->Write();
    
    hMpippimETA_PT2->Write();
    hMpippimETA_PT3->Write();
    hMpippimETA_Mix->Write();
    hMggpippim_PT2->Write();
    hMggpippim_PT3->Write();
    hMggpippim_Mix->Write();
    hinvM_gg_pippim_PT3->Write();
    //**************
    
    hinvM_gg_PT2->Write();    
    hinvM_gg_PT3->Write();    
    hinvM_ggMix->Write();

//////////////////////////
hEPhtoton1EMC_gen->Write();
    hEPhtoton2EMC_gen->Write();
    //before

    hEPhtoton1beforeRefit->Write();
    hEPhtoton2beforeRefit->Write();



    //refit
    hPullEInvPhoton1Conv_mass_mix->Write();
    hPullThetaPhoton1Conv_mass_mix->Write();
    hPullPhiPhoton1Conv_mass_mix->Write();

    hEPhtotonAfterRefit_mix->Write();
    hPPhtotonAfterRefit_mix->Write();
    hmassPhtotonAfterRefit_mix->Write();
    hthetaPhtotonAfterRefit_mix->Write();
    hphiPhtotonAfterRefit_mix->Write();


    hPullEInvPhoton1Conv_mass->Write();
    hPullThetaPhoton1Conv_mass->Write();
    hPullPhiPhoton1Conv_mass->Write();
    hPullRPhoton1Conv_mass->Write();
    hPullZPhoton1Conv_mass->Write();
    hPullEInvPhoton2Conv_mass->Write();
    hPullThetaPhoton2Conv_mass->Write();
    hPullPhiPhoton2Conv_mass->Write();
    hPullRPhoton2Conv_mass->Write();
    hPullZPhoton2Conv_mass->Write();
    hEPhtoton1AfterRefit->Write();
    hPPhtoton1AfterRefit->Write();
    hmassPhtoton1AfterRefit->Write();
    hthetaPhtoton1AfterRefit->Write();
    hphiPhtoton1AfterRefit->Write();


    hEPhtoton2AfterRefit->Write();
    hPPhtoton2AfterRefit->Write();
    hmassPhtoton2AfterRefit->Write();
    hthetaPhtoton2AfterRefit->Write();
    hphiPhtoton2AfterRefit->Write();

    hEETAAfterRefit->Write();
    hPETAAfterRefit->Write();
    hmassETAAfterRefit->Write();
    hthetaETAAfterRefit->Write();
    hphiETAAfterRefit->Write();
    hIterations->Write();
    hChi2->Write();
    hProb->Write();

    //residuals
    hResidua_Refit_Reco_photon1_E->Write();
    hResidua_Refit_Reco_photon1_theta->Write();
    hResidua_Refit_Reco_photon1_phi->Write();
    hResidua_Refit_Reco_photon1_R->Write();
    hResidua_Refit_Reco_photon1_Z->Write();

    
    hEETABeforeRefit_selection->Write();                     
    hPETABeforeRefit_selection->Write();
    hmassETABeforeRefit_selection->Write();
    hthetahtoton2BeforeRefit_selection->Write();
    hphiETABeforeRefit_selection->Write();

    hphiETABeforeRefit->Write();



    hEETABeforeRefit->Write();

    hResidua_True_Reco_photon1_E->Write();


    hEPhtotonbeforeRefitcomb->Write();
    htheta_reco->Write();
    htheta_reco150->Write();
    hthetaETABeforeRefit_selection->Write();
    hthetaETABeforeRefit->Write();

    Eta_m->Write();
    Eta_theta->Write();
    Eta_phi->Write();
       
    EtaP_m->Write();
    EtaP_theta->Write();
    EtaP_phi->Write();
 
    Eta_m_selection->Write();
    Eta_theta_selection->Write();
    Eta_phi_selection->Write();


    Eta_m_before->Write();
    Eta_theta_before->Write();
    Eta_phi_before->Write();



///////////////////////////
    
    output_file->Close();
    cout << "writing root tree done" << endl;
    cout<<endl<<endl<<endl;
    cout<<"total: "<<total<<"  refited: "<<refited<<"  perc:"<<(double)refited/(double)total<<endl<<endl<<endl;
    //*********************************
    //*******************************

    timer.Stop();
    timer.Print();

    return 0;
}
