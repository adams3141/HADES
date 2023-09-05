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
//#include "heventmixer.h"
//#include "TMultiEventMixer.h" 


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

//////////GEANT
#include "fwdet_res.h"
#include "hgeantkine.h"
#include "hparticlecandsim.h"
#include "hparticletool.h"
#include "hparticlegeant.h"
#include "TIterator.h"
#include "heventheader.h"
#include "hmdcdef.h"
#include "hmdcsetup.h"
#include "hmdccal1.h"
#include "hmdccal2.h"
#include "tofdef.h"
#include "htofhit.h"
#include "rpcdef.h"
#include "hrpcraw.h"
#include "hrpchit.h"
#include "hrpccluster.h"
#include "hloop.h"
#include "hcategory.h"
#include <TCanvas.h>
#include <TStyle.h>
#include <sstream>
#include <TGraph.h>
#include <TGraphErrors.h>
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

    Int_t getMotherIndex(HGeantKine* particle)
    {
      Int_t trackID=particle->getTrack();
      HGeantKine* particleParent=particle->getParent(trackID);
      Int_t parentID=0;
      if(particleParent!=0)
	parentID=particleParent->getID();
      return parentID;
      
      //return particle->getGeneratorInfo1();
    }


    bool isLepton(HParticleCand* particle)
    {
      double delta=0.05;
      double mquality=particle->getRichMatchingQuality();

      double dphi=particle->getDeltaPhi();
      double dtheta=particle->getDeltaTheta();
      
          
      if(particle->isFlagBit(kIsUsed))
	{
	  if(mquality==-1 || mquality>5)
	    return false;
	  if(particle->getBeta()<(1-delta) || particle->getBeta()>(1+delta))
	    return false;
	  // if(dtheta<-0.4 || dtheta>0.4)
	  // return false;
	  // if(dphi*TMath::Sin(particle->getTheta())>0.4 || dphi*TMath::Sin(particle->getTheta())<-0.4)
	  // return false;
	}
      else
	return false;
     
      return true;
    }

/////////GEANT END
#define PI 3.14159265

#define PR(x)                                                                                      \
    //std::cout << "++DEBUG: " << #x << " = |" << x << "| (" << __FILE__ << ", " << __LINE__ << ")\n";

using namespace std;

const double D2R = 1.74532925199432955e-02;
const double R2D = 57.2957795130823229;


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
  
    //LOADING CATEGORIES
    // Check if loop was properly initialized
    if (!loop->setInput("+HEmcNeutralCand,+HEmcCluster,+HGeantKine"))
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
    
    HCategory *catGeant = loop->getCategory("HGeantKine");//ADDED geant
    HCategory *catEmcClus = loop->getCategory("HEmcCluster");// EMC gen
     HCategory *catGeaKine = loop->getCategory("HGeantKine");
  if (!catGeant)
  {
    std::cout << "No kineCat in input!" << std::endl;
    exit(1);
  }
  //GEANT
 HCategory * fCatGeantKine = nullptr;
    fCatGeantKine = HCategoryManager::getCategory(catGeantKine, kTRUE, "catGeantKine");
    if (!fCatGeantKine)
    {
        cout << "No catGeantKine!" << endl;
        exit(EXIT_FAILURE);  // do you want a brute force exit ?
    }


  //GEANT END

    HEnergyLossCorrPar dEdxCorr;
    dEdxCorr.setDefaultPar("feb22");    

    
    //----------------------------------------------------------------------------------------------
    // Setting parameters for loop over events
    //----------------------------------------------------------------------------------------------

    Int_t entries = loop->getEntries(); // Number of entries in loop
    int limit_sta = anapars.start;      // Limit START - Where to start the loop
    int limit_sto = 100000;                  // Limit STOP - Where to stop the loop

    if (anapars.events >= 0)
        limit_sto = limit_sta + anapars.events;
    else
        limit_sto = entries;

    if (limit_sto > entries) limit_sto = entries;

    //----------------------------------------------------------------------------------------------
    // Specifying output file
    //----------------------------------------------------------------------------------------------

    TFile* output_file = TFile::Open(anapars.outfile, "RECREATE");
    output_file->cd();
    cout << "NEW ROOT TREE " << endl;

    //----------------------------------------------------------------------------------------------
  
    TH2F *hbeta_mom= new TH2F("hbeta_mom","hbeta_mom",1000,-4000.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_p= new TH2F("hbeta_mom_p","hbeta_mom_p",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pip= new TH2F("hbeta_mom_pip","hbeta_mom_pip",1000,0.,4000.,700,0.,1.4);
    TH2F *hbeta_mom_pim= new TH2F("hbeta_mom_pim","hbeta_mom_pim",1000,-4000.,0.,700,0.,1.4);

    TH2F* hmass_mom = new TH2F("hmass_mom","hmass_mom; p*q [MeV/c]; mass",1000,-2000,4000,1000,0,2000);
    TH1F *hmass=new TH1F("hmass","hmass",1000,-2000,2000);
    
    TH1F *hbeta=new TH1F("hbeta","hbeta",1000,0,2);
    TH1F *hbeta150=new TH1F("hbeta150","hbeta150",1000,0,2);
    TH1F *hbeta500=new TH1F("hbeta500","hbeta500",1000,0,2);
    TH1F *hbeta100_1=new TH1F("hbeta100_1","hbeta100_1",1000,0,2);
    TH1F *hbeta500_1=new TH1F("hbeta500_1","hbeta500_1",1000,0,2);
    TH1F *htime=new TH1F("htime","htime",2200,-200,2000);
    
    TH1F *htracklength=new TH1F("htracklength","htracklength",2000,2000,4000);
    TH1F *hcl_size=new TH1F("hcl_size","hcl_size",20,0,20);

    TH1F *hnumOfNeutral=new TH1F("hnumOfNeutral","hnumOfNeutral",100,0,100);
    TH1F *hnumOfNeutral_g100=new TH1F("hnumOfNeutral_g100","hnumOfNeutral_g100",100,0,100);
    TH1F *hg_energy=new TH1F("hg_energy","hg_energy",2000,0,2000);
    TH1F *hg_energy_cl1=new TH1F("hg_energy_cl1","hg_energy_cl1",2000,0,2000);


    TH1F *hOA_gg=new TH1F("hOA_gg","hOA_gg",180,0,180);
    TH2F *henergyvscell = new TH2F("henergyvscell", "henergyvscell; sec*255+cell; energy [MeV]", 1530,0,1530,500,0,1500);
    TH2F *htimevscell = new TH2F("htimevscell", "htimevscell; sec*255+cell; time [ns]", 1530,0,1530,110,-10,100);
    TH2F *henergyvscell1 = new TH2F("henergyvscell1", "henergyvscell1; sec*255+cell; energy [MeV]", 1530,0,1530,500,0,1500);
    TH2F *hbetavscell1 = new TH2F("hbetavscell1", "hbetavscell1; sec*255+cell; beta", 1530,0,1530,200,0,2);



    TH1F *hinvM_gg=new TH1F("hinvM_gg","hinvM_gg",700,0,700);
    char name[200];

/////////////////////////////////////////////////////////////HISTROGRAMS
//EMC gen
TH1F * hEPhton1EMC_gen= new TH1F("hEPhton1EMC_gen","E Phtoton1 EMC gen",100,0,3000);
TH1F * hEPhton2EMC_gen= new TH1F("hEPhton2EMC_gen","E Phtoton2 EMC gen",100,0,3000);
TH1F * hEPhtonEMC_gen= new TH1F("hEPhtonEMC_gen","E Phtoton comb EMC gen",100,0,3000);
//PULLS
 TH1F * hPullEInvPhoton1Conv_mass_mix= new TH1F("Pull_E_photon_1_mix","Photon pull 1/E",5000,-10,10);
    TH1F * hPullThetaPhoton1Conv_mass_mix= new TH1F("Pull_theta_photon_1_mix","Photon pull #theta",5000,-10,10);
    TH1F * hPullPhiPhoton1Conv_mass_mix= new TH1F("Pull_phi_photon_1_mix","Photon pull #phi",5000,-10,10);
 


    TH1F * hPullEInvPhoton1Conv_mass= new TH1F("Pull_E_photon_1","Pull E photon 1",500,-10,10);
    TH1F * hPullThetaPhoton1Conv_mass= new TH1F("Pull_theta_photon_1","Pull theta photon 1",500,-10,10);
    TH1F * hPullPhiPhoton1Conv_mass= new TH1F("Pull_phi_photon_1","Pull phiphoton 1",500,-10,10);
    TH1F * hPullRPhoton1Conv_mass= new TH1F("Pull_R_photon_1","Pull R photon 1",10000,-1,1);
    TH1F * hPullZPhoton1Conv_mass= new TH1F("Pull_Z_photon_1","Pull Z photon 1",10000,-1,1);

    TH1F * hPullEInvPhoton2Conv_mass= new TH1F("Pull_E_photon_2","Pull E photon 2",500,-10,10);
    TH1F * hPullThetaPhoton2Conv_mass= new TH1F("Pull_theta_photon_2","Pull theta photon 2",500,-10,10);
    TH1F * hPullPhiPhoton2Conv_mass= new TH1F("Pull_phi_photon_2","Pull phi photon 2",500,-10,10);
    TH1F * hPullRPhoton2Conv_mass= new TH1F("Pull_R_photon_2","Pull R photon 2",10000,-1,1);
    TH1F * hPullZPhoton2Conv_mass= new TH1F("Pull_Z_photon_2","Pull Z photon 2",10000,-1,1);

    //photons before refit
    TH1F * hEPhtoton1beforeRefit= new TH1F("Before Energy Phtoton1 ","E Phtoton1 Before Refit",100,0,3000);
    TH1F * hEPhtoton2beforeRefit= new TH1F("Before Energy Photon2","E Phtoton2 Before Refit",100,0,3000);
    TH1F * hEPhtotonbeforeRefitcomb= new TH1F("Before Energy Photon comb","E Phtoton comb Before Refit",100,0,3000);

    //photons after refit
    TH1F * hEPhtotonAfterRefit_mix= new TH1F("hEPhtotonAfterRefit_mix","Refit Photon E",100,0,3000);
    TH1F * hPPhtotonAfterRefit_mix= new TH1F("hPPhtotonAfterRefit_mix","Refit Photon P",100,0,3000);
    TH1F * hmassPhtotonAfterRefit_mix= new TH1F("hmassPhtotonAfterRefit_mix","mass Phtoton Refit_mix",100,-100,100);
    TH1F * hthetaPhtotonAfterRefit_mix= new TH1F("hthetaphtotonAfterRefit_mix","Refit Photon #theta",1000,-180,180);
    TH1F * hphiPhtotonAfterRefit_mix= new TH1F("hphiPhtotonAfterRefit_mix","Refit Photon #phi",300,-180,180);

     TH1F * hphiPhtotonBeforeRefit_mix= new TH1F("hphiPhtotonBeforeRefit_mix","Refit Photon #phi",300,-180,180);


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


    TH1F * hEPi0AfterRefit= new TH1F("hEPi0AfterRefit","hEPi0AfterRefit",1000,0,4000);
    TH1F * hPPi0AfterRefit= new TH1F("hPPi0AfterRefit","hPPi0AfterRefit",1000,0,4000);
    TH1F * hmassPi0AfterRefit= new TH1F("hmassPi0AfterRefit","hmassPi0AfterRefit",100000000,100,200);
    TH1F * hthetaPi0AfterRefit= new TH1F("hthetahtoton2AfterRefit","hthetaPi0AfterRefit",1000,-180,180);
    TH1F * hphiPi0AfterRefit= new TH1F("hphiPi0AfterRefit","hEphiPi0AfterRefit",1000,-180,180);


    //pi0 simple- from g+g without refit
     
    TH1F * hEPi0BeforeRefit_selection= new TH1F("hEPi0ABeforeRefit","hEPi0BeforeRefit Selection",1000,0,4000);
    TH1F * hEPi0BeforeRefit= new TH1F("hEPi0BeforeRefit","hEPi0BeforeRefit",1000,0,4000);

    TH1F * hPPi0BeforeRefit_selection= new TH1F("hPPi0BeforeRefit_selection","hPPi0BeforeRefit_selection",100,0,4000);
    TH1F * hPPi0BeforeRefit= new TH1F("hPPi0BeforeRefit","hPPi0BeforeRefit",100,0,4000);
    TH1F * hmassPi0BeforeRefit_selection= new TH1F("hmassPi0BeforeRefit_selection","hmassPi0BeforeRefit_selection",10000,100,200);
    TH1F * hmassPi0BeforeRefit= new TH1F("hmassPi0BeforeRefit","hmassPi0BeforeRefit",10000,100,200);
    TH1F * hthetahtoton2BeforeRefit_selection= new TH1F("hthetahtoton2BeforeRefit_selection","hthetaPi0BeforeRefit_selection",1000,-180,180);
    TH1F * hthetahtoton2BeforeRefit= new TH1F("hthetahtoton2BeforeRefit","hthetaPi0BeforeRefit",300,-180,180);

    TH1F * hphiPi0BeforeRefit_selection= new TH1F("hphiPi0BeforeRefit_selection","hEphiPi0BeforeRefit_selection",1000,-180,180);
    TH1F * hphiPi0BeforeRefit= new TH1F("hphiPi0BeforeReft","hEphiPi0BeforeRefit",1000,-180,180);

    TH1F * hthetaPi0BeforeRefit_selection= new TH1F("hthetaPi0BeforeRefit_selection","hEphiPi0BeforeRefit_selection",1000,-180,180);
    TH1F * hthetaPi0BeforeRefit= new TH1F("hthetaPi0BeforeReft","hEphiPi0BeforeRefit",1000,-180,180);
    //Residuals
    TH1F * hResidua_Refit_Reco_photon1_E= new TH1F("hResidua_Refit_Reco_photon1_E","Residua Refit-Reco E 1/E",10000,-100,100);
    TH1F * hResidua_Refit_Reco_photon1_theta= new TH1F("hResidua_Refit_Reco_photon1_theta","Residua Refit-Reco theta photon1",10000,-3*180,3*180);
    TH1F * hResidua_Refit_Reco_photon1_phi= new TH1F("hResidua_Refit_Reco_photon1_phi","Residua Refit-Reco phi photon1",10000,-3*180,3*180);
    TH1F * hResidua_Refit_Reco_photon1_R= new TH1F("Residua Refit-Reco R photon1","Residua Refit-Reco R photon1",10000,-100,100);
    TH1F * hResidua_Refit_Reco_photon1_Z= new TH1F("Residua Refit-Reco Z photon1","Residua Refit-Reco Z photon1",10000,-100,100);

     TH1F * hResidua_True_Reco_photon1_E= new TH1F("hResidua_True_Reco_photon1_E","Residua True-Reco E 1/E",10000,-100,100);
    TH1F * hResidua_True_Reco_photon1_theta= new TH1F("hResidua_True_Reco_photon1_theta","Residua True-Reco theta photon1",10000,-3*180,3*180);
    TH1F * hResidua_True_Reco_photon1_phi= new TH1F("hResidua_True_Reco_photon1_phi","Residua True-Reco phi photon1",10000,-3*180,3*180);
    TH1F * hResidua_True_Reco_photon1_R= new TH1F("Residua True-Reco R photon1","Residua True-Reco R photon1",10000,-100,100);
    TH1F * hResidua_True_Reco_photon1_Z= new TH1F("Residua True-Reco Z photon1","Residua True-Reco Z photon1",10000,-100,100);


    TH1F * hResidua_True_Refit_photon1_E= new TH1F("hResidua_True_Refit_photon1_E","Residua True-Refit E 1/E",10000,-100,100);
    TH1F * hResidua_True_Refit_photon1_theta= new TH1F("hResidua_True_Refit_photon1_theta","Residua True-Refit theta photon1",10000,-3*180,3*180);
    TH1F * hResidua_True_Refit_photon1_phi= new TH1F("hResidua_True_Refit_photon1_phi","Residua True-Refit phi photon1",10000,-3*180,3*180);
    TH1F * hResidua_True_Refit_photon1_R= new TH1F("Residua True-Refit R photon1","Residua True-Refit R photon1",10000,-100,100);
    TH1F * hResidua_True_Refit_photon1_Z= new TH1F("Residua True-Refit Z photon1","Residua True-Refit Z photon1",10000,-100,100);

  

    //Refit GEANT- performing refit on geant tracks
    //PULLS
    TH1F * hPullEInvPhoton1Conv_mass_geant= new TH1F("Pull_E_photon_1_GEANT","Pull E photon 1_GEANT",500,-10,10);
    TH1F * hPullThetaPhoton1Conv_mass_geant= new TH1F("Pull_theta_photon_1_GEANT","Pull theta photon 1_GEANT",500,-10,10);
    TH1F * hPullPhiPhoton1Conv_mass_geant= new TH1F("Pull_phi_photon_1_GEANT","Pull phiphoton 1_GEANT",500,-10,10);
    TH1F * hPullRPhoton1Conv_mass_geant= new TH1F("Pull_R_photon_1_GEANT","Pull R photon 1_GEANT",10000,-1,1);
    TH1F * hPullZPhoton1Conv_mass_geant= new TH1F("Pull_Z_photon_1_GEANT","Pull Z photon 1_GEANT",10000,-1,1);

    TH1F * hPullEInvPhoton2Conv_mass_geant= new TH1F("Pull_E_photon_2_GEANT","Pull E photon 2_GEANT",500,-10,10);
    TH1F * hPullThetaPhoton2Conv_mass_geant= new TH1F("Pull_theta_photon_2_GEANT","Pull theta photon 2_GEANT",500,-10,10);
    TH1F * hPullPhiPhoton2Conv_mass_geant= new TH1F("Pull_phi_photon_2_GEANT","Pull phi photon 2_GEANT",500,-10,10);
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

    TH1I * hIterations_GEANT = new TH1I("hIterations","hIterations",10,0,10);
    TH1F * hChi2_GEANT = new TH1F("hChi2_GEANT","hChi2_GEANT",100,0,1);
    TH1F * hProb_GEANT = new TH1F("hProb_GEANT","hProb_GEANT",100,0,1);


    TH1F * hEPi0AfterRefit_GEANT= new TH1F("hEPi0AfterRefit_GEANT","hEPi0AfterRefit_GEANT",100,0,4000);
    TH1F * hPPi0AfterRefit_GEANT= new TH1F("hPPi0AfterRefit_GEANT","hPPi0AfterRefit_GEANT",100,0,4000);
    TH1F * hmassPi0AfterRefit_GEANT= new TH1F("hmassPi0AfterRefit_GEANT","hmassPi0AfterRefit_GEANT",10000,100,200);
    TH1F * hthetaPi0AfterRefit_GEANT= new TH1F("hthetahtoton2AfterRefit_GEANT","hthetaPi0AfterRefit_GEANT",1000,-180,180);
    TH1F * hphiPi0AfterRefit_GEANT= new TH1F("hphiPi0AfterRefit_GEANT","hEphiPi0AfterRefit_GEANT",300,-180,180);


    TH1F * htheta_GEANT= new TH1F("htheta_GEANT","htheta_GEANT",1000,-180,180);
    TH1F * htheta_GEANT150= new TH1F("htheta_GEANT150","htheta_GEANT150",1000,-180,180);
    TH1F * htheta_reco= new TH1F("htheta_reco","htheta_reco",1000,-180,180);
    TH1F * htheta_reco150= new TH1F("htheta_reco150","htheta_reco150",1000,-180,180);






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
    vector<TLorentzVector> lv_gMix, lv_pi0;

    //*************************************    

    
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
    for (Int_t ev = limit_sta; ev < limit_sto; ev++) // event loop
    {
          if (ev % 10000 == 0)
        {
            printf("Event nr.: %d, progress: %.2f%%\n", ev, (double)(ev - limit_sta) / (limit_sto - limit_sta) * 100.);
        }

        /*Int_t nbytes =*/loop->nextEvent(ev); // get next event. categories will be cleared before
        /////////////////////////////////////////////
        //EMC GEN
        //hEPhton1EMC_gen
        ////
        //catGeaKine
 //HGeantKine* pKine2  = (HGeantKine*)catGeaKine->getObject(track2-1);
   //     id2   = pKine2->getID();
     //   mom2  = pKine2->getTotalMomentum();
        //  Int_t Ntracks = emc->getNTracks();
     

        ///
    vector<HGeantKine*> emc_gen_kine;
         Int_t nEmcCluster = catEmcClus->getEntries();
    
    for(Int_t e=0; e<nEmcCluster; e++) {
      HEmcClusterSim* emc = (HEmcClusterSim*)catEmcClus->getObject(e);
       Int_t track = emc->getTrack();
      HGeantKine* pKine  = (HGeantKine*)catGeaKine->getObject(track-1);
      emc_gen_kine.push_back(pKine);
      hEPhtonEMC_gen->Fill(pKine->getTotalMomentum());
      //cout<<pKine->getTotalMomentum()<<endl;;
      //hEPhtonEMC_gen->Fill(E2->getTotalMomentum());
     double mom2  = pKine->getTotalMomentum();
      //cout<<mom2<<endl;
/*
          pid.sector=emc->getSector();
	  pid.theta=emc->getTheta();
          pid.phi=emc->getPhi();
          pid.size=emc->getNCells();
          pid.cell=systemMap[emc->getCell()];
          pid.matchedCells = emc->getNMatchedCells(); // 0 : not matched (neutral particle)
	  //pid.matchedPart=emc_ModMult[emc->getSector()][systemMap[emc->getCell()]]++ ;
	  pid.matchedPart=-1 ;
          pid.flag=emc->isUsedInParticleCand();
          pid.energy=emc->getEnergy();
          pid.max_energy=emc->getMaxEnergy();
          pid.time=emc->getTime();
	  
          // Calculate Beta:
          Double_t trackLength = TMath::Sqrt(emc->getXLab()*emc->getXLab()+emc->getYLab()*emc->getYLab()+(emc->getZLab()+30)*(emc->getZLab()+30));  
          Double_t Beta        = (trackLength/1000.)/(pid.time*1.e-9) / TMath::C();
	  pid.distance=trackLength;
	  pid.beta=Beta;
	  
          //RPC variables
	  
	  pid.qualityDThDPh=emc->getQualDThDPh();
	  pid.qualityDTime=emc->getQualDTime();
          pid.tof_rpc=-1000;
          pid.theta_rpc=-1000;
          pid.phi_rpc=-1000;
	  /*
	  if(emc->getRpcIndex()>-1&&emc->getNMatchedCells()>0) {
            Int_t RPCInd=emc->getRpcIndex();
	    HRpcCluster* rpcClst   = HCategoryManager::getObject(rpcClst  ,fCatRpcCluster,RPCInd);
	    pid.tof_rpc=rpcClst->getTof();
	    pid.theta_rpc=rpcClst->getTheta();
	    pid.phi_rpc=rpcClst->getPhi();
	  }
	  */
      /*
      
          mult_emc++;
          pid.downflag=0;
          if( mult_emc>1) pid.downflag = -1; //only 1st tree entry in event
                                              //had downflag=0
                                              //downflag is not used anyway
	        
      Int_t Ntracks = emc->getNTracks();
      Int_t track = emc->getTrack();
      HGeantKine* pKine  = (HGeantKine*)catGeaKine->getObject(track-1);
      Float_t     mom  = pKine->getTotalMomentum();
      Int_t       id      = pKine->getID();
      HGeantKine* prnt    = pKine->getParent(track);
      Int_t       idp     = prnt!=NULL ? prnt->getID()    : 0; 
      Int_t       trkPi0I = prnt!=NULL ? prnt->getTrack() : 0;

      Int_t track2  = -1;
      Float_t mom2  = -1;
      Int_t   id2    = -1;
      if(Ntracks>1) {
	track2 = emc->getTrack(1);
        HGeantKine* pKine2  = (HGeantKine*)catGeaKine->getObject(track2-1);
        id2   = pKine2->getID();
        mom2  = pKine2->getTotalMomentum();
      }
      
      pid.ntracks= Ntracks;
      pid.id=id;
      pid.mom_geant   = mom;
      pid.parentid=idp;
      pid.id2=id2;
      pid.mom_geant2   = mom2;
	
      if(emc->getNMatchedCells()==0&&Beta>0.95) nGammas++; 

      //if(emc->getNMatchedCells()!=0) continue; 
      nt->Fill();
*/
    }  //loop over cluster



        ////////////////////////////////////////
        //geant
        HParticleCandSim* partH=nullptr;
	    //HFwDetCandSim*  partFT=nullptr;
         std::vector<HGeantKine *> geantkine_vec;

	    HGeantKine* kine=nullptr;

	
	    Int_t isDilepton=0;
	    HParticleTool tool;
	
	    HGeomVector base_Tg, dir_Tg;
			  
	    base_Tg.setX(0);
	    base_Tg.setY(0);
	    base_Tg.setZ(-25.);
	    dir_Tg.setX(0);
	    dir_Tg.setY(0);
	    dir_Tg.setZ(1);

	
	    TLorentzVector *lv_pp=new TLorentzVector;
	    TLorentzVector *lv_pppippim=new TLorentzVector;

        int knum=fCatGeantKine->getEntries();
        for(int p=0;p<knum;p++) 
        {
	        kine=HCategoryManager::getObject(kine, fCatGeantKine,p);
	        int kineID=kine->getID();
	        int mech=kine->getMechanism();
	        int kineparentID=getMotherIndex(kine);

            float weight=kine->getGeneratorWeight();
	        //if(kine->isInAcceptance()) cout<<"::::::::"<<endl;
	  
	        //hPRmotherindex->Fill(kineparentID);

	        int flep=0;
	        int flem=0;
            if(kineID==1&&kineparentID==7)//photon from pi0
            {
              //  cout<<"tak"<<endl;
            }
             if(kineID==1&&kineparentID==0)//photon from pi0
            {
              //  cout<<"tak"<<endl;
            }
             if(kineID==1)//photon from pi0
            {
                geantkine_vec.push_back(kine);
                htheta_GEANT->Fill(kine->getThetaDeg());
                //cout<<kine->getTotalMomentum()<<endl;
                if(kine->getTotalMomentum()>0.150)
                {
                  htheta_GEANT150->Fill(kine->getThetaDeg());
                }
                
              //  cout<<"tak"<<endl;
            }
        }
        if(geantkine_vec.size()==2)
        {   
            HGeantKine  phot_geant_cand_1=*geantkine_vec[0];
            HGeantKine  phot_geant_cand_2=*geantkine_vec[1];
             //      TH1F * htheta_GEANT= new TH1F("htheta_GEANT","htheta_GEANT",1000,-180,180);
            //cout<<phot_geant_cand_1.E()<<endl;
            //hEPhtoton1beforeRefit_GEANT->Fill()
           // hEPhtoton2beforeRefit_GEANT
        }

        //GEANT END
        //simulation reconstruction

      

	    lv_neutr.clear();
	    lv_neutr1.clear();

	

	
	
        HEventHeader* event_header = NULL;
        if (!(event_header = gHades->getCurrentEvent()->getHeader())) continue;

        Int_t TBit = (Int_t)event_header->getTBit();
        Double_t VertexX = event_header->getVertexReco().getX();
        Double_t VertexY = event_header->getVertexReco().getY();
        Double_t VertexZ = event_header->getVertexReco().getZ();

	
	
	    //****************************************************
       
	    Int_t nNeutral_ev = fEmcNeutralCand->getEntries();

	    hnumOfNeutral->Fill(nNeutral_ev);
	
        std::vector<KFitParticle *> kFit_gamma;
        // std::vector<HGeantKine *> geantkine_vec;
	    for (int j = 0; j < nNeutral_ev; ++j)
	    {
        //LOADING DATA FROM EVENT
            HEmcNeutralCand* neutr_cand = HCategoryManager::getObject(neutr_cand, fEmcNeutralCand, j);
            neutr_cand->calc4vectorProperties(0.0);//zmiana
            Float_t dist  = neutr_cand->getDistanceToEmc();
            Int_t ind=neutr_cand->getEmcClusterIndex();


            //geant
            HGeantKine *geantkine;
            geantkine = HCategoryManager::getObject(geantkine, catGeant, j);
            geantkine_vec.push_back(geantkine);
            //cout<<geantkine->getE()<<endl;


            //

            HEmcCluster *cl=nullptr;
            cl=HCategoryManager::getObject(cl, fEmcCluster, ind);
            Int_t cl_size = cl->getNCells();

            Int_t sec = cl->getSector();
            Int_t cel = cl->getCell();		
            if(cel<33) continue;
            

            Double_t energy  = cl->getEnergy();
            //HGeomVector trackVec(cl->getXLab()- VertexX, cl->getYLab()- VertexY, cl->getZLab() - VertexZ);
            //Double_t trackLength = trackVec.length();
            //trackVec  *= (energy/trackLength);

            Double_t tof  =cl->getTime();
            Float_t theta = cl->getTheta();
            Float_t phi = cl->getPhi();
            Double_t beta = neutr_cand->getBeta();
            int pid=neutr_cand->getPID();
            
            hcl_size->Fill(cl_size);
            hbeta->Fill(beta);
            htime->Fill(tof);
            //htracklength->Fill(trackLength);
            htracklength->Fill(dist);

            
            henergyvscell->Fill(sec*200+cel,energy);
            htimevscell->Fill(sec*200+cel,tof);
            
            if(energy>150 && cl_size==1)hbeta150->Fill(beta);
            if(energy>500 && cl_size==1)hbeta500->Fill(beta);
            
            if(energy>100 && cl_size==1 && tof>0 && tof<50){

            hbeta100_1->Fill(beta);
            henergyvscell1->Fill(sec*200+cel,energy);
            hbetavscell1->Fill(sec*200+cel,beta);
            }
            if(energy>500 && cl_size==1 && tof>0 && tof<50){
            hbeta500_1->Fill(beta);
            }
            hg_energy->Fill(energy);
            if(cl_size==1)hg_energy_cl1->Fill(energy);


            TLorentzVector lvg=* neutr_cand;
            //////////////////////////////////////////////////////
            HVertex primVertex = event_header->getVertexReco();



            /////////////////////////////////////////////////////////
            KFitParticle* candidate = new KFitParticle(&lvg,neutr_cand->getR(),neutr_cand->getZ());  //puste? 
             htheta_reco->Fill(candidate->getTheta()*TMath::RadToDeg());
             
            if ( pid==1&&energy>150){
              //cout<<candidate->getTheta()<<endl;
                lv_neutr1.push_back(lvg);	 //vector of TLorentz
                
                double errors[] = {0.055/(lvg.E()*std::sqrt(lvg.E())), 2.5*deg2rad,
                                    2.5*deg2rad,0.00001,0.00001};// need to be adjusted for gammmas /std::pow(lvg.E(),2)
                    FillData(neutr_cand, candidate, errors,0);//filling atributes of kfitparticles
                    
                kFit_gamma.push_back(candidate);
                gamma_vector.push_back(neutr_cand);
                    hEPhtotonbeforeRefitcomb->Fill(candidate->Energy());
              
  
            htheta_reco150->Fill(candidate->getTheta()*TMath::RadToDeg());
    //cout<<candidate->getTheta()<<endl;
               // cout<<lvg.E()<<" to  "<<candidate->Energy()<<endl;
            }



	    }
        //*********************************************************
        //**************   gg  *****************
    

        hnumOfNeutral_g100->Fill(lv_neutr1.size());
        int gamma_num=kFit_gamma.size();
        if(gamma_num>=2)//at least 2 photons  were selected
        {
            std::vector<KFitParticle> candsFit; //vector for photon candidates to refit
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
                    hEPi0BeforeRefit->Fill(lvg2.E());
                    if (oAngle1<oAngleCut) continue;
                    
                    //************** pi0 ************************
                    hinvM_gg->Fill(mass_gg);
                   // KFitParticle cand2Dec = *kFit_gamma[ii];
                    
                    KFitParticle *cand1DecPtr = kFit_gamma[ii];
                    KFitParticle cand1Dec = *cand1DecPtr;
                    candsFit.clear();
                    candsFit.push_back(cand1Dec);
                    KFitParticle *cand11DecPtr = kFit_gamma[jj];
                    KFitParticle  cand11Dec = *cand11DecPtr;
                    hEPhtoton1beforeRefit->Fill(cand1Dec.Energy());//energy before refit//gamma_vector[ii]->E());//
                    double EE1=cand1Dec.Energy(),EE2=cand11Dec.Energy();

                    hEPhtoton2beforeRefit->Fill(cand11Dec.Energy());//gamma_vector[jj]->E());//
                   // hEPhtotonbeforeRefitcomb->Fill(cand1Dec.Energy());
                   // hEPhtotonbeforeRefitcomb->Fill(cand11Dec.Energy());

                    //// emc gen
                    double E1=emc_gen_kine[ii]->getTotalMomentum();
                    double E2=emc_gen_kine[jj]->getTotalMomentum();
                    hEPhton1EMC_gen->Fill(E1);
                    hEPhton2EMC_gen->Fill(E2);
                  

                    ///
                    candsFit.push_back(cand11Dec);//filling vector with second photon cadidate
                    KinFitter FitterMass(candsFit);// fitter declaration
                    FitterMass.setVerbosity(0);//fitter does not print any information
                    FitterMass.addMassConstraint(134.9768);//mass constraint
                    FitterMass.setNumberOfIterations(10);
                    FitterMass.setConvergenceCriteria(1,9999.9,9999.9);//1 is for chi2 further ones are constraints for parameters better not to touch
                    //FitterMass.SetVerbosity(2);//fitter prints a lot of infomation
                    FitterMass.fit();
                    total++;

                  TLorentzVector pi0_simple=*kFit_gamma[ii]+*kFit_gamma[jj];
                            hEPi0BeforeRefit->Fill(pi0_simple.E());
                            
                            hPPi0BeforeRefit->Fill(pi0_simple.P());
                            //hPPi0BeforeRefit
                            hmassPi0BeforeRefit->Fill(pi0_simple.Mag());
                            //hmassPi0BeforeRefit
                            hthetahtoton2BeforeRefit->Fill(pi0_simple.Theta()*TMath::RadToDeg());
                           // cout<<pi0_simple.Theta()*TMath::RadToDeg()<<endl;
                            //hthetahtoton2BeforeRefit
                           // cout<<pi0_simple.Phi()*TMath::RadToDeg()<<endl<<endl;;
                            hphiPi0BeforeRefit->Fill(pi0_simple.Phi()*TMath::RadToDeg());
                            hthetaPi0BeforeRefit->Fill(pi0_simple.Theta()*TMath::RadToDeg());

                    if(FitterMass.isConverged())
                    {
                        hIterations->Fill(FitterMass.getIteration());
                        hChi2->Fill(FitterMass.getChi2());
                        hProb->Fill(FitterMass.getProb());
                            
                        // Get the pull distributions
                        if (FitterMass.getProb() > 0.01)
                        {
                          refited++;
                        }
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
                        
                       

                                

                            KFitParticle cand1mass = FitterMass.getDaughter(0); 
                            KFitParticle cand2mass = FitterMass.getDaughter(1);



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
                            hphiPhtotonBeforeRefit_mix->Fill(cand1mass.Phi()*180.0/M_PI);
                            hphiPhtotonBeforeRefit_mix->Fill(cand2mass.Phi()*180.0/M_PI);


                            //photon residuals
                            hResidua_Refit_Reco_photon1_E->Fill(cand1mass.Energy()-cand1Dec.Energy());
                            //cout<<cand1mass.Energy()-cand1Dec.Energy()<<"   "<<EE1<<"  "<<cand1mass.Energy()<<"   "<<cand11Dec.Energy()<<endl;
                           // hResidua_Refit_Reco_photon1_E->Fill(cand1mass.Energy()-cand1Dec.Energy());
                            hResidua_Refit_Reco_photon1_theta->Fill(cand1mass.Theta()*R2D-cand1Dec.Theta()*R2D);
                            hResidua_Refit_Reco_photon1_phi->Fill(cand1mass.Phi()*R2D-cand1Dec.Phi()*R2D);
                           // cout<<cand1mass.Phi()*R2D<<"  "<<cand1Dec.Phi()*R2D<<"  "<<cand1mass.getPhi()-cand1Dec.Phi()*R2D<<endl;
                            hResidua_Refit_Reco_photon1_R->Fill(cand1mass.getR()-cand1Dec.getR());
                            hResidua_Refit_Reco_photon1_Z->Fill(cand1mass.getZ()-cand1Dec.getZ());

                            hResidua_True_Refit_photon1_E->Fill(emc_gen_kine[ii]->getTotalMomentum()-cand1mass.Energy());
                            hResidua_True_Refit_photon1_theta->Fill(emc_gen_kine[ii]->getThetaDeg()-cand1mass.Theta()*R2D);
                            hResidua_True_Refit_photon1_phi->Fill(emc_gen_kine[ii]->getPhiDeg()-cand1mass.Phi()*R2D);
                          //  cout<<emc_gen_kine[ii]->getPhiDeg()<<" "<<cand1mass.Phi()*R2D<<endl;



                            hResidua_True_Reco_photon1_E->Fill(emc_gen_kine[ii]->getTotalMomentum()-cand1Dec.Energy());
                            hResidua_True_Reco_photon1_theta->Fill(emc_gen_kine[ii]->getThetaDeg()-cand1Dec.Theta()*R2D);
                            hResidua_True_Reco_photon1_phi->Fill(emc_gen_kine[ii]->getPhiDeg()-cand1Dec.Phi()*R2D);
                            //hResidua_True_Refit_photon1_R->Fill(cand1mass.getR()-cand1Dec.getR());
                            //hResidua_True_Refit_photon1_Z->Fill(cand1mass.getZ()-cand1Dec.getZ());
                            //cout<<emc_gen_kine[ii]->getTotalMomentum()<<"    "<<cand1mass.Energy()<<endl;
                            //hResidua_True_Refit_photon1_theta->Fill(geantkine_vec[ii]->getTheta()-cand1mass.Theta());
                            //hResidua_True_Refit_photon1_phi
                            //hResidua_True_Refit_photon1_R
                            //hResidua_True_Refit_photon1_Z



                            KFitParticle pi0_refit = FitterMass.getMother(); //obtaining pi0 mother partivle from refit
                            
                            
                            
                            hEPi0BeforeRefit_selection->Fill(pi0_simple.E());
                            
                            hPPi0BeforeRefit_selection->Fill(pi0_simple.P());
                            //hPPi0BeforeRefit
                            hmassPi0BeforeRefit_selection->Fill(pi0_simple.Mag());
                            //hmassPi0BeforeRefit
                            hthetahtoton2BeforeRefit_selection->Fill(pi0_simple.Theta()*TMath::RadToDeg());
                           // cout<<pi0_simple.Theta()*TMath::RadToDeg()<<endl;
                            //hthetahtoton2BeforeRefit
                           // cout<<pi0_simple.Phi()*TMath::RadToDeg()<<endl<<endl;;
                            hphiPi0BeforeRefit_selection->Fill(pi0_simple.Phi()*TMath::RadToDeg());
                            hthetaPi0BeforeRefit_selection->Fill(pi0_simple.Theta()*TMath::RadToDeg());
                           // hphiPi0BeforeRefit
                            
                            
                            
                            
                            
                            
                            
                            
                            
                            TLorentzVector pi0_simple_new=cand1mass+cand2mass;
                            
                            
                            
                            
                            
                            
                            hEPi0AfterRefit->Fill(pi0_simple_new.E());
                            hPPi0AfterRefit->Fill(std::sqrt(std::pow(pi0_simple_new.Px(),2)+std::pow(pi0_simple_new.Py(),2)+std::pow(pi0_simple_new.Pz(),2)));
                            double pi0_mass_refit=std::sqrt(std::pow(pi0_simple_new.E(),2)-std::pow(pi0_simple_new.Px(),2)-std::pow(pi0_simple_new.Py(),2)-std::pow(pi0_simple_new.Pz(),2));
                            hmassPi0AfterRefit->Fill(pi0_mass_refit);
                            hthetaPi0AfterRefit->Fill(pi0_simple_new.Theta()*180.0/M_PI);
                            hphiPi0AfterRefit->Fill(pi0_simple_new.Phi()*180.0/M_PI);
                       // }//prob
                    }
                    candsFit.pop_back();
                    candsFit.pop_back();
                }
                    
            }
	    }

	    //**********************************************************

    }//end event loop


    
	        //SAVING RESULTS
    output_file->cd();

    //hbeta_mom->Write();
    //hmass_mom->Write();
    //hmass->Write();
    //hbeta_mom_p->Write();
    //hbeta_mom_pip->Write();
    //hbeta_mom_pim->Write();

    /*
    hcl_size->Write();
    hnumOfNeutral->Write();  
    hnumOfNeutral_g100->Write();
    
    hbeta100->Write();
    hbeta100_1->Write();
    htime->Write();
    htracklength->Write();
    hg_energy->Write();
    hg_energy_cl1->Write();
    hOA_gg->Write();
    

    henergyvscell->Write();
    htimevscell->Write();
    henergyvscell1->Write();
    hbetavscell1->Write();
    */
    //*************************
    hbeta->Write();
    hbeta150->Write();
    //EMC gen
    
    hEPhton1EMC_gen->Write();
    hEPhton2EMC_gen->Write();
    hEPhtonEMC_gen->Write();
    //before

    hEPhtoton1beforeRefit->Write();
    hEPhtoton2beforeRefit->Write();
    hEPhtotonbeforeRefitcomb->Write();



    //refit
    hPullEInvPhoton1Conv_mass_mix->Write();
    hPullThetaPhoton1Conv_mass_mix->Write();
    hPullPhiPhoton1Conv_mass_mix->Write();

    hEPhtotonAfterRefit_mix->Write();
    hPPhtotonAfterRefit_mix->Write();
    hmassPhtotonAfterRefit_mix->Write();
    hthetaPhtotonAfterRefit_mix->Write();
    hphiPhtotonAfterRefit_mix->Write();

    hphiPhtotonBeforeRefit_mix->Write();


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

    hEPi0AfterRefit->Write();
    hPPi0AfterRefit->Write();
    hmassPi0AfterRefit->Write();
    hthetaPi0AfterRefit->Write();
    hphiPi0AfterRefit->Write();
    hIterations->Write();
    hChi2->Write();
    hProb->Write();

    //residuals
    hResidua_Refit_Reco_photon1_E->Write();
    hResidua_Refit_Reco_photon1_theta->Write();
    hResidua_Refit_Reco_photon1_phi->Write();
    hResidua_Refit_Reco_photon1_R->Write();
    hResidua_Refit_Reco_photon1_Z->Write();

    
    hEPi0BeforeRefit_selection->Write();                     
    hPPi0BeforeRefit_selection->Write();
    hmassPi0BeforeRefit_selection->Write();
    hthetahtoton2BeforeRefit_selection->Write();
    hphiPi0BeforeRefit_selection->Write();

 
                            
    hPPi0BeforeRefit->Write();

    hmassPi0BeforeRefit->Write();
    hthetahtoton2BeforeRefit->Write();
    hphiPi0BeforeRefit->Write();
    hthetaPi0BeforeRefit->Write();



    hEPi0BeforeRefit->Write();

    hResidua_True_Reco_photon1_E->Write();


    hResidua_True_Refit_photon1_E->Write();
    hResidua_True_Refit_photon1_theta->Write();
    hResidua_True_Refit_photon1_phi->Write();


    hResidua_True_Reco_photon1_E->Write();
    hResidua_True_Reco_photon1_theta->Write();
    hResidua_True_Reco_photon1_phi->Write();


    htheta_GEANT->Write();
    htheta_GEANT150->Write();
    htheta_reco->Write();
    htheta_reco150->Write();

    hinvM_gg->Write();

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
