#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>

using std::cout;
using std::endl;

void fit_hist(TH1* h, TF1 * f, float xmin, float xmax, TLatex * tex, float label_x, float mu_labeyl_y, float sigma_label_y)
{
	h->Fit(f, "", "S", xmin, xmax);
	TF1 *ff = (TF1*)h->GetListOfFunctions()->FindObject(f->GetName());
	ff->Print();
	double mean = ff->GetParameter(1);
	double stddev = ff->GetParameter(2);

	printf("[%s] #mu=%f  #sigma=%f\n", h->GetName(), mean, stddev);

	tex->DrawLatex(label_x, mu_labeyl_y, TString::Format("#mu = %.4f", mean));
	tex->DrawLatex(label_x, sigma_label_y, TString::Format("#sigma = %.2f", stddev));
}

void Histogram_adjust(TH1D *Hist,double title_size=-100,double label_size=-100,
double x_title_offset=-100,double y_title_offset=-100,
double x_label_offset=-100,double y_label_offset=-100,
 int x_divisions=-100,int y_divisions=-100)
{
    if(title_size!=-100)Hist->GetXaxis()->SetTitleSize(title_size);
    if(title_size!=-100)Hist->GetYaxis()->SetTitleSize(title_size);
    
	if(label_size!=-100)Hist->GetXaxis()->SetLabelSize(label_size);
    if(label_size!=-100)Hist->GetYaxis()->SetLabelSize(label_size);

	if(x_title_offset!=-100)Hist->GetXaxis()->SetTitleOffset(x_title_offset);
    if(y_title_offset!=-100)Hist->GetYaxis()->SetTitleOffset(y_title_offset);
   
    if(x_label_offset!=-100)Hist->GetXaxis()->SetLabelOffset(x_label_offset);
    if(y_label_offset!=-100)Hist->GetYaxis()->SetLabelOffset(y_label_offset);

   if(x_divisions!=-100)Hist->GetXaxis()->SetNdivisions(x_divisions);   
   if(y_divisions!=-100)Hist->GetYaxis()->SetNdivisions(y_divisions);
   Hist->GetXaxis()->CenterTitle(true);
    Hist->GetYaxis()->CenterTitle(true);
    Hist->GetXaxis()->SetNdivisions(5);
    Hist->SetLineWidth(2);
}
void Histogram_adjust(TH1F *Hist,double title_size=-100,double label_size=-100,
double x_title_offset=-100,double y_title_offset=-100,
double x_label_offset=-100,double y_label_offset=-100,
 int x_divisions=-100,int y_divisions=-100)
{
    if(title_size!=-100)Hist->GetXaxis()->SetTitleSize(title_size);
    if(title_size!=-100)Hist->GetYaxis()->SetTitleSize(title_size);
    
	if(label_size!=-100)Hist->GetXaxis()->SetLabelSize(label_size);
    if(label_size!=-100)Hist->GetYaxis()->SetLabelSize(label_size);

	if(x_title_offset!=-100)Hist->GetXaxis()->SetTitleOffset(x_title_offset);
    if(y_title_offset!=-100)Hist->GetYaxis()->SetTitleOffset(y_title_offset);
   
    if(x_label_offset!=-100)Hist->GetXaxis()->SetLabelOffset(x_label_offset);
    if(y_label_offset!=-100)Hist->GetYaxis()->SetLabelOffset(y_label_offset);

   if(x_divisions!=-100)Hist->GetXaxis()->SetNdivisions(x_divisions);   
   if(y_divisions!=-100)Hist->GetYaxis()->SetNdivisions(y_divisions);
   Hist->GetXaxis()->CenterTitle(true);
    Hist->GetYaxis()->CenterTitle(true);
    Hist->SetStats(kFALSE);
    Hist->GetXaxis()->SetNdivisions(5);
    Hist->SetLineWidth(2);
}
void Histogram_adjust(TH1I *Hist,double title_size=-100,double label_size=-100,
double x_title_offset=-100,double y_title_offset=-100,
double x_label_offset=-100,double y_label_offset=-100,
 int x_divisions=-100,int y_divisions=-100)
{
    if(title_size!=-100)Hist->GetXaxis()->SetTitleSize(title_size);
    if(title_size!=-100)Hist->GetYaxis()->SetTitleSize(title_size);
    
	if(label_size!=-100)Hist->GetXaxis()->SetLabelSize(label_size);
    if(label_size!=-100)Hist->GetYaxis()->SetLabelSize(label_size);

	if(x_title_offset!=-100)Hist->GetXaxis()->SetTitleOffset(x_title_offset);
    if(y_title_offset!=-100)Hist->GetYaxis()->SetTitleOffset(y_title_offset);
   
    if(x_label_offset!=-100)Hist->GetXaxis()->SetLabelOffset(x_label_offset);
    if(y_label_offset!=-100)Hist->GetYaxis()->SetLabelOffset(y_label_offset);

   if(x_divisions!=-100)Hist->GetXaxis()->SetNdivisions(x_divisions);   
   if(y_divisions!=-100)Hist->GetYaxis()->SetNdivisions(y_divisions);
   Hist->GetXaxis()->CenterTitle(true);
    Hist->GetYaxis()->CenterTitle(true);
    Hist->SetStats(kFALSE);
    Hist->GetXaxis()->SetNdivisions(5);
    Hist->SetLineWidth(2);
}
void analysis_result_adjustment_mom() {

	//TFile *f = new TFile("Lambda_kinfit_ForwardHad_test_fast.root");
    TLatex * tex = new TLatex;
	tex->SetNDC(true);
   
    TFile *f = new TFile("output_ecal.root");


	//f->ls();
    //cout<<endl<<endl<<"2"<<endl<<endl;
    //EMC gen

    TH1F *hEPhtoton1EMC_gen=(TH1F*)f->Get("hEPhtoton1EMC_gen");
    TH1F *hEPhtoton2EMC_gen=(TH1F*)f->Get("hEPhtoton2EMC_gen");


    //TH1D* Iterations=(TH1D*)f->Get("Iterations");
    TH1I * hIterations =(TH1I*)f->Get("hIterations");/////////////
    TH1F * hChi2 = (TH1F*)f->Get("hChi2");/////////////
    TH1F * hProb = (TH1F*)f->Get("hProb");///////////
  
    //TH1D *hXi2_Refit =(TH1D*)f->Get("hXi2_Refit"); 
    //TH1D *hprob_Refit =(TH1D*)f->Get("hprob_Refit");
    
    TH1F * hPullEInvPhoton1Conv_mass= (TH1F*)f->Get("Pull_E_photon_1");////////
    TH1F * hPullThetaPhoton1Conv_mass= (TH1F*)f->Get("Pull_theta_photon_1");//////
    TH1F * hPullPhiPhoton1Conv_mass= (TH1F*)f->Get("Pull_phi_photon_1");//////
  

    TH1F * hPullEInvPhoton2Conv_mass= (TH1F*)f->Get("Pull_E_photon_2");//////////
    TH1F * hPullThetaPhoton2Conv_mass= (TH1F*)f->Get("Pull_theta_photon_2");//////
    TH1F * hPullPhiPhoton2Conv_mass= (TH1F*)f->Get("Pull_phi_photon_2");///////


    //photons before refit
    TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");////////
    TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");//////////

    //photons after refit
    TH1F * hEPhtoton1AfterRefit= (TH1F*)f->Get("hEPhtoton1AfterRefit");///////
    TH1F * hPPhtoton1AfterRefit= (TH1F*)f->Get("hPPhtoton1AfterRefit");////////////
    TH1F * hmassPhtoton1AfterRefit= (TH1F*)f->Get("hmassPhtoton1AfterRefit");
    TH1F * hthetaPhtoton1AfterRefit= (TH1F*)f->Get("hthetahtoton1AfterRefit");
    TH1F * hphiPhtoton1AfterRefit= (TH1F*)f->Get("hphiPhtoton1AfterRefit");


    TH1F * hEPhtoton2AfterRefit= (TH1F*)f->Get("hEPhtoton2AfterRefit");//////
    TH1F * hPPhtoton2AfterRefit= (TH1F*)f->Get("hPPhtoton2AfterRefit");///////
    TH1F * hmassPhtoton2AfterRefit= (TH1F*)f->Get("hmassPhtoton2AfterRefit");
    TH1F * hthetaPhtoton2AfterRefit= (TH1F*)f->Get("hthetaphtoton2AfterRefit");
    TH1F * hphiPhtoton2AfterRefit= (TH1F*)f->Get("hphiPhtoton2AfterRefit");

   


    TH1F * hEETAAfterRefit= (TH1F*)f->Get("hEETAAfterRefit");
    TH1F * hPETAAfterRefit= (TH1F*)f->Get("hPETAAfterRefit");
    TH1F * hmassETAAfterRefit= (TH1F*)f->Get("hmassETAAfterRefit");
    TH1F * hthetaETAAfterRefit= (TH1F*)f->Get("hthetahtoton2AfterRefit");
    TH1F * hphiETAAfterRefit= (TH1F*)f->Get("hphiETAAfterRefit");
    TH1F * phiETABeforeRefit= (TH1F*)f->Get("hphiETABeforeRefit");



    TH1F * hEETABeforeRefit_selection=(TH1F*) f->Get("hEETABeforeRefit_selection");
    TH1F * hEETABeforeRefit= (TH1F*) f->Get("hEETABeforeRefit");



    //mix
    TH1F * hPullEInvPhoton1Conv_mass_mix=(TH1F*) f->Get("Pull_E_photon_1_mix");
    TH1F * hPullThetaPhoton1Conv_mass_mix=(TH1F*) f->Get("Pull_theta_photon_1_mix");
    TH1F * hPullPhiPhoton1Conv_mass_mix=(TH1F*) f->Get("Pull_phi_photon_1_mix");
    TH1F * hEPhtotonAfterRefit_mix=(TH1F*) f->Get("hEPhtotonAfterRefit_mix");
    TH1F * hPPhtotonAfterRefit_mix=(TH1F*) f->Get("hPPhtotonAfterRefit_mix");
    TH1F * hmassPhtotonAfterRefit_mix=(TH1F*) f->Get("hmassPhtotonAfterRefit_mix");
    TH1F * hthetaPhtotonAfterRefit_mix=(TH1F*) f->Get("hthetaphtotonAfterRefit_mix");
    TH1F * hphiPhtotonAfterRefit_mix=(TH1F*) f->Get("hphiPhtotonAfterRefit_mix");

    TH1F * hEPhtonEMC_gen=(TH1F*) f->Get("hEPhtonEMC_gen");
    TH1F * hEPhtotonbeforeRefitcomb=(TH1F*) f->Get("Before Energy Photon comb");

    TH1F * hbeta=(TH1F*) f->Get("hbeta");
    TH1F * hbeta150=(TH1F*) f->Get("hbeta150");


    TH1F * htheta_GEANT=(TH1F*) f->Get("htheta_GEANT");
    TH1F * htheta_GEANT150=(TH1F*) f->Get("htheta_GEANT150");
    
    TH1F * htheta_reco=(TH1F*) f->Get("htheta_reco");
    TH1F * htheta_reco150=(TH1F*) f->Get("htheta_reco150");

       
   
    TH1F * hResidua_Refit_Reco_photon1_E=(TH1F*) f->Get("hResidua_Refit_Reco_photon1_E");
    TH1F * hResidua_Refit_Reco_photon1_theta=(TH1F*) f->Get("hResidua_Refit_Reco_photon1_theta");
    TH1F * hResidua_Refit_Reco_photon1_phi=(TH1F*) f->Get("hResidua_Refit_Reco_photon1_phi");
   

    TH1F * hResidua_True_Refit_photon1_E=(TH1F*) f->Get("hResidua_True_Refit_photon1_E");
    TH1F * hResidua_True_Refit_photon1_theta=(TH1F*) f->Get("hResidua_True_Refit_photon1_theta");
    TH1F * hResidua_True_Refit_photon1_phi=(TH1F*) f->Get("hResidua_True_Refit_photon1_phi");

    TH1F * hResidua_True_Reco_photon1_E=(TH1F*) f->Get("hResidua_True_Reco_photon1_E");
    TH1F * hResidua_True_Reco_photon1_theta=(TH1F*) f->Get("hResidua_True_Reco_photon1_theta");
    TH1F * hResidua_True_Reco_photon1_phi=(TH1F*) f->Get("hResidua_True_Reco_photon1_phi");
   // TH1F * hEETABeforeRefit=(TH1F*) f->Get("hEETABeforeRefit");
    TH1F * hthetaETABeforeRefit=(TH1F*) f->Get("hthetaETABeforeReft");

    TH1F * hinvM_gg=(TH1F*) f->Get("hinvM_gg");
   
    TH1F * hphiETABeforeRefit=(TH1F*) f->Get("hphiETABeforeRefit");
    //TH1F *  hphiETAAfterRefit=(TH1F*) f->Get("hphiETAAfterRefit");


    TH1F * EtaP_m=(TH1F*) f->Get("EtaP_m");
    TH1F * EtaP_theta=(TH1F*) f->Get("EtaP_theta");
    TH1F * EtaP_phi=(TH1F*) f->Get("EtaP_phi");


    TH1F * Eta_m=(TH1F*) f->Get("Eta_m");
    TH1F * Eta_theta=(TH1F*) f->Get("Eta_theta");
    TH1F * Eta_phi=(TH1F*) f->Get("Eta_phi");


    TH1F * Eta_m_before=(TH1F*) f->Get("Eta_m_before");
    TH1F * Eta_theta_before=(TH1F*) f->Get("Eta_theta_before");
    TH1F * Eta_phi_before=(TH1F*) f->Get("Eta_phi_before");
    
















    TCanvas* a = new TCanvas("a", "a", 200,100,1920, 1080); //utworz kanwe
    a->Divide(3,1);	

	a->cd(1);
    gStyle->SetOptTitle(0);//gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hIterations->SetXTitle("Iterations");
	hIterations->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hIterations,0.11,0.09,0.4,0.7,-0.03,-100);
    hIterations->GetYaxis()->SetNdivisions(10);
    hIterations->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hIterations->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hIterations->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
    hIterations->Draw();
	a->Update();
   

    a->cd(2);
	gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.2);
	hChi2->SetXTitle("#chi^2");
	hChi2->SetYTitle("Counts");
    

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hChi2,0.11,0.09,0.4,0.7,-0.03,-100);
	gPad->SetLogy();
    	 hChi2->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);

    hChi2->Draw();
	a->Update();

    a->cd(3);
    gStyle->SetOptTitle(0);
    //gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.2);
	hProb->SetXTitle("Probability");
	hProb->SetYTitle("Counts");
    

	
	Histogram_adjust(hProb,0.11,0.09,0.4,0.7,-0.03,-100);
	
	gPad->SetLogy();
    hProb->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);
    hProb->Draw();
    TLine line;
    line.SetLineColor(kGreen);
    line.SetLineWidth(2);
    line.DrawLine(0.01,0,0.01,100000);
    TText *t = new TText();
    t->SetTextColor(kGreen+2);
    t->SetTextFont(43);
    
    t->SetTextSize(40);
    t->DrawText(0.05,15000,"Probability");
    t->DrawText(0.05,10000,"cut: P>0.01");
     

	a->Update();

a->SaveAs("Iterations.png");

    TCanvas* aa = new TCanvas("aa", "aa", 200,100,1920, 1080); //utworz kanwe
    aa->Divide(3,1);	

	aa->cd(1);
     

   
    gStyle->SetOptTitle(0);//gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hEPhtotonbeforeRefitcomb->SetXTitle("#gamma Energy [MeV]");
	hEPhtotonbeforeRefitcomb->SetYTitle("Counts");
Histogram_adjust(hEPhtotonbeforeRefitcomb,0.11,0.09,0.4,0.4,-0.03,-100);
Histogram_adjust(hEPhtotonAfterRefit_mix,0.11,0.09,0.4,0.4,-0.03,-100);

	//line.DrawLine(0,400,8000,400);
	//Histogram_adjust(hEPhtonEMC_gen,0.11,0.07,0.45,0.4,-0.02,-100);
   // hEPhtonEMC_gen->GetYaxis()->SetTitleOffset(0.25);
   // hIterations->GetYaxis()->SetNdivisions(10);
    hEPhtotonbeforeRefitcomb->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hEPhtotonbeforeRefitcomb->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
   // hEPhtonEMC_gen->Draw();
    hEPhtotonbeforeRefitcomb->Draw();
    gStyle->SetHistLineStyle(1);
    hEPhtotonAfterRefit_mix->Scale(2);
    hEPhtotonAfterRefit_mix->Draw("same,HIST");
    auto legendaa = new TLegend(0.4,0.7,0.9,0.9);
    legendaa->SetTextSize(0.05);
    //legendaa->AddEntry(hEPhtoton1beforeRefit,"All neutral-ideal","L");
    hEPhtotonbeforeRefitcomb->SetLineColor(kRed);
    legendaa->AddEntry(hEPhtotonbeforeRefitcomb,"Reconstruction","L");
    hEPhtotonAfterRefit_mix->SetLineColor(kGreen);
    legendaa->AddEntry(hEPhtotonAfterRefit_mix,"Refit photons x 2","L");



    legendaa->Draw();
	aa->Update();
   

    aa->cd(2);
	//gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.2);
	hbeta->SetXTitle("#beta");
	hbeta->SetYTitle("Counts");
    //=(TH1F*) f->Get("hbeta");
   // TH1F * hbeta150=

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hbeta,0.11,0.07,0.45,0.4,-0.02,-100);
    Histogram_adjust(hbeta150,0.11,0.07,0.45,0.4,-0.02,-100);
    hbeta->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
    hbeta->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hbeta->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
    hbeta->GetXaxis()->SetRangeUser(0.5,1.5);
    hbeta->Draw();
    hbeta150->Draw("Same");

       auto legendaa2 = new TLegend(0.6,0.7,0.9,0.9);
    legendaa2->SetTextSize(0.045);
    legendaa2->AddEntry(hbeta,"All neutral","L");
    hbeta->SetLineColor(kRed);
    legendaa2->AddEntry(hbeta150,"E>150MeV","L");
    hEPhtotonAfterRefit_mix->SetLineColor(kGreen);
   // legendaa2->AddEntry(hEPhtotonAfterRefit_mix,"Refit photons x 20","L");
    legendaa2->Draw();
    

	aa->Update();

    aa->cd(3);
     auto legendaa3 = new TLegend(0.50,0.7,0.9,0.9);
    
    legendaa3->SetTextSize(0.045);
 
 
    gStyle->SetOptTitle(0);
    //gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.2);
	htheta_reco->SetXTitle("#theta [#circ]");
	htheta_reco->SetYTitle("Counts");
    
   
    //TH1F * htheta_reco=(TH1F*) f->Get("htheta_reco");
    //TH1F * htheta_reco150=(TH1F*) f->Get("htheta_reco150");
	//line.DrawLine(0,400,8000,400);
    
	//Histogram_adjust(htheta_GEANT,0.11,0.07,0.45,0.4,-0.02,-100);
   // Histogram_adjust(htheta_GEANT150,0.07,0.05,-100,0.3,-100,-100);   
    Histogram_adjust(htheta_reco,0.11,0.07,0.45,0.4,-0.02,-100); 
     htheta_reco->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    htheta_reco->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   // htheta_GEANT->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
   // htheta_GEANT->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
 //   htheta_GEANT->GetXaxis()->SetNdivisions(5);
  //  htheta_GEANT->GetXaxis()->SetRangeUser(0,180);
 //  htheta_GEANT->Draw();
  //    legendaa3->AddEntry(htheta_GEANT,"ideal","L");
      //htheta_reco->Scale(5);
      htheta_reco->RebinX(4);
      htheta_reco->GetXaxis()->SetRangeUser(0,100);
       htheta_reco->Draw("same,hist");
    htheta_reco->SetLineColor(kRed);
    legendaa3->AddEntry(htheta_reco,"Reconstruction","L");
   // htheta_GEANT150->SetLineColor(kGreen);
    //htheta_GEANT150->Draw("same");
  //  legendaa3->AddEntry(htheta_GEANT150,"Geant E>150","L");
    legendaa3->Draw();
    

	aa->Update();

    aa->SaveAs("Photon_distr.png");



     TCanvas* ea = new TCanvas("ea", "ea", 200,100,1920, 1080);  //utworz kanwe
    ea->Divide(2,2);	
    ea->cd(1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    
    hEETABeforeRefit->SetXTitle("#ETA energy [MeV]");
    hEETABeforeRefit->SetYTitle("Counts");
    
	Histogram_adjust(hEETAAfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    hEETABeforeRefit->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hEETABeforeRefit->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hEETABeforeRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hEETABeforeRefit->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
   // hEETAAfterRefit->GetXaxis()->SetNdivisions(510,kTRUE);   
   Histogram_adjust(hEETABeforeRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
   
   hEETABeforeRefit->GetXaxis()->SetRangeUser(100,3000);

    hEETABeforeRefit->Draw();
    hEETAAfterRefit->Draw("same");
    hEETABeforeRefit->SetLineColor(kRed);
     auto legendmass1 = new TLegend(0.6,0.7,0.9,0.9);

    legendmass1->SetTextSize(0.05);
    
    legendmass1->AddEntry(hEETABeforeRefit,"before refit","L");
    legendmass1->AddEntry(hEETAAfterRefit,"Refit","L");
    
    legendmass1->Draw();
    ea->Update();
  
    ea->cd(2);

    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hthetaETABeforeRefit->SetXTitle("#ETA #theta [#circ]");
    hthetaETABeforeRefit->SetYTitle("Counts");
    hthetaETABeforeRefit->GetXaxis()->SetRangeUser(0,50);
	Histogram_adjust(hthetaETABeforeRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    hthetaETABeforeRefit->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hthetaETABeforeRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hthetaETABeforeRefit->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    //hthetaETABeforeRefit->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    //hthetaETABeforeRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hthetaETABeforeRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
   
    
    
    hthetaETABeforeRefit->Draw();
    hthetaETAAfterRefit->Draw("same");
    hthetaETABeforeRefit->SetLineColor(kRed);
     auto legendmass2 = new TLegend(0.7,0.7,0.9,0.9);

    legendmass2->SetTextSize(0.05);
    
    legendmass2->AddEntry(hthetaETAAfterRefit,"Refit","L");
    legendmass2->AddEntry(hthetaETABeforeRefit,"Before refit","L");
    
    legendmass2->Draw();

    ea->Update();


 
    ea->cd(3);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hmassETAAfterRefit->SetXTitle("#ETA mass [MeV/c^2]");
    hmassETAAfterRefit->SetYTitle("Counts");
    hmassETAAfterRefit->GetXaxis()->SetRangeUser(134.9765,134.9775);
	Histogram_adjust(hmassETAAfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    hmassETAAfterRefit->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hmassETAAfterRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hmassETAAfterRefit->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
   // hmassETAAfterRefit->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);
    hmassETAAfterRefit->GetXaxis()->SetNdivisions(3);  
  
    hmassETAAfterRefit->Draw();
    ea->Update();

    ea->cd(4);
   
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hinvM_gg->SetXTitle("#ETA mass [MeV/c^2]");
    hinvM_gg->SetYTitle("Counts");
    hinvM_gg->GetXaxis()->SetRangeUser(100,200);
    hinvM_gg->RebinX(50000);
    hinvM_gg->SetLineColor(kRed);
	Histogram_adjust(hmassETAAfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    Histogram_adjust(hinvM_gg,0.11,0.07,0.6,0.2,-0.001,-100);
    hinvM_gg->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1); 
    hinvM_gg->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hinvM_gg->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    auto legendmass = new TLegend(0.7,0.7,0.9,0.9);

    legendmass->SetTextSize(0.05);
    legendmass->AddEntry(hinvM_gg,"without Refit","L");
    legendmass->AddEntry(hmassETAAfterRefit,"Refit","L");
    
    hinvM_gg->Draw();
    hmassETAAfterRefit->Draw("Same");
    legendmass->Draw();
     ea->SaveAs("ETAmass_contr.png");


TCanvas* Photon_pulls = new TCanvas("Photon_pulls", "Photon_pulls", 200,100,1920, 1080);  //utworz kanwe
    Photon_pulls->Divide(3,1);	
//TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");
  //  TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");
	Photon_pulls->cd(1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);

    hPullEInvPhoton1Conv_mass_mix->SetXTitle("Pull 1/E");
    hPullEInvPhoton1Conv_mass_mix->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullEInvPhoton1Conv_mass_mix,0.11,0.09,0.4,0.7,-0.03,-100);
     hPullEInvPhoton1Conv_mass_mix->RebinX(5);
   // hPullEInvPhoton1Conv_mass_mix->GetXaxis()->SetRangeUser(-0.75,0.75);
     TF1 * f_res = new TF1("f_res", "gaus(0)", -30, 30);
    f_res->SetParameters(300, 0., 0.2, 3);
    hPullEInvPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hPullEInvPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hPullEInvPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hPullEInvPhoton1Conv_mass_mix->Draw();
	fit_hist(hPullEInvPhoton1Conv_mass_mix,f_res,-3,3,tex,0.65,0.8,0.75);
	//gPad->SetLogy();
    
	Photon_pulls->Update();
   

    Photon_pulls->cd(2);

    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullThetaPhoton1Conv_mass_mix->SetXTitle("Pull #theta");
    hPullThetaPhoton1Conv_mass_mix->SetYTitle("Counts");
    hPullThetaPhoton1Conv_mass_mix->RebinX(5);
    //hPullThetaPhoton1Conv_mass_mix->GetXaxis()->SetRangeUser(-0.75,0.75);
	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullThetaPhoton1Conv_mass_mix,0.11,0.09,0.4,0.7,-0.03,-100);

    hPullThetaPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hPullThetaPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hPullThetaPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    f_res->SetParameters(300, 0., 1.5, 3);
    
    hPullThetaPhoton1Conv_mass_mix->Draw();
	fit_hist(hPullThetaPhoton1Conv_mass_mix,f_res,-3.0,3.0,tex,0.66,0.8,0.75);
  
	Photon_pulls->Update();
    
    Photon_pulls->cd(3);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullPhiPhoton1Conv_mass_mix->RebinX(5);
    //hPullPhiPhoton1Conv_mass_mix->GetXaxis()->SetRangeUser(-0.75,0.75);
    hPullPhiPhoton1Conv_mass_mix->SetXTitle("Pull #phi");
    hPullPhiPhoton1Conv_mass_mix->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullPhiPhoton1Conv_mass_mix,0.11,0.09,0.4,0.7,-0.03,-100);
    hPullPhiPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hPullPhiPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hPullPhiPhoton1Conv_mass_mix->Draw();
    f_res->SetParameters(300, 0., 1.5, 3);
    	fit_hist(hPullPhiPhoton1Conv_mass_mix,f_res,-3.0,3.0,tex,0.66,0.8,0.75);
	Photon_pulls->Update();
Photon_pulls->SaveAs("Photon_pulls.png");


TCanvas* phi = new TCanvas("phi", "phi", 200,100,1920, 1080);  //utworz kanwe
phi->Divide(1,1);
phi->cd();	
gPad->SetBottomMargin(0.2);
	gPad->SetLeftMargin(0.25);
   // hPullPhiPhoton1Conv_mass_mix->RebinX(5);
    //hPullPhiPhoton1Conv_mass_mix->GetXaxis()->SetRangeUser(-0.75,0.75);
    hphiETABeforeRefit->SetXTitle("#phi [#circ]");
    hphiETABeforeRefit->SetYTitle("Counts");
     auto legendphi = new TLegend(0.6,0.7,0.9,0.9);
    hphiETABeforeRefit->SetLineColor(kRed);
    legendphi->SetTextSize(0.05);
    legendphi->AddEntry(hphiETABeforeRefit,"Before refit","L");
    legendphi->AddEntry(hphiETAAfterRefit,"After refit","L");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hphiETABeforeRefit,0.11,0.08,0.7,0.7,0.01,-100);
    Histogram_adjust(hphiETAAfterRefit,0.11,0.08,0.7,0.7,0.01,-100);
phiETABeforeRefit->RebinX(5);
hphiETAAfterRefit->RebinX(5);
hphiETABeforeRefit->GetYaxis()->SetRangeUser(0,2000);
hphiETABeforeRefit->Draw();
hphiETAAfterRefit->Draw("same,hist");
legendphi->Draw();
phi->SaveAs("phi.png");



TCanvas* Eta = new TCanvas("Eta", "Eta", 200,100,1920, 1080);  //utworz kanwe
Eta->Divide(3,2);
Eta->cd(1);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
hmassETAAfterRefit->SetXTitle("mass [#MeV]");
hmassETAAfterRefit->SetYTitle("Counts");
Histogram_adjust(hmassETAAfterRefit,0.11,0.08,0.7,0.7,0.01,-100);
hmassETAAfterRefit->RebinX(250);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
hmassETAAfterRefit->Draw();

Eta->cd(2);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
hthetaETABeforeRefit->SetXTitle("#theta [#circ]");
hthetaETABeforeRefit->SetYTitle("Counts");
Histogram_adjust(hthetaETABeforeRefit,0.11,0.08,0.7,0.7,0.01,-100);
hthetaETABeforeRefit->RebinX(5);
hthetaETABeforeRefit->GetXaxis()->SetRangeUser(0,80);
hthetaETABeforeRefit->Draw();

Eta->cd(3);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
Eta_phi->SetXTitle("#phi [#circ]");
Eta_phi->SetYTitle("Counts");
Histogram_adjust(Eta_phi,0.11,0.08,0.7,0.7,0.01,-100);
Eta_phi->RebinX(5);
Eta_phi->GetXaxis()->SetRangeUser(-200,200);
hphiETABeforeRefit->Draw();

Eta->cd(4);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
Eta_m_before->SetXTitle("mass without refit [#MeV]");
Eta_m_before->SetYTitle("Counts");
Histogram_adjust(Eta_m_before,0.11,0.08,0.7,0.7,0.01,-100);
Eta_m_before->RebinX(250);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
Eta_m_before->Draw();

Eta->cd(5);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
Eta_theta_before->SetXTitle("#theta without refit [#circ]");
Eta_theta_before->SetYTitle("Counts");
Histogram_adjust(Eta_theta_before,0.11,0.08,0.7,0.7,0.01,-100);
Eta_theta_before->RebinX(250);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
Eta_theta_before->Draw();

Eta->cd(6);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
Eta_phi_before->SetXTitle("#phi without refit [#circ]");
Eta_phi_before->SetYTitle("Counts");
Histogram_adjust(Eta_phi_before,0.11,0.08,0.7,0.7,0.01,-100);
Eta_phi_before->RebinX(250);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
Eta_phi_before->Draw();

Eta->SaveAs("Eta.png");
/*
Eta->cd(6);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
Eta_phi_before->SetXTitle("#phi without refit [#circ]");
Eta_phi_before->SetYTitle("Counts");
Histogram_adjust(Eta_phi_before,0.11,0.08,0.7,0.7,0.01,-100);
Eta_phi_before->RebinX(50);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
Eta_phi_before->Draw();
Eta->SaveAs("Eta.png");
*/



TCanvas* EtaP = new TCanvas("EtaP", "EtaP", 200,100,1920, 1080);  //utworz kanwe
EtaP->Divide(3,1);
EtaP->cd(1);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
EtaP_m->SetXTitle("mass [#MeV]");
EtaP_m->SetYTitle("Counts");
Histogram_adjust(Eta_m,0.11,0.08,0.7,0.7,0.01,-100);
EtaP_m->RebinX(2000);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
EtaP_m->Draw();

EtaP->cd(2);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
EtaP_theta->SetXTitle("#theta [#circ]");
EtaP_theta->SetYTitle("Counts");
Histogram_adjust(Eta_theta,0.11,0.08,0.7,0.7,0.01,-100);
EtaP_theta->RebinX(2000);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
EtaP_theta->Draw();

EtaP->cd(3);	
gPad->SetBottomMargin(0.2);
gPad->SetLeftMargin(0.25);
EtaP_phi->SetXTitle("#phi [#circ]");
EtaP_phi->SetYTitle("Counts");
Histogram_adjust(EtaP_phi,0.11,0.08,0.7,0.7,0.01,-100);
EtaP_phi->RebinX(2000);
//Eta_m->GetXaxis()->SetRangeUser(0,2000);
EtaP_phi->Draw();

EtaP->SaveAs("EtaP.png");

}
