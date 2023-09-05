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

   


    TH1F * hEPi0AfterRefit= (TH1F*)f->Get("hEPi0AfterRefit");
    TH1F * hPPi0AfterRefit= (TH1F*)f->Get("hPPi0AfterRefit");
    TH1F * hmassPi0AfterRefit= (TH1F*)f->Get("hmassPi0AfterRefit");
    TH1F * hthetaPi0AfterRefit= (TH1F*)f->Get("hthetahtoton2AfterRefit");
    TH1F * hphiPi0AfterRefit= (TH1F*)f->Get("hphiPi0AfterRefit");



    TH1F * hEPi0BeforeRefit_selection=(TH1F*) f->Get("hEPi0BeforeRefit_selection");
    TH1F * hEPi0BeforeRefit= (TH1F*) f->Get("hEPi0BeforeRefit");



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
   // TH1F * hEPi0BeforeRefit=(TH1F*) f->Get("hEPi0BeforeRefit");
    TH1F * hthetaPi0BeforeRefit=(TH1F*) f->Get("hthetaPi0BeforeReft");

    TH1F * hinvM_gg=(TH1F*) f->Get("hinvM_gg");
   
    TH1F * hphiPi0BeforeRefit=(TH1F*) f->Get("hphiPi0BeforeReft");


   // TH1F * hResidua_True_Reco_photon1_phi=(TH1F*) f->Get("hResidua_True_Reco_photon1_phi");
   // TH1F * hResidua_True_Refit_photon1_phi=(TH1F*) f->Get("hResidua_True_Refit_photon1_phi");
   // TH1F * hResidua_Refit_Reco_photon1_phi=(TH1F*) f->Get("hResidua_Refit_Reco_photon1_phi");
    //TH1F *  hphiPi0AfterRefit=(TH1F*) f->Get("hphiPi0AfterRefit");


    
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
    hIterations->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
     hIterations->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
    //hIterations->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
	
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
    
    t->SetTextSize(25);
    t->DrawText(0.05,15000,"Probability");
    t->DrawText(0.05,12000,"cut: P>0.01");
     

	a->Update();

a->SaveAs("Iterations.png");

    TCanvas* aa = new TCanvas("aa", "aa", 200,100,1920, 1080); //utworz kanwe
    aa->Divide(3,1);	

	aa->cd(1);
     

   
    gStyle->SetOptTitle(0);//gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hEPhtonEMC_gen->SetXTitle("#gamma Energy [MeV]");
	hEPhtonEMC_gen->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hEPhtonEMC_gen,0.11,0.07,0.45,0.4,-0.02,-100);
    hEPhtonEMC_gen->GetYaxis()->SetTitleOffset(0.25);
   // hIterations->GetYaxis()->SetNdivisions(10);
    hEPhtonEMC_gen->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hEPhtonEMC_gen->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
    hEPhtonEMC_gen->Draw();
    hEPhtotonbeforeRefitcomb->Draw("same");
    gStyle->SetHistLineStyle(1);
    hEPhtotonAfterRefit_mix->Scale(20);
    hEPhtotonAfterRefit_mix->Draw("same,HIST");
    auto legendaa = new TLegend(0.35,0.55,0.9,0.9);
    legendaa->SetTextSize(0.05);
    legendaa->AddEntry(hEPhtoton1beforeRefit,"All neutral-ideal","L");
    hEPhtotonbeforeRefitcomb->SetLineColor(kRed);
    legendaa->AddEntry(hEPhtotonbeforeRefitcomb,"Reconstruction","L");
    hEPhtotonAfterRefit_mix->SetLineColor(kGreen);
    legendaa->AddEntry(hEPhtotonAfterRefit_mix,"Refit photons x 20","L");



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
    hbeta->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
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
    //hEPhtotonAfterRefit_mix->SetLineColor(kGreen);
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
	htheta_GEANT->SetXTitle("#theta [#circ]");
	htheta_GEANT->SetYTitle("Counts");
    
   
    //TH1F * htheta_reco=(TH1F*) f->Get("htheta_reco");
    //TH1F * htheta_reco150=(TH1F*) f->Get("htheta_reco150");
	//line.DrawLine(0,400,8000,400);
    
	Histogram_adjust(htheta_GEANT,0.11,0.07,0.45,0.4,-0.02,-100);
   // Histogram_adjust(htheta_GEANT150,0.07,0.05,-100,0.3,-100,-100);   
    Histogram_adjust(htheta_reco,0.11,0.07,0.45,0.4,-0.02,-100); 
    htheta_GEANT->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    htheta_GEANT->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
    htheta_GEANT->GetXaxis()->SetNdivisions(5);
    htheta_GEANT->GetXaxis()->SetRangeUser(0,180);
   htheta_GEANT->Draw();
      legendaa3->AddEntry(htheta_GEANT,"ideal","L");
      htheta_reco->Scale(5);
      htheta_reco->RebinX(4);
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
    
    hEPi0BeforeRefit->SetXTitle("#pi0 energy [MeV]");
    hEPi0BeforeRefit->SetYTitle("Counts");
    hEPi0BeforeRefit->GetXaxis()->SetRangeUser(100,3000);
	Histogram_adjust(hEPi0AfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    hEPi0BeforeRefit->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hEPi0BeforeRefit->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hEPi0BeforeRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   // hEPi0AfterRefit->GetXaxis()->SetNdivisions(510,kTRUE);   
   Histogram_adjust(hEPi0BeforeRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
   
   

    hEPi0BeforeRefit->Draw();
    hEPi0AfterRefit->Draw("same");
    hEPi0BeforeRefit->SetLineColor(kRed);
     auto legendmass1 = new TLegend(0.6,0.7,0.9,0.9);

    legendmass1->SetTextSize(0.05);
    
    legendmass1->AddEntry(hEPi0BeforeRefit,"before refit","L");
    legendmass1->AddEntry(hEPi0AfterRefit,"Refit","L");
    
    legendmass1->Draw();
    ea->Update();
  
    ea->cd(2);

    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hthetaPi0BeforeRefit->SetXTitle("#pi0 #theta [#circ]");
    hthetaPi0BeforeRefit->SetYTitle("Counts");
    hthetaPi0AfterRefit->GetXaxis()->SetRangeUser(0,50);
	Histogram_adjust(hthetaPi0AfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    hthetaPi0BeforeRefit->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hthetaPi0BeforeRefit->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
   // hthetaPi0BeforeRefit->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
   // hthetaPi0BeforeRefit->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    hthetaPi0BeforeRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hthetaPi0BeforeRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
   
    
    hthetaPi0BeforeRefit->GetXaxis()->SetRangeUser(0,50);
    hthetaPi0BeforeRefit->Draw();
    hthetaPi0AfterRefit->Draw("same");
    hthetaPi0BeforeRefit->SetLineColor(kRed);
     auto legendmass2 = new TLegend(0.7,0.7,0.9,0.9);

    legendmass2->SetTextSize(0.05);
    
    legendmass2->AddEntry(hthetaPi0AfterRefit,"Refit","L");
    legendmass2->AddEntry(hthetaPi0BeforeRefit,"Before refit","L");
    
    legendmass2->Draw();

    ea->Update();


 
    ea->cd(3);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hmassPi0AfterRefit->SetXTitle("#pi0 mass [MeV/c^2]");
    hmassPi0AfterRefit->SetYTitle("Counts");
    hmassPi0AfterRefit->GetXaxis()->SetRangeUser(134.97675,134.97685);
	Histogram_adjust(hmassPi0AfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    hmassPi0AfterRefit->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hmassPi0AfterRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hmassPi0AfterRefit->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hmassPi0AfterRefit->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hmassPi0AfterRefit->GetXaxis()->SetNdivisions(3);  
  
    hmassPi0AfterRefit->Draw();
    ea->Update();

    ea->cd(4);
   
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hinvM_gg->SetXTitle("#pi0 mass [MeV/c^2]");
    hinvM_gg->SetYTitle("Counts");
    hinvM_gg->GetXaxis()->SetRangeUser(100,200);
    hinvM_gg->RebinX(50000);
    hinvM_gg->SetLineColor(kRed);
	Histogram_adjust(hmassPi0AfterRefit,0.11,0.07,0.6,0.2,-0.001,-100); 
    Histogram_adjust(hinvM_gg,0.11,0.07,0.6,0.2,-0.001,-100); 
    hinvM_gg->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hinvM_gg->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
    hinvM_gg->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);
    auto legendmass = new TLegend(0.6,0.7,0.9,0.9);

    legendmass->SetTextSize(0.05);
    legendmass->AddEntry(hinvM_gg,"without Refit","L");
    legendmass->AddEntry(hmassPi0AfterRefit,"Refit","L");
    
    hinvM_gg->Draw();
    hmassPi0AfterRefit->Draw("Same");
    legendmass->Draw();
     ea->SaveAs("pi0mass_contr.png");

/*

TCanvas* b = new TCanvas("b", "b", 200,100,800, 800); //utworz kanwe
    b->Divide(3,1);	
//TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");
  //  TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");
	b->cd(1);
    auto legend = new TLegend(0.6,0.7,0.9,0.9);

    legend->SetTextSize(0.05);
    gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hEPhtoton1beforeRefit->SetXTitle("Energy [MeV]");
	hEPhtoton1beforeRefit->SetYTitle("Counts");
    hEPhtoton2beforeRefit->SetXTitle("Energy [MeV]");
    hEPhtoton2beforeRefit->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hEPhtoton1beforeRefit,0.07,0.07,1,1.5,-100,-100);
	
	//gPad->SetLogy();
    hEPhtoton2beforeRefit->Draw();
    hEPhtoton2beforeRefit->SetLineColor(kRed);
    hEPhtoton1beforeRefit->Draw("same");
    legend->AddEntry(hEPhtoton1beforeRefit,"Photon 1","L");
    legend->AddEntry(hEPhtoton2beforeRefit,"Photon 2","L");
    legend->Draw();
	b->Update();
   

    b->cd(2);
	gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.5);
    
    auto legend1 = new TLegend(0.6,0.7,0.9,0.9);

    legend1->SetTextSize(0.05);

    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hEPhtoton2AfterRefit->SetLineColor(kRed);
    legend1->AddEntry(hEPhtoton1AfterRefit,"Photon 1","L");
    legend1->AddEntry(hEPhtoton2AfterRefit,"Photon 2","L");
    Histogram_adjust(hEPhtoton1AfterRefit,0.07,0.07,1,1,-100,-100);
    Histogram_adjust(hEPhtoton2AfterRefit,0.07,0.07,1,1,-100,-100);
    
    hEPhtoton2AfterRefit->Draw();
    hEPhtoton1AfterRefit->Draw("same");
    legend1->Draw();
	b->Update();
    
    b->cd(3);
    
    gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.5);
    
    auto legend2 = new TLegend(0.6,0.7,0.9,0.9);

    legend2->SetTextSize(0.05);

    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hEPhtoton2EMC_gen->SetLineColor(kRed);
    legend2->AddEntry(hEPhtoton1EMC_gen,"Photon 1","L");
    legend2->AddEntry(hEPhtoton2EMC_gen,"Photon 2","L");
    Histogram_adjust(hEPhtoton1EMC_gen,0.07,0.07,1,1,-100,-100);
    Histogram_adjust(hEPhtoton2EMC_gen,0.07,0.07,1,1,-100,-100);
    
    hEPhtoton2EMC_gen->Draw();
    hEPhtoton1EMC_gen->Draw("same");
    legend2->Draw();
	b->Update();
	
	
	b->Update();

b->SaveAs("Photon_energy.png");



*/
/*
TCanvas* b_mix = new TCanvas("b_mix", "b_mix", 200,100,800, 800); //utworz kanwe
    b_mix->Divide(3,3);	
//TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");
  //  TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");
	b_mix->cd(1);
  
    //gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullEInvPhoton1Conv_mass_mix->SetXTitle("Pull 1/E ");
    hPullEInvPhoton1Conv_mass_mix->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullEInvPhoton1Conv_mass_mix,0.07,0.07,1,1.5,-100,-100);
	
	//gPad->SetLogy();
    hPullEInvPhoton1Conv_mass_mix->Draw();
	b_mix->Update();
   

    b_mix->cd(2);
    //gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullThetaPhoton1Conv_mass_mix->SetXTitle("Pull photon #theta ");
    hPullThetaPhoton1Conv_mass_mix->SetYTitle("Counts");
    Histogram_adjust(hPullThetaPhoton1Conv_mass_mix,0.07,0.07,1,1.5,-100,-100);
    hPullThetaPhoton1Conv_mass_mix->Draw();
	b_mix->Update();
    
    b_mix->cd(3);
    // gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullPhiPhoton1Conv_mass_mix->SetXTitle("Pull photon #phi ");
    hPullPhiPhoton1Conv_mass_mix->SetYTitle("Counts");
    Histogram_adjust(hPullPhiPhoton1Conv_mass_mix,0.07,0.07,1,1.5,-100,-100);
    hPullPhiPhoton1Conv_mass_mix->Draw();
	b_mix->Update();
	b_mix->cd(4);
   //  gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPPhtotonAfterRefit_mix->SetXTitle("Momentum photon [Mev/c]");
    hPPhtotonAfterRefit_mix->SetYTitle("Counts");
    Histogram_adjust(hPPhtotonAfterRefit_mix,0.07,0.07,1,1.5,-100,-100);
    hPPhtotonAfterRefit_mix->Draw();
	b_mix->Update();
    b_mix->cd(5);
   // gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hmassPhtotonAfterRefit_mix->SetXTitle("mass photon [Mev/c^2]");
    hmassPhtotonAfterRefit_mix->SetYTitle("Counts");
    Histogram_adjust(hmassPhtotonAfterRefit_mix,0.07,0.07,1,1.5,-100,-100);
    hmassPhtotonAfterRefit_mix->Draw();
	b_mix->Update();
    b_mix->cd(6);
   //  gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hthetaPhtotonAfterRefit_mix->SetXTitle("#theta photon [#circ]");
    hthetaPhtotonAfterRefit_mix->SetYTitle("Counts");
    Histogram_adjust(hthetaPhtotonAfterRefit_mix,0.07,0.07,1,1.5,-100,-100);
    hthetaPhtotonAfterRefit_mix->Draw();
	b_mix->Update();
    b_mix->cd(7);
   //  gStyle->SetOptTitle(1);
    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hphiPhtotonAfterRefit_mix->SetXTitle("#phi photon [#circ]");
    hphiPhtotonAfterRefit_mix->SetYTitle("Counts");
    Histogram_adjust(hphiPhtotonAfterRefit_mix,0.07,0.07,1,1.5,-100,-100);
   
    hphiPhtotonAfterRefit_mix->Draw();
	b_mix->Update();


b_mix->SaveAs("Photon_mix.png");
*/
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
    hPullEInvPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
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
    hPullPhiPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);
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
    hphiPi0BeforeRefit->SetXTitle("#phi [#circ]");
    hphiPi0BeforeRefit->SetYTitle("Counts");
     auto legendphi = new TLegend(0.6,0.7,0.9,0.9);
    hphiPi0BeforeRefit->SetLineColor(kRed);
    legendphi->SetTextSize(0.05);
    legendphi->AddEntry(hphiPi0BeforeRefit,"Before refit","L");
    legendphi->AddEntry(hphiPi0AfterRefit,"After refit","L");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hphiPi0BeforeRefit,0.11,0.08,0.7,0.7,0.01,-100);
    Histogram_adjust(hphiPi0AfterRefit,0.11,0.08,0.7,0.7,0.01,-100);
   //hphiPi0BeforeRefit->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   // hphiPi0BeforeRefit->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
  //  hPullPhiPhoton1Conv_mass_mix->Draw();
    //f_res->SetParameters(300, 0., 1.5, 3);
    //	fit_hist(hPullPhiPhoton1Conv_mass_mix,f_res,-3.0,3.0,tex,0.66,0.8,0.75);
//	Photon_pulls->Update();
//Photon_pulls->SaveAs("Photon_pulls.png");
hphiPi0BeforeRefit->RebinX(5);
hphiPi0AfterRefit->RebinX(5);
hphiPi0BeforeRefit->GetYaxis()->SetRangeUser(0,2000);
hphiPi0BeforeRefit->Draw();
hphiPi0AfterRefit->Draw("same,hist");
legendphi->Draw();
phi->SaveAs("phi.png");


///////////////////////////////////////////
///////////////////////////////////////////

    TCanvas* Residua = new TCanvas("Residua", "Residua", 200,100,1920, 1080);  //utworz kanwe
    Residua->Divide(3,2);	
    //TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");
    //  TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");
	Residua->cd(3);
   // auto legendaRes = new TLegend(0.35,0.55,0.9,0.9);
   // legendaRes->SetTextSize(0.05);
    //legendaRes->AddEntry(hResidua_Refit_Reco_photon1_E,"Refit-Reo","L");
    //hResidua_True_Refit_photon1_E->SetLineColor(kRed);
    //legendaRes->AddEntry(hResidua_True_Refit_photon1_E,"Ideal-Refit","L");
    //hResidua_True_Reco_photon1_E->SetLineColor(kGreen);
    //legendaRes->AddEntry(hResidua_True_Reco_photon1_E,"Ideal-Reco","L");
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hResidua_Refit_Reco_photon1_E->SetXTitle("#gamma E Refit-Reco [MeV]");
    hResidua_Refit_Reco_photon1_E->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hResidua_Refit_Reco_photon1_E,0.09,0.07,0.6,0.4,-0.001,-100);
    Histogram_adjust(hResidua_True_Refit_photon1_E,0.09,0.07,0.6,0.4,-0.001,-100);
    Histogram_adjust(hResidua_True_Reco_photon1_E,0.09,0.07,0.6,0.4,-0.001,-100);
  



    hResidua_Refit_Reco_photon1_E->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hResidua_Refit_Reco_photon1_E->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
   hResidua_Refit_Reco_photon1_E->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
   hResidua_Refit_Reco_photon1_E->GetXaxis()->SetRangeUser(-0.2,0.2);
   hResidua_Refit_Reco_photon1_E->Draw();
    //hResidua_True_Refit_photon1_E->Draw("");
    //hResidua_True_Reco_photon1_E->Draw("");
	//legendaRes->Draw();
	//gPad->SetLogy();
    
	Residua->Update();
    Residua->cd(2);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hResidua_True_Refit_photon1_E->SetXTitle("#gamma E Ideal-Refit [MeV]");
   hResidua_True_Refit_photon1_E->SetYTitle("Counts");
   hResidua_True_Refit_photon1_E->RebinX(50);
   hResidua_True_Refit_photon1_E->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
   hResidua_True_Refit_photon1_E->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   hResidua_True_Refit_photon1_E->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
   hResidua_True_Refit_photon1_E->Draw();
   Residua->cd(1);
       gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hResidua_True_Reco_photon1_E->SetXTitle("#gamma E Ideal-Reco [MeV]");
   hResidua_True_Reco_photon1_E->SetYTitle("Counts");
   hResidua_True_Reco_photon1_E->RebinX(50);
   hResidua_True_Reco_photon1_E->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   hResidua_True_Reco_photon1_E->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
   hResidua_True_Reco_photon1_E->Draw();
   ///////////////////////////////////////
   ///////////////////////////////////////
   ///////////////////////////////////////
  
	Residua->cd(6);
   
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hResidua_Refit_Reco_photon1_theta->SetXTitle("#gamma #theta Refit-Reco [MeV]");
    hResidua_Refit_Reco_photon1_theta->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
    hResidua_Refit_Reco_photon1_theta->RebinX(50);
	Histogram_adjust(hResidua_Refit_Reco_photon1_theta,0.09,0.07,0.6,0.4,-0.001,-100);
    Histogram_adjust(hResidua_True_Refit_photon1_theta,0.09,0.07,0.6,0.4,-0.001,-100);
    Histogram_adjust(hResidua_True_Reco_photon1_theta,0.09,0.07,0.6,0.4,-0.001,-100);
  



   hResidua_Refit_Reco_photon1_theta->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
   hResidua_Refit_Reco_photon1_theta->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
   hResidua_Refit_Reco_photon1_theta->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
   hResidua_Refit_Reco_photon1_theta->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   hResidua_Refit_Reco_photon1_theta->GetXaxis()->SetRangeUser(-50,50);
   //hResidua_Refit_Reco_photon1_theta->RebinX(2);
   hResidua_Refit_Reco_photon1_theta->Draw();
   
    
	Residua->Update();
    Residua->cd(5);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hResidua_True_Refit_photon1_theta->SetXTitle("#gamma #theta Ideal-Refit [MeV]");
   hResidua_True_Refit_photon1_theta->SetYTitle("Counts");
   hResidua_True_Refit_photon1_theta->RebinX(50);
   hResidua_True_Refit_photon1_theta->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   hResidua_True_Refit_photon1_theta->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
   hResidua_True_Refit_photon1_theta->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
   hResidua_True_Refit_photon1_theta->GetXaxis()->SetRangeUser(-50,50);
   hResidua_True_Refit_photon1_theta->Draw();
   Residua->cd(4);
       gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
   hResidua_True_Reco_photon1_theta->SetXTitle("#gamma #theta Ideal-Reco [MeV]");
   hResidua_True_Reco_photon1_theta->SetYTitle("Counts");
   hResidua_True_Reco_photon1_theta->RebinX(50);
  // hResidua_True_Reco_photon1_theta->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
   hResidua_True_Reco_photon1_theta->GetYaxis()->ChangeLabel(3,-1,0,-1,-1,-1,-1);
   hResidua_True_Reco_photon1_theta->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
   hResidua_True_Reco_photon1_theta->GetXaxis()->SetRangeUser(-50,50);
   hResidua_True_Reco_photon1_theta->Draw();
  Residua->SaveAs("Residua.png");
   

   /////////////////

   ////////////////
   TCanvas* ab = new TCanvas("ab", "ab", 200,100,1920, 1080); //utworz kanwe
    ab->Divide(3,1);	

	ab->cd(1);
    gStyle->SetOptTitle(0);//gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hResidua_True_Reco_photon1_phi->SetXTitle("#phi Ideal-Reco [#circ]");
	hResidua_True_Reco_photon1_phi->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hResidua_True_Reco_photon1_phi,0.11,0.07,0.4,0.7,-0.01,-100);
  //  hIterations->GetYaxis()->SetNdivisions(10);
    hResidua_True_Reco_photon1_phi->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    hResidua_True_Reco_photon1_phi->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hResidua_True_Reco_photon1_phi->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
    //hIterations->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
    hResidua_True_Reco_photon1_phi->RebinX(10);
    hResidua_True_Reco_photon1_phi->GetXaxis()->SetRangeUser(-200,200);
    hResidua_True_Reco_photon1_phi->Draw();
	ab->Update();
   

    ab->cd(2);
	gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.2);
	hResidua_True_Refit_photon1_phi->SetXTitle("#phi Ideal-Refit [#circ]");
	hResidua_True_Refit_photon1_phi->SetYTitle("Counts");
    

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hResidua_True_Refit_photon1_phi,0.11,0.07,0.4,0.7,-0.01,-100);
    hResidua_True_Refit_photon1_phi->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hResidua_True_Refit_photon1_phi->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hResidua_True_Refit_photon1_phi->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
	//gPad->SetLogy();
   // hResidua_True_Refit_photon1_phi->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);
    hResidua_True_Refit_photon1_phi->GetXaxis()->SetRangeUser(-200,400);
    hResidua_True_Refit_photon1_phi->RebinX(100);
    //hResidua_True_Refit_photon1_phi->RebinX(10);
    hResidua_True_Refit_photon1_phi->Draw();
	ab->Update();

    ab->cd(3);
    gStyle->SetOptTitle(0);
    //gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.2);
    hResidua_Refit_Reco_photon1_phi->RebinX(10);
	hResidua_Refit_Reco_photon1_phi->SetXTitle("#phi Refit-Reco [#circ]");
	hResidua_Refit_Reco_photon1_phi->SetYTitle("Counts");
    

	
	Histogram_adjust(hResidua_Refit_Reco_photon1_phi,0.11,0.07,0.4,0.7,-0.01,-100);
   // hResidua_Refit_Reco_photon1_phi->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    hResidua_Refit_Reco_photon1_phi->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    hResidua_Refit_Reco_photon1_phi->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
	
	//gPad->SetLogy();
    //hResidua_Refit_Reco_photon1_phi->GetYaxis()->ChangeLabel(2,-1,0,-1,-1,-1,-1);
    hResidua_Refit_Reco_photon1_phi->GetXaxis()->SetRangeUser(-400,200);
    hResidua_Refit_Reco_photon1_phi->Draw();
    //TLine line;
   // line.SetLineColor(kGreen);
   // line.SetLineWidth(2);
    //line.DrawLine(0.01,0,0.01,100000);
   // TText *t = new TText();
   // t->SetTextColor(kGreen+2);
    //t->SetTextFont(43);
    
   // t->SetTextSize(25);
   // t->DrawText(0.05,15000,"Probability");
   // t->DrawText(0.05,12000,"cut: P>0.01");
     

	ab->Update();
    ab->SaveAs("phi_res.png");

    
   // ->"hResidua_True_Refit_photon1_phi",
   // hResidua_Refit_Reco_photon1_phi= new TH1F("hResidua_Refit_Reco_photon1_phi","Residua Refit-Reco phi photon1",10000,-100,100);


/*
    Photon_pulls->cd(2);

    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullThetaPhoton1Conv_mass_mix->SetXTitle("Pull #theta");
    hPullThetaPhoton1Conv_mass_mix->SetYTitle("Counts");
    hPullThetaPhoton1Conv_mass_mix->RebinX(5);
    hPullThetaPhoton1Conv_mass_mix->GetXaxis()->SetRangeUser(-0.75,0.75);
	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullThetaPhoton1Conv_mass_mix,0.09,0.09,0.4,0.4,-0.035,-100);

    hPullThetaPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    f_res->SetParameters(300, 0., 1.5, 3);
    
    hPullThetaPhoton1Conv_mass_mix->Draw();
	fit_hist(hPullThetaPhoton1Conv_mass_mix,f_res,-10.0,10.0,tex,0.66,0.8,0.75);
  
	Photon_pulls->Update();
    
    Photon_pulls->cd(3);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullPhiPhoton1Conv_mass_mix->RebinX(5);
    hPullPhiPhoton1Conv_mass_mix->GetXaxis()->SetRangeUser(-0.75,0.75);
    hPullPhiPhoton1Conv_mass_mix->SetXTitle("Pull #phi");
    hPullPhiPhoton1Conv_mass_mix->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullPhiPhoton1Conv_mass_mix,0.09,0.09,0.4,0.4,-0.035,-100);
    hPullPhiPhoton1Conv_mass_mix->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);
    
    hPullPhiPhoton1Conv_mass_mix->Draw();
    f_res->SetParameters(300, 0., 1.5, 3);
    	fit_hist(hPullPhiPhoton1Conv_mass_mix,f_res,-10.0,10.0,tex,0.66,0.8,0.75);
	Photon_pulls->Update();
    */

/*
 
   
    TH1F * hResidua_Refit_Reco_photon1_phi=(TH1F*) f->Get("hResidua_Refit_Reco_photon1_phi");
   

    
    
    TH1F * hResidua_True_Refit_photon1_phi=(TH1F*) f->Get("hResidua_True_Refit_photon1_phi");

   
    
    TH1F * hResidua_True_Reco_photon1_phi
/*

//PHOTON PULLS
TCanvas* c = new TCanvas("c", "c", 200,100,800, 800); //utworz kanwe
    c->Divide(3,1);	
    /////

//TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");
  //  TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");
	c->cd(1);
    auto legend3 = new TLegend(0.6,0.7,0.9,0.9);

    legend3->SetTextSize(0.05);

    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hPullEInvPhoton1Conv_mass->SetXTitle("");
	hPullEInvPhoton1Conv_mass->SetYTitle("Counts");
    hPullEInvPhoton2Conv_mass->SetXTitle("Energy [MeV]");
    hPullEInvPhoton2Conv_mass->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(hPullEInvPhoton1Conv_mass,0.07,0.07,1,1.5,-100,-100);
    //Histogram_adjust(hPullEInvPhoton2Conv_mass,0.07,0.07,1,1.5,-100,-100);
	
	//gPad->SetLogy();
    hPullEInvPhoton1Conv_mass->Draw();
    hPullEInvPhoton1Conv_mass->SetLineColor(kRed);
    hPullEInvPhoton1Conv_mass->Draw("same");
    legend3->AddEntry(hEPhtoton1beforeRefit,"Photon 1","L");
    legend3->AddEntry(hEPhtoton2beforeRefit,"Photon 2","L");
    legend3->Draw();
	c->Update();
   

    c->cd(2);
	gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.5);
    
    auto legend4 = new TLegend(0.6,0.7,0.9,0.9);

    legend4->SetTextSize(0.05);

    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullThetaPhoton1Conv_mass->SetLineColor(kRed);
    legend4->AddEntry(hPullThetaPhoton1Conv_mass,"Photon 1","L");
    legend4->AddEntry(hPullThetaPhoton2Conv_mass,"Photon 2","L");
    Histogram_adjust(hPullThetaPhoton1Conv_mass,0.07,0.07,1,1,-100,-100);
    Histogram_adjust(hPullThetaPhoton2Conv_mass,0.07,0.07,1,1,-100,-100);

    hPullThetaPhoton1Conv_mass->Draw();
    hPullThetaPhoton2Conv_mass->Draw("same");
    legend4->Draw();
	c->Update();
    

     c->cd(3);
	gStyle->SetTitleH(0.1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.5);
    
    auto legend5 = new TLegend(0.6,0.7,0.9,0.9);

    legend5->SetTextSize(0.05);

    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPullPhiPhoton2Conv_mass->SetLineColor(kRed);
    legend5->AddEntry(hPullPhiPhoton1Conv_mass,"Photon 1","L");
    legend5->AddEntry(hPullPhiPhoton2Conv_mass,"Photon 2","L");
    Histogram_adjust(hPullPhiPhoton1Conv_mass,0.07,0.07,1,1,-100,-100);
    Histogram_adjust(hPullPhiPhoton2Conv_mass,0.07,0.07,1,1,-100,-100);

    hPullPhiPhoton1Conv_mass->Draw();
    hPullPhiPhoton2Conv_mass->Draw("same");
    legend5->Draw();
	c->Update();
    


c->SaveAs("pull.png");
*/

/*


TCanvas* d = new TCanvas("d", "d", 200,100,800, 800); //utworz kanwe
    d->Divide(1,1);	
//TH1F * hEPhtoton1beforeRefit= (TH1F*)f->Get("Before Energy Phtoton1");
  //  TH1F * hEPhtoton2beforeRefit= (TH1F*)f->Get("Before Energy Photon2");
	d->cd(1);
    auto legend6 = new TLegend(0.6,0.7,0.9,0.9);

    legend6->SetTextSize(0.05);

    gStyle->SetTitleH(0.1);
	gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
	hEPi0BeforeRefit->SetXTitle("Pi0 mass [MeV]");
	hEPi0BeforeRefit->SetYTitle("Counts");

	//line.DrawLine(0,400,8000,400);
	Histogram_adjust(    hEPi0BeforeRefit,0.07,0.07,1,1.5,-100,-100);
	
	//gPad->SetLogy();
        hEPi0BeforeRefit->Draw();
    hEPi0BeforeRefit_selection->SetLineColor(kRed);
    hEPi0BeforeRefit_selection->Draw("same");

	
    hEPi0AfterRefit->SetLineColor(kGreen);

  
    d->Update();
    d->SaveAs("pi0mass.png");


*/
///////////////
/*
    TCanvas* e = new TCanvas("e", "e", 200,100,800, 800); //utworz kanwe
    e->Divide(3,2);	
    e->cd(1);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hEPi0AfterRefit->SetXTitle("Energy [MeV]");
    hEPi0AfterRefit->SetYTitle("Counts");
    hEPi0AfterRefit->GetXaxis()->SetRangeUser(100,3000);
	Histogram_adjust(hEPi0AfterRefit,0.07,0.04,1,1,-100,-100,10,-100);
    hEPi0AfterRefit->GetXaxis()->SetNdivisions(510,kTRUE);   

    hEPi0AfterRefit->Draw();
    e->Update();
   
    e->cd(2);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hPPi0AfterRefit->SetXTitle("Momentum [MeV/c]");
    hPPi0AfterRefit->SetYTitle("Counts");
    hPPi0AfterRefit->GetXaxis()->SetRangeUser(100,3000);
	Histogram_adjust(hPPi0AfterRefit,0.07,0.04,1,1,-100,-100);
    hPPi0AfterRefit->GetXaxis()->SetNdivisions(510,kTRUE);   
    hPPi0AfterRefit->Draw();
    e->Update();
    e->cd(3);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hmassPi0AfterRefit->SetXTitle("Pi0 mass [MeV/c^2]");
    hmassPi0AfterRefit->SetYTitle("Counts");
    hmassPi0AfterRefit->GetXaxis()->SetRangeUser(134.9,135.1);
	Histogram_adjust(hmassPi0AfterRefit,0.07,0.04,1,1,-100,-100);
   // hPPi0AfterRefit->GetXaxis()->SetNdivisions(510,kTRUE);   
  
    hmassPi0AfterRefit->Draw();
    e->Update();
    e->cd(4);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hthetaPi0AfterRefit->SetXTitle("Pi0 #theta [#circ]");
    hthetaPi0AfterRefit->SetYTitle("Counts");
    hthetaPi0AfterRefit->GetXaxis()->SetRangeUser(0,50);
	Histogram_adjust(hthetaPi0AfterRefit,0.07,0.04,1,1,-100,-100);
    hthetaPi0AfterRefit->Draw();
    e->Update();
    e->cd(5);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hphiPi0AfterRefit->SetXTitle("Pi0 #phi [#circ]");
    hphiPi0AfterRefit->SetYTitle("Counts");
	Histogram_adjust(hphiPi0AfterRefit,0.07,0.04,1,1,-100,-100);
    hphiPi0AfterRefit->Draw();
    e->Update();
    e->Draw();
    
    e->cd(6);
    gPad->SetBottomMargin(0.15);
	gPad->SetLeftMargin(0.25);
    hEPi0BeforeRefit->SetXTitle("Pi0 mass [MeV/c^2]");
    hEPi0BeforeRefit->SetYTitle("Counts");
    hEPi0BeforeRefit->GetXaxis()->SetRangeUser(100,200);
    hEPi0BeforeRefit->RebinX(50000);
    hEPi0BeforeRefit->SetLineColor(kRed);
	Histogram_adjust(hmassPi0AfterRefit,0.07,0.04,1,1,-100,-100);
    Histogram_adjust(hEPi0BeforeRefit,0.07,0.04,1,1,-100,-100);
    auto legendmass = new TLegend(0.6,0.7,0.9,0.9);

    legendmass->SetTextSize(0.05);
    legendmass->AddEntry(hEPi0BeforeRefit,"without Refit","L");
    legendmass->AddEntry(hmassPi0AfterRefit,"Refit","L");
    
    hEPi0BeforeRefit->Draw();
    hmassPi0AfterRefit->Draw("Same");
    legendmass->Draw();
     e->SaveAs("pi0mass_contr.png");

    /*
    TH1D *hpull_P_proton_prim=(TH1D*)f->Get("hpull_P_proton_prim");
    TH1D *hpull_Theta_proton_prim=(TH1D*)f->Get("hpull_Thetha_proton_prim");
    TH1D *hpull_Phi_proton_prim=(TH1D*)f->Get("hpull_Phi_proton_prim");
    TH1D *hpull_R_proton_prim =(TH1D*)f->Get("hpull_R_proton_prim");
    TH1D *hpull_Z_proton_prim=(TH1D*)f->Get("hpull_Z_proton_prim");

    TH1D *hpull_P_lambda_pim=(TH1D*)f->Get("hpull_P_lambda_pim");
    TH1D *hpull_Theta_lambda_pim=(TH1D*)f->Get("hpull_Thetha_lambda_pim");
    TH1D *hpull_Phi_lambda_pim=(TH1D*)f->Get("hpull_Phi_lambda_pim");
    TH1D *hpull_R_lambda_pim=(TH1D*)f->Get("hpull_R_lambda_pim");
    TH1D *hpull_Z_lambda_pim=(TH1D*)f->Get("hpull_Z_lambda_pim");

    TH1D *hpull_P_kp=(TH1D*)f->Get("hpull_P_kp");
    TH1D *hpull_Theta_kp=(TH1D*)f->Get("hpull_Thetha_kp");
    TH1D *hpull_Phi_kp=(TH1D*)f->Get("hpull_Phi_kp");
    TH1D *hpull_R_kp =(TH1D*)f->Get("hpull_R_kp");
    TH1D *hpull_Z_kp=(TH1D*)f->Get("hpull_Z_kp");



/////
    TH1D *hforward_prot_mom_refit=(TH1D*)f->Get("hforward_prot_mom_refit");
    TH1D *htotal_4mom_refit=(TH1D*)f->Get("htotal_4mom_refit");



    ////RESOLUTIONS
    TH1D *hPRes_prim_prot_before=(TH1D*)f->Get("hPRes_prim_prot_before");
    TH1D *hThetaRes_prim_prot_before=(TH1D*)f->Get("hThetaRes_prim_prot_before");
    TH1D *hPhiRes_prim_prot_before=(TH1D*)f->Get("hPhiRes_prim_prot_before");
    TH1D *hResR_prim_prot_before=(TH1D*)f->Get("hResR_prim_prot_before");
    TH1D *hResZ_prim_prot_before=(TH1D*)f->Get("hResZ_prim_prot_before");

    TH1D *hPRes_lambda_pim_before =(TH1D*)f->Get("hPRes_lambda_pim_before");
    TH1D *hThetaRes_lambda_pim_before=(TH1D*)f->Get("hThetaRes_lambda_pim_before");
    TH1D *hPhiRes_lambda_pim_before =(TH1D*)f->Get("hPhiRes_lambda_pim_before");
    TH1D *hResR_lambda_pim_before=(TH1D*)f->Get("hResR_lambda_pim_before");
    TH1D *hResZ_lambda_pim_before=(TH1D*)f->Get("hResZ_lambda_pim_before");

    TH1D *hPRes_lambda_pim=(TH1D*)f->Get("hPRes_lambda_pim");
    TH1D *hThetaRes_lambda_pim=(TH1D*)f->Get("hThetaRes_lambda_pim");
    TH1D *hPhiRes_lambda_pim=(TH1D*)f->Get("hPhiRes_lambda_pim");
    TH1D *hResR_lambda_pim=(TH1D*)f->Get("hResR_lambda_pim");
    TH1D *hResZ_lambda_pim=(TH1D*)f->Get("hResZ_lambda_pim");

    TH1D *hPRes_kaons_plus_before=(TH1D*)f->Get("hPRes_kaons_plus_before");
    TH1D *hThetaRes_kaons_plus_before=(TH1D*)f->Get("hThetaRes_kaons_plus_before");
    TH1D *hPhiRes_kaons_plus_before=(TH1D*)f->Get("hPhiRes_kaons_plus_before");
    TH1D *hResR_kaons_plus_before=(TH1D*)f->Get("hResR_kaons_plus_before");
    TH1D *hResZ_kaons_plus_before=(TH1D*)f->Get("hResZ_kaons_plus_before");

    TH1D *hPRes_kaons_plus=(TH1D*)f->Get("hPRes_kaons_plus_before");
    TH1D *hThetaRes_kaons_plus=(TH1D*)f->Get("hPRes_kaons_plus_before");
    TH1D *hPhiRes_kaons_plus=(TH1D*)f->Get("hPRes_kaons_plus_before");
    TH1D *hResR_kaons_plus=(TH1D*)f->Get("hPRes_kaons_plus_before");
    TH1D *hResZ_kaons_plus=(TH1D*)f->Get("hPRes_kaons_plus_before");

    hPRes_prim_prot_before->SetStats(kFALSE);
    hThetaRes_prim_prot_before->SetStats(kFALSE);
    hPhiRes_prim_prot_before->SetStats(kFALSE);
    hResR_prim_prot_before->SetStats(kFALSE);
    hResZ_prim_prot_before->SetStats(kFALSE);

    hPRes_lambda_pim_before->SetStats(kFALSE);
    hThetaRes_lambda_pim_before->SetStats(kFALSE);
    hPhiRes_lambda_pim_before->SetStats(kFALSE);
    hResR_lambda_pim_before->SetStats(kFALSE);
    hResZ_lambda_pim_before->SetStats(kFALSE);

    hPRes_kaons_plus_before->SetStats(kFALSE);
    hThetaRes_kaons_plus_before->SetStats(kFALSE);
    hPhiRes_kaons_plus_before->SetStats(kFALSE);
    hResR_kaons_plus_before->SetStats(kFALSE);
    hResZ_kaons_plus_before->SetStats(kFALSE);
    ////

    htotal_4mom_refit->SetXTitle("Reconstructed beam invariant mass [MeV]");
    htotal_4mom_refit->SetYTitle("Counts");
    hforward_prot_mom_refit->SetStats(kFALSE);
    htotal_4mom_refit->SetStats(kFALSE);
    hforward_prot_mom_refit->SetXTitle("Refited forward proton momentum [MeV]");
    hforward_prot_mom_refit->SetYTitle("Counts");
    ////////////////////////////////////////////

    hpull_P_proton_prim->SetStats(kFALSE);
    hpull_Theta_proton_prim->SetStats(kFALSE);
    hpull_Phi_proton_prim->SetStats(kFALSE);
    hpull_R_proton_prim->SetStats(kFALSE);
    hpull_Z_proton_prim->SetStats(kFALSE);

    hpull_P_lambda_pim->SetStats(kFALSE);
    hpull_Theta_lambda_pim->SetStats(kFALSE);
    hpull_Phi_lambda_pim->SetStats(kFALSE);
    hpull_R_lambda_pim->SetStats(kFALSE);
    hpull_Z_lambda_pim->SetStats(kFALSE);

    hpull_P_kp->SetStats(kFALSE);
    hpull_Theta_kp->SetStats(kFALSE);
    hpull_Phi_kp->SetStats(kFALSE);
    hpull_R_kp ->SetStats(kFALSE);
    hpull_Z_kp->SetStats(kFALSE);
    
    TH1D *hPRes_prim_prot=(TH1D*)f->Get("hPRes_prim_prot");
    TH1D *hThetaRes_prim_prot =(TH1D*)f->Get("hThetaRes_prim_prot");
    TH1D *hPhiRes_prim_prot=(TH1D*)f->Get("hPhiRes_prim_prot");
    TH1D *hResR_prim_prot=(TH1D*)f->Get("hResR_prim_prot");
    TH1D *hResZ_prim_prot=(TH1D*)f->Get("hResZ_prim_prot");

    TH1D *hPRes_fw_prot=(TH1D*)f->Get("hPRes_fw_prot");
    TH1D *hThetaRes_fw_prot=(TH1D*)f->Get("hThetaRes_fw_prot");
    TH1D *hPhiRes_fw_prot=(TH1D*)f->Get("hPhiRes_fw_prot");
    TH1D *hResR_fw_prot=(TH1D*)f->Get("hResR_fw_prot");
    TH1D *hResZ_fw_prot=(TH1D*)f->Get("hResZ_fw_prot");
   
  //Track parameters
    TH1D *hMom_lambda_pim =(TH1D*)f->Get("hMom_lambda_pim");
    TH1D *hTheta_lambda_pim =(TH1D*)f->Get("hTheta_lambda_pim");
    TH1D *hPhi_lambda_pim =(TH1D*)f->Get("hPhi_lambda_pim");
    TH1D *hR_lambda_pim =(TH1D*)f->Get("hR_lambda_pim");
    TH1D *hZ_lambda_pim=(TH1D*)f->Get("hZ_lambda_pim");
    
    TH1D *hTheta_lambda_pim_before =(TH1D*)f->Get("hTheta_lambda_pim_before");
    TH1D *hPhi_lambda_pim_before =(TH1D*)f->Get("hPhi_lambda_pim_before");
    TH1D *hMom_lambda_pim_before =(TH1D*)f->Get("hMom_lambda_pim_before");
    TH1D *hR_lambda_pim_before =(TH1D*)f->Get("hR_lambda_pim_before");
    TH1D *hZ_lambda_pim_before =(TH1D*)f->Get("hZ_lambda_pim_before");

    TH1D *hMom_fw_prot =(TH1D*)f->Get("hMom_fw_prot");
    TH1D *hTheta_fw_prot =(TH1D*)f->Get("hTheta_fw_prot");
    TH1D *hPhi_fw_prot =(TH1D*)f->Get("hPhi_fw_prot");
    TH1D *hR_fw_prot =(TH1D*)f->Get("hR_fw_prot");
    TH1D *hZ_fw_prot=(TH1D*)f->Get("hZ_fw_prot");
    
    TH1D *hMom_fw_prot_before =(TH1D*)f->Get("hMom_fw_prot_before");
    TH1D *hTheta_fw_prot_before=(TH1D*)f->Get("hTheta_fw_prot_before");
    TH1D *hPhi_fw_prot_before=(TH1D*)f->Get("hPhi_fw_prot_before");
    TH1D *hR_fw_prot_before =(TH1D*)f->Get("hR_fw_prot_before");
    TH1D *hZ_fw_prot_before =(TH1D*)f->Get("hZ_fw_prot_before");

    TH1D *hMom_kp =(TH1D*)f->Get("hMom_kp");
    TH1D *hTheta_kp=(TH1D*)f->Get("hTheta_kp");
    TH1D *hPhi_kp=(TH1D*)f->Get("hPhi_kp");
    TH1D *hR_kp =(TH1D*)f->Get("hR_kp");
    TH1D *hZ_kp =(TH1D*)f->Get("hZ_kp");

    TH1D *hPhi_kp_before =(TH1D*)f->Get("hPhi_kp_before");

    TH1D *hPhi_prim_prot =(TH1D*)f->Get("hPhi_prim_prot");
    TH1D *hPhi_prim_prot_before =(TH1D*)f->Get("hPhi_prim_prot_before");

    hPhi_kp_before->SetStats(kFALSE);
    hPhi_prim_prot->SetStats(kFALSE);
    hPhi_prim_prot_before->SetStats(kFALSE);
   
    ///////////////////
   
    hPRes_fw_prot->SetStats(kFALSE);
    hThetaRes_fw_prot->SetStats(kFALSE);
    hPhiRes_fw_prot->SetStats(kFALSE);
    hResR_fw_prot->SetStats(kFALSE);
    hResZ_fw_prot->SetStats(kFALSE);
   
  //parameters
    hMom_lambda_pim->SetStats(kFALSE);
    hTheta_lambda_pim->SetStats(kFALSE);
    hPhi_lambda_pim->SetStats(kFALSE);
    hR_lambda_pim->SetStats(kFALSE);
    hZ_lambda_pim->SetStats(kFALSE);
    
    hTheta_lambda_pim_before->SetStats(kFALSE);
    hPhi_lambda_pim_before->SetStats(kFALSE);
    hMom_lambda_pim_before->SetStats(kFALSE);
    hR_lambda_pim_before->SetStats(kFALSE);
    hZ_lambda_pim_before->SetStats(kFALSE);

    hMom_fw_prot->SetStats(kFALSE);
    hTheta_fw_prot->SetStats(kFALSE);
    hPhi_fw_prot->SetStats(kFALSE);
    hR_fw_prot->SetStats(kFALSE);
    hZ_fw_prot->SetStats(kFALSE);
    
    hMom_fw_prot_before->SetStats(kFALSE);
    hTheta_fw_prot_before->SetStats(kFALSE);
    hPhi_fw_prot_before->SetStats(kFALSE);
    hR_fw_prot_before->SetStats(kFALSE);
    hZ_fw_prot_before->SetStats(kFALSE);

    hMom_kp->SetStats(kFALSE);
    hTheta_kp->SetStats(kFALSE);
    hPhi_kp->SetStats(kFALSE);
    hR_kp->SetStats(kFALSE);
    hZ_kp->SetStats(kFALSE);

    //

//prob
    Iterations->SetStats(kFALSE);
    hXi2_Refit->SetStats(kFALSE);
    hprob_Refit->SetStats(kFALSE);
    Iterations->SetXTitle("Iterations");
    Iterations->SetYTitle("Counts");
    hXi2_Refit->SetXTitle("#chi^2");
    hprob_Refit->SetXTitle("Probability");
    
    hXi2_Refit->Rebin(10);
    hprob_Refit->Rebin(10);
    TCanvas *can=new TCanvas("can","can", 3000, 3000);
    can->Divide(3,1,0,0);
   
    can->cd(1);
     gPad->SetLogy();
    Histogram_adjust(Iterations, 0.08, 0.06,0.9,0.8,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    Iterations->Draw();
    can->cd(2);
    gPad->SetLogy();
    Histogram_adjust(hXi2_Refit, 0.08, 0.06,0.9,0.8,-100,-100,7,10);
     gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.2);
    hXi2_Refit->Draw();
    
    can->cd(3);
    gPad->SetLogy();
    Histogram_adjust(hprob_Refit, 0.08, 0.06,0.9,0.8,-100,-100,7,10);
     gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.1);
    gPad->SetRightMargin(0.05);
    hprob_Refit->Draw();
    can->Draw();
    can->SaveAs("Probability.png");
    
////////prob end

    hpull_P_proton_prim->SetYTitle("Counts");
    hpull_P_proton_prim->SetXTitle("HADES prot 1/P pull [1/MeV]");
    hpull_Theta_proton_prim->SetXTitle("HADES prot #theta pull [#circ]");
    hpull_Phi_proton_prim->SetXTitle("HADES prot #phi pull [#circ]");
    //
    hpull_P_lambda_pim->SetYTitle("Counts");
    hpull_P_lambda_pim->SetXTitle("HADES #pi- 1/P pull [1/MeV]");
    hpull_Theta_lambda_pim->SetXTitle("HADES #pi- #theta pull [#cir]");
    hpull_Phi_lambda_pim->SetXTitle("HADES #pi- #phi pull [#cir]");
    //
     hpull_P_kp->SetYTitle("Counts");
    hpull_P_kp->SetXTitle("HADES prot 1/P pull [1/MeV]");
    hpull_Theta_kp->SetXTitle("HADES prot #theta pull [#cir] ");
    hpull_Phi_kp->SetXTitle("HADES prot #phi pull [#cir]");
    //pulls
    //////////////////////////////
    hpull_P_proton_prim->Rebin(10);
    hpull_Theta_proton_prim->Rebin(10);
    hpull_Phi_proton_prim->Rebin(10);
    hpull_P_lambda_pim->Rebin(10);
    hpull_Theta_lambda_pim->Rebin(10);
    hpull_Phi_lambda_pim->Rebin(10);
    hpull_P_kp->Rebin(10);
    hpull_Theta_kp->Rebin(10);
    hpull_Phi_kp->Rebin(10);
    TCanvas *PULLS=new TCanvas("PULLS","PULLS", 3000, 3000);
    PULLS->Divide(3,3,0,0);
    PULLS->cd(1);
    gPad->SetLogy();
    Histogram_adjust(hpull_P_proton_prim, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hpull_P_proton_prim->Draw();

    PULLS->cd(2);
    gPad->SetLogy();
    Histogram_adjust(hpull_Theta_proton_prim, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.08);;
    gPad->SetRightMargin(0.05);
    hpull_Theta_proton_prim->Draw();

    PULLS->cd(3);
    gPad->SetLogy();
    Histogram_adjust(hpull_Phi_proton_prim, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.08);;
    gPad->SetRightMargin(0.05);
    hpull_Phi_proton_prim->Draw();

    //////////////////////////////
    PULLS->cd(4);
    gPad->SetLogy();
    Histogram_adjust(hpull_P_lambda_pim , 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hpull_P_lambda_pim->Draw();

    PULLS->cd(5);
    gPad->SetLogy();
    Histogram_adjust(hpull_Theta_lambda_pim, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.08);;
    gPad->SetRightMargin(0.05);
    hpull_Theta_lambda_pim->Draw();

    PULLS->cd(6);
    gPad->SetLogy();
    Histogram_adjust(hpull_Phi_lambda_pim, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.08);;
    gPad->SetRightMargin(0.05);
    hpull_Phi_lambda_pim->Draw();

    ///////////////////////////////
    PULLS->cd(7);
    gPad->SetLogy();
    Histogram_adjust(hpull_P_kp, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hpull_P_kp->Draw();

    PULLS->cd(8);
    gPad->SetLogy();
    Histogram_adjust(hpull_Theta_kp, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.08);;
    gPad->SetRightMargin(0.05);
    hpull_Theta_kp->Draw();

    PULLS->cd(9);
    gPad->SetLogy();
    Histogram_adjust(hpull_Phi_kp, 0.1, 0.1,0.9,0.55,-100,-100,10,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.08);;
    gPad->SetRightMargin(0.05);
    hpull_Phi_kp ->Draw();
    PULLS->SaveAs("PULLS.png");
    TCanvas *FWMom=new TCanvas("FWMom","FWMom", 3000, 3000);
    FWMom->cd();
     Histogram_adjust(hforward_prot_mom_refit, 0.07, 0.05,0.9,0.55,-100,-100,10,10);
    hforward_prot_mom_refit->Draw();
    gPad->SetBottomMargin(0.15);
   
    FWMom->SaveAs("hforward_prot_mom_refit.png");

    TCanvas *total_4mom_refit=new TCanvas("total_4mom_refit","total_4mom_refit", 3000, 3000);
    total_4mom_refit->cd();
    gPad->SetLogy();
    gPad->SetBottomMargin(0.15);
    htotal_4mom_refit->GetXaxis()->SetRangeUser(3440,3600);
    Histogram_adjust(htotal_4mom_refit, 0.07, 0.05,0.9,0.55,-100,-100,10,10);


    htotal_4mom_refit->Draw();
    total_4mom_refit->SaveAs("Total4mom.png");

/////////////////////////

/////////////////////////
    TCanvas *RES=new TCanvas("RES","RES", 3000, 3000);
    RES->Divide(3,3,0,0);
    RES->cd(1);
    //gPad->SetLogy();
    Histogram_adjust(hPRes_prim_prot_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPRes_prim_prot_before->SetXTitle("1/P HADES prot resolution");
    hPRes_prim_prot_before->Draw("Same  PMC");
    hPRes_prim_prot->SetLineColor(kRed);
    hPRes_prim_prot->Draw("Same  PMC");
    
    auto legend = new TLegend(0.15,0.6,0.48,0.9);
    legend->SetTextSize(0.15);
    legend->AddEntry(hPRes_prim_prot_before,"Before","L");
    legend->AddEntry(hPRes_prim_prot,"After","L");

   legend->Draw();
   //
RES->cd(2);
    //gPad->SetLogy();
    Histogram_adjust(hThetaRes_prim_prot_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hThetaRes_prim_prot_before->SetXTitle("#theta HADES prot resolution");
    hThetaRes_prim_prot_before->Draw("Same  PMC");
    hThetaRes_prim_prot->SetLineColor(kRed);
    hThetaRes_prim_prot->Draw("Same  PMC");
    
    auto legend2 = new TLegend(0.15,0.6,0.48,0.9);
    legend2->SetTextSize(0.15);
    legend2->AddEntry(hThetaRes_prim_prot_before,"Before","L");
    legend2->AddEntry(hThetaRes_prim_prot,"After","L");

   legend2->Draw();

   //
RES->cd(3);
    //gPad->SetLogy();
    Histogram_adjust(hPhiRes_prim_prot_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhiRes_prim_prot_before->SetXTitle("#phi HADES prot resolution");
    hPhiRes_prim_prot_before->Draw("Same  PMC");
    hPhiRes_prim_prot->SetLineColor(kRed);
    hPhiRes_prim_prot->Draw("Same  PMC");
    
    auto legend3 = new TLegend(0.15,0.6,0.48,0.9);
    legend3->SetTextSize(0.15);
    legend3->AddEntry(hPhiRes_prim_prot_before,"Before","L");
    legend3->AddEntry(hPhiRes_prim_prot,"After","L");

   legend->Draw();
//
RES->cd(4);
    //gPad->SetLogy();
    Histogram_adjust(hPRes_lambda_pim_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPRes_lambda_pim_before->SetXTitle("1/P #pi- resolution");
    hPRes_lambda_pim_before->Draw("Same  PMC");
    hPRes_lambda_pim->SetLineColor(kRed);
    hPRes_lambda_pim->Draw("Same  PMC");
    
    auto legend4 = new TLegend(0.15,0.6,0.48,0.9);
    legend4->SetTextSize(0.15);
    legend4->AddEntry(hPRes_lambda_pim_before,"Before","L");
    legend4->AddEntry(hPRes_lambda_pim,"After","L");

   legend4->Draw();
//

RES->cd(5);
    //gPad->SetLogy();
    Histogram_adjust(hThetaRes_lambda_pim_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hThetaRes_lambda_pim_before->SetXTitle("#theta #pi- resolution");
    hThetaRes_lambda_pim_before->Draw("Same  PMC");
    hThetaRes_lambda_pim->SetLineColor(kRed);
    hThetaRes_lambda_pim->Draw("Same  PMC");
    
    auto legend5 = new TLegend(0.15,0.6,0.48,0.9);
    legend5->SetTextSize(0.15);
    legend5->AddEntry(hThetaRes_lambda_pim_before,"Before","L");
    legend5->AddEntry(hThetaRes_lambda_pim,"After","L");

   legend5->Draw();

//
RES->cd(6);
    //gPad->SetLogy();
    Histogram_adjust(hPhiRes_lambda_pim_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhiRes_lambda_pim_before->SetXTitle("#phi #pi- resolution");
    hPhiRes_lambda_pim_before->Draw("Same  PMC");
    hPhiRes_lambda_pim->SetLineColor(kRed);
    hPhiRes_lambda_pim->Draw("Same  PMC");
    
    auto legend6 = new TLegend(0.15,0.6,0.48,0.9);
    legend6->SetTextSize(0.15);
    legend6->AddEntry(hPhiRes_lambda_pim_before,"Before","L");
    legend6->AddEntry(hPhiRes_lambda_pim,"After","L");

   legend6->Draw();
   //
   RES->cd(7);
    //gPad->SetLogy();
    Histogram_adjust(hPRes_kaons_plus_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPRes_kaons_plus_before->SetXTitle("1/P K+ resolution");
    hPRes_kaons_plus_before->Draw("Same  PMC");
    hPRes_kaons_plus->SetLineColor(kRed);
    hPRes_kaons_plus->Draw("Same  PMC");
    
    auto legend7 = new TLegend(0.15,0.6,0.48,0.9);
    legend7->SetTextSize(0.15);
    legend7->AddEntry(hPRes_kaons_plus_before,"Before","L");
    legend7->AddEntry(hPRes_kaons_plus,"After","L");

   legend7->Draw();
   //
    RES->cd(8);
    //gPad->SetLogy();
    Histogram_adjust(hThetaRes_kaons_plus_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hThetaRes_kaons_plus_before->SetXTitle("#theta K+ resolution");
    hThetaRes_kaons_plus_before->Draw("Same  PMC");
    hThetaRes_kaons_plus->SetLineColor(kRed);
    hThetaRes_kaons_plus->Draw("Same  PMC");
    
    auto legend8 = new TLegend(0.15,0.6,0.48,0.9);
    legend8->SetTextSize(0.15);
    legend8->AddEntry(hThetaRes_kaons_plus_before,"Before","L");
    legend8->AddEntry(hThetaRes_kaons_plus,"After","L");

   legend8->Draw();
   //
   RES->cd(9);

    //gPad->SetLogy();
    Histogram_adjust(hPhiRes_kaons_plus_before, 0.1, 0.1,0.9,0.55,-100,-100,5,10);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhiRes_kaons_plus_before->SetXTitle("#phi K+ resolution");
    hPhiRes_kaons_plus_before->Draw("Same  PMC");
    hPhiRes_kaons_plus->SetLineColor(kRed);
    hPhiRes_kaons_plus->Draw("Same  PMC");
    
    auto legend9 = new TLegend(0.15,0.6,0.48,0.9);
    legend9->SetTextSize(0.15);
    legend9->AddEntry(hPhiRes_kaons_plus_before,"Before","L");
    legend9->AddEntry(hPhiRes_kaons_plus,"After","L");

   legend9->Draw();
   RES->SaveAs("RES.png");
   //
    TCanvas *PHI=new TCanvas("PHI","PHI", 3000, 3000);
    PHI->Divide(3,2,0,0);
    PHI->cd(1);
    hPhi_prim_prot_before->SetXTitle("HADES prot #phi before refit");
    hPhi_prim_prot_before->SetYTitle("Counts");
    Histogram_adjust(hPhi_prim_prot_before, 0.1, 0.07,0.9,0.55,-100,-100,5,5);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhi_prim_prot_before->Draw();
     PHI->cd(2);
    hPhi_lambda_pim_before->SetXTitle("#pi- #phi before refit");
    hPhi_lambda_pim_before->SetYTitle("Counts");
    Histogram_adjust(hPhi_lambda_pim_before, 0.1, 0.07,0.9,0.55,-100,-100,5,5);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhi_lambda_pim_before->Draw();
    PHI->cd(3);
    hPhi_kp_before->SetXTitle("K+ #phi before refit");
    hPhi_kp_before->SetYTitle("Counts");
    Histogram_adjust(hPhi_kp_before, 0.1, 0.07,0.9,0.55,-100,-100,5,5);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhi_kp_before->Draw();
    
       PHI->cd(4);
    hPhi_prim_prot->SetXTitle("HADES prot #phi refit");
    hPhi_prim_prot->SetYTitle("Counts");
    Histogram_adjust(hPhi_prim_prot, 0.1, 0.07,0.9,0.55,-100,-100,5,5);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhi_prim_prot->Draw();
    PHI->cd(5);
    hPhi_lambda_pim->SetXTitle("#pi- #phi refit");
    hPhi_lambda_pim->SetYTitle("Counts");
    Histogram_adjust(hPhi_lambda_pim, 0.1, 0.07,0.9,0.55,-100,-100,5,5);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhi_lambda_pim->Draw();
    PHI->cd(6);
    hPhi_kp->SetXTitle("K+ #phi refit");
    hPhi_kp->SetYTitle("Counts");
    Histogram_adjust(hPhi_kp, 0.1, 0.07,0.9,0.55,-100,-100,5,5);
    gPad->SetBottomMargin(0.2);
    gPad->SetLeftMargin(0.15);
    gPad->SetRightMargin(0.05);
    hPhi_kp_before->Draw();
    PHI->SaveAs("PHI.png");

/*
    TH1D *hPRes_prim_prot = new TH1D("hPRes_prim_prot", "", 200, -0.001, 0.001);
    TH1D *hThetaRes_prim_prot = new TH1D("hThetaRes_prim_prot", "", 200, -10, 10);
    TH1D *hPhiRes_prim_prot = new TH1D("hPhiRes_prim_prot", "", 200, -10, 10);
    TH1D *hResR_prim_prot= new TH1D("hResR_prim_prot", "", 200, -35, 35);
    TH1D *hResZ_prim_prot = new TH1D("hResZ_prim_prot", "", 200, -100, 100);

H1D *hPRes_fw_prot=(TH1D*)f->Get("hPRes_fw_prot");
    TH1D *hThetaRes_fw_prot=(TH1D*)f->Get("hThetaRes_fw_prot");
    TH1D *hPhiRes_fw_prot=(TH1D*)f->Get("hPhiRes_fw_prot");
    TH1D *hResR_fw_prot=(TH1D*)f->Get("hResR_fw_prot");
    TH1D *hResZ_fw_prot=(TH1D*)f->Get("hResZ_fw_prot");
   
  //Track parameters
    TH1D *hMom_lambda_pim =(TH1D*)f->Get("hMom_lambda_pim");
    TH1D *hTheta_lambda_pim =(TH1D*)f->Get("hTheta_lambda_pim");
    TH1D *hPhi_lambda_pim =(TH1D*)f->Get("hPhi_lambda_pim");
    TH1D *hR_lambda_pim =(TH1D*)f->Get("hR_lambda_pim");
    TH1D *hZ_lambda_pim=(TH1D*)f->Get("hZ_lambda_pim");
    
    TH1D *hTheta_lambda_pim_before =(TH1D*)f->Get("hTheta_lambda_pim_before");
    TH1D *hPhi_lambda_pim_before =(TH1D*)f->Get("hPhi_lambda_pim_before");
    TH1D *hMom_lambda_pim_before =(TH1D*)f->Get("hMom_lambda_pim_before");
    TH1D *hR_lambda_pim_before =(TH1D*)f->Get("hR_lambda_pim_before");
    TH1D *hZ_lambda_pim_before =(TH1D*)f->Get("hZ_lambda_pim_before");

    TH1D *hMom_fw_prot =(TH1D*)f->Get("hMom_fw_prot");
    TH1D *hTheta_fw_prot =(TH1D*)f->Get("hTheta_fw_prot");
    TH1D *hPhi_fw_prot =(TH1D*)f->Get("hPhi_fw_prot");
    TH1D *hR_fw_prot =(TH1D*)f->Get("hR_fw_prot");
    TH1D *hZ_fw_prot=(TH1D*)f->Get("hZ_fw_prot");
    
    TH1D *hMom_fw_prot_before =(TH1D*)f->Get("hMom_fw_prot_before");
    TH1D *hTheta_fw_prot_before=(TH1D*)f->Get("hTheta_fw_prot_before");
    TH1D *hPhi_fw_prot_before=(TH1D*)f->Get("hPhi_fw_prot_before");
    TH1D *hR_fw_prot_before =(TH1D*)f->Get("hR_fw_prot_before");
    TH1D *hZ_fw_prot_before =(TH1D*)f->Get("hZ_fw_prot_before");

    TH1D *hMom_kp =(TH1D*)f->Get("hMom_kp");
    TH1D *hTheta_kp=(TH1D*)f->Get("hTheta_kp");
    TH1D *hPhi_kp=(TH1D*)f->Get("hPhi_kp");
    TH1D *hR_kp =(TH1D*)f->Get("hR_kp");
    TH1D *hZ_kp =(TH1D*)f->Get("hZ_kp");
    /*hPRes_prim_prot_before
    hThetaRes_prim_prot_before
    hPhiRes_prim_prot_before
   

    hPRes_lambda_pim_before
    hThetaRes_lambda_pim_before
    hPhiRes_lambda_pim_before
    

    hPRes_kaons_plus_before
    hThetaRes_kaons_plus_before
    hPhiRes_kaons_plus_before
    
    //FWMom->Divide(,3,0,0);


   /* TCanvas *IT_Chi2_prob=new TCanvas("IT_Chi2_prob","IT_Chi2_prob", 3000, 3000);
    IT_Chi2_prob->SetLeftMargin(0.15);
    IT_Chi2_prob->Divide(3,1,0,0);
    IT_Chi2_prob->cd(1);
    Iterations->Draw();
    IT_Chi2_prob->cd(2);
    hXi2_Refit->Draw();
    IT_Chi2_prob->cd(3);
    gPad->SetLogy();
    hprob_Refit->Draw();
    
    //angles K0S
/*    TH1D* hTheta_K0S_pim_before_Refit = (TH1D*)f->Get("hTheta_K0S_pim_before_Refit");
    TH1D* hPhi_K0S_pim_before_Refit = (TH1D*)f->Get("hPhi_K0S_pim_before_Refit");
    TH1D* hTheta_K0S_pim_Refit = (TH1D*)f->Get("hTheta_K0S_pim_Refit");
    TH1D* hPhi_K0S_pim_Refit = (TH1D*)f->Get("hPhi_K0S_pim_Refit");

    TH1D* hTheta_K0S_pip_before_Refit= (TH1D*)f->Get("hTheta_K0S_pip_before_Refit");
    TH1D* hPhi_K0S_pip_before_Refit = (TH1D*)f->Get("hPhi_K0S_pip_before_Refit");
    TH1D* hTheta_K0S_pip_Refit = (TH1D*)f->Get("hTheta_K0S_pip_Refit");
    TH1D* hPhi_K0S_pip_Refit = (TH1D*)f->Get("hPhi_K0S_pip_Refit");

    /////
    hTheta_K0S_pim_before_Refit->SetStats(kFALSE);
    hPhi_K0S_pim_before_Refit->SetStats(kFALSE);
    hTheta_K0S_pim_Refit->SetStats(kFALSE);
    hPhi_K0S_pim_Refit->SetStats(kFALSE);

    hTheta_K0S_pip_before_Refit->SetStats(kFALSE);
    hPhi_K0S_pip_before_Refit->SetStats(kFALSE);
    hTheta_K0S_pip_Refit->SetStats(kFALSE);
    hPhi_K0S_pip_Refit->SetStats(kFALSE);



    ////


    
    
     TCanvas *K0S_angles=new TCanvas("K0S_angles","K0S_angles", 3000, 3000);
    K0S_angles->SetLeftMargin(0.15);
    K0S_angles->Divide(4,2,0,0);
////
 K0S_angles->SetLeftMargin(0.13);
    K0S_angles->SetRightMargin(0.01);
    K0S_angles->SetTopMargin(0.01);
   // K0S_angles->SetLeftMargin(0.4);
   
   // FW_mom_prot_geand->GetXaxis()->SetRangeUser(0,4000);
    
    //Histogram_adjust(FW_mom_prot_geand, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);
    //FW_mom_prot_geand->SetXTitle("Forwad proton P_{gen} [MeV/c^2]");

    /// 
    
    K0S_angles->cd(1);
     TH1D* hTheta_K0S_pim_before_Refit_wide = new TH1D("hTheta_K0S_pim_before_Refit_wide", "", 90, -12.5, 100);
    hTheta_K0S_pim_before_Refit_wide->SetStats(kFALSE);
    hTheta_K0S_pim_before_Refit_wide->Add(hTheta_K0S_pim_before_Refit);
    //hTheta_K0S_pim_before_Refit_wide->GetXaxis()->SetRangeUser(-10,2000);

    hTheta_K0S_pim_before_Refit_wide->SetXTitle("#theta_{K^{0}_{S}#pi- #left[#circ#right]}");
    hTheta_K0S_pim_before_Refit_wide->GetXaxis()->SetRangeUser(-12,90);
    gPad->SetBottomMargin(0.3);
    hTheta_K0S_pim_before_Refit_wide->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pim_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pim_before_Refit_wide->GetXaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);

    Histogram_adjust(hTheta_K0S_pim_before_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hTheta_K0S_pim_before_Refit_wide->Draw();
    K0S_angles->cd(2);
    TH1D* hPhi_K0S_pim_before_Refit_wide = new TH1D("hPhi_K0S_pim_before_Refit_wide", "", 174, -75, 360);
    
    hPhi_K0S_pim_before_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_K0S_pim_before_Refit->GetBinContent(i);j++)
        {
            hPhi_K0S_pim_before_Refit_wide->Fill(i*2.5-0.1);
        }
    }
    //hPhi_K0S_pim_before_Refit_wide->Add(hPhi_K0S_pim_before_Refit);
    //hPhi_K0S_pim_before_Refit_wide->GetXaxis()->SetLimits(-25,360);


    hPhi_K0S_pim_before_Refit_wide->SetXTitle("#phi_{K^{0}_{S}#pi- #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    Histogram_adjust(hPhi_K0S_pim_before_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hPhi_K0S_pim_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
   // hPhi_K0S_pim_before_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
   // hPhi_K0S_pim_before_Refit->GetXaxis()->SetLimits(0,400);
        gPad->SetBottomMargin(0.3);
    hPhi_K0S_pim_before_Refit_wide->Draw();
    K0S_angles->cd(3);
     TH1D* hTheta_K0S_pip_before_Refit_wide = new TH1D("hTheta_K0S_pip_before_Refit_wide", "", 90, -12.5, 100);
    hTheta_K0S_pip_before_Refit_wide->SetStats(kFALSE);
    hTheta_K0S_pip_before_Refit_wide->Add(hTheta_K0S_pip_before_Refit);
    hTheta_K0S_pip_before_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);

    hTheta_K0S_pip_before_Refit_wide->SetXTitle("#theta_{K^{0}_{S}#pi+ #left[#circ#right]}");
    //hTheta_K0S_pip_before_Refit_wide->GetXaxis()->SetRangeUser(0,90);
    gPad->SetBottomMargin(0.3);
    hTheta_K0S_pip_before_Refit_wide->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pip_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pip_before_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);

    Histogram_adjust(hTheta_K0S_pip_before_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hTheta_K0S_pip_before_Refit_wide->Draw();
    K0S_angles->cd(4);
    TH1D* hPhi_K0S_pip_before_Refit_wide = new TH1D("hPhi_K0S_pip_before_Refit_wide", "", 174, -75, 360);
    
    hPhi_K0S_pip_before_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pip_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_K0S_pip_before_Refit->GetBinContent(i);j++)
        {
            hPhi_K0S_pip_before_Refit_wide->Fill(i*2.5-0.1);
        }
    }
    hPhi_K0S_pip_before_Refit_wide->SetXTitle("#phi_{K^{0}_{S}#pi+ #left[#circ#right]}");
    //hPhi_K0S_pip_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_K0S_pip_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_K0S_pip_before_Refit_wide->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    gPad->SetBottomMargin(0.3);
     gPad->SetRightMargin(0.01);
    Histogram_adjust(hPhi_K0S_pip_before_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hPhi_K0S_pip_before_Refit_wide->Draw();
    K0S_angles->cd(5);
     TH1D* hTheta_K0S_pim_Refit_wide = new TH1D("hTheta_K0S_pim_Refit_wide", "", 90, -12.5, 100);
    hTheta_K0S_pim_Refit_wide->SetStats(kFALSE);
    hTheta_K0S_pim_Refit_wide->Add(hTheta_K0S_pim_Refit);

    hTheta_K0S_pim_Refit_wide->SetXTitle("#theta_{K^{0}_{S}#pi- #left[#circ#right]}");
    
    hTheta_K0S_pim_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);
    gPad->SetBottomMargin(0.3);
    hTheta_K0S_pim_Refit_wide->GetYaxis()->ChangeLabel(10,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pim_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pim_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);

    Histogram_adjust(hTheta_K0S_pim_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hTheta_K0S_pim_Refit_wide->Draw();
    K0S_angles->cd(6);
    TH1D* hPhi_K0S_pim_Refit_wide = new TH1D("hPhi_K0S_pim_Refit_wide", "", 174, -75, 360);
    
    hPhi_K0S_pim_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_K0S_pim_Refit->GetBinContent(i);j++)
        {
            hPhi_K0S_pim_Refit_wide->Fill(i*2.5-0.1);
        }
    }
    //hPhi_K0S_pim_Refit_wide->Add(hPhi_K0S_pim_before_Refit);
    //hPhi_K0S_pim_Refit_wide->GetXaxis()->SetLimits(-25,360);
    hPhi_K0S_pim_Refit_wide->SetXTitle("#phi_{K^{0}_{S}#pi- #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    Histogram_adjust(hPhi_K0S_pim_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hPhi_K0S_pim_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
        hPhi_K0S_pim_Refit_wide->GetYaxis()->ChangeLabel(4,-1,0,-1,-1,-1,-1);

    //hPhi_K0S_pim_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_K0S_pim_Refit_wide->Draw();
    K0S_angles->cd(7);
     TH1D* hTheta_K0S_pip_Refit_wide = new TH1D("hTheta_K0S_pip_Refit_wide", "", 90, -12.5, 100);
    hTheta_K0S_pip_Refit_wide->SetStats(kFALSE);
    hTheta_K0S_pip_Refit_wide->Add(hTheta_K0S_pip_Refit);
    hTheta_K0S_pip_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);
    hTheta_K0S_pip_Refit_wide->SetXTitle("#theta_{K^{0}_{S}#pi+ #left[#circ#right]}");
    //hTheta_K0S_pip_Refit_wide->GetXaxis()->SetRangeUser(0,90);
    gPad->SetBottomMargin(0.3);
    Histogram_adjust(hTheta_K0S_pip_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.17,6,4);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hTheta_K0S_pip_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
   // hTheta_K0S_pip_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pip_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hTheta_K0S_pip_Refit_wide->Draw();
    K0S_angles->cd(8);
    TH1D* hPhi_K0S_pip_Refit_wide = new TH1D("hPhi_K0S_pip_Refit_wide", "", 174, -75, 360);
    
    hPhi_K0S_pip_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pip_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_K0S_pip_Refit->GetBinContent(i);j++)
        {
            hPhi_K0S_pip_Refit_wide->Fill(i*2.5-0.1);
        }
    }
    hPhi_K0S_pip_Refit_wide->SetXTitle("#phi_{K^{0}_{S}#pi+ #left[#circ#right]}");
    //hPhi_K0S_pip_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_K0S_pip_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    gPad->SetBottomMargin(0.3);
     gPad->SetRightMargin(0.01);
    Histogram_adjust(hPhi_K0S_pip_Refit_wide, 0.13, 0.055,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hPhi_K0S_pip_Refit_wide->Draw();
    K0S_angles->SaveAs("K0S_angless.png");


    //langles lambda
  
    TH1D* hTheta_Lambda_pim_before_Refit = (TH1D*)f->Get("hTheta_Lambda_pim_before_Refit");
    TH1D* hPhi_Lambda_pim_before_Refit = (TH1D*)f->Get("hPhi_Lambda_pim_before_Refit");
    TH1D* hTheta_Lambda_pim_Refit = (TH1D*)f->Get("hTheta_Lambda_pim_Refit");
    TH1D* hPhi_Lambda_pim_Refit = (TH1D*)f->Get("hPhi_Lambda_pim_Refit");

    TH1D* hTheta_Lambda_prot_before_Refit= (TH1D*)f->Get("hTheta_Lambda_prot_before_Refit");
    TH1D* hPhi_Lambda_prot_before_Refit = (TH1D*)f->Get("hPhi_Lambda_prot_before_Refit");
    TH1D* hTheta_Lambda_prot_Refit = (TH1D*)f->Get("hTheta_Lambda_prot_Refit");
    TH1D* hPhi_Lambda_prot_Refit = (TH1D*)f->Get("hPhi_Lambda_prot_Refit");

    /////

    hTheta_Lambda_pim_before_Refit->SetStats(kFALSE);
    hPhi_Lambda_pim_before_Refit->SetStats(kFALSE);
    hTheta_Lambda_pim_Refit->SetStats(kFALSE);
    hPhi_Lambda_pim_Refit->SetStats(kFALSE);

    hTheta_Lambda_prot_before_Refit->SetStats(kFALSE);
    hPhi_Lambda_prot_before_Refit->SetStats(kFALSE);
    hTheta_Lambda_prot_Refit->SetStats(kFALSE);
    hPhi_Lambda_prot_Refit->SetStats(kFALSE);

    /////
    
    TCanvas *Lambda_angles=new TCanvas("Lambda_angles","Lambda_angles", 3000, 3000);
    Lambda_angles->SetLeftMargin(0.15);
    Lambda_angles->Divide(4,2,0,0);
    Lambda_angles->SetRightMargin(0.01);
    Lambda_angles->SetTopMargin(0.01);
    
   
    Lambda_angles->cd(1);
    TH1D* hTheta_Lambda_pim_before_Refit_wide = new TH1D("hTheta_Lambda_pim_before_Refit_wide", "", 90, -12.5, 100);
    hTheta_Lambda_pim_before_Refit_wide->SetStats(kFALSE);
    hTheta_Lambda_pim_before_Refit_wide->Add(hTheta_Lambda_pim_before_Refit);
    hTheta_Lambda_pim_before_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);
    gPad->SetBottomMargin(0.3);
    hTheta_Lambda_pim_before_Refit_wide->SetXTitle("#theta_{#Lambda#pi- #left[#circ#right]}");
    hTheta_Lambda_pim_before_Refit_wide->GetYaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
    hTheta_Lambda_pim_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
     hTheta_Lambda_pim_before_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hTheta_Lambda_pim_before_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.17,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hTheta_Lambda_pim_before_Refit_wide->Draw();
    Lambda_angles->cd(2);
    
     TH1D* hPhi_Lambda_pim_before_Refit_wide = new TH1D("hPhi_Lambda_pim_before_Refit_wide", "", 174, -75, 360);
    
    hPhi_Lambda_pim_before_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_Lambda_pim_before_Refit->GetBinContent(i);j++)
        {
            hPhi_Lambda_pim_before_Refit_wide->Fill(i*2.5-0.1);
        }
    }
    //hPhi_Lambda_pim_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_pim_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_pim_before_Refit_wide->GetXaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);

    hPhi_Lambda_pim_before_Refit_wide->SetXTitle("#phi_{#Lambda#pi- #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    Histogram_adjust(hPhi_Lambda_pim_before_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hPhi_Lambda_pim_before_Refit_wide->Draw();
    Lambda_angles->cd(3);
    TH1D* hTheta_Lambda_prot_before_Refit_wide = new TH1D("hTheta_Lambda_prot_before_Refit_wide", "", 90, -12.5, 100);
    hTheta_Lambda_prot_before_Refit_wide->SetStats(kFALSE);
    hTheta_Lambda_prot_before_Refit_wide->Add(hTheta_Lambda_prot_before_Refit);
    hTheta_Lambda_prot_before_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);
    //hTheta_Lambda_prot_before_Refit_wide->GetXaxis()->SetRangeUser(0,100);
    hTheta_Lambda_prot_before_Refit_wide->SetXTitle("#theta_{#Lambda proton #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    //hTheta_Lambda_prot_before_Refit_wide->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hTheta_Lambda_prot_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
   // hTheta_Lambda_prot_before_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
   hTheta_Lambda_prot_before_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hTheta_Lambda_prot_before_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.19,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hTheta_Lambda_prot_before_Refit_wide->Draw();
    Lambda_angles->cd(4);
    TH1D* hPhi_Lambda_prot_before_Refit_wide = new TH1D("hPhi_Lambda_prot_before_Refit_wide", "", 174, -75, 360);
    
    hPhi_Lambda_prot_before_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_Lambda_prot_before_Refit->GetBinContent(i);j++)
        {
            hPhi_Lambda_prot_before_Refit_wide->Fill(i*2.5-0.1);
        }
    }
     gPad->SetRightMargin(0.01);
   // hPhi_Lambda_prot_before_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_prot_before_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_prot_before_Refit_wide->GetYaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);

    hPhi_Lambda_prot_before_Refit_wide->SetXTitle("#phi_{#Lambda proton #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    //hPhi_Lambda_prot_before_Refit_wide->SetXTitle("#phi before refit [#circ]");
    hPhi_Lambda_prot_before_Refit_wide->SetYTitle("Counts");
    Histogram_adjust(hPhi_Lambda_prot_before_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hPhi_Lambda_prot_before_Refit_wide->Draw();
    Lambda_angles->cd(5);
    TH1D* hTheta_Lambda_pim_Refit_wide = new TH1D("hTheta_Lambda_pim_Refit_wide", "", 90, -12.5, 100);
    hTheta_Lambda_pim_Refit_wide->SetStats(kFALSE);
    hTheta_Lambda_pim_Refit_wide->Add(hTheta_Lambda_pim_Refit);
    hTheta_Lambda_pim_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);
    hTheta_Lambda_pim_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    //hTheta_Lambda_pim_Refit_wide->GetXaxis()->SetRangeUser(0,100);
    hTheta_Lambda_pim_Refit_wide->SetXTitle("#theta_{#Lambda #pi- #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    hTheta_Lambda_pim_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hTheta_Lambda_pim_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.17,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hTheta_Lambda_pim_Refit_wide->Draw();
    
    Lambda_angles->cd(6);
       TH1D* hPhi_Lambda_pim_Refit_wide = new TH1D("hPhi_Lambda_pim_Refit_wide", "", 174, -75, 360);
    
    hPhi_Lambda_pim_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_Lambda_pim_Refit->GetBinContent(i);j++)
        {
            hPhi_Lambda_pim_Refit_wide->Fill(i*2.5-0.1);
        }
    }
    hPhi_Lambda_pim_Refit_wide->GetXaxis()->ChangeLabel(5,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_pim_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_pim_Refit_wide->SetXTitle("#phi_{#Lambda #pi- #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    Histogram_adjust(hPhi_Lambda_pim_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hPhi_Lambda_pim_Refit_wide->Draw();
    Lambda_angles->cd(7);
    TH1D* hTheta_Lambda_prot_Refit_wide = new TH1D("hTheta_Lambda_prot_Refit_wide", "", 90, -12.5, 100);
    hTheta_Lambda_prot_Refit_wide->SetStats(kFALSE);
    hTheta_Lambda_prot_Refit_wide->Add(hTheta_Lambda_prot_Refit);
    hTheta_Lambda_prot_Refit_wide->GetXaxis()->SetRangeUser(-12.5,100);
    hTheta_Lambda_prot_Refit_wide->GetXaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    hTheta_Lambda_prot_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_Lambda_prot_Refit_wide->SetXTitle("#theta_{#Lambda proton #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    //hTheta_Lambda_prot_Refit_wide->GetXaxis()->SetRangeUser(0,100);
    Histogram_adjust(hTheta_Lambda_prot_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.19,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
  //  hTheta_Lambda_prot_Refit_wide->GetYaxis()->ChangeLabel(6,-1,0,-1,-1,-1,-1);
    //hTheta_Lambda_prot_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_Lambda_prot_Refit_wide->Draw();
    Lambda_angles->cd(8);
    TH1D* hPhi_Lambda_prot_Refit_wide = new TH1D("hPhi_Lambda_prot_Refit_wide", "", 174, -75, 360);
    
    hPhi_Lambda_prot_Refit_wide->SetStats(kFALSE);
    for(int i=1;i<145;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPhi_Lambda_prot_Refit->GetBinContent(i);j++)
        {
            hPhi_Lambda_prot_Refit_wide->Fill(i*2.5-0.1);
        }
    }
     gPad->SetRightMargin(0.01);
    hPhi_Lambda_prot_Refit_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
   // hPhi_Lambda_prot_Refit_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Lambda_prot_Refit_wide->SetXTitle("#phi_{#Lambda proton #left[#circ#right]}");
    gPad->SetBottomMargin(0.3);
    Histogram_adjust(hPhi_Lambda_prot_Refit_wide, 0.13, 0.065,0.9,0.35,-100,-0.15,6,6);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);

    hPhi_Lambda_prot_Refit_wide->Draw();
    Lambda_angles->SaveAs("Lambda_angless.png");
    



    /// prob iterations
   
   
    TH1D* hLam_Prob = (TH1D*)f->Get("hLam_Prob");
    TH1D* hK0S_Prob = (TH1D*)f->Get("hK0S_Prob");
    TH1D* hPrim_Prob = (TH1D*)f->Get("hPrim_Prob");
    TH1D* hAll_Prob = (TH1D*)f->Get("hAll_Prob");

    TH1D* hLam_Iteration= (TH1D*)f->Get("hLam_Iteration");
    TH1D* hK0S_Iteration = (TH1D*)f->Get("hK0S_Iteration");
    TH1D* hPrim_Iteration = (TH1D*)f->Get("hPrim_Iteration");
    TH1D* hAll_Iteration = (TH1D*)f->Get("hAll_Iteration");
    TCanvas *prob_iteration=new TCanvas("prob_iteration","prob_iteration", 3000, 3000);


     hLam_Prob->SetStats(kFALSE);
    hK0S_Prob->SetStats(kFALSE);
    hPrim_Prob->SetStats(kFALSE);
    hAll_Prob->SetStats(kFALSE);

    hLam_Iteration->SetStats(kFALSE);
    hK0S_Iteration->SetStats(kFALSE);
    hPrim_Iteration->SetStats(kFALSE);
    hAll_Iteration->SetStats(kFALSE);
    //prob_iteration->SetLogy();
    gPad->SetLogy();
    prob_iteration->SetRightMargin(0.01);
    prob_iteration->SetTopMargin(0.01);
    prob_iteration->SetLeftMargin(0.1);
    prob_iteration->Divide(4,2,0,0);
   
    prob_iteration->cd(1);
    TH1D* hLam_Prob_wide = new TH1D("hLam_Prob_wide", "", 1250, -0.25, 1);
    
    hLam_Prob_wide->SetStats(kFALSE);
    for(int i=1;i<1250;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hLam_Prob->GetBinContent(i);j++)
        {
            hLam_Prob_wide->Fill(i*0.001-0.0001);
        }
    }


   //TH1D* hLam_Prob_wide = new TH1D("hLam_Prob_wide", "", 125, -0.25, 1);
    hLam_Prob_wide->SetStats(kFALSE);
    //hLam_Prob_wide->Add(hLam_Prob);
    //hLam_Prob_wide->GetXaxis()->SetRangeUser(-10,2000);

    Histogram_adjust(hLam_Prob_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,-100,-100);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    gPad->SetBottomMargin(0.3);
    hLam_Prob_wide->SetYTitle("Counts");
    hLam_Prob_wide->SetXTitle("P_{#Lambda}");//#Lambda refit probability");
    //prob_iteration->GetSelectedPad()->SetLogy();
    gPad->SetLogy();
    hLam_Prob_wide->Draw();
    gPad->SetLogy();
    hLam_Prob_wide->Rebin(10);
    hLam_Prob_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hLam_Prob_wide->GetXaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
    
    //hLam_Prob->SetYTitle("Counts /0.01");
    prob_iteration->cd(2);
      TH1D* hK0S_Prob_wide = new TH1D("hK0S_Prob_wide", "", 1250, -0.25, 1);
    
    hK0S_Prob_wide->SetStats(kFALSE);
    for(int i=1;i<1250;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hK0S_Prob->GetBinContent(i);j++)
        {
            hK0S_Prob_wide->Fill(i*0.001-0.0001);
        }
    }
    hK0S_Prob_wide->GetXaxis()->ChangeLabel(11,-1,0,-1,-1,-1,-1);
    hK0S_Prob_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hK0S_Prob_wide->GetXaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);

    gPad->SetBottomMargin(0.3);
    gPad->SetLogy();
    hK0S_Prob_wide->SetXTitle("P_{K^{0}_{S}}");//K^{0}_{S} refit probability");
    Histogram_adjust(hK0S_Prob_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,-100,-100);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hK0S_Prob_wide->Rebin(10);
    hK0S_Prob_wide->SetYTitle("Counts /0.01");
    hK0S_Prob_wide->Draw();
    prob_iteration->cd(3);
    TH1D* hPrim_Prob_wide = new TH1D("hPrim_Prob_wide", "", 1250, -0.25, 1);
    
    hPrim_Prob_wide->SetStats(kFALSE);
    for(int i=1;i<1250;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPrim_Prob->GetBinContent(i);j++)
        {
            hPrim_Prob_wide->Fill(i*0.001-0.0001);
        }
    }
   
    hPrim_Prob_wide->GetXaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);
    hPrim_Prob_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    

    hPrim_Prob_wide->SetXTitle("P_{prim}");//Primary refit probability");
    gPad->SetBottomMargin(0.3);
    gPad->SetLogy();
    Histogram_adjust(hPrim_Prob_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,-100,-100);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hPrim_Prob_wide->Rebin(10);
    hPrim_Prob_wide->SetYTitle("Counts /0.01");
    hPrim_Prob_wide->Draw();
    prob_iteration->cd(4);
     TH1D* hAll_Prob_wide = new TH1D("hAll_Prob_wide", "", 1250, -0.25, 1);
    
    hAll_Prob_wide->SetStats(kFALSE);
    for(int i=1;i<1250;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hAll_Prob->GetBinContent(i);j++)
        {
            hAll_Prob_wide->Fill(i*0.001-0.0001);
        }
    }
    hAll_Prob_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hAll_Prob_wide->GetXaxis()->ChangeLabel(7,-1,0,-1,-1,-1,-1);

    hAll_Prob_wide->SetXTitle("P_{Full}");
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.3);
    gPad->SetLogy();
    Histogram_adjust(hAll_Prob_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,-100,-100);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hAll_Prob_wide->Rebin(10);
    hAll_Prob_wide->SetYTitle("Counts /0.01");
    hAll_Prob_wide->Draw();
    hAll_Prob_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    prob_iteration->cd(5);
     TH1D* hLam_Iteration_wide = new TH1D("hLam_Iteration_wide", "", 11, -1, 10);
    
    hLam_Iteration_wide->SetStats(kFALSE);
    for(int i=1;i<11;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hLam_Iteration->GetBinContent(i);j++)
        {
            hLam_Iteration_wide->Fill(i-0.1);
        }
    }
    hLam_Iteration_wide->SetYTitle("Counts");
    hLam_Iteration_wide->SetXTitle("#Lambda refit iterations");
    gPad->SetLogy();
    gPad->SetBottomMargin(0.3);
   // hLam_Iteration_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hLam_Iteration_wide->GetXaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    hLam_Iteration_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hLam_Iteration_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,20,20);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    //hLam_Iteration_wide->GetXaxis()->SetRangeUser(0,10);
    hLam_Iteration_wide->Draw();
    prob_iteration->cd(6);
     TH1D* hK0S_Iteration_wide = new TH1D("hK0S_Iteration_wide", "", 11, -1, 10);
    
    hK0S_Iteration_wide->SetStats(kFALSE);
    for(int i=1;i<11;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hK0S_Iteration->GetBinContent(i);j++)
        {
            hK0S_Iteration_wide->Fill(i-0.1);
        }
    }
    hK0S_Iteration_wide->SetXTitle("K^{0}_{S} refit iterations");
   // hK0S_Iteration_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    //hK0S_Iteration_wide->GetXaxis()->ChangeLabel(11,-1,0,-1,-1,-1,-1);
    gPad->SetLogy();
    gPad->SetBottomMargin(0.3);
    hK0S_Iteration_wide->GetXaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    hK0S_Iteration_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hK0S_Iteration_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,20,20);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
   // hK0S_Iteration_wide->GetXaxis()->SetRangeUser(0,10);
    hK0S_Iteration_wide->Draw();
    prob_iteration->cd(7);
        TH1D* hPrim_Iteration_wide = new TH1D("hPrim_Iteration_wide", "", 11, -1, 10);
    
    hPrim_Iteration_wide->SetStats(kFALSE);
    for(int i=1;i<11;i++)
    {
       // cout<<hPhi_K0S_pim_before_Refit->GetBinContent(i)<<endl;
        for(int j=0;j<hPrim_Iteration->GetBinContent(i);j++)
        {
            hPrim_Iteration_wide->Fill(i-0.1);
        }
    }
    hPrim_Iteration_wide->SetXTitle("Primary refit iterations");
   // hPrim_Iteration_wide->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPrim_Iteration_wide->GetXaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    hPrim_Iteration_wide->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    gPad->SetLogy();
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.3);
  //  hPrim_Iteration_wide->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    Histogram_adjust(hPrim_Iteration_wide, 0.11, 0.08,0.9,0.35,-100,-0.15,20,20);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
//hPrim_Iteration_wide->GetXaxis()->SetRangeUser(0,10);

    hPrim_Iteration_wide->Draw();
    //prob_iteration->cd(8);
    //hAll_Iteration->Draw();
    prob_iteration->SaveAs("prob_iterations.png");

    
    
   
    //////////////////////////mom before and after
  

    TH1D* h4Mom_Lambda_before_Refit = (TH1D*)f->Get("h4Mom_Lambda_before_Refit");
    TH1D* h4Mom_Lambda_after_Refit = (TH1D*)f->Get("h4Mom_Lambda_after_Refit");
    TH1D* h4Mom_K0S_before_Refit = (TH1D*)f->Get("h4Mom_K0S_before_Refit");
    TH1D* h4Mom_K0S_after_Refit = (TH1D*)f->Get("h4Mom_K0S_after_Refit");
    


    h4Mom_Lambda_before_Refit->SetStats(kFALSE);
    h4Mom_Lambda_after_Refit->SetStats(kFALSE);
    h4Mom_K0S_before_Refit->SetStats(kFALSE);
    h4Mom_K0S_after_Refit->SetStats(kFALSE);
    TCanvas *can=new TCanvas("can","can", 3000, 3000);
    can->SetLeftMargin(0.13);
    can->SetRightMargin(0.01);
    can->SetTopMargin(0.01);
    can->SetLeftMargin(0.1);
   
    can->Divide(2,2,0,0);
    can->cd(1);
    gPad->SetBottomMargin(0.3);
	h4Mom_Lambda_before_Refit->Draw();
    
    //f_res->SetParameters(450, 0.000075, 0.0001, 3);
    h4Mom_Lambda_before_Refit->SetYTitle("Counts");
    h4Mom_Lambda_before_Refit->SetXTitle("M_{#Lambda} #left[MeV/c^2#right]");
    Histogram_adjust(h4Mom_Lambda_before_Refit, 0.16, 0.07,0.7,0.11,-100,-0.07,-100,9);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    h4Mom_Lambda_before_Refit->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_Lambda_before_Refit->GetXaxis()->ChangeLabel(11,-1,0,-1,-1,-1,-1);
    h4Mom_Lambda_before_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    //h4Mom_Lambda_before_Refit->GetYaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    h4Mom_Lambda_before_Refit->GetXaxis()->SetRangeUser(1000,1500);
    can->cd(2);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.3);
    h4Mom_K0S_before_Refit->SetXTitle("M_{K^{0}_{S}} #left[MeV/c^2#right]");
    h4Mom_K0S_before_Refit->SetYTitle("");
    h4Mom_K0S_before_Refit->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_before_Refit->GetXaxis()->ChangeLabel(9,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_before_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_before_Refit->GetYaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    Histogram_adjust( h4Mom_K0S_before_Refit, 0.16, 0.07,0.7,0.11,-100,-0.07,10,10);
    //h4Mom_K0S_before_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_before_Refit->GetXaxis()->SetRangeUser(200,1000);
    h4Mom_K0S_before_Refit->Draw();
    
    can->cd(3);
    gPad->SetBottomMargin(0.3);
    h4Mom_Lambda_after_Refit->SetXTitle("M_{#Lambda} #left[MeV/c^2#right]");
    h4Mom_Lambda_after_Refit->SetYTitle("Counts");
    h4Mom_Lambda_after_Refit->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_Lambda_after_Refit->GetXaxis()->ChangeLabel(11,-1,0,-1,-1,-1,-1);
    h4Mom_Lambda_after_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    Histogram_adjust(h4Mom_Lambda_after_Refit, 0.16, 0.07,0.7,0.11,-100,-0.07,10,10);
    //h4Mom_Lambda_after_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_Lambda_after_Refit->Draw();
     h4Mom_Lambda_after_Refit->GetXaxis()->SetRangeUser(1000,1500);
    can->cd(4);
    gPad->SetRightMargin(0.01);
    gPad->SetBottomMargin(0.3);
    h4Mom_K0S_after_Refit->GetXaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_after_Refit->GetXaxis()->ChangeLabel(9,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_after_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_after_Refit->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    
   // h4Mom_K0S_after_Refit->SetXTitle("K0S mass after refit [MeV/c^2]");
    h4Mom_K0S_after_Refit->SetXTitle("M_{K^{0}_{S}} #left[MeV/c^2#right]");
    h4Mom_K0S_after_Refit->SetYTitle("");
    Histogram_adjust( h4Mom_K0S_after_Refit, 0.16, 0.07,0.7,0.11,-100,-0.07,10,10);
   //h4Mom_K0S_after_Refit->Rebin(30);
    h4Mom_K0S_after_Refit->GetXaxis()->SetRangeUser(200,1000);
    //h4Mom_K0S_after_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    h4Mom_K0S_after_Refit->Draw();
    
	can->Draw();
    can->SaveAs("Reconstructed_mass.png");

    
/*    

	TH1D* hMom_main_before_Refit = (TH1D*)f->Get("hMom_main_before_Refit");

   	//TH1D* hPhi_Forward_before_Refit = (TH1D*)f->Get("hPhi_Forward_before_Refit");
	TH1D* hTheta_Main_before_Refit = (TH1D*)f->Get("hTheta_main_before_Refit");
	TH1D* hPhi_Main_before_Refit = (TH1D*)f->Get("hPhi_main_before_Refit");

	TH1D* hXi2_Refit= (TH1D*)f->Get("hXi2_Refit");
	TH1D* hMom_main_Refit = (TH1D*)f->Get("hMom_main_Refit");
	//TH1D* hPhi_Forward_Refit = (TH1D*)f->Get("hPhi_forward_Refit");
	TH1D* hTheta_Main_Refit = (TH1D*)f->Get("hTheta_main_Refit");
	TH1D* hPhi_Main_Refit = (TH1D*)f->Get("hPhi_main_Refit");
	TH1D* hProb_Refit = (TH1D*)f->Get("hPProb_Refit");

	
    hMom_main_before_Refit->SetStats(kFALSE);
   	//hPhi_Forward_before_Refit->SetStats(kFALSE);
	hTheta_Main_before_Refit->SetStats(kFALSE);
    hPhi_Main_before_Refit->SetStats(kFALSE);
	hXi2_Refit->SetStats(kFALSE);
	hMom_main_Refit->SetStats(kFALSE);
	//hPhi_Forward_Refit->SetStats(kFALSE);
	hTheta_Main_Refit->SetStats(kFALSE);
	hPhi_Main_Refit->SetStats(kFALSE);
    hProb_Refit->SetStats(kFALSE);

    
    TF1 * f_res = new TF1("f_res", "gaus(0)", -30, 30);
    f_res->SetParameters(300, 0., 1.5, 3);
    f_res->SetParameters(600, 0, 100, 10);
/*
    TCanvas *help=new TCanvas("help","help", 3000, 3000);
    help->Divide(4,2);
    help->cd(1);
    
	hMom_main_before_Refit->Draw();
    f_res->SetParameters(450, 0.000075, 0.0001, 3);
     
    Histogram_adjust(hMom_main_before_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);// x_title_offset=-100//, y_title_offset=-100,x_label_offset=-100, y_label_offset=-100,x_divisions=-100, y_divisions=-100);
    hMom_main_before_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hMom_main_before_Refit->GetYaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    
    help->cd(2);
    Histogram_adjust( hTheta_Main_before_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    hTheta_Main_before_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_Main_before_Refit->Draw();
    help->cd(3);
    Histogram_adjust( hPhi_Main_before_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    hPhi_Main_before_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Main_before_Refit->Draw();
    help->cd(4);
    Histogram_adjust( hXi2_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    hXi2_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hXi2_Refit->Draw();
    help->cd(5);
    Histogram_adjust(hMom_main_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    hMom_main_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hMom_main_Refit->GetYaxis()->ChangeLabel(12,-1,0,-1,-1,-1,-1);
    hMom_main_Refit->Draw();
    f_res->SetParameters(450, 0.000075, 0.00001, 0);
    help->cd(6);
    Histogram_adjust( hTheta_Main_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    hTheta_Main_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hTheta_Main_Refit->GetYaxis()->ChangeLabel(8,-1,0,-1,-1,-1,-1);
    hTheta_Main_Refit->Draw();
    f_res->SetParameters(250, 0., 1.5, 3);
    help->cd(7);
    Histogram_adjust(  hPhi_Main_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    hPhi_Main_Refit->GetYaxis()->SetTitle("Counts / 0.25mm");
    hPhi_Main_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    hPhi_Main_Refit->Draw();
    f_res->SetParameters(250, 0., 1.5, 3);
    help->cd(8);
    Histogram_adjust( hProb_Refit, 0.045, 0.035,0.9,0.9,-100,-100,-100,-100);
    
    //hProb_Refit->GetYaxis()->ChangeLabel(1,-1,0,-1,-1,-1,-1);
    //hProb_Refit->GetYaxis()->ChangeLabel(13,-1,0,-1,-1,-1,-1);
    hProb_Refit->Draw();

  
  
	help->Draw();
    help->SaveAs("Refit_results_sim_hadhad.png");
    */

}
