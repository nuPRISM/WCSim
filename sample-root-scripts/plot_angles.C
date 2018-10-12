#include <iostream>
#include <fstream>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <stdio.h>     
#include <stdlib.h>    
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <TRandom3.h>
#include <TSystem.h>
#include "TMatrixD.h"
#include "TVectorD.h"
#include <TMatrixDSym.h>


using namespace std;


#define NEVENTS 10000
#define NHITS 200
#define NPMTS 50000


double get_residual_position(double x0, double y0, double z0, double t0, 
			     double x,
			     double y,
			     double z,
			     double t,
			     double x2,
			     double y2,
			     double z2,
			     double t2,
			     double x3,
			     double y2x,
			     double z2x,
			     double t2x,
			     double xy,
			     double xz,
			     double xt
			     );

double get_residual_d(double x0, double y0, double z0, double t0, 
			     double x,
			     double y,
			     double z,
			     double t,
			     double x2,
			     double y2,
			     double z2,
			     double t2
			     );

double smear(double center, double sigma);

void smear_vtx();

void calculate_averages();

void calculate_reco_dir();

bool calculate_reco2_dir();

void get_range(TH1F * h, double *m, double *M);

double n = 1.336;
double c = 30/n; // cm / ns

double average_x=0;
double average_y=0;
double average_z=0;
double average_d=0;
double average_t=0;
double average_x2=0;
double average_y2=0;
double average_z2=0;
double average_t2=0;
double average_x3=0;
double average_y3=0;
double average_z3=0;
double average_t3=0;
double average_y2x=0;
double average_z2x=0;
double average_t2x=0;
double average_x2y=0;
double average_z2y=0;
double average_t2y=0;
double average_x2z=0;
double average_y2z=0;
double average_t2z=0;
double average_x2t=0;
double average_y2t=0;
double average_z2t=0;
double average_xy=0;
double average_xz=0;
double average_xt=0;
double average_yz=0;
double average_yt=0;
double average_zt=0;
double average_dx=0;
double average_dy=0;
double average_dz=0;
double average_dx2=0;
double average_dy2=0;
double average_dz2=0;
double average_dxy=0;
double average_dxz=0;
double average_dyz=0;

double vtx_sigma = 215; // cm

Float_t vtx_x, vtx_y, vtx_z;
Float_t dir_x, dir_y, dir_z;
Float_t kinetic;
Int_t n_hits;
Float_t hit_x[NHITS];
Float_t hit_y[NHITS];
Float_t hit_z[NHITS];
Float_t hit_t[NHITS];

double cos_theta;
double dist_x, dist_y, dist_z, dist;
double reco_dist_x, reco_dist_y, reco_dist_z, reco_dist, reco_dx, reco_dy, reco_dz;
double tof, tcorr;

double theta;
double pi = acos(-1.);
double cos_theta0 = 1./n;
double theta0 = acos(cos_theta0)*180./pi;
double tmin = 950.;
double tmax = 970.;
double vtx_t = (tmin + tmax)/2.;

double _cos_theta_min = 0.64; // reco
double _cos_theta_max = 1;

double _cos_theta_min2 = 0.61; // reco2
double _cos_theta_max2 = 0.97;

double cos_theta_min[5];
double cos_theta_max[5];

double cos_theta_min2[5];
double cos_theta_max2[5];

double residual_x0, residual_y0, residual_z0, residual_t0, residual_d;
double reco_dir_x, reco_dir_y, reco_dir_z;
double reco2_dir_x, reco2_dir_y, reco2_dir_z, reco2_dir_modulus;
double reco_cos_theta;
double reco2_cos_theta;
double reco_vtx_x, reco_vtx_y, reco_vtx_z;

double A[3][3];
double B[3];
double detA, detAA;
double Ai[3][3];
//double I[3][3];
TMatrixDSym AA(3);
TMatrixDSym AAi(3);
TVectorD AB(3);
TVectorD A_reco2_dir(3);

double n_hits_in_range;
double n_hits_in_range2;

TH1F * h_cos_theta_1;
TH1F * h_cos_theta_2;
TH1F * h_cos_theta_3;
TH1F * h_cos_theta_4;
TH1F * h_cos_theta_5;


int main(){

  clog << " cos_theta0 " << cos_theta0 << " theta0 " << theta0 << endl;

  gSystem->Load("libMatrix");

  TFile input("angles.root","READ");
  TFile * output = new TFile("angles_plot.root","RECREATE");

  TTree * angles_tree = (TTree*)input.Get("angles_tree");

  angles_tree->SetBranchAddress("vtx_x",&vtx_x);
  angles_tree->SetBranchAddress("vtx_y",&vtx_y);
  angles_tree->SetBranchAddress("vtx_z",&vtx_z);
  angles_tree->SetBranchAddress("dir_x",&dir_x);
  angles_tree->SetBranchAddress("dir_y",&dir_y);
  angles_tree->SetBranchAddress("dir_z",&dir_z);
  angles_tree->SetBranchAddress("kinetic",&kinetic);
  angles_tree->SetBranchAddress("n_hits",&n_hits);
  angles_tree->SetBranchAddress("hit_x",hit_x);
  angles_tree->SetBranchAddress("hit_y",hit_y);
  angles_tree->SetBranchAddress("hit_z",hit_z);
  angles_tree->SetBranchAddress("hit_t",hit_t);


  TH1F h_t("h_t","t;time [ns]",200,1,-1);
  TH1F h_tcorr("h_tcorr","tcorr;corrected time [ns]",200,1,-1);

  h_cos_theta_1 = new TH1F("h_cos_theta_1","1 MeV;cos#theta",200,-1,1);
  h_cos_theta_1->SetLineColor(kGray);
  h_cos_theta_1->SetLineWidth(2);
  TH1F h_theta_1("h_theta_1","1 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_1("h_theta_from_theta0_1","1 MeV;#theta - #theta_{0}",200,-360.,360.);
  h_cos_theta_2 = new TH1F("h_cos_theta_2","2 MeV;cos#theta",200,-1,1);
  h_cos_theta_2->SetLineColor(kBlack);
  h_cos_theta_2->SetLineWidth(2);
  TH1F h_theta_2("h_theta_2","2 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_2("h_theta_from_theta0_2","2 MeV;#theta - #theta_{0}",200,-360.,360.);
  h_cos_theta_3 = new TH1F("h_cos_theta_3","3 MeV;cos#theta",200,-1,1);
  h_cos_theta_3->SetLineColor(kRed);
  h_cos_theta_3->SetLineWidth(2);
  TH1F h_theta_3("h_theta_3","3 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_3("h_theta_from_theta0_3","3 MeV;#theta - #theta_{0}",200,-360.,360.);
  h_cos_theta_4 = new TH1F("h_cos_theta_4","4 MeV;cos#theta",200,-1,1);
  h_cos_theta_4->SetLineColor(kBlue);
  h_cos_theta_4->SetLineWidth(2);
  TH1F h_theta_4("h_theta_4","4 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_4("h_theta_from_theta0_4","4 MeV;#theta - #theta_{0}",200,-360.,360.);
  h_cos_theta_5 = new TH1F("h_cos_theta_5","5 MeV;cos#theta",200,-1,1);
  h_cos_theta_5->SetLineColor(kGreen);
  h_cos_theta_5->SetLineWidth(2);
  TH1F h_theta_5("h_theta_5","5 MeV;#theta",200,-360.,360.);
  TH1F h_theta_from_theta0_5("h_theta_from_theta0_5","5 MeV;#theta - #theta_{0}",200,-360.,360.);


  TH1F h_hits_fraction_1("h_hits_fraction_1",Form("1 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min,_cos_theta_max),50,0,1);
  h_hits_fraction_1.SetLineColor(kGray);
  h_hits_fraction_1.SetLineWidth(2);
  TH1F h_hits_fraction_2("h_hits_fraction_2",Form("2 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min,_cos_theta_max),50,0,1);
  h_hits_fraction_2.SetLineColor(kBlack);
  h_hits_fraction_2.SetLineWidth(2);
  TH1F h_hits_fraction_3("h_hits_fraction_3",Form("3 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min,_cos_theta_max),50,0,1);
  h_hits_fraction_3.SetLineColor(kRed);
  h_hits_fraction_3.SetLineWidth(2);
  TH1F h_hits_fraction_4("h_hits_fraction_4",Form("4 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min,_cos_theta_max),50,0,1);
  h_hits_fraction_4.SetLineColor(kBlue);
  h_hits_fraction_4.SetLineWidth(2);
  TH1F h_hits_fraction_5("h_hits_fraction_5",Form("5 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min,_cos_theta_max),50,0,1);
  h_hits_fraction_5.SetLineColor(kGreen);
  h_hits_fraction_5.SetLineWidth(2);

  TH1F h_hits_fraction2_1("h_hits_fraction2_1",Form("1 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min2,_cos_theta_max2),50,0,1);
  h_hits_fraction2_1.SetLineColor(kGray);
  h_hits_fraction2_1.SetLineWidth(2);
  TH1F h_hits_fraction2_2("h_hits_fraction2_2",Form("2 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min2,_cos_theta_max2),50,0,1);
  h_hits_fraction2_2.SetLineColor(kBlack);
  h_hits_fraction2_2.SetLineWidth(2);
  TH1F h_hits_fraction2_3("h_hits_fraction2_3",Form("3 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min2,_cos_theta_max2),50,0,1);
  h_hits_fraction2_3.SetLineColor(kRed);
  h_hits_fraction2_3.SetLineWidth(2);
  TH1F h_hits_fraction2_4("h_hits_fraction2_4",Form("4 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min2,_cos_theta_max2),50,0,1);
  h_hits_fraction2_4.SetLineColor(kBlue);
  h_hits_fraction2_4.SetLineWidth(2);
  TH1F h_hits_fraction2_5("h_hits_fraction2_5",Form("5 MeV;fraction of hits with cos#theta #in(%.2f, %.2f)",_cos_theta_min2,_cos_theta_max2),50,0,1);
  h_hits_fraction2_5.SetLineColor(kGreen);
  h_hits_fraction2_5.SetLineWidth(2);

  TH1F h_analytic_residual_x0("h_analytic_residual_x0","h_analytic_residual_x0; residual x0",100,-100,100);
  h_analytic_residual_x0.SetLineColor(kBlack);
  h_analytic_residual_x0.SetLineWidth(2);
  TH1F h_analytic_residual_y0("h_analytic_residual_y0","h_analytic_residual_y0; residual y0",100,1,-1);
  h_analytic_residual_y0.SetLineColor(kBlack);
  h_analytic_residual_y0.SetLineWidth(2);
  TH1F h_analytic_residual_z0("h_analytic_residual_z0","h_analytic_residual_z0; residual z0",100,1,-1);
  h_analytic_residual_z0.SetLineColor(kBlack);
  h_analytic_residual_z0.SetLineWidth(2);
  TH1F h_analytic_residual_t0("h_analytic_residual_t0","h_analytic_residual_t0; residual t0",100,1,-1);
  h_analytic_residual_t0.SetLineColor(kBlack);
  h_analytic_residual_t0.SetLineWidth(2);

  TH1F h_analytic_residual_d("h_analytic_residual_d","h_analytic_residual_d; residual d",100,1,-1);
  h_analytic_residual_d.SetLineColor(kBlack);
  h_analytic_residual_d.SetLineWidth(2);

  TH2F h_reco_dir_x__vs__dir_x("h_reco_dir_x__vs__dir_x","h_reco_dir_x__vs__dir_x; true dir x; reco dir x",100,-1,1,100,-1,1);
  TH2F h_reco_dir_y__vs__dir_y("h_reco_dir_y__vs__dir_y","h_reco_dir_y__vs__dir_y; true dir y; reco dir y",100,-1,1,100,-1,1);
  TH2F h_reco_dir_z__vs__dir_z("h_reco_dir_z__vs__dir_z","h_reco_dir_z__vs__dir_z; true dir z; reco dir z",100,-1,1,100,-1,1);

  TH1F h_residual_dir_x("h_residual_dir_x","h_residual_dir_x; reco dir x - true dir x",100,-2,2);
  TH1F h_residual_dir_y("h_residual_dir_y","h_residual_dir_y; reco dir y - true dir y",100,-2,2);
  TH1F h_residual_dir_z("h_residual_dir_z","h_residual_dir_z; reco dir z - true dir z",100,-2,2);

  TH1F h_reco_dir_x("h_reco_dir_x","h_reco_dir_x; reco dir x",100,-1,1);
  TH1F h_reco_dir_y("h_reco_dir_y","h_reco_dir_y; reco dir y",100,-1,1);
  TH1F h_reco_dir_z("h_reco_dir_z","h_reco_dir_z; reco dir z",100,-1,1);

  TH2F h_reco2_dir_x__vs__dir_x("h_reco2_dir_x__vs__dir_x","h_reco2_dir_x__vs__dir_x; true dir x; reco dir x",100,-1,1,100,-1,1);
  TH2F h_reco2_dir_y__vs__dir_y("h_reco2_dir_y__vs__dir_y","h_reco2_dir_y__vs__dir_y; true dir y; reco dir y",100,-1,1,100,-1,1);
  TH2F h_reco2_dir_z__vs__dir_z("h_reco2_dir_z__vs__dir_z","h_reco2_dir_z__vs__dir_z; true dir z; reco dir z",100,-1,1,100,-1,1);

  TH1F h_residual2_dir_x("h_residual2_dir_x","h_residual2_dir_x; reco dir x - true dir x",100,-2,2);
  TH1F h_residual2_dir_y("h_residual2_dir_y","h_residual2_dir_y; reco dir y - true dir y",100,-2,2);
  TH1F h_residual2_dir_z("h_residual2_dir_z","h_residual2_dir_z; reco dir z - true dir z",100,-2,2);

  TH1F h_reco2_dir_x("h_reco2_dir_x","h_reco2_dir_x; reco dir x",100,-1,1);
  TH1F h_reco2_dir_y("h_reco2_dir_y","h_reco2_dir_y; reco dir y",100,-1,1);
  TH1F h_reco2_dir_z("h_reco2_dir_z","h_reco2_dir_z; reco dir z",100,-1,1);

  TH1F h_reco_cos_theta_1("h_reco_cos_theta_1","1 MeV;reco_cos#theta",200,-1,1);
  h_reco_cos_theta_1.SetLineColor(kGray);
  h_reco_cos_theta_1.SetLineWidth(2);
  TH1F h_reco_cos_theta_2("h_reco_cos_theta_2","2 MeV;reco_cos#theta",200,-1,1);
  h_reco_cos_theta_2.SetLineColor(kBlack);
  h_reco_cos_theta_2.SetLineWidth(2);
  TH1F h_reco_cos_theta_3("h_reco_cos_theta_3","3 MeV;reco_cos#theta",200,-1,1);
  h_reco_cos_theta_3.SetLineColor(kRed);
  h_reco_cos_theta_3.SetLineWidth(2);
  TH1F h_reco_cos_theta_4("h_reco_cos_theta_4","4 MeV;reco_cos#theta",200,-1,1);
  h_reco_cos_theta_4.SetLineColor(kBlue);
  h_reco_cos_theta_4.SetLineWidth(2);
  TH1F h_reco_cos_theta_5("h_reco_cos_theta_5","5 MeV;reco_cos#theta",200,-1,1);
  h_reco_cos_theta_5.SetLineColor(kGreen);
  h_reco_cos_theta_5.SetLineWidth(2);

  TH1F h_reco2_cos_theta_1("h_reco2_cos_theta_1","1 MeV;reco_cos#theta",200,-1,1);
  h_reco2_cos_theta_1.SetLineColor(kGray);
  h_reco2_cos_theta_1.SetLineWidth(2);
  TH1F h_reco2_cos_theta_2("h_reco2_cos_theta_2","2 MeV;reco_cos#theta",200,-1,1);
  h_reco2_cos_theta_2.SetLineColor(kBlack);
  h_reco2_cos_theta_2.SetLineWidth(2);
  TH1F h_reco2_cos_theta_3("h_reco2_cos_theta_3","3 MeV;reco_cos#theta",200,-1,1);
  h_reco2_cos_theta_3.SetLineColor(kRed);
  h_reco2_cos_theta_3.SetLineWidth(2);
  TH1F h_reco2_cos_theta_4("h_reco2_cos_theta_4","4 MeV;reco_cos#theta",200,-1,1);
  h_reco2_cos_theta_4.SetLineColor(kBlue);
  h_reco2_cos_theta_4.SetLineWidth(2);
  TH1F h_reco2_cos_theta_5("h_reco2_cos_theta_5","5 MeV;reco_cos#theta",200,-1,1);
  h_reco2_cos_theta_5.SetLineColor(kGreen);
  h_reco2_cos_theta_5.SetLineWidth(2);

  TH2F h_reco_vtx_x__vs__vtx_x("h_reco_vtx_x__vs__vtx_x","h_reco_vtx_x__vs__vtx_x; true vtx x; reco vtx x",100,-6000,6000,100,-6000,6000);
  TH2F h_reco_vtx_y__vs__vtx_y("h_reco_vtx_y__vs__vtx_y","h_reco_vtx_y__vs__vtx_y; true vtx y; reco vtx y",100,-6000,6000,100,-6000,6000);
  TH2F h_reco_vtx_z__vs__vtx_z("h_reco_vtx_z__vs__vtx_z","h_reco_vtx_z__vs__vtx_z; true vtx z; reco vtx z",100,-6000,6000,100,-6000,6000);

  clog << " n entries " << angles_tree->GetEntries() << endl;

  for(int ientry = 0; ientry < angles_tree->GetEntries(); ientry++){
    angles_tree->GetEntry(ientry);

    if( n_hits == 0 ) continue;

    smear_vtx();

    h_reco_vtx_x__vs__vtx_x.Fill(reco_vtx_x, vtx_x);
    h_reco_vtx_y__vs__vtx_y.Fill(reco_vtx_y, vtx_y);
    h_reco_vtx_z__vs__vtx_z.Fill(reco_vtx_z, vtx_z);

    calculate_averages();

    calculate_reco_dir();

    h_reco_dir_x__vs__dir_x.Fill(dir_x, reco_dir_x);
    h_reco_dir_y__vs__dir_y.Fill(dir_y, reco_dir_y);
    h_reco_dir_z__vs__dir_z.Fill(dir_z, reco_dir_z);
    
    if( kinetic == 3. ){
      h_residual_dir_x.Fill(reco_dir_x - dir_x);
      h_residual_dir_y.Fill(reco_dir_y - dir_y);
      h_residual_dir_z.Fill(reco_dir_z - dir_z);
    }

    h_reco_dir_x.Fill(reco_dir_x);
    h_reco_dir_y.Fill(reco_dir_y);
    h_reco_dir_z.Fill(reco_dir_z);
    
    if( !calculate_reco2_dir() ) continue;
  
    h_reco2_dir_x__vs__dir_x.Fill(dir_x, reco2_dir_x);
    h_reco2_dir_y__vs__dir_y.Fill(dir_y, reco2_dir_y);
    h_reco2_dir_z__vs__dir_z.Fill(dir_z, reco2_dir_z);
    
    if( kinetic == 3. ){
      h_residual2_dir_x.Fill(reco2_dir_x - dir_x);
      h_residual2_dir_y.Fill(reco2_dir_y - dir_y);
      h_residual2_dir_z.Fill(reco2_dir_z - dir_z);
    }

    h_reco2_dir_x.Fill(reco2_dir_x);
    h_reco2_dir_y.Fill(reco2_dir_y);
    h_reco2_dir_z.Fill(reco2_dir_z);
    
  
    n_hits_in_range=0;
    n_hits_in_range2=0;

    for(int ihit=0; ihit < n_hits; ihit++){
      dist_x = hit_x[ihit] - reco_vtx_x;
      dist_y = hit_y[ihit] - reco_vtx_y;
      dist_z = hit_z[ihit] - reco_vtx_z;
      dist = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));

      tof = dist/c;
      tcorr = hit_t[ihit] - tof;
      h_t.Fill(hit_t[ihit]);
      h_tcorr.Fill(tcorr);

      cos_theta = (dist_x*dir_x + dist_y*dir_y + dist_z*dir_z)/dist;
      theta = acos(cos_theta)*180./pi;

      reco_cos_theta = (dist_x*reco_dir_x + dist_y*reco_dir_y + dist_z*reco_dir_z)/dist;

      reco2_cos_theta = (dist_x*reco2_dir_x + dist_y*reco2_dir_y + dist_z*reco2_dir_z)/dist;


      if( kinetic == 1. ){
	h_cos_theta_1->Fill(cos_theta);
	h_theta_1.Fill(theta);
	h_theta_from_theta0_1.Fill(theta - theta0);
	h_reco_cos_theta_1.Fill(reco_cos_theta);
	h_reco2_cos_theta_1.Fill(reco2_cos_theta);
      }else if( kinetic == 2. ){
	h_cos_theta_2->Fill(cos_theta);
	h_theta_2.Fill(theta);
	h_theta_from_theta0_2.Fill(theta - theta0);
	h_reco_cos_theta_2.Fill(reco_cos_theta);
	h_reco2_cos_theta_2.Fill(reco2_cos_theta);
      }else if( kinetic == 3. ){
	h_cos_theta_3->Fill(cos_theta);
	h_theta_3.Fill(theta);
	h_theta_from_theta0_3.Fill(theta - theta0);
	h_reco_cos_theta_3.Fill(reco_cos_theta);
	h_reco2_cos_theta_3.Fill(reco2_cos_theta);
      }else if( kinetic == 4. ){
	h_cos_theta_4->Fill(cos_theta);
	h_theta_4.Fill(theta);
	h_theta_from_theta0_4.Fill(theta - theta0);
	h_reco_cos_theta_4.Fill(reco_cos_theta);
	h_reco2_cos_theta_4.Fill(reco2_cos_theta);
      }else if( kinetic == 5. ){
	h_cos_theta_5->Fill(cos_theta);
	h_theta_5.Fill(theta);
	h_theta_from_theta0_5.Fill(theta - theta0);
	h_reco_cos_theta_5.Fill(reco_cos_theta);
	h_reco2_cos_theta_5.Fill(reco2_cos_theta);
      }else{
	clog << " unknown kinetic " << kinetic << endl;
      }

      if( reco_cos_theta >= _cos_theta_min && reco_cos_theta <= _cos_theta_max ){
	n_hits_in_range++;
      }
      
      if( reco2_cos_theta >= _cos_theta_min2 && reco2_cos_theta <= _cos_theta_max2 ){
	n_hits_in_range2++;
      }


    }

    residual_x0 = get_residual_position(vtx_x, vtx_y, vtx_z, vtx_t, 
					average_x,
					average_y,
					average_z,
					average_t,
					average_x2,
					average_y2,
					average_z2,
					average_t2,
					average_x3,
					average_y2x,
					average_z2x,
					average_t2x,
					average_xy,
					average_xz,
					average_xt
					);
    h_analytic_residual_x0.Fill(residual_x0);

    residual_d = get_residual_d(vtx_x, vtx_y, vtx_z, vtx_t, 
					average_x,
					average_y,
					average_z,
					average_t,
					average_x2,
					average_y2,
					average_z2,
					average_t2
					);
    h_analytic_residual_d.Fill(residual_d);
  }
  
  get_range(&h_reco_cos_theta_1, &cos_theta_min[0], &cos_theta_max[0]);
  get_range(&h_reco_cos_theta_2, &cos_theta_min[1], &cos_theta_max[1]);
  get_range(&h_reco_cos_theta_3, &cos_theta_min[2], &cos_theta_max[2]);
  get_range(&h_reco_cos_theta_4, &cos_theta_min[3], &cos_theta_max[3]);
  get_range(&h_reco_cos_theta_5, &cos_theta_min[4], &cos_theta_max[4]);
  
  get_range(&h_reco2_cos_theta_1, &cos_theta_min2[0], &cos_theta_max2[0]);
  get_range(&h_reco2_cos_theta_2, &cos_theta_min2[1], &cos_theta_max2[1]);
  get_range(&h_reco2_cos_theta_3, &cos_theta_min2[2], &cos_theta_max2[2]);
  get_range(&h_reco2_cos_theta_4, &cos_theta_min2[3], &cos_theta_max2[3]);
  get_range(&h_reco2_cos_theta_5, &cos_theta_min2[4], &cos_theta_max2[4]);

  h_hits_fraction_1.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min[0], cos_theta_max[0]));
  h_hits_fraction_2.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min[1], cos_theta_max[1]));
  h_hits_fraction_3.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min[2], cos_theta_max[2]));
  h_hits_fraction_4.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min[3], cos_theta_max[3]));
  h_hits_fraction_5.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min[4], cos_theta_max[4]));
  
  h_hits_fraction2_1.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min2[0], cos_theta_max2[0]));
  h_hits_fraction2_2.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min2[1], cos_theta_max2[1]));
  h_hits_fraction2_3.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min2[2], cos_theta_max2[2]));
  h_hits_fraction2_4.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min2[3], cos_theta_max2[3]));
  h_hits_fraction2_5.GetXaxis()->SetTitle(Form("fraction of hits with cos#theta #in(%.2f, %.2f)", cos_theta_min2[4], cos_theta_max2[4]));

  for(int ientry = 0; ientry < angles_tree->GetEntries(); ientry++){
    angles_tree->GetEntry(ientry);

    if( n_hits == 0 ) continue;

    smear_vtx();

    calculate_averages();

    calculate_reco_dir();

    if( !calculate_reco2_dir() ) continue;
  
    n_hits_in_range=0;
    n_hits_in_range2=0;

    int ikinetic = (int)(kinetic - 1.);

    for(int ihit=0; ihit < n_hits; ihit++){
      dist_x = hit_x[ihit] - reco_vtx_x;
      dist_y = hit_y[ihit] - reco_vtx_y;
      dist_z = hit_z[ihit] - reco_vtx_z;
      dist = sqrt(pow(dist_x,2) + pow(dist_y,2) + pow(dist_z,2));

      reco_cos_theta = (dist_x*reco_dir_x + dist_y*reco_dir_y + dist_z*reco_dir_z)/dist;

      reco2_cos_theta = (dist_x*reco2_dir_x + dist_y*reco2_dir_y + dist_z*reco2_dir_z)/dist;


      if( reco_cos_theta >= cos_theta_min[ikinetic] && reco_cos_theta <= cos_theta_max[ikinetic] ){
	n_hits_in_range++;
      }
      
      if( reco2_cos_theta >= cos_theta_min2[ikinetic] && reco2_cos_theta <= cos_theta_max2[ikinetic] ){
	n_hits_in_range2++;
      }
    }
    
    if( kinetic == 1. ){
      h_hits_fraction_1.Fill(n_hits_in_range/(double)n_hits);
      h_hits_fraction2_1.Fill(n_hits_in_range2/(double)n_hits);
    }else if( kinetic == 2. ){
      h_hits_fraction_2.Fill(n_hits_in_range/(double)n_hits);
      h_hits_fraction2_2.Fill(n_hits_in_range2/(double)n_hits);
    }else if( kinetic == 3. ){
      h_hits_fraction_3.Fill(n_hits_in_range/(double)n_hits);
      h_hits_fraction2_3.Fill(n_hits_in_range2/(double)n_hits);
    }else if( kinetic == 4. ){
      h_hits_fraction_4.Fill(n_hits_in_range/(double)n_hits);
      h_hits_fraction2_4.Fill(n_hits_in_range2/(double)n_hits);
    }else if( kinetic == 5. ){
      h_hits_fraction_5.Fill(n_hits_in_range/(double)n_hits);
      h_hits_fraction2_5.Fill(n_hits_in_range2/(double)n_hits);
    }

  }


  output->Write();

  return 0;

}

double get_residual_position(double x0, double y0, double z0, double t0, 
			     double x,
			     double y,
			     double z,
			     double t,
			     double x2,
			     double y2,
			     double z2,
			     double t2,
			     double x3,
			     double y2x,
			     double z2x,
			     double t2x,
			     double xy,
			     double xz,
			     double xt
			     ){

  double v2 = pow(c,2);

  double residual =
    - pow(x0,3)
    + x*(3.*pow(x0,2) + pow(y0,2) + pow(z0,2))
    - x0*(3.*x2 + y2 + z2 - v2*t2)
    + x3 + y2x + z2x - v2*t2x
    - 2.*y0*xy - 2.*z0*xz + 2.*x0*(y0*y + z0*z)
    - (pow(y0,2) + pow(z0,2))*x0
    + v2*(pow(t0,2)*(x0 - x) + 2.*t0*xt - 2.*t0*x0*t)
    ;

  //  clog << " x0 " << x0 << " y0 " << y0 << " z0 " << z0 << " t0 " << t0 << " x " << x << " y " << y << " z " << z << " t " << t << " x2 " << x2 << " y2 " << y2 << " z2 " << z2 << " t2 " << t2 << " x3 " << x3 << " y2x " << y2x << " z2x " << z2x << " t2x " << t2x << " xy " << xy << " xz " << xz << " xt " << xt << endl;

  //  clog << " residual " << residual << endl;

  return residual;

}


double get_residual_d(double x0, double y0, double z0, double t0, 
			     double x,
			     double y,
			     double z,
			     double t,
			     double x2,
			     double y2,
			     double z2,
			     double t2
			     ){

  double v2 = pow(c,2);

  double residual =
    x2 + pow(x0,2) - 2.*x0*x
    + y2 + pow(y0,2) - 2.*y0*y
    + z2 + pow(z0,2) - 2.*z0*z
    - v2*( t2 + pow(t0,2) - 2.*t0*t)
    ;

  //clog << " x0 " << x0 << " y0 " << y0 << " z0 " << z0 << " t0 " << t0 << " x " << x << " y " << y << " z " << z << " t " << t << " x2 " << x2 << " y2 " << y2 << " z2 " << z2 << " t2 " << t2 << endl;

  //clog << " residual " << residual << endl;

  return residual;

}


double smear(double center, double sigma){

  double val=gRandom->Gaus(center, sigma);
  return val;

}

void smear_vtx(){

  reco_vtx_x = smear(vtx_x, vtx_sigma);
  reco_vtx_y = smear(vtx_y, vtx_sigma);
  reco_vtx_z = smear(vtx_z, vtx_sigma);

  return;
}

void calculate_averages(){

  average_x=0;
  average_y=0;
  average_z=0;
  average_d=0;
  average_t=0;
  average_x2=0;
  average_y2=0;
  average_z2=0;
  average_t2=0;
  average_x3=0;
  average_y3=0;
  average_z3=0;
  average_t3=0;
  average_y2x=0;
  average_z2x=0;
  average_t2x=0;
  average_x2y=0;
  average_z2y=0;
  average_t2y=0;
  average_x2z=0;
  average_y2z=0;
  average_t2z=0;
  average_x2t=0;
  average_y2t=0;
  average_z2t=0;
  average_xy=0;
  average_xz=0;
  average_xt=0;
  average_yz=0;
  average_yt=0;
  average_zt=0;
  average_dx=0;
  average_dy=0;
  average_dz=0;
  average_dx2=0;
  average_dy2=0;
  average_dz2=0;
  average_dxy=0;
  average_dxz=0;
  average_dyz=0;
  
  for(int ihit=0; ihit < n_hits; ihit++){

    reco_dist_x = hit_x[ihit] - reco_vtx_x;
    reco_dist_y = hit_y[ihit] - reco_vtx_y;
    reco_dist_z = hit_z[ihit] - reco_vtx_z;
    reco_dist = sqrt(pow(reco_dist_x,2) + pow(reco_dist_y,2) + pow(reco_dist_z,2));
    reco_dx = reco_dist_x/reco_dist;
    reco_dy = reco_dist_y/reco_dist;
    reco_dz = reco_dist_z/reco_dist;

    average_dx += reco_dx;
    average_dy += reco_dy;
    average_dz += reco_dz;
    average_dx2 += pow(reco_dx,2);
    average_dy2 += pow(reco_dy,2);
    average_dz2 += pow(reco_dz,2);
    average_dxy += reco_dx*reco_dy;
    average_dxz += reco_dx*reco_dz;
    average_dyz += reco_dy*reco_dz;
    
    average_x += hit_x[ihit];
    average_y += hit_y[ihit];
    average_z += hit_z[ihit];
    average_t += hit_t[ihit];
    average_x2 += pow(hit_x[ihit],2);
    average_y2 += pow(hit_y[ihit],2);
    average_z2 += pow(hit_z[ihit],2);
    average_t2 += pow(hit_t[ihit],2);
    average_x3 += pow(hit_x[ihit],3);
    average_y3 += pow(hit_y[ihit],3);
    average_z3 += pow(hit_z[ihit],3);
    average_t3 += pow(hit_t[ihit],3);
    average_y2x += pow(hit_y[ihit],2)*hit_x[ihit];
    average_z2x += pow(hit_z[ihit],2)*hit_x[ihit];
    average_t2x += pow(hit_t[ihit],2)*hit_x[ihit];
    average_x2y += pow(hit_x[ihit],2)*hit_y[ihit];
    average_z2y += pow(hit_z[ihit],2)*hit_y[ihit];
    average_t2y += pow(hit_t[ihit],2)*hit_y[ihit];
    average_x2z += pow(hit_x[ihit],2)*hit_z[ihit];
    average_y2z += pow(hit_y[ihit],2)*hit_z[ihit];
    average_t2z += pow(hit_t[ihit],2)*hit_z[ihit];
    average_x2t += pow(hit_x[ihit],2)*hit_t[ihit];
    average_y2t += pow(hit_y[ihit],2)*hit_t[ihit];
    average_z2t += pow(hit_z[ihit],2)*hit_t[ihit];
    average_xy += hit_x[ihit]*hit_y[ihit];
    average_xz += hit_x[ihit]*hit_z[ihit];
    average_xt += hit_x[ihit]*hit_t[ihit];
    average_yz += hit_y[ihit]*hit_z[ihit];
    average_yt += hit_y[ihit]*hit_t[ihit];
    average_zt += hit_z[ihit]*hit_t[ihit];

  }
    
  average_x /= (double)n_hits;
  average_y /= (double)n_hits;
  average_z /= (double)n_hits;
  average_t /= (double)n_hits;
  average_x2 /= (double)n_hits;
  average_y2 /= (double)n_hits;
  average_z2 /= (double)n_hits;
  average_t2 /= (double)n_hits;
  average_x3 /= (double)n_hits;
  average_y3 /= (double)n_hits;
  average_z3 /= (double)n_hits;
  average_t3 /= (double)n_hits;
  average_y2x /= (double)n_hits;
  average_z2x /= (double)n_hits;
  average_t2x /= (double)n_hits;
  average_x2y /= (double)n_hits;
  average_z2y /= (double)n_hits;
  average_t2y /= (double)n_hits;
  average_x2z /= (double)n_hits;
  average_y2z /= (double)n_hits;
  average_t2z /= (double)n_hits;
  average_x2t /= (double)n_hits;
  average_y2t /= (double)n_hits;
  average_z2t /= (double)n_hits;
  average_xy /= (double)n_hits;
  average_xz /= (double)n_hits;
  average_xt /= (double)n_hits;
  average_yz /= (double)n_hits;
  average_yt /= (double)n_hits;
  average_zt /= (double)n_hits;
  average_dx /= (double)n_hits;
  average_dy /= (double)n_hits;
  average_dz /= (double)n_hits;
  average_dx2 /= (double)n_hits;
  average_dy2 /= (double)n_hits;
  average_dz2 /= (double)n_hits;
  average_dxy /= (double)n_hits;
  average_dxz /= (double)n_hits;
  average_dyz /= (double)n_hits;
  
  return;
}

void calculate_reco_dir(){

  average_d = sqrt(pow(average_x - reco_vtx_x,2) + pow(average_y - reco_vtx_y,2) + pow(average_z - reco_vtx_z,2));
  reco_dir_x = (average_x - reco_vtx_x)/average_d;
  reco_dir_y = (average_y - reco_vtx_y)/average_d;
  reco_dir_z = (average_z - reco_vtx_z)/average_d;
    

  return;
}

bool calculate_reco2_dir(){

  B[0] = cos_theta0*average_dx;
  B[1] = cos_theta0*average_dy;
  B[2] = cos_theta0*average_dz;

  for(int i=0; i<3; i++)
    AB[i] = B[i];
    
  A[0][0] = average_dx2;
  A[0][1] = average_dxy;
  A[0][2] = average_dxz;
  A[1][0] = A[0][1];
  A[1][1] = average_dy2;
  A[1][2] = average_dyz;
  A[2][0] = A[0][2];
  A[2][1] = A[1][2];
  A[2][2] = average_dz2;

  detA = A[0][0]*(A[1][1]*A[2][2] - pow(A[1][2],2)) -
    A[0][1]*(A[1][0]*A[2][2] - A[1][2]*A[2][0]) +
    A[0][2]*(A[1][0]*A[2][1] - A[1][1]*A[2][0]) ;
    
  for(int i=0; i<3; i++)
    for(int j=0; j<3; j++)
      AA[i][j] = A[i][j];

  detAA = AA.Determinant();

    
  // clog << " B " << endl;
  // clog << B[0] << " " << B[1] << " " << B[2] << " " << endl;
  // clog << " AB " << endl;
  // clog << AB[0] << " " << AB[1] << " " << AB[2] << " " << endl;
  // clog << " A " << endl;
  // clog << A[0][0] << " " << A[0][1] << " " << A[0][2] << " " << endl;
  // clog << A[1][0] << " " << A[1][1] << " " << A[1][2] << " " << endl;
  // clog << A[2][0] << " " << A[2][1] << " " << A[2][2] << " " << endl;
  // clog << " detA " << detA << endl;
  // clog << " AA " << endl;
  // clog << AA[0][0] << " " << AA[0][1] << " " << AA[0][2] << " " << endl;
  // clog << AA[1][0] << " " << AA[1][1] << " " << AA[1][2] << " " << endl;
  // clog << AA[2][0] << " " << AA[2][1] << " " << AA[2][2] << " " << endl;
  // clog << " detA " << detAA << endl;
    
  Ai[0][0] = (A[1][1]*A[2][2] - pow(A[1][2],2))/detA;
  Ai[0][1] = -(A[0][1]*A[2][2] - A[0][2]*A[2][1])/detA;
  Ai[0][2] = (A[0][1]*A[1][2] - A[0][2]*A[1][1])/detA;
  Ai[1][0] = Ai[0][1];
  Ai[1][1] = (A[0][0]*A[2][2] - pow(A[0][2],2))/detA;
  Ai[1][2] = -(A[0][0]*A[1][2] - A[0][2]*A[1][0])/detA;
  Ai[2][0] = Ai[0][2];
  Ai[2][1] = Ai[1][2];
  Ai[2][2] = (A[0][0]*A[1][1] - pow(A[0][1],2))/detA;
    
  AAi = AA.Invert(&detAA);

  if( fabs(detAA) < 1.e-10 ){
    clog << " not invertible; skip " << endl;
    return false;
  }

  // clog << " Ai " << endl;
  // clog << Ai[0][0] << " " << Ai[0][1] << " " << Ai[0][2] << " " << endl;
  // clog << Ai[1][0] << " " << Ai[1][1] << " " << Ai[1][2] << " " << endl;
  // clog << Ai[2][0] << " " << Ai[2][1] << " " << Ai[2][2] << " " << endl;
  // clog << " AAi " << endl;
  // clog << AAi[0][0] << " " << AAi[0][1] << " " << AAi[0][2] << " " << endl;
  // clog << AAi[1][0] << " " << AAi[1][1] << " " << AAi[1][2] << " " << endl;
  // clog << AAi[2][0] << " " << AAi[2][1] << " " << AAi[2][2] << " " << endl;
    
  // for(int i=0; i<3; i++){
  // 	for(int j=0; j<3; j++){
  // 	  I[i][j] = A[i][0]*Ai[0][j] + A[i][1]*Ai[1][j] + A[i][2]*Ai[2][j];
  // 	}
  // }
    
  // clog << " I " << endl;
  // clog << I[0][0] << " " << I[0][1] << " " << I[0][2] << " " << endl;
  // clog << I[1][0] << " " << I[1][1] << " " << I[1][2] << " " << endl;
  // clog << I[2][0] << " " << I[2][1] << " " << I[2][2] << " " << endl;

  reco2_dir_x = Ai[0][0]*B[0] + Ai[0][1]*B[1] + Ai[0][2]*B[2];
  reco2_dir_y = Ai[1][0]*B[0] + Ai[1][1]*B[1] + Ai[1][2]*B[2];
  reco2_dir_z = Ai[2][0]*B[0] + Ai[2][1]*B[1] + Ai[2][2]*B[2];
  // reco2_dir_modulus = sqrt(pow(reco2_dir_x,2) + pow(reco2_dir_y,2) + pow(reco2_dir_z,2));
  // reco2_dir_x /= reco2_dir_modulus;
  // reco2_dir_y /= reco2_dir_modulus;
  // reco2_dir_z /= reco2_dir_modulus;

  // clog << " reco2_dir_x " << reco2_dir_x << " reco2_dir_y " << reco2_dir_y << " reco2_dir_z " << reco2_dir_z << endl;

  A_reco2_dir = AAi*AB;

  // clog << " Areco2_dir_x " << A_reco2_dir[0] << " Areco2_dir_y " << A_reco2_dir[1] << " Areco2_dir_z " << A_reco2_dir[2] << endl;

  ////////////////////////////////////////////////
  ////    (A + eta) u = B                     ////
  ////    u = P - eta Q
  ////    P = Ai B
  ////    Q = Ai Ai B
  ////    1 = P2 - 2 eta PQ + eta2 Q2
  ////    eta2 Q2 - eta 2PQ + P2 - 1 = 0
  ////    Delta = (PQ)^2 + Q2(1 - P2)
  ////    eta = (PQ +- sq(Delta))/Q2
  
  

  TVectorD P(3);
  TVectorD Q(3);

  P = AAi*AB;
  Q = AAi*P;

  double P2 = P*P;
  double Q2 = Q*Q;
  double PQ = P*Q;

  double Delta = pow(PQ,2) + Q2*(1. - P2);

  if( Delta < 0. ){
    clog << " Delta " << Delta << endl;
    return false;
  }

  double eta1 = (PQ + sqrt(Delta))/Q2;
  double eta2 = (PQ - sqrt(Delta))/Q2;

  TVectorD u1(3);
  TVectorD u2(3);

  u1 = P - eta1*Q;
  u2 = P - eta2*Q;

  // clog << " u1 " << u1[0] << ", " << u1[1] << ", " << u1[2] << " norm " << u1*u1 << endl;
  // clog << " u2 " << u2[0] << ", " << u2[1] << ", " << u2[2] << " norm " << u2*u2 << endl;

  double chi2_1 = (AA*u1)*u1 - 2.* u1*AB;
  double chi2_2 = (AA*u2)*u2 - 2.* u2*AB;

  //  clog << " eta1 " << eta1 << " chi2_1 " << chi2_1 << " eta2 " << eta2 << " chi2_2 " << chi2_2 << endl;

  // clog << " Areco2_dir_x " << A_reco2_dir[0] << " Areco2_dir_y " << A_reco2_dir[1] << " Areco2_dir_z " << A_reco2_dir[2] << " norm " << A_reco2_dir*A_reco2_dir << endl;

  if( chi2_1 < chi2_2 ) A_reco2_dir = u1;
  else A_reco2_dir = u2;

  reco2_dir_x = A_reco2_dir[0];
  reco2_dir_y = A_reco2_dir[1];
  reco2_dir_z = A_reco2_dir[2];

  // clog << " reco2_dir_x " << reco2_dir_x << " reco2_dir_y " << reco2_dir_y << " reco2_dir_z " << reco2_dir_z << endl;

  // clog << " Areco2_dir_x " << A_reco2_dir[0] << " Areco2_dir_y " << A_reco2_dir[1] << " Areco2_dir_z " << A_reco2_dir[2] << " norm " << A_reco2_dir*A_reco2_dir << endl;

  return true;
}

void get_range(TH1F * h, double *m, double *M){

  int imax = h->GetMaximumBin();
  int i1 = imax;
  int i2 = imax;
  double total_integral = h->Integral();
  double partial_integral = h->Integral(i1,i2);
  double ratio = partial_integral/total_integral;
  while( ratio < 0.5 ){
    if( i1 > 0 )
      i1 --;
    if( i2 < h->GetXaxis()->GetNbins() )
      i2 ++;
    partial_integral = h->Integral(i1,i2);
    ratio = partial_integral/total_integral;
    
  }

  *m = h->GetBinCenter(i1);
  *M = h->GetBinCenter(i2);

  return;


}
