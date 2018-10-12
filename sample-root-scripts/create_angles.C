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

using namespace std;

struct vec
{
  double x;
  double y;
  double z;
};  


struct particle
{
  vec vertex;
  vec direction;
  double K;
};  

struct hit
{
  vec position;
  double t;
};  

struct pmt
{
  vec position;
};  

#define NEVENTS 10000
#define NHITS 200
#define NPMTS 50000

int main(){

  TFile * output = new TFile("angles.root","RECREATE");

  TTree angles_tree("angles_tree","angles tree");
  Float_t vtx_x, vtx_y, vtx_z;
  Float_t dir_x, dir_y, dir_z;
  Float_t kinetic;
  Int_t n_hits;
  Float_t hit_x[NHITS];
  Float_t hit_y[NHITS];
  Float_t hit_z[NHITS];
  Float_t hit_t[NHITS];
  angles_tree.Branch("vtx_x",&vtx_x,"vtx_x/F");
  angles_tree.Branch("vtx_y",&vtx_y,"vtx_y/F");
  angles_tree.Branch("vtx_z",&vtx_z,"vtx_z/F");
  angles_tree.Branch("dir_x",&dir_x,"dir_x/F");
  angles_tree.Branch("dir_y",&dir_y,"dir_y/F");
  angles_tree.Branch("dir_z",&dir_z,"dir_z/F");
  angles_tree.Branch("kinetic",&kinetic,"kinetic/F");
  angles_tree.Branch("n_hits",&n_hits,"n_hits/I");
  angles_tree.Branch("hit_x",hit_x,"hit_x[n_hits]/F");
  angles_tree.Branch("hit_y",hit_y,"hit_y[n_hits]/F");
  angles_tree.Branch("hit_z",hit_z,"hit_z[n_hits]/F");
  angles_tree.Branch("hit_t",hit_t,"hit_t[n_hits]/F");

  ifstream vertices_file;
  vertices_file.open("outputs/eminus_1MeV_random_4pi_000_40per/vertices.txt");

  ifstream directions_file;
  directions_file.open("outputs/eminus_1MeV_random_4pi_000_40per/directions.txt");

  particle particles[NEVENTS];

  vec vertex;
  vec direction;
  double K;
  particle p;
  int nparticles=0;


  while(!vertices_file.eof()){
    vertices_file >> vertex.x;
    vertices_file >> vertex.y;
    vertices_file >> vertex.z;
    p.vertex = vertex;
    directions_file >> K;
    p.K = K;
    directions_file >> direction.x;
    directions_file >> direction.y;
    directions_file >> direction.z;
    p.direction = direction;

    particles[nparticles] = p;

    nparticles++;

  } 
  nparticles --;
  vertices_file.close(); 
  directions_file.close(); 


  ifstream pmts_file;
  pmts_file.open("outputs/eminus_1MeV_random_4pi_000_40per/all_pmts.txt");

  pmt pmts[NPMTS];
  pmt apmt;

  vec position;
  int id;
  int npmts=0;


  while(!pmts_file.eof()){
    pmts_file >> id;
    pmts_file >> position.x;
    pmts_file >> position.y;
    pmts_file >> position.z;
    apmt.position = position;

    pmts[npmts] = apmt;

    npmts++;

  } 
  npmts--;
  pmts_file.close(); 

  double t;
  hit ahit;
  ifstream hits_file;
  int nhits;

  for(int ipart=0; ipart<nparticles; ipart++){

    hits_file.open(Form("outputs/eminus_1MeV_random_4pi_000_40per/all_hits_%d.txt", ipart+1));
    if (hits_file.is_open()){

      hit hits[NHITS];
      nhits=0;
      while(!hits_file.eof()){
	hits_file >> id;
	hits_file >> t;

	ahit.t = t;
	position = pmts[id-1].position;
	ahit.position = position;
	
	hits[nhits] = ahit;
	
	nhits++;
	
      } 
      nhits--;
      hits_file.close(); 

      vtx_x=particles[ipart].vertex.x;
      vtx_y=particles[ipart].vertex.y;
      vtx_z=particles[ipart].vertex.z;
      dir_x=particles[ipart].direction.x;
      dir_y=particles[ipart].direction.y;
      dir_z=particles[ipart].direction.z;
      kinetic=particles[ipart].K;
      n_hits=nhits;
      for(int ihit=0; ihit<nhits; ihit++){
	hit_x[ihit] = hits[ihit].position.x;
	hit_y[ihit] = hits[ihit].position.y;
	hit_z[ihit] = hits[ihit].position.z;
	hit_t[ihit] = hits[ihit].t;
      }

      angles_tree.Fill();

    }
  }


  output->Write();

  return 0;

}

