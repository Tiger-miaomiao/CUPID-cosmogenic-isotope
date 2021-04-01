#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4GeneralParticleSource.hh"

#include "G4RunManager.hh"
#include "G4ios.hh"

#include "statisticsLogger.hh"
#include "primaryGeneratorAction.hh"
#include "Randomize.hh"
#include "spectra.h"
#include "TH1D.h"
#include "TF1.h"

/*PrimaryGeneratorAction::PrimaryGeneratorAction(G4String source_name):
  G4VUserPrimaryGeneratorAction(),
  particleGun(0)
{

  sname=source_name;

  G4int nPrimaries = 1;
  particleGun = new G4ParticleGun(nPrimaries);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();

  G4String particleName="neutron";
  //if(sname=="gamma"){particleName = "gamma";}//MLCC "gamma" "neutron" "mu-"
  //if(sname=="neutron"){particleName = "neutron";}//MLCC "gamma" "neutron" "mu-"
  //if(sname=="muon"){particleName = "mu-";}//MLCC "gamma" "neutron" "mu-"

  G4ParticleDefinition* particle = particleTable->FindParticle(particleName);

  StatisticsLogger* logger = StatisticsLogger::GetInstance();
  logger->SetPrimaryPDGEncoding(particle->GetPDGEncoding());

  G4ThreeVector momentumDirection = G4ThreeVector(0, 0, 0);//MLCC
  G4ThreeVector position = G4ThreeVector(0,0,0);//MLCC

  particleGun->SetParticleEnergy(0);
  particleGun->SetParticlePosition(position);
  particleGun->SetParticleDefinition(particle);
  particleGun->SetParticleMomentumDirection(momentumDirection);

  }*/
PrimaryGeneratorAction::PrimaryGeneratorAction():
  G4VUserPrimaryGeneratorAction(),
  particleGun(0)
{
  particleGun = new G4GeneralParticleSource();
}

/*PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
  }*/
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete particleGun;
}

/*void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  TH1D* hrandom;
  TF1* fmrandom;

  if(sname=="gamma")
    {
      hrandom=new TH1D("hrandom","hrandom",10,Eg);
      for(int i=0;i<10;i++){hrandom->SetBinContent(i+1,Rg[i]);}
    }

  if(sname=="neutron")
    {
      hrandom=new TH1D("hrandom","hrandom",4,En);
      for(int i=0;i<4;i++){hrandom->SetBinContent(i+1,Rn[i]);}
    }
  if(sname=="muon")
    {
      fmrandom=new TF1("fmrandom","exp(-0.4*6.72*2.77)*pow((x+693.0*(1-exp(-0.4*6.72))),-3.77)",1.0,1000.0);
      //hrandom=new TH1D("hrandom","hrandom",10,Em);
      //for(int i=0;i<10;i++){hrandom->SetBinContent(i+1,Rm[i]);}
    }


 G4double Ep=0;
 G4double Gx=0;
 G4double Gy=0;
 G4double Gz=0;
 G4double Px=0;
 G4double Py=0;
 G4double Pz=0;

 //G4cout<<"sname="<<sname<<G4endl;

 if((sname=="gamma"||sname=="neutron")&&hrandom){SetParticle(Gx,Gy,Gz,Px,Py,Pz,Ep,hrandom);}
  if(sname=="muon"&&fmrandom){SetParticle(Gx,Gy,Gz,Px,Py,Pz,Ep,fmrandom);}


//if(Ep==0)
//{
//if(hrandom){hrandom->Delete();}
//if(fmrandom){fmrandom->Delete();}
//return;
//}

if(Ep==0&&hrandom)
{
hrandom->Delete();
return;
}

if(hrandom){hrandom->Delete();}
  //G4cout<<"energy="<<energy<<G4endl;

G4double energy = Ep*CLHEP::MeV;

//G4cout<<"primary energy="<<energy<<G4endl;

 G4ThreeVector position = G4ThreeVector(Gx*CLHEP::m,Gy*CLHEP::m,Gz*CLHEP::m);//MLCC
 G4ThreeVector  momentumDirection=G4ThreeVector(Px,Py,Pz);//MLCC

 //G4cout<<"Px="<<Px<<G4endl;

StatisticsLogger* logger = StatisticsLogger::GetInstance();
 logger->FillPrimarydistHistogram(Gx*100*CLHEP::cm,Gy*100*CLHEP::cm,Gz*100*CLHEP::cm);

  particleGun->SetParticleEnergy(energy);
  particleGun->SetParticlePosition(position);
  particleGun->SetParticleMomentumDirection(momentumDirection);
  particleGun->GeneratePrimaryVertex(anEvent);

  //if(hrandom){hrandom->Delete();}else if(fmrandom){fmrandom->Delete();}else{return;}

  //G4cout<<"Here3"<<G4endl;

}*/
void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  particleGun -> GeneratePrimaryVertex(anEvent);

  G4ThreeVector source_position=particleGun->GetParticlePosition();

StatisticsLogger* logger = StatisticsLogger::GetInstance();
 logger->FillPrimarydistHistogram(source_position.x(),source_position.y(),source_position.z());

}

/*
void PrimaryGeneratorAction::SetParticle(G4double &x,G4double &y,G4double &z,G4double &px,G4double &py,G4double &pz,G4double &e,TH1D* h)
{

//initialize source position in hemisphere, CJPL hole radius = 6.0m

  G4double u=2.0*G4UniformRand()-1.0;
  G4double v=2.0*G4UniformRand()-1.0;
  if((u*u+v*v)>1.0)return;
  x=6.0*(2.0*u*sqrt(1.0-(u*u+v*v))); 
  y=6.0*(2.0*v*sqrt(1.0-(u*u+v*v)));
  z=6.0*(1-2.0*(u*u+v*v));


  e=h->GetRandom(); 
  //e=e+(G4UniformRand()-0.5)*(h->GetXaxis()->GetBinWidth(h->GetXaxis()->FindBin(e)));

//initialize source momemtum direction
G4double dphi=6.28*G4UniformRand();
G4double dr=0.2*G4UniformRand();
G4double x0=sin(dphi)*dr;
G4double y0=cos(dphi)*dr;
G4double z0=(G4UniformRand()*0.6+0.4);//MLCC


//G4double x0=(G4UniformRand()-0.5)*0.4;//MLCC
//G4double y0=(G4UniformRand()-0.5)*0.4;//MLCC
//G4double z0=(G4UniformRand()*0.6+0.4);//MLCC
  //G4double x0=0;
  //G4double y0=0;
  //G4double z0=0;

 px=x0-x;
 py=y0-y;
 pz=z0-z;

 return;
}



void PrimaryGeneratorAction::SetParticle(G4double &x,G4double &y,G4double &z,G4double &px,G4double &py,G4double &pz,G4double &e,TF1* f)
{

//initialize source position in hemisphere, CJPL hole radius = 6.0m

  TF1* fm=new TF1("fm","(8.60*exp(-6.72*(1.0/cos(x))/0.45)+0.44*exp(-6.72*(1.0/cos(x))/0.87))*(1.0/cos(x))",0,1.57);
  G4double u=fm->GetRandom();
  G4double v=6.28*G4UniformRand();
  x=6.0*sin(u)*cos(v);
  y=6.0*sin(u)*sin(v);
  z=6.0*cos(u);
 
  e=1000.0*f->GetRandom();//MeV

//initialize source momemtum direction
G4double dphi=6.28*G4UniformRand();
G4double dr=0.2*G4UniformRand();
G4double x0=sin(dphi)*dr;
G4double y0=cos(dphi)*dr;
G4double z0=(G4UniformRand()*0.6+0.4);//MLCC

//G4double x0=(G4UniformRand()-0.5)*0.4;//MLCC
//G4double y0=(G4UniformRand()-0.5)*0.4;//MLCC
//G4double z0=(G4UniformRand()*0.6+0.4);//MLCC

//G4double x0=0;
//G4double y0=0;
//G4double z0=0;

 px=x0-x;
 py=y0-y;
 pz=z0-z;

 return;
}
*/
