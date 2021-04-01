//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file radioactivedecay/rdecay02/include/PhysicsList.hh
/// \brief Definition of the PhysicsList class
//

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef PhysicsList_h
#define PhysicsList_h 1

#include "G4VModularPhysicsList.hh"
#include "PhysListHadron.hh"
#include "G4HadronPhysicsQGSP_BERT_HP.hh"
#include "globals.hh"
#include <vector>

class G4VPhysicsConstructor;
class PhysicsListMessenger;
class G4ProductionCuts;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class PhysicsList: public G4VModularPhysicsList
{
public:
  PhysicsList(bool fActiveNeutrons, bool fActiveOptical, G4double fIonEner, G4double fIonStep)
	   : G4VModularPhysicsList(), aIonEner(fIonEner), aIonStep(fIonStep)
  {
  	PhysListInitialize(fActiveNeutrons,fActiveOptical);
  }
  
  PhysicsList(bool fActiveNeutrons, bool fActiveOptical)
	   : G4VModularPhysicsList(), aIonEner(0.), aIonStep(0.)
  {
  	PhysListInitialize(fActiveNeutrons,fActiveOptical);
  }
  
  PhysicsList() : G4VModularPhysicsList(), aIonEner(0.), aIonStep(0.)
  {
  	PhysListInitialize( true, true);
  }
  
  PhysicsList(G4double fIonEner, G4double fIonStep) : G4VModularPhysicsList(), aIonEner(fIonEner), aIonStep(fIonStep)
  {
  	PhysListInitialize( true, true);
  }

  virtual ~PhysicsList();

  void PhysListInitialize(bool fActiveNeutrons, bool fActiveOptical);
  virtual void ConstructParticle();

  virtual void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);
  void SetCutForProton(G4double);
  void SetCutForAlpha(G4double);
  void SetCutForGenericIon(G4double);

  void SelectPhysicsList(const G4String& name);
  virtual void ConstructProcess();

  void SetTargetCut(G4double val);
  void SetDetectorCut(G4double val);

  void SetCutScale(G4double val);
  void killRadioactiveDecays() { fradDecay = false; delete fRaddecayList; }
  void ConstructOp();

private:

  void AddExtraBuilders(G4bool flagHP);

  // hide assignment operator
  PhysicsList & operator=(const PhysicsList &right);
  PhysicsList(const PhysicsList&);

  bool  ActiveNeutrons;
  bool  ActiveOptical;
  G4double fCutForGamma;
  G4double fCutForElectron;
  G4double fCutForPositron;
  G4double fCutForProton;
  G4double fCutForAlpha;
  G4double fCutForGenericIon;
  G4double fCutForOpticalPhoton;
  G4double fSF;
  G4double aIonEner, aIonStep;
  bool fradDecay;

  G4VPhysicsConstructor*  fEmPhysicsList;
  G4VPhysicsConstructor*  fRaddecayList;
  G4VPhysicsConstructor*  fParticleList;
  G4VPhysicsConstructor*  fHadPhysicsList;
  PhysListHadron*  aHadPhysicsList;
  G4HadronPhysicsQGSP_BERT_HP *aBertPhysicsList;

  std::vector<G4VPhysicsConstructor*>  fHadronPhys;
  G4int fNhadcomp;  

  PhysicsListMessenger* fPMessenger;
  G4ProductionCuts* fDetectorCuts;
  G4ProductionCuts* fTargetCuts;

protected:
  virtual void AddTransportation();
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
