#include "G4Step.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "eventAction.hh"
#include "statisticsLogger.hh"
#include "steppingAction.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolume.hh"
#include "Randomize.hh"
#include <string>

SteppingAction::SteppingAction(EventAction* eventAction)
: G4UserSteppingAction(),
  fEventAction(eventAction),
  fScoringVolume(0)
{
  counter_primary = 0;
  counter_primary_aboveROILowerBound = 0;
  counter_gamma = 0;
  counter_gamma_aboveROILowerBound = 0;
  counter_beta = 0;
  counter_beta_aboveROILowerBound = 0;
  counter_alpha = 0;
  counter_alpha_aboveROILowerBound = 0;
  counter_neutron = 0;
  counter_neutron_aboveROILowerBound = 0;
  ROI_upperBound = 3.06*CLHEP::MeV;//MLCC
  ROI_lowerBound = 3.0*CLHEP::MeV;//MLCC
}

SteppingAction::~SteppingAction()
{}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
  G4Track* track = step->GetTrack();
  const G4StepPoint* preStep = step->GetPreStepPoint();
  assert(preStep);
  const G4StepPoint* postStep = step->GetPostStepPoint();
  assert(postStep);

  //energy deposition in detector


  G4String detector_id;
  char stt[20];

  //void AddEdep(G4double edep,int i) { fEdep[i] += edep; }

  for(int i=1;i<=36;i++)
    {
   sprintf(stt,"detector%d",i);
   detector_id=stt;

   //G4VPhysicalVolume* volume = step->GetPreStepPoint()->GetTouchableHandle()->GetVolume();
  if(preStep->GetPhysicalVolume()->GetName() == detector_id)
    {
G4double edepStep = step->GetTotalEnergyDeposit();

 if(edepStep>0){fEventAction->AddEdep(edepStep,i-1);}

    }
    }
    
  if (track->GetParentID() == 0) { // primary
    if (step->IsLastStepInVolume() && // last step in detector
        preStep->GetPhysicalVolume()->GetName() == "detector") { //MLCC
      const G4ThreeVector& momentumDirection = step->GetTrack()->GetMomentumDirection();
      if (momentumDirection.x() != 0) {
        G4double energy = postStep->GetKineticEnergy();

        counter_primary++;
        if (energy > ROI_lowerBound) counter_primary_aboveROILowerBound++;

        StatisticsLogger* logger = StatisticsLogger::GetInstance();
        logger->FillPrimaryHistogram(energy);
      } else {
        //G4cout << "Primary is going in opposite direction." << G4endl;
      }
    }
  }

  if (track->GetDefinition() == G4Gamma::GammaDefinition()) { // gamma
    if (step->IsLastStepInVolume() && // last step in detector
        preStep->GetPhysicalVolume()->GetName() == "detector") {
      const G4ThreeVector& momentumDirection = step->GetTrack()->GetMomentumDirection();
      if (momentumDirection.x() > 0) {
        G4double energy = postStep->GetKineticEnergy();

        counter_gamma++;
        if (energy > ROI_lowerBound) counter_gamma_aboveROILowerBound++;

        StatisticsLogger* logger = StatisticsLogger::GetInstance();
        logger->FillGammaHistogram(energy);
      } else {
        //G4cout << "Gamma is going in opposite direction." << G4endl;
      }
    }
  } else if (track->GetDefinition() == G4Electron::ElectronDefinition()) { // beta
    if (step->IsLastStepInVolume() && // last step in detector
        preStep->GetPhysicalVolume()->GetName() == "detector") {
      const G4ThreeVector& momentumDirection = step->GetTrack()->GetMomentumDirection();
      if (momentumDirection.x() > 0) {
        G4double energy = postStep->GetKineticEnergy();

        counter_beta++;
        if (energy > ROI_lowerBound) counter_beta_aboveROILowerBound++;

        StatisticsLogger* logger = StatisticsLogger::GetInstance();
        logger->FillBetaHistogram(energy);
      } else {
        //G4cout << "Beta is going in opposite direction." << G4endl;
      }
    }
  } else if (track->GetDefinition() == G4Alpha::AlphaDefinition()) { // alpha
    if (step->IsLastStepInVolume() && // last step in attenuator
        preStep->GetPhysicalVolume()->GetName() == "detector") {
      const G4ThreeVector& momentumDirection = step->GetTrack()->GetMomentumDirection();
      if (momentumDirection.x() > 0) {
        G4double energy = postStep->GetKineticEnergy();

        counter_alpha++;
        if (energy > ROI_lowerBound) counter_alpha_aboveROILowerBound++;

        StatisticsLogger* logger = StatisticsLogger::GetInstance();
        logger->FillAlphaHistogram(energy);
      } else {
        //G4cout << "Alpha is going in opposite direction." << G4endl;
      }
    }
  } else if (track->GetDefinition() == G4Neutron::NeutronDefinition()) { // Neutron
    if (step->IsLastStepInVolume() && // last step in attenuator
        preStep->GetPhysicalVolume()->GetName() == "detector") {
      const G4ThreeVector& momentumDirection = step->GetTrack()->GetMomentumDirection();
      if (momentumDirection.x() > 0) {
        G4double energy = postStep->GetKineticEnergy();

        counter_neutron++;
        if (energy > ROI_lowerBound) counter_neutron_aboveROILowerBound++;

        StatisticsLogger* logger = StatisticsLogger::GetInstance();
        logger->FillNeutronHistogram(energy);
      } else {
        //G4cout << "Neutron is going in opposite direction." << G4endl;
      }
    }
  }
}
