#include "G4RunManager.hh"

#include "statisticsLogger.hh"
#include "steppingAction.hh"
#include "runAction.hh"

#include "G4AccumulableManager.hh"

RunAction::RunAction()
: G4UserRunAction(),
  fEdep(0.)
{
G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->RegisterAccumulable(fEdep); 
}

RunAction::~RunAction()
{}

void RunAction::BeginOfRunAction(const G4Run*)
{
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Reset();
}

void RunAction::EndOfRunAction(const G4Run*)
{

  // Merge accumulables 
  G4AccumulableManager* accumulableManager = G4AccumulableManager::Instance();
  accumulableManager->Merge();

 G4double edep  = fEdep.GetValue();

 //G4cout<<"runaction energy depo:"<<edep<<G4endl;

  StatisticsLogger* logger = StatisticsLogger::GetInstance();
  G4RunManager* runManager = G4RunManager::GetRunManager();
  const SteppingAction* steppingAction = dynamic_cast<const SteppingAction*>(runManager->GetUserSteppingAction());
  logger->SetFinalNPrimaries(steppingAction->GetPrimaryCounter());
  logger->SetFinalNPrimariesAboveROILowerBound(steppingAction->GetPrimaryAboveROILowerBoundCounter());
  logger->SetFinalNGammas(steppingAction->GetGammaCounter());
  logger->SetFinalNGammasAboveROILowerBound(steppingAction->GetGammaAboveROILowerBoundCounter());
  logger->SetFinalNBetas(steppingAction->GetBetaCounter());
  logger->SetFinalNBetasAboveROILowerBound(steppingAction->GetBetaAboveROILowerBoundCounter());
  logger->SetFinalNAlphas(steppingAction->GetAlphaCounter());
  logger->SetFinalNAlphasAboveROILowerBound(steppingAction->GetAlphaAboveROILowerBoundCounter());
  logger->SetFinalNNeutrons(steppingAction->GetNeutronCounter());
  logger->SetFinalNNeutronsAboveROILowerBound(steppingAction->GetNeutronAboveROILowerBoundCounter());
  logger->AddNtupleRow();
  logger->WriteAndClose();
}

void RunAction::AddEdep(G4double edep)
{
  fEdep  += edep;
}
