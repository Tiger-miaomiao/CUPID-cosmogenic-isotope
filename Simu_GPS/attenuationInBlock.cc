#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "PhysicsList.hh"
#include "statisticsLogger.hh"
#include "detectorConstruction.hh"
#include "primaryGeneratorAction.hh"
#include "runAction.hh"
#include "steppingAction.hh"
#include "eventAction.hh"

#include "G4PhysListFactory.hh"
#include "G4VUserPhysicsList.hh"

int main(int argc, char** argv)
{

  //if(argc!=2) return 0;

  //G4String source_name;

 //if(argc==2){source_name=argv[1];}

///////setup random seed
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  G4long seed = time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);

  G4RunManager* runManager = new G4RunManager();

  DetectorConstruction* detectorConstruction = new DetectorConstruction();
  runManager->SetUserInitialization(detectorConstruction);

/*
  //runManager->SetUserInitialization(new PhysicsList());
  PhysicsList *aPhysicsList = new PhysicsList();
  runManager->SetUserInitialization(aPhysicsList);
  //aPhysicsList->SelectPhysicsList("Livermore_EM");//MLCC
  //aPhysicsList->SelectPhysicsList("Hadron");
  //aPhysicsList->SelectPhysicsList("QGSP_BERT_HP");
*/

G4PhysListFactory *physListFactory = new G4PhysListFactory(); 
G4VUserPhysicsList *aphysicsList =  physListFactory->GetReferencePhysList("Shielding"); 
runManager->SetUserInitialization(aphysicsList);
 

  runManager->SetUserAction(new PrimaryGeneratorAction());

  RunAction* runAction = new RunAction();
  runManager->SetUserAction(runAction);

  EventAction* eventAction = new EventAction(runAction);
  runManager->SetUserAction(eventAction);

  SteppingAction* steppingAction = new SteppingAction(eventAction);
  runManager->SetUserAction(steppingAction);

  runManager->Initialize();

  //Visualization

  G4VisManager* visManager = new G4VisExecutive();//LM
  visManager->Initialize();//LM

  G4UImanager* UImanager = G4UImanager::GetUIpointer();//LM
  UImanager->ApplyCommand("/run/verbose 1");//LM
  //UImanager->ApplyCommand("/event/verbose 1");//LM
  //UImanager->ApplyCommand("/tracking/verbose 1");//LM
  UImanager->ApplyCommand("/control/execute vis.mac");//LM

  //  const int nEvents = 1000;
  //const int nEvents = 10000;//set No. Events

  //StatisticsLogger* logger = StatisticsLogger::GetInstance();
  //logger->SetInitialNPrimaries(nEvents);
  //runManager->BeamOn(nEvents);

  G4UIExecutive* ui = new G4UIExecutive(argc, argv);//LM
  ui->SessionStart();//LM

  delete ui;//LM
  delete visManager;//LM
  delete runManager;

  return 0;
}
