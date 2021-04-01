#include "Randomize.hh"

#include "eventAction.hh"
#include "runAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"

#include "statisticsLogger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(RunAction* runAction)
: G4UserEventAction(),
  fRunAction(runAction)
{
  fEventNo = 0;
} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{    
  //G4cout << "Event number: " << fEventNo << G4endl;
  G4Random::setTheSeed(fEventNo);

for(int i=0;i<=35;i++){fEdep[i] = 0;}
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event*)
{   
  fEventNo++;
  //fRunAction->AddEdep(fEdep);

StatisticsLogger* logger = StatisticsLogger::GetInstance();

 G4int count_M1=0;
 G4double rec_edep=0;

 for(int i=0;i<=35;i++)
   {
 if(fEdep[i]>0)
   {

     //if(fEdep[i]>0){count_M1++;}

     rec_edep=G4RandGauss::shoot(fEdep[i],(fEdep[i]*0.001+2.25*0.001)/2.25);
count_M1++;
logger->FillOriginalHistogram(rec_edep);

//G4cout<<"crystal:"<<i+1<<" edep:"<<fEdep[i]<<" rec_edep:"<<rec_edep<<G4endl;
   }
   }

 if(count_M1==1&&rec_edep>0)
{
logger->FillOriginalM1Histogram(rec_edep);
}

 if(fEventNo%500==0){G4cout<<"fEventNo="<<fEventNo<<G4endl;}
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
