#ifndef EventAction_h
#define EventAction_h 1

#include "globals.hh"

#include "G4UserEventAction.hh"

class RunAction;
/// Event action class
///

class EventAction : public G4UserEventAction
{
public:
  EventAction(RunAction* runAction);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event* event);
  virtual void EndOfEventAction(const G4Event* event);
  void AddEdep(G4double edep,int i) { fEdep[i] += edep; }
private:
  G4int fEventNo;
    RunAction* fRunAction;
    G4double     fEdep[36];
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
