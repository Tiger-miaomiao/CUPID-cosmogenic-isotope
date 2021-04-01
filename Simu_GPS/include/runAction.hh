#include "G4UserRunAction.hh"
#include "G4Accumulable.hh"
#include "globals.hh"

class G4Run;

class RunAction : public G4UserRunAction {
public:
  RunAction();
  virtual ~RunAction();

public:
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);
  void AddEdep (G4double edep);

private:
    G4Accumulable<G4double> fEdep;
};
