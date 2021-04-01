#include "G4UserSteppingAction.hh"
#include "G4ios.hh"

class SteppingAction : public G4UserSteppingAction
{
public:
	SteppingAction();
	virtual ~SteppingAction();

public:
	virtual void UserSteppingAction(const G4Step*);
	G4int GetCounter();

private:
	int counter;
};

inline G4int SteppingAction::GetCounter()
{
	return counter;
}
