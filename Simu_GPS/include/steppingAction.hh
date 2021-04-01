#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include "eventAction.hh"

class EventAction;
class G4LogicalVolume;

class SteppingAction : public G4UserSteppingAction {
public:
  SteppingAction(EventAction* eventAction);
  virtual ~SteppingAction();

public:
  virtual void UserSteppingAction(const G4Step*);
  inline G4int GetPrimaryCounter() const { return counter_primary; };
  inline G4int GetPrimaryAboveROILowerBoundCounter() const { return counter_primary_aboveROILowerBound; };
  inline G4int GetGammaCounter() const { return counter_gamma; };
  inline G4int GetGammaAboveROILowerBoundCounter() const { return counter_gamma_aboveROILowerBound; };
  inline G4int GetBetaCounter() const { return counter_beta; };
  inline G4int GetBetaAboveROILowerBoundCounter() const { return counter_beta_aboveROILowerBound; };
  inline G4int GetAlphaCounter() const { return counter_alpha; };
  inline G4int GetAlphaAboveROILowerBoundCounter() const { return counter_alpha_aboveROILowerBound; };
  inline G4int GetNeutronCounter() const { return counter_neutron; };
  inline G4int GetNeutronAboveROILowerBoundCounter() const { return counter_neutron_aboveROILowerBound; };
  inline G4double GetROILowerBound() const { return ROI_lowerBound; };

private:

  EventAction*  fEventAction;
  G4LogicalVolume* fScoringVolume;
  int counter_primary;
  int counter_primary_aboveROILowerBound;
  int counter_gamma;
  int counter_gamma_aboveROILowerBound;
  int counter_beta;
  int counter_beta_aboveROILowerBound;
  int counter_alpha;
  int counter_alpha_aboveROILowerBound;
  int counter_neutron;
  int counter_neutron_aboveROILowerBound;
  G4double ROI_lowerBound;
  G4double ROI_upperBound;
};
