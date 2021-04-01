#include "G4VUserDetectorConstruction.hh"

class G4VPhysicalVolume;
class G4LogicalVolume;
class G4Region;

class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
  DetectorConstruction();
  virtual ~DetectorConstruction();

G4LogicalVolume* GetScoringVolume() const { return fScoringVolume; }

public:
  virtual G4VPhysicalVolume* Construct();

G4LogicalVolume*  fScoringVolume;
G4LogicalVolume**  logicBox_LMO;

G4Region*          fCopperShield;
G4Region*          fLeadShield;
G4Region*          fSSShield;

};
