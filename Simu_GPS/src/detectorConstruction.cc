#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4ThreeVector.hh"
#include "G4VSolid.hh"
#include "G4UnionSolid.hh"
#include "G4SubtractionSolid.hh"

#include "statisticsLogger.hh"
#include "detectorConstruction.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4Region.hh"
#include "G4RegionStore.hh"
#include "G4ProductionCuts.hh"

DetectorConstruction::DetectorConstruction()
  : fScoringVolume(0),logicBox_LMO(NULL)
{
  logicBox_LMO = new G4LogicalVolume*[36];
}

DetectorConstruction::~DetectorConstruction()
{
  delete [] logicBox_LMO;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  G4String name; // name of material
  G4double density; // density of material
  G4int nComponents, nAtoms; // number of composing atom types, number of atoms composing molecule of material

  G4NistManager* nist = G4NistManager::Instance();

  G4Material* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
  G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
  G4Material* Cu = nist->FindOrBuildMaterial("G4_Cu");
  G4Material* Pb = nist->FindOrBuildMaterial("G4_Pb");
  G4Material* Al = nist->FindOrBuildMaterial("G4_Al");
  G4Material* SS =nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
  G4Element* Te = nist->FindOrBuildElement("Te");
  G4Element* O = nist->FindOrBuildElement("O");
  G4Element* H = nist->FindOrBuildElement("H");
  G4Element* C = nist->FindOrBuildElement("C");
  G4Element* B = nist->FindOrBuildElement("B");
  G4Element* Li = nist->FindOrBuildElement("Li");
  G4Element* Mo = nist->FindOrBuildElement("Mo");

  // Li2MoO4
  G4Material* LMO = new G4Material(name = "LMO", density = 3.07*CLHEP::g/CLHEP::cm3, nComponents = 3);
  LMO->AddElement(Mo, nAtoms = 1);
  LMO->AddElement(O, nAtoms = 4);
  LMO->AddElement(Li, nAtoms = 2);

  // C2H4B
  G4Material* CHB = new G4Material(name = "CHB", density= 1.09*CLHEP::g/CLHEP::cm3, nComponents=3);
 CHB->AddElement(B, 0.05);
 CHB->AddElement(H, 0.136);
 CHB->AddElement(C, 0.814);

  G4double worldX = 6.0*CLHEP::m;
  G4double worldY = 6.0*CLHEP::m;
  G4double worldZ = 6.0*CLHEP::m;
  G4Box* solidWorld = new G4Box("World",worldX,worldY,worldZ);

  G4LogicalVolume* logicWorld = 
    new G4LogicalVolume(solidWorld,vacuum,"World");

  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,G4ThreeVector(),logicWorld,"World",0,false,0,true);

  //detector 
  const G4double thickness = 4.5*CLHEP::cm;//晶体尺寸
  StatisticsLogger* logger = StatisticsLogger::GetInstance();
  logger->SetAttenuatorThickness(thickness/CLHEP::cm);

  G4double boxX1 = 2.25*CLHEP::cm;
  G4double boxY1 = 2.25*CLHEP::cm;
  G4double boxZ1 = 2.25*CLHEP::cm;

  G4ThreeVector GVDetector;
  G4String detector;
  char st[20];


  //detector LMO array

for(int i=1;i<=36;i++)
{
  sprintf(st,"detector%d",i);
  detector=st;
  if(i>=1&&i<=9){GVDetector = G4ThreeVector(3.0*CLHEP::cm,3.0*CLHEP::cm,(110.0-i*6.0)*CLHEP::cm);}
  if(i>=10&&i<=18){GVDetector = G4ThreeVector(-3.0*CLHEP::cm,3.0*CLHEP::cm,(110.0-(i-9)*6.0)*CLHEP::cm);}
  if(i>=19&&i<=27){GVDetector = G4ThreeVector(3.0*CLHEP::cm,-3.0*CLHEP::cm,(110.0-(i-18)*6.0)*CLHEP::cm);}
  if(i>=28&&i<=36){GVDetector = G4ThreeVector(-3.0*CLHEP::cm,-3.0*CLHEP::cm,(110.0-(i-27)*6.0)*CLHEP::cm);}
  G4Box* solidBox = new G4Box(detector,boxX1,boxY1,boxZ1);
  logicBox_LMO[i-1] = new G4LogicalVolume(solidBox,LMO,detector);
  G4VPhysicalVolume* physBox_LMO = new G4PVPlacement(0,GVDetector,logicBox_LMO[i-1],detector,logicWorld,false,0,true);
}

 
//Frame (Cu)

G4double RminTPlate=0*CLHEP::cm, RmaxTPlate=20*CLHEP::cm, deltaZTPlate= 0.5*CLHEP::cm, PhiminTPlate=0., deltaPhiTPlate=360*CLHEP::degree;
G4double RminTScrew=0*CLHEP::cm, RmaxTScrew=0.25*CLHEP::cm, deltaZTScrew= 2.75*CLHEP::cm, PhiminTScrew=0., deltaPhiTScrew=360*CLHEP::degree;
G4VSolid* box_frame = new G4Box("box_frame",6*CLHEP::cm,6*CLHEP::cm,0.25*CLHEP::cm);
G4VSolid* box_sub = new G4Box("box_sub",2.6*CLHEP::cm,2.6*CLHEP::cm,0.25*CLHEP::cm);

G4VSolid* solidPlate = new G4Tubs("Plate",RminTPlate,RmaxTPlate,deltaZTPlate,PhiminTPlate,deltaPhiTPlate);
G4VSolid* solidScrew = new G4Tubs("Screw",RminTScrew,RmaxTScrew,deltaZTScrew,PhiminTScrew,deltaPhiTScrew);

G4VSolid* box_frame_sub = new G4SubtractionSolid("box_frame_sub", box_frame, box_sub,0, G4ThreeVector(2.9*CLHEP::cm,2.9*CLHEP::cm,0.));
box_frame_sub = new G4SubtractionSolid("box_frame_sub", box_frame_sub, box_sub,0, G4ThreeVector(2.9*CLHEP::cm,-2.9*CLHEP::cm,0.));
box_frame_sub = new G4SubtractionSolid("box_frame_sub", box_frame_sub, box_sub,0, G4ThreeVector(-2.9*CLHEP::cm,2.9*CLHEP::cm,0.));
box_frame_sub = new G4SubtractionSolid("box_frame_sub", box_frame_sub, box_sub,0, G4ThreeVector(-2.9*CLHEP::cm,-2.9*CLHEP::cm,0.));

box_frame_sub = new G4UnionSolid("box_frame_sub", box_frame_sub, solidScrew,0,G4ThreeVector(5.75*CLHEP::cm,5.75*CLHEP::cm,-3.0*CLHEP::cm));
box_frame_sub = new G4UnionSolid("box_frame_sub", box_frame_sub, solidScrew,0,G4ThreeVector(-5.75*CLHEP::cm,5.75*CLHEP::cm,-3.0*CLHEP::cm));
box_frame_sub = new G4UnionSolid("box_frame_sub", box_frame_sub, solidScrew,0,G4ThreeVector(5.75*CLHEP::cm,-5.75*CLHEP::cm,-3.0*CLHEP::cm));
box_frame_sub = new G4UnionSolid("box_frame_sub", box_frame_sub, solidScrew,0,G4ThreeVector(-5.75*CLHEP::cm,-5.75*CLHEP::cm,-3.0*CLHEP::cm));
box_frame_sub = new G4UnionSolid("box_frame_sub", box_frame_sub, solidScrew,0,G4ThreeVector(0,0,-3.0*CLHEP::cm));

 for(int i=0;i<=8;i++)
   {
     solidPlate = new G4UnionSolid("Frame", solidPlate, box_frame_sub,0,G4ThreeVector(0,0,(-1.0-6.0*i)*CLHEP::cm));
   }

G4LogicalVolume* logicFrame = new G4LogicalVolume(solidPlate,Cu,"Frame");
G4VPhysicalVolume* physBox_Frame = new G4PVPlacement(0,G4ThreeVector(0,0,108*CLHEP::cm),logicFrame,"Frame",logicWorld,false,0,true);





  //top shield (Al)
  G4ThreeVector GVAlTop = G4ThreeVector(0,0,170.0*CLHEP::cm);
  G4double boxXAlTop = 50.0*CLHEP::cm;
  G4double boxYAlTop = 50.0*CLHEP::cm;
  G4double boxZAlTop = 7.5*CLHEP::cm;
  G4Box* solidBoxAlTop = new G4Box("BoxAlTop",boxXAlTop,boxYAlTop,boxZAlTop);
  G4LogicalVolume* logicBox_AlTop = new G4LogicalVolume(solidBoxAlTop,Al,"BoxAlTop");
  G4VPhysicalVolume* physBox_AlTop = new G4PVPlacement(0,GVAlTop,logicBox_AlTop,"BoxAlTop",logicWorld,false,0,true);


//top ss plate

G4ThreeVector GVSSTopPlate = G4ThreeVector(0,0,150.0*CLHEP::cm);
G4double RminSSTopPlate=0*CLHEP::cm, RmaxSSTopPlate=40.0*CLHEP::cm, deltaZSSTopPlate= 1.0*CLHEP::cm, PhiminSSTopPlate=0., deltaPhiSSTopPlate=360*CLHEP::degree;
G4Tubs* solidSSTopPlate = new G4Tubs("SSTopPlate",RminSSTopPlate,RmaxSSTopPlate,deltaZSSTopPlate,PhiminSSTopPlate,deltaPhiSSTopPlate);
G4LogicalVolume* logicSSTopPlate = new G4LogicalVolume(solidSSTopPlate,SS,"SSTopPlate");  
G4VPhysicalVolume* physSSTopPlate = new G4PVPlacement(0,GVSSTopPlate,logicSSTopPlate,"SSTopPlate",logicWorld,false,0); 


//top copper plate (Cu)  

 G4ThreeVector GVTCuTopPlate = G4ThreeVector(0,0,140.0*CLHEP::cm);
 G4double RminTCuTopPlate=0*CLHEP::cm, RmaxTCuTopPlate=35.0*CLHEP::cm, deltaZTCuTopPlate= 1.0*CLHEP::cm, PhiminTCuTopPlate=0., deltaPhiTCuTopPlate=360*CLHEP::degree;
 G4Tubs* solidTCuTopPlate = new G4Tubs("CuTopPlate",RminTCuTopPlate,RmaxTCuTopPlate,deltaZTCuTopPlate,PhiminTCuTopPlate,deltaPhiTCuTopPlate);   
 G4LogicalVolume* logicTCuTopPlate = new G4LogicalVolume(solidTCuTopPlate,Cu,"CuTopPlate");
 G4VPhysicalVolume* physTCuTopPlate = new G4PVPlacement(0,GVTCuTopPlate,logicTCuTopPlate,"CuTopPlate",logicWorld,false,0);


//top copper shield (Cu)
  G4ThreeVector GVCuTopShield = G4ThreeVector(0,0,120.0*CLHEP::cm);
 G4Tubs* solidTCuTopShield = new G4Tubs("CuTopShield",0*CLHEP::cm,25.0*CLHEP::cm,6.0*CLHEP::cm,0,360*CLHEP::degree);   
  G4LogicalVolume* logicTCuTopShield = new G4LogicalVolume(solidTCuTopShield,Cu,"CuTopShield");
  G4VPhysicalVolume* physTCuTopShield = new G4PVPlacement(0,GVCuTopShield,logicTCuTopShield,"CuTopShield",logicWorld,false,0,true);

  

//thermal shield1 (Cu)
 G4ThreeVector GVTub1 = G4ThreeVector(0,0,65.0*CLHEP::cm);
 G4double RminTub1=25.0*CLHEP::cm, RmaxTub1=25.2*CLHEP::cm, deltaZTub1= 35.0*CLHEP::cm, PhiminTub1=0., deltaPhiTub1=360*CLHEP::degree;
 G4Tubs* solidTub1 = new G4Tubs("ThermalTub1",RminTub1,RmaxTub1,deltaZTub1,PhiminTub1,deltaPhiTub1);
 G4LogicalVolume* logicTub1 = new G4LogicalVolume(solidTub1,Cu,"ThermalTub1");  
 G4VPhysicalVolume* physTub1 = new G4PVPlacement(0,GVTub1,logicTub1,"ThermalTub1",logicWorld,false,0); 
  


//thermal shield2 (Cu)
 G4ThreeVector GVTub2 = G4ThreeVector(0,0,65.0*CLHEP::cm);
 G4double RminTub2=27.0*CLHEP::cm, RmaxTub2=27.2*CLHEP::cm, deltaZTub2= 35.0*CLHEP::cm, PhiminTub2=0., deltaPhiTub2=360*CLHEP::degree;
 G4Tubs* solidTub2 = new G4Tubs("ThermalTub2",RminTub2,RmaxTub2,deltaZTub2,PhiminTub2,deltaPhiTub2);
 G4LogicalVolume* logicTub2 = new G4LogicalVolume(solidTub2,Cu,"ThermalTub2");  
 G4VPhysicalVolume* physTub2 = new G4PVPlacement(0,GVTub2,logicTub2,"ThermalTub2",logicWorld,false,0); 


		
//bottom shield (Cu)

  G4ThreeVector GVCuBottom = G4ThreeVector(0,0,6.0*CLHEP::cm);
  G4double boxXCuBottom = 25.0*CLHEP::cm;
  G4double boxYCuBottom = 25.0*CLHEP::cm;
  G4double boxZCuBottom = 6.0*CLHEP::cm;
  //G4ThreeVector position4 = G4ThreeVector(3.0*CLHEP::cm,0,0);
  G4Box* solidBoxCuBottom = new G4Box("CuBottomShield",boxXCuBottom,boxYCuBottom,boxZCuBottom);
  G4LogicalVolume* logicTCuBottomShield = new G4LogicalVolume(solidBoxCuBottom,Cu,"CuBottomShield");
  G4VPhysicalVolume* physTCuBottom = new G4PVPlacement(0,GVCuBottom,logicTCuBottomShield,"CuBottomShield",logicWorld,false,0,true);





  ////////////////////////SetCuts threshold

  G4ProductionCuts* fCopperShieldCuts   = new G4ProductionCuts();
  fCopperShieldCuts->SetProductionCut(0.5*CLHEP::mm,"gamma");
  fCopperShieldCuts->SetProductionCut(0.5*CLHEP::mm,"e-");
  fCopperShieldCuts->SetProductionCut(0.5*CLHEP::mm,"e+");

  G4ProductionCuts* fSSShieldCuts   = new G4ProductionCuts();
  fSSShieldCuts->SetProductionCut(0.5*CLHEP::mm,"gamma");
  fSSShieldCuts->SetProductionCut(0.5*CLHEP::mm,"e-");
  fSSShieldCuts->SetProductionCut(0.5*CLHEP::mm,"e+");

  //G4ProductionCuts* fLeadShieldCuts   = new G4ProductionCuts();
  //fLeadShieldCuts->SetProductionCut(50*CLHEP::mm,"gamma");
  //fLeadShieldCuts->SetProductionCut(50*CLHEP::mm,"e-");
  //fLeadShieldCuts->SetProductionCut(50*CLHEP::mm,"e+");

  fSSShield = new G4Region("SSShield");
  fSSShield->SetProductionCuts(fSSShieldCuts);
  fSSShield->AddRootLogicalVolume(logicSSTopPlate);

  fCopperShield = new G4Region("CopperShield");
  fCopperShield->SetProductionCuts(fCopperShieldCuts);
  fCopperShield->AddRootLogicalVolume(logicTCuTopPlate);
  fCopperShield->AddRootLogicalVolume(logicTCuTopShield);
  fCopperShield->AddRootLogicalVolume(logicTub1);
  fCopperShield->AddRootLogicalVolume(logicTub2);
  fCopperShield->AddRootLogicalVolume(logicTCuBottomShield);
  fCopperShield->AddRootLogicalVolume(logicFrame);

//fLeadShield = new G4Region("LeadShield");
//fLeadShield->SetProductionCuts(fLeadShieldCuts);
//fLeadShield->AddRootLogicalVolume(logicBox_PbL1);
//fLeadShield->AddRootLogicalVolume(logicBox_PbL2);
//fLeadShield->AddRootLogicalVolume(logicBox_PbL3);
//fLeadShield->AddRootLogicalVolume(logicBox_PbL4);

  
//Set Visualization 
//(r,g,b,a) a: Alpha is opacity. Its default value 1 means "opaque".
logicWorld->SetVisAttributes(G4Colour(1.0,1.0,1.0,1.0));
//logicBox_detector->SetVisAttributes(G4Colour(0.5,0.5,0.5,0.5)); 
logicTCuBottomShield->SetVisAttributes(G4Colour(1.0,1.0,0,0.5));
logicTub1->SetVisAttributes(G4Colour(0,1.0,0,0.5));
logicTub2->SetVisAttributes(G4Colour(0,1.0,0,0.5));
logicTCuTopShield->SetVisAttributes(G4Colour(1.0,1.0,0,0.5));

logicBox_AlTop->SetVisAttributes(G4Colour(1.0,1.0,1.0,0.5));
//fScoringVolume = logicBox_detector;

  return physWorld;
}




