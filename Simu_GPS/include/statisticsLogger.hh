#include "globals.hh"

class StatisticsLogger {
public:
  static StatisticsLogger* GetInstance();

public:
  void SetPrimaryPDGEncoding(G4int primary);
  void SetAttenuatorThickness(G4double thickness);
  void SetInitialNPrimaries(G4int n);
  void SetFinalNPrimaries(G4int n);
  void SetFinalNPrimariesAboveROILowerBound(G4int n);
  void SetFinalNGammas(G4int n);
  void SetFinalNGammasAboveROILowerBound(G4int n);
  void SetFinalNBetas(G4int n);
  void SetFinalNBetasAboveROILowerBound(G4int n);
  void SetFinalNAlphas(G4int n);
  void SetFinalNAlphasAboveROILowerBound(G4int n);
  void SetFinalNNeutrons(G4int n);
  void SetFinalNNeutronsAboveROILowerBound(G4int n);
  G4bool AddNtupleRow();
  void FillPrimaryHistogram(G4double energy);
  void FillGammaHistogram(G4double energy);
  void FillBetaHistogram(G4double energy);
  void FillAlphaHistogram(G4double energy);
  void FillNeutronHistogram(G4double energy);
  void FillOriginalHistogram(G4double energy);
  void FillOriginalM1Histogram(G4double energy);
  void FillPrimarydistHistogram(G4double x,G4double y,G4double z);
  void WriteAndClose();

private:
  StatisticsLogger();
  virtual ~StatisticsLogger();

private:
  static StatisticsLogger* fInstance;
};
