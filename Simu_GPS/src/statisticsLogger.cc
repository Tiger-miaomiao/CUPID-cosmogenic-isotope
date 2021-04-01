#include "g4root.hh"

#include "statisticsLogger.hh"

StatisticsLogger* StatisticsLogger::fInstance = NULL;

StatisticsLogger::StatisticsLogger()
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->SetVerboseLevel(1);

  analysis->OpenFile("result");

  analysis->SetFirstHistoId(0);
  analysis->CreateH1("primaries", "Energy of primaries after attenuation (MeV)", 1000, 0, 10,"MeV"); // histogram 0
  analysis->CreateH1("gammas", "Energy of gammas after attenuation (MeV)", 350, 0, 3.5,"MeV"); // histogram 1
  analysis->CreateH1("electrons", "Energy of electrons after attenuation (MeV)", 350, 0, 3.5,"MeV"); // histogram 2
  analysis->CreateH1("alphas", "Energy of alphas after attenuation (MeV)", 350, 0, 3.5,"MeV"); // histogram 3
  analysis->CreateH1("neutrons", "Energy of neutrons after attenuation (MeV)", 350, 0, 3.5,"MeV"); // histogram 4
  analysis->CreateH1("originals", "Energy of originals after attenuation (MeV)", 1000, 0, 5.0,"MeV"); // histogram 5
  analysis->CreateH1("originals_M1", "Energy of originals M1 after attenuation (MeV)", 1000, 0, 5.0,"MeV"); // histogram 6
  analysis->CreateH3("primarydist", "primarydist", 500,-250,250,500,-250,250,200,0,2000,"cm","cm","cm"); // histogram 0

  analysis->CreateNtuple("statistics", "statistics");
  analysis->SetFirstNtupleColumnId(0);
  analysis->CreateNtupleIColumn("primaryPDGEncoding"); // column 0
  analysis->CreateNtupleIColumn("attenuatorThickness_cm"); // column 1
  analysis->CreateNtupleIColumn("nPrimaries_initial"); // column 2
  analysis->CreateNtupleIColumn("nPrimaries_final"); // column 3
  analysis->CreateNtupleIColumn("nPrimaries_final_aboveROILowerBound"); // column 4
  analysis->CreateNtupleIColumn("nGammas_final"); // column 5
  analysis->CreateNtupleIColumn("nGammas_final_aboveROILowerBound"); // column 6
  analysis->CreateNtupleIColumn("nBetas_final"); // column 7 
  analysis->CreateNtupleIColumn("nBetas_final_aboveROILowerBound"); // column 8
  analysis->CreateNtupleIColumn("nAlphas_final"); // column 9
  analysis->CreateNtupleIColumn("nAlphas_final_aboveROILowerBound"); // column 10
  analysis->CreateNtupleIColumn("nNeutrons_final"); // column 11
  analysis->CreateNtupleIColumn("nNeutrons_final_aboveROILowerBound"); // column 12
  analysis->FinishNtuple();
}

StatisticsLogger::~StatisticsLogger()
{;}

StatisticsLogger* StatisticsLogger::GetInstance()
{
  if (fInstance == NULL) {
    fInstance = new StatisticsLogger();
  }
  return fInstance;
}

void StatisticsLogger::WriteAndClose()
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->Write();
  analysis->CloseFile();
}

void StatisticsLogger::SetPrimaryPDGEncoding(G4int primary)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(0, primary);
}

void StatisticsLogger::SetAttenuatorThickness(G4double thickness)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(1, thickness);
}

void StatisticsLogger::SetInitialNPrimaries(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(2, n);
}

void StatisticsLogger::SetFinalNPrimaries(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(3, n);
}

void StatisticsLogger::SetFinalNPrimariesAboveROILowerBound(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(4, n);
}

void StatisticsLogger::SetFinalNGammas(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(5, n);
}

void StatisticsLogger::SetFinalNGammasAboveROILowerBound(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(6, n);
}

void StatisticsLogger::SetFinalNBetas(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(7, n);
}

void StatisticsLogger::SetFinalNBetasAboveROILowerBound(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(8, n);
}

void StatisticsLogger::SetFinalNAlphas(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(9, n);
}

void StatisticsLogger::SetFinalNAlphasAboveROILowerBound(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(10, n);
}

void StatisticsLogger::SetFinalNNeutrons(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(11, n);
}

void StatisticsLogger::SetFinalNNeutronsAboveROILowerBound(G4int n)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  analysis->FillNtupleIColumn(12, n);
}

G4bool StatisticsLogger::AddNtupleRow()
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  return analysis->AddNtupleRow();
}

void StatisticsLogger::FillPrimaryHistogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(0, energy);
}

void StatisticsLogger::FillGammaHistogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(1, energy);
}

void StatisticsLogger::FillBetaHistogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(2, energy);
}

void StatisticsLogger::FillAlphaHistogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(3, energy);
}

void StatisticsLogger::FillNeutronHistogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(4, energy);
}

void StatisticsLogger::FillOriginalHistogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(5, energy);
}

void StatisticsLogger::FillOriginalM1Histogram(G4double energy)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  analysis->FillH1(6, energy);
}

void StatisticsLogger::FillPrimarydistHistogram(G4double x,G4double y,G4double z)
{
  G4AnalysisManager* analysis = G4AnalysisManager::Instance();
  //std::cout << "analysis->Fill()" << std::endl;
  //G4cout<<"id:x:y:z = "<<analysis->GetH3Id("primarydist")<<" "<<x<<" "<<y<<" "<<z<<G4endl;
  analysis->FillH3(0,x,y,z);
}

