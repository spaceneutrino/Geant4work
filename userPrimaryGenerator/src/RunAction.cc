#include "RunAction.hh"
#include "G4AnalysisManager.hh"
#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::RunAction()
{
  // set printing event number per each event
  G4RunManager::GetRunManager()->SetPrintProgress(1);

  // Create analysis manager
  // The choice of the output format is done via the specified
  // file extension.
  auto analysisManager = G4AnalysisManager::Instance();

  // See Event Action
  // Create directories
  //analysisManager->SetHistoDirectoryName("histograms");
  //analysisManager->SetNtupleDirectoryName("ntuple");
  analysisManager->SetVerboseLevel(1);
  analysisManager->SetNtupleMerging(true);
    // Note: merging ntuples is available only with Root output

  // Book histograms, ntuple
  //

  // Creating histograms
  // analysisManager->CreateH1("N","Number", 10, 0., 20);
  

  // Creating ntuple
  //
 // analysisManager->CreateNtuple("B4", "Edep and TrackL");
 // analysisManager->CreateNtupleDColumn("Eabs");
 // analysisManager->CreateNtupleDColumn("Egap");
 // analysisManager->CreateNtupleDColumn("Labs");
 // analysisManager->CreateNtupleDColumn("Lgap");
 // analysisManager->FinishNtuple();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

RunAction::~RunAction()
{
   G4RunManager::GetRunManager()->SetPrintProgress(1000);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::BeginOfRunAction(const G4Run* /*run*/)
{
  //inform the runManager to save random number seed
  //G4RunManager::GetRunManager()->SetRandomNumberStore(true);

  // Get analysis manager
  auto analysisManager = G4AnalysisManager::Instance();
  

  // Open an output file
  //
  //G4String fileName = "B3Test.root";
  // Other supported output types:
  // G4String fileName = "B4.csv";
  // G4String fileName = "B4.hdf5";
//  // G4String fileName = "B4.xml";
//  analysisManager->OpenFile(fileName);
//  G4cout << "Using " << analysisManager->GetType() << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void RunAction::EndOfRunAction(const G4Run* /*run*/)
{}
  // print histogram statistics
  //
 
 // analysisManager->Write();
 // analysisManager->CloseFile();


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

