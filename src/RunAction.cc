#include "RunAction.hh"

#include "G4AnalysisManager.hh"
#include "HitNtuple.hh"

RunAction::RunAction() {
  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->SetDefaultFileType("root");
  analysisManager->SetVerboseLevel(0);

  // "Bgo": one row per BGO hit (BgoHit, see BgoHit.hh) -- column order
  // here must match HitNtuple.hh's kBgo* indices exactly.
  analysisManager->CreateNtuple("Bgo", "BGO crystal (SCNT) hits");
  analysisManager->CreateNtupleIColumn("eventID");   // kBgoEventID
  analysisManager->CreateNtupleIColumn("crystalID");  // kBgoCrystalID, 1-30
  analysisManager->CreateNtupleDColumn("edepMeV");    // kBgoEdepMeV
  analysisManager->CreateNtupleDColumn("x_cm");        // kBgoXCm
  analysisManager->CreateNtupleDColumn("y_cm");        // kBgoYCm
  analysisManager->CreateNtupleDColumn("z_cm");        // kBgoZCm
  analysisManager->CreateNtupleDColumn("t_ns");        // kBgoTimeNs
  analysisManager->CreateNtupleIColumn("trackID");    // kBgoTrackID
  analysisManager->CreateNtupleSColumn("particle");   // kBgoParticle
  analysisManager->FinishNtuple();  // -> HitNtuple::kBgoNtupleID (0)

  // "Dsssd": one row per DSSSD hit (DsssdHit, see DsssdHit.hh) -- no
  // per-strip channel yet (see README), just the whole slab's own hits.
  analysisManager->CreateNtuple("Dsssd", "DSSSD hits");
  analysisManager->CreateNtupleIColumn("eventID");  // kDsssdEventID
  analysisManager->CreateNtupleDColumn("edepMeV");   // kDsssdEdepMeV
  analysisManager->CreateNtupleDColumn("x_cm");       // kDsssdXCm
  analysisManager->CreateNtupleDColumn("y_cm");       // kDsssdYCm
  analysisManager->CreateNtupleDColumn("z_cm");       // kDsssdZCm
  analysisManager->CreateNtupleDColumn("t_ns");       // kDsssdTimeNs
  analysisManager->CreateNtupleIColumn("trackID");   // kDsssdTrackID
  analysisManager->CreateNtupleSColumn("particle");  // kDsssdParticle
  analysisManager->FinishNtuple();  // -> HitNtuple::kDsssdNtupleID (1)
}

void RunAction::BeginOfRunAction(const G4Run*) {
  G4AnalysisManager::Instance()->OpenFile("dragon_hits.root");
}

void RunAction::EndOfRunAction(const G4Run*) {
  auto* analysisManager = G4AnalysisManager::Instance();
  analysisManager->Write();
  analysisManager->CloseFile();
}
