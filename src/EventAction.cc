#include "EventAction.hh"

#include <cstdio>

#include "BgoHit.hh"
#include "DsssdHit.hh"
#include "G4AnalysisManager.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
#include "G4SystemOfUnits.hh"
#include "HitNtuple.hh"

void EventAction::EndOfEventAction(const G4Event* event) {
  G4HCofThisEvent* hce = event->GetHCofThisEvent();
  if (!hce) return;
  G4SDManager* sdManager = G4SDManager::GetSDMpointer();
  auto* analysisManager = G4AnalysisManager::Instance();
  const G4int eventID = event->GetEventID();

  const G4int bgoHcID = sdManager->GetCollectionID("BgoSD/BgoHitsCollection");
  if (bgoHcID >= 0) {
    auto* hits = static_cast<BgoHitsCollection*>(hce->GetHC(bgoHcID));
    if (hits) {
      double totalEdepMeV = 0.0;
      std::printf("HITS BGO   n=%zu crystals=[", hits->entries());
      for (std::size_t i = 0; i < hits->entries(); ++i) {
        const BgoHit* hit = (*hits)[i];
        totalEdepMeV += hit->edepMeV;
        std::printf("%s%d", i == 0 ? "" : ",", hit->crystalID);

        using namespace HitNtuple;
        analysisManager->FillNtupleIColumn(kBgoNtupleID, kBgoEventID, eventID);
        analysisManager->FillNtupleIColumn(kBgoNtupleID, kBgoCrystalID, hit->crystalID);
        analysisManager->FillNtupleDColumn(kBgoNtupleID, kBgoEdepMeV, hit->edepMeV);
        analysisManager->FillNtupleDColumn(kBgoNtupleID, kBgoXCm, hit->position.x() / cm);
        analysisManager->FillNtupleDColumn(kBgoNtupleID, kBgoYCm, hit->position.y() / cm);
        analysisManager->FillNtupleDColumn(kBgoNtupleID, kBgoZCm, hit->position.z() / cm);
        analysisManager->FillNtupleDColumn(kBgoNtupleID, kBgoTimeNs, hit->timeNs);
        analysisManager->FillNtupleIColumn(kBgoNtupleID, kBgoTrackID, hit->trackID);
        analysisManager->FillNtupleSColumn(kBgoNtupleID, kBgoParticle, hit->particleName);
        analysisManager->AddNtupleRow(kBgoNtupleID);
      }
      std::printf("] edep_total=%.6f MeV\n", totalEdepMeV);
    }
  }

  const G4int dsssdHcID = sdManager->GetCollectionID("DsssdSD/DsssdHitsCollection");
  if (dsssdHcID >= 0) {
    auto* hits = static_cast<DsssdHitsCollection*>(hce->GetHC(dsssdHcID));
    if (hits) {
      double totalEdepMeV = 0.0;
      for (std::size_t i = 0; i < hits->entries(); ++i) {
        const DsssdHit* hit = (*hits)[i];
        totalEdepMeV += hit->edepMeV;

        using namespace HitNtuple;
        analysisManager->FillNtupleIColumn(kDsssdNtupleID, kDsssdEventID, eventID);
        analysisManager->FillNtupleDColumn(kDsssdNtupleID, kDsssdEdepMeV, hit->edepMeV);
        analysisManager->FillNtupleDColumn(kDsssdNtupleID, kDsssdXCm, hit->position.x() / cm);
        analysisManager->FillNtupleDColumn(kDsssdNtupleID, kDsssdYCm, hit->position.y() / cm);
        analysisManager->FillNtupleDColumn(kDsssdNtupleID, kDsssdZCm, hit->position.z() / cm);
        analysisManager->FillNtupleDColumn(kDsssdNtupleID, kDsssdTimeNs, hit->timeNs);
        analysisManager->FillNtupleIColumn(kDsssdNtupleID, kDsssdTrackID, hit->trackID);
        analysisManager->FillNtupleSColumn(kDsssdNtupleID, kDsssdParticle, hit->particleName);
        analysisManager->AddNtupleRow(kDsssdNtupleID);
      }
      std::printf("HITS DSSSD n=%zu edep_total=%.6f MeV\n", hits->entries(), totalEdepMeV);
    }
  }
}
