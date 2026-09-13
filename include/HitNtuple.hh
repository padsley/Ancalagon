// Shared ntuple/column layout for the "Bgo" and "Dsssd" ROOT TTrees --
// RunAction.cc creates columns in exactly this order (G4AnalysisManager
// assigns ntuple/column IDs by creation order, not by name), EventAction.cc
// fills them; kept in one place so the two can't drift out of sync.
#pragma once

namespace HitNtuple {

constexpr int kBgoNtupleID = 0;
constexpr int kDsssdNtupleID = 1;

// Bgo columns.
constexpr int kBgoEventID = 0;
constexpr int kBgoCrystalID = 1;
constexpr int kBgoEdepMeV = 2;
constexpr int kBgoXCm = 3;
constexpr int kBgoYCm = 4;
constexpr int kBgoZCm = 5;
constexpr int kBgoTimeNs = 6;
constexpr int kBgoTrackID = 7;
constexpr int kBgoParticle = 8;

// Dsssd columns.
constexpr int kDsssdEventID = 0;
constexpr int kDsssdEdepMeV = 1;
constexpr int kDsssdXCm = 2;
constexpr int kDsssdYCm = 3;
constexpr int kDsssdZCm = 4;
constexpr int kDsssdTimeNs = 5;
constexpr int kDsssdTrackID = 6;
constexpr int kDsssdParticle = 7;

}  // namespace HitNtuple
