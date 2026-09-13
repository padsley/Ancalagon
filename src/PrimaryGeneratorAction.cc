#include "PrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "ReactionConfig.hh"

PrimaryGeneratorAction::PrimaryGeneratorAction(double x0Cm, double y0Cm, double momentumMeV,
                                                double z0Cm) {
  fGun = new G4ParticleGun(1);
  G4ParticleDefinition* proton = G4ParticleTable::GetParticleTable()->FindParticle("proton");
  fGun->SetParticleDefinition(proton);
  fGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  fGun->SetParticleMomentum(momentumMeV * MeV);
  fGun->SetParticlePosition(G4ThreeVector(x0Cm * cm, y0Cm * cm, z0Cm * cm));
}

PrimaryGeneratorAction::PrimaryGeneratorAction(double x0Cm, double y0Cm, double momentumMeV,
                                                double z0Cm, int ionZ, int ionA,
                                                int ionChargeState)
    : fNeedsIonLookup(true), fIonZ(ionZ), fIonA(ionA), fIonChargeState(ionChargeState) {
  fGun = new G4ParticleGun(1);
  fGun->SetParticleMomentumDirection(G4ThreeVector(0.0, 0.0, 1.0));
  fGun->SetParticleMomentum(momentumMeV * MeV);
  fGun->SetParticlePosition(G4ThreeVector(x0Cm * cm, y0Cm * cm, z0Cm * cm));
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const std::string& reactionFilePath, double x0Cm,
                                               double y0Cm, double z0Cm)
    : fIsReaction(true), fReactionX0Cm(x0Cm), fReactionY0Cm(y0Cm), fReactionZ0Cm(z0Cm) {
  const ReactionConfig config = ReactionConfig::Load(reactionFilePath);
  fReaction = new ReactionKinematics(config);
  fRecoilZ = config.recoil.Z;
  fRecoilA = config.recoil.A;
  fRecoilChargeState = config.recoilChargeState;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
  delete fReaction;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* event) {
  if (fIsReaction) {
    GenerateReactionPrimaries(event);
    return;
  }
  if (fNeedsIonLookup) {
    G4ParticleDefinition* ion = G4IonTable::GetIonTable()->GetIon(fIonZ, fIonA, 0.0);
    fGun->SetParticleDefinition(ion);
    fGun->SetParticleCharge(fIonChargeState * eplus);
    fNeedsIonLookup = false;
  }
  fGun->GeneratePrimaryVertex(event);
}

void PrimaryGeneratorAction::GenerateReactionPrimaries(G4Event* event) {
  const ReactionKinematics::Event ev = fReaction->GenerateEvent();

  const G4ThreeVector vertexPos(fReactionX0Cm * cm, fReactionY0Cm * cm, fReactionZ0Cm * cm);
  auto* vertex = new G4PrimaryVertex(vertexPos, 0.0);

  G4ParticleDefinition* recoilIon = G4IonTable::GetIonTable()->GetIon(fRecoilZ, fRecoilA, 0.0);
  auto* recoil = new G4PrimaryParticle(recoilIon, ev.recoilMomentumMeV.x() * MeV,
                                        ev.recoilMomentumMeV.y() * MeV,
                                        ev.recoilMomentumMeV.z() * MeV);
  recoil->SetCharge(fRecoilChargeState * eplus);
  vertex->SetPrimary(recoil);

  // One gamma per cascade step actually sampled this event (see
  // ReactionKinematics::GenerateEvent) -- usually 1, but 2 for this
  // pilot's own reaction whenever the resonance decays through its
  // level-1 intermediate state first (see reactions/o15ag_19ne.reaction).
  for (const G4ThreeVector& pGamma : ev.gammaMomentaMeV) {
    auto* gamma = new G4PrimaryParticle(G4Gamma::Gamma(), pGamma.x() * MeV, pGamma.y() * MeV,
                                         pGamma.z() * MeV);
    vertex->SetPrimary(gamma);
  }

  event->AddPrimaryVertex(vertex);
}
