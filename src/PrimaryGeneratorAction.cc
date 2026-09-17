#include "PrimaryGeneratorAction.hh"

#include <algorithm>
#include <cmath>

#include "G4Event.hh"
#include "G4Gamma.hh"
#include "G4IonTable.hh"
#include "G4Material.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTable.hh"
#include "G4PhysicalConstants.hh"
#include "G4PrimaryParticle.hh"
#include "G4PrimaryVertex.hh"
#include "G4SystemOfUnits.hh"
#include "G4ThreeVector.hh"
#include "GasStoppingPower.hh"
#include "Randomize.hh"
#include "ReactionConfig.hh"
#include "TargetChamber.hh"

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
                                                int ionChargeState, double angleXDeg)
    : fNeedsIonLookup(true), fIonZ(ionZ), fIonA(ionA), fIonChargeState(ionChargeState) {
  fGun = new G4ParticleGun(1);
  const double angleXRad = angleXDeg * CLHEP::pi / 180.0;
  fGun->SetParticleMomentumDirection(
      G4ThreeVector(std::sin(angleXRad), 0.0, std::cos(angleXRad)));
  fGun->SetParticleMomentum(momentumMeV * MeV);
  fGun->SetParticlePosition(G4ThreeVector(x0Cm * cm, y0Cm * cm, z0Cm * cm));
}

PrimaryGeneratorAction::PrimaryGeneratorAction(const std::string& reactionFilePath, double x0Cm,
                                               double y0Cm, double z0Cm)
    : fIsReaction(true),
      fConfig(ReactionConfig::Load(reactionFilePath)),
      fReactionX0Cm(x0Cm),
      fReactionY0Cm(y0Cm),
      fReactionZ0Cm(z0Cm) {
  fReaction = new ReactionKinematics(fConfig);
  fRecoilZ = fConfig.recoil.Z;
  fRecoilA = fConfig.recoil.A;
  fRecoilChargeState = fConfig.recoilChargeState;
}

PrimaryGeneratorAction::~PrimaryGeneratorAction() {
  delete fGun;
  delete fReaction;
  delete fStoppingPower;
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
  double vertexZCm = fReactionZ0Cm;
  double beamKineticEnergyMeV = fReaction->BeamKineticEnergyMeV();

  // BKIN given: track the beam's own real energy loss through the target
  // gas (see ReactionConfig.hh/GasStoppingPower.hh) instead of firing
  // every event from the same fixed point at the same fixed (energy-
  // loss-free) beam energy.
  if (fConfig.beamEntranceKineticEnergyMeV > 0.0) {
    if (!fStoppingPowerReady) {
      G4Material* targetGas = G4Material::GetMaterial(TargetChamber::TargetGasMaterialName());
      fStoppingPower =
          new GasStoppingPower(fConfig.beam.Z, fConfig.beam.A, targetGas,
                                fConfig.beamEntranceKineticEnergyMeV,
                                2.0 * TargetChamber::BeamPathHalfLengthCm());
      fStoppingPowerReady = true;
    }
    beamKineticEnergyMeV = fReaction->SampleBeamKineticEnergyMeV(
        fStoppingPower->ExitKineticEnergyMeV(), fStoppingPower->EntranceKineticEnergyMeV());
    const double meanDepthCm = fStoppingPower->DepthForKineticEnergyCm(beamKineticEnergyMeV);

    // Straggling (see GasStoppingPower's own header comment): the mean
    // dE/dx curve alone says every ion of the same entrance energy
    // reaches beamKineticEnergyMeV at exactly meanDepthCm -- real energy
    // loss is a statistical process, so a real ion actually reaches it
    // at a depth that fluctuates around that mean. Model this by
    // perturbing the *depth-crossing* energy (not the reaction's own
    // kinematics, which stay at the Breit-Wigner-sampled
    // beamKineticEnergyMeV above -- that's still the energy the reaction
    // itself happens at) by a Gaussian sample with Bohr's straggling
    // sigma at that depth, then re-finding the depth a real ion reaches
    // *that* energy at.
    const double sigmaMeV = fStoppingPower->StragglingSigmaMeV(meanDepthCm);
    double straggledKeMeV = beamKineticEnergyMeV + G4RandGauss::shoot(0.0, sigmaMeV);
    straggledKeMeV = std::min(std::max(straggledKeMeV, fStoppingPower->ExitKineticEnergyMeV()),
                               fStoppingPower->EntranceKineticEnergyMeV());
    const double depthCm = fStoppingPower->DepthForKineticEnergyCm(straggledKeMeV);
    vertexZCm = -TargetChamber::BeamPathHalfLengthCm() + depthCm;
  }

  const ReactionKinematics::Event ev = fReaction->GenerateEvent(beamKineticEnergyMeV);

  const G4ThreeVector vertexPos(fReactionX0Cm * cm, fReactionY0Cm * cm, vertexZCm * cm);
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
