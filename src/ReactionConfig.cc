#include "ReactionConfig.hh"

#include <fstream>
#include <sstream>
#include <stdexcept>

namespace {

// Reads whitespace-separated tokens off a card's own line (the card
// keyword itself already consumed by the caller).
std::vector<std::string> Tokenize(const std::string& line) {
  std::istringstream iss(line);
  std::vector<std::string> tokens;
  std::string tok;
  while (iss >> tok) tokens.push_back(tok);
  return tokens;
}

}  // namespace

ReactionConfig ReactionConfig::Load(const std::string& path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("ReactionConfig::Load: cannot open '" + path + "'");
  }

  ReactionConfig cfg;
  bool haveBeam = false, haveTarg = false, haveRecl = false, haveEres = false;

  std::string line;
  while (std::getline(in, line)) {
    std::istringstream lineStream(line);
    std::string card;
    lineStream >> card;
    if (card.empty() || card == "COMM") continue;
    if (card == "SENT") break;

    const std::string rest = line.substr(line.find(card) + card.size());
    const std::vector<std::string> tok = Tokenize(rest);

    if (card == "BEAM" || card == "TARG" || card == "RECL") {
      if (tok.size() < 3) {
        throw std::runtime_error("ReactionConfig::Load: '" + path + "': " + card +
                                  " needs Z A massExcessMeV" +
                                  (card == "RECL" ? " chargeState" : ""));
      }
      ParticleSpec spec;
      spec.Z = std::stoi(tok[0]);
      spec.A = std::stoi(tok[1]);
      spec.massExcessMeV = std::stod(tok[2]);
      if (card == "BEAM") {
        cfg.beam = spec;
        haveBeam = true;
      } else if (card == "TARG") {
        cfg.target = spec;
        haveTarg = true;
      } else {
        if (tok.size() < 4) {
          throw std::runtime_error("ReactionConfig::Load: '" + path +
                                    "': RECL needs Z A massExcessMeV chargeState");
        }
        cfg.recoil = spec;
        cfg.recoilChargeState = std::stoi(tok[3]);
        haveRecl = true;
      }
    } else if (card == "ERES") {
      if (tok.empty()) {
        throw std::runtime_error("ReactionConfig::Load: '" + path + "': ERES needs a value");
      }
      cfg.resonanceEnergyMeV = std::stod(tok[0]);
      haveEres = true;
    } else if (card == "RWID") {
      if (tok.empty()) {
        throw std::runtime_error("ReactionConfig::Load: '" + path + "': RWID needs a value");
      }
      cfg.resonanceWidthMeV = std::stod(tok[0]);
    } else if (card == "BKIN") {
      if (tok.empty()) {
        throw std::runtime_error("ReactionConfig::Load: '" + path + "': BKIN needs a value");
      }
      cfg.beamEntranceKineticEnergyMeV = std::stod(tok[0]);
    } else if (card == "RTUN") {
      if (tok.empty()) {
        throw std::runtime_error("ReactionConfig::Load: '" + path + "': RTUN needs a value");
      }
      cfg.magneticFieldRetuneScale = std::stod(tok[0]);
    } else if (card == "LEVL") {
      if (tok.size() < 3) {
        throw std::runtime_error("ReactionConfig::Load: '" + path +
                                  "': LEVL needs level excitationMeV lifetimeSeconds");
      }
      const int level = std::stoi(tok[0]);
      Level lvl;
      lvl.excitationMeV = std::stod(tok[1]);
      lvl.lifetimeS = std::stod(tok[2]);
      cfg.levels[level] = lvl;
    } else if (card == "BRAT") {
      if (tok.size() < 3) {
        throw std::runtime_error("ReactionConfig::Load: '" + path +
                                  "': BRAT needs fromLevel percent toLevel");
      }
      Branch b;
      b.fromLevel = std::stoi(tok[0]);
      b.percent = std::stod(tok[1]);
      b.toLevel = std::stoi(tok[2]);
      cfg.branches.push_back(b);
    } else {
      throw std::runtime_error("ReactionConfig::Load: '" + path + "': unknown card '" + card + "'");
    }
  }

  if (!haveBeam || !haveTarg || !haveRecl || !haveEres) {
    throw std::runtime_error("ReactionConfig::Load: '" + path +
                              "' is missing one of BEAM/TARG/RECL/ERES");
  }
  if (cfg.branches.empty()) {
    throw std::runtime_error("ReactionConfig::Load: '" + path +
                              "' has no BRAT cards -- the resonance (-1) needs at least one");
  }
  return cfg;
}
