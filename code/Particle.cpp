#include "Particle.hpp"

std::vector<ParticleType*> Particle::fParticleType(1, nullptr);

int Particle::fNParticleType = 0;

int Particle::FindParticle(std::string name) {
  for(int i = 0; i != static_cast<int>(fParticleType.size()) - 1; ++i) {
    if(fParticleType[i]->getName() == name) {
      return i;
    }
  }
  return -1;
}

void Particle::Boost(double bx, double by, double bz) {
  double energy = Energy();

  double b2 = bx * bx + by * by + bz * bz;
  double gamma = 1.0 / std::sqrt(1.0 - b2);
  double bp = bx * fP.fPx + by * fP.fPy + bz * fP.fPz;
  double gamma2 = b2 > 0 ? (gamma - 1.0) / b2 : 0.0;

  fP.fPx += gamma2 * bp * bx + gamma * bx * energy;
  fP.fPy += gamma2 * bp * by + gamma * by * energy;
  fP.fPz += gamma2 * bp * bz + gamma * bz * energy;
}

Particle::Particle(std::string name, P p) : fP{p} {
  fIParticle = FindParticle(name);
} 

int Particle::getIParticle() { return fIParticle; }

void Particle::AddParticle(std::string pt_name, double pt_mass, int pt_charge, double pt_width) {
  ++fNParticleType;
  fParticleType[fNParticleType - 1] = new ParticleType{pt_name, pt_mass, pt_charge};
  if(pt_width != 0.) {
    fParticleType[fNParticleType - 1] = new ResonanceType{*fParticleType[fNParticleType -1], pt_width};
  }
  fParticleType.push_back(nullptr);
}

void Particle::Print() {
  for(int i = 0; i != static_cast<int>(fParticleType.size()) - 1; ++i) {
    std::cout << "Name: " << fParticleType[i]->getName() << '\n';
    std::cout << "Mass: " << fParticleType[i]->getMass() << '\n';
    std::cout << "Charge: " << fParticleType[i]->getCharge() << '\n';
    std::cout << "Width: " << fParticleType[i]->getWidth() << "\n\n";
  }
}

void Particle::PrintParticle() {
  std::cout << "Index: " << fIParticle << '\n';
  std::cout << "Name: " << fParticleType[fIParticle]->getName() << '\n';
  std::cout << "Px: " << fP.fPx << '\n';
  std::cout << "Py: " << fP.fPy << '\n';
  std::cout << "Pz: " << fP.fPz << "\n\n";
}

P Particle::getP() { return fP; }

std::string Particle::getName() {
  return fParticleType[fIParticle]->getName();
}

double Particle::getMass() { 
  return fParticleType[fIParticle]->getMass();
}

int Particle::getCharge() {
  return fParticleType[fIParticle]->getCharge();
}

double Particle::Energy() {
  double m = getMass();
  double p_x = fP.fPx;
  double p_y = fP.fPy;
  double p_z = fP.fPz;
  return std::sqrt((m * m) + (p_x * p_x) + (p_y * p_y) + (p_z * p_z));
}

double Particle::invMass(Particle& p) {
  double sumE = Energy() + p.Energy();

  double sumPx = getP().fPx + p.getP().fPx;
  double sumPy = getP().fPy + p.getP().fPy;
  double sumPz = getP().fPz + p.getP().fPz;

  return std::sqrt(sumE * sumE - sumPx * sumPx - sumPy * sumPy - sumPz * sumPz);
}

void Particle::setP(double px, double py, double pz) {
  fP.fPx = px;
  fP.fPy = py;
  fP.fPz = pz;
}

int Particle::Decay2Body(Particle &dau1, Particle &dau2) {
  if (getMass() == 0.0) {
    std::cout << "Decayment cannot be preformed if mass is zero\n";
    return 1;
  }

  double massMot = getMass();
  double massDau1 = dau1.getMass();
  double massDau2 = dau2.getMass();

  if (fIParticle > -1) {
    float x1, x2, w, y1, y2;

    double invnum = 1. / RAND_MAX;

    do {
      x1 = 2.0 * rand() * invnum - 1.0;
      x2 = 2.0 * rand() * invnum - 1.0;
      w = x1 * x1 + x2 * x2;
    } while (w >= 1.0);

    w = std::sqrt((-2.0 * std::log(w)) / w);
    y1 = x1 * w;
    y2 = x2 * w;
    massMot += fParticleType[fIParticle]->getWidth() * y1;
  }

  if (massMot < massDau1 + massDau2) {
    std::cout << "Decayment cannot be preformed because mass is too low in this channel\n";
    return 2;
  }

  double pout = std::sqrt((massMot * massMot - (massDau1 + massDau2) * (massDau1 + massDau2)) * (massMot * massMot - (massDau1 - massDau2) * (massDau1 - massDau2))) / massMot * 0.5;

  double norm = 2 * M_PI / RAND_MAX;

  double phi = rand() * norm;
  double theta = rand() * norm * 0.5 - M_PI / 2.;
  dau1.setP(pout * std::sin(theta) * std::cos(phi), pout * std::sin(theta) * std::sin(phi), pout * std::cos(theta));
  dau2.setP(-pout * std::sin(theta) * std::cos(phi), -pout * std::sin(theta) * std::sin(phi), -pout * std::cos(theta));

  double energy = std::sqrt(fP.fPx * fP.fPx + fP.fPy * fP.fPy + fP.fPz * fP.fPz + massMot * massMot);

  double bx = fP.fPx / energy;
  double by = fP.fPy / energy;
  double bz = fP.fPz / energy;

  dau1.Boost(bx, by, bz);
  dau2.Boost(bx, by, bz);

  return 0;
}