#ifndef PARTICLE_HPP
#define PARTICLE_HPP

#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <random>

struct P {
  double fPx = 0;
  double fPy = 0;
  double fPz = 0;
};

class Particle {
private:
  static std::vector<ParticleType*> fParticleType;
  std::string fName;
  static int fNParticleType; 
  int fIParticle; 
  P fP;
  static int FindParticle(std::string name);
  void Boost(double bx, double by, double bz);

public:
  Particle(std::string name, P p);
  int getIParticle();
  static void AddParticle(std::string pt_name, double pt_mass, int pt_charge, double pt_width);
  static void Print();
  void PrintParticle();
  P getP();
  std::string getName();
  double getMass();
  int getCharge();
  double Energy();
  double invMass(Particle& p);
  void setP(double px, double py, double pz);
  int Decay2Body(Particle& dau1, Particle& dau2);
};

#endif