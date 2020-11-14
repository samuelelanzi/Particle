#ifndef PARTICLETYPE_HPP
#define PARTICLETYPE_HPP
#include <iostream>

class ParticleType {
protected:
  std::string fName;
  double fMass;
  int fCharge;

public:
  ParticleType(std::string name, double mass, int charge);
  std::string getName();
  double getMass();
  int getCharge();
  virtual double getWidth();
  virtual void Print();
};

#endif