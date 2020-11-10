// header file
#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "invMass.hpp"
#include "rndmCharge.hpp"

// ROOT's header file
#include "TCanvas.h"
#include "TH1.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TStyle.h"

// STL c++
#include <random>

int main() {
  using namespace std;

  ParticleType* pion = new ParticleType{"pion", 0.13957, 0}; // 0 is default
  ParticleType* kaon = new ParticleType{"kaon", 0.49367, 0};
  ParticleType* proton = new ParticleType{"proton", 0.93827, 0};
  ParticleType* K_s = new ParticleType{"K*", 0.89166, 0};
  ResonanceType* K_resonance = new ResonanceType{*K_s, 0.050};

  vector<Particle> particle_v{};

  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(111);

  TH1F* hPhi = new TH1F("hPhi", "Phi Distribution", 100, 0., 2 * M_PI);
  TH1F* hTheta = new TH1F("hTheta", "Theta Distribution", 100, 0., M_PI);
  TH1F* hP = new TH1F("hP", "Momentum Distribution", 1000, 0, 7);
  TH1F* hPtr = new TH1F("hPtr", "Trasversal Momentum Distribution", 1000, 0, 5);
  TH1F* hE = new TH1F("hE", "Energy Distribution", 1000, 0, 6);
  TH1F* hPT = new TH1F("hPT", "Particle Types Distribution", 4, 0, 5);
  TH1F* hMass = new TH1F("hMass", "Mass Invariant Distribution", 100, 0, 100);

  random_device rd;
  mt19937 gen(rd());
  uniform_int_distribution<> distrib(1, 100);
  int seed = distrib(gen);

  gRandom->SetSeed(seed);

  for (int i = 0; i != 1e5; ++i) {
    for (int j = 0; j != 1e2; ++j) {
      int prob_type = distrib(gen);
      int charge = rndmCharge(prob_type);

      double phi = gRandom->Uniform(0., 2 * M_PI);
      hPhi->Fill(phi);
      double theta = gRandom->Uniform(0., M_PI);
      hTheta->Fill(theta);
      double p_ = gRandom->Exp(1);
      hP->Fill(p_);

      P pi_linearMomentum;
      P ka_linearMomentum;
      P pr_linearMomentum;
      P ks_linearMomentum;

      if (prob_type <= 80) {
        hPT->Fill(1);
        pion->setCharge(charge);
        Particle pi{pion, "pion", pi_linearMomentum};
        pi.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pi.Energy());
        particle_v.push_back(pi);
      }

      else if (prob_type > 80 && prob_type <= 90) {
        hPT->Fill(2);
        kaon->setCharge(charge);
        Particle ka{kaon, "kaon", ka_linearMomentum};
        ka.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(ka.Energy());
        particle_v.push_back(ka);
      }

      else if (prob_type > 90 && prob_type <= 99) {
        hPT->Fill(3);
        proton->setCharge(charge);
        Particle pr{proton, "proton", pr_linearMomentum};
        pr.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pr.Energy());
        particle_v.push_back(pr);
      }

      else {
        hPT->Fill(4);
        Particle ks{K_s, "K*", ka_linearMomentum};
        ks.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(ks.Energy());
        particle_v.push_back(ks);
      }
    }
  }

  TCanvas* cPT =
      new TCanvas("cPT", "Particle Types Distribution", 100, 100, 1100, 700);

  hPT->Draw();

  TCanvas* cAngles =
      new TCanvas("cAngles", "Angles Distribution", 200, 100, 1100, 700);

  {
    cAngles->Divide(1, 2);

    cAngles->cd(1);
    hPhi->Fit("pol0");
    hPhi->Draw();

    cAngles->cd(2);
    hTheta->Fit("pol0");
    hTheta->Draw();
  }

  TCanvas* cPE = new TCanvas("cPE", "Momentum & Energy", 300, 100, 1100, 700);

  {
    cPE->Divide(1, 3);

    cPE->cd(1);
    hP->Fit("expo");
    hP->Draw();

    cPE->cd(2);
    hPtr->Draw();

    cPE->cd(3);
    hE->Draw();
  }
}