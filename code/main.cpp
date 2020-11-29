// header file
#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"
#include "invMass.hpp"

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

  ParticleType* pion = new ParticleType{"pion+", 0.13957, 1};
  ParticleType* pion_ = new ParticleType{"pion-", 0.13957, -1};

  ParticleType* kaon = new ParticleType{"kaon+", 0.49367, 1};
  ParticleType* kaon_ = new ParticleType{"kaon-", 0.49367, -1};

  ParticleType* proton = new ParticleType{"proton+", 0.93827, 1};
  ParticleType* proton_ = new ParticleType{"proton-", 0.93827, -1};

  ParticleType* K_s = new ParticleType{"K*", 0.89166, 0};
  ResonanceType* K_resonance = new ResonanceType{*K_s, 0.050};

  vector<Particle> particle_v{};
  particle_v.reserve(1e7);

  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(111);

  TH1F* hPhi = new TH1F("hPhi", "Phi Distribution", 100, 0., 2 * M_PI);
  TH1F* hTheta = new TH1F("hTheta", "Theta Distribution", 100, 0., M_PI);
  TH1F* hP = new TH1F("hP", "Momentum Distribution", 1000, 0, 7);
  TH1F* hPtr = new TH1F("hPtr", "Trasversal Momentum Distribution", 1000, 0, 5);
  TH1F* hE = new TH1F("hE", "Energy Distribution", 1000, 0, 6);
  TH1F* hPT = new TH1F("hPT", "Particle Types Distribution", 7, 0, 8);
  TH1F* hMass = new TH1F("hMass", "Mass Invariant Distribution", 100, 0, 100);

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> distrib(0, 1);

  gRandom->SetSeed();

  for (int i = 0; i != 1e5; ++i) {
    for (int j = 0; j != 1e2; ++j) {
      double prob_type = distrib(gen);

      double phi = gRandom->Uniform(0., 2 * M_PI);
      hPhi->Fill(phi);

      double theta = gRandom->Uniform(0., M_PI);
      hTheta->Fill(theta);

      double p_ = gRandom->Exp(1);
      hP->Fill(p_);

      P linearMomentum;

      if (prob_type <= 0.4) {
        hPT->Fill(1);
        Particle pi{pion, "pion+", linearMomentum};
        pi.setIParticle(particle_v.size());
        pi.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pi.Energy());
        particle_v.push_back(pi);
      }

      else if (prob_type > 0.4 && prob_type <= 0.8) {
        hPT->Fill(2);
        Particle pi_{pion_, "pion-", linearMomentum};
        pi_.setIParticle(particle_v.size());
        pi_.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pi_.Energy());
        particle_v.push_back(pi_);
      }

      else if (prob_type > 0.8 && prob_type <= 0.85) {
        hPT->Fill(3);
        Particle ka{kaon, "kaon+", linearMomentum};
        ka.setIParticle(particle_v.size());
        ka.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(ka.Energy());
        particle_v.push_back(ka);
      }

      else if (prob_type > 0.85 && prob_type <= 0.90) {
        hPT->Fill(4);
        Particle ka_{kaon_, "kaon-", linearMomentum};
        ka_.setIParticle(particle_v.size());
        ka_.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(ka_.Energy());
        particle_v.push_back(ka_);
      }

      else if (prob_type > 0.9 && prob_type <= 0.945) {
        hPT->Fill(5);
        Particle pr{proton, "proton+", linearMomentum};
        pr.setIParticle(particle_v.size());
        pr.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pr.Energy());
        particle_v.push_back(pr);
      }

      else if (prob_type > 0.945 && prob_type <= 0.99) {
        hPT->Fill(6);
        Particle pr_{proton_, "proton-", linearMomentum};
        pr_.setIParticle(particle_v.size());
        pr_.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pr_.Energy());
        particle_v.push_back(pr_);
      }

      else {
        hPT->Fill(7);
        Particle ks{K_s, "K*", linearMomentum};
        ks.setIParticle(particle_v.size());
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
