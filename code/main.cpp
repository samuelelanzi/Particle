// header file
#include "Particle.hpp"
#include "ParticleType.hpp"
#include "ResonanceType.hpp"

// ROOT's header file
#include "TCanvas.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1F.h"
#include "TRandom.h"
#include "TStyle.h"

int main() {
  using namespace std;

  Particle::AddParticle("pion+", 0.13957, 1, 0.);
  Particle::AddParticle("pion-", 0.13957, -1, 0.);

  Particle::AddParticle("kaon+", 0.49367, 1, 0.);
  Particle::AddParticle("kaon-", 0.49367, -1, 0.);

  Particle::AddParticle("proton+", 0.93827, 1, 0.);
  Particle::AddParticle("proton-", 0.93827, -1, 0.);

  Particle::AddParticle("K*", 0.89166, 0, 0.050);

  vector<Particle> particle_v{};

  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(111);

  TH1F* hPt = new TH1F("hPT", "Particle Types Distribution", 7, 0, 7);
  TH1F* hPhi = new TH1F("hPhi", "Phi Distribution", 100, 0., 2 * M_PI);
  TH1F* hTheta = new TH1F("hTheta", "Theta Distribution", 100, 0., M_PI);

  TH1F* hP = new TH1F("hP", "Momentum Distribution", 80, 0, 7);
  TH1F* hPtr = new TH1F("hPtr", "Trasversal Momentum Distribution", 80, 0, 5);
  TH1F* hE = new TH1F("hE", "Energy Distribution", 80, 0, 6);
  TH1F* hMass = new TH1F("hMass", "Mass Invariant Distribution", 320, 0, 4);

  TH1F* hMass_dc = new TH1F("hMass_dc", "Mass Invariant Distribution with Discordant Charges", 320, 0, 4);
  TH1F* hMass_sc = new TH1F("hMass_sc", "Mass Invariant Distribution with Same Charges", 320, 0, 4);
  TH1F* hMass_pkD = new TH1F("hMass_pkD", "Mass Invariant Distribution Pion+/Kaon- Pion-/Kaon+", 160, 0, 2);
  TH1F* hMass_pkC = new TH1F("hMass_pkC", "Mass Invariant Distribution Pion+/Kaon+ Pion-/Kaon-", 160, 0, 2);

  TH1F* hMass_pkDecay = new TH1F("hMass_pkDecay", "Mass Invariant Distribution Decay K* in Pion+/Kaon- Pion-/Kaon+", 320, 0, 4);
  TH1F* pkD_pkC = new TH1F("pkD_pkC", "Difference between Pion+/Kaon- Pion-/Kaon+ Minus Pion+/Kaon+ Pion-/Kaon-", 160, 0, 2);
  TH1F* dc_sc = new TH1F("dc_sc", "Difference between Discordant Charges and Same Charges", 320, 0, 4);

  random_device rd;
  mt19937 gen(rd());
  uniform_real_distribution<> distrib(0, 1);
  uniform_int_distribution<> int_distrib(0, 1);

  gRandom->SetSeed();

  for(int i = 0; i != 1e5; ++i) {
    for(int j = 0; j != 1e2; ++j) {
      int prob_type_decay = int_distrib(gen);
      double prob_type = distrib(gen);

      double phi = gRandom->Uniform(0., 2 * M_PI);
      hPhi->Fill(phi);

      double theta = gRandom->Uniform(0., M_PI);
      hTheta->Fill(theta);

      double p_ = gRandom->Exp(1);
      hP->Fill(p_);

      P linearMomentum;

      if(prob_type <= 0.4) {
        Particle pion {"pion+", linearMomentum};
        pion.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(pion);
        hPt->Fill(pion.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pion.Energy());
      }

      else if(prob_type > 0.4 && prob_type <= 0.8) {
        Particle pion_ {"pion-", linearMomentum};
        pion_.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(pion_);
        hPt->Fill(pion_.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(pion_.Energy());
      }

      else if (prob_type > 0.8 && prob_type <= 0.85) {
        Particle kaon{"kaon+", linearMomentum};
        kaon.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(kaon);
        hPt->Fill(kaon.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(kaon.Energy());
      }

      else if (prob_type > 0.85 && prob_type <= 0.90) {
        Particle kaon_{"kaon-", linearMomentum};
        kaon_.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(kaon_);
        hPt->Fill(kaon_.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(kaon_.Energy());
      }

      else if (prob_type > 0.9 && prob_type <= 0.945) {
        Particle proton{"proton+", linearMomentum};
        proton.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(proton);
        hPt->Fill(proton.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(proton.Energy());
      }

      else if (prob_type > 0.945 && prob_type <= 0.99) {
        Particle proton_{"proton-", linearMomentum};
        proton_.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(proton_);
        hPt->Fill(proton_.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(proton_.Energy());
      }

      else {
        Particle k_star{"K*", linearMomentum};
        k_star.setP(p_ * sin(theta) * cos(phi), p_ * sin(theta) * sin(phi), p_ * cos(theta));
        particle_v.push_back(k_star);
        hPt->Fill(k_star.getIParticle());
        hPtr->Fill(sqrt(pow(p_ * sin(theta) * cos(phi), 2) + pow(p_ * sin(theta) * sin(phi), 2)));
        hE->Fill(k_star.Energy());

        if(prob_type_decay != 0) {
          Particle pionD {"pion+", linearMomentum};
          Particle kaonD_ {"kaon-", linearMomentum};
          k_star.Decay2Body(pionD, kaonD_);
          particle_v.push_back(pionD);
          particle_v.push_back(kaonD_);
          hMass_pkDecay->Fill(pionD.invMass(kaonD_));
        } else {
          Particle pionD_ {"pion-", linearMomentum};
          Particle kaonD {"kaon+", linearMomentum};
          k_star.Decay2Body(pionD_, kaonD);
          particle_v.push_back(pionD_);
          particle_v.push_back(kaonD);
          hMass_pkDecay->Fill(pionD_.invMass(kaonD));
        }
      }
    }

    for(int i = 0; i != static_cast<int>(particle_v.size()); ++i) {
      for(int j = i + 1; j != static_cast<int>(particle_v.size()); ++j) {
        hMass->Fill(particle_v[i].invMass(particle_v[j]));

        if(particle_v[i].getCharge() == particle_v[j].getCharge()) {
          hMass_sc->Fill(particle_v[i].invMass(particle_v[j]));
        } else {
          hMass_dc->Fill(particle_v[i].invMass(particle_v[j]));
        }

        if((particle_v[i].getName() == "pion+" && particle_v[j].getName() == "kaon-") || (particle_v[i].getName() == "kaon-" && particle_v[j].getName() == "pion+")) {
          hMass_pkD->Fill(particle_v[i].invMass(particle_v[j]));
        } else if((particle_v[i].getName() == "pion-" && particle_v[j].getName() == "kaon+") || (particle_v[i].getName() == "kaon+" && particle_v[j].getName() == "pion-")) {
          hMass_pkD->Fill(particle_v[i].invMass(particle_v[j]));
        }

        if((particle_v[i].getName() == "pion+" && particle_v[j].getName() == "kaon+") || (particle_v[i].getName() == "kaon+" && particle_v[j].getName() == "pion+")) {
          hMass_pkC->Fill(particle_v[i].invMass(particle_v[j]));
        } else if((particle_v[i].getName() == "pion-" && particle_v[j].getName() == "kaon-") || (particle_v[i].getName() == "kaon-" && particle_v[j].getName() == "pion-")) {
          hMass_pkC->Fill(particle_v[i].invMass(particle_v[j]));
        }
      }
    } 
    particle_v.clear(); 
  }

  TFile* particle = new TFile("particle.root", "CREATE");
  {
    particle->Write();
    hPt->Write();
    hPhi->Write();
    hTheta->Write();
    hP->Write();
    hPtr->Write();
    hE->Write();
    hMass->Write();
    hMass_dc->Write();
    hMass_sc->Write();
    hMass_pkD->Write();
    hMass_pkC->Write();
    hMass_pkDecay->Write();
    particle->Close();
  }

  TCanvas* c1 = new TCanvas("c1", "Particle Types & Angles Distribution", 100, 100, 1100, 700);
  {
    c1->Divide(1, 2);

    c1->cd(1);
    hPt->Draw("E");
    hPt->Draw("SAME");

    auto c1_2 = c1->cd(2);
    c1_2->Divide(2, 1);

    c1_2->cd(1);
    hPhi->Fit("pol0");
    hPhi->Draw("E");
    hPhi->Draw("SAME");

    c1_2->cd(2);
    hTheta->Fit("pol0");
    hTheta->Draw("E");
    hTheta->Draw("SAME");
  }

  TCanvas* c2 = new TCanvas("c2", "Momentum, Energy and Mass Invariant", 200, 100, 1100, 700);
  {
    c2->Divide(2, 2);

    c2->cd(1);
    hP->Fit("expo");
    hP->Draw("E");
    hP->Draw("SAME");

    c2->cd(2);
    hPtr->Draw("E");
    hPtr->Draw("SAME");

    c2->cd(3);
    hE->Draw("E");
    hE->Draw("SAME");

    c2->cd(4);
    hMass->Draw();
  }

  TCanvas* c3 = new TCanvas("c3", "Mass Invariant between certain Particle Types", 300, 100, 1100, 700);
  {
    c3->Divide(2, 2);

    c3->cd(1);
    hMass_dc->Draw();

    c3->cd(2);
    hMass_sc->Draw();

    c3->cd(3);
    hMass_pkD->Draw();

    c3->cd(4);
    hMass_pkC->Draw();
  }

  TCanvas* c4 = new TCanvas("c4", "Mass Invariant Decay", 400, 100, 1100, 700);
  {
    c4->Divide(1, 3);

    c4->cd(1);
    hMass_pkDecay->Fit("gaus");
    hMass_pkDecay->Draw();

    c4->cd(2);
    pkD_pkC->Sumw2();
    pkD_pkC->Add(hMass_pkD, hMass_pkC, 1, -1);
    pkD_pkC->SetEntries(pkD_pkC->Integral());
    pkD_pkC->Fit("gaus");
    pkD_pkC->Draw("HIST, SAME");

    c4->cd(3);
    dc_sc->Add(hMass_dc, hMass_sc, 1, -1);
    dc_sc->Draw();
  }
}