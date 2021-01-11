void setStyle() {
  // Cosmetics
  gROOT->SetStyle("Plain");
  gStyle->SetOptStat(112210);
  gStyle->SetOptFit(111);
  gStyle->SetOptTitle(0);
}

void analysis() {
  // Open Root File
  TFile* f = new TFile("particle.root");

  // Find Histograms
  auto hPt = (TH1F*)f->Get("hPt");
  auto hPhi = (TH1F*)f->Get("hPhi");
  auto hTheta = (TH1F*)f->Get("hTheta");
  auto hP = (TH1F*)f->Get("hP");
  auto hPtr = (TH1F*)f->Get("hPtr");
  auto hE = (TH1F*)f->Get("hE");
  auto hMass = (TH1F*)f->Get("hMass");
  auto hMass_dc = (TH1F*)f->Get("hMass_dc");
  auto hMass_sc = (TH1F*)f->Get("hMass_sc");
  auto hMass_pkD = (TH1F*)f->Get("hMass_pkD");
  auto hMass_pkC = (TH1F*)f->Get("hMass_pkC");
  auto hMass_pkDecay = (TH1F*)f->Get("hMass_pkDecay");

  TH1F *pkD_pkC = new TH1F("pkD_pkC", "Difference between Pion+/Kaon- Pion-/Kaon+ Minus Pion+/Kaon+ Pion-/Kaon-", 160, 0, 2);
  TH1F *dc_sc = new TH1F("dc_sc", "Difference between Discordant Charges and Same Charges", 160, 0, 2);

  // First Canvas
  TCanvas *c1 = new TCanvas("c1", "Particle Types, Angles and Momentum Distributioins", 100, 100, 1100, 700);
  {
    c1->Divide(2, 2);

    c1->cd(1);
    hPt->GetXaxis()->SetTitle("Types of Particles");
    hPt->GetYaxis()->SetTitle("Entries");
    hPt->GetXaxis()->CenterTitle();
    hPt->GetXaxis()->CenterTitle();
    hPt->SetLineColor(4);
    hPt->SetFillColor(4);
    hPt->SetFillStyle(3002);
    hPt->Draw("E");
    hPt->Draw("SAME");

    c1->cd(2);
    hPhi->GetXaxis()->SetTitle("Phi Angle");
    hPhi->GetYaxis()->SetTitle("Entries");
    hPhi->GetXaxis()->CenterTitle();
    hPhi->GetXaxis()->CenterTitle();
    hPhi->SetLineColor(4);
    hPhi->SetFillColor(4);
    hPhi->SetFillStyle(3002);
    hPhi->Fit("pol0");
    hPhi->Draw("E");
    hPhi->Draw("SAME");

    c1->cd(3);
    hTheta->GetXaxis()->SetTitle("Theta Angle");
    hTheta->GetYaxis()->SetTitle("Entries");
    hTheta->GetXaxis()->CenterTitle();
    hTheta->GetXaxis()->CenterTitle();
    hTheta->SetLineColor(4);
    hTheta->SetFillColor(4);
    hTheta->SetFillStyle(3002);
    hTheta->Fit("pol0");
    hTheta->Draw("E");
    hTheta->Draw("SAME");

    c1->cd(4);
    hP->GetXaxis()->SetTitle("Momentum");
    hP->GetYaxis()->SetTitle("Entries");
    hP->GetXaxis()->CenterTitle();
    hP->GetXaxis()->CenterTitle();
    hP->SetLineColor(4);
    hP->SetFillColor(4);
    hP->SetFillStyle(3002);
    hP->Fit("expo");
    hP->Draw("E");
    hP->Draw("SAME");

    c1->Print("c1.gif");
  }

  //Second Canvas
  TCanvas *c2 = new TCanvas("c2", "Mass Invariant", 400, 100, 1100, 700);
  {
    c2->Divide(1, 3);

    c2->cd(1);
    hMass_pkDecay->GetXaxis()->SetTitle("Mass Invariant #pi K Decay");
    hMass_pkDecay->GetYaxis()->SetTitle("Entries");
    hMass_pkDecay->GetXaxis()->CenterTitle();
    hMass_pkDecay->GetXaxis()->CenterTitle();
    hMass_pkDecay->SetLineColor(4);
    hMass_pkDecay->SetFillColor(4);
    hMass_pkDecay->SetFillStyle(3002);
    hMass_pkDecay->Fit("gaus");
    hMass_pkDecay->Draw("E, SAME");

    c2->cd(2);
    pkD_pkC->GetXaxis()->SetTitle("Mass Invariant #pi K Discordant - #pi K Concordant");
    pkD_pkC->GetYaxis()->SetTitle("Entries");
    pkD_pkC->GetXaxis()->CenterTitle();
    pkD_pkC->GetXaxis()->CenterTitle();
    pkD_pkC->SetLineColor(4);
    pkD_pkC->SetFillColor(4);
    pkD_pkC->SetFillStyle(3002);
    pkD_pkC->Sumw2();
    pkD_pkC->Add(hMass_pkD, hMass_pkC, 1, -1);
    pkD_pkC->SetEntries(pkD_pkC->Integral());
    pkD_pkC->Fit("gaus");
    pkD_pkC->Draw("HIST, SAME");

    c2->cd(3);
    dc_sc->GetXaxis()->SetTitle("Mass Invariant Discordant - Concordant");
    dc_sc->GetYaxis()->SetTitle("Entries");
    dc_sc->GetXaxis()->CenterTitle();
    dc_sc->GetXaxis()->CenterTitle();
    dc_sc->SetLineColor(4);
    dc_sc->SetFillColor(4);
    dc_sc->SetFillStyle(3002);
    dc_sc->Sumw2();
    dc_sc->Add(hMass_dc, hMass_sc, 1, -1);
    dc_sc->SetEntries(dc_sc->Integral());
    dc_sc->Fit("gaus");
    dc_sc->Draw("HIST, SAME");

    c2->Print("c2.gif");
  }
}