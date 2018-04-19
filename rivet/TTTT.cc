// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {


  /// @brief Add a short analysis description here
  class TTTT : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(TTTT);


    /// @name Analysis methods
    //@{

    /// Book histograms and initialise projections before the run
    void init() {

      // Initialise and register projections
      PromptFinalState pls(ChargedLeptons(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV), true);
      declare(pls, "PromptLeptons");

      // build jets out of all stable particles that aren't invisible
      // and aren't prompt leptons
      VetoedFinalState vfs(FinalState(Cuts::abseta < 5.0 && Cuts::pT > 100*MeV));
      vfs.addVetoOnThisFinalState(pls);

      declare(FastJets(vfs, FastJets::ANTIKT, 0.4), "Jets");
      declare(FastJets(vfs, FastJets::ANTIKT, 1.0), "FatJets");


      njets = bookH("njets", 21, -0.5, 20.5, "number of jets in acceptance");
      nfatjets = bookH("nfatjets", 4, -0.5, 3.5, "number of fatjets in acceptance");
      nleps = bookH("nleps", 5, -0.5, 4.5, "number of prompt leptons in acceptance");

      njets_lJ = bookH("njets_lJ", 21, -0.5, 20.5, "number of jets in acceptance (l+J)");
      njets_lJJ = bookH("njets_lJJ", 21, -0.5, 20.5, "number of jets in acceptance (l+JJ)");
      njets_ssJ = bookH("njets_ssJ", 21, -0.5, 20.5, "number of jets in acceptance (ss+J)");

      nfatjets_lJ = bookH("nfatjets_lJ", 4, -0.5, 3.5, "number of fatjets in acceptance (l+J)");
      nfatjets_lJJ = bookH("nfatjets_lJJ", 4, -0.5, 3.5, "number of fatjets in acceptance (l+JJ)");
      nfatjets_ssJ = bookH("nfatjets_ssJ", 4, -0.5, 3.5, "number of fatjets in acceptance (ss+J)");

      naddjets_lJ = bookH("naddjets_lJ", 21, -0.5, 20.5, "number of additional jets in acceptance (l+J)");
      naddjets_lJJ = bookH("naddjets_lJJ", 21, -0.5, 20.5, "number of additional jets in acceptance (l+JJ)");
      naddjets_ssJ = bookH("naddjets_ssJ", 21, -0.5, 20.5, "number of additional jets in acceptance (ss+J)");

      ptth_lJ = bookH("ptth_lJ", 25, 0, 2*TeV, "hadronic top $p_\\mathrm{T}$ (l+J)");
      pttl_lJ = bookH("pttl_lJ", 25, 0, 2*TeV, "leptonic top $p_\\mathrm{T}$ (l+J)");
      ptth1_lJJ = bookH("ptth1_lJJ", 25, 0, 2*TeV, "leading hadronic top $p_\\mathrm{T}$ (l+JJ)");
      ptth2_lJJ = bookH("ptth2_lJJ", 25, 0, 2*TeV, "subleading hadronic top $p_\\mathrm{T}$ (l+JJ)");
      ptth_ssJ = bookH("ptth_ssJ", 25, 0, 2*TeV, "hadronic top $p_\\mathrm{T}$ (ss+J)");
      pttl_ssJ = bookH("pttl_ssJ", 25, 0, 2*TeV, "leptonic top $p_\\mathrm{T}$ (ss+J)");

      ptl1_lJ = bookH("ptl1_lJ", 25, 0, 1*TeV, "leading lepton $p_\\mathrm{T}$ (l+J)");
      ptl1_lJJ = bookH("ptl1_lJJ", 25, 0, 1*TeV, "leading lepton $p_\\mathrm{T}$ (l+JJ)");
      ptl1_ssJ = bookH("ptl1_ssJ", 25, 0, 1*TeV, "leading lepton $p_\\mathrm{T}$ (ss+J)");
      ptl2_ssJ = bookH("ptl2_ssJ", 25, 0, 1*TeV, "subleading lepton $p_\\mathrm{T}$ (ss+J)");

      mth_lJ = bookH("mth_lJ", 25, 0, 500*GeV, "hadronic top mass (l+J)");
      mtl_lJ = bookH("mtl_lJ", 25, 0, 500*GeV, "leptonic top mass (l+J)");
      mth1_lJJ = bookH("mth1_lJJ", 25, 0, 500*GeV, "leading hadronic top mass (l+JJ)");
      mth2_lJJ = bookH("mth2_lJJ", 25, 0, 500*GeV, "leading hadronic top mass (l+JJ)");
      mth_ssJ = bookH("mth_ssJ", 25, 0, 500*GeV, "hadronic top mass (ss+J)");
      mtl_ssJ = bookH("mtl_ssJ", 25, 0, 500*GeV, "leptonic top mass (ss+J)");

      dphitt_lJ = bookH("dphitt_lJ", 20, 0, 4, "$d\\phi(tt)$ (l+J)");
      dphitt_lJJ = bookH("dphitt_lJJ", 20, 0, 4, "$d\\phi(tt)$ (l+JJ)");
      dphitt_ssJ = bookH("dphitt_ssJ", 20, 0, 4, "$d\\phi(tt)$ (ss+J)");

      pttt_lJ = bookH("pttt_lJ", 25, 0, 500*GeV, "$tt$ $p_\\mathrm{T}$ (l+J)");
      pttt_lJJ = bookH("pttt_lJJ", 25, 0, 500*GeV, "$tt$ $p_\\mathrm{T}$ (l+JJ)");
      pttt_ssJ = bookH("pttt_ssJ", 25, 0, 500*GeV, "$tt$ $p_\\mathrm{T}$ (ss+J)");

      mtt_lJ = bookH("mtt_lJ", 15, 0, 3*TeV, "$tt$ invariant mass (l+J)");
      mtt_lJJ = bookH("mtt_lJJ", 15, 0, 3*TeV, "$tt$ invariant mass (l+JJ)");
      mtt_ssJ = bookH("mtt_ssJ", 15, 0, 3*TeV, "$tt$ invariant mass (ss+J)");

    }

    Histo1DPtr bookH(const string& path, double nb, double bmin, double bmax, const string& title) {
      Histo1DPtr h = bookHisto1D(path, nb, bmin, bmax, title);
      allHists.push_back(h);
      return h;
    }

    vector<const Jet*> additionalJets(const Jets& jets, const Jets& fatjets) {
      vector<const Jet*> addjets;

      for (const Jet& j : jets) {
        const FourMomentum mom = j.mom();
        bool pass = true;

        for (const Jet& fj : fatjets) {
          if (deltaR(fj.mom(), mom) > 1.2)
            continue;

          pass = false;
          break;
        }

        if (pass)
          addjets.push_back(&j);
      }

      return addjets;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& leps = apply<PromptFinalState>(event, "PromptLeptons").particles();
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV);
      const Jets& fatjets =
        apply<FastJets>(event, "FatJets").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.0 && Cuts::mass > 100*GeV);

      double weight = event.weight();

      nleps->fill(leps.size(), weight);
      njets->fill(jets.size(), weight);
      nfatjets->fill(fatjets.size(), weight);


      if (leps.size() == 1 && fatjets.size() == 1) {
        njets_lJ->fill(jets.size(), weight);
        nfatjets_lJ->fill(fatjets.size(), weight);
        ptl1_lJ->fill(leps[0].pt(), weight);

        Jets goodfatjets = fatjets;
        vector<const Jet*> goodjets = additionalJets(jets, goodfatjets);
        naddjets_lJ->fill(goodjets.size(), weight);

        const Jet* bestjet = NULL;
        double drmin = -1;
        for (const Jet* j : goodjets) {
          double dr = deltaR(leps[0].mom(), j->mom());
          if (drmin < 0 || dr < drmin) {
            bestjet = j;
            drmin = dr;
          }
        }

        if (drmin >= 0) {
          // colinear approximation
          FourMomentum tl = leps[0].mom() + leps[0].mom() + bestjet->mom();
          FourMomentum th = goodfatjets[0].mom();
          FourMomentum tt = tl + th;

          ptth_lJ->fill(th.pt(), weight);
          pttl_lJ->fill(tl.pt(), weight);

          mtl_lJ->fill(tl.mass(), weight);
          mth_lJ->fill(th.mass(), weight);

          dphitt_lJ->fill(abs(deltaPhi(tl, th)), weight);

          pttt_lJ->fill(tt.pt(), weight);
          mtt_lJ->fill(tt.mass(), weight);
        }

      } else if (leps.size() == 1 && fatjets.size() >= 2) {
        njets_lJJ->fill(jets.size(), weight);
        nfatjets_lJJ->fill(fatjets.size(), weight);
        ptth1_lJJ->fill(fatjets[0].pt(), weight);
        ptth2_lJJ->fill(fatjets[1].pt(), weight);
        mth1_lJJ->fill(fatjets[0].mass(), weight);
        mth2_lJJ->fill(fatjets[1].mass(), weight);
        ptl1_lJJ->fill(leps[0].pt(), weight);

        Jets goodfatjets;
        goodfatjets.push_back(fatjets[0]);
        goodfatjets.push_back(fatjets[1]);
        vector<const Jet*> goodjets = additionalJets(jets, goodfatjets);
        naddjets_lJJ->fill(goodjets.size(), weight);

        FourMomentum t1 = goodfatjets[0].mom();
        FourMomentum t2 = goodfatjets[1].mom();
        FourMomentum tt = t1 + t2;

        mth1_lJJ->fill(t1.mass(), weight);
        mth2_lJJ->fill(t2.mass(), weight);

        dphitt_lJJ->fill(abs(deltaPhi(t1, t2)), weight);

        pttt_lJJ->fill(tt.pt(), weight);
        mtt_lJJ->fill(tt.mass(), weight);

      } else if (leps.size() == 2 && fatjets.size() >= 1 && leps[0].threeCharge()*leps[1].threeCharge() > 0 ) {
        njets_ssJ->fill(jets.size(), weight);
        nfatjets_ssJ->fill(fatjets.size(), weight);
        ptl1_ssJ->fill(leps[0].pt(), weight);
        ptl2_ssJ->fill(leps[1].pt(), weight);

        Jets goodfatjets;
        goodfatjets.push_back(fatjets[0]);
        vector<const Jet*> goodjets = additionalJets(jets, goodfatjets);
        naddjets_ssJ->fill(goodjets.size(), weight);

        const Jet* bestjet = NULL;
        double drmin = -1;
        for (const Jet* j : goodjets) {
          double dr = deltaR(leps[0].mom(), j->mom());
          if (drmin < 0 || dr < drmin) {
            bestjet = j;
            drmin = dr;
          }
        }

        if (drmin >= 0) {
          // colinear approximation
          FourMomentum tl = leps[0].mom() + leps[0].mom() + bestjet->mom();
          FourMomentum th = goodfatjets[0].mom();
          FourMomentum tt = tl + th;

          pttl_ssJ->fill(tl.pt(), weight);
          ptth_ssJ->fill(th.pt(), weight);

          mtl_ssJ->fill(tl.mass(), weight);
          mth_ssJ->fill(th.mass(), weight);

          dphitt_ssJ->fill(abs(deltaPhi(tl, th)), weight);

          pttt_ssJ->fill(tt.pt(), weight);
          mtt_ssJ->fill(tt.mass(), weight);
        }

      }

      return;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (Histo1DPtr& h : allHists) {
        Histo1DPtr hnorm = make_shared<Histo1D>(*h);
        normalize(hnorm);
        hnorm->setPath(hnorm->path() + "_norm");
        addAnalysisObject(hnorm);

        scale(h, crossSection()/sumOfWeights());
      }

      return;
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr njets;
    Histo1DPtr nfatjets;
    Histo1DPtr nleps;

    Histo1DPtr njets_lJ;
    Histo1DPtr njets_lJJ;
    Histo1DPtr njets_ssJ;

    Histo1DPtr nfatjets_lJ;
    Histo1DPtr nfatjets_lJJ;
    Histo1DPtr nfatjets_ssJ;

    Histo1DPtr naddjets_lJ;
    Histo1DPtr naddjets_lJJ;
    Histo1DPtr naddjets_ssJ;

    Histo1DPtr ptth_lJ;
    Histo1DPtr pttl_lJ;
    Histo1DPtr ptth1_lJJ;
    Histo1DPtr ptth2_lJJ;
    Histo1DPtr ptth_ssJ;
    Histo1DPtr pttl_ssJ;

    Histo1DPtr ptl1_lJ;
    Histo1DPtr ptl1_lJJ;
    Histo1DPtr ptl1_ssJ;
    Histo1DPtr ptl2_ssJ;

    Histo1DPtr mth_lJ;
    Histo1DPtr mtl_lJ;
    Histo1DPtr mth1_lJJ;
    Histo1DPtr mth2_lJJ;
    Histo1DPtr mth_ssJ;
    Histo1DPtr mtl_ssJ;

    Histo1DPtr dphitt_lJ;
    Histo1DPtr dphitt_lJJ;
    Histo1DPtr dphitt_ssJ;

    Histo1DPtr pttt_lJ;
    Histo1DPtr pttt_lJJ;
    Histo1DPtr pttt_ssJ;

    Histo1DPtr mtt_lJ;
    Histo1DPtr mtt_lJJ;
    Histo1DPtr mtt_ssJ;

    vector<Histo1DPtr> allHists;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TTTT);


}
