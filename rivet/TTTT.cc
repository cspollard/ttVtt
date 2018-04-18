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


      njets = bookHisto1D("njets", 15, -0.5, 14.5, "number of jets in acceptance");
      nfatjets = bookHisto1D("nfatjets", 4, -0.5, 3.5, "number of fatjets in acceptance");
      nleps = bookHisto1D("nleps", 5, -0.5, 4.5, "number of prompt leptons in acceptance");

      njets_lJ = bookHisto1D("njets_lJ", 21, -0.5, 20.5, "number of jets in acceptance (l+J)");
      njets_lJJ = bookHisto1D("njets_lJJ", 21, -0.5, 20.5, "number of jets in acceptance (l+JJ)");
      njets_ssJ = bookHisto1D("njets_ssJ", 21, -0.5, 20.5, "number of jets in acceptance (ss+J)");

      nfatjets_lJ = bookHisto1D("nfatjets_lJ", 4, -0.5, 3.5, "number of fatjets in acceptance (l+J)");
      nfatjets_lJJ = bookHisto1D("nfatjets_lJJ", 4, -0.5, 3.5, "number of fatjets in acceptance (l+JJ)");
      nfatjets_ssJ = bookHisto1D("nfatjets_ssJ", 4, -0.5, 3.5, "number of fatjets in acceptance (ss+J)");

      ptJ1_lJ = bookHisto1D("ptJ1_lJ", 25, 0, 1*TeV, "leading fatjet $p_\\mathrm{T}$ (l+J)");
      ptJ1_lJJ = bookHisto1D("ptJ1_lJJ", 25, 0, 1*TeV, "leading fatjet $p_\\mathrm{T}$ (l+JJ)");
      ptJ2_lJJ = bookHisto1D("ptJ2_lJJ", 25, 0, 1*TeV, "subleading fatjet $p_\\mathrm{T}$ (l+JJ)");
      ptJ1_ssJ = bookHisto1D("ptJ1_ssll", 25, 0, 1*TeV, "leading fatjet $p_\\mathrm{T}$ (ss+J)");

      ptl1_lJ = bookHisto1D("ptl1_lJ", 25, 0, 1*TeV, "leading lepton $p_\\mathrm{T}$ (l+J)");
      ptl1_lJJ = bookHisto1D("ptl1_lJJ", 25, 0, 1*TeV, "leading lepton $p_\\mathrm{T}$ (l+JJ)");
      ptl1_ssJ = bookHisto1D("ptl1_ssJ", 25, 0, 1*TeV, "leading lepton $p_\\mathrm{T}$ (ss+J)");
      ptl2_ssJ = bookHisto1D("ptl2_ssJ", 25, 0, 1*TeV, "subleading lepton $p_\\mathrm{T}$ (ss+J)");

      mJ1_lJ = bookHisto1D("mJ_lJ", 25, 0, 500*GeV, "fatjet mass (l+J)");
      mJ1_lJJ = bookHisto1D("mJ_lJJ", 25, 0, 500*GeV, "fatjet mass (l+JJ)");
      mJ2_lJJ = bookHisto1D("mJ_lJJ", 25, 0, 500*GeV, "fatjet mass (l+JJ)");
      mJ1_ssJ = bookHisto1D("mJ_ssJ", 25, 0, 500*GeV, "fatjet mass (ss+J)");

      mlJ_lJ = bookHisto1D("mlJ_lJ", 15, 0, 3*TeV, "lepton+fatjet invariant mass (l+J)");
      mJJ_lJJ = bookHisto1D("mJJ_lJJ", 15, 0, 3*TeV, "fatjet+fatjet invariant mass (l+JJ)");
      mlJ_ssJ = bookHisto1D("mlJ_ssJ", 15, 0, 3*TeV, "lepton+fatjet invariant mass (ss+J)");

      allHists =
        { njets, nfatjets, nleps
        , njets_lJ, nfatjets_lJ
        , njets_lJJ, nfatjets_lJJ
        , njets_ssJ, nfatjets_ssJ
        , ptJ1_lJ, ptJ1_lJJ, ptJ2_lJJ, ptJ1_ssJ
        , ptl1_lJ, ptl1_lJJ, ptl1_ssJ, ptl2_ssJ
        , mJ1_lJ, mJ1_lJJ, mJ2_lJJ, mJ1_ssJ 
        , mlJ_lJ, mJJ_lJJ, mlJ_ssJ
        };
      
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& leps = apply<PromptFinalState>(event, "PromptLeptons").particles();
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV);
      const Jets& fatjets = apply<FastJets>(event, "FatJets").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.5);

      double weight = event.weight();

      nleps->fill(leps.size(), weight);
      njets->fill(jets.size(), weight);
      nfatjets->fill(fatjets.size(), weight);


      if (leps.size() == 1 && fatjets.size() == 1) {
        njets_lJ->fill(jets.size(), weight);
        nfatjets_lJ->fill(fatjets.size(), weight);
        ptJ1_lJ->fill(fatjets[0].pt(), weight);
        mJ1_lJ->fill(fatjets[0].mass(), weight);
        ptl1_lJ->fill(leps[0].pt(), weight);
        mlJ_lJ->fill((fatjets[0].mom() + leps[0].mom()).mass());
      } else if (leps.size() == 1 && fatjets.size() >= 2) {
        njets_lJJ->fill(jets.size(), weight);
        nfatjets_lJJ->fill(fatjets.size(), weight);
        ptJ1_lJJ->fill(fatjets[0].pt(), weight);
        ptJ2_lJJ->fill(fatjets[1].pt(), weight);
        mJ1_lJJ->fill(fatjets[0].mass(), weight);
        mJ2_lJJ->fill(fatjets[1].mass(), weight);
        ptl1_lJJ->fill(leps[0].pt(), weight);
        mJJ_lJJ->fill((fatjets[0].mom() + fatjets[1].mom()).mass());
      } else if (leps.size() == 2 && fatjets.size() >= 1 && leps[0].threeCharge()*leps[1].threeCharge() > 0 ) {
        njets_ssJ->fill(jets.size(), weight);
        nfatjets_ssJ->fill(fatjets.size(), weight);
        ptJ1_ssJ->fill(fatjets[0].pt(), weight);
        mJ1_ssJ->fill(fatjets[0].mass(), weight);
        ptl1_ssJ->fill(leps[0].pt(), weight);
        ptl2_ssJ->fill(leps[1].pt(), weight);
        mlJ_ssJ->fill((fatjets[0].mom() + leps[0].mom()).mass(), weight);
      }

      return;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      for (Histo1DPtr& h : allHists)
        scale(h, crossSection()/sumOfWeights());

      return;
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr njets, nfatjets, nleps
      , njets_lJ, nfatjets_lJ
      , njets_lJJ, nfatjets_lJJ
      , njets_ssJ, nfatjets_ssJ
      , ptJ1_lJ, ptJ1_lJJ, ptJ2_lJJ, ptJ1_ssJ
      , ptl1_lJ, ptl1_lJJ, ptl1_ssJ, ptl2_ssJ
      , mJ1_lJ, mJ1_lJJ, mJ2_lJJ, mJ1_ssJ 
      , mlJ_lJ, mJJ_lJJ, mlJ_ssJ;

    vector<Histo1DPtr> allHists;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TTTT);


}
