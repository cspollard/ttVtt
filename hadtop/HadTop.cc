// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "Rivet/Projections/FastJets.hh"

namespace Rivet {

  /// @brief Add a short analysis description here
  class HadTop : public Analysis {
  public:

    /// Constructor
    DEFAULT_RIVET_ANALYSIS_CTOR(HadTop);


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

      // Book histograms
      hPDF = bookHisto2D("ttPDF", 50, 0, 400, 5, 0, 5, "ttPDF", "mass", "njets", "probability");
      hPDF0b = bookHisto2D("ttPDF0b", 50, 0, 400, 5, 0, 5, "ttPDF0b", "mass", "njets", "probability");
      hPDF1b = bookHisto2D("ttPDF1b", 50, 0, 400, 5, 0, 5, "ttPDF1b", "mass", "njets", "probability");
      hPDF2j0b = bookHisto1D("ttPDF2j0b", 50, 0, 400, "ttPDF2j0b", "mass", "probability");
      hPDF2j1b = bookHisto1D("ttPDF2j1b", 50, 0, 400, "ttPDF2j1b", "mass", "probability");
      hPDF3j0b = bookHisto1D("ttPDF3j0b", 50, 0, 400, "ttPDF3j0b", "mass", "probability");
      hPDF3j1b = bookHisto1D("ttPDF3j1b", 50, 0, 400, "ttPDF3j1b", "mass", "probability");
      hPDF4j0b = bookHisto1D("ttPDF4j0b", 50, 0, 400, "ttPDF4j0b", "mass", "probability");
      hPDF4j1b = bookHisto1D("ttPDF4j1b", 50, 0, 400, "ttPDF4j1b", "mass", "probability");

      hTopPtEta = bookHisto2D("TopPtEta", 50, 0, 500*GeV, 50, 0, 5, "TopPtEta", "pt", "eta", "probability");

      return;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {
      for (const Particle& t : event.allParticles()) {
        if (t.abspid() != 6)
          continue;
        else {
          hTopPtEta->fill(t.pt(), t.abseta(), event.weight());
          break;
        }
      }



      const Jets& jets =
        apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);

      size_t njets = jets.size();

      size_t nbjets = 0;
      for (const Jet& j : jets)
        if (j.bTagged(Cuts::pT > 5*GeV))
          nbjets++;


      FourMomentum alljets;
      for (const Jet& j : jets)
        alljets += j;

      double mass = alljets.mass();
      hPDF->fill(mass, njets, event.weight());

      if (nbjets == 0) {
        hPDF0b->fill(mass, njets, event.weight());
        if (njets == 2)
          hPDF2j0b->fill(mass, event.weight());
        else if (njets == 3)
          hPDF3j0b->fill(mass, event.weight());
        else if (njets == 4)
          hPDF4j0b->fill(mass, event.weight());
      } else if (nbjets == 1) {
        hPDF1b->fill(mass, njets, event.weight());
        if (njets == 2)
          hPDF2j1b->fill(mass, event.weight());
        else if (njets == 3)
          hPDF3j1b->fill(mass, event.weight());
        else if (njets == 4)
          hPDF4j1b->fill(mass, event.weight());
      }

      return;
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      vector<Histo2DPtr> h2ds = { hPDF, hPDF0b, hPDF1b, hTopPtEta };
      for (Histo2DPtr& h : h2ds)
        scale(h, 1.0/sumOfWeights());

      vector<Histo1DPtr> h1ds = { hPDF2j0b, hPDF2j1b, hPDF3j0b, hPDF3j1b, hPDF4j0b, hPDF4j1b };
      for (Histo1DPtr& h : h1ds)
        scale(h, 1.0/sumOfWeights());

      return;
    }

    //@}


    /// @name Histograms
    //@{
    Histo2DPtr hPDF, hPDF0b, hPDF1b, hTopPtEta;
    Histo1DPtr hPDF2j0b, hPDF2j1b, hPDF3j0b, hPDF3j1b, hPDF4j0b, hPDF4j1b;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(HadTop);
}
