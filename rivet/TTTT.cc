// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"

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
      ChargedLeptons cls(Cuts::abseta < 2.5 && Cuts::pT > 25*GeV);
      declare(PromptFinalState(cls, true), "PromptLeptons");

      nleps = bookHisto1D("nleps", 5, -0.5, 4.5, "number of leptons in acceptance");
      nss = bookCounter("nss");

    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& leps = apply<PromptFinalState>(event, "PromptLeptons").particles();

      nleps->fill(leps.size(), event.weight());

      // check for SS dilepton
      if (leps.size() == 2 && leps[0].threeCharge()*leps[1].threeCharge() > 0)
        nss->fill(event.weight());
    }


    /// Normalise histograms etc., after the run
    void finalize() {

      normalize(nleps); // normalize to unity

    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr nleps;
    CounterPtr nss;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TTTT);


}
