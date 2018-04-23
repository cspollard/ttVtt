// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"

namespace Rivet {

  // from https://stackoverflow.com/questions/3418231/replace-part-of-a-string-with-another-string
  string replace(const string& str, const string& from, const string& to) {
    string ret = str;
    size_t start_pos = str.find(from);
    if(start_pos == std::string::npos)
      return ret;
    ret.replace(start_pos, from.length(), to);
    return ret;
  }


  const string ptstr = "\\ensuremath{p_\\mathrm{T}}";
  const string ptttstr = "\\ensuremath{p_\\mathrm{T}(tt)}";
  const string mttstr = "\\ensuremath{m_{tt}}";
  const string mstr = "\\ensuremath{m}";
  const string nstr = "\\ensuremath{n}";
  const string dphittstr = "\\ensuremath{\\Delta\\phi(tt)}";

  string dxdy(const string& x, const string& y, const string& xunit, const string& yunit) {
    return "\\ensuremath{\\frac{\\mathrm{d}" + x + "}{\\mathrm{d}" + y + "} \\Big[ \\frac{" + xunit + "}{" + yunit + "} \\Big]}";
  }

  string dsigdy(const string& x, const string& xunit) {
    return dxdy(x, "\\sigma", xunit, "\\mathrm{pb}");
  }

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


      njets = bookH("njets", 21, -0.5, 20.5, "njets", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets = bookH("ncentjets", 21, -0.5, 20.5, "ncentjets", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets = bookH("nfwdjets", 11, -0.5, 10.5, "nfwdjets", "forward jet multiplicity", dsigdy(nstr, "1"));
      nfatjets = bookH("nfatjets", 4, -0.5, 3.5, "nfatjets", "large-$R$ jet multiplicity", dsigdy(nstr, "1"));
      nleps = bookH("nleps", 5, -0.5, 4.5, "nleps", "prompt lepton multiplicity", dsigdy(nstr, "1"));

      njets_JJ = bookH("njets_JJ", 21, -0.5, 20.5, "njets_JJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_JJ = bookH("ncentjets_JJ", 21, -0.5, 20.5, "ncentjets_JJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_JJ = bookH("nfwdjets_JJ", 11, -0.5, 10.5, "nfwdjets_JJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_JJ = bookH("naddjets_JJ", 21, -0.5, 20.5, "naddjets_JJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      nfatjets_JJ = bookH("nfatjets_JJ", 4, -0.5, 3.5, "nfatjets_JJ", "large-$R$ jet multiplicity", dsigdy(nstr, "1"));

      ptth1_JJ = bookH("ptth1_JJ", 25, 0, 2, "ptth1_JJ", "leading hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      ptth2_JJ = bookH("ptth2_JJ", 25, 0, 2, "ptth2_JJ", "subleading top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      mth1_JJ = bookH("mth1_JJ", 25, 0, 500, "mth1_JJ", "leading hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      mth2_JJ = bookH("mth2_JJ", 25, 0, 500, "mth2_JJ", "leading hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      dphitt_JJ = bookH("dphitt_JJ", 20, 0, 4, "dphitt_JJ", dphittstr, dsigdy(dphittstr, "\\mathrm{rad}"));
      pttt_JJ = bookH("pttt_JJ", 25, 0, 1, "pttt_JJ", ptttstr + " [TeV]", dsigdy(ptttstr, "\\mathrm{TeV}"));
      mtt_JJ = bookH("mtt_JJ", 15, 0, 3, "mtt_JJ", "$tt$ invariant mass [TeV]", dsigdy(mttstr, "\\mathrm{TeV}"));

      njets_lJ = bookH("njets_lJ", 21, -0.5, 20.5, "njets_lJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_lJ = bookH("ncentjets_lJ", 21, -0.5, 20.5, "ncentjets_lJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_lJ = bookH("nfwdjets_lJ", 11, -0.5, 10.5, "nfwdjets_lJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_lJ = bookH("naddjets_lJ", 21, -0.5, 20.5, "naddjets_lJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      nfatjets_lJ = bookH("nfatjets_lJ", 4, -0.5, 3.5, "nfatjets_lJ", "large-$R$ jet multiplicity", dsigdy(nstr, "1"));

      ptth_lJ = bookH("ptth_lJ", 25, 0, 2, "ptth_lJ", "hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      pttl_lJ = bookH("pttl_lJ", 25, 0, 2, "pttl_lJ", "leptonic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      ptl1_lJ = bookH("ptl1_lJ", 25, 0, 1e3, "ptl1_lJ", "leading lepton $p_\\mathrm{T}$ [GeV]", dsigdy(ptstr, "\\mathrm{GeV}"));
      mth_lJ = bookH("mth_lJ", 25, 0, 500, "mth_lJ", "hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      mtl_lJ = bookH("mtl_lJ", 25, 0, 500, "mtl_lJ", "leptonic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      dphitt_lJ = bookH("dphitt_lJ", 20, 0, 4, "dphitt_lJ", dphittstr, dsigdy(dphittstr, "\\mathrm{rad}"));
      pttt_lJ = bookH("pttt_lJ", 25, 0, 1, "pttt_lJ", ptttstr + " [TeV]", dsigdy(ptttstr, "\\mathrm{TeV}"));
      mtt_lJ = bookH("mtt_lJ", 15, 0, 3, "mtt_lJ", "$tt$ invariant mass [TeV]", dsigdy(mttstr, "\\mathrm{TeV}"));

      njets_lJJ = bookH("njets_lJJ", 21, -0.5, 20.5, "njets_lJJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_lJJ = bookH("ncentjets_lJJ", 21, -0.5, 20.5, "ncentjets_lJJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_lJJ = bookH("nfwdjets_lJJ", 11, -0.5, 10.5, "nfwdjets_lJJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_lJJ = bookH("naddjets_lJJ", 21, -0.5, 20.5, "naddjets_lJJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      nfatjets_lJJ = bookH("nfatjets_lJJ", 4, -0.5, 3.5, "nfatjets_lJJ", "large-$R$ jet multiplicity", dsigdy(nstr, "1"));

      ptl1_lJJ = bookH("ptl1_lJJ", 25, 0, 1e3, "ptl1_lJJ", "leading lepton $p_\\mathrm{T}$ [GeV]", dsigdy(ptstr, "\\mathrm{GeV}"));
      ptth1_lJJ = bookH("ptth1_lJJ", 25, 0, 2, "ptth1_lJJ", "leading hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      ptth2_lJJ = bookH("ptth2_lJJ", 25, 0, 2, "ptth2_lJJ", "subleading hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      mth1_lJJ = bookH("mth1_lJJ", 25, 0, 500, "mth1_lJJ", "leading hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      mth2_lJJ = bookH("mth2_lJJ", 25, 0, 500, "mth2_lJJ", "subleading hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      dphitt_lJJ = bookH("dphitt_lJJ", 20, 0, 4, "dphitt_lJJ", dphittstr, dsigdy(dphittstr, "\\mathrm{rad}"));
      pttt_lJJ = bookH("pttt_lJJ", 25, 0, 1, "pttt_lJJ", ptttstr + " [TeV]", dsigdy(ptttstr, "\\mathrm{TeV}"));
      mtt_lJJ = bookH("mtt_lJJ", 15, 0, 3, "mtt_lJJ", "$tt$ invariant mass [TeV]", dsigdy(mttstr, "\\mathrm{TeV}"));

      njets_ssJ = bookH("njets_ssJ", 21, -0.5, 20.5, "njets_ssJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_ssJ = bookH("ncentjets_ssJ", 21, -0.5, 20.5, "ncentjets_ssJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_ssJ = bookH("nfwdjets_ssJ", 11, -0.5, 10.5, "nfwdjets_ssJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_ssJ = bookH("naddjets_lJ", 21, -0.5, 20.5, "naddjets_lJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      nfatjets_ssJ = bookH("nfatjets_ssJ", 4, -0.5, 3.5, "nfatjets_ssJ", "large-$R$ jet multiplicity", dsigdy(nstr, "1"));

      ptth_ssJ = bookH("ptth_ssJ", 25, 0, 2, "ptth_ssJ", "hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      pttl_ssJ = bookH("pttl_ssJ", 25, 0, 2, "pttl_ssJ", "leptonic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      ptl1_ssJ = bookH("ptl1_ssJ", 25, 0, 1e3, "ptl1_lJJ", "leading lepton $p_\\mathrm{T}$ [GeV]", dsigdy(ptstr, "\\mathrm{GeV}"));
      ptl2_ssJ = bookH("ptl2_ssJ", 25, 0, 1e3, "ptl2_ssJ", "subleading lepton $p_\\mathrm{T}$ [GeV]", dsigdy(ptstr, "\\mathrm{GeV}"));
      mth_ssJ = bookH("mth_ssJ", 25, 0, 500, "mth_ssJ", "hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      mtl_ssJ = bookH("mtl_ssJ", 25, 0, 500, "mtl_ssJ", "leptonic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      dphitt_ssJ = bookH("dphitt_ssJ", 20, 0, 4, "dphitt_ssJ", dphittstr, dsigdy(dphittstr, "\\mathrm{rad}"));
      pttt_ssJ = bookH("pttt_ssJ", 25, 0, 1, "pttt_ssJ", ptttstr + " [TeV]", dsigdy(ptttstr, "\\mathrm{TeV}"));
      mtt_ssJ = bookH("mtt_ssJ", 15, 0, 3, "mtt_ssJ", "$tt$ invariant mass [TeV]", dsigdy(mttstr, "\\mathrm{TeV}"));

    }

    Histo1DPtr bookH(const string& path, double nb, double bmin, double bmax
        , const string& title, const string& xlabel, const string& ylabel) {
      Histo1DPtr h = bookHisto1D(path, nb, bmin, bmax, title, xlabel, ylabel);
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
      const Jets& centjets =
        apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      const Jets& fwdjets =
        apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta > 2.5);
      const Jets& fatjets =
        apply<FastJets>(event, "FatJets").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.0 && Cuts::mass > 100*GeV);

      double weight = event.weight();

      nleps->fill(leps.size(), weight);
      njets->fill(jets.size(), weight);
      ncentjets->fill(centjets.size(), weight);
      nfwdjets->fill(fwdjets.size(), weight);
      nfatjets->fill(fatjets.size(), weight);


      if (leps.size() == 0 && fatjets.size() >= 2) {
        njets_JJ->fill(jets.size(), weight);
        ncentjets_JJ->fill(centjets.size(), weight);
        nfwdjets_JJ->fill(fwdjets.size(), weight);
        nfatjets_JJ->fill(fatjets.size(), weight);

        ptth1_JJ->fill(fatjets[0].pt()/TeV, weight);
        ptth2_JJ->fill(fatjets[1].pt()/TeV, weight);
        mth1_JJ->fill(fatjets[0].mass()/GeV, weight);
        mth2_JJ->fill(fatjets[1].mass()/GeV, weight);

        Jets goodfatjets;
        goodfatjets.push_back(fatjets[0]);
        goodfatjets.push_back(fatjets[1]);
        vector<const Jet*> goodjets = additionalJets(jets, goodfatjets);
        naddjets_JJ->fill(goodjets.size(), weight);

        FourMomentum t1 = goodfatjets[0].mom();
        FourMomentum t2 = goodfatjets[1].mom();
        FourMomentum tt = t1 + t2;

        mth1_JJ->fill(t1.mass()/GeV, weight);
        mth2_JJ->fill(t2.mass()/GeV, weight);

        dphitt_JJ->fill(abs(deltaPhi(t1, t2)), weight);

        pttt_JJ->fill(tt.pt()/TeV, weight);
        mtt_JJ->fill(tt.mass()/TeV, weight);

      } if (leps.size() == 1 && fatjets.size() == 1) {
        njets_lJ->fill(jets.size(), weight);
        ncentjets_lJ->fill(centjets.size(), weight);
        nfwdjets_lJ->fill(fwdjets.size(), weight);
        nfatjets_lJ->fill(fatjets.size(), weight);
        ptl1_lJ->fill(leps[0].pt()/GeV, weight);

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

          ptth_lJ->fill(th.pt()/TeV, weight);
          pttl_lJ->fill(tl.pt()/TeV, weight);

          mtl_lJ->fill(tl.mass()/GeV, weight);
          mth_lJ->fill(th.mass()/GeV, weight);

          dphitt_lJ->fill(abs(deltaPhi(tl, th)), weight);

          pttt_lJ->fill(tt.pt()/TeV, weight);
          mtt_lJ->fill(tt.mass()/TeV, weight);
        }

      } else if (leps.size() == 1 && fatjets.size() >= 2) {
        njets_lJJ->fill(jets.size(), weight);
        ncentjets_lJJ->fill(centjets.size(), weight);
        nfwdjets_lJJ->fill(fwdjets.size(), weight);
        nfatjets_lJJ->fill(fatjets.size(), weight);

        ptth1_lJJ->fill(fatjets[0].pt()/TeV, weight);
        ptth2_lJJ->fill(fatjets[1].pt()/TeV, weight);
        mth1_lJJ->fill(fatjets[0].mass()/GeV, weight);
        mth2_lJJ->fill(fatjets[1].mass()/GeV, weight);
        ptl1_lJJ->fill(leps[0].pt()/GeV, weight);

        Jets goodfatjets;
        goodfatjets.push_back(fatjets[0]);
        goodfatjets.push_back(fatjets[1]);
        vector<const Jet*> goodjets = additionalJets(jets, goodfatjets);
        naddjets_lJJ->fill(goodjets.size(), weight);

        FourMomentum t1 = goodfatjets[0].mom();
        FourMomentum t2 = goodfatjets[1].mom();
        FourMomentum tt = t1 + t2;

        mth1_lJJ->fill(t1.mass()/GeV, weight);
        mth2_lJJ->fill(t2.mass()/GeV, weight);

        dphitt_lJJ->fill(abs(deltaPhi(t1, t2)), weight);

        pttt_lJJ->fill(tt.pt()/TeV, weight);
        mtt_lJJ->fill(tt.mass()/TeV, weight);

      } else if (leps.size() == 2 && fatjets.size() >= 1 && leps[0].threeCharge()*leps[1].threeCharge() > 0 ) {
        njets_ssJ->fill(jets.size(), weight);
        ncentjets_ssJ->fill(centjets.size(), weight);
        nfwdjets_ssJ->fill(fwdjets.size(), weight);
        nfatjets_ssJ->fill(fatjets.size(), weight);
        ptl1_ssJ->fill(leps[0].pt()/GeV, weight);
        ptl2_ssJ->fill(leps[1].pt()/GeV, weight);

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

          pttl_ssJ->fill(tl.pt()/TeV, weight);
          ptth_ssJ->fill(th.pt()/TeV, weight);

          mtl_ssJ->fill(tl.mass()/GeV, weight);
          mth_ssJ->fill(th.mass()/GeV, weight);

          dphitt_ssJ->fill(abs(deltaPhi(tl, th)), weight);

          pttt_ssJ->fill(tt.pt()/TeV, weight);
          mtt_ssJ->fill(tt.mass()/TeV, weight);
        }

      }

      return;
    }


    /// Normalise histograms etc., after the run
    void finalize() {
      const string onebysig = "\\ensuremath{\\frac{1}{\\sigma}}";


      for (Histo1DPtr& h : allHists) {
        Histo1DPtr hnorm = make_shared<Histo1D>(*h);
        normalize(hnorm);
        hnorm->setPath(hnorm->path() + "_norm");
        hnorm->setAnnotation("YLabel", onebysig + replace(hnorm->annotation("YLabel"), "pb", "1"));
        addAnalysisObject(hnorm);

        scale(h, crossSection()/sumOfWeights());
      }

      return;
    }

    //@}


    /// @name Histograms
    //@{
    Histo1DPtr njets;
    Histo1DPtr ncentjets;
    Histo1DPtr nfwdjets;
    Histo1DPtr nfatjets;
    Histo1DPtr nleps;

    Histo1DPtr njets_JJ;
    Histo1DPtr ncentjets_JJ;
    Histo1DPtr nfwdjets_JJ;
    Histo1DPtr naddjets_JJ;
    Histo1DPtr nfatjets_JJ;

    Histo1DPtr ptth1_JJ;
    Histo1DPtr ptth2_JJ;
    Histo1DPtr mth1_JJ;
    Histo1DPtr mth2_JJ;
    Histo1DPtr dphitt_JJ;
    Histo1DPtr pttt_JJ;
    Histo1DPtr mtt_JJ;

    Histo1DPtr njets_lJ;
    Histo1DPtr ncentjets_lJ;
    Histo1DPtr nfwdjets_lJ;
    Histo1DPtr naddjets_lJ;
    Histo1DPtr nfatjets_lJ;

    Histo1DPtr ptth_lJ;
    Histo1DPtr pttl_lJ;
    Histo1DPtr ptl1_lJ;
    Histo1DPtr mth_lJ;
    Histo1DPtr mtl_lJ;
    Histo1DPtr dphitt_lJ;
    Histo1DPtr pttt_lJ;
    Histo1DPtr mtt_lJ;

    Histo1DPtr njets_lJJ;
    Histo1DPtr ncentjets_lJJ;
    Histo1DPtr nfwdjets_lJJ;
    Histo1DPtr naddjets_lJJ;
    Histo1DPtr nfatjets_lJJ;

    Histo1DPtr ptth1_lJJ;
    Histo1DPtr ptth2_lJJ;
    Histo1DPtr ptl1_lJJ;
    Histo1DPtr mth1_lJJ;
    Histo1DPtr mth2_lJJ;
    Histo1DPtr dphitt_lJJ;
    Histo1DPtr pttt_lJJ;
    Histo1DPtr mtt_lJJ;

    Histo1DPtr njets_ssJ;
    Histo1DPtr ncentjets_ssJ;
    Histo1DPtr nfwdjets_ssJ;
    Histo1DPtr naddjets_ssJ;
    Histo1DPtr nfatjets_ssJ;

    Histo1DPtr ptth_ssJ;
    Histo1DPtr pttl_ssJ;
    Histo1DPtr ptl1_ssJ;
    Histo1DPtr ptl2_ssJ;
    Histo1DPtr mth_ssJ;
    Histo1DPtr mtl_ssJ;
    Histo1DPtr dphitt_ssJ;
    Histo1DPtr pttt_ssJ;
    Histo1DPtr mtt_ssJ;

    vector<Histo1DPtr> allHists;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TTTT);


}
