// -*- C++ -*-
#include "Rivet/Analysis.hh"
#include "Rivet/Projections/FinalState.hh"
#include "Rivet/Projections/PromptFinalState.hh"
#include "Rivet/Projections/ChargedLeptons.hh"
#include "Rivet/Projections/FastJets.hh"
#include "Rivet/Projections/VetoedFinalState.hh"
#include "YODA/ReaderYODA.h"

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
  const string chi2str = "\\ensuremath{\\chi^2}";
  const string logttprobstr = "\\ensuremath{log(\\text{Prob}(tt)})";
  const string ptttstr = "\\ensuremath{p_\\mathrm{T}(tt)}";
  const string mttstr = "\\ensuremath{m_{tt}}";
  const string mstr = "\\ensuremath{m}";
  const string nstr = "\\ensuremath{n}";
  const string dphittstr = "\\ensuremath{\\Delta\\phi(tt)}";

  string dxdy(const string& x, const string& y, const string& xunit, const string& yunit) {
    return "\\ensuremath{\\frac{\\mathrm{d}" + x + "}{\\mathrm{d}" + y + "} \\Big[ \\frac{" + xunit + "}{" + yunit + "} \\Big]}";
  }

  string dsigdy(const string& y, const string& yunit) {
    return dxdy("\\sigma", y, "\\mathrm{pb}", yunit);
  }

  bool cmp_pt(const Jet& j1, const Jet& j2) {
    return j1.pt() > j2.pt();
  }

  template<class T>
  pair<vector<T>, vector<T>> splitAt(const vector<T>& v, size_t n) {
    if (v.size() < n)
      return make_pair(v, vector<T>());
    else
      return make_pair(vector<T>(v.begin(), v.begin()+n), vector<T>(v.begin()+n, v.end()));
  }


  double mw = 80.51*GeV;
  double mt_mw = 85.17*GeV;
  double sigw = 12.07*GeV;
  double sig_mt_mw = 16.05*GeV;

  double chi2_hadhad(const Jets& jets) {
    double minchi2 = 1e9;
    if (jets.size() < 6)
      return minchi2;

    // TODO
    // only look at up to the first eight jets for now...
    Jets js;
    size_t mx = jets.size();
    if (mx > 8)
      mx = 8;

    for (size_t i = 0; i < mx; i++)
      js.push_back(jets[i]);

    // TODO
    // this will be extremely inefficient for large jet multiplicities
    // because we will double-check some permutations.
    // oh well.
    size_t count = 0;
    do {
      double mjj1 = (js[0].mom() + js[1].mom()).mass();
      double mt_w1 = (js[0].mom() + js[1].mom() + js[2].mom()).mass();
      double mjj2 = (js[3].mom() + js[4].mom()).mass();
      double mt_w2 = (js[3].mom() + js[4].mom() + js[5].mom()).mass();

      double w1term = (mjj1 - mw) / sigw;
      double t1term = (mt_w1 - mt_mw) / sig_mt_mw;
      double w2term = (mjj2 - mw) / sigw;
      double t2term = (mt_w2 - mt_mw) / sig_mt_mw;

      double chi2 = w1term*w1term + t1term*t1term + w2term*w2term + t2term*t2term;

      if (chi2 < minchi2)
        minchi2 = chi2;

      count++;
    } while (next_permutation(js.begin(), js.end(), cmp_pt));

    return minchi2;
  }

  double topProb(const Histo2D& topPDF0b, const Histo2D& topPDF1b, const Jets& jets) {
    size_t nj = jets.size();
    if (nj != 2 && nj != 3)
      return 0.0;

    size_t nb = 0;
    FourMomentum alljets;
    for (const Jet& jet : jets) {
      if (jet.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5))
        nb++;

      alljets += jet.mom();
    }

    double mass = alljets.mass();
    if (nb > 1 || mass > 400)
      return 0.0;

    // cout << "mass: " << mass << endl;
    // cout << "njets: " << nj << endl;
    // cout << "nbjets: " << nb << endl;

    if (nb)
      return topPDF1b.binAt(alljets.mass(), nj + 0.1).volume();
    else
      return topPDF0b.binAt(alljets.mass(), nj + 0.1).volume();
  }

  double ttProb(const Histo2D& topPDF0b, const Histo2D& topPDF1b, const Jets& jets) {
    size_t nj = jets.size();

    if (nj < 4)
      return 0.0;

    // only look at the first 8 jets for now.
    Jets js;
    for (size_t i = 0; i < jets.size() && i < 8; i++)
      js.push_back(jets[i]);

    double bestprob = 1e-50;
    do {
      double prob;
      pair<Jets, Jets> top1_rest;
      pair<Jets, Jets> top2_rest;
 
        
      top1_rest = splitAt(js, 2);
      top2_rest = splitAt(top1_rest.second, 2);

      // cout << "top1 length: " << top1_rest.first.size() << endl;
      // cout << "top2 length: " << top2_rest.first.size() << endl;

      prob = topProb(topPDF0b, topPDF1b, top1_rest.first) * topProb(topPDF0b, topPDF1b, top2_rest.first);

      if (prob > bestprob)
        bestprob = prob;


      if (nj < 5)
        continue;

      top1_rest = splitAt(js, 2);
      top2_rest = splitAt(top1_rest.second, 3);

      // cout << "top1 length: " << top1_rest.first.size() << endl;
      // cout << "top2 length: " << top2_rest.first.size() << endl;

      prob = topProb(topPDF0b, topPDF1b, top1_rest.first) * topProb(topPDF0b, topPDF1b, top2_rest.first);

      if (prob > bestprob)
        bestprob = prob;


      top1_rest = splitAt(js, 3);
      top2_rest = splitAt(top1_rest.second, 2);

      // cout << "top1 length: " << top1_rest.first.size() << endl;
      // cout << "top2 length: " << top2_rest.first.size() << endl;

      prob = topProb(topPDF0b, topPDF1b, top1_rest.first) * topProb(topPDF0b, topPDF1b, top2_rest.first);

      if (prob > bestprob)
        bestprob = prob;


      if (nj < 6)
        continue;

      top1_rest = splitAt(js, 3);
      top2_rest = splitAt(top1_rest.second, 3);

      // cout << "top1 length: " << top1_rest.first.size() << endl;
      // cout << "top2 length: " << top2_rest.first.size() << endl;

      prob = topProb(topPDF0b, topPDF1b, top1_rest.first) * topProb(topPDF0b, topPDF1b, top2_rest.first);

      if (prob > bestprob)
        bestprob = prob;

    } while (next_permutation(js.begin(), js.end(), cmp_pt));

    cout << "bestprob: " << bestprob << endl;
    return bestprob;
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

      // read in control histograms for top tagging

      YODA::Reader& r = YODA::ReaderYODA::create();

      vector<AnalysisObject*> inputHists = r.read("toptemplate.yoda");

      for (AnalysisObject* aoptr : inputHists) {
        if (aoptr->path() == "/HadTop/ttPDF0b") {
          topPDF0b = * ((Histo2D*) aoptr);
          cout << "found /HadTop/ttPDF0b" << endl;
        } else if (aoptr->path() == "/HadTop/ttPDF1b") {
          topPDF1b = * ((Histo2D*) aoptr);
          cout << "found /HadTop/ttPDF1b" << endl;
        }

        // clean up after ourselves.
        delete aoptr;
        continue;
      }

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
      nbjets = bookH("nbjets", 21, -0.5, 20.5, "nbjets", "$b$-jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets = bookH("nfwdjets", 11, -0.5, 10.5, "nfwdjets", "forward jet multiplicity", dsigdy(nstr, "1"));
      ntopjets = bookH("ntopjets", 4, -0.5, 3.5, "ntopjets", "top-tagged jet multiplicity", dsigdy(nstr, "1"));
      ntopbjets = bookH("ntopbjets", 4, -0.5, 3.5, "ntopbjets", "$b$ + top-tagged jet multiplicity", dsigdy(nstr, "1"));
      nleps = bookH("nleps", 5, -0.5, 4.5, "nleps", "prompt lepton multiplicity", dsigdy(nstr, "1"));

      njets_JJ = bookH("njets_JJ", 21, -0.5, 20.5, "njets_JJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_JJ = bookH("ncentjets_JJ", 21, -0.5, 20.5, "ncentjets_JJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_JJ = bookH("nfwdjets_JJ", 11, -0.5, 10.5, "nfwdjets_JJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_JJ = bookH("naddjets_JJ", 21, -0.5, 20.5, "naddjets_JJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      ntopjets_JJ = bookH("ntopjets_JJ", 4, -0.5, 3.5, "ntopjets_JJ", "top-tagged jet multiplicity", dsigdy(nstr, "1"));
      naddbjets_JJ = bookH("naddbjets_JJ", 21, -0.5, 20.5, "naddbjets_JJ", "additional $b$-jet multiplicity", dsigdy(nstr, "1"));
      naddljets_JJ = bookH("naddljets_JJ", 21, -0.5, 20.5, "naddljets_JJ", "additional light-jet multiplicity", dsigdy(nstr, "1"));
      ntopbjets_JJ = bookH("ntopbjets_JJ", 4, -0.5, 3.5, "ntopbjets_JJ", "$b$ + top-tagged jet multiplicity", dsigdy(nstr, "1"));

      ptth1_JJ = bookH("ptth1_JJ", 25, 0, 2, "ptth1_JJ", "leading hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      ptth2_JJ = bookH("ptth2_JJ", 25, 0, 2, "ptth2_JJ", "subleading top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      mth1_JJ = bookH("mth1_JJ", 25, 0, 500, "mth1_JJ", "leading hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      mth2_JJ = bookH("mth2_JJ", 25, 0, 500, "mth2_JJ", "leading hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      dphitt_JJ = bookH("dphitt_JJ", 20, 0, 4, "dphitt_JJ", dphittstr, dsigdy(dphittstr, "\\mathrm{rad}"));
      pttt_JJ = bookH("pttt_JJ", 25, 0, 1, "pttt_JJ", ptttstr + " [TeV]", dsigdy(ptttstr, "\\mathrm{TeV}"));
      mtt_JJ = bookH("mtt_JJ", 15, 0, 3, "mtt_JJ", "$tt$ invariant mass [TeV]", dsigdy(mttstr, "\\mathrm{TeV}"));
      chi2_JJ = bookH("chi2_JJ", 20, 0, 1000, "chi2_JJ", chi2str, dsigdy(chi2str, "1"));
      logttprob_JJ = bookH("logttprob_JJ", 20, -20, 0, "logttprob_JJ", logttprobstr, dsigdy(logttprobstr, "1"));

      njets_lJ = bookH("njets_lJ", 21, -0.5, 20.5, "njets_lJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_lJ = bookH("ncentjets_lJ", 21, -0.5, 20.5, "ncentjets_lJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_lJ = bookH("nfwdjets_lJ", 11, -0.5, 10.5, "nfwdjets_lJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_lJ = bookH("naddjets_lJ", 21, -0.5, 20.5, "naddjets_lJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      ntopjets_lJ = bookH("ntopjets_lJ", 4, -0.5, 3.5, "ntopjets_lJ", "top-tagged jet multiplicity", dsigdy(nstr, "1"));
      naddbjets_lJ = bookH("naddbjets_lJ", 21, -0.5, 20.5, "naddbjets_lJ", "additional $b$-jet multiplicity", dsigdy(nstr, "1"));
      naddljets_lJ = bookH("naddljets_lJ", 21, -0.5, 20.5, "naddljets_lJ", "additional light-jet multiplicity", dsigdy(nstr, "1"));
      ntopbjets_lJ = bookH("ntopbjets_lJ", 4, -0.5, 3.5, "ntopbjets_lJ", "$b$ + top-tagged jet multiplicity", dsigdy(nstr, "1"));

      ptth_lJ = bookH("ptth_lJ", 25, 0, 2, "ptth_lJ", "hadronic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      pttl_lJ = bookH("pttl_lJ", 25, 0, 2, "pttl_lJ", "leptonic top $p_\\mathrm{T}$ [TeV]", dsigdy(ptstr, "\\mathrm{TeV}"));
      ptl1_lJ = bookH("ptl1_lJ", 25, 0, 1e3, "ptl1_lJ", "leading lepton $p_\\mathrm{T}$ [GeV]", dsigdy(ptstr, "\\mathrm{GeV}"));
      mth_lJ = bookH("mth_lJ", 25, 0, 500, "mth_lJ", "hadronic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      mtl_lJ = bookH("mtl_lJ", 25, 0, 500, "mtl_lJ", "leptonic top mass [GeV]", dsigdy(mstr, "\\mathrm{GeV}"));
      dphitt_lJ = bookH("dphitt_lJ", 20, 0, 4, "dphitt_lJ", dphittstr, dsigdy(dphittstr, "\\mathrm{rad}"));
      pttt_lJ = bookH("pttt_lJ", 25, 0, 1, "pttt_lJ", ptttstr + " [TeV]", dsigdy(ptttstr, "\\mathrm{TeV}"));
      mtt_lJ = bookH("mtt_lJ", 15, 0, 3, "mtt_lJ", "$tt$ invariant mass [TeV]", dsigdy(mttstr, "\\mathrm{TeV}"));
      chi2_lJ = bookH("chi2_lJ", 20, 0, 1000, "chi2_lJ", chi2str, dsigdy(chi2str, "1"));
      logttprob_lJ = bookH("logttprob_lJ", 20, -20, 0, "logttprob_lJ", logttprobstr, dsigdy(logttprobstr, "1"));

      njets_lJJ = bookH("njets_lJJ", 21, -0.5, 20.5, "njets_lJJ", "jet multiplicity", dsigdy(nstr, "1"));
      ncentjets_lJJ = bookH("ncentjets_lJJ", 21, -0.5, 20.5, "ncentjets_lJJ", "central jet multiplicity", dsigdy(nstr, "1"));
      nfwdjets_lJJ = bookH("nfwdjets_lJJ", 11, -0.5, 10.5, "nfwdjets_lJJ", "forward jet multiplicity", dsigdy(nstr, "1"));
      naddjets_lJJ = bookH("naddjets_lJJ", 21, -0.5, 20.5, "naddjets_lJJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      ntopjets_lJJ = bookH("ntopjets_lJJ", 4, -0.5, 3.5, "ntopjets_lJJ", "top-tagged jet multiplicity", dsigdy(nstr, "1"));
      naddbjets_lJJ = bookH("naddbjets_lJJ", 21, -0.5, 20.5, "naddbjets_lJJ", "additional $b$-jet multiplicity", dsigdy(nstr, "1"));
      naddljets_lJJ = bookH("naddljets_lJJ", 21, -0.5, 20.5, "naddljets_lJJ", "additional light-jet multiplicity", dsigdy(nstr, "1"));
      ntopbjets_lJJ = bookH("ntopbjets_lJJ", 4, -0.5, 3.5, "ntopbjets_lJJ", "$b$ + top-tagged jet multiplicity", dsigdy(nstr, "1"));

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
      naddjets_ssJ = bookH("naddjets_ssJ", 21, -0.5, 20.5, "naddjets_ssJ", "additional jet multiplicity", dsigdy(nstr, "1"));
      ntopjets_ssJ = bookH("ntopjets_ssJ", 4, -0.5, 3.5, "ntopjets_ssJ", "top-tagged jet multiplicity", dsigdy(nstr, "1"));
      naddbjets_ssJ = bookH("naddbjets_ssJ", 21, -0.5, 20.5, "naddbjets_ssJ", "additional $b$-jet multiplicity", dsigdy(nstr, "1"));
      naddljets_ssJ = bookH("naddljets_ssJ", 21, -0.5, 20.5, "naddljets_ssJ", "additional light-jet multiplicity", dsigdy(nstr, "1"));
      ntopbjets_ssJ = bookH("ntopbjets_ssJ", 4, -0.5, 3.5, "ntopbjets_ssJ", "$b$ + top-tagged jet multiplicity", dsigdy(nstr, "1"));

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

    Jets additionalJets(const Jets& jets, const Jets& topjets) {
      Jets addjets;

      for (const Jet& j : jets) {
        const FourMomentum mom = j.mom();
        bool pass = true;

        // soft jets should not be removed, since they are likely to
        // come from the spectators.
        if (j.pt() < 60*GeV) {
          addjets.push_back(j);
          continue;
        }

        for (const Jet& fj : topjets) {
          if (deltaR(fj.mom(), mom) > 1.2)
            continue;

          pass = false;
          break;
        }

        if (pass)
          addjets.push_back(j);
      }

      return addjets;
    }


    Jets btaggedJets(const Jets& jets) {
      Jets bjets;
      for (const Jet& j : jets) {
        // only can tag b-hadrons in the tracker fiducial volume.
        if (j.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5))
          bjets.push_back(j);
      }

      return bjets;
    }


    /// Perform the per-event analysis
    void analyze(const Event& event) {

      const Particles& leps = apply<PromptFinalState>(event, "PromptLeptons").particles();
      const Jets& jets = apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV);
      const Jets& centjets =
        apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta < 2.5);
      const Jets& fwdjets =
        apply<FastJets>(event, "Jets").jetsByPt(Cuts::pT > 25*GeV && Cuts::abseta > 2.5);
      const Jets& topjets =
        apply<FastJets>(event, "FatJets").jetsByPt(Cuts::pT > 300*GeV && Cuts::abseta < 2.0 && Cuts::mass > 100*GeV);

      double weight = event.weight();

      nleps->fill(leps.size(), weight);
      njets->fill(jets.size(), weight);
      nbjets->fill(btaggedJets(jets).size(), weight);
      ncentjets->fill(centjets.size(), weight);
      nfwdjets->fill(fwdjets.size(), weight);
      ntopjets->fill(topjets.size(), weight);
      ntopbjets->fill(btaggedJets(topjets).size(), weight);


      if (leps.size() == 0 && topjets.size() >= 2) {
        njets_JJ->fill(jets.size(), weight);
        ncentjets_JJ->fill(centjets.size(), weight);
        nfwdjets_JJ->fill(fwdjets.size(), weight);
        ntopjets_JJ->fill(topjets.size(), weight);

        ptth1_JJ->fill(topjets[0].pt()/TeV, weight);
        ptth2_JJ->fill(topjets[1].pt()/TeV, weight);
        mth1_JJ->fill(topjets[0].mass()/GeV, weight);
        mth2_JJ->fill(topjets[1].mass()/GeV, weight);

        Jets goodtopjets;
        goodtopjets.push_back(topjets[0]);
        goodtopjets.push_back(topjets[1]);
        ntopbjets_JJ->fill(btaggedJets(goodtopjets).size(), weight);

        Jets addjets = additionalJets(jets, goodtopjets);
        naddjets_JJ->fill(addjets.size(), weight);
        Jets addbjets = btaggedJets(addjets);
        naddbjets_JJ->fill(addbjets.size(), weight);
        naddljets_JJ->fill(addjets.size()-addbjets.size(), weight);

        FourMomentum t1 = goodtopjets[0].mom();
        FourMomentum t2 = goodtopjets[1].mom();
        FourMomentum tt = t1 + t2;

        mth1_JJ->fill(t1.mass()/GeV, weight);
        mth2_JJ->fill(t2.mass()/GeV, weight);

        dphitt_JJ->fill(abs(deltaPhi(t1, t2)), weight);

        pttt_JJ->fill(tt.pt()/TeV, weight);
        mtt_JJ->fill(tt.mass()/TeV, weight);


        if (addjets.size() >= 6)
          chi2_JJ->fill(chi2_hadhad(addjets), weight);

        if (addjets.size() >= 4)
          logttprob_JJ->fill(log(ttProb(topPDF0b, topPDF1b, addjets)), weight);

      } if (leps.size() == 1 && topjets.size() == 1) {
        njets_lJ->fill(jets.size(), weight);
        ncentjets_lJ->fill(centjets.size(), weight);
        nfwdjets_lJ->fill(fwdjets.size(), weight);
        ntopjets_lJ->fill(topjets.size(), weight);
        ptl1_lJ->fill(leps[0].pt()/GeV, weight);

        Jets goodtopjets = topjets;
        ntopbjets_lJ->fill(btaggedJets(goodtopjets).size(), weight);

        Jets addjetstmp = additionalJets(jets, goodtopjets);

        // look for the closest b-tagged jet to the lepton.
        // assume this is coming from the leptonically decaying top
        // quark from the resonance.
        double drmin = -1;
        int drmin_idx = -1;
        for (size_t i = 0; i < addjetstmp.size(); i++) {
          const Jet& j = addjetstmp[i];

          if (!j.bTagged(Cuts::pT > 5*GeV && Cuts::abseta < 2.5))
            continue;

          double dr = deltaR(leps[0].mom(), j.mom());
          if (drmin_idx < 0 || dr < drmin) {
            drmin = dr;
            drmin_idx = i;
          }
        }

        if (drmin_idx >= 0) {
          const Jet& bestjet = addjetstmp[drmin_idx].mom();

          Jets addjets;
          for (size_t i = 0; i < addjetstmp.size(); i++)
            if (i != drmin_idx)
              addjets.push_back(addjetstmp[i]);


          naddjets_lJ->fill(addjets.size(), weight);

          Jets addbjets = btaggedJets(addjets);
          naddbjets_lJ->fill(addbjets.size(), weight);
          naddljets_lJ->fill(addjets.size()-addbjets.size(), weight);


          // colinear approximation
          FourMomentum tl = leps[0].mom() + leps[0].mom() + bestjet;
          FourMomentum th = goodtopjets[0].mom();
          FourMomentum tt = tl + th;

          ptth_lJ->fill(th.pt()/TeV, weight);
          pttl_lJ->fill(tl.pt()/TeV, weight);

          mtl_lJ->fill(tl.mass()/GeV, weight);
          mth_lJ->fill(th.mass()/GeV, weight);

          dphitt_lJ->fill(abs(deltaPhi(tl, th)), weight);

          pttt_lJ->fill(tt.pt()/TeV, weight);
          mtt_lJ->fill(tt.mass()/TeV, weight);

          if (addjets.size() >= 6)
            chi2_lJ->fill(chi2_hadhad(addjets), weight);

          if (addjets.size() >= 4)
            logttprob_lJ->fill(log(ttProb(topPDF0b, topPDF1b, addjets)), weight);
        }

      } else if (leps.size() == 1 && topjets.size() >= 2) {
        njets_lJJ->fill(jets.size(), weight);
        ncentjets_lJJ->fill(centjets.size(), weight);
        nfwdjets_lJJ->fill(fwdjets.size(), weight);
        ntopjets_lJJ->fill(topjets.size(), weight);

        ptth1_lJJ->fill(topjets[0].pt()/TeV, weight);
        ptth2_lJJ->fill(topjets[1].pt()/TeV, weight);
        mth1_lJJ->fill(topjets[0].mass()/GeV, weight);
        mth2_lJJ->fill(topjets[1].mass()/GeV, weight);
        ptl1_lJJ->fill(leps[0].pt()/GeV, weight);

        Jets goodtopjets;
        goodtopjets.push_back(topjets[0]);
        goodtopjets.push_back(topjets[1]);
        ntopbjets_lJJ->fill(btaggedJets(goodtopjets).size(), weight);

        Jets addjets = additionalJets(jets, goodtopjets);
        naddjets_lJJ->fill(addjets.size(), weight);

        Jets addbjets = btaggedJets(addjets);
        naddbjets_lJJ->fill(addbjets.size(), weight);
        naddljets_lJJ->fill(addjets.size()-addbjets.size(), weight);

        FourMomentum t1 = goodtopjets[0].mom();
        FourMomentum t2 = goodtopjets[1].mom();
        FourMomentum tt = t1 + t2;

        mth1_lJJ->fill(t1.mass()/GeV, weight);
        mth2_lJJ->fill(t2.mass()/GeV, weight);

        dphitt_lJJ->fill(abs(deltaPhi(t1, t2)), weight);

        pttt_lJJ->fill(tt.pt()/TeV, weight);
        mtt_lJJ->fill(tt.mass()/TeV, weight);

      } else if (leps.size() == 2 && topjets.size() >= 1 && leps[0].threeCharge()*leps[1].threeCharge() > 0 ) {
        njets_ssJ->fill(jets.size(), weight);
        ncentjets_ssJ->fill(centjets.size(), weight);
        nfwdjets_ssJ->fill(fwdjets.size(), weight);
        ntopjets_ssJ->fill(topjets.size(), weight);
        ptl1_ssJ->fill(leps[0].pt()/GeV, weight);
        ptl2_ssJ->fill(leps[1].pt()/GeV, weight);

        Jets goodtopjets;
        goodtopjets.push_back(topjets[0]);
        ntopbjets_ssJ->fill(btaggedJets(goodtopjets).size(), weight);

        Jets addjets = additionalJets(jets, goodtopjets);
        naddjets_ssJ->fill(addjets.size(), weight);

        Jets addbjets = btaggedJets(addjets);
        naddbjets_ssJ->fill(addbjets.size(), weight);
        naddljets_ssJ->fill(addjets.size()-addbjets.size(), weight);

        const Jet* bestjet = NULL;
        double drmin = -1;
        for (const Jet& j : addbjets) {
          double dr = deltaR(leps[0].mom(), j.mom());
          if (drmin < 0 || dr < drmin) {
            bestjet = &j;
            drmin = dr;
          }
        }

        if (drmin >= 0) {
          // colinear approximation
          FourMomentum tl = leps[0].mom() + leps[0].mom() + bestjet->mom();
          FourMomentum th = goodtopjets[0].mom();
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
    Histo1DPtr nbjets;
    Histo1DPtr ncentjets;
    Histo1DPtr nfwdjets;
    Histo1DPtr ntopjets;
    Histo1DPtr ntopbjets;
    Histo1DPtr nleps;

    Histo1DPtr njets_JJ;
    Histo1DPtr ncentjets_JJ;
    Histo1DPtr nfwdjets_JJ;
    Histo1DPtr naddjets_JJ;
    Histo1DPtr naddbjets_JJ;
    Histo1DPtr naddljets_JJ;
    Histo1DPtr ntopjets_JJ;
    Histo1DPtr ntopbjets_JJ;

    Histo1DPtr ptth1_JJ;
    Histo1DPtr ptth2_JJ;
    Histo1DPtr mth1_JJ;
    Histo1DPtr mth2_JJ;
    Histo1DPtr dphitt_JJ;
    Histo1DPtr pttt_JJ;
    Histo1DPtr mtt_JJ;
    Histo1DPtr chi2_JJ;
    Histo1DPtr logttprob_JJ;

    Histo1DPtr njets_lJ;
    Histo1DPtr ncentjets_lJ;
    Histo1DPtr nfwdjets_lJ;
    Histo1DPtr naddjets_lJ;
    Histo1DPtr naddbjets_lJ;
    Histo1DPtr naddljets_lJ;
    Histo1DPtr ntopjets_lJ;
    Histo1DPtr ntopbjets_lJ;

    Histo1DPtr ptth_lJ;
    Histo1DPtr pttl_lJ;
    Histo1DPtr ptl1_lJ;
    Histo1DPtr mth_lJ;
    Histo1DPtr mtl_lJ;
    Histo1DPtr dphitt_lJ;
    Histo1DPtr pttt_lJ;
    Histo1DPtr mtt_lJ;
    Histo1DPtr chi2_lJ;
    Histo1DPtr logttprob_lJ;

    Histo1DPtr njets_lJJ;
    Histo1DPtr ncentjets_lJJ;
    Histo1DPtr nfwdjets_lJJ;
    Histo1DPtr naddjets_lJJ;
    Histo1DPtr naddbjets_lJJ;
    Histo1DPtr naddljets_lJJ;
    Histo1DPtr ntopjets_lJJ;
    Histo1DPtr ntopbjets_lJJ;

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
    Histo1DPtr naddbjets_ssJ;
    Histo1DPtr naddljets_ssJ;
    Histo1DPtr ntopjets_ssJ;
    Histo1DPtr ntopbjets_ssJ;

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

    Histo2D topPDF0b;
    Histo2D topPDF1b;
    //@}


  };


  // The hook for the plugin system
  DECLARE_RIVET_PLUGIN(TTTT);


}
