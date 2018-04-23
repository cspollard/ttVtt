GeV = 1
TeV = GeV*1e3
MeV = GeV*1e-3

# from
# https://stackoverflow.com/questions/8793772/how-to-split-a-sequence-according-to-a-predicate
def partition(pred, iterable):
  from itertools import tee

  def notpred(x):
    return not pred(x)

  'Use a predicate to partition entries into false entries and true entries'
  # partition(is_odd, range(10)) --> 0 2 4 6 8   and  1 3 5 7 9
  t1, t2 = tee(iterable)
  return filter(pred, t1,), filter(notpred, t2)


def hasAncestor(p, pid):
  if p.id == pid:
    return True
  else:
    for q in p.mothers():
      if hasAncestor(q, pid):
        return True
      else:
        continue

  return False


def main(infname, outfname):
  import ROOT
  import pylhe
  import itertools as it

  outf = ROOT.TFile(outfname, "CREATE")

  hpz_g1 = ROOT.TH1F("hpz_g1", "hpz_g1", 100, 0, 4*TeV)
  hpz_g2 = ROOT.TH1F("hpz_g2", "hpz_g2", 100, 0, 4*TeV)
  hmgg = ROOT.TH1F("hmgg", "hmgg", 100, 0, 4*TeV)
  hpzgg = ROOT.TH1F("hpzgg", "hpzgg", 100, 0, 4*TeV)

  hptv1 = ROOT.TH1F("hptv1", "hptv1", 100, 0, 2*TeV)
  hpzv1 = ROOT.TH1F("hpzv1", "hpzv1", 100, 0, 2*TeV)

  hptt1 = ROOT.TH1F("hptt1", "hptt1", 100, 0, 1*TeV)
  hptt2 = ROOT.TH1F("hptt2", "hptt2", 100, 0, 1*TeV)
  hptt3 = ROOT.TH1F("hptt3", "hptt3", 100, 0, 1*TeV)
  hptt4 = ROOT.TH1F("hptt4", "hptt4", 100, 0, 1*TeV)
  habsetat1 = ROOT.TH1F("habsetat1", "habsetat1", 100, 0, 5)
  habsetat2 = ROOT.TH1F("habsetat2", "habsetat2", 100, 0, 5)
  habsetat3 = ROOT.TH1F("habsetat3", "habsetat3", 100, 0, 5)
  habsetat4 = ROOT.TH1F("habsetat4", "habsetat4", 100, 0, 5)

  hmtt_centraltops = ROOT.TH1F("hmtt_centraltops", "hmtt_centraltops", 100, 0, 2*TeV)
  hmtt_leadingtops = ROOT.TH1F("hmtt_leadingtops", "hmtt_leadingtops", 100, 0, 2*TeV)

  hmtt_fromV = ROOT.TH1F("hmtt_fromV", "hmtt_fromV", 100, 0, 2*TeV)

  hpt_t1_fromV = ROOT.TH1F("hpt_t1_fromV", "hpt_t1_fromV", 100, 0, 1*TeV)
  hpt_t2_fromV = ROOT.TH1F("hpt_t2_fromV", "hpt_t2_fromV", 100, 0, 1*TeV)
  hpt_t1_notFromV = ROOT.TH1F("hpt_t1_notFromV", "hpt_t1_notFromV", 100, 0, 1*TeV)
  hpt_t2_notFromV = ROOT.TH1F("hpt_t2_notFromV", "hpt_t2_notFromV", 100, 0, 1*TeV)

  hpz_t1_fromV = ROOT.TH1F("hpz_t1_fromV", "hpz_t1_fromV", 100, 0, 1*TeV)
  hpz_t2_fromV = ROOT.TH1F("hpz_t2_fromV", "hpz_t2_fromV", 100, 0, 1*TeV)
  hpz_t1_notFromV = ROOT.TH1F("hpz_t1_notFromV", "hpz_t1_notFromV", 100, 0, 1*TeV)
  hpz_t2_notFromV = ROOT.TH1F("hpz_t2_notFromV", "hpz_t2_notFromV", 100, 0, 1*TeV)

  habseta_t1_fromV = ROOT.TH1F("habseta_t1_fromV", "habseta_t1_fromV", 100, 0, 5)
  habseta_t2_fromV = ROOT.TH1F("habseta_t2_fromV", "habseta_t2_fromV", 100, 0, 5)
  habseta_t1_notFromV = ROOT.TH1F("habseta_t1_notFromV", "habseta_t1_notFromV", 100, 0, 5)
  habseta_t2_notFromV = ROOT.TH1F("habseta_t2_notFromV", "habseta_t2_notFromV", 100, 0, 5)

  hpt_t3_notFromV = ROOT.TH1F("hpt_t3_notFromV", "hpt_t3_notFromV", 100, 0, 1*TeV)
  hpt_t4_notFromV = ROOT.TH1F("hpt_t4_notFromV", "hpt_t4_notFromV", 100, 0, 1*TeV)
  habseta_t3_notFromV = ROOT.TH1F("habseta_t3_notFromV", "habseta_t3_notFromV", 100, 0, 5)
  habseta_t4_notFromV = ROOT.TH1F("habseta_t4_notFromV", "habseta_t4_notFromV", 100, 0, 5)

  allHists = \
    [ hpz_g1, hpz_g2, hmgg
    , hptv1, hpzv1
    , hptt1, hptt2, hptt3, hptt4
    , habsetat1, habsetat2 , habsetat3, habsetat4
    , hmtt_centraltops, hmtt_leadingtops, hmtt_fromV 
    , hpt_t1_fromV, hpt_t2_fromV
    , hpt_t1_notFromV, hpt_t2_notFromV
    , hpz_t1_fromV, hpz_t2_fromV
    , hpz_t1_notFromV, hpz_t2_notFromV
    , habseta_t1_fromV, habseta_t2_fromV
    , habseta_t1_notFromV, habseta_t2_notFromV
    , hpt_t3_notFromV, hpt_t4_notFromV
    , habseta_t3_notFromV, habseta_t4_notFromV
    ]

  for evt in pylhe.readLHE(infname):
    glus = [ p for p in evt.particles if abs(p.id) == 21 ]
    v1s = [ p for p in evt.particles if abs(p.id) == 6000055 ]
    tops = [ p for p in evt.particles if abs(p.id) == 6 ]

    # attach the TLV to the particles
    for p in glus:
      p.tlv = ROOT.TLorentzVector(p.px, p.py, p.pz, p.e)
    for p in v1s:
      p.tlv = ROOT.TLorentzVector(p.px, p.py, p.pz, p.e)
    for p in tops:
      p.tlv = ROOT.TLorentzVector(p.px, p.py, p.pz, p.e)

    tops.sort(key = lambda t: t.tlv.Pt(), reverse=True)
    wgt = evt.eventinfo.weight

    if len(v1s) == 1:
      hptv1.Fill(v1s[0].tlv.Pt(), wgt)
      hpzv1.Fill(abs(v1s[0].tlv.Pz()), wgt)

    if len(glus) == 2:
      hpz_g1.Fill(abs(glus[0].tlv.Pz()), wgt)
      hpz_g2.Fill(abs(glus[1].tlv.Pz()), wgt)
      gg = glus[0].tlv + glus[1].tlv
      hmgg.Fill(gg.M(), wgt)
      hpzgg.Fill(gg.Pz(), wgt)

    if len(tops) != 4:
      print "woops!"
      continue

    hptt1.Fill(tops[0].tlv.Pt(), wgt)
    hptt2.Fill(tops[1].tlv.Pt(), wgt)
    hptt3.Fill(tops[2].tlv.Pt(), wgt)
    hptt4.Fill(tops[3].tlv.Pt(), wgt)

    habsetat1.Fill(abs(tops[0].tlv.Eta()), wgt)
    habsetat2.Fill(abs(tops[1].tlv.Eta()), wgt)
    habsetat3.Fill(abs(tops[2].tlv.Eta()), wgt)
    habsetat4.Fill(abs(tops[3].tlv.Eta()), wgt)

    hmtt_leadingtops.Fill((tops[0].tlv + tops[1].tlv).M(), wgt)

    eta_ordered_tops = sorted(map(lambda t: (abs(t.tlv.Eta()), t), tops))
    hmtt_centraltops.Fill(
        (eta_ordered_tops[0][1].tlv + eta_ordered_tops[1][1].tlv).M(), wgt)

    (topsFromV, topsNotFromV) = list(partition(lambda p: hasAncestor(p, 6000055), tops))

    if len(topsFromV) == 2:
      hpt_t1_fromV.Fill(topsFromV[0].tlv.Pt(), wgt)
      hpt_t2_fromV.Fill(topsFromV[1].tlv.Pt(), wgt)
      hpz_t1_fromV.Fill(topsFromV[0].tlv.Pz(), wgt)
      hpz_t2_fromV.Fill(topsFromV[1].tlv.Pz(), wgt)
      habseta_t1_fromV.Fill(abs(topsFromV[0].tlv.Eta()), wgt)
      habseta_t2_fromV.Fill(abs(topsFromV[1].tlv.Eta()), wgt)
      hmtt_fromV.Fill(
          (topsFromV[0].tlv + topsFromV[1].tlv).M(), wgt)

    if len(topsNotFromV) >= 2:
      hpt_t1_notFromV.Fill(topsNotFromV[0].tlv.Pt(), wgt)
      hpt_t2_notFromV.Fill(topsNotFromV[1].tlv.Pt(), wgt)
      hpz_t1_notFromV.Fill(topsNotFromV[0].tlv.Pz(), wgt)
      hpz_t2_notFromV.Fill(topsNotFromV[1].tlv.Pz(), wgt)
      habseta_t1_notFromV.Fill(abs(topsNotFromV[0].tlv.Eta()), wgt)
      habseta_t2_notFromV.Fill(abs(topsNotFromV[1].tlv.Eta()), wgt)

    if len(topsNotFromV) >= 4:
      hpt_t3_notFromV.Fill(topsNotFromV[2].tlv.Pt(), wgt)
      hpt_t4_notFromV.Fill(topsNotFromV[3].tlv.Pt(), wgt)
      habseta_t3_notFromV.Fill(abs(topsNotFromV[2].tlv.Eta()), wgt)
      habseta_t4_notFromV.Fill(abs(topsNotFromV[3].tlv.Eta()), wgt)

  # scale by 100/fb
  for h in allHists:
    h.Scale(100e3)

  outf.Write()
  outf.Close()

  return


if __name__ == "__main__":
  from sys import argv
  main(argv[1], argv[2])
