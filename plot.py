GeV = 1
TeV = GeV*1e3
MeV = GeV*1e-3


def main(infname, outfname):
    import ROOT
    import pylhe

    outf = ROOT.TFile(outfname, "CREATE")

    hmtt = ROOT.TH1F("hmtt", "hmtt", 100, 0, 2*TeV)
    hptt1 = ROOT.TH1F("hptt1", "hptt1", 100, 0, 1*TeV)
    hptt2 = ROOT.TH1F("hptt2", "hptt2", 100, 0, 1*TeV)
    hptt3 = ROOT.TH1F("hptt3", "hptt3", 100, 0, 1*TeV)
    hptt4 = ROOT.TH1F("hptt4", "hptt4", 100, 0, 1*TeV)
    habsetat1 = ROOT.TH1F("habsetat1", "habsetat1", 100, 0, 5)
    habsetat2 = ROOT.TH1F("habsetat2", "habsetat2", 100, 0, 5)
    habsetat3 = ROOT.TH1F("habsetat3", "habsetat3", 100, 0, 5)
    habsetat4 = ROOT.TH1F("habsetat4", "habsetat4", 100, 0, 5)

    for evt in pylhe.readLHE(infname):
        tops = \
            [ ROOT.TLorentzVector(t.px, t.py, t.pz, t.e)
              for t in evt.particles if abs(t.id) == 6
            ]

        tops.sort(key = lambda t: t.Pt(), reverse=True)

        if len(tops) != 4:
            print "woops!"
            continue
        else:
            hmtt.Fill((tops[0] + tops[1]).M(), evt.eventinfo.weight)
            hptt1.Fill(tops[0].Pt(), evt.eventinfo.weight)
            hptt2.Fill(tops[1].Pt(), evt.eventinfo.weight)
            hptt3.Fill(tops[2].Pt(), evt.eventinfo.weight)
            hptt4.Fill(tops[3].Pt(), evt.eventinfo.weight)

            habsetat1.Fill(abs(tops[0].Eta()), evt.eventinfo.weight)
            habsetat2.Fill(abs(tops[1].Eta()), evt.eventinfo.weight)
            habsetat3.Fill(abs(tops[2].Eta()), evt.eventinfo.weight)
            habsetat4.Fill(abs(tops[3].Eta()), evt.eventinfo.weight)

    # scale by 100/fb
    hmtt.Scale(100e3)
    hptt1.Scale(100e3)
    hptt2.Scale(100e3)
    hptt3.Scale(100e3)
    hptt4.Scale(100e3)
    habsetat1.Scale(100e3)
    habsetat2.Scale(100e3)
    habsetat3.Scale(100e3)
    habsetat4.Scale(100e3)

    outf.Write()
    outf.Close()

    return



def fromV(p):
    for m in p.mothers():
        if m.id == 6000055:
            return True

    return False


if __name__ == "__main__":
    from sys import argv
    main(argv[1], argv[2])
