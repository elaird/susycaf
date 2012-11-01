import ROOT as r
import os, subprocess
import glob



canvas = r.TCanvas()
canvas.SetRightMargin(0.2)
canvas.SetTickx()
canvas.SetTicky()

models = ["T2bb","T1bbbb"]
HTbins = ["375","875"]
pdfSets = ["genMRST2006nnlo", "gencteq66"]

for model in models :

    for ht in HTbins :

            epsFileName = "%s_%s_%s_acc_ratio.eps"%(model,ht,"cteq66OverMrst")
            numFile = r.TFile("%s_%s_calo_ge2_%s.root"%(model,ht,pdfSets[0]),"READ")
            numHist = numFile.Get("nEvents")
    
            denFile = r.TFile("%s_%s_calo_ge2_%s.root"%(model,ht,pdfSets[1]),"READ")
            denHist = denFile.Get("nEvents")


            result = numHist.Clone()
            result.Divide(denHist)
    
            result.SetTitle(";m_{parent} (GeV);m_{LSP} (GeV);ratio") 
            result.SetMarkerStyle(20)
            result.SetStats(False)

            if "T1bbbb" in model :
                result.SetMaximum(1.15)
                result.SetMinimum(0.65)
                line = r.TLine(300,125,2025,1850)
                #line = r.TLine(300,125,1200,1025)
                line2 = r.TLine(300,50,300,125)
                lineDiag = r.TLine(100,100,2025,2025)
                
            if "T2bb" in model :
                result.GetXaxis().SetRangeUser(0,1200)
                result.GetYaxis().SetRangeUser(0,1200)
                result.SetMaximum(1.06)
                result.SetMinimum(0.85)
                line = r.TLine(300,125,1225,1050)
                #line = r.TLine(300,125,1200,1025)
                line2 = r.TLine(300,50,300,125)
                lineDiag = r.TLine(100,100,1225,1225)

            result.Draw("colz")

            lineDiag.SetLineStyle(2)
           
            line.SetLineWidth(2)
            line2.SetLineWidth(2)
            lineDiag.SetLineWidth(2)
           
            line.Draw("lsame")
            line2.Draw("lsame")
            lineDiag.Draw("lsame")
           
            canvas.Print(epsFileName)
            result.Write()


           
            os.system("epstopdf "+ epsFileName)
            os.remove(epsFileName)



