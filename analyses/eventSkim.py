import supy
import steps,samples

class eventSkim(supy.analysis) :
    def listOfSteps(self,_) :
        return [ steps.printer.progressPrinter(2,300),
                 steps.other.runLsEventFilter("/home/hep/elaird1/84_darrens_event/list.txt"),
                 #steps.other.runLsEventFilter("/home/hep/elaird1/75_rob_sync/v1/robs_events/ht300.txt"),
                 #steps.other.runLsEventFilter("/home/hep/elaird1/58_wpol_events/v2/38_misRun_event_ls.txt"),
                 #steps.other.runLsEventFilter("/home/hep/elaird1/58_wpol_events/v2/39_Run_event_ls.txt"),
                 #steps.other.runLsEventFilter("/home/hep/elaird1/58_wpol_events/v2/39_misRun_event_ls.txt"),
                 steps.other.skimmer(),
                 ]

    def listOfCalculables(self,_) :
        return supy.calculables.zeroArgs()

    def listOfSamples(self,_) :
        return supy.samples.specify(names = ["qcd_mg_ht_100_250", "qcd_mg_ht_250_500", "qcd_mg_ht_500_1000", "qcd_mg_ht_1000_inf"] )

    def listOfSampleDictionaries(self) :
        return [samples.jetmet,samples.mc]
