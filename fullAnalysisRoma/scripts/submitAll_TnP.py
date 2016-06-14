#
# usage: %prog [opts] --cfg cmssw.py dataset doPUreweighting(0/1) sampleIndex PUweightsFile x-section kFactor
#
# Backgrounds: sampleID>0 && sampleID<100
# Signal (DY): sampleID>100 && sampleID<10000
# Data:        sampleID>10000

# MC DY signal
#./submitBatchTnP.py --cfg tnpAnaBATCH.py DYLL 1 101 pileupWeights___processedAndGolden_2016B_june10__69mb.root 2008  1 

# Data
##./submitBatchTnP.py --cfg tnpAnaBATCH.py singleEle2016B_p3v1  0 10001 pippo 1  1 
./submitBatchTnP.py --cfg tnpAnaBATCH.py singleEle2016B_p3v2  0 10001 pippo 1  1 
./submitBatchTnP.py --cfg tnpAnaBATCH.py singleEle2016B_p5v2  0 10001 pippo 1  1 
./submitBatchTnP.py --cfg tnpAnaBATCH.py singleEle2016B_p6v2  0 10001 pippo 1  1 

