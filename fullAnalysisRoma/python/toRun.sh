# in locale
# fggRunJobs.py --load jobs_singleelectron_80.json -x cmsRun tnpAnaFlashgg.py maxEvents=10 -d pippo

# in batch
fggRunJobs.py --load jobs_singleelectron_80.json -x cmsRun tnpAnaFlashgg.py -d pippo -n 50 -q 8nh  &
