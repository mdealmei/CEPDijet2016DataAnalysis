from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.General.requestName = 'test_RunB'
config.General.workArea = 'crab_projects'
config.General.transferOutputs = True
config.General.transferLogs = True

config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'miniAOD-prod_PAT_filtered.py'

config.Data.inputDataset = '/JetHT/Run2016B-23Sep2016-v3/AOD'
config.Data.inputDBS = 'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 200000#20
#config.Data.totalUnits = 1 ## how many files to analyze, just for test
config.Data.lumiMask = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt'
config.Data.outLFNDirBase = '/store/user/mdealmei'
config.Data.publication = False
config.Data.outputDatasetTag = 'Run2016B-miniAOD-custom'

config.Site.storageSite = 'T2_CH_CERN'
