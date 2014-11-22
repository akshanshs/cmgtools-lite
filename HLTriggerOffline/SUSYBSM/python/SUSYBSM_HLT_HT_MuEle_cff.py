import FWCore.ParameterSet.Config as cms

SUSY_HLT_HT_MuEle = cms.EDAnalyzer("SUSY_HLT_MuEle_Hadronic",
  trigSummary = cms.InputTag("hltTriggerSummaryAOD",'', 'HLT'), #to use with test sample
  #trigSummary = cms.InputTag("hltTriggerSummaryAOD"),
  MuonCollection = cms.InputTag("muons"),
  ElectronCollection = cms.InputTag("gedGsfElectrons"),
  pfJetCollection = cms.InputTag("ak4PFJetsCHS"),
  caloJetCollection = cms.InputTag("ak4CaloJets"),
  TriggerResults = cms.InputTag('TriggerResults','','HLT'), #to use with test sample
  #TriggerResults = cms.InputTag('TriggerResults','','HLT'),
  HLTProcess = cms.string('HLT'),
  TriggerPath = cms.string('HLT_Mu8_Ele8_CaloIdL_TrkIdVL_Mass8_PFNoPUHT300_v'),
  TriggerPathAuxiliaryForMuEle = cms.string('HLT_PFHT900_v'),
  TriggerPathAuxiliaryForHadronic = cms.string('HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v'),
  TriggerFilter = cms.InputTag('HLTElectronMuonInvMassFilter', '', 'HLT'), #the last filter in the path 
  PtThrJet = cms.untracked.double(40.0),
  EtaThrJet = cms.untracked.double(3.0)
)


SUSY_HLT_HT_MuEle_FASTSIM = cms.EDAnalyzer("SUSY_HLT_MuEle_Hadronic",
  trigSummary = cms.InputTag("hltTriggerSummaryAOD",'', 'HLT'), #to use with test sample
  #trigSummary = cms.InputTag("hltTriggerSummaryAOD"),
  MuonCollection = cms.InputTag("muons"),
  ElectronCollection = cms.InputTag("gedGsfElectrons"),
  pfJetCollection = cms.InputTag("ak4PFJetsCHS"),
  caloJetCollection = cms.InputTag("ak4CaloJets"),
  TriggerResults = cms.InputTag('TriggerResults','','HLT'), #to use with test sample
  #TriggerResults = cms.InputTag('TriggerResults','','HLT'),
  HLTProcess = cms.string('HLT'),
  TriggerPath = cms.string('HLT_Mu8_Ele8_CaloIdL_TrkIdVL_Mass8_PFNoPUHT300_v'),
  TriggerPathAuxiliaryForMuEle = cms.string('HLT_PFHT900_v'),
  TriggerPathAuxiliaryForHadronic = cms.string('HLT_Mu23_TrkIsoVVL_Ele12_Gsf_CaloId_TrackId_Iso_MediumWP_v'),
  TriggerFilter = cms.InputTag('HLTElectronMuonInvMassFilter', '', 'HLT'), #the last filter in the path
  PtThrJet = cms.untracked.double(40.0),
  EtaThrJet = cms.untracked.double(3.0)
)


SUSY_HLT_HT_MuEle_POSTPROCESSING = cms.EDAnalyzer("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/SUSYBSM/HLT_Mu8_Ele8_CaloId_TrkIdVL_Mass8_PFHTT300"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(""),
    efficiency     = cms.vstring(
       "pfHTTurnOn_eff 'Turn-on vs HT; PFHT (GeV); #epsilon' pfHTTurnOn_num pfHTTurnOn_den",
       "MuTurnOn_eff 'Turn-on vs Mu pT; pT (GeV); #epsilon' MuTurnOn_num MuTurnOn_den",
       "EleTurnOn_eff 'Turn-on vs Ele pT; pT (GeV); #epsilon' EleTurnOn_num EleTurnOn_den",
    )
)


SUSY_HLT_HT_MuEle_FASTSIM_POSTPROCESSING = cms.EDAnalyzer("DQMGenericClient",
    subDirs        = cms.untracked.vstring("HLT/SUSYBSM/HLT_Mu8_Ele8_CaloId_TrkIdVL_Mass8_PFHTT300"),
    verbose        = cms.untracked.uint32(2), # Set to 2 for all messages
    resolution     = cms.vstring(""),
    efficiency     = cms.vstring(
       "pfHTTurnOn_eff 'Turn-on vs HT; PFHT (GeV); #epsilon' pfHTTurnOn_num pfHTTurnOn_den",
       "MuTurnOn_eff 'Turn-on vs Mu pT; pT (GeV); #epsilon' MuTurnOn_num MuTurnOn_den",
       "EleTurnOn_eff 'Turn-on vs Ele pT; pT (GeV); #epsilon' EleTurnOn_num EleTurnOn_den",
    )
)



