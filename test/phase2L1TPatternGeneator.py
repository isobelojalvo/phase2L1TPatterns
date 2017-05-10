import FWCore.ParameterSet.Config as cms
import os
process = cms.Process("patternGeneration")



# import standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4Reco_cff')
process.load('Configuration.Geometry.GeometryExtended2023D4_cff')

process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
# ---- Global Tag :

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:upgradePLS3', '')

process.load('Configuration.Geometry.GeometryExtended2023D4_cff')

# input
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(1000) )
process.source = cms.Source("PoolSource",
   fileNames = cms.untracked.vstring(
        'file:singleE-gen-sim-raw.root'
        #'file:event8021-pi0.root'
        )
)


############################################################
# L1 tracking
############################################################

process.load("L1Trigger.TrackFindingTracklet.L1TrackletTracks_cff")

# run only the tracking (no MC truth associators)
process.TTTracks = cms.Path(process.L1TrackletTracks)

# run the tracking AND MC truth associators)
process.TTTracksWithTruth = cms.Path(process.L1TrackletTracksWithAssociators)


# L1 Cluster Producer
process.load("L1Trigger.phase2Demonstrator.L1CaloClusterProducer_cff")
process.L1Clusters = cms.Path(process.L1CaloClusterProducer)


# ----------------------------------------------------------------------------------------------
# 
# Analyzer starts here

process.patterns = cms.EDAnalyzer('Phase2L1TPatternGenerator',
                                  L1TrackInputTag = cms.InputTag("TTTracksFromTracklet", "Level1TTTracks"),               ## TTTrack input
                                  L1TrackPrimaryVertexTag = cms.InputTag("L1TkPrimaryVertex"),
                                  summaryCardInputFileName  = cms.untracked.string("inputPatterns.txt"),
                                  summaryCardOutputFileName = cms.untracked.string("outputPatterns.txt"),
                                  L1Clusters = cms.InputTag("L1CaloClusterProducer","L1Phase2CaloClusters"),
                                  ecalTPGsBarrel = cms.InputTag("simEcalEBTriggerPrimitiveDigis","","HLT"),
                                  hcalTPGs       = cms.InputTag("simHcalTriggerPrimitiveDigis","","HLT")
)


# output module
process.out = cms.OutputModule( "PoolOutputModule",
                                fileName = cms.untracked.string("CaloClusters.root"),
                                fastCloning = cms.untracked.bool( False ),
                                outputCommands = cms.untracked.vstring('keep *',
                                                                       #'keep *_*_L1Phase2CaloClusters_*', 
                                                                       )
)
process.FEVToutput_step = cms.EndPath(process.out)

process.p = cms.Path(process.patterns)

process.schedule = cms.Schedule(process.TTTracksWithTruth,process.L1Clusters,process.p)
