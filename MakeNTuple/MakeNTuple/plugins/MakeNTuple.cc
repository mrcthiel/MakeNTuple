// -*- C++ -*-
//
// Package:    MakeNTuple/MakeNTuple
// Class:      MakeNTuple
// 
/**\class MakeNTuple MakeNTuple.cc MakeNTuple/MakeNTuple/plugins/MakeNTuple.cc

Description: [one line class summary]

Implementation:
[Notes on implementation]
*/
//
// Original Author:  Mauricio Thiel
//         Created:  Tue, 09 Apr 2019 16:31:41 GMT
//
//


// system include files
#include <memory>
#include <iostream> 
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "DataFormats/VertexReco/interface/Vertex.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "TRandom.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "TLorentzVector.h" 
#include <vector>
#include <TMath.h>

//HLT
#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"


//PAT
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/JetReco/interface/GenJet.h"
#include "DataFormats/JetReco/interface/GenJetCollection.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/MET.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

//ROOT
#include "TTree.h" 
#include <TGraph.h>

// P P S
/*
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastRecHit.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastRecHitContainer.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastTrack.h"
#include "FastSimDataFormats/CTPPSFastSim/interface/CTPPSFastTrackContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
*/
#include "DataFormats/CTPPSDetId/interface/CTPPSDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSLocalTrackLite.h"
#include "DataFormats/ProtonReco/interface/ForwardProton.h"

#include "CondFormats/RunInfo/interface/LHCInfo.h"
#include "CondFormats/DataRecord/interface/LHCInfoRcd.h"

using namespace std; 
using namespace reco;
using namespace edm;

//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.

class MakeNTuple : public edm::one::EDAnalyzer<edm::one::SharedResources>  {
	public:
		explicit MakeNTuple(const edm::ParameterSet&);
		~MakeNTuple();

		static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


	private:
		virtual void beginJob() override;
		virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
		virtual void endJob() override;

		// ----------member data ---------------------------
		edm::Handle<pat::JetCollection> jetsAK8;
		edm::EDGetTokenT<pat::JetCollection> jetsAK8Token;
		edm::Handle<pat::JetCollection> jetsAK4;
		edm::EDGetTokenT<pat::JetCollection> jetsAK4Token;
		edm::Handle<reco::GenJetCollection> genjets;
		edm::EDGetTokenT<reco::GenJetCollection> genJetsToken;
		edm::Handle<reco::GenJetCollection> genjetsAK8;
		edm::EDGetTokenT<reco::GenJetCollection> genJetsTokenAK8;
		edm::Handle<pat::MuonCollection> muons;
		edm::EDGetTokenT<pat::MuonCollection> muonsToken;
		edm::Handle<pat::ElectronCollection> electrons;
		edm::EDGetTokenT<pat::ElectronCollection> electronsToken;
		edm::Handle<pat::METCollection> MET;
		edm::EDGetTokenT<pat::METCollection> MetToken;
		edm::Handle<std::vector<pat::PackedCandidate>> PFCand;
		edm::EDGetTokenT<std::vector<pat::PackedCandidate>> PFCandToken;
		edm::Handle<std::vector<PileupSummaryInfo> > PupInfo;
		edm::EDGetTokenT<std::vector<PileupSummaryInfo> > PileupSumInfoInputTag;
		edm::Handle<reco::VertexCollection> vertices;
		edm::EDGetTokenT<reco::VertexCollection> verticesToken;
//		edm::Handle<reco::GenParticleCollection> genP;
//		edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

		edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsSingleRPToken_;
		edm::EDGetTokenT<std::vector<reco::ForwardProton> > recoProtonsMultiRPToken_;

		edm::LumiReWeighting *LumiWeights_;

		bool MC, MC_Signal, DATA;

		edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;



		//FOR PU REWGT

		TTree* EventBranchs;
		int BX, Run, LumiSection, EventNum, xangle;
		double PUWeight;
		double nPU;

		//VERTEX INFO
		int nVtx;
		bool vtx_isValid, vtx_isFake;
		double vtx_z;



		HLTConfigProvider hltConfig_;
		HLTPrescaleProvider hltPrescaleProvider_;

		std::vector<double> *ArmF_xi_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_xi_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_t_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_t_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_thx_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_thx_gen        = new std::vector<double> ();
		std::vector<double> *ArmF_thy_gen        = new std::vector<double> ();
		std::vector<double> *ArmB_thy_gen        = new std::vector<double> ();

		std::vector<double> *ProtCand_xi		= new std::vector<double> ();
		std::vector<double> *ProtCand_t			= new std::vector<double> ();
		std::vector<double> *ProtCand_ThX		= new std::vector<double> ();
		std::vector<double> *ProtCand_ThY		= new std::vector<double> ();
		std::vector<double> *ProtCand_rpid		= new std::vector<double> ();
		std::vector<double> *ProtCand_arm       	= new std::vector<double> ();
		std::vector<double> *ProtCand_ismultirp		= new std::vector<double> ();


		std::vector<double> *muon_px        = new std::vector<double> ();
		std::vector<double> *muon_py        = new std::vector<double> ();
		std::vector<double> *muon_pz        = new std::vector<double> ();
		std::vector<double> *muon_pt        = new std::vector<double> ();
		std::vector<double> *muon_E         = new std::vector<double> ();
		std::vector<double> *muon_vtxZ      = new std::vector<double> ();
		std::vector<double> *muon_phi       = new std::vector<double> ();
		std::vector<double> *muon_eta       = new std::vector<double> ();
		std::vector<double> *muon_PFBasedIso     = new std::vector<double> ();
		std::vector<double> *muon_TrackBasedIso     = new std::vector<double> ();
		std::vector<bool> *muon_isTightMuon     = new std::vector<bool> ();
		std::vector<bool> *muon_isMediumMuon     = new std::vector<bool> ();
		std::vector<bool> *muon_isLooseMuon     = new std::vector<bool> ();
		std::vector<bool> *muon_isHighPtMuon     = new std::vector<bool> ();


		std::vector<double> *electron_px        = new std::vector<double> ();
		std::vector<double> *electron_py        = new std::vector<double> ();
		std::vector<double> *electron_pz        = new std::vector<double> ();
		std::vector<double> *electron_pt        = new std::vector<double> ();
		std::vector<double> *electron_E         = new std::vector<double> ();
		std::vector<double> *electron_vtxZ      = new std::vector<double> ();
		std::vector<double> *electron_phi       = new std::vector<double> ();
		std::vector<double> *electron_eta       = new std::vector<double> ();
		std::vector<bool> *electron_isTightElectron     = new std::vector<bool> ();
		std::vector<bool> *electron_isMediumElectron     = new std::vector<bool> ();
		std::vector<bool> *electron_isLooseElectron     = new std::vector<bool> ();
		std::vector<bool> *electron_isVetoElectron     = new std::vector<bool> ();

		std::vector<double> *jetAK4_px        = new std::vector<double> ();
		std::vector<double> *jetAK4_py        = new std::vector<double> ();
		std::vector<double> *jetAK4_pz        = new std::vector<double> ();
		std::vector<double> *jetAK4_pt        = new std::vector<double> ();
		std::vector<double> *jetAK4_E         = new std::vector<double> ();
		std::vector<double> *jetAK4_phi       = new std::vector<double> ();
		std::vector<double> *jetAK4_eta       = new std::vector<double> ();
		std::vector<double> *jetAK4_btag       = new std::vector<double> ();
		std::vector<bool> *jetAK4_isLoose       = new std::vector<bool> ();
		std::vector<bool> *jetAK4_isTight       = new std::vector<bool> ();

		std::vector<double> *genJets_px        = new std::vector<double> ();
		std::vector<double> *genJets_py        = new std::vector<double> ();
		std::vector<double> *genJets_pz        = new std::vector<double> ();
		std::vector<double> *genJets_pt        = new std::vector<double> ();
		std::vector<double> *genJets_E         = new std::vector<double> ();
		std::vector<double> *genJets_phi       = new std::vector<double> ();
		std::vector<double> *genJets_eta       = new std::vector<double> ();


		std::vector<double> *jetAK8_px        = new std::vector<double> ();
		std::vector<double> *jetAK8_py        = new std::vector<double> ();
		std::vector<double> *jetAK8_pz        = new std::vector<double> ();
		std::vector<double> *jetAK8_pt        = new std::vector<double> ();
		std::vector<double> *jetAK8_E         = new std::vector<double> ();
		std::vector<double> *jetAK8_phi       = new std::vector<double> ();
		std::vector<double> *jetAK8_eta       = new std::vector<double> ();
		std::vector<double> *jetAK8_btag       = new std::vector<double> ();
		std::vector<bool> *jetAK8_isLoose       = new std::vector<bool> ();
		std::vector<bool> *jetAK8_isTight       = new std::vector<bool> ();
		std::vector<double> *jetAK8_prunedMass       = new std::vector<double> ();
		std::vector<double> *jetAK8_tau21       = new std::vector<double> ();

		std::vector<double> *genJetsAK8_px        = new std::vector<double> ();
		std::vector<double> *genJetsAK8_py        = new std::vector<double> ();
		std::vector<double> *genJetsAK8_pz        = new std::vector<double> ();
		std::vector<double> *genJetsAK8_pt        = new std::vector<double> ();
		std::vector<double> *genJetsAK8_E         = new std::vector<double> ();
		std::vector<double> *genJetsAK8_phi       = new std::vector<double> ();
		std::vector<double> *genJetsAK8_eta       = new std::vector<double> ();


		double METPx, METPy, METPt, METphi;

		std::vector<double> *pfphi      = new std::vector<double> ();
		std::vector<double> *pfeta      = new std::vector<double> ();

		std::vector<int> *pffromPV      = new std::vector<int> ();
		std::vector<double> *pfdz      = new std::vector<double> ();
		std::vector<double> *pfpt      = new std::vector<double> ();

		std::vector<std::string> *HLT_name        = new std::vector<std::string> ();
		std::vector<bool> *HLT_pass        = new std::vector<bool> ();
		std::vector<int> *HLT_prescale        = new std::vector<int> ();


		std::vector<double> *muon_gen_px        = new std::vector<double> ();
		std::vector<double> *muon_gen_py       = new std::vector<double> ();
		std::vector<double> *muon_gen_pz        = new std::vector<double> ();
		std::vector<double> *muon_gen_pt        = new std::vector<double> ();
		std::vector<double> *muon_gen_E        = new std::vector<double> ();
		std::vector<double> *muon_gen_phi        = new std::vector<double> ();
		std::vector<double> *muon_gen_eta        = new std::vector<double> ();

		std::vector<double> *neut_gen_px        = new std::vector<double> ();
		std::vector<double> *neut_gen_py       = new std::vector<double> ();
		std::vector<double> *neut_gen_pz        = new std::vector<double> ();
		std::vector<double> *neut_gen_pt        = new std::vector<double> ();
		std::vector<double> *neut_gen_E        = new std::vector<double> ();
		std::vector<double> *neut_gen_phi        = new std::vector<double> ();
		std::vector<double> *neut_gen_eta        = new std::vector<double> ();

		std::vector<double> *ele_gen_px        = new std::vector<double> ();
		std::vector<double> *ele_gen_py       = new std::vector<double> ();
		std::vector<double> *ele_gen_pz        = new std::vector<double> ();
		std::vector<double> *ele_gen_pt        = new std::vector<double> ();
		std::vector<double> *ele_gen_E        = new std::vector<double> ();
		std::vector<double> *ele_gen_phi        = new std::vector<double> ();
		std::vector<double> *ele_gen_eta        = new std::vector<double> ();

		std::vector<double> *qrk_gen_px        = new std::vector<double> ();
		std::vector<double> *qrk_gen_py       = new std::vector<double> ();
		std::vector<double> *qrk_gen_pz        = new std::vector<double> ();
		std::vector<double> *qrk_gen_pt        = new std::vector<double> ();
		std::vector<double> *qrk_gen_E        = new std::vector<double> ();
		std::vector<double> *qrk_gen_phi        = new std::vector<double> ();
		std::vector<double> *qrk_gen_eta        = new std::vector<double> ();

		std::vector<double> *Wlep_gen_px        = new std::vector<double> ();
		std::vector<double> *Wlep_gen_py       = new std::vector<double> ();
		std::vector<double> *Wlep_gen_pz        = new std::vector<double> ();
		std::vector<double> *Wlep_gen_pt        = new std::vector<double> ();
		std::vector<double> *Wlep_gen_E        = new std::vector<double> ();
		std::vector<double> *Wlep_gen_phi        = new std::vector<double> ();
		std::vector<double> *Wlep_gen_eta        = new std::vector<double> ();

		std::vector<double> *Whad_gen_px        = new std::vector<double> ();
		std::vector<double> *Whad_gen_py       = new std::vector<double> ();
		std::vector<double> *Whad_gen_pz        = new std::vector<double> ();
		std::vector<double> *Whad_gen_pt        = new std::vector<double> ();
		std::vector<double> *Whad_gen_E        = new std::vector<double> ();
		std::vector<double> *Whad_gen_phi        = new std::vector<double> ();
		std::vector<double> *Whad_gen_eta        = new std::vector<double> ();

		std::vector<double> *gen_particle_pt	= new std::vector<double> ();
		std::vector<double> *gen_particle_Id	= new std::vector<double> ();
		std::vector<double> *gen_particle_phi	= new std::vector<double> ();
		std::vector<double> *gen_particle_eta	= new std::vector<double> ();
		std::vector<double> *gen_particle_status= new std::vector<double> ();


		std::vector<std::string> HLT_list = {"HLT_IsoMu24_v", "HLT_Mu50_v", "HLT_Ele27_WPTight_Gsf_v", "HLT_DoubleMu38NoFiltersNoVtx_v", "HLT_DoubleEle33_CaloIdL_MW_v"};

		std::vector<double> *genP_pass		= new std::vector<double> ();

		int count_gen = 0;
                int count_gen_lep = 0;
                int count_gen_q_b = 0;
                int count_gen_lq_g = 0;
                int count_gen_lq_w = 0;
                int count_gen_g = 0;
                int count_gen_others = 0;
};



MakeNTuple::MakeNTuple(const edm::ParameterSet& iConfig):
	jetsAK8Token (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK8")))
	, jetsAK4Token (consumes<pat::JetCollection>(iConfig.getParameter<edm::InputTag>("jetsAK4")))
	, genJetsToken (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJets")))
	, genJetsTokenAK8 (consumes<reco::GenJetCollection>(iConfig.getParameter<edm::InputTag>("genJetsAK8")))
	, muonsToken (consumes<pat::MuonCollection>(iConfig.getParameter<edm::InputTag>("muons")))
	, electronsToken (consumes<pat::ElectronCollection>(iConfig.getParameter<edm::InputTag>("electrons")))
	, MetToken (consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("MET")))
	, PFCandToken (consumes<std::vector<pat::PackedCandidate>>(iConfig.getParameter<edm::InputTag>("PFCand")))
	, PileupSumInfoInputTag (consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupSumInfoInputTag")))
	, verticesToken (consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices")))
//	, genParticlesToken_(consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("genP")))
	, recoProtonsSingleRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonSingleRPTag" ) ) )
	, recoProtonsMultiRPToken_   ( consumes<std::vector<reco::ForwardProton> >      ( iConfig.getParameter<edm::InputTag>( "ppsRecoProtonMultiRPTag" ) ) )
	, LumiWeights_(0)
	, MC(iConfig.getParameter<bool>("MC"))
	, MC_Signal(iConfig.getParameter<bool>("MC_Signal"))
	, DATA(iConfig.getParameter<bool>("DATA"))
	, triggerResultsToken_(consumes<edm::TriggerResults>              (iConfig.getParameter<edm::InputTag>("TriggerResults")))
	, hltPrescaleProvider_(iConfig, consumesCollector(), *this)
{
	//now do what ever initialization is needed
	usesResource("TFileService");
	edm::Service<TFileService> fs;
	EventBranchs                            = fs->make<TTree>( "Events","Events" );
	EventBranchs->Branch("Run", &Run, "Run/I");
	EventBranchs->Branch("LumiSection", &LumiSection, "LumiSection/I");
	EventBranchs->Branch("BX", &BX, "BX/I");
	EventBranchs->Branch("xangle", &xangle, "xangle/I");

	EventBranchs->Branch("EventNum", &EventNum, "EventNum/I");
	EventBranchs->Branch("nVtx", &nVtx, "nVtx/I");

	EventBranchs->Branch("PUWeight",&PUWeight,"PUWeight/D");
	EventBranchs->Branch("nPU",&nPU,"nPU/D");


	EventBranchs->Branch("vtx_z",&vtx_z,"vtx_z/D");
	EventBranchs->Branch("vtx_isValid",&vtx_isValid,"vtx_isValid/B");
	EventBranchs->Branch("vtx_isFake",&vtx_isFake,"vtx_isFake/B");


	EventBranchs->Branch("PUWeight",&PUWeight,"PUWeight/D");
	EventBranchs->Branch("nPU",&nPU,"nPU/D");

	const edm::InputTag& PileupSumInfoInputTags = edm::InputTag( "slimmedAddPileupInfo" ) ;

	LumiWeights_ = new edm::LumiReWeighting("PileupMC.root", "MyDataPileupHistogram.root", "input_Event/N_TrueInteractions", "pileup", PileupSumInfoInputTags);

	EventBranchs->Branch("ArmF_xi_gen","std::vector<double>",&ArmF_xi_gen);
	EventBranchs->Branch("ArmB_xi_gen","std::vector<double>",&ArmB_xi_gen);
	EventBranchs->Branch("ArmF_t_gen","std::vector<double>",&ArmF_t_gen);
	EventBranchs->Branch("ArmB_t_gen","std::vector<double>",&ArmB_t_gen);
	EventBranchs->Branch("ArmF_thx_gen","std::vector<double>",&ArmF_thx_gen);
	EventBranchs->Branch("ArmB_thx_gen","std::vector<double>",&ArmB_thx_gen);
	EventBranchs->Branch("ArmF_thy_gen","std::vector<double>",&ArmF_thy_gen);
	EventBranchs->Branch("ArmB_thy_gen","std::vector<double>",&ArmB_thy_gen);

	EventBranchs->Branch("ProtCand_xi","std::vector<double>",&ProtCand_xi);
	EventBranchs->Branch("ProtCand_t","std::vector<double>",&ProtCand_t);
	EventBranchs->Branch("ProtCand_ThX","std::vector<double>",&ProtCand_ThX);
	EventBranchs->Branch("ProtCand_ThY","std::vector<double>",&ProtCand_ThY);
	EventBranchs->Branch("ProtCand_rpid","std::vector<double>",&ProtCand_rpid);
	EventBranchs->Branch("ProtCand_arm","std::vector<double>",&ProtCand_arm);
	EventBranchs->Branch("ProtCand_ismultirp","std::vector<double>",&ProtCand_ismultirp);

	EventBranchs->Branch("muon_px","std::vector<double>",&muon_px);
	EventBranchs->Branch("muon_py","std::vector<double>",&muon_py);
	EventBranchs->Branch("muon_pz","std::vector<double>",&muon_pz);
	EventBranchs->Branch("muon_pt","std::vector<double>",&muon_pt);
	EventBranchs->Branch("muon_E","std::vector<double>",&muon_E);
	EventBranchs->Branch("muon_vtxZ","std::vector<double>",&muon_vtxZ);
	EventBranchs->Branch("muon_phi","std::vector<double>",&muon_phi);
	EventBranchs->Branch("muon_eta","std::vector<double>",&muon_eta);
	EventBranchs->Branch("muon_PFBasedIso","std::vector<double>",&muon_PFBasedIso);
	EventBranchs->Branch("muon_TrackBasedIso","std::vector<double>",&muon_TrackBasedIso);
	EventBranchs->Branch("muon_isTightMuon","std::vector<bool>",&muon_isTightMuon);
	EventBranchs->Branch("muon_isMediumMuon","std::vector<bool>",&muon_isMediumMuon);
	EventBranchs->Branch("muon_isLooseMuon","std::vector<bool>",&muon_isLooseMuon);
	EventBranchs->Branch("muon_isHighPtMuon","std::vector<bool>",&muon_isHighPtMuon);

	EventBranchs->Branch("electron_px","std::vector<double>",&electron_px);
	EventBranchs->Branch("electron_py","std::vector<double>",&electron_py);
	EventBranchs->Branch("electron_pz","std::vector<double>",&electron_pz);
	EventBranchs->Branch("electron_pt","std::vector<double>",&electron_pt);
	EventBranchs->Branch("electron_E","std::vector<double>",&electron_E);
	EventBranchs->Branch("electron_vtxZ","std::vector<double>",&electron_vtxZ);
	EventBranchs->Branch("electron_phi","std::vector<double>",&electron_phi);
	EventBranchs->Branch("electron_eta","std::vector<double>",&electron_eta);
	EventBranchs->Branch("electron_isTightElectron","std::vector<bool>",&electron_isTightElectron);
	EventBranchs->Branch("electron_isMediumElectron","std::vector<bool>",&electron_isMediumElectron);
	EventBranchs->Branch("electron_isLooseElectron","std::vector<bool>",&electron_isLooseElectron);
	EventBranchs->Branch("electron_isVetoElectron","std::vector<bool>",&electron_isVetoElectron);


	EventBranchs->Branch("jetAK4_px","std::vector<double>",&jetAK4_px);
	EventBranchs->Branch("jetAK4_py","std::vector<double>",&jetAK4_py);
	EventBranchs->Branch("jetAK4_pz","std::vector<double>",&jetAK4_pz);
	EventBranchs->Branch("jetAK4_pt","std::vector<double>",&jetAK4_pt);
	EventBranchs->Branch("jetAK4_E","std::vector<double>",&jetAK4_E);
	EventBranchs->Branch("jetAK4_phi","std::vector<double>",&jetAK4_phi);
	EventBranchs->Branch("jetAK4_eta","std::vector<double>",&jetAK4_eta);
	EventBranchs->Branch("jetAK4_btag","std::vector<double>",&jetAK4_btag);
	EventBranchs->Branch("jetAK4_isLoose","std::vector<bool>",&jetAK4_isLoose);
	EventBranchs->Branch("jetAK4_isTight","std::vector<bool>",&jetAK4_isTight);

	EventBranchs->Branch("genJets_px","std::vector<double>",&genJets_px);
	EventBranchs->Branch("genJets_py","std::vector<double>",&genJets_py);
	EventBranchs->Branch("genJets_pz","std::vector<double>",&genJets_pz);
	EventBranchs->Branch("genJets_pt","std::vector<double>",&genJets_pt);
	EventBranchs->Branch("genJets_E","std::vector<double>",&genJets_E);
	EventBranchs->Branch("genJets_phi","std::vector<double>",&genJets_phi);
	EventBranchs->Branch("genJets_eta","std::vector<double>",&genJets_eta);


	EventBranchs->Branch("jetAK8_px","std::vector<double>",&jetAK8_px);
	EventBranchs->Branch("jetAK8_py","std::vector<double>",&jetAK8_py);
	EventBranchs->Branch("jetAK8_pz","std::vector<double>",&jetAK8_pz);
	EventBranchs->Branch("jetAK8_pt","std::vector<double>",&jetAK8_pt);
	EventBranchs->Branch("jetAK8_E","std::vector<double>",&jetAK8_E);
	EventBranchs->Branch("jetAK8_phi","std::vector<double>",&jetAK8_phi);
	EventBranchs->Branch("jetAK8_eta","std::vector<double>",&jetAK8_eta);
	EventBranchs->Branch("jetAK8_btag","std::vector<double>",&jetAK8_btag);
	EventBranchs->Branch("jetAK8_isLoose","std::vector<bool>",&jetAK8_isLoose);
	EventBranchs->Branch("jetAK8_isTight","std::vector<bool>",&jetAK8_isTight);
	EventBranchs->Branch("jetAK8_prunedMass","std::vector<double>",&jetAK8_prunedMass);
	EventBranchs->Branch("jetAK8_tau21","std::vector<double>",&jetAK8_tau21);

	EventBranchs->Branch("pfphi","std::vector<double>",&pfphi);
	EventBranchs->Branch("pfeta","std::vector<double>",&pfeta);

	EventBranchs->Branch("pffromPV","std::vector<int>",&pffromPV);
	EventBranchs->Branch("pfdz","std::vector<double>",&pfdz);
	EventBranchs->Branch("pfpt","std::vector<double>",&pfpt);


	EventBranchs->Branch("METPx", &METPx, "METPx/D");
	EventBranchs->Branch("METPy", &METPy, "METPy/D");
	EventBranchs->Branch("METPt", &METPt, "METPt/D");
	EventBranchs->Branch("METphi", &METphi, "METphi/D");

	EventBranchs->Branch("HLT_name","std::vector<std::string>",&HLT_name);
	EventBranchs->Branch("HLT_pass","std::vector<bool>",&HLT_pass);
	EventBranchs->Branch("HLT_prescale","std::vector<int>",&HLT_prescale);

	EventBranchs->Branch("muon_gen_px","std::vector<double>",&muon_gen_px);
	EventBranchs->Branch("muon_gen_py","std::vector<double>",&muon_gen_py);
	EventBranchs->Branch("muon_gen_pz","std::vector<double>",&muon_gen_pz);
	EventBranchs->Branch("muon_gen_pt","std::vector<double>",&muon_gen_pt);
	EventBranchs->Branch("muon_gen_E","std::vector<double>",&muon_gen_E);
	EventBranchs->Branch("muon_gen_phi","std::vector<double>",&muon_gen_phi);
	EventBranchs->Branch("muon_gen_eta","std::vector<double>",&muon_gen_eta);

	EventBranchs->Branch("neut_gen_px","std::vector<double>",&neut_gen_px);
	EventBranchs->Branch("neut_gen_py","std::vector<double>",&neut_gen_py);
	EventBranchs->Branch("neut_gen_pz","std::vector<double>",&neut_gen_pz);
	EventBranchs->Branch("neut_gen_pt","std::vector<double>",&neut_gen_pt);
	EventBranchs->Branch("neut_gen_E","std::vector<double>",&neut_gen_E);
	EventBranchs->Branch("neut_gen_phi","std::vector<double>",&neut_gen_phi);
	EventBranchs->Branch("neut_gen_eta","std::vector<double>",&neut_gen_eta);

	EventBranchs->Branch("ele_gen_px","std::vector<double>",&ele_gen_px);
	EventBranchs->Branch("ele_gen_py","std::vector<double>",&ele_gen_py);
	EventBranchs->Branch("ele_gen_pz","std::vector<double>",&ele_gen_pz);
	EventBranchs->Branch("ele_gen_pt","std::vector<double>",&ele_gen_pt);
	EventBranchs->Branch("ele_gen_E","std::vector<double>",&ele_gen_E);
	EventBranchs->Branch("ele_gen_phi","std::vector<double>",&ele_gen_phi);
	EventBranchs->Branch("ele_gen_eta","std::vector<double>",&ele_gen_eta);

	EventBranchs->Branch("qrk_gen_px","std::vector<double>",&qrk_gen_px);
	EventBranchs->Branch("qrk_gen_py","std::vector<double>",&qrk_gen_py);
	EventBranchs->Branch("qrk_gen_pz","std::vector<double>",&qrk_gen_pz);
	EventBranchs->Branch("qrk_gen_pt","std::vector<double>",&qrk_gen_pt);
	EventBranchs->Branch("qrk_gen_E","std::vector<double>",&qrk_gen_E);
	EventBranchs->Branch("qrk_gen_phi","std::vector<double>",&qrk_gen_phi);
	EventBranchs->Branch("qrk_gen_eta","std::vector<double>",&qrk_gen_eta);

	EventBranchs->Branch("Wlep_gen_px","std::vector<double>",&Wlep_gen_px);
	EventBranchs->Branch("Wlep_gen_py","std::vector<double>",&Wlep_gen_py);
	EventBranchs->Branch("Wlep_gen_pz","std::vector<double>",&Wlep_gen_pz);
	EventBranchs->Branch("Wlep_gen_pt","std::vector<double>",&Wlep_gen_pt);
	EventBranchs->Branch("Wlep_gen_E","std::vector<double>",&Wlep_gen_E);
	EventBranchs->Branch("Wlep_gen_phi","std::vector<double>",&Wlep_gen_phi);
	EventBranchs->Branch("Wlep_gen_eta","std::vector<double>",&Wlep_gen_eta);

	EventBranchs->Branch("Whad_gen_px","std::vector<double>",&Whad_gen_px);
	EventBranchs->Branch("Whad_gen_py","std::vector<double>",&Whad_gen_py);
	EventBranchs->Branch("Whad_gen_pz","std::vector<double>",&Whad_gen_pz);
	EventBranchs->Branch("Whad_gen_pt","std::vector<double>",&Whad_gen_pt);
	EventBranchs->Branch("Whad_gen_E","std::vector<double>",&Whad_gen_E);
	EventBranchs->Branch("Whad_gen_phi","std::vector<double>",&Whad_gen_phi);
	EventBranchs->Branch("Whad_gen_eta","std::vector<double>",&Whad_gen_eta);

	EventBranchs->Branch("gen_particle_pt","std::vector<double>",&gen_particle_pt);
	EventBranchs->Branch("gen_particle_Id","std::vector<double>",&gen_particle_Id);
	EventBranchs->Branch("gen_particle_phi","std::vector<double>",&gen_particle_phi);
	EventBranchs->Branch("gen_particle_eta","std::vector<double>",&gen_particle_eta);
	EventBranchs->Branch("gen_particle_status","std::vector<double>",&gen_particle_status);

}


MakeNTuple::~MakeNTuple()
{

	// do anything here that needs to be done at desctruction time
	// (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
	void
MakeNTuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
	ArmF_xi_gen->clear();
	ArmB_xi_gen->clear();
	ArmF_t_gen->clear();
	ArmB_t_gen->clear();
	ArmF_thx_gen->clear();
	ArmB_thx_gen->clear();
	ArmF_thy_gen->clear();
	ArmB_thy_gen->clear();

	ProtCand_xi->clear();
	ProtCand_t->clear();
	ProtCand_ThX->clear();
	ProtCand_ThY->clear();
	ProtCand_rpid->clear();
	ProtCand_arm->clear();
	ProtCand_ismultirp->clear();

	muon_px->clear();
	muon_py->clear();
	muon_pz->clear();
	muon_pt->clear();
	muon_E->clear();
	muon_vtxZ->clear();
	muon_phi->clear();
	muon_eta->clear();
	muon_PFBasedIso->clear();
	muon_TrackBasedIso->clear();
	muon_isTightMuon->clear();
	muon_isMediumMuon->clear();
	muon_isLooseMuon->clear();
	muon_isHighPtMuon->clear();
	electron_px->clear();
	electron_py->clear();
	electron_pz->clear();
	electron_pt->clear();
	electron_E->clear();
	electron_vtxZ->clear();
	electron_phi->clear();
	electron_eta->clear();
	electron_isTightElectron->clear();
	electron_isMediumElectron->clear();
	electron_isLooseElectron->clear();
	electron_isVetoElectron->clear();
	jetAK4_px->clear();
	jetAK4_py->clear();
	jetAK4_pz->clear();
	jetAK4_pt->clear();
	jetAK4_E->clear();
	jetAK4_phi->clear();
	jetAK4_eta->clear();
	jetAK4_btag->clear();
	jetAK4_isLoose->clear();
	jetAK4_isTight->clear();
	genJets_px->clear();
	genJets_py->clear();
	genJets_pz->clear();
	genJets_pt->clear();
	genJets_E->clear();
	genJets_phi->clear();
	genJets_eta->clear();
	jetAK8_px->clear();
	jetAK8_py->clear();
	jetAK8_pz->clear();
	jetAK8_pt->clear();
	jetAK8_E->clear();
	jetAK8_phi->clear();
	jetAK8_eta->clear();
	jetAK8_btag->clear();
	jetAK8_isLoose->clear();
	jetAK8_isTight->clear();
	jetAK8_prunedMass->clear();
	jetAK8_tau21->clear();
	genJetsAK8_px->clear();
	genJetsAK8_py->clear();
	genJetsAK8_pz->clear();
	genJetsAK8_pt->clear();
	genJetsAK8_E->clear();
	genJetsAK8_phi->clear();
	genJetsAK8_eta->clear();
	pfphi->clear();
	pfeta->clear();
	pffromPV->clear();
	pfdz->clear();
	pfpt->clear();
	HLT_name->clear();
	HLT_pass->clear();
	HLT_prescale->clear();
	muon_gen_px->clear();
	muon_gen_py->clear();
	muon_gen_pz->clear();
	muon_gen_pt->clear();
	muon_gen_E->clear();
	muon_gen_phi->clear();
	muon_gen_eta->clear();
	neut_gen_px->clear();
	neut_gen_py->clear();
	neut_gen_pz->clear();
	neut_gen_pt->clear();
	neut_gen_E->clear();
	neut_gen_phi->clear();
	neut_gen_eta->clear();
	ele_gen_px->clear();
	ele_gen_py->clear();
	ele_gen_pz->clear();
	ele_gen_pt->clear();
	ele_gen_E->clear();
	ele_gen_phi->clear();
	ele_gen_eta->clear();
	qrk_gen_px->clear();
	qrk_gen_py->clear();
	qrk_gen_pz->clear();
	qrk_gen_pt->clear();
	qrk_gen_E->clear();
	qrk_gen_phi->clear();
	qrk_gen_eta->clear();
	Wlep_gen_px->clear();
	Wlep_gen_py->clear();
	Wlep_gen_pz->clear();
	Wlep_gen_pt->clear();
	Wlep_gen_E->clear();
	Wlep_gen_phi->clear();
	Wlep_gen_eta->clear();
	Whad_gen_px->clear();
	Whad_gen_py->clear();
	Whad_gen_pz->clear();
	Whad_gen_pt->clear();
	Whad_gen_E->clear();
	Whad_gen_phi->clear();
	Whad_gen_eta->clear();

	gen_particle_pt->clear();
	gen_particle_Id->clear();
	gen_particle_phi->clear();
	gen_particle_eta->clear();
	gen_particle_status->clear();

	using namespace edm;
	using namespace std;

	iEvent.getByToken(jetsAK8Token, jetsAK8);
	iEvent.getByToken(jetsAK4Token, jetsAK4);
	if(MC){
		iEvent.getByToken(genJetsToken, genjets);
		iEvent.getByToken(genJetsTokenAK8, genjetsAK8);
	}

	iEvent.getByToken(muonsToken, muons);
	iEvent.getByToken(electronsToken, electrons);
	iEvent.getByToken(MetToken, MET);
	iEvent.getByToken(PFCandToken, PFCand);
	iEvent.getByToken(verticesToken, vertices);


	if(MC){
		iEvent.getByToken(PileupSumInfoInputTag, PupInfo);
	}

	//	if(MC_Signal || DATA){
	edm::Handle<vector<reco::ForwardProton>> recoMultiRPProtons;
	iEvent.getByToken(recoProtonsMultiRPToken_, recoMultiRPProtons);
	edm::Handle<vector<reco::ForwardProton>> recoSingleRPProtons;
	iEvent.getByToken(recoProtonsSingleRPToken_, recoSingleRPProtons);
	//	}

	// Event Info
	BX = iEvent.bunchCrossing();
	Run = iEvent.id().run();
	LumiSection = iEvent.luminosityBlock();
	EventNum = iEvent.id().event();

	nVtx = vertices->size();

	vtx_isValid = false;
	vtx_isFake = true;
	vtx_z = -999;

	if (nVtx>0){
		vtx_isValid = vertices->at(0).isValid();
		vtx_isFake = vertices->at(0).isFake();
		vtx_z = vertices->at(0).z();
	}

	// X A N G L E   I N F O
	/*	if(DATA){
		edm::ESHandle<LHCInfo> pSetup;
		const string label = "";
		iSetup.get<LHCInfoRcd>().get(label, pSetup);
		const LHCInfo* pInfo = pSetup.product();
		std::cout << pInfo->crossingAngle() << std::endl;
		xangle = pInfo->crossingAngle();
		}
		*/


	// F O R  P I L E U P  R E W E I G T I N G
	if(MC){
		const edm::EventBase* iEventB = dynamic_cast<const edm::EventBase*>(&iEvent);
		if (LumiWeights_) PUWeight = LumiWeights_->weight( (*iEventB) );
		std::vector<PileupSummaryInfo>::const_iterator PVI;
		int npv = -1;
		for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
			int BXs = PVI->getBunchCrossing();
			if(BXs == 0) {
				npv = PVI->getTrueNumInteractions();
				continue;
			}
		}
		nPU = npv;
	}

	// F O R  T R I G G E R  I N F O
	if(!MC_Signal){
		edm::Handle<edm::TriggerResults> hltResults;
		iEvent.getByToken(triggerResultsToken_, hltResults);
		const edm::TriggerNames& trigNames = iEvent.triggerNames(*hltResults);


		for (size_t i=0; i<trigNames.size(); i++) {
			for (size_t j =0; j<HLT_list.size(); j++) {
				if (trigNames.triggerNames().at(i).find(HLT_list.at(j)) != std::string::npos){
					HLT_name->push_back(trigNames.triggerNames().at(i));
//					std::cout << "trigNames.triggerNames().at(i): " << trigNames.triggerNames().at(i) << std::endl;
					HLT_pass->push_back(hltResults->accept(i));
					HLT_prescale->push_back(hltPrescaleProvider_.prescaleValue(iEvent, iSetup,trigNames.triggerName(i)));
				}
			}
		}
	}



	if (DATA || MC_Signal){
		for(size_t i = 0;i<recoMultiRPProtons->size();i++){
			if(recoMultiRPProtons->at(i).validFit()){
				CTPPSDetId rpId((*recoMultiRPProtons->at(i).contributingLocalTracks().begin())->getRPId());
				//			int decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
				int armId = rpId.arm();
				ProtCand_xi->push_back(recoMultiRPProtons->at(i).xi());
				//			ProtCand_t->push_back(recoMultiRPProtons->at(i).t());
				ProtCand_ThX->push_back(recoMultiRPProtons->at(i).thetaX());
				ProtCand_ThY->push_back(recoMultiRPProtons->at(i).thetaY());
				ProtCand_rpid->push_back(-999);
				ProtCand_arm->push_back(armId);
				ProtCand_ismultirp->push_back(1);
			}
		}

		for(size_t i = 0;i<recoSingleRPProtons->size();i++){
			if(recoSingleRPProtons->at(i).validFit()){
				CTPPSDetId rpId((*recoSingleRPProtons->at(i).contributingLocalTracks().begin())->getRPId());
				int decRPId = rpId.arm()*100 + rpId.station()*10 + rpId.rp();
				//                      int armId = rpId.arm();
				ProtCand_xi->push_back(recoSingleRPProtons->at(i).xi());
				//                      ProtCand_t->push_back(recoSingleRPProtons->at(i).t());
				ProtCand_ThX->push_back(recoSingleRPProtons->at(i).thetaX());
				ProtCand_ThY->push_back(recoSingleRPProtons->at(i).thetaY());
				ProtCand_rpid->push_back(decRPId);
				ProtCand_arm->push_back(-999);
				ProtCand_ismultirp->push_back(0);
			}
		}
	}


	/*
	   if (MC && !MC_Signal){
	   double strips_xi_arm0_N;
	   double strips_xi_arm0_F;
	   double strips_xi_arm1_N;
	   double strips_xi_arm1_F;
	   float a, b, c, d;


	   Int_t ncols;
	   Int_t nlines = 0;
	   FILE *fp = fopen("RPs_xi.txt","r");

	   int rand_no = rand() % 557673;

	   while (nlines<rand_no+1) {

	   ncols = fscanf(fp,"%f %f %f %f",&a, &b, &c, &d);


	   if (ncols < 0) break;

	   if (nlines==rand_no){
	   strips_xi_arm0_N = a;
	   strips_xi_arm0_F = b;
	   strips_xi_arm1_N = c;
	   strips_xi_arm1_F = d;
	   }

	   nlines++;

	   }

	   fclose(fp);
	   if (strips_xi_arm0_N>0.005) {ArmF_N_xi->push_back(strips_xi_arm0_N);}
	   if (strips_xi_arm0_F>0.005) {ArmF_F_xi->push_back(strips_xi_arm0_F);}
	   if (strips_xi_arm1_N>0.005) {ArmB_N_xi->push_back(strips_xi_arm1_N);}
	   if (strips_xi_arm1_F>0.005) {ArmB_F_xi->push_back(strips_xi_arm1_F);}

	   }
	   */
	double AK8NHF, AK8NEMF, AK8NumConst, AK8CHF, AK8CHM, AK8CEMF;
	bool AK8tightJetID, AK8looseJetID;
	for(size_t i = 0;i<jetsAK8->size();i++){
		if(jetsAK8->at(i).pt()>170 && abs(jetsAK8->at(i).eta())<2.4){
			AK8NHF  = jetsAK8->at(i).neutralHadronEnergyFraction();
			AK8NEMF = jetsAK8->at(i).neutralEmEnergyFraction();
			AK8CHF  = jetsAK8->at(i).chargedHadronEnergyFraction();
			AK8CEMF = jetsAK8->at(i).chargedEmEnergyFraction();
			AK8NumConst = jetsAK8->at(i).chargedMultiplicity()+jetsAK8->at(i).neutralMultiplicity();
			AK8CHM      = jetsAK8->at(i).chargedMultiplicity();
			AK8tightJetID = (AK8NHF<0.90 && AK8NEMF<0.90 && AK8NumConst>1) && ((abs(jetsAK8->at(i).eta())<=2.4 && AK8CHF>0 && AK8CHM>0 && AK8CEMF<0.99) || abs(jetsAK8->at(i).eta())>2.4) && abs(jetsAK8->at(i).eta())<=2.7;
			AK8looseJetID = (AK8NHF<0.99 && AK8NEMF<0.99 && AK8NumConst>1) && ((abs(jetsAK8->at(i).eta())<=2.4 && AK8CHF>0 && AK8CHM>0 && AK8CEMF<0.99) || abs(jetsAK8->at(i).eta())>2.4) && abs(jetsAK8->at(i).eta())<=2.7 ;
			jetAK8_px->push_back(jetsAK8->at(i).px());
			jetAK8_py->push_back(jetsAK8->at(i).py());
			jetAK8_pz->push_back(jetsAK8->at(i).pz());
			jetAK8_pt->push_back(jetsAK8->at(i).pt());
			jetAK8_E->push_back(jetsAK8->at(i).energy());
			jetAK8_phi->push_back(jetsAK8->at(i).phi());
			jetAK8_eta->push_back(jetsAK8->at(i).eta());
			jetAK8_btag->push_back(jetsAK8->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
			jetAK8_isLoose->push_back(AK8tightJetID);
			jetAK8_isTight->push_back(AK8looseJetID);
			jetAK8_prunedMass->push_back(jetsAK8->at(i).userFloat("ak8PFJetsCHSPrunedMass"));
			jetAK8_tau21->push_back((jetsAK8->at(i).userFloat("NjettinessAK8CHS:tau2"))/(jetsAK8->at(i).userFloat("NjettinessAK8CHS:tau1")));
			//			std::cout << "jetsAK8->at(i).pt(): " << jetsAK8->at(i).pt() << "; jetsAK8->at(i).phi(): " << jetsAK8->at(i).phi() << "; jetsAK8->at(i).eta(): " << jetsAK8->at(i).eta() << std::endl;
		}
	}

	double NHF, NEMF, NumConst, CHF, CHM, CEMF;
	bool tightJetID, looseJetID;
	for(size_t i = 0;i<jetsAK4->size();i++){
		if(jetsAK4->at(i).pt()>20 && abs(jetsAK4->at(i).eta())<2.4){
			NHF  = jetsAK4->at(i).neutralHadronEnergyFraction();
			NEMF = jetsAK4->at(i).neutralEmEnergyFraction();
			CHF  = jetsAK4->at(i).chargedHadronEnergyFraction();
			CEMF = jetsAK4->at(i).chargedEmEnergyFraction();
			NumConst = jetsAK4->at(i).chargedMultiplicity()+jetsAK4->at(i).neutralMultiplicity();
			CHM      = jetsAK4->at(i).chargedMultiplicity(); 
			tightJetID = (NHF<0.90 && NEMF<0.90 && NumConst>1) && ((abs(jetsAK4->at(i).eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jetsAK4->at(i).eta())>2.4) && abs(jetsAK4->at(i).eta())<=2.7;
			looseJetID = (NHF<0.99 && NEMF<0.99 && NumConst>1) && ((abs(jetsAK4->at(i).eta())<=2.4 && CHF>0 && CHM>0 && CEMF<0.99) || abs(jetsAK4->at(i).eta())>2.4) && abs(jetsAK4->at(i).eta())<=2.7 ;
			jetAK4_px->push_back(jetsAK4->at(i).px());
			jetAK4_py->push_back(jetsAK4->at(i).py());
			jetAK4_pz->push_back(jetsAK4->at(i).pz());
			jetAK4_pt->push_back(jetsAK4->at(i).pt());
			jetAK4_E->push_back(jetsAK4->at(i).energy());
			jetAK4_phi->push_back(jetsAK4->at(i).phi());
			jetAK4_eta->push_back(jetsAK4->at(i).eta());
			jetAK4_btag->push_back(jetsAK4->at(i).bDiscriminator("pfCombinedInclusiveSecondaryVertexV2BJetTags"));
			jetAK4_isLoose->push_back(tightJetID);
			jetAK4_isTight->push_back(looseJetID);
		}
	}

	if(MC){
		for(size_t i = 0;i<genjets->size();i++){
			if(genjets->at(i).pt()>20 && abs(genjets->at(i).eta())<2.4){
				genJets_px->push_back(genjets->at(i).px());
				genJets_py->push_back(genjets->at(i).py());
				genJets_pz->push_back(genjets->at(i).pz());
				genJets_pt->push_back(genjets->at(i).pt());
				genJets_E->push_back(genjets->at(i).energy());
				genJets_phi->push_back(genjets->at(i).phi());
				genJets_eta->push_back(genjets->at(i).eta());
			}
		}
		for(size_t i = 0;i<genjetsAK8->size();i++){
			if(genjetsAK8->at(i).pt()>170 && abs(genjetsAK8->at(i).eta())<2.4){
				genJets_px->push_back(genjetsAK8->at(i).px());
				genJets_py->push_back(genjetsAK8->at(i).py());
				genJets_pz->push_back(genjetsAK8->at(i).pz());
				genJets_pt->push_back(genjetsAK8->at(i).pt());
				genJets_E->push_back(genjetsAK8->at(i).energy());
				genJets_phi->push_back(genjetsAK8->at(i).phi());
				genJets_eta->push_back(genjetsAK8->at(i).eta());
			}
		}
	}


	for(size_t i = 0;i<muons->size();i++){
		if(muons->at(i).pt()>20 && abs(muons->at(i).eta())<2.4){
			muon_px->push_back(muons->at(i).px());
			muon_py->push_back(muons->at(i).py());
			muon_pz->push_back(muons->at(i).pz());
			muon_pt->push_back(muons->at(i).pt());
			muon_E->push_back(muons->at(i).energy());
			muon_vtxZ->push_back(muons->at(i).vertex().z());
			muon_phi->push_back(muons->at(i).phi());
			muon_eta->push_back(muons->at(i).eta());
			muon_PFBasedIso->push_back((muons->at(i).pfIsolationR04().sumChargedHadronPt + max(0., muons->at(i).pfIsolationR04().sumNeutralHadronEt + muons->at(i).pfIsolationR04().sumPhotonEt - 0.5*muons->at(i).pfIsolationR04().sumPUPt))/muons->at(i).pt());
			muon_TrackBasedIso->push_back((muons->at(i).isolationR03().sumPt)/(muons->at(i).pt()));
			muon_isTightMuon->push_back(muon::isTightMuon(muons->at(i), vertices->at(0)));
			muon_isMediumMuon->push_back(muon::isMediumMuon(muons->at(i)));
			muon_isLooseMuon->push_back(muon::isLooseMuon(muons->at(i)));
			muon_isHighPtMuon->push_back(muon::isHighPtMuon(muons->at(i), vertices->at(0)));
		}
	}

	for(size_t i = 0;i<electrons->size();i++){
		if(electrons->at(i).pt()>20 && abs(electrons->at(i).eta())<2.4){
			electron_px->push_back(electrons->at(i).px());
			electron_py->push_back(electrons->at(i).py());
			electron_pz->push_back(electrons->at(i).pz());
			electron_pt->push_back(electrons->at(i).pt());
			electron_E->push_back(electrons->at(i).energy());
			electron_vtxZ->push_back(electrons->at(i).vertex().z());
			electron_phi->push_back(electrons->at(i).phi());
			electron_eta->push_back(electrons->at(i).eta());
			electron_isTightElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Summer16-80X-V1-tight"));
			electron_isMediumElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Summer16-80X-V1-medium"));
			electron_isLooseElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Summer16-80X-V1-loose"));
			electron_isVetoElectron->push_back(electrons->at(i).electronID("cutBasedElectronID-Summer16-80X-V1-veto"));
		}
	}

	METPx = (MET->front()).px();
	METPy = (MET->front()).py();
	METPt = (MET->front()).pt();
	METphi = (MET->front()).phi();

	int npfVtx = 0;
	for (size_t pf = 0; pf < PFCand->size(); pf++) {
		if (abs(PFCand->at(pf).pdgId()) == 211 || abs(PFCand->at(pf).pdgId()) == 11 || abs(PFCand->at(pf).pdgId()) == 13){
//			if (PFCand->at(pf).fromPV(0)>1) {
//				if (abs(PFCand->at(pf).dz())<0.15){
					npfVtx++;
					pfphi->push_back(PFCand->at(pf).phiAtVtx());
					pfeta->push_back(PFCand->at(pf).eta());
					pffromPV->push_back(PFCand->at(pf).fromPV(0));
					pfdz->push_back(PFCand->at(pf).dz());
					pfpt->push_back(PFCand->at(pf).pt());
//				}
//			}
		}
	}

	// G E N   I N F O
/*	iEvent.getByToken(genParticlesToken_, genP);
	if(MC){
		for(size_t i=0; i< genP->size(); i++){
			const GenParticle & p = (*genP)[i];
			//			if(p.pt()>150) std::cout << "p.pdgId(): " << p.pdgId() << ";  p.pt(): " << p.pt() << ";  p.status(): " << p.status() << ";  p.phi(): " << p.phi() << ";  p.eta(): " << p.eta() << std::endl;

			//			const Candidate * mom = p.mother();
			//			std::cout << mom->pdgId() << std::endl;
			if(p.pt()>130){
				gen_particle_pt->push_back(p.pt());
				gen_particle_Id->push_back(p.pdgId());
				gen_particle_phi->push_back(p.phi());
				gen_particle_eta->push_back(p.eta());
				gen_particle_status->push_back(p.status());
			}

			if(fabs(p.pdgId()) == 24 && p.numberOfDaughters()==2){
				if(fabs(p.daughter(0)->pdgId())<7 && fabs(p.daughter(0)->pdgId()>0)){
					Whad_gen_px->push_back(p.px());
					Whad_gen_py->push_back(p.py());
					Whad_gen_pz->push_back(p.pz());
					Whad_gen_pt->push_back(p.pt());
					Whad_gen_E->push_back(p.energy());
					Whad_gen_phi->push_back(p.phi());
					Whad_gen_eta->push_back(p.eta());
				}
				if(fabs(p.daughter(0)->pdgId())<15 && fabs(p.daughter(0)->pdgId()>10)){
					Wlep_gen_px->push_back(p.px());
					Wlep_gen_py->push_back(p.py());
					Wlep_gen_pz->push_back(p.pz());
					Wlep_gen_pt->push_back(p.pt());
					Wlep_gen_E->push_back(p.energy());
					Wlep_gen_phi->push_back(p.phi());
					Wlep_gen_eta->push_back(p.eta());
				}
				if(fabs(p.daughter(0)->pdgId())<7 && fabs(p.daughter(0)->pdgId()>0)){
					qrk_gen_px->push_back(p.daughter(0)->px());
					qrk_gen_py->push_back(p.daughter(0)->py());
					qrk_gen_pz->push_back(p.daughter(0)->pz());
					qrk_gen_pt->push_back(p.daughter(0)->pt());
					qrk_gen_E->push_back(p.daughter(0)->energy());
					qrk_gen_phi->push_back(p.daughter(0)->phi());
					qrk_gen_eta->push_back(p.daughter(0)->eta());
				}
				if(fabs(p.daughter(1)->pdgId())<7 && fabs(p.daughter(1)->pdgId()>0)){
					qrk_gen_px->push_back(p.daughter(1)->px());
					qrk_gen_py->push_back(p.daughter(1)->py());
					qrk_gen_pz->push_back(p.daughter(1)->pz());
					qrk_gen_pt->push_back(p.daughter(1)->pt());
					qrk_gen_E->push_back(p.daughter(1)->energy());
					qrk_gen_phi->push_back(p.daughter(1)->phi());
					qrk_gen_eta->push_back(p.daughter(1)->eta());
				}
				if(fabs(p.daughter(0)->pdgId())==12 || fabs(p.daughter(0)->pdgId()==14)){
					neut_gen_px->push_back(p.daughter(0)->px());
					neut_gen_py->push_back(p.daughter(0)->py());
					neut_gen_pz->push_back(p.daughter(0)->pz());
					neut_gen_pt->push_back(p.daughter(0)->pt());
					neut_gen_E->push_back(p.daughter(0)->energy());
					neut_gen_phi->push_back(p.daughter(0)->phi());
					neut_gen_eta->push_back(p.daughter(0)->eta());
				}
				if(fabs(p.daughter(1)->pdgId())==12 || fabs(p.daughter(1)->pdgId()==14)){
					neut_gen_px->push_back(p.daughter(1)->px());
					neut_gen_py->push_back(p.daughter(1)->py());
					neut_gen_pz->push_back(p.daughter(1)->pz());
					neut_gen_pt->push_back(p.daughter(1)->pt());
					neut_gen_E->push_back(p.daughter(1)->energy());
					neut_gen_phi->push_back(p.daughter(1)->phi());
					neut_gen_eta->push_back(p.daughter(1)->eta());
				}
				if(fabs(p.daughter(0)->pdgId())==13){
					muon_gen_px->push_back(p.daughter(0)->px());
					muon_gen_py->push_back(p.daughter(0)->py());
					muon_gen_pz->push_back(p.daughter(0)->pz());
					muon_gen_pt->push_back(p.daughter(0)->pt());
					muon_gen_E->push_back(p.daughter(0)->energy());
					muon_gen_phi->push_back(p.daughter(0)->phi());
					muon_gen_eta->push_back(p.daughter(0)->eta());
				}
				if(fabs(p.daughter(0)->pdgId())==11){
					ele_gen_px->push_back(p.daughter(0)->px());
					ele_gen_py->push_back(p.daughter(0)->py());
					ele_gen_pz->push_back(p.daughter(0)->pz());
					ele_gen_pt->push_back(p.daughter(0)->pt());
					ele_gen_E->push_back(p.daughter(0)->energy());
					ele_gen_phi->push_back(p.daughter(0)->phi());
					ele_gen_eta->push_back(p.daughter(0)->eta());
				}
				if(fabs(p.daughter(1)->pdgId())==13){
					muon_gen_px->push_back(p.daughter(1)->px());
					muon_gen_py->push_back(p.daughter(1)->py());
					muon_gen_pz->push_back(p.daughter(1)->pz());
					muon_gen_pt->push_back(p.daughter(1)->pt());
					muon_gen_E->push_back(p.daughter(1)->energy());
					muon_gen_phi->push_back(p.daughter(1)->phi());
					muon_gen_eta->push_back(p.daughter(1)->eta());
				}
				if(fabs(p.daughter(1)->pdgId())==11){
					ele_gen_px->push_back(p.daughter(1)->px());
					ele_gen_py->push_back(p.daughter(1)->py());
					ele_gen_pz->push_back(p.daughter(1)->pz());
					ele_gen_pt->push_back(p.daughter(1)->pt());
					ele_gen_E->push_back(p.daughter(1)->energy());
					ele_gen_phi->push_back(p.daughter(1)->phi());
					ele_gen_eta->push_back(p.daughter(1)->eta());
				}
			}
		}
	}
*/

	int BTag_ak4 = 0;
	if(jetAK8_pt->size()==1){
		for (size_t k = 0; k<jetAK4_phi->size(); k++) {
			if(!(jetAK4_pt->at(k)>20 && abs(jetAK4_eta->at(k))<2.4 && jetAK4_isTight->at(k))) continue;
			if (pow(pow(deltaPhi(jetAK4_phi->at(k),(jetAK8_phi->at(0))),2)+pow(( jetAK4_eta->at(k)-(jetAK8_eta->at(0))),2),0.5)>0.8) {
				if(jetAK4_btag->at(k)>0.9535){BTag_ak4++;}
			}
		}
	}


	double id_temp = -50.;
	double pt_max_temp = -999.;
	int i_max = -999;
	//	std::vector<double> *genP_pass;
/*	genP_pass->clear();
	if(Whad_gen_pt->size()==1 && jetAK8_pt->size()==1 && BTag_ak4<1){
		if(jetAK8_pt->at(0)>200){
			if((sqrt(pow(deltaPhi(Whad_gen_phi->at(0),jetAK8_phi->at(0)),2) + pow((Whad_gen_eta->at(0)-jetAK8_eta->at(0)),2))>0.5)){
				std::cout << "=================================================" << std::endl;
				std::cout << "passou1" << std::endl;
				for (size_t i=0; i< genP->size(); i++){
					const GenParticle & p = (*genP)[i];
					if((sqrt(pow(deltaPhi(p.phi(),jetAK8_phi->at(0)),2) + pow((p.eta()-jetAK8_eta->at(0)),2))>0.5) && abs(p.pdgId())!=6 && abs(p.pdgId())!=24) {
						if(p.pt()>pt_max_temp) {
							pt_max_temp = p.pt();
							id_temp = p.pdgId();
							i_max = i;
						}
					}
				}
				if(i_max>-1){
					std::cout << "MAIOR PT: " << id_temp << std::endl;
					const GenParticle & p = (*genP)[i_max];				
					const Candidate * mom = p.mother();
					std::cout << "-> p.pdgId(): " << p.pdgId() << "; mom->pdgId(): " << mom->pdgId() << "; (mom->mother())->pdgId(): " << (mom->mother())->pdgId();
					if(abs((mom->mother())->pdgId())!=2212){
						std::cout << "; (mom->mother()->mother())->pdgId(): " << (mom->mother()->mother())->pdgId()  <<  std::endl;
					} else {
						std::cout << std::endl;
					}
					std::cout << "-> p.status(): " << p.status() << "; mom->status(): " << mom->status() << "; (mom->mother())->status(): " << (mom->mother())->status();
                                        if(abs((mom->mother())->pdgId())!=2212){  
                                                std::cout << "; (mom->mother()->mother())->status(): " << (mom->mother()->mother())->status()  <<  std::endl;
                                        } else {
                                                std::cout << std::endl;
                                        }
					if(abs(id_temp)>10 && abs(id_temp)<19) count_gen_lep++;
                                        if(abs(id_temp)==21 || abs(id_temp)==9) count_gen_g++;
                                        if(abs(id_temp)==5) count_gen_q_b++;
                                        if(abs(id_temp)>0 && abs(id_temp)<5){
						if(abs((mom->mother())->pdgId())==2212 || abs((mom->mother())->pdgId())==21 || (abs((mom->mother())->pdgId())<5 && abs((mom->mother())->pdgId())>0)){
							if(abs((mom->mother())->pdgId())==2212) count_gen_lq_g++;
							if(abs((mom->mother())->pdgId())==21){
								if(abs((mom->mother()->mother())->pdgId())==2212) count_gen_lq_g++;
							}
							if(abs((mom->mother())->pdgId())<5 && abs((mom->mother())->pdgId())>0){
								if(abs((mom->mother()->mother())->pdgId())==2212) count_gen_lq_g++;
							}
						}
						if(abs(mom->pdgId())==24 || abs((mom->mother())->pdgId())==24) count_gen_lq_w++;
 					}
					count_gen++;
				} else count_gen_others++;
        			std::cout << "-> count_gen: " << count_gen << std::endl;
                                std::cout << "-> count_gen_lep: " << count_gen_lep << std::endl;
                                std::cout << "-> count_gen_q_b: " << count_gen_q_b << std::endl;
                                std::cout << "-> count_gen_lq_g: " << count_gen_lq_g << std::endl;
                                std::cout << "-> count_gen_lq_w: " << count_gen_lq_w << std::endl;
                                std::cout << "-> count_gen_g: " << count_gen_g << std::endl;
                                std::cout << "-> count_gen_others: " << count_gen_others << std::endl;
			}
		}
	}
*/

	EventBranchs->Fill();


#ifdef THIS_IS_AN_EVENT_EXAMPLE
	Handle<ExampleData> pIn;
	iEvent.getByLabel("example",pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
	ESHandle<SetupData> pSetup;
	iSetup.get<SetupRecord>().get(pSetup);
#endif
}


// ------------ method called once each job just before starting event loop  ------------
	void 
MakeNTuple::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
	void 
MakeNTuple::endJob() 
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
MakeNTuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
	//The following says we do not know what parameters are allowed so do no validation
	// Please change this to state exactly what you do use, even if it is no parameters
	edm::ParameterSetDescription desc;
	desc.setUnknown();
	descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(MakeNTuple);
