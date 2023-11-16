
////////////////////////////////////////////////////////////////////////
// Class:       AxionAnalyzer
// Plugin Type: analyzer (Unknown Unknown)
// File:        Atmospheric_module.cc
//
// Generated at Mon Jan 17 19:47:24 2022 by Leonardo Peres using cetskelgen
// from  version .
////////////////////////////////////////////////////////////////////////

#include "AxionAnalyzer.h"

// For a linearly correction of the charge deposited by an electron shower
double correctionGradient = 0.985;
double correctionIntercept = -0.02;
double fRecombFactor = 0.63;

constexpr int kDefInt = -9999;
constexpr int kDefMaxNRecoTracks = 1000;
constexpr int kDefDoub = (double)(kDefInt);

// Detector Limits =========================
float fFidVolXmin = 0;
float fFidVolXmax = 0;
float fFidVolYmin = 0;
float fFidVolYmax = 0;
float fFidVolZmin = 0;
float fFidVolZmax = 0;

axion::AxionAnalyzer::AxionAnalyzer(fhicl::ParameterSet const &p)
    : EDAnalyzer{p}, // More initializers here.
      fGeneratorModuleLabel(p.get<std::string>("GeneratorModuleLabel")),
      fGeantModuleLabel(p.get<std::string>("GeantModuleLabel")),
      fShowerModuleLabel(p.get<std::string>("ShowerModuleLabel")),
      fTrackModuleLabel(p.get<std::string>("TrackModuleLabel")),
      fPFParticleModuleLabel(p.get<std::string>("PFParticleModuleLabel")),
      fCaloModuleLabel(p.get<std::string>("CaloModuleLabel")),
      fPIDModuleLabel(p.get<std::string>("PIDModuleLabel")),
      fSpacePointModuleLabel(p.get<std::string>("SpacePointModuleLabel")),
      fHitModuleLabel(p.get<std::string>("HitModuleLabel")),
      fCVNModuleLabel(p.get<std::string>("CVNModuleLabel")),
      fMVAPIDModuleLabel(p.get<std::string>("MVAPIDModuleLabel")),
      fEnergyRecNCLabel(p.get<std::string>("EnergyRecNCLabel")),
      fSaveGeantInfo(p.get<bool>("SaveGeantInfo")),
      fCheatVertex(p.get<bool>("CheatVertex")),
      fShowerRecoSave(p.get<bool>("ShowerRecoSave")),
      fWireModuleLabel(p.get<std::string>("WireModuleLabel")),
      fNeutrinoEnergyRecoAlg(p.get<fhicl::ParameterSet>("NeutrinoEnergyRecoAlg"), fTrackModuleLabel, fShowerModuleLabel,
                             fHitModuleLabel, fWireModuleLabel, fTrackModuleLabel, fShowerModuleLabel, fPFParticleModuleLabel),
      fCalorimetryAlg(p.get<fhicl::ParameterSet>("CalorimetryAlg"))
{
    // Call appropriate consumes<>() for any products to be retrieved by this module.
}

void axion::AxionAnalyzer::ResetCounters()
{
    // std::cout << "Reseting counters... " << std::endl;

    fCosThetaDetTotalMom = -2;
    fCosPhiDetTotalMom = -2;
    fnTracks = 0;
    fnShowers = 0;
    fnSpacePoints = 0;
    fTotalMomentumP = 0;
    fPIDALongestTrack = 0;
    fHighestTrackSummedADC = 0;
    fLongestTrack = 0;
    fHighestShowerSummedADC = 0;
    fLargeShowerOpenAngle = -1;
    fLongestShower = -1;
    fCVN_NCScore = -1;
    fFracTotalChargeLongTrack = -1;
    fAvarageTrackLength = 0;
    fEventRecoEnergy = -1;
    //fEventRecoEnergy_numu = -1;
    fEventRecoEnergy_nue = -1;

    /*
      fRecoMethodUsed_NCEnergy = -2;
      fRecoVertex_NCEnergy.clear();
      fNuLorentzVector_NCEnergy.clear();
      fLepLorentzVector_NCEnergy.clear();
      fHadLorentzVector_NCEnergy.clear();
      flongestTrackContained_NCEnergy = -2;
      ftrackMomMethod_NCEnergy = -2;
    */

    fIsNC_CVNPred = false;

    fnGeantParticles_Primaries = 0;
    fnGenParticles = 0;

    fCVN_0protonsProbability = -1;
    fCVN_1protonsProbability = -1;
    fCVN_2protonsProbability = -1;
    fCVN_NprotonsProbability = -1;
    fCVN_0pionsProbability = -1;
    fCVN_0pizerosProbability = -1;
    fCVN_0neutronsProbability = -1;
    fCVN_IsAntineutrinoProbability = -1;
    fCVN_NumuProbability = -1;
    fCVN_NueProbability = -1;
    fCVN_NutauProbability = -1;
    fCVN_FlavourPred = -1;
    fCVN_NProtonsPred = -1;
    fCVN_AntineutrinoPred = -1;
    fCVN_NPionsPred = -1;
    fCVN_NPions0Pred = -1;

    fDistVertex = -2;
    fDiffCosAngleTotalMom = -2;
    fDiffCosAngleLongestTrack = -2;

    fNHits = 0;

    fDiffCosAngleTotalMom_AllProtons = -2;
    fTotalMomentumP_AllProtons = 0;
    fCosThetaDetTotalMom_AllProtons = -2;
    fCosPhiDetTotalMom_AllProtons = -2;
    fTotalMomRecoRangeUnitVect_AllProtons.clear();

    fDiffCosAngleTotalMom_AllMuons = -2;
    fTotalMomentumP_AllMuons = 0;
    fCosThetaDetTotalMom_AllMuons = -2;
    fCosPhiDetTotalMom_AllMuons = -2;
    fTotalMomRecoRangeUnitVect_AllMuons.clear();

    fDiffCosAngleTotalMom_MCS = -2;
    fTotalMomentumP_MCS = 0;
    fCosThetaDetTotalMom_MCS = -2;
    fCosPhiDetTotalMom_MCS = -2;
    fTotalMomMCSUnitVect.clear();

    fCCNC.clear();
    fNPrimaryDaughters.clear();
    fPrimaryPDGReco.clear();
    fNPrimaries = 0;
    InvertTrack = false;
    fMCGammaE = 0;
    fTotalMomentumTrueMag = 0;

    fDaughterTrackID.clear();
    fTrackStartX.clear();
    fTrackStartY.clear();
    fTrackStartZ.clear();
    fTrackEndX.clear();
    fTrackEndY.clear();
    fTrackEndZ.clear();
    fIsTrackFlipped.clear();
    fPIDAwithFlipMaybe.clear();
    fPIDANoFlip.clear();
    fPrimaryRecoVertex.clear();
    /*
    //MVA bits

      fRecoTrackMVAEvalRatio = 0;
      fRecoTrackMVAConcentration = 0;
      fRecoTrackMVACoreHaloRatio = 0;
      fRecoTrackMVAConicalness = 0;
      fRecoTrackMVAdEdxStart = 0;
      fRecoTrackMVAdEdxEnd = 0;
      fRecoTrackMVAdEdxEndRatio = 0;
      fRecoTrackMVAElectron = 0;
      fRecoTrackMVAPion = 0;
      fRecoTrackMVAMuon = 0;
      fRecoTrackMVAProton = 0;
      fRecoTrackMVAPhoton = 0;

      fRecoTrackTruePDG.clear();
      fRecoTrackTruePrimary.clear();
      fRecoTrackTrueMomX.clear();
      fRecoTrackTrueMomY.clear();
      fRecoTrackTrueMomZ.clear();
      fRecoTrackTrueMomT.clear();
      fRecoTrackTrueStartX.clear();
      fRecoTrackTrueStartY.clear();
      fRecoTrackTrueStartZ.clear();
      fRecoTrackTrueStartT.clear();
      fRecoTrackTrueEndX.clear();
      fRecoTrackTrueEndY.clear();
      fRecoTrackTrueEndZ.clear();
      fRecoTrackTrueEndT.clear();
    */

    fminChi2PDG.clear();
    fminChi2value.clear();
    fTotalMomRecoRangeUnitVect.clear();
    fTotalMomRecoCalVectUnit.clear();

    fTrksMom_AllProtons.clear();
    fTrksMom_AllMuons.clear();
    fTrksMom_BestFit.clear();
    fTrksMom_MCS.clear();

    fMCPartGenPDG.clear();
    fMCGammaGenMomentum.clear();
    fMCGammaGenEndMomentum.clear();
    fMCPartGenStatusCode.clear();

    fnGeantParticles = 0;
   // fThetaNuLepton.clear();
   // fMCNuMomentum.clear();
    fMCPrimaryGammaPDG.clear();
    fMCInitialPositionGamma.clear();
    fMCCosAzimuthGamma = -2;
    fMCTrackId.clear();
    fMCPdgCode.clear();
    fMCProcess.clear();
    fMCMomentum.clear();
    fMCStartEnergy.clear();
    fMCStatusCode.clear();
    fisMCinside.clear();
    fTotalMomentumUnitVect.clear();
    fMCKineticEnergy.clear();
    fTopologyNProton = -1;

    fShowerID.clear();
    fEnergyShowerLinearlyCorrected.clear();
    fShowerDirectionX.clear();
    fShowerDirectionY.clear();
    fShowerDirectionZ.clear();
    fShowerDirectionErr.clear();
    fShowerShowerStart.clear();
    fShowerShowerStartErr.clear();
    fShowerbest_plane.clear();
    fShowerLength.clear();
    fShowerOpenAngle.clear();
    fShowerEnergy.clear();
}

void axion::AxionAnalyzer::analyze(art::Event const &evt)
{
    // pid::Chi2PIDAlg fChiAlg;
    ResetCounters();

    TVector3 TotalMomentumRecoRange;
    TVector3 TotalMomentumRecoRange_AllProtons;
    TVector3 TotalMomentumRecoRange_AllMuons;
    TVector3 TotalMomentumMCS;
    TVector3 TotalMomentumRecoCal;
    TVector3 TrueEventDirection;
    TLorentzVector TotalMomentumTrue;
    // Get Geometry
    art::ServiceHandle<geo::Geometry const> geom;
    GeoLimits(geom, 10, 10, 10); // Exclude 10cm in each border

    // and the clock data for event
    auto const clockdata = art::ServiceHandle<detinfo::DetectorClocksService const>()->DataFor(evt);

    // Detector properties
    auto const detProp = art::ServiceHandle<detinfo::DetectorPropertiesService const>()->DataFor(evt, clockdata);

    // channel quality
    // lariov::ChannelStatusProvider const& channelStatus = art::ServiceHandle<lariov::ChannelStatusService const>()->GetProvider();

    size_t nplanes = geom->Nplanes();
    std::vector<std::vector<unsigned int>> hits(nplanes);
    m_run = evt.run();
    m_event = evt.id().event();

    //std::cout << "I'm here 1 " << std::endl;
    std::cout << "  Run: " << m_run << std::endl;
    std::cout << "  Event: " << m_event << std::endl;
    art::Ptr<recob::Shower> Shower_Longest;
    art::Ptr<recob::Track> Longest_trk;

    // Get CVN results
    art::Handle<std::vector<cvn::Result>> cvnResult;
    evt.getByLabel(fCVNModuleLabel, "cvnresult", cvnResult); // not sure if i need an instance name??

    if (!cvnResult->empty())
    {
        fCVN_NCScore = (*cvnResult)[0].GetNCProbability();
        fCVN_0protonsProbability = (*cvnResult)[0].Get0protonsProbability();
        fCVN_1protonsProbability = (*cvnResult)[0].Get1protonsProbability();
        fCVN_2protonsProbability = (*cvnResult)[0].Get2protonsProbability();
        fCVN_NprotonsProbability = (*cvnResult)[0].GetNprotonsProbability();
        fCVN_0pionsProbability = (*cvnResult)[0].Get0pionsProbability();
        fCVN_0pizerosProbability = (*cvnResult)[0].Get0pizerosProbability();
        fCVN_0neutronsProbability = (*cvnResult)[0].Get0neutronsProbability();
        fCVN_IsAntineutrinoProbability = (*cvnResult)[0].GetIsAntineutrinoProbability();
        fCVN_NumuProbability = (*cvnResult)[0].GetNumuProbability();
        fCVN_NueProbability = (*cvnResult)[0].GetNueProbability();
        fCVN_NutauProbability = (*cvnResult)[0].GetNutauProbability();
        fCVN_FlavourPred = (*cvnResult)[0].PredictedFlavour();
        fCVN_NProtonsPred = (*cvnResult)[0].PredictedProtons();
        fCVN_AntineutrinoPred = (*cvnResult)[0].PredictedIsAntineutrino();
        fCVN_NPionsPred = (*cvnResult)[0].PredictedPions();
        fCVN_NPions0Pred = (*cvnResult)[0].PredictedPizeros();
        //std::cout << "I'm here 2" << std::endl;
        if ((*cvnResult)[0].PredictedFlavour() == cvn::kFlavNC)
        {
            fIsNC_CVNPred = true;
        }
    }

    std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(evt)));
    fEventRecoEnergy = energyRecoHandle->fNuLorentzVector.E();
    /*
      art::Handle<dune::EnergyRecoOutput> EnergyNCOutput;
      evt.getByLabel(fEnergyRecNCLabel, "energyrecnc", EnergyNCOutput);

      fRecoMethodUsed_NCEnergy = 	(*EnergyNCOutput).recoMethodUsed;
      fRecoVertex_NCEnergy = {(*EnergyNCOutput).fRecoVertex.X(), (*EnergyNCOutput).fRecoVertex.Y(), (*EnergyNCOutput).fRecoVertex.Z()};
      fNuLorentzVector_NCEnergy ={(*EnergyNCOutput).fNuLorentzVector.X(),(*EnergyNCOutput).fNuLorentzVector.Y(), (*EnergyNCOutput).fNuLorentzVector.Z(), (*EnergyNCOutput).fNuLorentzVector.E()};
      fLepLorentzVector_NCEnergy = {(*EnergyNCOutput).fLepLorentzVector.X(),(*EnergyNCOutput).fLepLorentzVector.Y(), (*EnergyNCOutput).fLepLorentzVector.Z(), (*EnergyNCOutput).fLepLorentzVector.E()};
      fHadLorentzVector_NCEnergy = {(*EnergyNCOutput).fHadLorentzVector.X(),(*EnergyNCOutput).fHadLorentzVector.Y(), (*EnergyNCOutput).fHadLorentzVector.Z(), (*EnergyNCOutput).fHadLorentzVector.E()};
      flongestTrackContained_NCEnergy = (*EnergyNCOutput).longestTrackContained;
      ftrackMomMethod_NCEnergy = (*EnergyNCOutput).trackMomMethod;
    */
    std::map<int, float> PDGtoMass;
    PDGtoMass.insert(std::pair<int, float>(2212, 0.938272));
    PDGtoMass.insert(std::pair<int, float>(211, 0.13957));
    PDGtoMass.insert(std::pair<int, float>(321, 0.493677));
    PDGtoMass.insert(std::pair<int, float>(13, 0.105658));

    TVector3 vertical(0, 1, 0);

    auto const &MCTruthHandle = evt.getValidHandle<std::vector<simb::MCTruth>>(fGeneratorModuleLabel);
    auto const &MCTruthObjs = *MCTruthHandle;

    if (MCTruthObjs.size() > 0)
    {
        //std::cout << "I'm here 3 " << std::endl;
        for (size_t i = 0; i < MCTruthObjs.size(); i++)
        {
            simb::MCTruth MCTruthObj = MCTruthObjs[i];
            fnGenParticles = MCTruthObj.NParticles();
            //std::cout << "fnGenParticles = " << fnGenParticles << std::endl;

            const simb::MCNeutrino &MCNeutrino = MCTruthObj.GetNeutrino();
            fCCNC.push_back(MCNeutrino.CCNC());

            for (int j = 0; j < fnGenParticles; j++)
            {
                simb::MCParticle mc_part = MCTruthObj.GetParticle(j);
                fMCPartGenPDG.push_back(mc_part.PdgCode());
                fMCPartGenStatusCode.push_back(mc_part.StatusCode());

                //int part_cond = IsGammaAxion(mc_part.PdgCode(), mc_part.StatusCode());
                std::vector<double> tmp_mcMomentum = {mc_part.Momentum().X(), mc_part.Momentum().Y(), mc_part.Momentum().Z(), mc_part.Momentum().E()};
                fMCGammaGenMomentum.push_back(tmp_mcMomentum);
                fMCPrimaryGammaPDG.push_back(mc_part.PdgCode());
                TrueEventDirection = mc_part.Momentum().Vect().Unit();
                std::vector<double> tmp_fMCInitialPositionGamma = {mc_part.Vx(), mc_part.Vy(), mc_part.Vz()};
                fMCInitialPositionGamma.push_back(tmp_fMCInitialPositionGamma);
                fMCCosAzimuthGamma = vertical * mc_part.Momentum().Vect().Unit();
                //std::cout << "tmp_fMCInitialPositionGamma[0] = " << tmp_fMCInitialPositionGamma[0] << std::endl; 
                fMCGammaE = mc_part.Momentum().E();
            }
        }
    }

    // std::cout << "Number of CC and/or NC interactions: " << fCCNC.size() << std::endl;
    // std::cout << "Nu Interaction (0=CC 1=NC): " << fCCNC[0] << std::endl;

    // Truth information (Save geant4 info)
    if (fSaveGeantInfo)
    {

        auto const &MCParticleHandle = evt.getValidHandle<std::vector<simb::MCParticle>>(fGeantModuleLabel);
        auto const &MCParticleObjs = *MCParticleHandle;

        fnGeantParticles = MCParticleObjs.size();
        fnGeantParticles_Primaries = 0;
        //std::cout << "I'm here 4 " << std::endl;
        if (fnGeantParticles > 0)
        {
            for (size_t i = 0; i < MCParticleObjs.size(); i++)
            {
                //std::cout << "I'm here 5 " << std::endl;
                const simb::MCParticle MCParticleObj = MCParticleObjs[i];
                if (MCParticleObj.Process() != "primary")
                {
                    continue;
                }
                fnGeantParticles_Primaries++;
                fMCTrackId.push_back(MCParticleObj.TrackId());
                fMCPdgCode.push_back(MCParticleObj.PdgCode());
                fMCProcess.push_back(MCParticleObj.Process());
                fMCMomentum.push_back({MCParticleObj.Px(), MCParticleObj.Py(), MCParticleObj.Pz()});
                fMCStartEnergy.push_back(MCParticleObj.E(0));
                fMCKineticEnergy.push_back(MCParticleObj.E(0) - MCParticleObj.Mass());
                fMCStatusCode.push_back(MCParticleObj.StatusCode());
                geo::Point_t EndMCPos(MCParticleObj.EndPosition().X(), MCParticleObj.EndPosition().Y(), MCParticleObj.EndPosition().Z());
                bool isMCinside_tmp = insideFV(EndMCPos);
                fisMCinside.push_back(isMCinside_tmp);

                // True total momentum -> it has to be a visible particle,
                //               it has to be stopped inside the detector,
                //               it has to be a primary particle,
                //               it has to be a stable final state particle.
                //               so, we do not considered neutrons in the true total momentum.

                if (IsVisibleParticle(MCParticleObj.PdgCode(), MCParticleObj.Process()) && isMCinside_tmp && MCParticleObj.StatusCode() == 1)
                {
                    TotalMomentumTrue += MCParticleObj.Momentum();
                }
            }
        }

        TVector3 TotalMomTrue = TotalMomentumTrue.Vect().Unit();
        std::vector<double> TotalMomTrueXYZ = {TotalMomTrue.X(), TotalMomTrue.Y(), TotalMomTrue.Z()};
        fTotalMomentumUnitVect = TotalMomTrueXYZ;
        fTotalMomentumTrueMag = TotalMomentumTrue.Mag();
        fTopologyNProton = Topology(fMCPdgCode);
    }
    //std::cout << "I'm here 6 " << std::endl;
    m_AllEvents->Fill();

    // Collect the PFParticles from the event
    PFParticleHandle pfParticleHandle;
    evt.getByLabel(fPFParticleModuleLabel, pfParticleHandle);

    TrackHandle trackHandle;
    evt.getByLabel(fTrackModuleLabel, trackHandle);

    ShowerHandle showerHandle;
    evt.getByLabel(fShowerModuleLabel, showerHandle);

    SpacePointHandle spacepointHandle;
    std::vector<art::Ptr<recob::SpacePoint>> SpacePointVect;
    if (evt.getByLabel(fPFParticleModuleLabel, spacepointHandle))
        art::fill_ptr_vector(SpacePointVect, spacepointHandle);

    HitHandle hitHandle;
    std::vector<art::Ptr<recob::Hit>> HitVect;
    if (evt.getByLabel(fHitModuleLabel, hitHandle))
        art::fill_ptr_vector(HitVect, hitHandle);
    //std::cout << "I'm here 7 " << std::endl;

    if (!pfParticleHandle.isValid())
    {
        mf::LogDebug("AxionAnalyzer") << "  Failed to find the PFParticles." << std::endl;
        return;
    }

    // Get the associations between PFParticles and tracks/showers from the event and track from calo
    art::FindManyP<recob::Track> pfPartToTrackAssoc(pfParticleHandle, evt, fTrackModuleLabel);
    art::FindManyP<recob::Shower> pfPartToShowerAssoc(pfParticleHandle, evt, fShowerModuleLabel);
    art::FindManyP<anab::Calorimetry> TrackToCaloAssoc(trackHandle, evt, fCaloModuleLabel);
    art::FindManyP<recob::Hit> TrackToHitsAssoc(trackHandle, evt, fTrackModuleLabel);
    art::FindManyP<anab::ParticleID> TrackToPIDAssoc(trackHandle, evt, fPIDModuleLabel);
    art::FindManyP<recob::SpacePoint> pfpToSpacePoint(pfParticleHandle, evt, fSpacePointModuleLabel);
    art::FindManyP<recob::Vertex> pfPartToVertex(pfParticleHandle, evt, fPFParticleModuleLabel);
    art::FindManyP<recob::Hit> ShowerToHitAssoc(showerHandle, evt, fShowerModuleLabel);
    art::FindManyP<recob::SpacePoint> ShowerToSpacePoint(showerHandle, evt, fShowerModuleLabel);
    art::FindManyP<dune::EnergyRecoOutput> TrackToEnergyNCAssoc(trackHandle, evt, fEnergyRecNCLabel);
    // std::cout << "I'm in the association part, dad! " << std::endl;
    art::FindManyP<anab::MVAPIDResult> fmpidt(trackHandle, evt, fMVAPIDModuleLabel);
    art::FindManyP<anab::MVAPIDResult> fmpids(showerHandle, evt, fMVAPIDModuleLabel);
    // art::FindManyP<recob::Hit> TrackToHitsAssoc(trackHandle, evt, fHitModuleLabel);

    const std::vector<art::Ptr<recob::PFParticle>> pfparticleVect = dune_ana::DUNEAnaEventUtils::GetPFParticles(evt, fPFParticleModuleLabel);
    //std::cout << "I'm here 8 " << std::endl;
    size_t pfpID = 99999;
    bool isInsideFD = false; // flag for fidicial cut in the reconstructed vertex
                             // Double_t VertexXYZ[3] = {};

    // Block adapted from ConsolidatedPFParticleAnalysisTemplate should work fine!
    for (size_t iPfp = 0; iPfp < pfparticleVect.size(); iPfp++)
    {

        // const std::vector<art::Ptr<recob::SpacePoint>> &associatedSpacePoints = pfpToSpacePoint.at(iPfp);
        const std::vector<art::Ptr<recob::Vertex>> &associatedVertex = pfPartToVertex.at(iPfp);

        if (!(pfparticleVect[iPfp]->IsPrimary() /*&& (pfparticleVect[iPfp]->PdgCode() == 14 || pfparticleVect[iPfp]->PdgCode() == 12)*/ ))
            continue;

        // const int pdg(pParticle->PdgCode());
        // std::cout << "We have Primary Particles! Yep" << std::endl;
        // std::cout << "associatedVertex.size() = " << associatedVertex.size() << std::endl;
        //std::cout << "I'm here 9 " << std::endl;
        for (const art::Ptr<recob::Vertex> &vtx : associatedVertex)
        {
            if (!(vtx->isValid()))
                continue; // Valid Vertice
            //std::cout << "I'm here 10 " << std::endl;
            vtxPos vertex = vtx->position();
            if (!insideFV(vertex))
                continue;                                                       // Fiducial Reco Cut 20 cm away from the border
            fPrimaryRecoVertex.push_back({vertex.X(), vertex.Y(), vertex.Z()}); // Save vertex position
            isInsideFD = true;
            fDistVertex = pow(pow((fMCInitialPositionGamma.at(0).at(0) - vertex.X()), 2) + pow((fMCInitialPositionGamma.at(0).at(1) - vertex.Y()), 2) + pow((fMCInitialPositionGamma.at(0).at(2) - vertex.Z()), 2), 0.5);
        }

        pfpID = pfparticleVect[iPfp]->Self();
        fNPrimaryDaughters.push_back(pfparticleVect[iPfp]->NumDaughters());
        fPrimaryPDGReco.push_back(pfparticleVect[iPfp]->PdgCode());

        fNPrimaries++;
    }
    std::cout << "I'm here 11 " << std::endl;
    if (pfpID == 99999 || !isInsideFD)
        return; // neutrino candidate and reconstructed inside fiducial cut

    fnSpacePoints += SpacePointVect.size();

    double TotalHitsADC = 0;
    std::cout << "I'm here 12 " << std::endl;
    for (size_t iHit = 0; iHit < HitVect.size(); iHit++)
    {

        TotalHitsADC += HitVect[iHit]->SummedADC();
    }

    fNHits = HitVect.size();

    double LongestTrack = 1e-10;
    double HighestTrackSummedADC = 1e-10;
    long unsigned int LongestTrackID = -2;
    double PIDALongestTrack = 0;
    double AllLengthTracksSummed = 0;
    std::cout << "I'm here 13 " << std::endl;

    for (size_t iPfp = 0; iPfp < pfparticleVect.size(); iPfp++)
    {

        //if (pfparticleVect[iPfp]->Parent() != neutrinoID) continue;
        const std::vector<art::Ptr<recob::Track>> &associatedTracks = pfPartToTrackAssoc.at(iPfp);

        float minChi2value;
        int minChi2PDG;

        double trk_P_proton = 0;
        double trk_P_muon = 0;
        double trk_MSC = 0;

        if (!associatedTracks.empty())
        {

            // std::cout << "We have Tracks! Yep" << std::endl;
            for (const art::Ptr<recob::Track> &trk : associatedTracks)
            {
                std::vector<art::Ptr<recob::Hit>> trackHits = TrackToHitsAssoc.at(trk.key());

                if (trk->NumberTrajectoryPoints() < 3)
                    continue;

                if (!insideFV(trk->End()))
                    continue;
                if (!insideFV(trk->Start()))
                    continue;

                trk_P_proton = trkm.GetTrackMomentum(trk->Length(), 2212);
                trk_P_muon = trkm.GetTrackMomentum(trk->Length(), 13);
                trk_MSC = trkm.GetMomentumMultiScatterLLHD(trk);

                // int trkidtruth = TruthMatchUtils::TrueParticleIDFromTotalTrueEnergy(clockData, trackHits, true);
                // const simb::MCParticle *particle = pi_serv->TrackIdToParticle_P(trkidtruth);
                /* int g4id = TruthMatchUtils::TrueParticleIDFromTotalRecoHits(clockdata, trackHits, 1);
                if (TruthMatchUtils::Valid(g4id)){
                  std::cout << "TruthMatchUtils::Valid" << std::endl;
                art::ServiceHandle<cheat::ParticleInventoryService> pi_serv;
                const simb::MCParticle* matched_mcparticle = pi_serv->ParticleList().at(g4id);
                  if (matched_mcparticle){
                    //Fill variables
                    std::cout << "Fill variables" << std::endl;
                    fRecoTrackTruePDG.push_back(matched_mcparticle->PdgCode());
                    if (matched_mcparticle->Mother()==0) fRecoTrackTruePrimary.push_back(true);
                    else fRecoTrackTruePrimary.push_back(false);
                    fRecoTrackTrueMomX.push_back(matched_mcparticle->Momentum().X());
                    fRecoTrackTrueMomY.push_back(matched_mcparticle->Momentum().Y());
                    fRecoTrackTrueMomZ.push_back(matched_mcparticle->Momentum().Z());
                    fRecoTrackTrueMomT.push_back(matched_mcparticle->Momentum().T());
                    fRecoTrackTrueStartX.push_back(matched_mcparticle->Position(0).X());
                    fRecoTrackTrueStartY.push_back(matched_mcparticle->Position(0).Y());
                    fRecoTrackTrueStartZ.push_back(matched_mcparticle->Position(0).Z());
                    fRecoTrackTrueStartT.push_back(matched_mcparticle->Position(0).T());
                    fRecoTrackTrueEndX.push_back(matched_mcparticle->EndPosition().X());
                    fRecoTrackTrueEndY.push_back(matched_mcparticle->EndPosition().Y());
                    fRecoTrackTrueEndZ.push_back(matched_mcparticle->EndPosition().Z());
                    fRecoTrackTrueEndT.push_back(matched_mcparticle->EndPosition().T());
                  }
                }*/

                fDaughterTrackID.push_back(trk.key());

                float trackADC = 0;

                for (const art::Ptr<recob::Hit> &hit : trackHits)
                {

                    trackADC += hit->SummedADC();
                }
                if (trackADC > fHighestTrackSummedADC)
                    HighestTrackSummedADC = trackADC;

                TVector3 TrackDirectionLongestTrack(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());

                AllLengthTracksSummed += trk->Length();

                /*
                const std::vector<art::Ptr<dune::EnergyRecoOutput>> &energyNCfromTracks = TrackToEnergyNCAssoc.at(trk.key());

                for(const art::Ptr<dune::EnergyRecoOutput> &energy_nc : energyNCfromTracks){
                  fRecoMethodUsed_NCEnergy = 	(*energy_nc).recoMethodUsed;
                fRecoVertex_NCEnergy = {(*energy_nc).fRecoVertex.X(), (*energy_nc).fRecoVertex.Y(), (*energy_nc).fRecoVertex.Z()};
                fNuLorentzVector_NCEnergy ={(*energy_nc).fNuLorentzVector.X(),(*energy_nc).fNuLorentzVector.Y(), (*energy_nc).fNuLorentzVector.Z(), (*energy_nc).fNuLorentzVector.E()};
                fLepLorentzVector_NCEnergy = {(*energy_nc).fLepLorentzVector.X(),(*energy_nc).fLepLorentzVector.Y(), (*energy_nc).fLepLorentzVector.Z(), (*energy_nc).fLepLorentzVector.E()};
                fHadLorentzVector_NCEnergy = {(*energy_nc).fHadLorentzVector.X(),(*energy_nc).fHadLorentzVector.Y(), (*energy_nc).fHadLorentzVector.Z(), (*energy_nc).fHadLorentzVector.E()};
                flongestTrackContained_NCEnergy = (*energy_nc).longestTrackContained;
                ftrackMomMethod_NCEnergy = (*energy_nc).trackMomMethod;
                }
                */

                // Variables for the longest track
                if (trk->Length() > LongestTrack)
                {
                    LongestTrack = trk->Length();
                    LongestTrackID = trk.key();
                    fDiffCosAngleLongestTrack = TrackDirectionLongestTrack * TrueEventDirection;
                    Longest_trk = trk;

                    // std::cout << "I'm accessing the MMVAPID, dad! " << std::endl;
                    /*std::vector<art::Ptr<anab::MVAPIDResult> > pids = fmpidt.at(trk.key());

                    if (pids.at(0).isAvailable()){
                      fRecoTrackMVAEvalRatio = pids.at(0)->evalRatio;
                      fRecoTrackMVAConcentration = pids.at(0)->concentration;
                      fRecoTrackMVACoreHaloRatio = pids.at(0)->coreHaloRatio;
                      fRecoTrackMVAConicalness = pids.at(0)->conicalness;
                      fRecoTrackMVAdEdxStart = pids.at(0)->dEdxStart;
                      fRecoTrackMVAdEdxEnd = pids.at(0)->dEdxEnd;
                      fRecoTrackMVAdEdxEndRatio = pids.at(0)->dEdxEndRatio;
                      std::map<std::string,double> mvaOutMap = pids.at(0)->mvaOutput;
                    if (!(mvaOutMap.empty())){
                        //Get the PIDs
                        fRecoTrackMVAElectron = mvaOutMap["electron"];
                        fRecoTrackMVAPion = mvaOutMap["pich"];
                        fRecoTrackMVAMuon = mvaOutMap["muon"];
                        fRecoTrackMVAProton = mvaOutMap["proton"];
                        fRecoTrackMVAPhoton = mvaOutMap["photon"];
                    }
                    }*/
                }

                std::vector<art::Ptr<anab::ParticleID>> trackPID = TrackToPIDAssoc.at(trk.key());
                std::map<int, float> PDGtoChi2;

                for (size_t i = 0; i < trackPID.size(); i++)
                {

                    std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(i)->ParticleIDAlgScores();
                    minChi2value = 9999;
                    minChi2PDG = 9999;

                    // Loop through AlgScoresVec and find the variables we want
                    for (size_t i_algscore = 0; i_algscore < AlgScoresVec.size(); i_algscore++)
                    {
                        anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
                        // std::cout << "AlgScore.fAlgName = " << AlgScore.fAlgName << std::endl;
                        //  if (AlgScore.fPlaneMask[2] != 1)  continue; //Only collection plane
                        if (AlgScore.fAlgName == "Chi2")
                        {
                            // Sum Chi2 hyphotesis all planes
                            // if(AlgScore.fAssumedPdg == 13) std::cout <<  "AlgScore.fAssumedPdg = " << AlgScore.fAssumedPdg << std::endl;
                            // if(!(AlgScore.fAssumedPdg == 13 || AlgScore.fAssumedPdg == 211 || AlgScore.fAssumedPdg == 2212 || AlgScore.fAssumedPdg == 321)) std::cout << "AlgScore.fAssumedPdg = " << AlgScore.fAssumedPdg << std::endl;
                            PDGtoChi2[AlgScore.fAssumedPdg] = AlgScore.fValue;
                        }
                    }
                }

                // std::cout << "\nThe map PDGtoChi2 is : \n";
                // std::cout << "\tKEY\tELEMENT\n";
                std::map<int, float>::iterator it;
                for (it = PDGtoChi2.begin(); it != PDGtoChi2.end(); ++it)
                {
                    // std::cout << '\t' << it->first << '\t' << it->second << '\n';
                    if (it->second < minChi2value)
                    {
                        minChi2value = it->second;
                        minChi2PDG = it->first;
                    }
                }
                if (minChi2PDG == 9999)
                    continue;
                fnTracks++;

                // std::cout << "min Ch2 hypothesis : " << minChi2PDG << " with " << minChi2value << ". \n";
                fminChi2value.push_back(minChi2value);
                fminChi2PDG.push_back(minChi2PDG);
                // std::min(chi2_pion_backward, std::min(chi2_pion_forward, std::min(chi2_proton_backward, chi2_proton_forward)));

                const std::vector<art::Ptr<anab::Calorimetry>> associatedCalo = TrackToCaloAssoc.at(trk.key());
                if (associatedCalo.empty())
                    continue;
                float KE = 0;

                for (const art::Ptr<anab::Calorimetry> &cal : associatedCalo)
                {
                    if (!cal->PlaneID().isValid)
                        continue;
                    int planenum = cal->PlaneID().Plane;
                    // std::cout << "pid: " << pidout.ParticleIDAlgScores.at(0) << std::endl;
                    std::vector<float> temp_dEdx = cal->dEdx();

                    if (planenum == 2)
                    {
                        KE = cal->KineticEnergy();
                    }

                    temp_dEdx.clear();
                } // Calo

                TVector3 TrackDirection(trk->StartDirection().X(), trk->StartDirection().Y(), trk->StartDirection().Z());
                TVector3 DaughterStartPoint(trk->Start().X(), trk->Start().Y(), trk->Start().Z());
                TVector3 DaughterEndPoint(trk->End().X(), trk->End().Y(), trk->End().Z());

                fTrackStartX.push_back(DaughterStartPoint.X());
                fTrackStartY.push_back(DaughterStartPoint.Y());
                fTrackStartZ.push_back(DaughterStartPoint.Z());
                fTrackEndX.push_back(DaughterEndPoint.X());
                fTrackEndY.push_back(DaughterEndPoint.Y());
                fTrackEndZ.push_back(DaughterEndPoint.Z());

                // Cheat Vertex, all track directions based on the true vertex
                double PIDANoFlip = CalcPIDA(associatedCalo, InvertTrack);
                fPIDANoFlip.push_back(PIDANoFlip);

                // IsTrackBackwardsAndFlip(associatedCalo, DaughterStartPoint, DaughterEndPoint, InvertTrack);

                if (InvertTrack)
                    TrackDirection = -TrackDirection;

                double PIDAwithFlipMaybe = CalcPIDA(associatedCalo, InvertTrack);
                // std::cout << "InvertTrack = " << InvertTrack << "." << std::endl;
                //   std::cout << "PIDA = " << PIDA << std::endl;
                fPIDAwithFlipMaybe.push_back(PIDAwithFlipMaybe);
                fIsTrackFlipped.push_back(InvertTrack);

                if (trk.key() == LongestTrackID)
                    PIDALongestTrack = PIDAwithFlipMaybe;

                fTrksMom_AllProtons.push_back({trk_P_proton * TrackDirection.X(), trk_P_proton * TrackDirection.Y(), trk_P_proton * TrackDirection.Z()});
                fTrksMom_AllMuons.push_back({trk_P_muon * TrackDirection.X(), trk_P_muon * TrackDirection.Y(), trk_P_muon * TrackDirection.Z()});
                fTrksMom_MCS.push_back({trk_MSC * TrackDirection.X(), trk_MSC * TrackDirection.Y(), trk_MSC * TrackDirection.Z()});

                if (minChi2PDG == 13 || minChi2PDG == 211)
                {
                    TotalMomentumRecoRange += trk_P_muon * TrackDirection;
                    fTrksMom_BestFit.push_back({trk_P_muon * TrackDirection.X(), trk_P_muon * TrackDirection.Y(), trk_P_muon * TrackDirection.Z(), PDGtoMass[minChi2PDG]});
                }

                if (minChi2PDG == 2212 || minChi2PDG == 321)
                {
                    TotalMomentumRecoRange += trk_P_proton * TrackDirection;
                    fTrksMom_BestFit.push_back({trk_P_proton * TrackDirection.X(), trk_P_proton * TrackDirection.Y(), trk_P_proton * TrackDirection.Z(), PDGtoMass[minChi2PDG]});
                }
                TotalMomentumRecoRange_AllProtons += trk_P_proton * TrackDirection;
                TotalMomentumRecoRange_AllMuons += trk_P_muon * TrackDirection;
                TotalMomentumMCS += trk_MSC * TrackDirection;

                double Pcal = sqrt((KE + PDGtoMass[minChi2PDG]) * (KE + PDGtoMass[minChi2PDG]) - PDGtoMass[minChi2PDG] * PDGtoMass[minChi2PDG]);
                TotalMomentumRecoCal += Pcal * TrackDirection;

            } // Tracks- loop

            fAvarageTrackLength = AllLengthTracksSummed / fnTracks;

        } // if there is a track
    }
    std::cout << "I'm here 14 " << std::endl;

    for (size_t iPfp = 0; iPfp < pfparticleVect.size(); iPfp++)
    {
        if (fShowerRecoSave)
        {
            std::cout << "I'm here 15 " << std::endl;

            //if (pfparticleVect[iPfp]->Parent() != neutrinoID) continue;
            // SHOWERS RECO INFO ====================================================================
            const std::vector<art::Ptr<recob::Shower>> &associatedShowers = pfPartToShowerAssoc.at(iPfp);
            fnShowers += associatedShowers.size();
            // std::cout << "associatedShowers.size() = " << associatedShowers.size() << std::endl;

            if (!associatedShowers.empty())
            {
                for (const art::Ptr<recob::Shower> &Shower : associatedShowers)
                {
                  std::cout << "I'm here 16" << std::endl;
                    std::vector<art::Ptr<recob::SpacePoint>> showersp = ShowerToSpacePoint.at(Shower.key());
                    // std::cout << "showersp.size() = " << showersp.size() << std::endl;
                    if (showersp.size() == 0)
                        continue;

                    if (Shower->Direction().X() == -999)
                        continue;

                    const std::vector<art::Ptr<recob::Hit>> electronHits(dune_ana::DUNEAnaHitUtils::GetHitsOnPlane(dune_ana::DUNEAnaShowerUtils::GetHits(Shower, evt, fShowerModuleLabel), 2));
                    const double electronObservedCharge(dune_ana::DUNEAnaHitUtils::LifetimeCorrectedTotalHitCharge(clockdata, detProp, electronHits));
                    const double uncorrectedElectronEnergy = fCalorimetryAlg.ElectronsFromADCArea(electronObservedCharge, 2) * 1. / fRecombFactor / util::kGeVToElectrons;
                    double Showerenergy = (uncorrectedElectronEnergy - correctionIntercept) / correctionGradient;
                    fEnergyShowerLinearlyCorrected.push_back(Showerenergy);
                    fShowerID.push_back(Shower->ID());
                    fShowerDirectionX.push_back(Shower->Direction().X());
                    fShowerDirectionY.push_back(Shower->Direction().Y());
                    fShowerDirectionZ.push_back(Shower->Direction().Z());
                    fShowerDirectionErr.push_back({Shower->DirectionErr().X(), Shower->DirectionErr().Y(), Shower->DirectionErr().Z()});
                    fShowerShowerStart.push_back({Shower->ShowerStart().X(), Shower->ShowerStart().Y(), Shower->ShowerStart().Z()});
                    fShowerShowerStartErr.push_back({Shower->ShowerStartErr().X(), Shower->ShowerStartErr().Y(), Shower->ShowerStartErr().Z()});
                    fShowerbest_plane.push_back(Shower->best_plane());
                    fShowerEnergy.push_back(Shower->Energy());
                    fShowerLength.push_back(Shower->Length());
                    fShowerOpenAngle.push_back(Shower->OpenAngle());

                    TVector3 showerDirection(Shower->Direction().X(), Shower->Direction().Y(), Shower->Direction().Z());
                    std::vector<art::Ptr<recob::Hit>> ShowerHits = ShowerToHitAssoc.at(Shower.key());

                    TotalMomentumRecoRange += Showerenergy * showerDirection;
                    TotalMomentumRecoRange_AllProtons += Showerenergy * showerDirection;
                    TotalMomentumRecoRange_AllMuons += Showerenergy * showerDirection;
                    TotalMomentumMCS += Showerenergy * showerDirection;
                    if (Shower->Length() > fLongestShower)
                    {
                        Shower_Longest = Shower;
                        fLongestShower = Shower->Length();
                    }
                    if (Shower->OpenAngle() > fLargeShowerOpenAngle)
                        fLargeShowerOpenAngle = Shower->OpenAngle();

                    float showerADC = 0;
                    for (const art::Ptr<recob::Hit> &hit : ShowerHits)
                    {
                        showerADC += hit->SummedADC();
                    }

                    if (showerADC > fHighestShowerSummedADC)
                        fHighestShowerSummedADC = showerADC;
                }
            }
        }
    }

    // TVector3 z(0,0,1);
    // std::cout << "TotalMomentumRecoRange.Mag() = " << TotalMomentumRecoRange.Mag() << std::endl;
    if (TotalMomentumRecoRange.Mag() > 0.0)
    {
        std::cout << "I'm here 17 " << std::endl;
        std::cout << "Longest_trk.Length= " << Longest_trk->Length() << std::endl;
        //std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle_numu(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(Longest_trk, evt)));
        std::cout << "I'm here 18-1 " << std::endl;
        //fEventRecoEnergy_numu = energyRecoHandle_numu->fNuLorentzVector.E();
        
        std::cout << "I'm here 18-2 " << std::endl;

        std::unique_ptr<dune::EnergyRecoOutput> energyRecoHandle_nue(std::make_unique<dune::EnergyRecoOutput>(fNeutrinoEnergyRecoAlg.CalculateNeutrinoEnergy(Shower_Longest, evt)));
        fEventRecoEnergy_nue = energyRecoHandle_nue->fNuLorentzVector.E();

        fHighestTrackSummedADC = HighestTrackSummedADC;
        fLongestTrack = LongestTrack;
        fPIDALongestTrack = PIDALongestTrack;
        fFracTotalChargeLongTrack = HighestTrackSummedADC / TotalHitsADC;

        std::cout << "I'm here 19 " << std::endl;
        fTotalMomentumP = TotalMomentumRecoRange.Mag();
        fCosThetaDetTotalMom = TotalMomentumRecoRange.Unit().CosTheta();
        fCosPhiDetTotalMom = cos(TotalMomentumRecoRange.Unit().Phi());
        fTotalMomRecoRangeUnitVect = {TotalMomentumRecoRange.Unit().X(), TotalMomentumRecoRange.Unit().Y(), TotalMomentumRecoRange.Unit().Z()};

        fTotalMomentumP_AllProtons = TotalMomentumRecoRange_AllProtons.Mag();
        fCosThetaDetTotalMom_AllProtons = TotalMomentumRecoRange_AllProtons.Unit().CosTheta();
        fCosPhiDetTotalMom_AllProtons = cos(TotalMomentumRecoRange_AllProtons.Unit().Phi());
        fTotalMomRecoRangeUnitVect_AllProtons = {TotalMomentumRecoRange_AllProtons.Unit().X(), TotalMomentumRecoRange_AllProtons.Unit().Y(), TotalMomentumRecoRange_AllProtons.Unit().Z()};

        fTotalMomentumP_AllMuons = TotalMomentumRecoRange_AllMuons.Mag();
        fCosThetaDetTotalMom_AllMuons = TotalMomentumRecoRange_AllMuons.Unit().CosTheta();
        fCosPhiDetTotalMom_AllMuons = cos(TotalMomentumRecoRange_AllMuons.Unit().Phi());
        fTotalMomRecoRangeUnitVect_AllMuons = {TotalMomentumRecoRange_AllMuons.Unit().X(), TotalMomentumRecoRange_AllMuons.Unit().Y(), TotalMomentumRecoRange_AllMuons.Unit().Z()};

        fTotalMomentumP_MCS = TotalMomentumMCS.Mag();
        fCosThetaDetTotalMom_MCS = TotalMomentumMCS.Unit().CosTheta();
        fCosPhiDetTotalMom_MCS = cos(TotalMomentumMCS.Unit().Phi());
        fTotalMomMCSUnitVect = {TotalMomentumMCS.Unit().X(), TotalMomentumMCS.Unit().Y(), TotalMomentumMCS.Unit().Z()};

        // std::cout << "fTotalMomRecoRangeUnitVect =" << fTotalMomRecoRangeUnitVect <<std::endl;
        fTotalMomRecoCalVectUnit = {TotalMomentumRecoCal.Unit().X(), TotalMomentumRecoCal.Unit().Y(), TotalMomentumRecoCal.Unit().Z()};

        fDiffCosAngleTotalMom = TotalMomentumRecoRange.Unit() * TrueEventDirection;
        fDiffCosAngleTotalMom_AllProtons = TotalMomentumRecoRange_AllProtons.Unit() * TrueEventDirection;
        fDiffCosAngleTotalMom_AllMuons = TotalMomentumRecoRange_AllMuons.Unit() * TrueEventDirection;
        fDiffCosAngleTotalMom_MCS = TotalMomentumMCS.Unit() * TrueEventDirection;
        // std::cout << "fTotalMomRecoCalVectUnit =" << fTotalMomRecoCalVectUnit <<std::endl;
    }
    std::cout << "I'm here 20 " << std::endl;
    // Fill the Tree just for a NC event, and if there is at least one track per event in the BDT variables
    if (TotalMomentumRecoRange.Mag() > 0 && isInsideFD)
        m_AtmTree->Fill();
}

// Evaluate the charge deposition in the beggining and at the end of the track
void axion::AxionAnalyzer::IsTrackBackwardsAndFlip(std::vector<art::Ptr<anab::Calorimetry>> calos, TVector3 trkStart, TVector3 trkEnd, bool IsFlipped)
{

    // I need to decide which way to go through the hits...
    double SumdEdx_St = 0., SumdEdx_En = 0.; // ResRng_St = 0, ResRng_En = 0;
    unsigned int AvHits = 5;
    // if (calos[2]->dEdx().size()) {
    //   ResRng_St  = calos[2]->ResidualRange()[0];
    //   ResRng_En  = calos[2]->ResidualRange()[calos[2]->dEdx().size()-1];
    // }
    // If don't have 2*AvHits (10) hits on the collection plane then use mid point + 1 hit
    if (calos[2]->dEdx().size() < (2 * AvHits))
    {
        AvHits = 1 + (0.5 * calos[2]->dEdx().size());
    }
    for (unsigned int PlHit = 0; PlHit < calos[2]->dEdx().size(); ++PlHit)
    { // loop through hits on the collection plane
        if (PlHit <= AvHits)
        {
            SumdEdx_St += calos[2]->dEdx()[PlHit];
        }
        if (calos[2]->dEdx().size() - PlHit <= AvHits)
        {
            SumdEdx_En += calos[2]->dEdx()[PlHit];
        }
        // std::cout << "Looking at hit " << PlHit << " of " << (int)calos[2]->dEdx().size() << "...SumdEdx_St = " << SumdEdx_St << ", and SumdEdx_En = " << SumdEdx_En << std::endl;
    }
    double AvdEdx_St = SumdEdx_St / AvHits;
    double AvdEdx_En = SumdEdx_En / AvHits;
    // The dEdx at the start of the track should be less than that at the end...
    bool LowdEdxSt = false;
    std::cout << "AvdEdx_St = " << AvdEdx_St << "\t"
              << "AvdEdx_En = " << AvdEdx_En << "." << std::endl;
    if (AvdEdx_St < AvdEdx_En)
        LowdEdxSt = true;

    if (!LowdEdxSt)
    {
        std::cout << "Track backwards, energy deposition at the start higher than at the end!" << std::endl;
        double TmpX, TmpY, TmpZ;
        TmpX = trkStart.X();
        TmpY = trkStart.Y();
        TmpZ = trkStart.Z();
        trkStart.SetX(trkEnd.X());
        trkStart.SetY(trkEnd.Y());
        trkStart.SetZ(trkEnd.Z());
        trkEnd.SetX(TmpX);
        trkEnd.SetY(TmpY);
        trkEnd.SetZ(TmpZ);
    }

    InvertTrack = !LowdEdxSt;
    std::cout << "InvertTrack = " << InvertTrack << std::endl;
}

double axion::AxionAnalyzer::CalcPIDA(std::vector<art::Ptr<anab::Calorimetry>> calos, bool IsInverted)
{

    double PIDA = 0;
    int UsedHits = 0, TotHits = 0;
    double dEdxSum = 0;

    // *********** How do I deal with backwards tracks....What do I do if only one is backwards?.....
    //  If the dEdx is the wrong way around then without truth I would assume that the track is backwards.
    //  This means that I should use whether the MC is correct as a later cut.
    //  So if MCCorOrient == false then get PIDA using start of 'track'
    for (int Plane = 0; Plane < (int)calos.size(); ++Plane)
    { // Loop through planes
        double PlanePIDA = 0;
        int PlaneHits = 0;
        for (int PlaneHit = 0; PlaneHit < (int)calos[Plane]->dEdx().size(); ++PlaneHit)
        { // loop through hits on each plane
            double ThisdEdx = 0;
            double ThisResR = 0;
            if (IsInverted)
            {
                ThisdEdx = calos[Plane]->dEdx()[PlaneHit];
                ThisResR = calos[Plane]->ResidualRange()[(int)calos[Plane]->dEdx().size() - 1 - PlaneHit];
            }
            else
            {
                ThisdEdx = calos[Plane]->dEdx()[PlaneHit];
                ThisResR = calos[Plane]->ResidualRange()[PlaneHit];
            }
            // Increment TotHits and dEdx sum
            dEdxSum += ThisdEdx;
            ++TotHits;
            // ==== If MCCorOrient == true
            // Work out PIDA if ResRange < 30 cm
            if (ThisResR < 30)
            { // Only want PIDA for last 30 cm
                PlanePIDA += PIDAFunc(ThisdEdx, ThisResR);
                ++PlaneHits;
            } // If ResRange < 30 cm
              // ===== This is where I need to do things...
        }     // Loop over hits on each plane
              // Increment whole track PIDA.
        PIDA += PlanePIDA;
        UsedHits += PlaneHits;
        // Work out PIDA for this plane
        PlanePIDA = PlanePIDA / PlaneHits;
    } // Loop over planes

    if (UsedHits) // If had any hits, work out PIDA and calculate
        PIDA = PIDA / UsedHits;
    if (PIDA > 60)
        PIDA = 60;
    // AvdEdx = dEdxSum / TotHits;
    return PIDA;
} // CalcPIDA

double axion::AxionAnalyzer::PIDAFunc(double dedx, double resrng)
{
    double Va = dedx * pow(resrng, 0.42);
    return Va;
}

void axion::AxionAnalyzer::beginJob()
{

    std::cout << " axion::AxionAnalyzer::beginJob() - initializing..." << std::endl;

    art::ServiceHandle<art::TFileService const> tfs;
    m_AtmTree = tfs->make<TTree>("Atm", "AtmosphericAnalysis");
    m_AllEvents = tfs->make<TTree>("AllEvents", "AllEvents");

    m_AtmTree->Branch("CCNC", &fCCNC);
    m_AtmTree->Branch("run", &m_run, "run/I");
    m_AtmTree->Branch("event", &m_event, "event/I");
    m_AtmTree->Branch("LongestTrack", &fLongestTrack);

    m_AtmTree->Branch("NPrimaryDaughters", &fNPrimaryDaughters);
    m_AtmTree->Branch("NPrimaries", &fNPrimaries);
    m_AtmTree->Branch("DaughterTrackID", &fDaughterTrackID);
    m_AtmTree->Branch("CVN_NCScore", &fCVN_NCScore);
    m_AtmTree->Branch("CVN_0protonsProbability", &fCVN_0protonsProbability);
    m_AtmTree->Branch("CVN_1protonsProbability", &fCVN_1protonsProbability);
    m_AtmTree->Branch("CVN_2protonsProbability", &fCVN_2protonsProbability);
    m_AtmTree->Branch("CVN_NprotonsProbability", &fCVN_NprotonsProbability);
    m_AtmTree->Branch("CVN_0pionsProbability", &fCVN_0pionsProbability);
    m_AtmTree->Branch("CVN_0pizerosProbability", &fCVN_0pizerosProbability);
    m_AtmTree->Branch("CVN_0neutronsProbability", &fCVN_0neutronsProbability);
    m_AtmTree->Branch("CVN_IsAntineutrinoProbability", &fCVN_IsAntineutrinoProbability);
    m_AtmTree->Branch("CVN_NumuProbability", &fCVN_NumuProbability);
    m_AtmTree->Branch("CVN_NueProbability", &fCVN_NueProbability);
    m_AtmTree->Branch("CVN_NutauProbability", &fCVN_NutauProbability);
    m_AtmTree->Branch("IsNC_CVNPred", &fIsNC_CVNPred);
    m_AtmTree->Branch("CVN_FlavourPred", &fCVN_FlavourPred);
    m_AtmTree->Branch("CVN_NProtonsPred", &fCVN_NProtonsPred);
    m_AtmTree->Branch("CVN_AntineutrinoPred", &fCVN_AntineutrinoPred);
    m_AtmTree->Branch("CVN_NPionsPred", &fCVN_NPionsPred);
    m_AtmTree->Branch("CVN_NPions0Pred", &fCVN_NPions0Pred);
    m_AtmTree->Branch("MCGammaE", &fMCGammaE);

    m_AtmTree->Branch("MCKineticEnergy", &fMCKineticEnergy);
    m_AtmTree->Branch("TrackStartX", &fTrackStartX);
    m_AtmTree->Branch("TrackStartY", &fTrackStartY);
    m_AtmTree->Branch("TrackStartZ", &fTrackStartZ);
    m_AtmTree->Branch("TrackEndX", &fTrackEndX);
    m_AtmTree->Branch("TrackEndY", &fTrackEndY);
    m_AtmTree->Branch("TrackEndZ", &fTrackEndZ);

    /*
      m_AtmTree->Branch("RecoMethodUsed_NCEnergy", &fRecoMethodUsed_NCEnergy);
      m_AtmTree->Branch("RecoVertex_NCEnergy", &fRecoVertex_NCEnergy);
      m_AtmTree->Branch("NuLorentzVector_NCEnergy", &fNuLorentzVector_NCEnergy);
      m_AtmTree->Branch("LepLorentzVector_NCEnergy", &fLepLorentzVector_NCEnergy);
      m_AtmTree->Branch("HadLorentzVector_NCEnergy", &fHadLorentzVector_NCEnergy);
      m_AtmTree->Branch("longestTrackContained_NCEnergy", &flongestTrackContained_NCEnergy);
      m_AtmTree->Branch("trackMomMethod_NCEnergy", &ftrackMomMethod_NCEnergy);
    */

    m_AtmTree->Branch("EventRecoEnergy_Charge", &fEventRecoEnergy);
    //m_AtmTree->Branch("EventRecoEnergy_numu", &fEventRecoEnergy_numu);
    m_AtmTree->Branch("EventRecoEnergy_nue", &fEventRecoEnergy_nue);

    m_AtmTree->Branch("PrimaryRecoVertex", &fPrimaryRecoVertex);
    m_AtmTree->Branch("PIDA_NoFlip", &fPIDANoFlip);
    m_AtmTree->Branch("PIDAwithFlipMaybe", &fPIDAwithFlipMaybe);
    m_AtmTree->Branch("HighestTrackSummedADC", &fHighestTrackSummedADC);
    m_AtmTree->Branch("HighestShowerSummedADC", &fHighestShowerSummedADC);
    m_AtmTree->Branch("PIDALongestTrack", &fPIDALongestTrack);
    m_AtmTree->Branch("IsTrackFlipped", &fIsTrackFlipped);
    m_AtmTree->Branch("PrimaryPDGReco", &fPrimaryPDGReco);
    m_AtmTree->Branch("LargeShowerOpenAngle", &fLargeShowerOpenAngle);
    m_AtmTree->Branch("LongestShower", &fLongestShower);

    m_AtmTree->Branch("DiffCosAngleTotalMom", &fDiffCosAngleTotalMom);
    m_AtmTree->Branch("DiffCosAngleLongestTrack", &fDiffCosAngleLongestTrack);
    m_AtmTree->Branch("FracTotalChargeLongTrack", &fFracTotalChargeLongTrack);
    m_AtmTree->Branch("AvarageTrackLength", &fAvarageTrackLength);

    m_AtmTree->Branch("ShowerID", &fShowerID);
    m_AtmTree->Branch("EnergyShowerLinearlyCorrected", &fEnergyShowerLinearlyCorrected);
    m_AtmTree->Branch("ShowerDirectionX", &fShowerDirectionX);
    m_AtmTree->Branch("ShowerDirectionY", &fShowerDirectionY);
    m_AtmTree->Branch("ShowerDirectionZ", &fShowerDirectionZ);
    m_AtmTree->Branch("ShowerDirectionErr", &fShowerDirectionErr);
    m_AtmTree->Branch("ShowerShowerStart", &fShowerShowerStart);
    m_AtmTree->Branch("ShowerShowerStartErr", &fShowerShowerStartErr);
    m_AtmTree->Branch("Showerbest_plane", &fShowerbest_plane);
    m_AtmTree->Branch("ShowerEnergy", &fShowerEnergy);
    m_AtmTree->Branch("ShowerLength", &fShowerLength);
    m_AtmTree->Branch("ShowerOpenAngle", &fShowerOpenAngle);

    m_AtmTree->Branch("TotalMomRecoRangeUnitVect", &fTotalMomRecoRangeUnitVect);
    m_AtmTree->Branch("TotalMomRecoCalVectUnit", &fTotalMomRecoCalVectUnit);
    m_AtmTree->Branch("CosThetaDetTotalMom", &fCosThetaDetTotalMom);
    m_AtmTree->Branch("CosPhiDetTotalMom", &fCosPhiDetTotalMom);
    m_AtmTree->Branch("nTracks", &fnTracks);
    m_AtmTree->Branch("nShowers", &fnShowers);
    m_AtmTree->Branch("TotalMomentumP", &fTotalMomentumP);
    m_AtmTree->Branch("nSpacePoints", &fnSpacePoints);
    m_AtmTree->Branch("minChi2value", &fminChi2value);
    m_AtmTree->Branch("minChi2PDG", &fminChi2PDG);
    m_AtmTree->Branch("DistVertex", &fDistVertex);
    m_AtmTree->Branch("NHits", &fNHits);
    /*
      m_AtmTree->Branch("RecoTrackTruePDG", &fRecoTrackTruePDG);
      m_AtmTree->Branch("RecoTrackTruePrimary", &fRecoTrackTruePrimary);
      m_AtmTree->Branch("RecoTrackTrueMomX", &fRecoTrackTrueMomX);
      m_AtmTree->Branch("RecoTrackTrueMomY", &fRecoTrackTrueMomY);
      m_AtmTree->Branch("RecoTrackTrueMomZ", &fRecoTrackTrueMomZ);
      m_AtmTree->Branch("RecoTrackTrueMomT", &fRecoTrackTrueMomT);
      m_AtmTree->Branch("RecoTrackTrueStartX", &fRecoTrackTrueStartX);
      m_AtmTree->Branch("RecoTrackTrueStartY", &fRecoTrackTrueStartY);
      m_AtmTree->Branch("RecoTrackTrueStartZ", &fRecoTrackTrueStartZ);
      m_AtmTree->Branch("RecoTrackTrueStartT", &fRecoTrackTrueStartT);
      m_AtmTree->Branch("RecoTrackTrueEndX", &fRecoTrackTrueEndX);
      m_AtmTree->Branch("RecoTrackTrueEndY", &fRecoTrackTrueEndY);
      m_AtmTree->Branch("RecoTrackTrueEndZ", &fRecoTrackTrueEndZ);
      m_AtmTree->Branch("RecoTrackTrueEndT", &fRecoTrackTrueEndT);

      m_AtmTree->Branch("RecoTrackMVAEvalRatio", &fRecoTrackMVAEvalRatio);
      m_AtmTree->Branch("RecoTrackMVAConcentration", &fRecoTrackMVAConcentration);
      m_AtmTree->Branch("RecoTrackMVACoreHaloRatio", &fRecoTrackMVACoreHaloRatio);
      m_AtmTree->Branch("RecoTrackMVAConicalness", &fRecoTrackMVAConicalness);
      m_AtmTree->Branch("RecoTrackMVAdEdxStart", &fRecoTrackMVAdEdxStart);
      m_AtmTree->Branch("RecoTrackMVAdEdxEnd", &fRecoTrackMVAdEdxEnd);
      m_AtmTree->Branch("RecoTrackMVAdEdxEndRatio", &fRecoTrackMVAdEdxEndRatio);
      m_AtmTree->Branch("RecoTrackMVAElectron", &fRecoTrackMVAElectron);
      m_AtmTree->Branch("RecoTrackMVAPion", &fRecoTrackMVAPion);
      m_AtmTree->Branch("RecoTrackMVAMuon", &fRecoTrackMVAMuon);
      m_AtmTree->Branch("RecoTrackMVAProton", &fRecoTrackMVAProton);
      m_AtmTree->Branch("RecoTrackMVAPhoton", &fRecoTrackMVAPhoton);
    */

    // m_AtmTree->Branch("SunDirectionFromTrueBDM", &fSunDirectionFromTrueBDM);
    // m_AtmTree->Branch("TruePrimaryBDMVertex", &fPrimaryBDMVertex);

    m_AtmTree->Branch("MCPartGenPDG", &fMCPartGenPDG);
    m_AtmTree->Branch("MCGammaGenMomentum", &fMCGammaGenMomentum);
    m_AtmTree->Branch("MCGammaGenEndMomentum", &fMCGammaGenEndMomentum);
    m_AtmTree->Branch("MCPartGenStatusCode", &fMCPartGenStatusCode);
    //m_AtmTree->Branch("ThetaNuLepton", &fThetaNuLepton);
    m_AtmTree->Branch("MCPrimaryGammaPDG", &fMCPrimaryGammaPDG);
    //m_AtmTree->Branch("MCNuMomentum", &fMCNuMomentum);
    m_AtmTree->Branch("MCInitialPositionGamma", &fMCInitialPositionGamma);
    m_AtmTree->Branch("MCCosAzimuthGamma", &fMCCosAzimuthGamma);
    m_AtmTree->Branch("nGeantParticles", &fnGeantParticles);
    m_AtmTree->Branch("MCTrackId", &fMCTrackId);
    m_AtmTree->Branch("MCPdgCode", &fMCPdgCode);
    m_AtmTree->Branch("MCProcess", &fMCProcess);
    m_AtmTree->Branch("MCMomentum", &fMCMomentum);
    m_AtmTree->Branch("MCStartEnergy", &fMCStartEnergy);
    m_AtmTree->Branch("MCStatusCode", &fMCStatusCode);
    m_AtmTree->Branch("isMCinside", &fisMCinside);
    m_AtmTree->Branch("TotalMomentumUnitVect", &fTotalMomentumUnitVect);
    m_AtmTree->Branch("TopologyNProton", &fTopologyNProton);
    m_AtmTree->Branch("TotalMomentumTrueMag", &fTotalMomentumTrueMag);

    m_AtmTree->Branch("TrksMom_AllProtons", &fTrksMom_AllProtons);
    m_AtmTree->Branch("TrksMom_AllMuons", &fTrksMom_AllMuons);
    m_AtmTree->Branch("TrksMom_BestFit", &fTrksMom_BestFit);
    m_AtmTree->Branch("TrksMom_MCS", &fTrksMom_MCS);

    m_AtmTree->Branch("DiffCosAngleTotalMom_AllProtons", &fDiffCosAngleTotalMom_AllProtons);
    m_AtmTree->Branch("TotalMomentumP_AllProtons", &fTotalMomentumP_AllProtons);
    m_AtmTree->Branch("CosThetaDetTotalMom_AllProtons", &fCosThetaDetTotalMom_AllProtons);
    m_AtmTree->Branch("CosPhiDetTotalMom_AllProtons", &fCosPhiDetTotalMom_AllProtons);
    m_AtmTree->Branch("TotalMomRecoRangeUnitVect_AllProtons", &fTotalMomRecoRangeUnitVect_AllProtons);

    m_AtmTree->Branch("DiffCosAngleTotalMom_AllMuons", &fDiffCosAngleTotalMom_AllMuons);
    m_AtmTree->Branch("TotalMomentumP_AllMuons", &fTotalMomentumP_AllMuons);
    m_AtmTree->Branch("CosThetaDetTotalMom_AllMuons", &fCosThetaDetTotalMom_AllMuons);
    m_AtmTree->Branch("CosPhiDetTotalMom_AllMuons", &fCosPhiDetTotalMom_AllMuons);
    m_AtmTree->Branch("TotalMomRecoRangeUnitVect_AllMuons", &fTotalMomRecoRangeUnitVect_AllMuons);

    m_AtmTree->Branch("DiffCosAngleTotalMom_MCS", &fDiffCosAngleTotalMom_MCS);
    m_AtmTree->Branch("TotalMomentumP_MCS", &fTotalMomentumP_MCS);
    m_AtmTree->Branch("CosThetaDetTotalMom_MCS", &fCosThetaDetTotalMom_MCS);
    m_AtmTree->Branch("CosPhiDetTotalMom_MCS", &fCosPhiDetTotalMom_MCS);
    m_AtmTree->Branch("TotalMomMCSUnitVect", &fTotalMomMCSUnitVect);

    m_AtmTree->Branch("nGeantParticles_Primaries", &fnGeantParticles_Primaries);
    m_AtmTree->Branch("nGenParticles", &fnGenParticles);

    m_AllEvents->Branch("event", &m_event, "event/I");
    m_AllEvents->Branch("nGeantParticles_Primaries", &fnGeantParticles_Primaries);
    m_AllEvents->Branch("nGenParticles", &fnGenParticles);
    m_AllEvents->Branch("CCNC", &fCCNC);
   // m_AllEvents->Branch("ThetaNuLepton", &fThetaNuLepton);
    m_AllEvents->Branch("MCPrimaryNuPDG", &fMCPrimaryGammaPDG);
   // m_AllEvents->Branch("MCNuMomentum", &fMCNuMomentum);
    m_AllEvents->Branch("MCInitialPositionGamma", &fMCInitialPositionGamma);
    m_AllEvents->Branch("MCTrackId", &fMCTrackId);
    m_AllEvents->Branch("MCPdgCode", &fMCPdgCode);
    m_AllEvents->Branch("MCProcess", &fMCProcess);
    m_AllEvents->Branch("MCMomentum", &fMCMomentum);
    m_AllEvents->Branch("MCCosAzimuthGamma", &fMCCosAzimuthGamma);
    m_AllEvents->Branch("MCStartEnergy", &fMCStartEnergy);
    m_AllEvents->Branch("MCKineticEnergy", &fMCKineticEnergy);
    m_AllEvents->Branch("isMCinside", &fisMCinside);
    m_AllEvents->Branch("TopologyNProton", &fTopologyNProton);
    m_AllEvents->Branch("MCKineticEnergy", &fMCKineticEnergy);
    m_AllEvents->Branch("TotalMomentumTrueMag", &fTotalMomentumTrueMag);
    m_AllEvents->Branch("NHits", &fNHits);

    //  std::cout << " axion::AxionAnalyzer::beginJob() - End" << std::endl;
}

void axion::AxionAnalyzer::endJob()
{
    // Implementation of optional member function here.
}

bool axion::AxionAnalyzer::IsVisibleParticle(int Pdg, std::string process)
{

    Pdg = abs(Pdg);
    bool condition = false;

    if ((Pdg == 130 || Pdg == 130 || Pdg == 211 || (Pdg > 300 && Pdg < 400) || Pdg == 2212 || (Pdg > 3000 && Pdg < 4000) || Pdg == 22 || Pdg == 13 || Pdg == 11) && process == "primary")
    {
        condition = true;
    }

    return condition;
}

int axion::AxionAnalyzer::IsGammaAxion(int Pdg, int status_code)
{
    Pdg = abs(Pdg);
    int condition = 0;
    if (!(Pdg == 22))
        return condition;

    if (status_code == 0)
        condition = 1;
    if (status_code == 1 || status_code == 15)
        condition = 2;

    return condition;
}

void axion::AxionAnalyzer::GeoLimits(art::ServiceHandle<geo::Geometry const> &geom, float fFidVolCutX, float fFidVolCutY, float fFidVolCutZ)
{

    // Define histogram boundaries (cm).
    // For now only draw cryostat=0.
    /*double minx = 1e9;
    double maxx = -1e9;
    double miny = 1e9;
    double maxy = -1e9;
    double minz = 1e9;
    double maxz = -1e9;*/

    const geo::BoxBoundedGeo &GeoLimitsSize = geom->WorldBox();

    fFidVolXmin = GeoLimitsSize.MinX() + fFidVolCutX;
    fFidVolXmax = GeoLimitsSize.MaxX() - fFidVolCutX;
    fFidVolYmin = GeoLimitsSize.MinY() + fFidVolCutY;
    fFidVolYmax = GeoLimitsSize.MaxY() - fFidVolCutY;
    fFidVolZmin = GeoLimitsSize.MinZ() + fFidVolCutZ;
    fFidVolZmax = GeoLimitsSize.MaxX() - fFidVolCutZ;
    /*for (unsigned int i = 0; i < geom->NTPC(); ++i)
    {
      double local[3] = {0., 0., 0.};
      double world[3] = {0., 0., 0.};
      const geo::TPCGeo &tpc = geom->TPC(i,j);
      tpc.LocalToWorld(local, world);
      if (minx > world[0] - geom->DetHalfWidth(i))
        minx = world[0] - geom->DetHalfWidth(i);
      if (maxx < world[0] + geom->DetHalfWidth(i))
        maxx = world[0] + geom->DetHalfWidth(i);
      if (miny > world[1] - geom->DetHalfHeight(i))
        miny = world[1] - geom->DetHalfHeight(i);
      if (maxy < world[1] + geom->DetHalfHeight(i))
        maxy = world[1] + geom->DetHalfHeight(i);
      if (minz > world[2] - geom->DetLength(i) / 2.)
        minz = world[2] - geom->DetLength(i) / 2.;
      if (maxz < world[2] + geom->DetLength(i) / 2.)
        maxz = world[2] + geom->DetLength(i) / 2.;
    }*/
    /*
        fFidVolXmin = minx + fFidVolCutX;
        fFidVolXmax = maxx - fFidVolCutX;
        fFidVolYmin = miny + fFidVolCutY;
        fFidVolYmax = maxy - fFidVolCutY;
        fFidVolZmin = minz + fFidVolCutZ;
        fFidVolZmax = maxz - fFidVolCutZ;*/

} // GeoLimits()

bool axion::AxionAnalyzer::insideFV(geo::Point_t const &vertex)
{

    double const x = vertex.X();
    double const y = vertex.Y();
    double const z = vertex.Z();

    return x > fFidVolXmin && x < fFidVolXmax &&
           y > fFidVolYmin && y < fFidVolYmax &&
           z > fFidVolZmin && z < fFidVolZmax;

} // insideFV()

int axion::AxionAnalyzer::Topology(std::vector<int> fMCPdgCode)
{
    int ProtonCnt = -1;
    int PionPCnt = -1;
    int PionNCnt = -1;
    int Pion0Cnt = -1;
    ProtonCnt = std::count(fMCPdgCode.begin(), fMCPdgCode.end(), 2212);
    PionPCnt = std::count(fMCPdgCode.begin(), fMCPdgCode.end(), 211);
    PionNCnt = std::count(fMCPdgCode.begin(), fMCPdgCode.end(), -211);
    Pion0Cnt = std::count(fMCPdgCode.begin(), fMCPdgCode.end(), 111);
    if (PionPCnt == 0 && Pion0Cnt == 0 && PionNCnt == 0)
    {
        return ProtonCnt;
    }
    else
    {
        return -1;
    }
} // Topology

DEFINE_ART_MODULE(axion::AxionAnalyzer)