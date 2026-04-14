// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

/// \file correlationsEse.cxx
/// \brief Same- and mixed-event two-particle correlations on derived CF data with spherocity as CorrelationContainer user axis
/// \author Neelkamal Mallick (neelkamal.mallick@cern.ch)

#include "PWGCF/Core/CorrelationContainer.h"
#include "PWGCF/Core/PairCuts.h"
#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/CFSpherocity.h"

#include "CCDB/BasicCCDBManager.h"
#include "Common/Core/RecoDecay.h"
#include "Common/DataModel/TrackSelectionTables.h"
#include "CommonConstants/MathConstants.h"
#include "DataFormatsParameters/GRPObject.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <THn.h>
#include <TFile.h>

#include <chrono>
#include <memory>
#include <type_traits>
#include <vector>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;
using namespace constants::math;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct CorrelationEseTask {
  SliceCache cache;

  O2_DEFINE_CONFIGURABLE(cfgCutVertex, float, 7.0f, "Accepted z-vertex range")
  O2_DEFINE_CONFIGURABLE(cfgCutPt, float, 0.5f, "Minimal pT for tracks")
  O2_DEFINE_CONFIGURABLE(cfgCutEta, float, 0.8f, "Eta range for tracks")

  O2_DEFINE_CONFIGURABLE(cfgPtOrder, int, 1, "Only consider pairs for which pT,1 < pT,2 (0 = OFF, 1 = ON)");
  O2_DEFINE_CONFIGURABLE(cfgTriggerCharge, int, 0, "Select on charge of trigger particle: 0 = all; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgAssociatedCharge, int, 0, "Select on charge of associated particle: 0 = all charged; 1 = positive; -1 = negative");
  O2_DEFINE_CONFIGURABLE(cfgPairCharge, int, 0, "Select on charge of particle pair: 0 = all; 1 = like sign; -1 = unlike sign");

  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCut, float, -1, "Two track cut: -1 = off; >0 otherwise distance value (suggested: 0.02)");
  O2_DEFINE_CONFIGURABLE(cfgTwoTrackCutMinRadius, float, 0.8f, "Two track cut: radius in m from which two track cuts are applied");
  O2_DEFINE_CONFIGURABLE(cfgLocalEfficiency, int, 0, "0 = OFF and 1 = ON for local efficiency");
  O2_DEFINE_CONFIGURABLE(cfgDropStepRECO, bool, false, "choice to drop step RECO if efficiency correction is used")
  O2_DEFINE_CONFIGURABLE(cfgTrackBitMask, uint16_t, 0, "BitMask for track selection systematics; refer to the enum TrackSelectionCuts in filtering task");

  O2_DEFINE_CONFIGURABLE(cfgEfficiencyTrigger, std::string, "", "CCDB path to efficiency object for trigger particles")
  O2_DEFINE_CONFIGURABLE(cfgEfficiencyAssociated, std::string, "", "CCDB path to efficiency object for associated particles")

  O2_DEFINE_CONFIGURABLE(cfgNoMixedEvents, int, 5, "Number of mixed events per event")

  O2_DEFINE_CONFIGURABLE(cfgVerbosity, int, 1, "Verbosity level (0 = major, 1 = per collision)")

  ConfigurableAxis axisVertex{"axisVertex", {7, -7, 7}, "vertex axis for histograms"};
  ConfigurableAxis axisDeltaPhi{"axisDeltaPhi", {72, -PIHalf, PIHalf * 3}, "delta phi axis for histograms"};
  ConfigurableAxis axisDeltaEta{"axisDeltaEta", {40, -2, 2}, "delta eta axis for histograms"};
  ConfigurableAxis axisPtTrigger{"axisPtTrigger", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 10.0}, "pt trigger axis for histograms"};
  ConfigurableAxis axisPtAssoc{"axisPtAssoc", {VARIABLE_WIDTH, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0}, "pt associated axis for histograms"};
  ConfigurableAxis axisMultiplicity{"axisMultiplicity", {VARIABLE_WIDTH, 0, 5, 10, 20, 30, 40, 50, 100.1}, "multiplicity / centrality axis for histograms"};
  ConfigurableAxis axisSpherocity{"axisSpherocity", {20, 0.f, 1.f}, "spherocity user axis for CorrelationContainer"};

  ConfigurableAxis axisVertexEfficiency{"axisVertexEfficiency", {10, -10, 10}, "vertex axis for efficiency histograms"};
  ConfigurableAxis axisEtaEfficiency{"axisEtaEfficiency", {20, -1.0, 1.0}, "eta axis for efficiency histograms"};
  ConfigurableAxis axisPtEfficiency{"axisPtEfficiency", {VARIABLE_WIDTH, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.25, 1.5, 1.75, 2.0, 2.25, 2.5, 2.75, 3.0, 3.25, 3.5, 3.75, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0}, "pt axis for efficiency histograms"};

  static constexpr float kCfgPairCutDefaults[1][5] = {{-1, -1, -1, -1, -1}};
  Configurable<LabeledArray<float>> cfgPairCut{"cfgPairCut", {kCfgPairCutDefaults[0], 5, {"Photon", "K0", "Lambda", "Phi", "Rho"}}, "Pair cuts on various particles"};

  Filter collisionZVtxFilter = nabs(aod::collision::posZ) < cfgCutVertex;
  Filter cfTrackFilter = (nabs(aod::cftrack::eta) < cfgCutEta) && (aod::cftrack::pt > cfgCutPt) && ncheckbit(aod::track::trackType, as<uint8_t>(cfgTrackBitMask));

  OutputObj<CorrelationContainer> same{"sameEvent"};
  OutputObj<CorrelationContainer> mixed{"mixedEvent"};

  std::vector<float> efficiencyAssociatedCache;

  struct Config {
    bool mPairCuts = false;
    THn* mEfficiencyTrigger = nullptr;
    THn* mEfficiencyAssociated = nullptr;
    bool efficiencyLoaded = false;
  } cfg;

  HistogramRegistry registry{"registry"};
  PairCuts mPairCuts;

  Service<o2::ccdb::BasicCCDBManager> ccdb;

  using DerivedCollisionsEse = soa::Filtered<soa::Join<aod::CFCollisions, aod::CFSpherocity>>;
  using DerivedTracks = soa::Filtered<aod::CFTracks>;

  void init(o2::framework::InitContext&)
  {
    registry.add("multiplicity", "event multiplicity", {HistType::kTH1F, {{1000, 0, 100, "/multiplicity/centrality"}}});
    registry.add("yields", "multiplicity/centrality vs pT vs eta", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {40, 0, 20, "p_{T}"}, {100, -2, 2, "#eta"}}});
    registry.add("etaphi", "multiplicity/centrality vs eta vs phi", {HistType::kTH3F, {{100, 0, 100, "/multiplicity/centrality"}, {100, -2, 2, "#eta"}, {200, 0, o2::constants::math::TwoPI, "#varphi"}}});

    const int maxMixBin = AxisSpec(axisMultiplicity).getNbins() * AxisSpec(axisVertex).getNbins();
    registry.add("eventcount_same", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("eventcount_mixed", "bin", {HistType::kTH1F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}}});
    registry.add("trackcount_same", "bin", {HistType::kTH2F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {10, -0.5, 9.5}}});
    registry.add("trackcount_mixed", "bin", {HistType::kTH3F, {{maxMixBin + 2, -2.5, -0.5 + maxMixBin, "bin"}, {10, -0.5, 9.5}, {10, -0.5, 9.5}}});

    mPairCuts.SetHistogramRegistry(&registry);

    if (cfgPairCut->get("Photon") > 0 || cfgPairCut->get("K0") > 0 || cfgPairCut->get("Lambda") > 0 || cfgPairCut->get("Phi") > 0 || cfgPairCut->get("Rho") > 0) {
      mPairCuts.SetPairCut(PairCuts::Photon, cfgPairCut->get("Photon"));
      mPairCuts.SetPairCut(PairCuts::K0, cfgPairCut->get("K0"));
      mPairCuts.SetPairCut(PairCuts::Lambda, cfgPairCut->get("Lambda"));
      mPairCuts.SetPairCut(PairCuts::Phi, cfgPairCut->get("Phi"));
      mPairCuts.SetPairCut(PairCuts::Rho, cfgPairCut->get("Rho"));
      cfg.mPairCuts = true;
    }

    if (cfgTwoTrackCut > 0) {
      mPairCuts.SetTwoTrackCuts(cfgTwoTrackCut, cfgTwoTrackCutMinRadius);
    }

    std::vector<AxisSpec> corrAxis = {{axisDeltaEta, "#Delta#eta"},
                                      {axisPtAssoc, "p_{T} (GeV/c)"},
                                      {axisPtTrigger, "p_{T} (GeV/c)"},
                                      {axisMultiplicity, "multiplicity / centrality"},
                                      {axisDeltaPhi, "#Delta#varphi (rad)"},
                                      {axisVertex, "z-vtx (cm)"}};
    std::vector<AxisSpec> effAxis = {{axisEtaEfficiency, "#eta"},
                                     {axisPtEfficiency, "p_{T} (GeV/c)"},
                                     {axisVertexEfficiency, "z-vtx (cm)"}};
    std::vector<AxisSpec> userAxis = {{axisSpherocity, "S_{0}"}};

    same.setObject(new CorrelationContainer("sameEvent", "sameEvent", corrAxis, effAxis, userAxis));
    mixed.setObject(new CorrelationContainer("mixedEvent", "mixedEvent", corrAxis, effAxis, userAxis));

    same->setTrackEtaCut(cfgCutEta);
    mixed->setTrackEtaCut(cfgCutEta);

    if (!cfgEfficiencyAssociated.value.empty()) {
      efficiencyAssociatedCache.reserve(512);
    }

    ccdb->setURL("http://alice-ccdb.cern.ch");
    ccdb->setCaching(true);
    ccdb->setLocalObjectValidityChecking();

    auto now = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
    ccdb->setCreatedNotAfter(now);
  }

  int getMagneticField(uint64_t timestamp)
  {
    static o2::parameters::GRPObject* grpo = nullptr;
    if (grpo == nullptr) {
      grpo = ccdb->getForTimeStamp<o2::parameters::GRPObject>("GLO/GRP/GRP", timestamp);
      if (grpo == nullptr) {
        LOGF(fatal, "GRP object not found for timestamp %llu", timestamp);
        return 0;
      }
      LOGF(info, "Retrieved GRP for timestamp %llu with magnetic field of %d kG", timestamp, grpo->getNominalL3Field());
    }
    return grpo->getNominalL3Field();
  }

  template <typename TTracks>
  void fillQA(float multiplicity, const TTracks& tracks)
  {
    registry.fill(HIST("multiplicity"), multiplicity);
    for (const auto& track1 : tracks) {
      registry.fill(HIST("yields"), multiplicity, track1.pt(), track1.eta());
      registry.fill(HIST("etaphi"), multiplicity, track1.eta(), track1.phi());
    }
  }

  void loadEfficiency(uint64_t timestamp)
  {
    if (cfg.efficiencyLoaded) {
      return;
    }
    if (cfgEfficiencyTrigger.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyTrigger = TFile::Open(cfgEfficiencyTrigger.value.c_str(), "READ");
        cfg.mEfficiencyTrigger = reinterpret_cast<THn*>(fEfficiencyTrigger->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyTrigger = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyTrigger, timestamp);
      }
      if (cfg.mEfficiencyTrigger == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for trigger particles from %s", cfgEfficiencyTrigger.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for trigger particles from %s (%p)", cfgEfficiencyTrigger.value.c_str(), (void*)cfg.mEfficiencyTrigger);
    }
    if (cfgEfficiencyAssociated.value.empty() == false) {
      if (cfgLocalEfficiency > 0) {
        TFile* fEfficiencyAssociated = TFile::Open(cfgEfficiencyAssociated.value.c_str(), "READ");
        cfg.mEfficiencyAssociated = reinterpret_cast<THn*>(fEfficiencyAssociated->Get("ccdb_object"));
      } else {
        cfg.mEfficiencyAssociated = ccdb->getForTimeStamp<THnT<float>>(cfgEfficiencyAssociated, timestamp);
      }
      if (cfg.mEfficiencyAssociated == nullptr) {
        LOGF(fatal, "Could not load efficiency histogram for associated particles from %s", cfgEfficiencyAssociated.value.c_str());
      }
      LOGF(info, "Loaded efficiency histogram for associated particles from %s (%p)", cfgEfficiencyAssociated.value.c_str(), (void*)cfg.mEfficiencyAssociated);
    }
    cfg.efficiencyLoaded = true;
  }

  double getEfficiencyCorrection(THn* eff, float eta, float pt, float multiplicity, float posZ)
  {
    int effVars[4];
    effVars[0] = eff->GetAxis(0)->FindBin(eta);
    effVars[1] = eff->GetAxis(1)->FindBin(pt);
    effVars[2] = eff->GetAxis(2)->FindBin(multiplicity);
    effVars[3] = eff->GetAxis(3)->FindBin(posZ);
    return eff->GetBinContent(effVars);
  }

  template <CorrelationContainer::CFStep step, typename TTarget, typename TTracks1, typename TTracks2>
  void fillCorrelationsEse(TTarget target, TTracks1 const& tracks1, TTracks2 const& tracks2, float multiplicity, float posZ, int magField, float eventWeight, float spherocity)
  {
    if constexpr (step == CorrelationContainer::kCFStepCorrected) {
      if (cfg.mEfficiencyAssociated) {
        efficiencyAssociatedCache.clear();
        efficiencyAssociatedCache.reserve(tracks2.size());
        for (const auto& track : tracks2) {
          efficiencyAssociatedCache.push_back(getEfficiencyCorrection(cfg.mEfficiencyAssociated, track.eta(), track.pt(), multiplicity, posZ));
        }
      }
    }

    for (const auto& track1 : tracks1) {
      if (cfgTriggerCharge != 0 && cfgTriggerCharge * track1.sign() < 0) {
        continue;
      }

      float triggerWeight = eventWeight;
      if constexpr (step == CorrelationContainer::kCFStepCorrected) {
        if (cfg.mEfficiencyTrigger) {
          triggerWeight *= getEfficiencyCorrection(cfg.mEfficiencyTrigger, track1.eta(), track1.pt(), multiplicity, posZ);
        }
      }

      target->getTriggerHist()->Fill(step, track1.pt(), multiplicity, posZ, spherocity, triggerWeight);

      for (const auto& track2 : tracks2) {
        if constexpr (std::is_same_v<TTracks1, TTracks2>) {
          if (track1.globalIndex() == track2.globalIndex()) {
            continue;
          }
        }

        if (cfgPtOrder != 0 && track2.pt() >= track1.pt()) {
          continue;
        }

        if (cfgAssociatedCharge != 0) {
          if (cfgAssociatedCharge * track2.sign() < 0) {
            continue;
          }
        } else if (track2.sign() == 0) {
          continue;
        }

        if (cfgPairCharge != 0 && cfgPairCharge * track1.sign() * track2.sign() < 0) {
          continue;
        }

        if constexpr (std::is_same_v<TTracks1, TTracks2>) {
          if constexpr (step >= CorrelationContainer::kCFStepReconstructed) {
            if (cfg.mPairCuts && mPairCuts.conversionCuts(track1, track2)) {
              continue;
            }
            if (cfgTwoTrackCut > 0 && mPairCuts.twoTrackCut(track1, track2, magField)) {
              continue;
            }
          }
        }

        float associatedWeight = triggerWeight;
        if constexpr (step == CorrelationContainer::kCFStepCorrected) {
          if (cfg.mEfficiencyAssociated) {
            associatedWeight *= efficiencyAssociatedCache[track2.filteredIndex()];
          }
        }

        float deltaPhi = RecoDecay::constrainAngle(track1.phi() - track2.phi(), -o2::constants::math::PIHalf);

        target->getPairHist()->Fill(step, track1.eta() - track2.eta(), track2.pt(), track1.pt(), multiplicity, deltaPhi, posZ, spherocity, associatedWeight);
      }
    }
  }

  template <class CollType, class TTracks1, class TTracks2>
  void processSameDerivedEseT(CollType const& collision, TTracks1 const& tracks1, TTracks2 const& tracks2)
  {
    using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
    BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true};

    const float spherocity = collision.spherocity();

    if (cfgVerbosity > 0) {
      LOGF(info, "processSameDerivedEse: Tracks for collision: %d/%d | Vertex: %.1f | Multiplicity: %.1f | Spherocity: %.4f", tracks1.size(), tracks2.size(), collision.posZ(), collision.multiplicity(), spherocity);
    }
    loadEfficiency(collision.timestamp());

    const auto multiplicity = collision.multiplicity();
    int field = 0;
    if (cfgTwoTrackCut > 0) {
      field = getMagneticField(collision.timestamp());
    }

    int bin = configurableBinningDerived.getBin({collision.posZ(), collision.multiplicity()});
    registry.fill(HIST("eventcount_same"), bin);
    registry.fill(HIST("trackcount_same"), bin, tracks1.size());
    fillQA(multiplicity, tracks1);

    const bool hasEfficiency = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
    const bool fillReco = !(cfgDropStepRECO && hasEfficiency);

    if (fillReco) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepReconstructed);
      fillCorrelationsEse<CorrelationContainer::kCFStepReconstructed>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, spherocity);
    }
    if (hasEfficiency) {
      same->fillEvent(multiplicity, CorrelationContainer::kCFStepCorrected);
      fillCorrelationsEse<CorrelationContainer::kCFStepCorrected>(same, tracks1, tracks2, multiplicity, collision.posZ(), field, 1.0f, spherocity);
    }
  }

  void processSameDerivedEse(DerivedCollisionsEse::iterator const& collision, DerivedTracks const& tracks)
  {
    processSameDerivedEseT(collision, tracks, tracks);
  }
  PROCESS_SWITCH(CorrelationEseTask, processSameDerivedEse, "Process same event on derived data with spherocity user axis", true);

  template <class CollType, typename... TrackTypes>
  void processMixedDerivedEseT(CollType const& collisions, TrackTypes&&... tracks)
  {
    using BinningTypeDerived = ColumnBinningPolicy<aod::collision::PosZ, aod::cfcollision::Multiplicity>;
    BinningTypeDerived configurableBinningDerived{{axisVertex, axisMultiplicity}, true};

    auto tracksTuple = std::make_tuple(std::forward<TrackTypes>(tracks)...);
    using TA = std::tuple_element<0, decltype(tracksTuple)>::type;
    using TB = std::tuple_element<std::tuple_size_v<decltype(tracksTuple)> - 1, decltype(tracksTuple)>::type;
    Pair<CollType, TA, TB, BinningTypeDerived> pairs{configurableBinningDerived, cfgNoMixedEvents, -1, collisions, tracksTuple, &cache};

    for (auto it = pairs.begin(); it != pairs.end(); it++) {
      auto& [collision1, tracks1, collision2, tracks2] = *it;
      const float spherocity = collision1.spherocity();
      float eventWeight = 1.0f / it.currentWindowNeighbours();
      int field = 0;
      if (cfgTwoTrackCut > 0) {
        field = getMagneticField(collision1.timestamp());
      }

      if (cfgVerbosity > 0) {
        int bin = configurableBinningDerived.getBin({collision1.posZ(), collision1.multiplicity()});
        LOGF(info, "processMixedDerivedEse: Mixed collisions bin: %d pair: [%d, %d] %d (%.3f, %.3f, S0=%.4f), %d (%.3f, %.3f)", bin, it.isNewWindow(), it.currentWindowNeighbours(), collision1.globalIndex(), collision1.posZ(), collision1.multiplicity(), spherocity, collision2.globalIndex(), collision2.posZ(), collision2.multiplicity());
      }

      bool hasEfficiencyMixed = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
      bool fillRecoMixed = !(cfgDropStepRECO && hasEfficiencyMixed);

      if (it.isNewWindow()) {
        loadEfficiency(collision1.timestamp());
        hasEfficiencyMixed = (cfg.mEfficiencyAssociated != nullptr || cfg.mEfficiencyTrigger != nullptr);
        fillRecoMixed = !(cfgDropStepRECO && hasEfficiencyMixed);

        if (fillRecoMixed) {
          mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepReconstructed);
        }
      }

      int bin = configurableBinningDerived.getBin({collision1.posZ(), collision1.multiplicity()});
      registry.fill(HIST("eventcount_mixed"), bin);
      registry.fill(HIST("trackcount_mixed"), bin, tracks1.size(), tracks2.size());

      if (fillRecoMixed) {
        fillCorrelationsEse<CorrelationContainer::kCFStepReconstructed>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight, spherocity);
      }

      if (hasEfficiencyMixed) {
        if (it.isNewWindow()) {
          mixed->fillEvent(collision1.multiplicity(), CorrelationContainer::kCFStepCorrected);
        }
        fillCorrelationsEse<CorrelationContainer::kCFStepCorrected>(mixed, tracks1, tracks2, collision1.multiplicity(), collision1.posZ(), field, eventWeight, spherocity);
      }
    }
  }

  void processMixedDerivedEse(DerivedCollisionsEse const& collisions, DerivedTracks const& tracks)
  {
    processMixedDerivedEseT(collisions, tracks);
  }
  PROCESS_SWITCH(CorrelationEseTask, processMixedDerivedEse, "Process mixed events on derived data with spherocity user axis", true);
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<CorrelationEseTask>(cfgc)};
}
