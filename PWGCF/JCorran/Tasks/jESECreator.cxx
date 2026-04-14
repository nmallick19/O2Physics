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
///
/// \file jESECreator.cxx
/// \brief Task to create spherocity columns joinable to CFCollisions (reco) and CFMcCollisions (MC + reco on MC axis)
/// \author Neelkamal Mallick (neelkamal.mallick@cern.ch)
/// \since Apr 2026
///

#include "PWGCF/DataModel/CorrelationsDerived.h"
#include "PWGCF/JCorran/DataModel/CFSpherocity.h"

#include "Common/DataModel/Centrality.h"
#include "Common/DataModel/EventSelection.h"
#include "Common/Core/TrackSelection.h"
#include "Common/DataModel/TrackSelectionTables.h"

#include <array>
#include <cmath>
#include <vector>

#include "CommonConstants/MathConstants.h"
#include "Framework/ASoAHelpers.h"
#include "Framework/AnalysisTask.h"
#include "Framework/HistogramRegistry.h"
#include "Framework/RunningWorkflowInfo.h"
#include "Framework/runDataProcessing.h"

#include <TH1.h>

using namespace o2;
using namespace o2::framework;
using namespace o2::framework::expressions;

#define O2_DEFINE_CONFIGURABLE(NAME, TYPE, DEFAULT, HELP) Configurable<TYPE> NAME{#NAME, DEFAULT, HELP};

struct EventShapeProducer {
  Produces<aod::CFSpherocity> outputSpherocityColl;
  Produces<aod::CFSpherocityMc> outputSpherocityMc;

  O2_DEFINE_CONFIGURABLE(cfgZvtxMax, float, 10.0, "Maximum z-vertex range for collisions");
  O2_DEFINE_CONFIGURABLE(minTracks, int, 10, "Minimum number of tracks/particles per collision required to estimate spherocity");
  O2_DEFINE_CONFIGURABLE(minPtTrack, float, 0.2, "Minimum pT of tracks/particles required to estimate spherocity");
  O2_DEFINE_CONFIGURABLE(AbsEtaMax, float, 0.8, "Maximum |eta| of tracks/particles required to estimate spherocity");
  O2_DEFINE_CONFIGURABLE(usePtWeighted, bool, false, "Use pt-weighted spherocity");
  O2_DEFINE_CONFIGURABLE(cfgTrackBitMask, uint16_t, 0, "Track selection bitmask to use as defined in the filterCorrelations.cxx task");
  O2_DEFINE_CONFIGURABLE(outputQAHistos, bool, false, "Whether to output QA histograms for spherocity");
  O2_DEFINE_CONFIGURABLE(mcUsePrimariesOnly, bool, true, "For MC spherocity: use only charged physical-primary particles");
  O2_DEFINE_CONFIGURABLE(cfgEseMinEntries, int, 200, "Minimum entries in per-mult s0 reference before applying 20/80%% spherocity cuts in processSpherocityEse");

  static constexpr float kInvalidRecoOnMcAxis = -2.f;
  static constexpr int kMultEseNBins = 8;
  static constexpr std::array<float, 9> kMultEseEdges{0.f, 10.f, 20.f, 30.f, 40.f, 50.f, 60.f, 70.f, 100.f};

  HistogramRegistry QAhistos{"QAhistos", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAhistosCent{"QAhistosCent", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAhistosSpherocity{"QAhistosSpherocity", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};
  HistogramRegistry QAhistosEse{"EseQA", {}, OutputObjHandlingPolicy::AnalysisObject, false, true};

  std::vector<float> nx, ny, px, py;

  float minimizeSpherocity(size_t nTracks, float sumPt, bool usePtWeightedFlag)
  {
    if (nTracks < static_cast<size_t>(minTracks))
      return -2.0f;

    if (sumPt == 0.0f)
      return -1.5f;

    if (!usePtWeightedFlag && std::abs(static_cast<float>(nTracks) - sumPt) > 0.01f) {
      LOGF(info, "Spherocity calculation: number of tracks (%zu) does not match sum of weights (%f)", nTracks, sumPt);
      return -1.0f;
    }

    float retval = 999.0f;

    for (size_t i = 0; i < nTracks; i++) {
      float numerator = 0.0f;
      const float nyi = ny[i];
      const float nxi = nx[i];

      for (size_t j = 0; j < nTracks; j++) {
        numerator += std::abs(nyi * px[j] - nxi * py[j]);
      }

      const float sFull = std::pow(numerator / sumPt, 2);
      if (sFull < retval)
        retval = sFull;
    }

    retval = retval * M_PI * M_PI / 4.0f;

    if (retval < 0.0f || retval > 1.0f) {
      LOGF(info, "Spherocity value is out of range: %f", retval);
      return -0.5f;
    }

    return retval;
  }

  template <typename TTrack>
  bool trackPassesCfSelection(TTrack const& track)
  {
    const float minPt = minPtTrack;
    const float etaMax = AbsEtaMax;
    const uint16_t maskBits = cfgTrackBitMask;
    if (track.pt() <= minPt || std::abs(track.eta()) >= etaMax)
      return false;
    const uint8_t mask = static_cast<uint8_t>(maskBits);
    if (mask > 0 && (static_cast<uint8_t>(track.trackType()) & mask) != mask)
      return false;
    return true;
  }

  /// Spherocity from reconstructed CFTracks (per CFCollision).
  template <typename TTracks>
  float calculateSpherocityReco(TTracks const& tracks, bool usePtWeightedFlag, bool outputQAHistosFlag, bool applyCfCutsInLoop = false,
                                int32_t filterCfCollisionId = -1)
  {
    if (applyCfCutsInLoop) {
      size_t nAcc = 0;
      for (auto const& track : tracks) {
        if (filterCfCollisionId >= 0 && track.cfCollisionId() != filterCfCollisionId)
          continue;
        if (trackPassesCfSelection(track))
          nAcc++;
      }
      if (nAcc < static_cast<size_t>(minTracks))
        return -2.0f;

      if (nx.size() < nAcc) {
        nx.resize(nAcc);
        ny.resize(nAcc);
        px.resize(nAcc);
        py.resize(nAcc);
      }

      float sumPt = 0.0f;
      size_t idx = 0;
      for (auto const& track : tracks) {
        if (filterCfCollisionId >= 0 && track.cfCollisionId() != filterCfCollisionId)
          continue;
        if (!trackPassesCfSelection(track))
          continue;

        const float cosPhi = std::cos(track.phi());
        const float sinPhi = std::sin(track.phi());

        nx[idx] = cosPhi;
        ny[idx] = sinPhi;

        if (usePtWeightedFlag) {
          const float pt = track.pt();
          px[idx] = pt * cosPhi;
          py[idx] = pt * sinPhi;
          sumPt += pt;
        } else {
          px[idx] = cosPhi;
          py[idx] = sinPhi;
          sumPt += 1.0f;
        }
        if (outputQAHistosFlag) {
          QAhistos.fill(HIST("pt"), track.pt());
          QAhistos.fill(HIST("eta"), track.eta());
          QAhistos.fill(HIST("phi"), track.phi());
        }
        idx++;
      }
      return minimizeSpherocity(nAcc, sumPt, usePtWeightedFlag);
    }

    const size_t nTracks = tracks.size();
    if (nTracks < static_cast<size_t>(minTracks))
      return -2.0f;

    if (nx.size() < nTracks) {
      nx.resize(nTracks);
      ny.resize(nTracks);
      px.resize(nTracks);
      py.resize(nTracks);
    }

    float sumPt = 0.0f;
    size_t idx = 0;
    for (auto& track : tracks) {
      const float cosPhi = std::cos(track.phi());
      const float sinPhi = std::sin(track.phi());

      nx[idx] = cosPhi;
      ny[idx] = sinPhi;

      if (usePtWeightedFlag) {
        const float pt = track.pt();
        px[idx] = pt * cosPhi;
        py[idx] = pt * sinPhi;
        sumPt += pt;
      } else {
        px[idx] = cosPhi;
        py[idx] = sinPhi;
        sumPt += 1.0f;
      }
      if (outputQAHistosFlag) {
        QAhistos.fill(HIST("pt"), track.pt());
        QAhistos.fill(HIST("eta"), track.eta());
        QAhistos.fill(HIST("phi"), track.phi());
      }
      idx++;
    }

    return minimizeSpherocity(nTracks, sumPt, usePtWeightedFlag);
  }

  /// Spherocity from MC CFMcParticles for one MC collision (charged primaries optional).
  template <typename TTracksMc>
  float calculateSpherocityMc(TTracksMc const& mcParts, bool usePtWeightedFlag, bool outputQAHistosFlag, size_t* nMcAcceptedOut = nullptr,
                              int32_t filterCfMcCollisionId = -1)
  {
    size_t nAcc = 0;
    for (const auto& p : mcParts) {
      if (filterCfMcCollisionId >= 0 && p.cfMcCollisionId() != filterCfMcCollisionId)
        continue;
      if (mcUsePrimariesOnly && !p.isPhysicalPrimary())
        continue;
      if (p.sign() == 0)
        continue;
      if (p.pt() <= minPtTrack || std::abs(p.eta()) >= AbsEtaMax)
        continue;
      nAcc++;
    }

    if (nMcAcceptedOut != nullptr) {
      *nMcAcceptedOut = nAcc;
    }

    if (nAcc < static_cast<size_t>(minTracks))
      return -2.0f;

    if (nx.size() < nAcc) {
      nx.resize(nAcc);
      ny.resize(nAcc);
      px.resize(nAcc);
      py.resize(nAcc);
    }

    float sumPt = 0.0f;
    size_t idx = 0;
    for (const auto& p : mcParts) {
      if (filterCfMcCollisionId >= 0 && p.cfMcCollisionId() != filterCfMcCollisionId)
        continue;
      if (mcUsePrimariesOnly && !p.isPhysicalPrimary())
        continue;
      if (p.sign() == 0)
        continue;
      if (p.pt() <= minPtTrack || std::abs(p.eta()) >= AbsEtaMax)
        continue;

      const float cosPhi = std::cos(p.phi());
      const float sinPhi = std::sin(p.phi());
      nx[idx] = cosPhi;
      ny[idx] = sinPhi;
      if (usePtWeightedFlag) {
        const float pt = p.pt();
        px[idx] = pt * cosPhi;
        py[idx] = pt * sinPhi;
        sumPt += pt;
      } else {
        px[idx] = cosPhi;
        py[idx] = sinPhi;
        sumPt += 1.0f;
      }
      if (outputQAHistosFlag) {
        QAhistos.fill(HIST("ptMC"), p.pt());
        QAhistos.fill(HIST("etaMC"), p.eta());
        QAhistos.fill(HIST("phiMC"), p.phi());
      }
      idx++;
    }

    return minimizeSpherocity(nAcc, sumPt, usePtWeightedFlag);
  }

  Filter collisionFilter = (nabs(aod::collision::posZ) < cfgZvtxMax);
  Filter mcCollisionFilter = (nabs(aod::mccollision::posZ) < cfgZvtxMax);
  Filter cftrackFilter = (aod::cftrack::pt > minPtTrack) && (nabs(aod::cftrack::eta) < AbsEtaMax) && ncheckbit(aod::track::trackType, as<uint8_t>(cfgTrackBitMask));

  using myCollisions = soa::Filtered<aod::CFCollisions>; // this already contains multiplicity computed at CFFilter level
  using myMcCollisions = soa::Filtered<aod::CFMcCollisions>;
  using myFilteredTracks = soa::Filtered<aod::CFTracks>;

  /// One CFSpherocity row per CFCollision: use unfiltered aod::CFCollisions (not soa::Filtered<>).
  /// Filtered collisions would skip rows and break Join<CFCollisions, CFSpherocity> (e.g. 1096 vs 1058).
  SliceCache cache;
  Preslice<aod::CFTracks> perCfCollision = aod::cftrack::cfCollisionId;

  inline int getCentralityBin(float multiplicity)
  {
    if (multiplicity < 20.0f)
      return 0;
    else if (multiplicity < 40.0f)
      return 1;
    else if (multiplicity < 60.0f)
      return 2;
    else if (multiplicity < 100.0f)
      return 3;
    return -1;
  }

  /// Multiplicity bins for ESE: edges 0,10,...,100 → [0,10), …, [70,100].
  static int getMultEseBinIndex(float mult)
  {
    if (mult < kMultEseEdges[0] || mult > kMultEseEdges[8])
      return -1;
    for (int i = 0; i < kMultEseNBins; ++i) {
      const bool last = (i == kMultEseNBins - 1);
      if (mult >= kMultEseEdges[i] && (last ? mult <= kMultEseEdges[i + 1] : mult < kMultEseEdges[i + 1]))
        return i;
    }
    return -1;
  }

  /// qProb in (0,1): value of x such that the cumulative integral up to that bin reaches qProb * total.
  static float quantileFromTH1(TH1 const* h, double qProb)
  {
    if (h == nullptr || h->GetEntries() < 1)
      return -1.f;
    const double sum = h->Integral();
    if (sum <= 0.)
      return -1.f;
    const double target = qProb * sum;
    double cum = 0.;
    const int nb = h->GetNbinsX();
    for (int i = 1; i <= nb; ++i) {
      cum += h->GetBinContent(i);
      if (cum >= target - 1e-12)
        return static_cast<float>(h->GetBinCenter(i));
    }
    return static_cast<float>(h->GetBinCenter(nb));
  }

  void fillEseS0Dist(int bin, float s0)
  {
    switch (bin) {
      case 0:
        QAhistosEse.fill(HIST("s0_dist_0"), s0);
        break;
      case 1:
        QAhistosEse.fill(HIST("s0_dist_1"), s0);
        break;
      case 2:
        QAhistosEse.fill(HIST("s0_dist_2"), s0);
        break;
      case 3:
        QAhistosEse.fill(HIST("s0_dist_3"), s0);
        break;
      case 4:
        QAhistosEse.fill(HIST("s0_dist_4"), s0);
        break;
      case 5:
        QAhistosEse.fill(HIST("s0_dist_5"), s0);
        break;
      case 6:
        QAhistosEse.fill(HIST("s0_dist_6"), s0);
        break;
      case 7:
        QAhistosEse.fill(HIST("s0_dist_7"), s0);
        break;
      default:
        break;
    }
  }

  void fillEsePtSplow(int bin, float pt)
  {
    switch (bin) {
      case 0:
        QAhistosEse.fill(HIST("pt_0_splow"), pt);
        break;
      case 1:
        QAhistosEse.fill(HIST("pt_1_splow"), pt);
        break;
      case 2:
        QAhistosEse.fill(HIST("pt_2_splow"), pt);
        break;
      case 3:
        QAhistosEse.fill(HIST("pt_3_splow"), pt);
        break;
      case 4:
        QAhistosEse.fill(HIST("pt_4_splow"), pt);
        break;
      case 5:
        QAhistosEse.fill(HIST("pt_5_splow"), pt);
        break;
      case 6:
        QAhistosEse.fill(HIST("pt_6_splow"), pt);
        break;
      case 7:
        QAhistosEse.fill(HIST("pt_7_splow"), pt);
        break;
      default:
        break;
    }
  }

  void fillEsePtInt(int multBin, float pt)
  {
    switch (multBin) {
      case 0:
        QAhistosEse.fill(HIST("pt_0_int"), pt);
        break;
      case 1:
        QAhistosEse.fill(HIST("pt_1_int"), pt);
        break;
      case 2:
        QAhistosEse.fill(HIST("pt_2_int"), pt);
        break;
      case 3:
        QAhistosEse.fill(HIST("pt_3_int"), pt);
        break;
      case 4:
        QAhistosEse.fill(HIST("pt_4_int"), pt);
        break;
      case 5:
        QAhistosEse.fill(HIST("pt_5_int"), pt);
        break;
      case 6:
        QAhistosEse.fill(HIST("pt_6_int"), pt);
        break;
      case 7:
        QAhistosEse.fill(HIST("pt_7_int"), pt);
        break;
      default:
        break;
    }
  }

  void fillEsePtSphigh(int bin, float pt)
  {
    switch (bin) {
      case 0:
        QAhistosEse.fill(HIST("pt_0_sphigh"), pt);
        break;
      case 1:
        QAhistosEse.fill(HIST("pt_1_sphigh"), pt);
        break;
      case 2:
        QAhistosEse.fill(HIST("pt_2_sphigh"), pt);
        break;
      case 3:
        QAhistosEse.fill(HIST("pt_3_sphigh"), pt);
        break;
      case 4:
        QAhistosEse.fill(HIST("pt_4_sphigh"), pt);
        break;
      case 5:
        QAhistosEse.fill(HIST("pt_5_sphigh"), pt);
        break;
      case 6:
        QAhistosEse.fill(HIST("pt_6_sphigh"), pt);
        break;
      case 7:
        QAhistosEse.fill(HIST("pt_7_sphigh"), pt);
        break;
      default:
        break;
    }
  }

  std::shared_ptr<TH1> getEseS0Hist(int bin)
  {
    switch (bin) {
      case 0:
        return QAhistosEse.get<TH1>(HIST("s0_dist_0"));
      case 1:
        return QAhistosEse.get<TH1>(HIST("s0_dist_1"));
      case 2:
        return QAhistosEse.get<TH1>(HIST("s0_dist_2"));
      case 3:
        return QAhistosEse.get<TH1>(HIST("s0_dist_3"));
      case 4:
        return QAhistosEse.get<TH1>(HIST("s0_dist_4"));
      case 5:
        return QAhistosEse.get<TH1>(HIST("s0_dist_5"));
      case 6:
        return QAhistosEse.get<TH1>(HIST("s0_dist_6"));
      case 7:
        return QAhistosEse.get<TH1>(HIST("s0_dist_7"));
      default:
        return nullptr;
    }
  }

  void processSpherocityData(aod::CFCollisions const& collisions, myFilteredTracks const& allTracks)
  {
    for (auto const& coll : collisions) {
      if (std::abs(coll.posZ()) >= cfgZvtxMax) {
        outputSpherocityColl(-2.f);
        continue;
      }
      auto tracks = allTracks.sliceBy(perCfCollision, coll.globalIndex());
      float spherocity = calculateSpherocityReco(tracks, usePtWeighted, outputQAHistos);
      outputSpherocityColl(spherocity);
      const float multiplicity = coll.multiplicity();
      const float posZ = coll.posZ();
      QAhistos.fill(HIST("multS0"), multiplicity, spherocity);
      QAhistos.fill(HIST("spherocity"), spherocity);
      QAhistos.fill(HIST("posZ_coll"), posZ);

      if (!outputQAHistos) {
        continue;
      }

      const size_t nTracks = tracks.size();
      const int centralityBin = getCentralityBin(multiplicity);

      QAhistos.fill(HIST("ntracks"), nTracks);
      QAhistos.fill(HIST("mult"), multiplicity);

      if (centralityBin == 0) {
        QAhistosCent.fill(HIST("ntracks_0_20"), nTracks);
        QAhistosSpherocity.fill(HIST("spherocity_0_20"), spherocity);
      } else if (centralityBin == 1) {
        QAhistosCent.fill(HIST("ntracks_20_40"), nTracks);
        QAhistosSpherocity.fill(HIST("spherocity_20_40"), spherocity);
      } else if (centralityBin == 2) {
        QAhistosCent.fill(HIST("ntracks_40_60"), nTracks);
        QAhistosSpherocity.fill(HIST("spherocity_40_60"), spherocity);
      } else if (centralityBin == 3) {
        QAhistosCent.fill(HIST("ntracks_60_100"), nTracks);
        QAhistosSpherocity.fill(HIST("spherocity_60_100"), spherocity);
      }
    }
  }
  PROCESS_SWITCH(EventShapeProducer, processSpherocityData, "Spherocity from reco tracks; output joinable to CFCollisions", true);

  /// Per bin of kMultEseEdges (8 bins): s0 reference and 20%% / 80%% cuts → splow/sphigh pT.
  /// Integrated pT (int): full spherocity 0–1 (no s0 cuts); same mult bin index as splow/sphigh (kMultEseEdges).
  void processSpherocityEse(aod::CFCollisions const& collisions, myFilteredTracks const& allTracks)
  {
    for (auto const& coll : collisions) {
      if (std::abs(coll.posZ()) >= cfgZvtxMax) {
        continue;
      }
      auto tracks = allTracks.sliceBy(perCfCollision, coll.globalIndex());
      const float spherocity = calculateSpherocityReco(tracks, usePtWeighted, false);
      if (spherocity < 0.f || spherocity > 1.f)
        continue;

      const float multiplicity = coll.multiplicity();
      const int multBin = getMultEseBinIndex(multiplicity);
      if (multBin < 0)
        continue;

      QAhistosEse.fill(HIST("ese_nevent_int"), static_cast<float>(multBin) + 0.5f);

      fillEseS0Dist(multBin, spherocity);

      auto hS0 = getEseS0Hist(multBin);
      if (hS0 == nullptr)
        continue;

      const double q20 = quantileFromTH1(hS0.get(), 0.20);
      const double q80 = quantileFromTH1(hS0.get(), 0.80);
      const bool enoughEntries = hS0->GetEntries() >= cfgEseMinEntries;
      const bool useCuts = enoughEntries && q20 >= 0.f && q80 >= 0.f;

      if (useCuts) {
        if (spherocity <= q20) {
          QAhistosEse.fill(HIST("ese_nevent_splow"), static_cast<float>(multBin) + 0.5f);
        }
        if (spherocity >= q80) {
          QAhistosEse.fill(HIST("ese_nevent_sphigh"), static_cast<float>(multBin) + 0.5f);
        }
      }

      for (auto const& track : tracks) {
        const float pt = track.pt();
        fillEsePtInt(multBin, pt);
        if (useCuts) {
          if (spherocity <= q20)
            fillEsePtSplow(multBin, pt);
          if (spherocity >= q80)
            fillEsePtSphigh(multBin, pt);
        }
      }
    }
  }
  PROCESS_SWITCH(EventShapeProducer, processSpherocityEse, "ESE QA: s0 percentiles per mult bin and pT spectra in EseQA (requires reco spherocity)", false);

  void processSpherocityMC(myMcCollisions::iterator const& mcCollision,
                           aod::CFMcParticles const& mcParticles,
                           soa::SmallGroups<aod::CFCollisionsWithLabel> const& collisions,
                           aod::CFTracks const& tracks)
  {
    const int32_t mcCollId = static_cast<int32_t>(mcCollision.globalIndex());
    size_t nMcAccepted = 0;
    float spherocityMc = calculateSpherocityMc(mcParticles, usePtWeighted, outputQAHistos, &nMcAccepted, mcCollId);
    const float posZ = mcCollision.posZ();
    QAhistos.fill(HIST("spherocityMc"), spherocityMc);
    QAhistos.fill(HIST("posZ_mc"), posZ);

    float spherocityReco = kInvalidRecoOnMcAxis;
    bool hasRecoCollision = false;
    float multColl = 0.f;
    float posZColl = 0.f;
    if (collisions.size() > 0) {
      const float zvtxMax = cfgZvtxMax;
      for (const auto& recoColl : collisions) {
        if (std::abs(recoColl.posZ()) >= zvtxMax)
          continue;
        hasRecoCollision = true;
        multColl = recoColl.multiplicity();
        posZColl = recoColl.posZ();
        const int32_t recoCollId = static_cast<int32_t>(recoColl.globalIndex());
        spherocityReco = calculateSpherocityReco(tracks, usePtWeighted, outputQAHistos, true, recoCollId);
        QAhistos.fill(HIST("spherocityRecoOnMc"), spherocityReco);
        QAhistos.fill(HIST("spherocityMc_vs_spherocityReco"), spherocityMc, spherocityReco);
        QAhistos.fill(HIST("posZ_mc_vs_posZ_reco"), posZ, recoColl.posZ());
        break;
      }
    }

    if (outputQAHistos) {
      QAhistos.fill(HIST("ntracksMC"), static_cast<float>(nMcAccepted));
      QAhistos.fill(HIST("ntracks"), static_cast<float>(nMcAccepted));
    }

    if (hasRecoCollision) {
      QAhistos.fill(HIST("spherocity"), spherocityReco);
      QAhistos.fill(HIST("multS0"), multColl, spherocityReco);
      QAhistos.fill(HIST("posZ_coll"), posZColl);
      if (outputQAHistos) {
        QAhistos.fill(HIST("mult"), multColl);
        const int centralityBin = getCentralityBin(multColl);
        if (centralityBin == 0) {
          QAhistosCent.fill(HIST("ntracks_0_20"), static_cast<float>(nMcAccepted));
          QAhistosSpherocity.fill(HIST("spherocity_0_20"), spherocityReco);
        } else if (centralityBin == 1) {
          QAhistosCent.fill(HIST("ntracks_20_40"), static_cast<float>(nMcAccepted));
          QAhistosSpherocity.fill(HIST("spherocity_20_40"), spherocityReco);
        } else if (centralityBin == 2) {
          QAhistosCent.fill(HIST("ntracks_40_60"), static_cast<float>(nMcAccepted));
          QAhistosSpherocity.fill(HIST("spherocity_40_60"), spherocityReco);
        } else if (centralityBin == 3) {
          QAhistosCent.fill(HIST("ntracks_60_100"), static_cast<float>(nMcAccepted));
          QAhistosSpherocity.fill(HIST("spherocity_60_100"), spherocityReco);
        }
      }
    }

    outputSpherocityMc(spherocityMc, spherocityReco);
  }
  PROCESS_SWITCH(EventShapeProducer, processSpherocityMC, "Spherocity from MC particles + reco on MC axis; output joinable to CFMcCollisions", false);

  void init(InitContext const&)
  {
    if (doprocessSpherocityData && doprocessSpherocityMC) {
      LOGF(fatal, "Enable only one of processSpherocityData or processSpherocityMC, not both.");
    }
    if (doprocessSpherocityEse && !doprocessSpherocityData) {
      LOGF(fatal, "processSpherocityEse requires processSpherocityData (reco spherocity and CFTracks).");
    }

    AxisSpec axisPt = {100, 0.0, 50.0};
    AxisSpec axisEta = {100, -1.2f, 1.2f};
    AxisSpec axisPhi = {100, 0.f, o2::constants::math::TwoPI};
    AxisSpec axisNtracks = {100, 0.0, 200.0};
    AxisSpec axisMult = {100, 0.0, 100.0};
    AxisSpec axisSpherocity = {100, 0.0, 1.0};
    AxisSpec axisPosZ = {100, -20., 20.};

    QAhistos.add("multS0", "", {HistType::kTH2F, {axisMult, axisSpherocity}});
    QAhistos.add("spherocity", "", {HistType::kTH1F, {axisSpherocity}});
    QAhistos.add("posZ_coll", "", {HistType::kTH1F, {axisPosZ}});
    if (doprocessSpherocityMC) {
      QAhistos.add("spherocityMc", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistos.add("spherocityRecoOnMc", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistos.add("posZ_mc", "", {HistType::kTH1F, {axisPosZ}});
      QAhistos.add("spherocityMc_vs_spherocityReco", "", {HistType::kTH2F, {axisSpherocity, axisSpherocity}});
      QAhistos.add("posZ_mc_vs_posZ_reco", "", {HistType::kTH2F, {axisPosZ, axisPosZ}});
    }

    if (doprocessSpherocityEse) {
      // ESE pT: 0.1 (0.2–1.0), 0.2 (1.0–4.0), then 4.0–4.5, 4.5–5, 5–6, 6–8 GeV/c
      AxisSpec axisPtEse{std::vector<double>{0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0, 3.2, 3.4, 3.6, 3.8, 4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4, 5.6, 5.8, 6.0, 6.2, 6.4, 6.6, 6.8, 7.0, 7.2, 7.4, 7.6, 7.8, 8.0, 8.5, 9.0, 10.0},
                         "p_{T} (GeV/c)"};
      AxisSpec axisEseMultBin = {kMultEseNBins, 0., static_cast<float>(kMultEseNBins)};
      QAhistosEse.add("ese_nevent_splow", "", {HistType::kTH1F, {axisEseMultBin}});
      QAhistosEse.add("ese_nevent_sphigh", "", {HistType::kTH1F, {axisEseMultBin}});
      QAhistosEse.add("ese_nevent_int", "", {HistType::kTH1F, {axisEseMultBin}});
      QAhistosEse.add("s0_dist_0", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_1", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_2", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_3", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_4", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_5", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_6", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("s0_dist_7", "", {HistType::kTH1F, {axisSpherocity}});
      QAhistosEse.add("pt_0_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_0_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_1_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_1_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_2_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_2_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_3_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_3_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_4_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_4_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_5_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_5_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_6_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_6_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_7_splow", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_7_sphigh", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_0_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_1_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_2_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_3_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_4_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_5_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_6_int", "", {HistType::kTH1F, {axisPtEse}});
      QAhistosEse.add("pt_7_int", "", {HistType::kTH1F, {axisPtEse}});
    }

    if (!outputQAHistos)
      return;
    if (!doprocessSpherocityMC) {
      QAhistos.add("pt", "", {HistType::kTH1F, {axisPt}});
      QAhistos.add("eta", "", {HistType::kTH1F, {axisEta}});
      QAhistos.add("phi", "", {HistType::kTH1F, {axisPhi}});
    }
    QAhistos.add("ntracks", "", {HistType::kTH1F, {axisNtracks}});
    QAhistos.add("mult", "", {HistType::kTH1F, {axisMult}});
    if (doprocessSpherocityMC) {
      QAhistos.add("pt", "", {HistType::kTH1F, {axisPt}});
      QAhistos.add("eta", "", {HistType::kTH1F, {axisEta}});
      QAhistos.add("phi", "", {HistType::kTH1F, {axisPhi}});
      QAhistos.add("ptMC", "", {HistType::kTH1F, {axisPt}});
      QAhistos.add("etaMC", "", {HistType::kTH1F, {axisEta}});
      QAhistos.add("phiMC", "", {HistType::kTH1F, {axisPhi}});
      QAhistos.add("ntracksMC", "", {HistType::kTH1F, {axisNtracks}});
    }

    QAhistosCent.add("ntracks_0_20", "", {HistType::kTH1F, {axisNtracks}});
    QAhistosCent.add("ntracks_20_40", "", {HistType::kTH1F, {axisNtracks}});
    QAhistosCent.add("ntracks_40_60", "", {HistType::kTH1F, {axisNtracks}});
    QAhistosCent.add("ntracks_60_100", "", {HistType::kTH1F, {axisNtracks}});

    QAhistosSpherocity.add("spherocity_0_20", "", {HistType::kTH1F, {axisSpherocity}});
    QAhistosSpherocity.add("spherocity_20_40", "", {HistType::kTH1F, {axisSpherocity}});
    QAhistosSpherocity.add("spherocity_40_60", "", {HistType::kTH1F, {axisSpherocity}});
    QAhistosSpherocity.add("spherocity_60_100", "", {HistType::kTH1F, {axisSpherocity}});
  }
};

WorkflowSpec defineDataProcessing(ConfigContext const& cfgc)
{
  return WorkflowSpec{
    adaptAnalysisTask<EventShapeProducer>(cfgc),
  };
}
