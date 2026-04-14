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

#ifndef PWGCF_JCORRAN_DATAMODEL_CFSPHEROCITY_H_
#define PWGCF_JCORRAN_DATAMODEL_CFSPHEROCITY_H_

#include "PWGCF/DataModel/CorrelationsDerived.h"

#include "Framework/ASoA.h"

namespace o2::aod
{
namespace cfeventshape
{
DECLARE_SOA_COLUMN(Spherocity, spherocity, float);           //! Reco-track spherocity (per CFCollision)
DECLARE_SOA_COLUMN(SpherocityMc, spherocityMc, float);       //! MC-particle spherocity (per CFMcCollision)
DECLARE_SOA_COLUMN(SpherocityReco, spherocityReco, float);   //! Reco spherocity on MC axis (per CFMcCollision)
} // namespace cfeventshape

DECLARE_SOA_TABLE(CFSpherocity, "AOD", "CFSPHEROCOLL", cfeventshape::Spherocity);
DECLARE_SOA_TABLE(CFSpherocityMc, "AOD", "CFSPHEROMC", cfeventshape::SpherocityMc, cfeventshape::SpherocityReco);

using CFCollisionsWithSpherocity = soa::Join<CFCollisions, CFSpherocity>;
using CFCollisionWithSpherocity = CFCollisionsWithSpherocity::iterator;
using CFMcCollisionsWithSpherocity = soa::Join<CFMcCollisions, CFSpherocityMc>;
using CFMcCollisionWithSpherocity = CFMcCollisionsWithSpherocity::iterator;
} // namespace o2::aod

#endif // PWGCF_JCORRAN_DATAMODEL_CFSPHEROCITY_H_
