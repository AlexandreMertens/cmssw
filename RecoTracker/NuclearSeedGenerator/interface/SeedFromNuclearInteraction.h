#ifndef SeedFromNuclearInteraction_H
#define SeedFromNuclearInteraction_H

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"

#include "RecoTracker/NuclearSeedGenerator/interface/TangentHelix.h"

#include <boost/shared_ptr.hpp>

class FreeTrajectoryState;

class SeedFromNuclearInteraction {
private :
  typedef TrajectoryMeasurement                       TM;
  typedef TrajectoryStateOnSurface                    TSOS;
  typedef edm::OwnVector<TrackingRecHit>              recHitContainer;
  typedef TransientTrackingRecHit::ConstRecHitPointer ConstRecHitPointer;
  typedef std::vector<ConstRecHitPointer>             ConstRecHitContainer;

public :
  SeedFromNuclearInteraction(const edm::EventSetup& es, const edm::ParameterSet& iConfig);
 
  /// Fill all data members from 2 TM's where the first one is supposed to be at the interaction point
  void setMeasurements(const TSOS& tsosAtInteractionPoint, ConstRecHitPointer ihit, ConstRecHitPointer ohit);

  /// Fill all data members from 1 TSOS and 2 rec Hits and using the circle associated to the primary track as constraint
  void setMeasurements(TangentHelix& primHelix, const TSOS& inner_TSOS, ConstRecHitPointer ihit, ConstRecHitPointer ohit);

  PTrajectoryStateOnDet trajectoryState() const { return *pTraj; }

  FreeTrajectoryState* stateWithError() const;

  FreeTrajectoryState* stateWithError(TangentHelix& helix) const;

  PropagationDirection direction() const { return alongMomentum; }

  recHitContainer hits() const; 

  TrajectorySeed TrajSeed() const { return TrajectorySeed(trajectoryState(),hits(),direction()); }
 
  bool isValid() const { return isValid_; }

  const TSOS& updatedTSOS() const { return *updatedTSOS_; }

  const TSOS& initialTSOS() const { return *initialTSOS_; }

  GlobalPoint outerHitPosition() const {
          return pDD->idToDet(outerHitDetId())->surface().toGlobal(outerHit_->localPosition()); 
  }

  DetId outerHitDetId() const { return outerHit_->geographicalId(); }

  ConstRecHitPointer outerHit() const { return outerHit_; }

private :
  bool                                     isValid_;        /**< check if the seed is valid */

  ConstRecHitContainer                     theHits;         /**< all the hits to be used to update the */
                                                            /*   initial freeTS and to be fitted       */

  ConstRecHitPointer                       innerHit_;        /**< Pointer to the hit of the inner TM */
  ConstRecHitPointer                       outerHit_;        /**< Pointer to the outer hit */

  boost::shared_ptr<TSOS>                      updatedTSOS_;   /**< Final TSOS */

  boost::shared_ptr<TSOS>                      initialTSOS_;     /**< Initial TSOS used as input */
  
  boost::shared_ptr<FreeTrajectoryState>       freeTS_;          /**< Initial FreeTrajectoryState */

  boost::shared_ptr<PTrajectoryStateOnDet> pTraj;           /**< the final persistent TSOS */


  const Propagator*                        thePropagator;
  edm::ESHandle<TrackerGeometry>           pDD;

  bool construct();

  double rescaleDirectionFactor; /**< Rescale the direction error */
  double rescalePositionFactor;  /**< Rescale the position error */
  double rescaleCurvatureFactor; /**< Rescale the curvature error */
  double ptMin;                  /**< Minimum transverse momentum of the seed */
};
#endif
