///////////////////////////////////////////////////////////////////
// TrackingVolumeArrayCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYTOOLS_TRACKINGVOLUMEARRAYCREATOR_H
#define ACTS_GEOMETRYTOOLS_TRACKINGVOLUMEARRAYCREATOR_H

// Core module
#include "ACTS/Tools/ITrackingVolumeArrayCreator.hpp"
#include "ACTS/Utilities/BinnedArray.hpp"
#include "ACTS/Utilities/Logger.hpp"
// STL
#include <algorithm>
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

    class Layer;
    class TrackingVolume;
    
    typedef std::pair< TrackingVolumePtr, Vector3D>   TrackingVolumeOrderPosition;
    
    /** @class TrackingVolumeArrayCreator

      The TrackingVolumeArrayCreator is a simple Tool that helps to construct
      binned arrays of TrackingVolumes for both, confinement in another volume 
      and navigation issues.
     
      @author Andreas.Salzburger@cern.ch   
     */

    class TrackingVolumeArrayCreator : public ITrackingVolumeArrayCreator {

      public:
      struct Config
      {
	std::shared_ptr<Logger>                 logger;                      //!< logging instance

	Config():
	  logger(getDefaultLogger("LayerArrayCreator",Logging::INFO))
	  {}
      };
      
        /** Constructor */
      TrackingVolumeArrayCreator(const Config& c):
	m_config(c)
      {}

      /** Destructor */
      virtual ~TrackingVolumeArrayCreator() = default;

      /** create a tracking volume array */
      std::shared_ptr<const TrackingVolumeArray> trackingVolumeArray(const TrackingVolumeVector& vols, BinningValue bVal) const;

      void setConfiguration(const Config& c) {m_config = c;}

       /** Get configuration method */
      Config getConfiguration() const {return m_config;}
      
      private:
        /** Configuration struct */
        Config m_config;
        const Logger& logger() const {return *m_config.logger;}
      
    };
}

#endif

