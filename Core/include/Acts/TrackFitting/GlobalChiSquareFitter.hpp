// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/SourceLink.hpp"
#include "Acts/EventData/TrackProxyConcept.hpp"
#include "Acts/EventData/Types.hpp"
#include "Acts/EventData/VectorMultiTrajectory.hpp"
#include "Acts/EventData/VectorTrackContainer.hpp"
#include "Acts/EventData/detail/CorrectedTransformationFreeToBound.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/Interactions.hpp"
#include "Acts/Propagator/DirectNavigator.hpp"
#include "Acts/Propagator/PropagatorOptions.hpp"
#include "Acts/Propagator/StandardAborters.hpp"
#include "Acts/Propagator/detail/PointwiseMaterialInteraction.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitterEnergyLossMode.hpp"
#include "Acts/TrackFitting/GlobalChiSquareFitterError.hpp"
#include "Acts/TrackFitting/detail/VoidFitterComponents.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/TrackHelpers.hpp"

#include <functional>
#include <limits>
#include <memory>
#include <type_traits>
#include <unordered_map>

namespace Acts::Experimental {

/// @addtogroup track_fitting
/// @{

namespace Gx2fConstants {
constexpr std::string_view gx2fnUpdateColumn = "Gx2fnUpdateColumn";

// Mask for the track states. We don't need Predicted and Filtered
constexpr TrackStatePropMask trackStateMask = TrackStatePropMask::Smoothed |
                                              TrackStatePropMask::Jacobian |
                                              TrackStatePropMask::Calibrated;

// A projector used for scattering. By using Jacobian * phiThetaProjector one
// gets only the derivatives for the variables phi and theta.
const Eigen::Matrix<double, eBoundSize, 2> phiThetaProjector = [] {
  Eigen::Matrix<double, eBoundSize, 2> m =
      Eigen::Matrix<double, eBoundSize, 2>::Zero();
  m(eBoundPhi, 0) = 1.0;
  m(eBoundTheta, 1) = 1.0;
  return m;
}();
}  // namespace Gx2fConstants

/// Extension struct which holds delegates to customise the GX2F behaviour
template <typename traj_t>
struct Gx2FitterExtensions {
  /// Type alias for mutable track state proxy from multi-trajectory
  using TrackStateProxy = typename MultiTrajectory<traj_t>::TrackStateProxy;
  /// Type alias for const track state proxy from multi-trajectory
  using ConstTrackStateProxy =
      typename MultiTrajectory<traj_t>::ConstTrackStateProxy;
  /// Type alias for track parameters from track state proxy
  using Parameters = typename TrackStateProxy::Parameters;

  /// Type alias for calibrator delegate to process measurements
  using Calibrator =
      Delegate<void(const GeometryContext&, const CalibrationContext&,
                    const SourceLink&, TrackStateProxy)>;

  /// Type alias for updater delegate to incorporate measurements into track
  /// parameters
  using Updater = Delegate<Result<void>(const GeometryContext&, TrackStateProxy,
                                        const Logger&)>;

  /// Type alias for outlier finder delegate to identify measurement outliers
  using OutlierFinder = Delegate<bool(ConstTrackStateProxy)>;

  /// The Calibrator is a dedicated calibration algorithm that allows
  /// to calibrate measurements using track information, this could be
  /// e.g. sagging for wires, module deformations, etc.
  Calibrator calibrator;

  /// The updater incorporates measurement information into the track parameters
  Updater updater;

  /// Determines whether a measurement is supposed to be considered as an
  /// outlier
  OutlierFinder outlierFinder;

  /// Retrieves the associated surface from a source link
  SourceLinkSurfaceAccessor surfaceAccessor;

  /// Default constructor which connects the default void components
  Gx2FitterExtensions() {
    calibrator.template connect<&Acts::detail::voidFitterCalibrator<traj_t>>();
    updater.template connect<&Acts::detail::voidFitterUpdater<traj_t>>();
    outlierFinder.template connect<&Acts::detail::voidOutlierFinder<traj_t>>();
    surfaceAccessor.connect<&Acts::detail::voidSurfaceAccessor>();
  }
};

/// Combined options for the Global-Chi-Square fitter.
///
/// @tparam traj_t The trajectory type
template <typename traj_t>
struct Gx2FitterOptions {
  /// PropagatorOptions with context.
  ///
  /// @param gctx The geometry context for this fit
  /// @param mctx The magnetic context for this fit
  /// @param cctx The calibration context for this fit
  /// @param extensions_ The KF extensions
  /// @param pOptions The plain propagator options
  /// @param rSurface The reference surface for the fit to be expressed at
  /// @param mScattering Whether to include multiple scattering
  /// @param eLoss Whether to include energy loss
  /// @param freeToBoundCorrection_ Correction for non-linearity effect during transform from free to bound
  /// @param nUpdateMax_ Max number of iterations for updating the parameters
  /// @param relChi2changeCutOff_ Check for convergence (abort condition). Set to 0 to skip.
  /// @param eLossMode Whether to use the mean or the most probable energy loss
  /// @param nMaterialUpdateMax_ Number of iterations of the material fit
  Gx2FitterOptions(const GeometryContext& gctx,
                   const MagneticFieldContext& mctx,
                   std::reference_wrapper<const CalibrationContext> cctx,
                   Gx2FitterExtensions<traj_t> extensions_,
                   const PropagatorPlainOptions& pOptions,
                   const Surface* rSurface = nullptr, bool mScattering = false,
                   bool eLoss = false,
                   const FreeToBoundCorrection& freeToBoundCorrection_ =
                       FreeToBoundCorrection(false),
                   const std::size_t nUpdateMax_ = 5,
                   double relChi2changeCutOff_ = 1e-5,
                   Gx2fEnergyLossMode eLossMode = Gx2fEnergyLossMode::Mean,
                   const std::size_t nMaterialUpdateMax_ = 1)
      : geoContext(gctx),
        magFieldContext(mctx),
        calibrationContext(cctx),
        extensions(extensions_),
        propagatorPlainOptions(pOptions),
        referenceSurface(rSurface),
        multipleScattering(mScattering),
        energyLoss(eLoss),
        energyLossMode(eLossMode),
        freeToBoundCorrection(freeToBoundCorrection_),
        nUpdateMax(nUpdateMax_),
        nMaterialUpdateMax(nMaterialUpdateMax_),
        relChi2changeCutOff(relChi2changeCutOff_) {}

  /// Contexts are required and the options must not be default-constructible.
  Gx2FitterOptions() = delete;

  /// Context object for the geometry
  std::reference_wrapper<const GeometryContext> geoContext;
  /// Context object for the magnetic field
  std::reference_wrapper<const MagneticFieldContext> magFieldContext;
  /// context object for the calibration
  std::reference_wrapper<const CalibrationContext> calibrationContext;

  /// Extensions for calibration and outlier finding
  Gx2FitterExtensions<traj_t> extensions;

  /// The trivial propagator options
  PropagatorPlainOptions propagatorPlainOptions;

  /// The reference Surface
  const Surface* referenceSurface = nullptr;

  /// Whether to consider multiple scattering
  bool multipleScattering = false;

  /// Whether to consider energy loss
  bool energyLoss = false;

  /// Whether to use the mean or the most probable energy loss.
  /// Defaults to the mean, matching the KF and CKF.
  Gx2fEnergyLossMode energyLossMode = Gx2fEnergyLossMode::Mean;

  /// Whether to include non-linear correction during global to local
  /// transformation
  FreeToBoundCorrection freeToBoundCorrection;

  /// Max number of iterations during the fit (abort condition)
  std::size_t nUpdateMax = 5;

  /// Number of iterations of the material fit. The material parameters are
  /// fitted after the main loop converged. One iteration reproduces the
  /// historical behaviour.
  std::size_t nMaterialUpdateMax = 1;

  /// Check for convergence (abort condition). Set to 0 to skip.
  double relChi2changeCutOff = 1e-7;
};

/// Result container for a global chi-square fit.
template <typename traj_t>
struct Gx2FitterResult {
  /// Fitted states that the actor has handled.
  traj_t* fittedStates{nullptr};

  /// This is the index of the 'tip' of the track stored in multitrajectory.
  /// This corresponds to the last measurement state in the multitrajectory.
  /// Since this KF only stores one trajectory, it is unambiguous.
  /// Acts::TrackTraits::kInvalid is the start of a trajectory.
  std::size_t lastMeasurementIndex = Acts::kTrackIndexInvalid;

  /// This is the index of the 'tip' of the states stored in multitrajectory.
  /// This corresponds to the last state in the multitrajectory.
  /// Since this KF only stores one trajectory, it is unambiguous.
  /// Acts::TrackTraits::kInvalid is the start of a trajectory.
  std::size_t lastTrackIndex = Acts::kTrackIndexInvalid;

  /// The optional Parameters at the provided surface
  std::optional<BoundTrackParameters> fittedParameters;

  /// Counter for states with non-outlier measurements
  std::size_t measurementStates = 0;

  /// Counter for measurements holes
  /// A hole correspond to a surface with an associated detector element with no
  /// associated measurement. Holes are only taken into account if they are
  /// between the first and last measurements.
  std::size_t measurementHoles = 0;

  /// Counter for handled states
  std::size_t processedStates = 0;

  /// Counter for handled measurements
  std::size_t processedMeasurements = 0;

  /// Indicator if track fitting has been done
  bool finished = false;

  /// Measurement surfaces without hits
  std::vector<const Surface*> missedActiveSurfaces;

  /// Measurement surfaces handled in both forward and
  /// backward filtering
  std::vector<const Surface*> passedAgainSurfaces;

  /// Count how many surfaces have been hit
  std::size_t surfaceCount = 0;
};

/// @brief A container to store the material properties of a surface
///
/// This struct holds the scattering angles, the inverse covariance of the
/// material, and a validity flag indicating whether the material is valid for
/// the scattering process.
struct Gx2fMaterialProperties {
 public:
  /// @brief Constructor to initialize the material properties.
  ///
  /// @param scatteringAngles_ The vector of scattering angles.
  /// @param invCovarianceMaterial_ The inverse covariance of the material.
  /// @param materialIsValid_ A boolean flag indicating whether the material is valid.
  Gx2fMaterialProperties(const BoundVector& scatteringAngles_,
                         const double invCovarianceMaterial_,
                         const bool materialIsValid_)
      : m_scatteringAngles(scatteringAngles_),
        m_invCovarianceMaterial(invCovarianceMaterial_),
        m_materialIsValid(materialIsValid_) {}

  /// @brief Accessor for the scattering angles (const version)
  /// @return Const reference to the vector of scattering angles
  const BoundVector& scatteringAngles() const { return m_scatteringAngles; }

  /// @brief Accessor for the scattering angles (mutable version)
  /// @return Mutable reference to the vector of scattering angles for modification
  BoundVector& scatteringAngles() { return m_scatteringAngles; }

  /// @brief Accessor for the inverse covariance of the material
  /// @return Inverse covariance value computed from material properties (e.g., Highland formula)
  double invCovarianceMaterial() const { return m_invCovarianceMaterial; }

  /// @brief Accessor for the material validity flag
  /// @return True if material is valid for scattering calculations, false for vacuum or zero thickness
  bool materialIsValid() const { return m_materialIsValid; }

  /// @brief Invalidate this material surface, so it is ignored everywhere
  void invalidateMaterial() { m_materialIsValid = false; }

  /// @brief Accessor for the deterministic q/p change from energy loss (const version)
  /// @return The expected q/p change, evaluated at the local q/p. Not a free parameter.
  double expectedQOverPOffset() const { return m_expectedQOverPOffset; }

  /// @brief Accessor for the deterministic q/p change from energy loss (mutable version)
  /// @return Mutable reference for refreshing the expectation during propagation
  double& expectedQOverPOffset() { return m_expectedQOverPOffset; }

  /// @brief Accessor for the fitted deviation from the expected q/p change (const version)
  /// @return The free parameter of the energy loss fit. Starts at 0.
  double qOverPOffset() const { return m_qOverPOffset; }

  /// @brief Accessor for the fitted deviation from the expected q/p change (mutable version)
  /// @return Mutable reference to the free parameter for the fit update
  double& qOverPOffset() { return m_qOverPOffset; }

  /// @brief The total q/p change applied at this surface
  /// @return Sum of the deterministic expectation and the fitted deviation
  double totalQOverPOffset() const {
    return m_expectedQOverPOffset + m_qOverPOffset;
  }

  /// @brief Accessor for the inverse covariance of the q/p straggling (const version)
  /// @return 1/sigma^2 of the energy loss straggling, converted to q/p units
  double invCovarianceQOverP() const { return m_invCovarianceQOverP; }

  /// @brief Accessor for the inverse covariance of the q/p straggling (mutable version)
  /// @return Mutable reference for refreshing the straggling during propagation
  double& invCovarianceQOverP() { return m_invCovarianceQOverP; }

 private:
  /// Vector of scattering angles. The vector is usually all zeros except for
  /// eBoundPhi and eBoundTheta.
  BoundVector m_scatteringAngles;

  /// Inverse covariance of the material. Compute with e.g. the Highland
  /// formula.
  double m_invCovarianceMaterial;

  /// Flag indicating whether the material is valid. Commonly vacuum and zero
  /// thickness material will be ignored.
  bool m_materialIsValid;

  /// Deterministic q/p change from the mean/mode energy loss at this surface.
  /// Refreshed from the local q/p while the fit iterates and frozen for the
  /// final propagation.
  double m_expectedQOverPOffset = 0.;

  /// Fitted deviation from the expected q/p change. This is the free parameter.
  double m_qOverPOffset = 0.;

  /// Inverse covariance of the q/p straggling. Compute with e.g. the Landau
  /// sigma converted to q/p units.
  double m_invCovarianceQOverP = 0.;
};

/// @brief Describes which per-material-surface parameters are fitted, and where
///        they live in the extended parameter vector
///
/// The extended parameter vector is laid out as
///   [0, eBoundSize)                      the bound track parameters
///   [eBoundSize + stride() * k, ...)     the parameters of material surface k
/// where each fitted material surface contributes stride() parameters:
///   phi, theta      if fitScattering()
///   delta(q/p)      if fitEnergyLoss()
class Gx2fParameterLayout {
 public:
  /// @brief Default constructor for a system without any material parameters
  Gx2fParameterLayout() = default;

  /// @brief Constructor
  ///
  /// @param fitScattering_ Fit two scattering angles per material surface
  /// @param fitEnergyLoss_ Fit one q/p deviation per material surface
  /// @param nMaterialSurfaces_ Number of fitted material surfaces
  Gx2fParameterLayout(bool fitScattering_, bool fitEnergyLoss_,
                      std::size_t nMaterialSurfaces_)
      : m_fitScattering{fitScattering_},
        m_fitEnergyLoss{fitEnergyLoss_},
        m_nMaterialSurfaces{
            fitScattering_ || fitEnergyLoss_ ? nMaterialSurfaces_ : 0u} {}

  /// @brief Whether the scattering angles are fitted
  /// @return True if each material surface contributes two angle parameters
  bool fitScattering() const { return m_fitScattering; }

  /// @brief Whether the energy loss deviation is fitted
  /// @return True if each material surface contributes one q/p parameter
  bool fitEnergyLoss() const { return m_fitEnergyLoss; }

  /// @brief Whether any material parameters are fitted at all
  /// @return True if material surfaces contribute to the extended system
  bool fitMaterial() const { return m_fitScattering || m_fitEnergyLoss; }

  /// @brief Number of extended parameters contributed by one material surface
  /// @return The stride between two consecutive material surface blocks
  std::size_t stride() const {
    return (m_fitScattering ? 2u : 0u) + (m_fitEnergyLoss ? 1u : 0u);
  }

  /// @brief Accessor for the number of fitted material surfaces
  /// @return Number of material surfaces contributing to the extended system
  std::size_t nMaterialSurfaces() const { return m_nMaterialSurfaces; }

  /// @brief Total dimension of the extended system
  /// @return eBoundSize plus the parameters of all fitted material surfaces
  std::size_t nDims() const {
    return eBoundSize + stride() * m_nMaterialSurfaces;
  }

  /// @brief Index of the phi column of the k-th fitted material surface
  ///
  /// The theta column is at scatteringOffset(k) + 1.
  ///
  /// @param k Index of the fitted material surface
  /// @return Position of the phi parameter in the extended system
  std::size_t scatteringOffset(std::size_t k) const {
    assert(m_fitScattering && k < m_nMaterialSurfaces);
    return eBoundSize + stride() * k;
  }

  /// @brief Index of the delta(q/p) column of the k-th fitted material surface
  ///
  /// @param k Index of the fitted material surface
  /// @return Position of the energy loss parameter in the extended system
  std::size_t energyLossOffset(std::size_t k) const {
    assert(m_fitEnergyLoss && k < m_nMaterialSurfaces);
    return eBoundSize + stride() * k + (m_fitScattering ? 2u : 0u);
  }

 private:
  /// Whether two scattering angles per material surface are fitted
  bool m_fitScattering = false;

  /// Whether one q/p deviation per material surface is fitted
  bool m_fitEnergyLoss = false;

  /// Number of fitted material surfaces
  std::size_t m_nMaterialSurfaces = 0;
};

/// @brief A container to manage all properties of a gx2f system
///
/// This struct manages the mathematical infrastructure for the gx2f. It
/// initializes and maintains the extended aMatrix and extended bVector.
struct Gx2fSystem {
 public:
  /// @brief Constructor to initialize matrices and vectors to zero based on the parameter layout.
  ///
  /// @param layout Describes which material parameters are fitted and where they live.
  explicit Gx2fSystem(const Gx2fParameterLayout& layout)
      : m_layout{layout},
        m_nDims{layout.nDims()},
        m_aMatrix{Eigen::MatrixXd::Zero(m_nDims, m_nDims)},
        m_bVector{Eigen::VectorXd::Zero(m_nDims)} {}

  /// @brief Accessor for the parameter layout of the extended system
  /// @return The layout describing which material parameters are fitted
  const Gx2fParameterLayout& layout() const { return m_layout; }

  /// @brief Accessor for the number of dimensions of the extended system
  /// @return Number of dimensions for the aMatrix and bVector (bound parameters + material parameters)
  std::size_t nDims() const { return m_nDims; }

  /// @brief Accessor for the accumulated chi-squared value (const version)
  /// @return Current sum of chi-squared contributions from measurements and material
  double chi2() const { return m_chi2; }

  /// @brief Accessor for the accumulated chi-squared value (mutable version)
  /// @return Mutable reference to chi-squared sum for modification during fitting
  double& chi2() { return m_chi2; }

  /// @brief Accessor for the extended system matrix (const version)
  /// @return Const reference to the aMatrix containing measurement and material contributions
  const Eigen::MatrixXd& aMatrix() const { return m_aMatrix; }

  /// @brief Accessor for the extended system matrix (mutable version)
  /// @return Mutable reference to the aMatrix for adding measurement and material contributions
  Eigen::MatrixXd& aMatrix() { return m_aMatrix; }

  /// @brief Accessor for the extended system vector (const version)
  /// @return Const reference to the bVector containing measurement and material contributions
  const Eigen::VectorXd& bVector() const { return m_bVector; }

  /// @brief Accessor for the extended system vector (mutable version)
  /// @return Mutable reference to the bVector for adding measurement and material contributions
  Eigen::VectorXd& bVector() { return m_bVector; }

  /// @brief Accessor for the number of degrees of freedom (const version)
  /// @return Current number of degrees of freedom from processed measurements
  std::size_t ndf() const { return m_ndf; }

  /// @brief Accessor for the number of degrees of freedom (mutable version)
  /// @return Mutable reference to NDF counter for incrementing during measurement processing
  std::size_t& ndf() { return m_ndf; }

  /// @brief Determines the minimum number of degrees of freedom required for the fit
  ///
  /// Automatically deduces the required NDF based on the system configuration.
  /// We have only 3 cases, because we always have l0, l1, phi, theta:
  /// - 4: no magnetic field -> q/p is empty
  /// - 5: no time measurement -> time is not fittable
  /// - 6: full fit with all parameters
  ///
  /// @return Required NDF based on which parameters can be fitted
  std::size_t findRequiredNdf() {
    std::size_t ndfSystem = 0;
    if (m_aMatrix(4, 4) == 0) {
      ndfSystem = 4;
    } else if (m_aMatrix(5, 5) == 0) {
      ndfSystem = 5;
    } else {
      ndfSystem = 6;
    }

    return ndfSystem;
  }

  /// @brief Checks if the system has sufficient degrees of freedom for fitting
  /// @return True if NDF exceeds the minimum required for the parameter configuration
  bool isWellDefined() { return m_ndf > findRequiredNdf(); }

 private:
  /// Layout of the material parameters within the extended system
  Gx2fParameterLayout m_layout;

  /// Number of dimensions of the (extended) system
  std::size_t m_nDims;

  /// Sum of chi-squared values.
  double m_chi2 = 0.;

  /// Extended matrix for accumulation.
  Eigen::MatrixXd m_aMatrix;

  /// Extended vector for accumulation.
  Eigen::VectorXd m_bVector;

  /// Number of degrees of freedom of the system
  std::size_t m_ndf = 0u;
};

/// @brief Adds a measurement to the GX2F equation system in a modular backend function.
///
/// This function processes measurement data and integrates it into the GX2F
/// system.
///
/// @param extendedSystem All parameters of the current equation system to update.
/// @param jacobianFromStart The Jacobian matrix from the start to the current state.
/// @param covarianceMeasurement The covariance matrix of the measurement.
/// @param predicted The predicted state vector based on the track state.
/// @param measurement The measurement vector.
/// @param projector The projection matrix.
/// @param logger A logger instance.
///
/// @note The dynamic Eigen matrices are suboptimal. We could think of
/// templating again in the future on kMeasDims. We currently use dynamic
/// matrices to reduce the memory during compile time.
void addMeasurementToGx2fSumsBackend(
    Gx2fSystem& extendedSystem,
    const std::vector<BoundMatrix>& jacobianFromStart,
    const Eigen::MatrixXd& covarianceMeasurement, const BoundVector& predicted,
    const Eigen::VectorXd& measurement, const Eigen::MatrixXd& projector,
    const Logger& logger);

/// @brief Process measurements and fill the aMatrix and bVector
///
/// The function processes each measurement for the GX2F Actor fitting process.
/// It extracts the information from the track state and adds it to aMatrix,
/// bVector, and chi2sum.
///
/// @tparam kMeasDim Number of dimensions of the measurement
/// @tparam track_state_t The type of the track state
///
/// @param extendedSystem All parameters of the current equation system to update
/// @param jacobianFromStart The Jacobian matrix from start to the current state
/// @param trackState The track state to analyse
/// @param logger A logger instance
template <std::size_t kMeasDim, typename track_state_t>
void addMeasurementToGx2fSums(Gx2fSystem& extendedSystem,
                              const std::vector<BoundMatrix>& jacobianFromStart,
                              const track_state_t& trackState,
                              const Logger& logger) {
  const SquareMatrix<kMeasDim> covarianceMeasurement =
      trackState.template calibratedCovariance<kMeasDim>();

  const BoundVector predicted = trackState.smoothed();

  const Vector<kMeasDim> measurement =
      trackState.template calibrated<kMeasDim>();

  const Matrix<kMeasDim, eBoundSize> projector =
      trackState.template projectorSubspaceHelper<kMeasDim>().projector();

  addMeasurementToGx2fSumsBackend(extendedSystem, jacobianFromStart,
                                  covarianceMeasurement, predicted, measurement,
                                  projector, logger);
}

/// @brief Process material and fill the aMatrix and bVector
///
/// The function processes each material for the GX2F Actor fitting process.
/// It extracts the information from the track state and adds it to aMatrix,
/// bVector, and chi2sum.
///
/// @tparam track_state_t The type of the track state
///
/// @param extendedSystem All parameters of the current equation system
/// @param nMaterialsHandled How many materials we already handled. Used for the offset.
/// @param materialMap The material map, containing all material properties
/// @param trackState The track state to analyse
/// @param logger A logger instance
template <typename track_state_t>
void addMaterialToGx2fSums(
    Gx2fSystem& extendedSystem, const std::size_t nMaterialsHandled,
    const std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>&
        materialMap,
    const track_state_t& trackState, const Logger& logger) {
  // Get and store geoId for the current material surface
  const GeometryIdentifier geoId = trackState.referenceSurface().geometryId();
  const auto materialMapId = materialMap.find(geoId);
  if (materialMapId == materialMap.end()) {
    ACTS_ERROR("No material properties found for material surface " << geoId);
    throw std::runtime_error(
        "No material properties found for material surface.");
  }

  const Gx2fParameterLayout& layout = extendedSystem.layout();

  if (layout.fitEnergyLoss()) {
    // The position, where we need to insert the values in aMatrix and bVector
    const std::size_t elPosition = layout.energyLossOffset(nMaterialsHandled);

    // The expected q/p change is already part of the model, so only the fitted
    // deviation from it enters the penalty. Its expectation is 0, exactly like
    // for the scattering angles.
    const double deltaQOverP = materialMapId->second.qOverPOffset();
    const double invCovQOverP = materialMapId->second.invCovarianceQOverP();

    extendedSystem.aMatrix()(elPosition, elPosition) += invCovQOverP;
    extendedSystem.bVector()(elPosition, 0) -= invCovQOverP * deltaQOverP;
    extendedSystem.chi2() += invCovQOverP * deltaQOverP * deltaQOverP;

    ACTS_VERBOSE("Energy loss contribution in addMaterialToGx2fSums:\n"
                 << "    invCov:               " << invCovQOverP << "\n"
                 << "    elPosition:           " << elPosition << "\n"
                 << "    delta(q/p):           " << deltaQOverP << "\n"
                 << "    aMatrix contribution: " << invCovQOverP << "\n"
                 << "    bVector contribution: " << invCovQOverP * deltaQOverP
                 << "\n"
                 << "    chi2sum contribution: "
                 << invCovQOverP * deltaQOverP * deltaQOverP << "\n");
  }

  if (!layout.fitScattering()) {
    return;
  }

  // The position, where we need to insert the values in aMatrix and bVector
  const std::size_t deltaPosition = layout.scatteringOffset(nMaterialsHandled);

  const BoundVector& scatteringAngles =
      materialMapId->second.scatteringAngles();

  const double invCovTheta = materialMapId->second.invCovarianceMaterial();
  const double sinThetaLoc = std::sin(trackState.smoothed()[eBoundTheta]);
  const double invCovPhi = invCovTheta * sinThetaLoc * sinThetaLoc;

  // Phi contribution
  extendedSystem.aMatrix()(deltaPosition, deltaPosition) += invCovPhi;
  extendedSystem.bVector()(deltaPosition, 0) -=
      invCovPhi * scatteringAngles[eBoundPhi];
  extendedSystem.chi2() +=
      invCovPhi * scatteringAngles[eBoundPhi] * scatteringAngles[eBoundPhi];

  // Theta Contribution
  extendedSystem.aMatrix()(deltaPosition + 1, deltaPosition + 1) += invCovTheta;
  extendedSystem.bVector()(deltaPosition + 1, 0) -=
      invCovTheta * scatteringAngles[eBoundTheta];
  extendedSystem.chi2() += invCovTheta * scatteringAngles[eBoundTheta] *
                           scatteringAngles[eBoundTheta];

  ACTS_VERBOSE(
      "Contributions in addMaterialToGx2fSums:\n"
      << "    invCov:        " << invCovPhi << "\n"
      << "    sinThetaLoc:   " << sinThetaLoc << "\n"
      << "    deltaPosition: " << deltaPosition << "\n"
      << "    Phi:\n"
      << "        scattering angle:     " << scatteringAngles[eBoundPhi] << "\n"
      << "        aMatrix contribution: " << invCovPhi << "\n"
      << "        bVector contribution: "
      << invCovPhi * scatteringAngles[eBoundPhi] << "\n"
      << "        chi2sum contribution: "
      << invCovPhi * scatteringAngles[eBoundPhi] * scatteringAngles[eBoundPhi]
      << "\n"
      << "    Theta:\n"
      << "        scattering angle:     " << scatteringAngles[eBoundTheta]
      << "\n"
      << "        aMatrix contribution: " << invCovTheta << "\n"
      << "        bVector contribution: "
      << invCovTheta * scatteringAngles[eBoundTheta] << "\n"
      << "        chi2sum contribution: "
      << invCovTheta * scatteringAngles[eBoundTheta] *
             scatteringAngles[eBoundTheta]
      << "\n");

  return;
}

/// @brief Fill the GX2F system with data from a track
///
/// This function processes a track proxy and updates the aMatrix, bVector, and
/// chi2 values for the GX2F fitting system.
///
/// @note @p handleMaterial must match what the actor did during the
/// propagation that produced @p track, and is therefore independent of whether
/// material parameters are actually fitted (which is described by the layout of
/// @p extendedSystem). The actor resets the transport Jacobian on every surface
/// where it handles material, so such states must be chained here even when
/// they contribute no free parameters. If the two ever disagree, the Jacobians
/// of all downstream measurements are silently wrong.
///
/// @tparam track_proxy_t The type of the track proxy
///
/// @param track A constant track proxy to inspect
/// @param extendedSystem All parameters of the current equation system
/// @param handleMaterial Whether the actor handled material during the propagation
/// @param materialMap Map of geometry identifiers to material properties,
///        containing scattering angles and validation status
/// @param geoIdVector A vector to store geometry identifiers for tracking processed elements
/// @param logger A logger instance
template <TrackProxyConcept track_proxy_t>
void fillGx2fSystem(
    const track_proxy_t track, Gx2fSystem& extendedSystem,
    const bool handleMaterial,
    const std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>&
        materialMap,
    std::vector<GeometryIdentifier>& geoIdVector, const Logger& logger) {
  assert((handleMaterial || !extendedSystem.layout().fitMaterial()) &&
         "Cannot fit material parameters on surfaces the actor did not "
         "handle.");

  std::vector<BoundMatrix> jacobianFromStart;
  jacobianFromStart.emplace_back(BoundMatrix::Identity());

  for (const auto& trackState : track.trackStates()) {
    // Get and store geoId for the current surface
    const GeometryIdentifier geoId = trackState.referenceSurface().geometryId();
    ACTS_DEBUG("Start to investigate trackState on surface " << geoId);
    const auto typeFlags = trackState.typeFlags();
    const bool stateHasMeasurement = typeFlags.hasMeasurement();
    const bool stateHasMaterial = typeFlags.hasMaterial();

    // First we figure out, if we would need to look into material
    // surfaces at all. Later, we also check, if the material slab is
    // valid, otherwise we modify this flag to ignore the material
    // completely.
    bool doMaterial = handleMaterial && stateHasMaterial;
    if (doMaterial) {
      const auto materialMapId = materialMap.find(geoId);
      assert(materialMapId != materialMap.end() &&
             "No material properties found for material surface.");
      doMaterial = doMaterial && materialMapId->second.materialIsValid();
    }

    // We only consider states with a measurement (and/or material)
    if (!stateHasMeasurement && !doMaterial) {
      ACTS_DEBUG("    Skip state.");
      continue;
    }

    // update all Jacobians from start
    for (auto& jac : jacobianFromStart) {
      jac = trackState.jacobian() * jac;
    }

    // Handle measurement
    if (stateHasMeasurement) {
      ACTS_DEBUG("    Handle measurement.");

      const auto measDim = trackState.calibratedSize();

      if (measDim < 1 || 6 < measDim) {
        ACTS_ERROR("Can not process state with measurement with "
                   << measDim << " dimensions.");
        throw std::domain_error(
            "Found measurement with less than 1 or more than 6 dimension(s).");
      }

      extendedSystem.ndf() += measDim;

      visit_measurement(measDim, [&](auto N) {
        addMeasurementToGx2fSums<N>(extendedSystem, jacobianFromStart,
                                    trackState, logger);
      });
    }

    // Handle material. Only surfaces that contribute free parameters open a
    // new Jacobian and a new block of columns.
    if (doMaterial && extendedSystem.layout().fitMaterial()) {
      ACTS_DEBUG("    Handle material");
      // Add for this material a new Jacobian, starting from this surface.
      jacobianFromStart.emplace_back(BoundMatrix::Identity());

      // Add the material contribution to the system
      addMaterialToGx2fSums(extendedSystem, geoIdVector.size(), materialMap,
                            trackState, logger);

      geoIdVector.emplace_back(geoId);
    }
  }
}

/// @brief Count the valid material states in a track for material calculations.
///
/// This function counts the valid material surfaces encountered in a track
/// by examining each track state. The count is based on the presence of
/// material flags and the availability of material information for each
/// surface.
///
/// @tparam track_proxy_t The type of the track proxy
///
/// @param track A constant track proxy to inspect
/// @param materialMap Map of geometry identifiers to material properties,
///        containing scattering angles and validation status
/// @param logger A logger instance
/// @return Number of valid material states in the track
template <TrackProxyConcept track_proxy_t>
std::size_t countMaterialStates(
    const track_proxy_t track,
    const std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>&
        materialMap,
    const Logger& logger) {
  std::size_t nMaterialSurfaces = 0;
  ACTS_DEBUG("Count the valid material surfaces.");
  for (const auto& trackState : track.trackStates()) {
    const auto typeFlags = trackState.typeFlags();
    const bool stateHasMaterial = typeFlags.hasMaterial();

    if (!stateHasMaterial) {
      continue;
    }

    // Get and store geoId for the current material surface
    const GeometryIdentifier geoId = trackState.referenceSurface().geometryId();

    const auto materialMapId = materialMap.find(geoId);
    assert(materialMapId != materialMap.end() &&
           "No material properties found for material surface.");
    if (!materialMapId->second.materialIsValid()) {
      continue;
    }

    nMaterialSurfaces++;
  }

  return nMaterialSurfaces;
}

/// @brief Solve the gx2f system to get the delta parameters for the update
///
/// This function computes the delta parameters for the GX2F Actor fitting
/// process by solving the linear equation system [a] * delta = b. It uses the
/// column-pivoting Householder QR decomposition for numerical stability.
///
/// @param extendedSystem All parameters of the current equation system
/// @return Delta parameters for the GX2F update
Eigen::VectorXd computeGx2fDeltaParams(const Gx2fSystem& extendedSystem);

/// @brief Deterministic q/p change caused by the energy loss in a material slab
///
/// Mirrors the energy loss update of
/// @ref Acts::detail::performMaterialInteraction, so that the GX2F applies the
/// same correction as the KF.
///
/// The @c computeEnergyLoss* functions return the positive magnitude of the
/// loss. Following the same convention, the energy after the slab is
/// `E - eLoss * direction`. For forward propagation the energy decreases,
/// therefore the momentum decreases and the magnitude of q/p increases. The
/// returned offset thus has the same sign as @p qOverP.
///
/// @param slab Material slab, already corrected for the path length
/// @param particleHypothesis The particle hypothesis
/// @param qOverP Local q/p in front of the slab
/// @param direction The propagation direction
/// @param mode Whether to use the mean or the most probable energy loss
/// @return The q/p offset to be added to the local q/p
double computeGx2fQOverPOffset(const MaterialSlab& slab,
                               const ParticleHypothesis& particleHypothesis,
                               double qOverP, Direction direction,
                               Gx2fEnergyLossMode mode);

/// @brief Update parameters (and scattering angles if applicable)
///
/// @param params Parameters to be updated
/// @param deltaParamsExtended Delta parameters for bound parameter and material parameters
/// @param layout Layout of the material parameters in the extended system
/// @param materialMap Map of geometry identifiers to material properties,
///        containing all scattering angles and covariances
/// @param geoIdVector Vector of geometry identifiers corresponding to material surfaces
void updateGx2fParams(
    BoundTrackParameters& params, const Eigen::VectorXd& deltaParamsExtended,
    const Gx2fParameterLayout& layout,
    std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>& materialMap,
    const std::vector<GeometryIdentifier>& geoIdVector);

/// @brief Calculate and update the covariance of the fitted parameters
///
/// This function calculates the covariance of the fitted parameters using
/// cov = inv([a])
/// It then updates the first square block of size ndfSystem. This ensures,
/// that we only update the covariance for fitted parameters. (In case of
/// no qop/time fit)
///
/// @param fullCovariancePredicted The covariance matrix to update
/// @param extendedSystem All parameters of the current equation system
void updateGx2fCovarianceParams(BoundMatrix& fullCovariancePredicted,
                                Gx2fSystem& extendedSystem);

/// Global Chi Square fitter (GX2F) implementation.
///
/// @tparam propagator_t Type of the propagation class
///
/// TODO Write description
template <typename propagator_t, typename traj_t>
class Gx2Fitter {
  /// The navigator type
  using Gx2fNavigator = typename propagator_t::Navigator;

  /// The navigator has DirectNavigator type or not
  static constexpr bool isDirectNavigator =
      std::is_same_v<Gx2fNavigator, DirectNavigator>;

  static constexpr auto kInvalid = kTrackIndexInvalid;

 public:
  /// @brief Constructor for the Global Chi-Square Fitter
  ///
  /// Initializes the fitter with a propagator and optional logger.
  /// The fitter uses iterative fitting with a linear equation system
  /// to minimize chi-squared including multiple scattering effects.
  ///
  /// @param pPropagator The propagator instance for track propagation
  /// @param _logger Logger instance for debugging output (optional)
  explicit Gx2Fitter(const propagator_t& pPropagator,
                     std::unique_ptr<const Logger> _logger =
                         getDefaultLogger("Gx2Fitter", Logging::INFO))
      : m_propagator(pPropagator),
        m_logger{std::move(_logger)},
        m_actorLogger{m_logger->cloneWithSuffix("Actor")},
        m_addToSumLogger{m_logger->cloneWithSuffix("AddToSum")} {}

 private:
  /// The propagator for the transport and material update
  propagator_t m_propagator;

  /// The logger instance
  std::unique_ptr<const Logger> m_logger;
  std::unique_ptr<const Logger> m_actorLogger;
  std::unique_ptr<const Logger> m_addToSumLogger;

  const Logger& logger() const { return *m_logger; }

  /// @brief Propagator Actor plugin for the GX2F
  ///
  /// @tparam parameters_t The type of parameters used for "local" parameters.
  /// @tparam calibrator_t The type of calibrator
  /// @tparam outlier_finder_t Type of the outlier finder class
  ///
  /// The GX2F Actor does not rely on the measurements to be sorted along the
  /// track.
  class Actor {
   public:
    /// Broadcast the result_type
    using result_type = Gx2FitterResult<traj_t>;

    /// The target surface
    const Surface* targetSurface = nullptr;

    /// Allows retrieving measurements for a surface
    const std::unordered_map<const Surface*, SourceLink>* inputMeasurements{};

    /// Whether to consider multiple scattering.
    bool multipleScattering = false;

    /// Whether to consider energy loss.
    bool energyLoss = false;

    /// Whether to use the mean or the most probable energy loss.
    Gx2fEnergyLossMode energyLossMode = Gx2fEnergyLossMode::Mean;

    /// Whether to re-evaluate the deterministic q/p correction from the local
    /// q/p of this propagation. Must be disabled for the final propagation,
    /// where the applied model has to stay exactly the fitted one.
    bool refreshEnergyLossExpectation = true;

    /// Whether to include non-linear correction during global to local
    /// transformation
    FreeToBoundCorrection freeToBoundCorrection;

    /// Input MultiTrajectory
    std::shared_ptr<MultiTrajectory<traj_t>> outputStates;

    /// The logger instance
    const Logger* actorLogger{nullptr};

    /// Logger helper
    const Logger& logger() const { return *actorLogger; }

    Gx2FitterExtensions<traj_t> extensions;

    /// The Surface being
    SurfaceReached targetReached;

    /// Calibration context for the fit
    const CalibrationContext* calibrationContext{nullptr};

    /// The particle hypothesis is needed for estimating scattering angles
    const BoundTrackParameters* parametersWithHypothesis = nullptr;

    /// The materialMap stores for each visited surface their material
    /// properties
    std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties>*
        materialMap = nullptr;

    /// @brief Gx2f actor operation
    ///
    /// @tparam propagator_state_t is the type of Propagator state
    /// @tparam stepper_t Type of the stepper
    /// @tparam navigator_t Type of the navigator
    ///
    /// @param state is the mutable propagator state object
    /// @param stepper The stepper in use
    /// @param navigator The navigator in use
    /// @param result is the mutable result state object
    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t>
    Result<void> act(propagator_state_t& state, const stepper_t& stepper,
                     const navigator_t& navigator, result_type& result,
                     const Logger& /*logger*/) const {
      assert(result.fittedStates && "No MultiTrajectory set");

      // Check if we can stop to propagate
      if (result.measurementStates == inputMeasurements->size()) {
        ACTS_DEBUG("Actor: finish: All measurements have been found.");
        result.finished = true;
      } else if (state.navigation.navigationBreak) {
        ACTS_DEBUG("Actor: finish: navigationBreak.");
        result.finished = true;
      }

      // End the propagation and return to the fitter
      if (result.finished) {
        // Remove the missing surfaces that occur after the last measurement
        if (result.measurementStates > 0) {
          result.missedActiveSurfaces.resize(result.measurementHoles);
        }

        return Result<void>::success();
      }

      // We are only interested in surfaces. If we are not on a surface, we
      // continue the navigation
      auto surface = navigator.currentSurface(state.navigation);
      if (surface == nullptr) {
        return Result<void>::success();
      }

      ++result.surfaceCount;
      const GeometryIdentifier geoId = surface->geometryId();
      ACTS_DEBUG("Surface " << geoId << " detected.");

      const bool surfaceIsSensitive = surface->isSensitive();
      const bool surfaceHasMaterial = surface->hasMaterial();
      // First we figure out, if we would need to look into material surfaces at
      // all. Later, we also check, if the material slab is valid, otherwise we
      // modify this flag to ignore the material completely.
      const bool handleMaterial = multipleScattering || energyLoss;
      bool doMaterial = handleMaterial && surfaceHasMaterial;

      // Found material - add a material entry if not done yet.
      // Handling will happen later
      if (doMaterial) {
        ACTS_DEBUG("    The surface contains material, ...");

        auto materialMapId = materialMap->find(geoId);
        const bool isNewEntry = (materialMapId == materialMap->end());
        // Unlike the Highland sigma, which is only a weight, the expected
        // energy loss is a systematic shift of the model. Keeping it at the
        // seed's q/p would re-inject exactly the bias we want to remove, so we
        // re-evaluate it from the local q/p on every propagation while the fit
        // iterates.
        const bool refreshEnergyLoss =
            energyLoss && refreshEnergyLossExpectation;

        if (isNewEntry || refreshEnergyLoss) {
          const Result<MaterialSlab> slabResult =
              Acts::detail::evaluateMaterialSlab(
                  state, stepper, *surface,
                  Acts::detail::determineMaterialUpdateMode(
                      state, navigator, MaterialUpdateMode::FullUpdate));
          if (!slabResult.ok()) {
            ACTS_DEBUG("GlobalChiSquareFitter | "
                       << "Failed to evaluate material slab: "
                       << slabResult.error());
            return Result<void>::failure(slabResult.error());
          }
          const MaterialSlab& slab = *slabResult;
          const bool slabIsValid = !slab.isVacuum();

          if (isNewEntry) {
            ACTS_DEBUG("    ... create entry in material map.");

            double invSigma2 = 0.;
            if (slabIsValid) {
              const auto& particle =
                  parametersWithHypothesis->particleHypothesis();

              const double sigma =
                  static_cast<double>(Acts::computeMultipleScatteringTheta0(
                      slab, particle.absolutePdg(), particle.mass(),
                      static_cast<float>(
                          parametersWithHypothesis->parameters()[eBoundQOverP]),
                      particle.absoluteCharge()));
              ACTS_VERBOSE(
                  "        The Highland formula gives sigma = " << sigma);
              invSigma2 = 1. / std::pow(sigma, 2);
            } else {
              ACTS_VERBOSE("        Material slab is not valid.");
            }

            materialMap->emplace(
                geoId, Gx2fMaterialProperties{BoundVector::Zero(), invSigma2,
                                              slabIsValid});
            materialMapId = materialMap->find(geoId);
          } else {
            ACTS_DEBUG("    ... found entry in material map.");
          }

          if (refreshEnergyLoss && materialMapId->second.materialIsValid()) {
            const auto& particle = stepper.particleHypothesis(state.stepping);
            const double qOverPLocal = stepper.qOverP(state.stepping);

            materialMapId->second.expectedQOverPOffset() =
                computeGx2fQOverPOffset(slab, particle, qOverPLocal,
                                        state.options.direction,
                                        energyLossMode);

            const double sigmaQOverP =
                static_cast<double>(Acts::computeEnergyLossLandauSigmaQOverP(
                    slab, static_cast<float>(particle.mass()),
                    static_cast<float>(qOverPLocal),
                    static_cast<float>(particle.absoluteCharge())));

            if (!std::isfinite(sigmaQOverP) || sigmaQOverP <= 0.) {
              // Without a straggling width the energy loss parameter would be
              // unconstrained. Drop the surface entirely, so that the actor,
              // fillGx2fSystem and countMaterialStates stay consistent.
              ACTS_WARNING(
                  "Material surface "
                  << geoId
                  << " has no valid energy loss straggling. Ignoring it.");
              materialMapId->second.invalidateMaterial();
            } else {
              materialMapId->second.invCovarianceQOverP() =
                  1. / (sigmaQOverP * sigmaQOverP);
              ACTS_VERBOSE("        Expected delta(q/p) = "
                           << materialMapId->second.expectedQOverPOffset()
                           << ", sigma(q/p) = " << sigmaQOverP);
            }
          }
        } else {
          ACTS_DEBUG("    ... found entry in material map.");
        }

        doMaterial = doMaterial && materialMapId->second.materialIsValid();
      }

      // Here we handle all measurements
      if (auto sourceLinkIt = inputMeasurements->find(surface);
          sourceLinkIt != inputMeasurements->end()) {
        ACTS_DEBUG("    The surface contains a measurement.");

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // TODO generalize the update of the currentTrackIndex
        auto& fittedStates = *result.fittedStates;

        // Add a <trackStateMask> TrackState entry multi trajectory. This
        // allocates storage for all components, which we will set later.
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                        result.lastTrackIndex);
        const std::size_t currentTrackIndex = trackStateProxy.index();

        // Set the trackStateProxy components with the state from the ongoing
        // propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            return res.error();
          }
          // Not const since, we might need to update with scattering angles
          auto& [boundParams, jacobian, pathLength] = *res;

          // For material surfaces, we also update the parameters with the
          // available material information
          if (doMaterial) {
            const auto materialMapId = materialMap->find(geoId);
            ACTS_VERBOSE("        boundParams before the update: "
                         << boundParams.parameters().transpose());
            if (multipleScattering) {
              ACTS_DEBUG("    Update parameters with scattering angles.");
              ACTS_VERBOSE(
                  "        scatteringAngles: "
                  << materialMapId->second.scatteringAngles().transpose());
              boundParams.parameters() +=
                  materialMapId->second.scatteringAngles();
            }
            if (energyLoss) {
              ACTS_DEBUG("    Update q/p with the energy loss.");
              ACTS_VERBOSE("        delta(q/p): "
                           << materialMapId->second.totalQOverPOffset());
              boundParams.parameters()[eBoundQOverP] +=
                  materialMapId->second.totalQOverPOffset();
            }
            ACTS_VERBOSE("        boundParams after the update: "
                         << boundParams.parameters().transpose());
          }

          // Fill the track state
          trackStateProxy.smoothed() = boundParams.parameters();
          trackStateProxy.smoothedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;

          if (doMaterial) {
            stepper.update(state.stepping,
                           transformBoundToFreeParameters(
                               trackStateProxy.referenceSurface(),
                               state.geoContext, trackStateProxy.smoothed()),
                           trackStateProxy.smoothed(),
                           trackStateProxy.smoothedCovariance(), *surface);
          }
        }

        // We have smoothed parameters, so calibrate the uncalibrated input
        // measurement
        extensions.calibrator(state.geoContext, *calibrationContext,
                              sourceLinkIt->second, trackStateProxy);

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.setHasParameters();
        if (surfaceHasMaterial) {
          typeFlags.setHasMaterial();
        }

        // Set the measurement type flag
        typeFlags.setIsMeasurement();
        // We count the processed measurement
        ++result.processedMeasurements;

        result.lastMeasurementIndex = currentTrackIndex;
        result.lastTrackIndex = currentTrackIndex;

        // TODO check for outlier first
        // We count the state with measurement
        ++result.measurementStates;

        // We count the processed state
        ++result.processedStates;

        // Update the number of holes count only when encountering a
        // measurement
        result.measurementHoles = result.missedActiveSurfaces.size();

        return Result<void>::success();
      }

      if (doMaterial) {
        // Here we handle material for multipleScattering. If holes exist, we
        // also handle them already. We create a full trackstate (unlike for
        // simple holes), since we need to evaluate the material later
        ACTS_DEBUG(
            "    The surface contains no measurement, but material and maybe "
            "a hole.");

        // Transport the covariance to the surface
        stepper.transportCovarianceToBound(state.stepping, *surface,
                                           freeToBoundCorrection);

        // TODO generalize the update of the currentTrackIndex
        auto& fittedStates = *result.fittedStates;

        // Add a <trackStateMask> TrackState entry multi trajectory. This
        // allocates storage for all components, which we will set later.
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                        result.lastTrackIndex);
        const std::size_t currentTrackIndex = trackStateProxy.index();

        // Set the trackStateProxy components with the state from the ongoing
        // propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            return res.error();
          }
          // Not const since, we might need to update with scattering angles
          auto& [boundParams, jacobian, pathLength] = *res;

          // For material surfaces, we also update the angles with the
          // available scattering information
          // We know already, that we handle material here, but we still need to
          // distinguish which of the two effects are enabled.
          const auto materialMapId = materialMap->find(geoId);
          ACTS_VERBOSE("        boundParams before the update: "
                       << boundParams.parameters().transpose());
          if (multipleScattering) {
            ACTS_DEBUG("    Update parameters with scattering angles.");
            ACTS_VERBOSE(
                "        scatteringAngles: "
                << materialMapId->second.scatteringAngles().transpose());
            boundParams.parameters() +=
                materialMapId->second.scatteringAngles();
          }
          if (energyLoss) {
            ACTS_DEBUG("    Update q/p with the energy loss.");
            ACTS_VERBOSE("        delta(q/p): "
                         << materialMapId->second.totalQOverPOffset());
            boundParams.parameters()[eBoundQOverP] +=
                materialMapId->second.totalQOverPOffset();
          }
          ACTS_VERBOSE("        boundParams after the update: "
                       << boundParams.parameters().transpose());

          // Fill the track state
          trackStateProxy.smoothed() = boundParams.parameters();
          trackStateProxy.smoothedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;

          stepper.update(state.stepping,
                         transformBoundToFreeParameters(
                             trackStateProxy.referenceSurface(),
                             state.geoContext, trackStateProxy.smoothed()),
                         trackStateProxy.smoothed(),
                         trackStateProxy.smoothedCovariance(), *surface);
        }

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.setHasParameters();
        typeFlags.setHasMaterial();

        // Set hole only, if we are on a sensitive surface and had
        // measurements before (no holes before the first measurement)
        const bool precedingMeasurementExists = (result.measurementStates > 0);
        if (surfaceIsSensitive && precedingMeasurementExists) {
          ACTS_DEBUG("    Surface is also sensitive. Marked as hole.");
          typeFlags.setIsHole();

          // Count the missed surface
          result.missedActiveSurfaces.push_back(surface);
        }

        result.lastTrackIndex = currentTrackIndex;

        ++result.processedStates;

        return Result<void>::success();
      }

      if (surfaceIsSensitive || surfaceHasMaterial) {
        // Here we handle holes. If material hasn't been handled before
        // (because both material effects are turned off), we will also handle
        // it here
        if (handleMaterial) {
          ACTS_DEBUG(
              "    The surface contains no measurement, but maybe a hole.");
        } else {
          ACTS_DEBUG(
              "    The surface contains no measurement, but maybe a hole "
              "and/or material.");
        }

        // We only create track states here if there is already a measurement
        // detected (no holes before the first measurement) or if we encounter
        // material
        const bool precedingMeasurementExists = (result.measurementStates > 0);
        if (!precedingMeasurementExists && !surfaceHasMaterial) {
          ACTS_DEBUG(
              "    Ignoring hole, because there are no preceding "
              "measurements.");
          return Result<void>::success();
        }

        auto& fittedStates = *result.fittedStates;

        // Add a <trackStateMask> TrackState entry multi trajectory. This
        // allocates storage for all components, which we will set later.
        typename traj_t::TrackStateProxy trackStateProxy =
            fittedStates.makeTrackState(Gx2fConstants::trackStateMask,
                                        result.lastTrackIndex);
        const std::size_t currentTrackIndex = trackStateProxy.index();

        // Set the trackStateProxy components with the state from the
        // ongoing propagation
        {
          trackStateProxy.setReferenceSurface(surface->getSharedPtr());
          // Bind the transported state to the current surface
          auto res = stepper.boundState(state.stepping, *surface, false,
                                        freeToBoundCorrection);
          if (!res.ok()) {
            return res.error();
          }
          const auto& [boundParams, jacobian, pathLength] = *res;

          // Fill the track state
          trackStateProxy.smoothed() = boundParams.parameters();
          trackStateProxy.smoothedCovariance() = state.stepping.cov;

          trackStateProxy.jacobian() = jacobian;
          trackStateProxy.pathLength() = pathLength;
        }

        // Get and set the type flags
        auto typeFlags = trackStateProxy.typeFlags();
        typeFlags.setHasParameters();
        if (surfaceHasMaterial) {
          ACTS_DEBUG("    It is material.");
          typeFlags.setHasMaterial();
        }

        // Set hole only, if we are on a sensitive surface
        if (surfaceIsSensitive && precedingMeasurementExists) {
          ACTS_DEBUG("    It is a hole.");
          typeFlags.setIsHole();
          // Count the missed surface
          result.missedActiveSurfaces.push_back(surface);
        }

        result.lastTrackIndex = currentTrackIndex;

        ++result.processedStates;

        return Result<void>::success();
      }

      ACTS_DEBUG("    The surface contains no measurement/material/hole.");
      return Result<void>::success();
    }

    template <typename propagator_state_t, typename stepper_t,
              typename navigator_t, typename result_t>
    bool checkAbort(propagator_state_t& /*state*/, const stepper_t& /*stepper*/,
                    const navigator_t& /*navigator*/, const result_t& result,
                    const Logger& /*logger*/) const {
      if (result.finished) {
        return true;
      }
      return false;
    }
  };

 public:
  /// Fit implementation
  ///
  /// @tparam source_link_iterator_t Iterator type used to pass source links
  /// @tparam track_container_t Type of the track container backend
  /// @tparam holder_t Type defining track container backend ownership
  ///
  /// @param it Begin iterator for the fittable uncalibrated measurements
  /// @param end End iterator for the fittable uncalibrated measurements
  /// @param sParameters The initial track parameters
  /// @param gx2fOptions Gx2FitterOptions steering the fit
  /// @param trackContainer Input track container storage to append into
  /// @note The input measurements are given in the form of @c SourceLink s.
  /// It's the calibrators job to turn them into calibrated measurements used in
  /// the fit.
  ///
  /// @return the output as an output track
  template <typename source_link_iterator_t,
            TrackContainerFrontend track_container_t>
  Result<typename track_container_t::TrackProxy> fit(
      source_link_iterator_t it, source_link_iterator_t end,
      const BoundTrackParameters& sParameters,
      const Gx2FitterOptions<traj_t>& gx2fOptions,
      track_container_t& trackContainer) const
    requires(!isDirectNavigator)
  {
    // Preprocess Measurements (SourceLinks -> map)
    // To be able to find measurements later, we put them into a map.
    // We need to copy input SourceLinks anyway, so the map can own them.
    ACTS_VERBOSE("Preparing " << std::distance(it, end)
                              << " input measurements");
    std::unordered_map<const Surface*, SourceLink> inputMeasurements{};

    for (; it != end; ++it) {
      inputMeasurements.try_emplace(gx2fOptions.extensions.surfaceAccessor(*it),
                                    *it);
    }

    // Store, which material effects we want to consider. We still need to pass
    // these options to the Actor.
    const bool multipleScattering = gx2fOptions.multipleScattering;
    const bool energyLoss = gx2fOptions.energyLoss;

    // Create the ActorList
    using GX2FActor = Actor;

    using GX2FResult = typename GX2FActor::result_type;
    using Actors = Acts::ActorList<GX2FActor>;

    using PropagatorOptions = typename propagator_t::template Options<Actors>;

    BoundTrackParameters params = sParameters;
    double chi2sum = 0;
    double oldChi2sum = std::numeric_limits<double>::max();

    // We need to create a temporary track container. We create several times a
    // new track and delete it after updating the parameters. However, if we
    // would work on the externally provided track container, it would be
    // difficult to remove the correct track, if it contains more than one.
    typename track_container_t::TrackContainerBackend trackContainerTempBackend;
    traj_t trajectoryTempBackend;
    TrackContainer trackContainerTemp{trackContainerTempBackend,
                                      trajectoryTempBackend};

    // Create an index of the 'tip' of the track stored in multitrajectory. It
    // is needed outside the update loop. It will be updated with each iteration
    // and used for the final track
    std::size_t tipIndex = kInvalid;

    // The materialMap stores for each visited surface their material
    // properties
    std::unordered_map<GeometryIdentifier, Gx2fMaterialProperties> materialMap;

    // This will be filled during the updates with the final covariance of the
    // track parameters.
    BoundMatrix fullCovariancePredicted = BoundMatrix::Identity();

    ACTS_VERBOSE("Initial parameters: " << params.parameters().transpose());

    /// Actual Fitting /////////////////////////////////////////////////////////
    ACTS_DEBUG("Start to iterate");

    // Iterate the fit and improve result. Abort after n steps or after
    // convergence.
    // nUpdate is initialized outside to save its state for the track.
    std::size_t nUpdate = 0;
    for (nUpdate = 0; nUpdate < gx2fOptions.nUpdateMax; nUpdate++) {
      ACTS_DEBUG("nUpdate = " << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax);

      // set up propagator and co
      PropagatorOptions propagatorOptions{gx2fOptions.propagatorPlainOptions};

      // Add the measurement surface as external surface to the navigator.
      // We will try to hit those surface by ignoring boundary checks.
      for (const auto& [surface, _] : inputMeasurements) {
        propagatorOptions.navigation.appendExternalSurface(*surface);
      }

      auto& gx2fActor = propagatorOptions.actorList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.multipleScattering = false;
      // Unlike the scattering angles, which start at 0 and therefore leave the
      // material-free main loop untouched, the energy loss is a non-zero
      // systematic bias from the first iteration. We apply its deterministic
      // part here, but do not yet fit the deviation from it.
      gx2fActor.energyLoss = energyLoss;
      gx2fActor.energyLossMode = gx2fOptions.energyLossMode;
      gx2fActor.refreshEnergyLossExpectation = true;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.materialMap = &materialMap;
      gx2fActor.parametersWithHypothesis = &params;

      auto propagatorState = m_propagator.makeState(propagatorOptions);

      auto propagatorInitResult =
          m_propagator.initialize(propagatorState, params);
      if (!propagatorInitResult.ok()) {
        ACTS_DEBUG("Propagation initialization failed: "
                   << propagatorInitResult.error());
        return propagatorInitResult.error();
      }

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trajectoryTempBackend;

      // Clear the track container. It could be more performant to update the
      // existing states, but this needs some more thinking.
      trackContainerTemp.clear();

      // Run the fitter
      auto propagationResult = m_propagator.propagate(propagatorState);

      auto result =
          m_propagator.makeResult(std::move(propagatorState), propagationResult,
                                  propagatorOptions, false);

      if (!result.ok()) {
        ACTS_DEBUG("Propagation failed: " << result.error());
        return result.error();
      }

      // TODO Improve Propagator + Actor [allocate before loop], rewrite
      // makeMeasurements
      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      auto track = trackContainerTemp.makeTrack();
      tipIndex = gx2fResult.lastMeasurementIndex;

      // It could happen, that no measurements were found. Then the track would
      // be empty and the following operations would be invalid. Usually, this
      // only happens during the first iteration, due to bad initial parameters.
      if (tipIndex == kInvalid) {
        ACTS_INFO("Did not find any measurements in nUpdate "
                  << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax);
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      track.tipIndex() = tipIndex;
      track.linkForward();

      // No material parameters are fitted in the main loop: the scattering
      // angles are held at 0.
      const Gx2fParameterLayout layout{/*fitScattering_=*/false,
                                       /*fitEnergyLoss_=*/false,
                                       /*nMaterialSurfaces_=*/0u};

      // System that we fill with the information gathered by the actor and
      // evaluate later
      Gx2fSystem extendedSystem{layout};

      // This vector stores the IDs for each visited material in order. We use
      // it later for updating the scattering angles. We cannot use
      // materialMap directly, since we cannot guarantee, that we will visit
      // all stored material in each propagation.
      std::vector<GeometryIdentifier> geoIdVector;

      // The actor handles material whenever energy loss is enabled, and resets
      // the transport Jacobian on every surface where it does. Those states
      // must be chained here even though they contribute no free parameters.
      fillGx2fSystem(track, extendedSystem, /*handleMaterial=*/energyLoss,
                     materialMap, geoIdVector, *m_addToSumLogger);

      chi2sum = extendedSystem.chi2();

      // This check takes into account the evaluated dimensions of the
      // measurements. To fit, we need at least NDF+1 measurements. However, we
      // count n-dimensional measurements for n measurements, reducing the
      // effective number of needed measurements. We might encounter the case,
      // where we cannot use some (parts of a) measurements, maybe if we do not
      // support that kind of measurement. This is also taken into account here.
      // We skip the check during the first iteration, since we cannot guarantee
      // to hit all/enough measurement surfaces with the initial parameter
      // guess.
      // We skip the check during the first iteration, since we cannot guarantee
      // to hit all/enough measurement surfaces with the initial parameter
      // guess.
      if ((nUpdate > 0) && !extendedSystem.isWellDefined()) {
        ACTS_INFO("Not enough measurements. Require "
                  << extendedSystem.findRequiredNdf() + 1 << ", but only "
                  << extendedSystem.ndf() << " could be used.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      Eigen::VectorXd deltaParamsExtended =
          computeGx2fDeltaParams(extendedSystem);

      ACTS_VERBOSE("aMatrix:\n"
                   << extendedSystem.aMatrix() << "\n"
                   << "bVector:\n"
                   << extendedSystem.bVector() << "\n"
                   << "deltaParamsExtended:\n"
                   << deltaParamsExtended << "\n"
                   << "oldChi2sum = " << oldChi2sum << "\n"
                   << "chi2sum = " << extendedSystem.chi2());

      if ((gx2fOptions.relChi2changeCutOff != 0) && (nUpdate > 0) &&
          (std::abs(extendedSystem.chi2() / oldChi2sum - 1) <
           gx2fOptions.relChi2changeCutOff)) {
        ACTS_DEBUG("Abort with relChi2changeCutOff after "
                   << nUpdate + 1 << "/" << gx2fOptions.nUpdateMax
                   << " iterations.");
        updateGx2fCovarianceParams(fullCovariancePredicted, extendedSystem);
        break;
      }

      if (extendedSystem.chi2() > oldChi2sum + 1e-5) {
        ACTS_DEBUG("chi2 not converging monotonically in update " << nUpdate);
      }

      // If this is the final iteration, update the covariance and break.
      // Otherwise, we would update the scattering angles too much.
      if (nUpdate == gx2fOptions.nUpdateMax - 1) {
        // Since currently most of our tracks converge in 4-5 updates, we want
        // to set nUpdateMax higher than that to guarantee convergence for most
        // tracks. In cases, where we set a smaller nUpdateMax, it's because we
        // want to investigate the behaviour of the fitter before it converges,
        // like in some unit-tests.
        if (gx2fOptions.nUpdateMax > 5) {
          ACTS_INFO("Did not converge in " << gx2fOptions.nUpdateMax
                                           << " updates.");
          return Experimental::GlobalChiSquareFitterError::DidNotConverge;
        }

        updateGx2fCovarianceParams(fullCovariancePredicted, extendedSystem);
        break;
      }

      updateGx2fParams(params, deltaParamsExtended, layout, materialMap,
                       geoIdVector);
      ACTS_VERBOSE("Updated parameters: " << params.parameters().transpose());

      oldChi2sum = extendedSystem.chi2();
    }
    ACTS_DEBUG("Finished to iterate");
    ACTS_VERBOSE("Final parameters: " << params.parameters().transpose());
    /// Finish Fitting /////////////////////////////////////////////////////////

    /// Actual MATERIAL Fitting ////////////////////////////////////////////////
    ACTS_DEBUG("Start to evaluate material");
    for (std::size_t nMaterialUpdate = 0;
         (multipleScattering || energyLoss) &&
         nMaterialUpdate < gx2fOptions.nMaterialUpdateMax;
         nMaterialUpdate++) {
      ACTS_DEBUG("nMaterialUpdate = " << nMaterialUpdate + 1 << "/"
                                      << gx2fOptions.nMaterialUpdateMax);

      // Set up the propagator
      PropagatorOptions propagatorOptions{gx2fOptions.propagatorPlainOptions};

      // Add the measurement surface as external surface to the navigator.
      // We will try to hit those surface by ignoring boundary checks.
      for (const auto& [surface, _] : inputMeasurements) {
        propagatorOptions.navigation.appendExternalSurface(*surface);
      }

      auto& gx2fActor = propagatorOptions.actorList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.multipleScattering = multipleScattering;
      gx2fActor.energyLoss = energyLoss;
      gx2fActor.energyLossMode = gx2fOptions.energyLossMode;
      // Relinearise the deterministic energy loss around the parameters of this
      // iteration
      gx2fActor.refreshEnergyLossExpectation = true;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.materialMap = &materialMap;
      gx2fActor.parametersWithHypothesis = &params;

      auto propagatorState = m_propagator.makeState(propagatorOptions);

      auto propagatorInitResult =
          m_propagator.initialize(propagatorState, params);
      if (!propagatorInitResult.ok()) {
        ACTS_DEBUG("Propagation initialization failed: "
                   << propagatorInitResult.error());
        return propagatorInitResult.error();
      }

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trajectoryTempBackend;

      // Clear the track container. It could be more performant to update the
      // existing states, but this needs some more thinking.
      trackContainerTemp.clear();

      // Run the fitter
      auto propagationResult = m_propagator.propagate(propagatorState);

      auto result =
          m_propagator.makeResult(std::move(propagatorState), propagationResult,
                                  propagatorOptions, false);

      if (!result.ok()) {
        ACTS_DEBUG("Propagation failed: " << result.error());
        return result.error();
      }

      // TODO Improve Propagator + Actor [allocate before loop], rewrite
      // makeMeasurements
      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      auto track = trackContainerTemp.makeTrack();
      tipIndex = gx2fResult.lastMeasurementIndex;

      // It could happen, that no measurements were found. Then the track would
      // be empty and the following operations would be invalid. Usually, this
      // only happens during the first iteration, due to bad initial parameters.
      if (tipIndex == kInvalid) {
        ACTS_INFO("Did not find any measurements in material fit.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      track.tipIndex() = tipIndex;
      track.linkForward();

      // Count the material surfaces, to set up the system. In the multiple
      // scattering case, we need to extend our system.
      const std::size_t nMaterialSurfaces =
          countMaterialStates(track, materialMap, *m_addToSumLogger);

      // We need 6 dimensions for the bound parameters and additional dimensions
      // for the parameters of each material surface.
      const Gx2fParameterLayout layout{/*fitScattering_=*/multipleScattering,
                                       /*fitEnergyLoss_=*/energyLoss,
                                       nMaterialSurfaces};

      // System that we fill with the information gathered by the actor and
      // evaluate later
      Gx2fSystem extendedSystem{layout};

      // This vector stores the IDs for each visited material in order. We use
      // it later for updating the scattering angles. We cannot use
      // materialMap directly, since we cannot guarantee, that we will visit
      // all stored material in each propagation.
      std::vector<GeometryIdentifier> geoIdVector;

      fillGx2fSystem(track, extendedSystem, /*handleMaterial=*/true,
                     materialMap, geoIdVector, *m_addToSumLogger);

      chi2sum = extendedSystem.chi2();

      // This check takes into account the evaluated dimensions of the
      // measurements. To fit, we need at least NDF+1 measurements. However, we
      // count n-dimensional measurements for n measurements, reducing the
      // effective number of needed measurements. We might encounter the case,
      // where we cannot use some (parts of a) measurements, maybe if we do not
      // support that kind of measurement. This is also taken into account here.
      // We skip the check during the first iteration, since we cannot guarantee
      // to hit all/enough measurement surfaces with the initial parameter
      // guess.
      if ((nUpdate > 0) && !extendedSystem.isWellDefined()) {
        ACTS_INFO("Not enough measurements. Require "
                  << extendedSystem.findRequiredNdf() + 1 << ", but only "
                  << extendedSystem.ndf() << " could be used.");
        return Experimental::GlobalChiSquareFitterError::NotEnoughMeasurements;
      }

      Eigen::VectorXd deltaParamsExtended =
          computeGx2fDeltaParams(extendedSystem);

      ACTS_VERBOSE("aMatrix:\n"
                   << extendedSystem.aMatrix() << "\n"
                   << "bVector:\n"
                   << extendedSystem.bVector() << "\n"
                   << "deltaParamsExtended:\n"
                   << deltaParamsExtended << "\n"
                   << "oldChi2sum = " << oldChi2sum << "\n"
                   << "chi2sum = " << extendedSystem.chi2());

      chi2sum = extendedSystem.chi2();

      updateGx2fParams(params, deltaParamsExtended, layout, materialMap,
                       geoIdVector);
      ACTS_VERBOSE("Updated parameters: " << params.parameters().transpose());

      updateGx2fCovarianceParams(fullCovariancePredicted, extendedSystem);
    }
    ACTS_DEBUG("Finished to evaluate material");
    ACTS_VERBOSE(
        "Final parameters after material: " << params.parameters().transpose());
    /// Finish MATERIAL Fitting ////////////////////////////////////////////////

    ACTS_VERBOSE("Final material parameters:");
    for (const auto& [key, value] : materialMap) {
      if (!value.materialIsValid()) {
        continue;
      }
      const auto& angles = value.scatteringAngles();
      ACTS_VERBOSE("    theta = "
                   << angles[eBoundTheta] << " | phi = " << angles[eBoundPhi]
                   << " | expected delta(q/p) = "
                   << value.expectedQOverPOffset()
                   << " | fitted delta(q/p) = " << value.qOverPOffset());
    }

    ACTS_VERBOSE("Final covariance:\n" << fullCovariancePredicted);

    // Propagate again with the final covariance matrix. This is necessary to
    // obtain the propagated covariance for each state.
    // We also need to recheck the result and find the tipIndex, because at this
    // step, we will not ignore the boundary checks for measurement surfaces. We
    // want to create trackstates only on surfaces, that we actually hit.
    if (gx2fOptions.nUpdateMax > 0) {
      ACTS_VERBOSE("Propagate with the final covariance.");
      // update covariance
      params.covariance() = fullCovariancePredicted;

      // set up the propagator
      PropagatorOptions propagatorOptions{gx2fOptions.propagatorPlainOptions};
      auto& gx2fActor = propagatorOptions.actorList.template get<GX2FActor>();
      gx2fActor.inputMeasurements = &inputMeasurements;
      gx2fActor.multipleScattering = multipleScattering;
      gx2fActor.energyLoss = energyLoss;
      gx2fActor.energyLossMode = gx2fOptions.energyLossMode;
      // Freeze the fitted model. If we re-evaluated the expectation here, the
      // output track states would carry a different q/p kink than the one the
      // reported parameters and covariance were computed with.
      gx2fActor.refreshEnergyLossExpectation = false;
      gx2fActor.extensions = gx2fOptions.extensions;
      gx2fActor.calibrationContext = &gx2fOptions.calibrationContext.get();
      gx2fActor.actorLogger = m_actorLogger.get();
      gx2fActor.materialMap = &materialMap;
      gx2fActor.parametersWithHypothesis = &params;

      auto propagatorState = m_propagator.makeState(propagatorOptions);

      auto propagatorInitResult =
          m_propagator.initialize(propagatorState, params);
      if (!propagatorInitResult.ok()) {
        ACTS_DEBUG("Propagation initialization failed: "
                   << propagatorInitResult.error());
        return propagatorInitResult.error();
      }

      auto& r = propagatorState.template get<Gx2FitterResult<traj_t>>();
      r.fittedStates = &trackContainer.trackStateContainer();

      // Run the fitter
      auto propagationResult = m_propagator.propagate(propagatorState);

      auto result =
          m_propagator.makeResult(std::move(propagatorState), propagationResult,
                                  propagatorOptions, false);

      if (!result.ok()) {
        ACTS_DEBUG("Propagation failed: " << result.error());
        return result.error();
      }

      auto& propRes = *result;
      GX2FResult gx2fResult = std::move(propRes.template get<GX2FResult>());

      if (tipIndex != gx2fResult.lastMeasurementIndex) {
        ACTS_INFO("Final fit used unreachable measurements.");
        tipIndex = gx2fResult.lastMeasurementIndex;

        // It could happen, that no measurements were found. Then the track
        // would be empty and the following operations would be invalid.
        if (tipIndex == kInvalid) {
          ACTS_INFO("Did not find any measurements in final propagation.");
          return Experimental::GlobalChiSquareFitterError::
              NotEnoughMeasurements;
        }
      }
    }

    if (!trackContainer.hasColumn(
            Acts::hashString(Gx2fConstants::gx2fnUpdateColumn))) {
      trackContainer.template addColumn<std::uint32_t>("Gx2fnUpdateColumn");
    }

    // Prepare track for return
    auto track = trackContainer.makeTrack();
    track.tipIndex() = tipIndex;
    track.parameters() = params.parameters();
    track.covariance() = fullCovariancePredicted;
    track.setReferenceSurface(params.referenceSurface().getSharedPtr());

    if (trackContainer.hasColumn(
            Acts::hashString(Gx2fConstants::gx2fnUpdateColumn))) {
      ACTS_DEBUG("Add nUpdate to track");
      track.template component<std::uint32_t>("Gx2fnUpdateColumn") =
          static_cast<std::uint32_t>(nUpdate);
    }

    // TODO write test for calculateTrackQuantities
    calculateTrackQuantities(track);

    // Set the chi2sum for the track summary manually, since we don't calculate
    // it for each state
    track.chi2() = chi2sum;

    track.linkForward();

    // Return the converted Track
    return track;
  }
};

/// @}

}  // namespace Acts::Experimental
