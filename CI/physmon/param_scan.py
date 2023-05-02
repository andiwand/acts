#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tempfile
import sys
import uproot
from benderopt import minimize

sys.path += [
    str(Path(__file__).parent.parent.parent / "Examples/Scripts/Python/"),
]

# this has to happen before we import the ACTS module
import acts.examples

# @TODO: Fix failure in gain matrix smoothing
# See https://github.com/acts-project/acts/issues/1215
acts.logging.setFailureThreshold(acts.logging.FATAL)

from truth_tracking_kalman import runTruthTrackingKalman
from truth_tracking_gsf import runTruthTrackingGsf
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    PhiConfig,
    ParticleConfig,
    addFatras,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    ParticleSmearingSigmas,
    SeedFinderConfigArg,
    SeedFinderOptionsArg,
    SeedingAlgorithm,
    TruthEstimatedSeedingAlgorithmConfigArg,
    addCKFTracks,
    CKFPerformanceConfig,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addVertexFitting,
    VertexFinder,
    TrackSelectorRanges,
)

srcdir = Path(__file__).resolve().parent.parent.parent

u = acts.UnitConstants

matDeco = acts.IMaterialDecorator.fromFile(
    srcdir / "thirdparty/OpenDataDetector/data/odd-material-maps.root",
    level=acts.logging.INFO,
)
detector, trackingGeometry, decorators = getOpenDataDetector(
    getOpenDataDetectorDirectory(), matDeco
)
digiConfig = srcdir / "thirdparty/OpenDataDetector/config/odd-digi-smearing-config.json"
geoSel = srcdir / "thirdparty/OpenDataDetector/config/odd-seeding-config.json"

field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

def run_vertexing(minWeight):
    maxIterations = 1
    fitter = VertexFinder.AMVF
    mu = 3
    events = 1
    with tempfile.TemporaryDirectory() as temp:
        s = acts.examples.Sequencer(
            events=events, numThreads=-1, logLevel=acts.logging.INFO
        )

        tp = Path(temp)

        for d in decorators:
            s.addContextDecorator(d)

        rnd = acts.examples.RandomNumbers(seed=43)

        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
            EtaConfig(-3.0, 3.0),
            ParticleConfig(4, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 0),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=mu,
            rnd=rnd,
            outputDirRoot="/home/frusso/hep/acts/Examples/Scripts/out/vertex_mu_scan/",
        )

        addFatras(
            s,
            trackingGeometry,
            field,
            rnd=rnd,
        )

        addDigitization(
            s,
            trackingGeometry,
            field,
            digiConfigFile=digiConfig,
            rnd=rnd,
        )

        addSeeding(
            s,
            trackingGeometry,
            field,
            SeedFinderConfigArg(
                r=(None, 200 * u.mm),  # rMin=default, 33mm
                deltaR=(1 * u.mm, 60 * u.mm),
                collisionRegion=(-250 * u.mm, 250 * u.mm),
                z=(-2000 * u.mm, 2000 * u.mm),
                maxSeedsPerSpM=1,
                sigmaScattering=5,
                radLengthPerSeed=0.1,
                minPt=500 * u.MeV,
                impactMax=3 * u.mm,
            ),
            SeedFinderOptionsArg(bFieldInZ=1.99724 * u.T),
            seedingAlgorithm=SeedingAlgorithm.Default,
            geoSelectionConfigFile=geoSel,
            rnd=rnd,  # only used by SeedingAlgorithm.TruthSmeared
        )

        addCKFTracks(
            s,
            trackingGeometry,
            field,
            CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
            TrackSelectorRanges(
                pt=(500 * u.MeV, None),
                loc0=(-4.0 * u.mm, 4.0 * u.mm),
                absEta=(None, 3.0),
            ),
        )

        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(maximumSharedHits=3),
            CKFPerformanceConfig(ptMin=400.0 * u.MeV, nMeasurementsMin=6),
        )

        addVertexFitting(
            s,
            field,
            vertexFinder=fitter,
            outputDirRoot=tp,
            minWeight=minWeight,
            maxIterations=maxIterations,
        )

        s.run()

        del s

        perf_file = tp / f"performance_vertexing.root"
        assert perf_file.exists(), "Performance file not found"
        #print("\n\n\n\n"+f"{perf_file}:vertexing"+"\n\n\n\n")
        rf = uproot.open(f"{perf_file}:vertexing")
        nreco = rf["nRecoVtx"].array(library="np")
        print("nreco: ", np.mean(nreco))
        return -np.mean(nreco)


with acts.FpeMonitor():

    '''
    ### VERTEX MU SCAN
    events = 5
    mus = (100,)#(1, 10, 25, 50, 75, 100, 125, 150, 175, 200)
    fitters = (VertexFinder.AMVF,)#(VertexFinder.Iterative, VertexFinder.AMVF)
    nrecos = np.zeros((len(fitters), len(mus), events))
    for fitter_ind, fitter in enumerate(fitters):
        for mu_ind, mu in enumerate(mus):
            start = datetime.datetime.now()

            nrecos[fitter_ind, mu_ind, :] = run_vertexing(fitter, mu, events)
    
    for mu_ind, mu in enumerate(mus):
        print(f'for mu = {mu}:')
        for fitter_ind, fitter in enumerate(fitters):
            print(f'avg reconstructed vertices for {fitter.name}: ', np.mean(nrecos[fitter_ind, mu_ind]))
    '''
    '''
    # We define the parameters we want to optimize:
    optimization_problem_parameters = [
        {
            "name": "minWeight", 
            "category": "uniform",
            "search_space": {
                "low": 0.0001,
                "high": 0.1,
            }
        }
    ]

    # We launch the optimization
    best_sample = minimize(run_vertexing, optimization_problem_parameters, number_of_evaluation=1)
    #opt = -run_vertexing(best_sample["minWeight"])
    #print(best_sample["minWeight"], opt)
    '''
    run_vertexing(0.0001)
