#!/usr/bin/env python3
from pathlib import Path
import numpy as np
import tempfile
import sys
import csv
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

def run_vertexing(mu, parts_per_vtx):
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
            ParticleConfig(parts_per_vtx, acts.PdgParticle.eMuon, randomizeCharge=True),
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(10 * u.um, 10 * u.um, 50 * u.mm, 0),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            multiplicity=mu,
            rnd=rnd,
            outputDirCsv=outputDir,
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

        s.addWriter(acts.examples.CsvTrackParameterWriter(level=acts.logging.INFO,
                                                              inputTrajectories="trajectories", 
                                                              outputDir=outputDir, 
                                                              outputStem="tracks.csv"))
        """
        addVertexFitting(
            s,
            field,
            vertexFinder=fitter,
            outputDirRoot=tp,
            minWeight=minWeight,
            maxIterations=maxIterations,
        )
        """
        s.run()

        del s


with acts.FpeMonitor():
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

    outputDir = "/home/frusso/hep/acts/Examples/Scripts/out/unit_test_AMVF/"
    mu = 3
    parts_per_vtx = 4
    run_vertexing(mu=mu, parts_per_vtx=parts_per_vtx)
    newOutputDir = "/home/frusso/hep/acts/Tests/Data/"

    vtx_file_old = outputDir + "event000000000-particles.csv"
    vtx_file_new = newOutputDir + "vertexing_event_mu"+str(mu)+"_vertices_truth.csv"
    with open(vtx_file_old) as oldcsvfile, open(vtx_file_new, mode="w") as newcsvfile:
        csv_reader = csv.reader(oldcsvfile, delimiter=',')
        csv_writer = csv.writer(newcsvfile, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count==0:
                csv_writer.writerow(["posX", "posY", "posZ", 
                                     "covXX", "covXY", "covXZ", 
                                     "covYX", "covYY", "covYZ", 
                                     "covZX", "covZY", "covZZ", 
                                     "nTracks", "trkWeight(1)", 
                                     "trkComp(1)", "trkChi2(1)"])
            elif ((line_count-1)%parts_per_vtx==0):
                csv_writer.writerow([row[3], row[4], row[5], 
                                     0.0, 0.0, 0.0, 
                                     0.0, 0.0, 0.0, 
                                     0.0, 0.0, 0.0, 
                                     parts_per_vtx, 0.0, 
                                     0.0, 0.0])
            line_count += 1
        
    trk_file_old = outputDir + "event000000000-tracks.csv"
    trk_file_new = newOutputDir + "vertexing_event_mu"+str(mu)+"_tracks.csv"
    conv_fac = 1e-3
    with open(trk_file_old) as oldcsvfile, open(trk_file_new, mode="w") as newcsvfile:
        csv_reader = csv.reader(oldcsvfile, delimiter=',')
        csv_writer = csv.writer(newcsvfile, delimiter=',')
        line_count = 0
        for row in csv_reader:
            if line_count == 0:
                csv_writer.writerow([row[0], row[1], row[2], row[3], row[4], row[5],
                                    row[6], row[12], row[13], row[14], row[15], row[16],
                                    row[7], row[18], row[19], row[20], row[21], 
                                    row[8], row[24], row[25], row[26], 
                                    row[9], row[30], row[31], 
                                    row[10], row[36], 
                                    row[11]])
            else:
                csv_writer.writerow([row[0], row[1], row[2], row[3], str(float(row[4])*conv_fac), row[5],
                                    row[6], row[12], row[13], row[14], str(float(row[15])*conv_fac), row[16],
                                    row[7], row[18], row[19], str(float(row[20])*conv_fac), row[21], 
                                    row[8], row[24], str(float(row[25])*conv_fac), row[26], 
                                    row[9], str(float(row[30])*conv_fac), row[31], 
                                    str(float(row[10])*conv_fac**2), str(float(row[36])*conv_fac), 
                                    row[11]])
            line_count += 1
    