#!/usr/bin/env python3
import os, argparse, pathlib, contextlib, acts, acts.examples
from acts.examples.simulation import (
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    addFatras,
    addGeant4,
    ParticleSelectorConfig,
    addDigitization,
)
from acts.examples.reconstruction import (
    addSeeding,
    TruthSeedRanges,
    addCKFTracks,
    CKFPerformanceConfig,
    TrackSelectorRanges,
    addAmbiguityResolution,
    AmbiguityResolutionConfig,
    addAmbiguityResolutionML,
    AmbiguityResolutionMLConfig,
    addVertexFitting,
    VertexFinder,
)
from common import getOpenDataDetectorDirectory
from acts.examples.odd import getOpenDataDetector

#odd = open data detector

parser = argparse.ArgumentParser(description="Full chain with the OpenDataDetector") #for treating the command line arguments

parser.add_argument("--events", "-n", help="Number of events", type=int, default=100) #can specify number of events using --events or n flag
parser.add_argument(
    "--geant4", help="Use Geant4 instead of fatras", action="store_true"
) #
parser.add_argument(
    "--ttbar",
    help="Use Pythia8 (ttbar, pile-up 200) instead of particle gun",
    action="store_true",
) #if ttbar flag is set, args["ttbar"] is set to True, otherwise to False
parser.add_argument(
    "--MLSolver",
    help="Use the Ml Ambiguity Solver instead of the classical one",
    action="store_true",
)

args = vars(parser.parse_args()) #parse arguments

ttbar_pu200 = args["ttbar"] #top antitop produktion (then Pythia simulation)
g4_simulation = args["geant4"]
ambiguity_MLSolver = args["MLSolver"]
u = acts.UnitConstants #todo
geoDir = getOpenDataDetectorDirectory() #path to detector directory
outputDir = pathlib.Path.cwd() / "odd_output" #path to output directory

# acts.examples.dump_args_calls(locals())  # show python binding calls

oddMaterialMap = geoDir / "data/odd-material-maps.root"
oddDigiConfig = geoDir / "config/odd-digi-smearing-config.json"
oddSeedingSel = geoDir / "config/odd-seeding-config.json"
oddMaterialDeco = acts.IMaterialDecorator.fromFile(oddMaterialMap) #instead of simulating the entire detector, there are predefined surfaces to which we map the material -- allows to compute energy loss

detector, trackingGeometry, decorators = getOpenDataDetector( #we don´t know detector, only describe it by material mapping, tracking geometry is built bu acts modules on the other hand, decorator is what defines material mapping
    geoDir, mdecorator=oddMaterialDeco
) 

field = acts.ConstantBField(acts.Vector3(0.0, 0.0, 2.0 * u.T))
rnd = acts.examples.RandomNumbers(seed=42) #random number generator

# TODO Geant4 currently crashes with FPE monitoring
with acts.FpeMonitor() if not g4_simulation else contextlib.nullcontext(): #fpe = floating point exception, to find position in code where fpe occurs
    s = acts.examples.Sequencer( #sequencer is the entire workflow, avoids python loops to improve performance
        events=args["events"],
        numThreads=1,
        outputDir=str(outputDir),
    )

    #decide initial conditions for particles
    #difference between addPythia and addParticleGun?
    if not ttbar_pu200:
        addParticleGun(
            s,
            MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True), #p_T uniform zw 1 und 10 GeV
            EtaConfig(-3.0, 3.0, uniform=True), #p_z aus eta
            ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True), #generate (anti)muon pairs, direction are uncorrelated, but same vertex, allows vertexing
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ), #vertices are generated from a multivariate Gaussian distribution with 0 mean, time coordinate is only built-in for future use, its value should not impact outcome
            multiplicity=50, #in total, we will create multiplicity*n*2 particles, where 2 is because we create pairs of (anti)particles and n is the number of events
            rnd=rnd, 
        )
    else:
        addPythia8(
            s,
            hardProcess=["Top:qqbar2ttbar=on"], #top antitop production
            npileup=200, #(avg?) number of pairs per event
            vtxGen=acts.examples.GaussianVertexGenerator(
                stddev=acts.Vector4(
                    0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns
                ),
                mean=acts.Vector4(0, 0, 0, 0),
            ),
            rnd=rnd,
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
        )
    #decide how to simulate particle movement
    #Fatras = fast track simulation, acts software, not as elaborate (e.g., no secondary particles)
    if g4_simulation:
        if s.config.numThreads != 1:
            raise ValueError("Geant 4 simulation does not support multi-threading")

        # Pythia can sometime simulate particles outside the world volume, a cut on the Z of the track help mitigate this effect
        # Older version of G4 might not work, this as has been tested on version `geant4-11-00-patch-03`
        # For more detail see issue #1578
        addGeant4(
            s,
            detector,
            trackingGeometry,
            field,
            preSelectParticles=ParticleSelectorConfig(
                eta=(-3.0, 3.0),
                absZ=(0, 1e4),
                rho=(0, 1e3),
                pt=(150 * u.MeV, None),
                removeNeutral=True,
            ), #throw away particles generated by Pythia whose inital conditions are not in the ranges specified above
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            rnd=rnd, #why do we need randomness in the simulation?
        )
    else:
        addFatras(
            s,
            trackingGeometry,
            field,
            preSelectParticles=ParticleSelectorConfig(
                eta=(-3.0, 3.0),
                pt=(150 * u.MeV, None),
                removeNeutral=True,
            )
            if ttbar_pu200
            else ParticleSelectorConfig(),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            rnd=rnd,
        )
    #to blurr particle tracks (finite resolution of position due to finite cell size of the tracker)
    addDigitization(
        s,
        trackingGeometry,
        field,
        digiConfigFile=oddDigiConfig,
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
        rnd=rnd,
    )
    #seeding = take part of tracking system and try connecting triplets, look which initial conditions would explain these triplets and apply cuts in momentum and initial position
    # -> gives track candidates
    #Truthseeding (flag not activated here): only take three best hits whose initial conditions are in a certain range, throw away other hits, and then start seeding
    addSeeding(
        s,
        trackingGeometry,
        field,
        TruthSeedRanges(pt=(1.0 * u.GeV, None), eta=(-3.0, 3.0), nHits=(9, None)) #conditions for truth seeding (not applied here)
        if ttbar_pu200
        else TruthSeedRanges(), 
        geoSelectionConfigFile=oddSeedingSel, #specifies subpart of tracking system where we perform seeding
        outputDirRoot=outputDir,
    )
    #CKF = Combinatorial Kalman Filter
    #from seeding we get rought esitimate of momentum start at nearest hit to beam pipe and start propagating, try out which hits fit the best combinatorically and sequetially update
    #momentum
    addCKFTracks(
        s,
        trackingGeometry,
        field,
        CKFPerformanceConfig(
            ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0,
            nMeasurementsMin=7,
        ), #conditions for performance output file
        TrackSelectorRanges(
            pt=(1.0 * u.GeV, None),
            absEta=(None, 3.0),
            loc0=(-4.0 * u.mm, 4.0 * u.mm),
        ), #conditions for keeping the tracks
        outputDirRoot=outputDir,
        # outputDirCsv=outputDir,
    )

    #during seeding and CKF several tracks corresponding to the same particle are found - these are removed here
    if ambiguity_MLSolver:
        addAmbiguityResolutionML(
            s,
            AmbiguityResolutionMLConfig(nMeasurementsMin=7),
            CKFPerformanceConfig(
                ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0, nMeasurementsMin=7
            ),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
            onnxModelFile=os.path.dirname(__file__)
            + "/MLAmbiguityResolution/duplicateClassifier.onnx",
        )
    else:
        addAmbiguityResolution(
            s,
            AmbiguityResolutionConfig(
                maximumSharedHits=3, maximumIterations=10000, nMeasurementsMin=7
            ),
            CKFPerformanceConfig(
                ptMin=1.0 * u.GeV if ttbar_pu200 else 0.0,
                nMeasurementsMin=7,
            ),
            outputDirRoot=outputDir,
            # outputDirCsv=outputDir,
        )
    
    #clustering of tracks that belong to the same vertex
    addVertexFitting(
        s,
        field,
        vertexFinder=VertexFinder.Iterative,
        outputDirRoot=outputDir,
    )

    s.run()
