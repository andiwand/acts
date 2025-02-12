#!/usr/bin/env python3
from pathlib import Path

import acts
import acts.examples
from acts.examples.simulation import (
    addParticleGun,
    addFatras,
    EtaConfig,
    MomentumConfig,
    ParticleConfig,
)

u = acts.UnitConstants


def runFatras(trackingGeometry, field, outputDir, s: acts.examples.Sequencer = None):
    s = s or acts.examples.Sequencer(events=1, numThreads=1)
    s.config.logLevel = acts.logging.INFO
    rnd = acts.examples.RandomNumbers()
    addParticleGun(
        s,
        EtaConfig(-3.0, 3.0),
        MomentumConfig(1 * u.GeV, 1 * u.GeV, transverse=True),
        ParticleConfig(1, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )
    outputDir = Path(outputDir)
    addFatras(
        s,
        trackingGeometry,
        field,
        outputDirCsv=outputDir / "csv",
        outputDirRoot=outputDir,
        rnd=rnd,
    )
    return s


if "__main__" == __name__:
    detector = acts.examples.GenericDetector()
    trackingGeometry = detector.trackingGeometry()
    decorators = detector.contextDecorators()

    field = acts.ConstantBField(acts.Vector3(0, 0, 2 * u.T))

    runFatras(trackingGeometry, field, Path.cwd()).run()
