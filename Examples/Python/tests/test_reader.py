import pytest
from pathlib import Path

from helpers import (
    geant4Enabled,
    AssertCollectionExistsAlg,
)

import acts
from acts import UnitConstants as u
from acts.examples import (
    RootParticleWriter,
    RootParticleReader,
    RootMaterialTrackReader,
    RootTrackSummaryReader,
    CsvParticleWriter,
    CsvParticleReader,
    CsvMeasurementWriter,
    CsvMeasurementReader,
    CsvSimHitWriter,
    CsvSimHitReader,
    Sequencer,
)
from acts.examples.odd import getOpenDataDetector, getOpenDataDetectorDirectory


@pytest.mark.root
def test_root_particle_reader(tmp_path, conf_const, ptcl_gun):
    # need to write out some particles first
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    _, h3conv = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.WARNING,
            inputParticles=h3conv.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    # reset sequencer for reading

    s2 = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s2.addReader(
        conf_const(
            RootParticleReader,
            acts.logging.WARNING,
            outputParticles="particles_generated",
            filePath=str(file),
        )
    )

    alg = AssertCollectionExistsAlg(
        "particles_generated", "check_alg", acts.logging.WARNING
    )
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_particle_reader(tmp_path, conf_const, ptcl_gun):
    s = Sequencer(numThreads=1, events=10, logLevel=acts.logging.WARNING)
    _, h3conv = ptcl_gun(s)

    out = tmp_path / "csv"

    out.mkdir()

    s.addWriter(
        conf_const(
            CsvParticleWriter,
            acts.logging.WARNING,
            inputParticles=h3conv.config.outputParticles,
            outputStem="particle",
            outputDir=str(out),
        )
    )

    s.run()

    # reset the seeder
    s = Sequencer(numThreads=1, logLevel=acts.logging.WARNING)

    s.addReader(
        conf_const(
            CsvParticleReader,
            acts.logging.WARNING,
            inputDir=str(out),
            inputStem="particle",
            outputParticles="input_particles",
        )
    )

    alg = AssertCollectionExistsAlg(
        "input_particles", "check_alg", acts.logging.WARNING
    )

    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.parametrize(
    "reader",
    [RootParticleReader, RootTrackSummaryReader],
)
@pytest.mark.root
def test_root_reader_interface(reader, conf_const, tmp_path):
    assert hasattr(reader, "Config")

    config = reader.Config

    assert hasattr(config, "filePath")

    kw = {"level": acts.logging.INFO, "filePath": str(tmp_path / "file.root")}

    assert conf_const(reader, **kw)


@pytest.mark.slow
@pytest.mark.root
@pytest.mark.odd
@pytest.mark.skipif(not geant4Enabled, reason="Geant4 not set up")
def test_root_material_track_reader(material_recording):
    input_tracks = material_recording / "geant4_material_tracks.root"
    assert input_tracks.exists()

    s = Sequencer(numThreads=1)

    s.addReader(
        RootMaterialTrackReader(
            level=acts.logging.INFO,
            fileList=[str(input_tracks)],
            outputMaterialTracks="material-tracks",
        )
    )

    alg = AssertCollectionExistsAlg(
        "material-tracks", "check_alg", acts.logging.WARNING
    )
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 2


@pytest.mark.csv
def test_csv_meas_reader(tmp_path, fatras, trk_geo, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        CsvMeasurementWriter(
            level=acts.logging.INFO,
            inputMeasurements=digiAlg.config.outputMeasurements,
            inputClusters=digiAlg.config.outputClusters,
            inputMeasurementSimHitsMap=digiAlg.config.outputMeasurementSimHitsMap,
            outputDir=str(out),
        )
    )

    # Write hits, so we can later construct the measurement-particles-map
    s.addWriter(
        CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputDir=str(out),
            outputStem="hits",
        )
    )

    s.run()

    # read back in
    s = Sequencer(numThreads=1)

    s.addReader(
        CsvSimHitReader(
            level=acts.logging.INFO,
            outputSimHits=simAlg.config.outputSimHits,
            inputDir=str(out),
            inputStem="hits",
        )
    )

    s.addReader(
        conf_const(
            CsvMeasurementReader,
            level=acts.logging.WARNING,
            outputMeasurements="measurements",
            outputMeasurementSimHitsMap="simhitsmap",
            outputMeasurementParticlesMap="meas_ptcl_map",
            inputSimHits=simAlg.config.outputSimHits,
            inputDir=str(out),
        )
    )

    algs = [
        AssertCollectionExistsAlg(k, f"check_alg_{k}", acts.logging.WARNING)
        for k in ("measurements", "simhitsmap", "meas_ptcl_map")
    ]
    for alg in algs:
        s.addAlgorithm(alg)

    s.run()

    for alg in algs:
        assert alg.events_seen == 10


@pytest.mark.csv
def test_csv_simhits_reader(tmp_path, fatras, conf_const):
    s = Sequencer(numThreads=1, events=10)
    evGen, simAlg, digiAlg = fatras(s)

    out = tmp_path / "csv"
    out.mkdir()

    s.addWriter(
        CsvSimHitWriter(
            level=acts.logging.INFO,
            inputSimHits=simAlg.config.outputSimHits,
            outputDir=str(out),
            outputStem="hits",
        )
    )

    s.run()

    s = Sequencer(numThreads=1)

    s.addReader(
        conf_const(
            CsvSimHitReader,
            level=acts.logging.INFO,
            inputDir=str(out),
            inputStem="hits",
            outputSimHits="simhits",
        )
    )

    alg = AssertCollectionExistsAlg("simhits", "check_alg", acts.logging.WARNING)
    s.addAlgorithm(alg)

    s.run()

    assert alg.events_seen == 10


@pytest.mark.root
def test_buffered_reader(tmp_path, conf_const, ptcl_gun):
    # Test the buffered reader with the ROOT particle reader
    # need to write out some particles first
    eventsInBuffer = 5
    eventsToProcess = 10

    s = Sequencer(numThreads=1, events=eventsInBuffer, logLevel=acts.logging.WARNING)
    _, h3conv = ptcl_gun(s)

    file = tmp_path / "particles.root"
    s.addWriter(
        conf_const(
            RootParticleWriter,
            acts.logging.WARNING,
            inputParticles=h3conv.config.outputParticles,
            filePath=str(file),
        )
    )

    s.run()

    # reset sequencer for reading
    s2 = Sequencer(events=eventsToProcess, numThreads=1, logLevel=acts.logging.WARNING)

    reader = acts.examples.RootParticleReader(
        level=acts.logging.WARNING,
        outputParticles="particles_input",
        filePath=str(file),
    )

    s2.addReader(
        acts.examples.BufferedReader(
            level=acts.logging.WARNING,
            upstreamReader=reader,
            bufferSize=eventsInBuffer,
        )
    )

    alg = AssertCollectionExistsAlg(
        "particles_input", "check_alg", acts.logging.WARNING
    )
    s2.addAlgorithm(alg)

    s2.run()

    assert alg.events_seen == eventsToProcess
