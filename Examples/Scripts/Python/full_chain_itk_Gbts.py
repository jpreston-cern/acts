#!/usr/bin/env python3
import pathlib, acts, acts.examples, acts.examples.itk
from acts.examples.simulation import (#import all the config and functions we need from reconstruction.py
    addParticleGun,
    MomentumConfig,
    EtaConfig,
    ParticleConfig,
    addPythia8,
    ParticleSelectorConfig,
    addGenParticleSelection,
    addFatras,
    addDigitization,
    addDigiParticleSelection,
)
from acts.examples.reconstruction import (#same here 
    addSeeding, #important function 
    SeedingAlgorithm,
    addCKFTracks,
    TrackSelectorConfig,
)
#settings we can change for whne we run 
ttbar_pu200 = False
u = acts.UnitConstants #units used in acts 
geo_dir = pathlib.Path("acts-itk") #itk files 
outputDir = pathlib.Path.cwd() / "itk_output"
# acts.examples.dump_args_calls(locals())  # show acts.examples python binding calls

detector = acts.examples.itk.buildITkGeometry(geo_dir) #heavy itk object that holds everything
trackingGeometry = detector.trackingGeometry()# the actual volume 
field = acts.examples.MagneticFieldMapXyz(str(geo_dir / "bfield/ATLAS-BField-xyz.root"))
rnd = acts.examples.RandomNumbers(seed=42)

s = acts.examples.Sequencer(events=1, numThreads=1, outputDir=str(outputDir))

if not ttbar_pu200:
    addParticleGun(
        s,
        MomentumConfig(1.0 * u.GeV, 10.0 * u.GeV, transverse=True),
        EtaConfig(-4.0, 4.0, uniform=True),
        ParticleConfig(2, acts.PdgParticle.eMuon, randomizeCharge=True),
        rnd=rnd,
    )
else:
    addPythia8(
        s,
        hardProcess=["Top:qqbar2ttbar=on"],
        npileup=200,
        vtxGen=acts.examples.GaussianVertexGenerator(
            stddev=acts.Vector4(0.0125 * u.mm, 0.0125 * u.mm, 55.5 * u.mm, 5.0 * u.ns),
            mean=acts.Vector4(0, 0, 0, 0),
        ),
        rnd=rnd,
        outputDirRoot=outputDir,
    )

    addGenParticleSelection(
        s,
        ParticleSelectorConfig(
            rho=(0.0 * u.mm, 28.0 * u.mm),
            absZ=(0.0 * u.mm, 1.0 * u.m),
            eta=(-4.0, 4.0),
            pt=(150 * u.MeV, None),
        ),
    )

addFatras(
    s,
    trackingGeometry,
    field,
    rnd=rnd,
    outputDirRoot=outputDir,
)

addDigitization(
    s,
    trackingGeometry,
    field,
    digiConfigFile=geo_dir
    / "itk-hgtd/itk-smearing-config.json",  # change this file to make it do digitization
    outputDirRoot=outputDir,
    rnd=rnd,
)

addDigiParticleSelection(
    s,
    ParticleSelectorConfig(
        pt=(1.0 * u.GeV, None),
        eta=(-4.0, 4.0),
        measurements=(9, None),
        removeNeutral=True,
    ),
)

addSeeding(
    s, #seequencer (what actually runs everything after making algorithm choices)
    trackingGeometry, #ACTS detector geometry
    field, #magnetic field 
    seedingAlgorithm=SeedingAlgorithm.Gbts, #naming the algorithm  so add seeding can choose it in its definition
    *acts.examples.itk.itkSeedingAlgConfig(#this filled in the itk.py file and initialises and fills python config files (not in c++ yet) using a helper function
        acts.examples.itk.InputSpacePointsType.PixelSpacePoints
    ),
    geoSelectionConfigFile=geo_dir / "itk-hgtd/geoSelection-ITk.json", #specififes which layers are actually to be looked at (means you dont accidently load suppoet structure hits ect)
    layerMappingConfigFile=geo_dir / "itk-hgtd/GbtsMapping.csv", #file for chaning between acts id's to athena/gbts id's
    ConnectorInputConfigFile=geo_dir / "itk-hgtd/GbtsBinTable.txt", #connection table 
    outputDirRoot=outputDir,
)

addCKFTracks(
    s,
    trackingGeometry,
    field,
    TrackSelectorConfig(
        pt=(1.0 * u.GeV if ttbar_pu200 else 0.0, None),
        absEta=(None, 4.0),
        nMeasurementsMin=6,
    ),
    outputDirRoot=outputDir,
)


s.run()
