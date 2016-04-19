
dbLoadDatabase("flame.dbd")
flame_registerRecordDeviceDriver(pdbbase)

flamePrepare("THESIM", "input.lat")
dbLoadRecords("core.db", "P=TST:,SIM=THESIM")
dbLoadRecords("orbit.db", "P=TST:,SIM=THESIM,NELM=%(nmarkers)d")
dbLoadTemplate("lattice.substitutions", "PREF=TST:,SIM=THESIM")

iocInit()
