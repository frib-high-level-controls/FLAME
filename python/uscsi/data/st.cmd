
epicsEnvSet("EPICS_DB_INCLUDE_PATH", "../vaccel:.")

dbLoadDatabase("flame.dbd")
flame_registerRecordDeviceDriver(pdbbase)

< vconfig

flamePrepare("THESIM", "$(LATTICE=input.lat)")
dbLoadRecords("core.db", "P=TST:,SIM=THESIM")
dbLoadRecords("orbit.db", "P=TST:,SIM=THESIM,NELM=%(nmarkers)d")
dbLoadTemplate("lattice.substitutions", "PREF=TST:,SIM=THESIM")

iocInit()
