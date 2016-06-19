
epicsEnvSet("EPICS_DB_INCLUDE_PATH", "../vaccel:.")

dbLoadDatabase("flamedemo.dbd")
flamedemo_registerRecordDeviceDriver(pdbbase)

< vconfig

flamePrepare("THESIM", "$(LATTICE=input.lat)")
dbLoadRecords("core.db", "P=TST:,SIM=THESIM")
dbLoadRecords("orbit.db", "P=TST:,SIM=THESIM,NELM=%(nmarkers)d")
dbLoadTemplate("lattice.substitutions", "PREF=TST:,SIM=THESIM")

iocInit()
