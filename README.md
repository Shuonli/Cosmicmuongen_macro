# Cosmic muon generator analysis macros

CondorRun.sh and run_condor.sh are used to run the cosmic muon generator on condor, several millions of generated events are needed to make all towers to have a well-defined MIP peak.

G4_HcalOut_ref.C is modified so that the simulated ADC output is corresponding to the high gain mode.

G4Setup_sPHENIX.C is modified to increase the boundary size.

oHCalsimulationanalysis.C and oHCalsimulationanalysisalltower.C takes the TTree output and convert that to the MIP energy histogram and preform the rate calculation. 
