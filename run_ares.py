import ares 

#This file should be able to read parameters from either command line or a file.
#Then at the end it should output the dTb(z) to some file which can be read in by
# the c++ code.

sim = ares.simulations.Global21cm(verbose = False, ThreadNum = 4)
sim.run()
