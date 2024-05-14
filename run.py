import os
import subprocess
import numpy as np
import glob

compile = "gfortran -ffree-line-length-none Step1_rz.f90"
subprocess.run(compile, shell=True)
subprocess.run(['./a.out'], shell=True)


flist = glob.glob('*.csv')
for fname in flist:
    data = np.genfromtxt(fname)
    np.savetxt(fname, data, delimiter=',', fmt='%.8f')
