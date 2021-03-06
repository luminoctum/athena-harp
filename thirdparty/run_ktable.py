#! /usr/bin/env python2.7
from multiprocessing import Pool
import os, re, subprocess

# file architecture
hitbin  = "./hitbin"
rfm     = "./rfm"
hitran  = "./HITRAN2012.par"
hitfile = "./12.5-10000.hit"
run_rfm = "./run_rfm.py"
kcoeff  = "./kcoeff"

# ktable specifics
atm     = "new_jup_lat30_modelc.txt"
mol     = "CH4 C2H2 C2H4 C2H6"

# temperature grid
temp    = "-5 5 3"

# flags
generate_tab = True
generate_nc  = True
cia     = "H2-H2-eq H2-He-eq"
dist    = "5.03"

# spectral bands 
band= []
band.append("10 200 0.01")
band.append("200 400 0.01")
band.append("400 600 0.01")
band.append("600 700 0.01")
band.append("700 860 0.01")
band.append("860 960 0.01")
band.append("960 1200 0.01")
band.append("1200 1400 0.01")
band.append("1400 1700 0.01")
band.append("1700 2100 0.01")
band.append("2100 3450 0.01")
band.append("3450 4800 0.01")
band.append("4800 6300 0.01")
band.append("6300 7800 0.01")
band.append("7800 9300 0.01")

# number of parallel threads
nthreads = len(band)

# run ktable in single thread
def RunSingleKtable(wave):
  wname = re.sub(' ', '-', wave)
  print "working on band %s ..." % wname
  wmin, wmax, dw = map(float, wave.split())
  tab_folder = str(wmin) + "-" + str(wmax)
  kinp = 'kcoeff.inp-'+wname
  kncfile = "kcoeff." + wname +".nc"
  # create tab files and kcoeff.inp
  if generate_tab:
    script = [run_rfm,'--hitbin',hitbin,'--rfm',rfm,'--par',hitran,'--hit',hitfile,
                      '--atm',atm,'--wave',wave,'--temp',temp,
                      '--molecule',mol,'--output',kinp,'--rundir',tab_folder]
    out, err = subprocess.Popen(script,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE).communicate()
    print out, err
  # run kcoeff
  if generate_nc:
    script = [kcoeff,'-i',kinp,'-o','kcoeff.%s.nc' % wname]
    out, err = subprocess.Popen(script,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE).communicate()
    print out, err
  print "band %s finishes." % wave

# parallel on spectral bands
pool = Pool(nthreads)
pool.map(RunSingleKtable, band)
#RunSingleKtable(band[0])
