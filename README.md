# High-performance Atmospheric Radiation Package (HARP)

## How to make a line by line radiative transfer calculation

1. git clone the source code
```
git clone https://github.com/luminoctum/athena-harp harp-1.0
```


2. cd into third party software folder
```
cd harp-1.0/thirdparty
```

3. untar RFM and DISORT softwares
```
tar -xvf rfm-4.33.tar.gz
tar -xvf cdisort-2.1.3a.tar.gz
```

4. compile third party softwares
```
make
```

5. make a experiment folder
```
mkdir 1d-rad-jupiter
```

6. create or copy model input file
```
cp input/athinput.radjup ./1d-rad-jupiter/
```

7. copy the configuration of the input file and configure
```
head 1d-rad-jupiter/athinput.radjup
./configure --comp=a5 --prob=radconv --eos=heterogeneous --main=radiation -netcdf -disort
make -j8
```

8. cd to experiment folder
```
cd 1d-rad-jupiter
```

9. link excutables
```
ln -s ../bin/athena ./
ln -s YOUR_HITRAN_PARFILE ./
ln -s ../thirdparty/hitbin ./
cp ../thirdparty/run_ktable.csh ./
ln -s ../thirdparty/run_rfm.py ./
ln -s ../thirdparty/kcoeff ./
```

10. run hibtin
```
./hitbin
```

11. run ktable
```
./run_ktable.csh
```

12. run radiative transfer code
```
./athena -i athinput
```
