<comment>
problem   = One layer shallow water model
reference = 
configure = --prob=sw --eos=shallow_water -netcdf -mpi

<job>
problem_id = sw        # problem ID: basename of output filenames

<output1>
file_type  = hst       # History data dump
dt         = 0.1      # time increment between outputs

<output2>
file_type  = netcdf    # netcdf data dump
variable   = prim      # variables to be output
dt         = 0.1       # time increment between outputs

<time>
cfl_number = 0.4       # The Courant, Friedrichs, & Lewy (CFL) Number
nlim       = 100000    # cycle limit
tlim       = 40.0      # time limit

<mesh>
nx1        = 256       # Number of zones in X1-direction
x1min      = -10.0      # minimum value of X1
x1max      =  10.0      # maximum value of X1
ix1_bc     = periodic  # inner-X1 boundary flag
ox1_bc     = periodic  # inner-X1 boundary flag

nx2        = 128       # Number of zones in X2-direction
x2min      = -5.     # minimum value of X2
x2max      = 5.      # maximum value of X2
ix2_bc     = periodic  # inner-X2 boundary flag
ox2_bc     = periodic  # inner-X2 boundary flag

nx3        = 1         # Number of zones in X3-direction
x3min      = -0.5      # minimum value of X3
x3max      = 0.5       # maximum value of X3
ix3_bc     = periodic # inner-X3 boundary flag
ox3_bc     = periodic # inner-X3 boundary flag

<meshblock>
nx1       = 64
nx2       = 64
nx3       = 1

<problem>
iprob = 1
amp   = 0.01
drat  = 2.0
vflow = 0.5
