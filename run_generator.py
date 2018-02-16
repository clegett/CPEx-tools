import mmm
import math
import os
import errno
import sys

RHOST = 500
DRRIM = 100
RINCL = 10
NINCL = 100
mypack = mmm.Pack.from_file('example_input/10.dat')
spheres_in_pack = 10
ocfile = 'example_input/rerun_oc.csv'

queue = 'devel'
cpumodel = 'bro'
nodes = 100
tpn = 28
walltime = '1:59:00'

runname = str(spheres_in_pack) + 'x' + str(RHOST) + '-' + str(DRRIM) + '-' \
          + str(RINCL) + '-' + str(NINCL)
runoutputdir = ('/nobackupp8/clegett/' + runname + '/')
outdir = runname + '/'

try:
    with open(ocfile) as f:
        oclines = f.read().splitlines()
except IOError as e:
    sys.exit('I/O error: file {}: {}'.format(ocfile, e))

ocarray = []
for line in oclines:
    ocarray.append(line.split(','))
ocarray.pop(0) 

mypack.rescale_pack(RHOST + DRRIM)
mypack.center_pack()
mycluster = mmm.Cluster()
mycluster.grainlist = []

for sphere in mypack.sphere_coords:
    mycluster.grainlist.append(mmm.Grain.new_grain(sphere[0], sphere[1],
                               sphere[2], RHOST, DRRIM, RINCL, NINCL))

myrun = mmm.ModelRun(name=runname, fixed_or_random=mmm.RunType.FIXED)
myrun.set_option('number_spheres', (1 + 1 + NINCL) * len(mycluster.grainlist))
myrun.set_option('max_number_iterations', 5000)

try:
    os.makedirs(outdir)
except OSError as e:
    if e.errno != errno.EEXIST:
        raise

outfile = outdir + runname + '.inp'
try:
    with open(outfile, 'w') as out:
        for line in ocarray:
            for grain in mycluster.grainlist:
                grain.set_grain_oc(line[1], line[2], line[3], line[4],
                                   line[5], line[6])

            myrun.set_option('length_scale_factor', (2*math.pi/float(line[0])))
            myrun.set_option('output_file', runoutputdir + 'dat/' + line[0]
                             + 'nm.dat')
            myrun.set_option('scattering_coefficient_file',
                             (runoutputdir + 'sc/' + line[0] + 'nm.sc.dat'))
            myrun.set_option('run_print_file', runoutputdir + 'run_print.dat')
            myrun.set_option('azimuth_average_scattering_matrix', 1)
            myrun.set_option('delta_scattering_angle_deg', 1)

            runoutput = myrun.formatted_options()
            print(runoutput, file=out)
            print('sphere_sizes_and_positions', file=out)
            for grain in mycluster.grainlist:
                grain.print_rxyznk(out)
            if line[0] is not ocarray[len(ocarray) - 1][0]:
                print('new_run', file=out)
except IOError as e:
    sys.exit('I/O error: file {}: {}'.format(outfile, e))

pbsfile = outdir + runname + '.pbs'
try:
    with open(pbsfile, 'w') as pbs:
        print('#PBS -N ' + runname, file=pbs)
        print('#PBS -q ' + queue, file=pbs)
        print('#PBS -l select=' + str(nodes) + ':ncpus=' + str(tpn) + 
              ':mpiprocs=' + str(tpn) + ':model=' + cpumodel, file=pbs)
        print('#PBS -l walltime=' + walltime, file=pbs)
        print('#PBS -e ' + runoutputdir + runname + '.err', file=pbs)
        print('#PBS -o ' + runoutputdir + runname + '.out', file=pbs)
        print('#PBS -M carey.legett@stonybrook.edu', file=pbs)
        print('#PBS -m abe', file=pbs)
        print('', file=pbs)
        print('module load comp-intel/2016.2.181 mpi-sgi/mpt', file=pbs)
        print('', file=pbs)
        print('mkdir ' + runoutputdir, file=pbs)
        print('mkdir ' + runoutputdir + 'dat', file=pbs)
        print('mkdir ' + runoutputdir + 'sc', file=pbs)
        print('', file=pbs)
        print('cd $PBS_O_WORKDIR', file=pbs)
        print('', file=pbs)
        print('mpiexec -np ' + str(nodes * tpn) + ' ./mstm_ttv2.3.exe ' +
              runname + '.inp', file=pbs)
        print('mv ~/' + runname + '* ' + runoutputdir, file=pbs)
except IOError as e:
    sys.exit('I/O error: file {}: {}'.format(pbsfile, e))

try:
    with open(outdir + runname + 'stats.txt', 'w') as stats:
        print('runname: ' + runname, file=stats)
        print('spheres in pack: ' + str(spheres_in_pack), file=stats)
        print('cluster packing fraction: '
              + str(mycluster.get_packing_fraction()), file=stats)
        print('cluster x-section: ' + str(mycluster.get_geom_xsection()),
              file=stats)
        print('cluster bounding sphere radius: '
              + str(mycluster.get_bounding_sphere().r), file=stats)
except IOError as e:
    sys.exit('IOError while writing {}: {}'.format((outdir + runname
                                                    + 'stats.txt'), e))
