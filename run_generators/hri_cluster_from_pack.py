#!/usr/bin/env python3
"""This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
"""

import mmm
import math
import os
import errno
import sys

__author__ = 'Carey Legett'
__contact__ = 'carey.legett@stonybrook.edu'
__copyright__ = 'Copyright 2018, Stony Brook University'
__credits__ = ['Carey Legett', 'Freenode #r/linux', 'Freenode #python']
__date__ = '2018/02/19'
__deprecated__ = False
__email__ = 'carey.legett@stonybrook.edu'
__license__ = 'GPLv3'
__maintainer = 'Carey Legett'
__status__ = 'Development'
__version__ = '1.0'


def main():
    """
    INPUT VARIABLES
    ############################################################################
    Edit this section to change the parameters of the MSTM run
    """
    rhost = 500         # Radius of the host particle
    drrim = 100         # Thickness of the rim
    rincl = 10          # Radius of the inclusions in the rim
    nincl = 100         # Number of inclusions in the rim
    mypack = mmm.Pack.from_file('example_input/10.dat')  # Location of pack file
    spheres_in_pack = 10    # Number of spheres in the pack
    ocfile = 'example_input/rerun_oc.csv'   # Optical constants file

    queue = 'devel'     # Pleiades queue name
    cpumodel = 'bro'    # Pleiades processor model
    nodes = 100         # Pleiades number of nodes to request
    tpn = 28            # Pleiades threads per node
    walltime = '1:59:00'  # Pleiades walltime to request

    """
    End of input variables section. You probably do NOT want to edit anything 
    below here.
    ############################################################################
    """

    runname = '{!s}x{!s}-{!s}-{!s}-{!s}'.format(spheres_in_pack, rhost, drrim,
                                                rincl, nincl)
    runoutputdir = '/nobackupp8/clegett/{}/'.format(runname)
    outdir = '{}/'.format(runname)

    try:
        with open(ocfile) as f:
            oclines = f.read().splitlines()
    except IOError as e:
        sys.exit('I/O error: file {}: {}'.format(ocfile, e))

    ocarray = []
    for line in oclines:
        ocarray.append(line.split(','))
    ocarray.pop(0)

    mypack.rescale_pack(rhost + drrim)
    mypack.center_pack()
    mycluster = mmm.Cluster()
    mycluster.grainlist = []

    for sphere in mypack.sphere_coords:
        mycluster.grainlist.append(mmm.Grain.new_grain(sphere[0], sphere[1],
                                                       sphere[2], rhost, drrim,
                                                       rincl, nincl))

    myrun = mmm.ModelRun(name=runname, fixed_or_random=mmm.RunType.FIXED)
    myrun.set_option('number_spheres', (1 + 1 + nincl) *
                     len(mycluster.grainlist))
    myrun.set_option('max_number_iterations', 5000)

    try:
        os.makedirs(outdir)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

    outfile = '{}{}.inp'.format(outdir, runname)
    try:
        with open(outfile, 'w') as out:
            for line in ocarray:
                for grain in mycluster.grainlist:
                    grain.set_grain_oc(line[1], line[2], line[3], line[4],
                                       line[5], line[6])

                myrun.set_option('length_scale_factor',
                                 (2 * math.pi / float(line[0])))
                myrun.set_option('output_file',
                                 '{}dat/{}.nm.dat'.format(runoutputdir,
                                                          line[0]))
                myrun.set_option('scattering_coefficient_file',
                                 '{}sc/{}.nm.sc.dat'.format(runoutputdir,
                                                            line[0]))
                myrun.set_option('run_print_file',
                                 '{}run_print.dat'.format(runoutputdir))
                myrun.set_option('azimuth_average_scattering_matrix', 1)
                myrun.set_option('delta_scattering_angle_deg', 1)

                runoutput = myrun.formatted_options()
                print(runoutput, file=out)
                print('sphere_sizes_and_positions', file=out)
                for grain in mycluster.grainlist:
                    print(grain.get_rxyznk(), file=out)
                if line[0] is not ocarray[len(ocarray) - 1][0]:
                    print('new_run', file=out)
    except IOError as e:
        sys.exit('I/O error: file {}: {}'.format(outfile, e))

    pbsfile = '{}{}.pbs'.format(outdir, runname)
    try:
        with open(pbsfile, 'w') as pbs:
            print('#PBS -N {}'.format(runname), file=pbs)
            print('#PBS -q {}'.format(queue), file=pbs)
            print('#PBS -l select={!s}:ncpus={!s}:mpiprocs={!s}:model={}'
                  ''.format(nodes, tpn, tpn, cpumodel), file=pbs)
            print('#PBS -l walltime={}'.format(walltime), file=pbs)
            print('#PBS -e {}{}.err'.format(runoutputdir, runname), file=pbs)
            print('#PBS -o {}{}.out'.format(runoutputdir, runname), file=pbs)
            print('#PBS -M carey.legett@stonybrook.edu', file=pbs)
            print('#PBS -m abe', file=pbs)
            print('', file=pbs)
            print('module load comp-intel/2016.2.181 mpi-sgi/mpt', file=pbs)
            print('', file=pbs)
            print('mkdir {}'.format(runoutputdir), file=pbs)
            print('mkdir {}dat'.format(runoutputdir), file=pbs)
            print('mkdir {}sc'.format(runoutputdir), file=pbs)
            print('', file=pbs)
            print('cd $PBS_O_WORKDIR', file=pbs)
            print('', file=pbs)
            print('mpiexec -np {!s} ./mstm_ttv2.3.exe {}.inp'
                  ''.format((nodes * tpn), runname), file=pbs)
            print('mv ~/{}* {}'.format(runname, runoutputdir), file=pbs)
    except IOError as e:
        sys.exit('I/O error: file {}: {}'.format(pbsfile, e))

    try:
        with open('{}{}stats.txt'.format(outdir, runname), 'w') as stats:
            print('runname: {}'.format(runname), file=stats)
            print('spheres in pack: {!s}'.format(spheres_in_pack), file=stats)
            print('cluster packing fraction: {!s}'
                  ''.format(mycluster.get_packing_fraction()), file=stats)
            print('cluster x-section: {!s}'.format(
                mycluster.get_geom_xsection()), file=stats)
            print('cluster bounding sphere radius: {}'
                  ''.format(mycluster.get_bounding_sphere().r), file=stats)
    except IOError as e:
        sys.exit('IOError while writing {}: {}'.format(
            '{}{}stats.txt'.format(outdir, runname), e))


if __name__ == '__main__':
    main()
