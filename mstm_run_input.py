#!/usr/bin/env python3
"""Control class for the MSTM Model Manager

"""

import sys
import numpy as np

__author__ = 'Carey Legett'
__copyright__ = 'Copyright 2018, Stony Brook University'
__credits__ = ['Carey Legett']

__version__ = '1.0'
__maintainer__ = 'Carey Legett'
__email__ = 'carey.legett@stonybrook.edu'
__status__ = 'Production'

class Pack:
    """A class containing the data generated by a PackLSD run.

    This class holds the metadata and x,y,z coordinates generated by a PackLSD
    run.

    Attributes:


    """
    def __init__(self, dims=None, num_of_particles=None, dispersity=None,
                 sphere_radius=None, sphere_coords=[[None, None, None]]):
        self.dims = dims
        self.num_of_particles = num_of_particles
        self.dispersity = dispersity
        self.sphere_radius = sphere_radius
        self.sphere_coords = sphere_coords

    @classmethod
    def from_file(cls, a_filename):

        # Check for the version of python being used and use appropriate flags
        # for opening the input file as necessary
        if sys.version_info[0] == 2:
            raccess = 'rb'
            kwargs = {}
        else:
            raccess = 'rt'
            kwargs = {'newline':''}

        try:
            with open(a_filename, raccess, **kwargs) as infile:
                rawfile = infile.read().splitlines()
            if not rawfile:
                sys.exit('No data in file ' + a_filename)
        except IOError as e:
            sys.exit('I/O error: file {}: {}'.format(a_filename, e))
        except:
            sys.exit('Unexpected error: {}'.format(sys.exc_info()))

        lines = []
        for line in rawfile:
            lines.append(line.split())

        return(cls(lines[0][0], lines[1][0], lines[1][1]), lines[3][0],
                   lines[6:])

    def center_pack(self):
        xtotal = 0
        ytotal = 0
        ztotal = 0
        for line in self.sphere_coords:
            xtotal += line[0]
            ytotal += line[1]
            ztotal += line[2]
        length = len(self.sphere_coords)
        avgx = xtotal/length
        avgy = ytotal/length
        avgz = ztotal/length
        for line in self.sphere_coords:
            line[0] = line[0] - avgx
            line[1] = line[1] - avgy
            line[2] = line[2] - avgz

    def rescale_pack(self, new_sphere_radius):
        multiplier = new_sphere_radius / self.sphere_radius
        self.sphere_coords = [[coord * multiplier for coord in row] for row in
            self.sphere_coords]