# -*- coding: utf-8 -*-
""" Classes for MSTM Model Manager's data model

This module provides the classes necessary for the MVC "Model" portion of the
MSTM Model Manager program.

Classes:
    Sphere
    ModelOption
    ModelOptionValue
    RunType
    ModelRun
    Grain
    Cluster
    Pack

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
"""

import math
import attr
import random
import sys
from enum import Enum
import re
import glob
import json
import pickle

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


class Sphere(object):
    """A class describing the spheres used in the MSTM model

    This class holds the physical and optical properties of the spheres used in
    the MSTM model. The physical properties are required variables while the
    optical properties are optional. The units for x, y, z, and r should be
    identical.

    Attributes:
        x (float): the x coordinate for the center of the sphere
        y (float): the y coordinate for the center of the sphere
        z (float): the z coordinate for the center of the sphere
        r (float): the radius of the sphere
        n (float): the real portion of the index of refraction
        k (float): the imaginary portion of the index of refraction (extinction
            coefficient)
        real_chiral (float): the real portion of the chiral factor beta
        imag_chiral (float): the imaginary portion of the chiral factor beta
        tmatrix_file (Str): the name of the file containing the previously
        calculated T-matrix for this sphere

    """
    def __init__(self, x, y, z, r, n=None, k=None, real_chiral=None,
                 imag_chiral=None, tmatrix_file=None):
        """Constructor

        Note: x, y, z, and r must be in identical units

        Args:
            x (float): the x coordinate for the center of the sphere
            y (float): the y coordinate for the center of the sphere
            z (float): the z coordinate for the center of the sphere
            r (float): the radius of the sphere
            n (float, optional): the real portion of the index of refraction
              (default=None)
            k (float, optional): the imaginary portion of the index of
              refraction a.k.a. extinction coefficient (default=None)
            real_chiral (float, optional): the real portion of the chiral
              factor "beta" (default=None)
            imag_chiral (float, optional): the imaginary portion of the chiral
              factor "beta" (Default=None)
            tmatrix_file (Str, optional): the name of the file containing the
              previously calculated T-matrix for this sphere (default=None)

        Returns: a Sphere object

        """
        self.x = x
        self.y = y
        self.z = z
        self.r = r
        self.n = n
        self.k = k
        self.real_chiral = real_chiral
        self.imag_chiral = imag_chiral
        self.tmatrix_file = tmatrix_file

    def distance_from(self, another_sphere=None, x1=None, y1=None, z1=None):
        """Return the distance of the center of this sphere from the provided
        coordinates.

        Args:
            another_sphere -- a second sphere to get the distance from
            x1 -- the x coord of the point from which distance is calculated
            y1 -- the y coord of the point from which distance is calculated
            z1 -- the z coord of the point from which distance is calculated

        Returns:
            A float representing the distance from the provided coords to the
            center of this sphere.

        """
        if another_sphere is None:
            return math.sqrt(
                (self.x - x1) ** 2 + (self.y - y1) ** 2 + (self.z - z1) ** 2)
        else:
            return (math.sqrt((self.x - another_sphere.x) ** 2 + (self.y -
                    another_sphere.y) ** 2 + (self.z - another_sphere.z) ** 2))

    def geom_x_sec(self):
        """Return the geometric cross sectional area of this sphere."""
        return math.pi * self.r ** 2

    def get_tabbed_rxyz(self):
        return '{r!s}\t{x!s}\t{y!s}\t{z!s}'.format(**self.__dict__)

    def get_tabbed_rxyznk(self):
        if self.n is None:
            raise ValueError('attribute n cannot be None')
        elif self.k is None:
            raise ValueError('attribute k cannot be None')
        else:
            return('{r!s}\t{x!s}\t{y!s}\t{z!s}\t{n!s}\t{k!s}'
                   ''.format(**self.__dict__))

    def get_tabbed_chiral(self):
        if self.n is None:
            raise ValueError('attribute n cannot be None')
        elif self.k is None:
            raise ValueError('attribute k cannot be None')
        elif self.real_chiral is None:
            raise ValueError('attribute real_chiral cannot be None')
        elif self.imag_chiral is None:
            raise ValueError('attribute imag_chiral cannot be None')
        else:
            return('{r!s}\t{x!s}\t{y!s}\t{z!s}\t{n!s}\t{k!s}\t'
                   '{real_chiral!s}\t{imag_chiral!s}'.format(**self.__dict__))

    def get_tabbed_tmatrix(self):
        if self.tmatrix_file is None:
            raise ValueError('attribute tmatrix_file cannot be None')
        else:
            return('{r!s}\t{x!s}\t{y!s}\t{z!s}\t{tmatrix_file!s}'
                   ''.format(**self.__dict__))


class ModelOption(object):
    """A class to store a name & default value of an option for the MSTM model.

    This class holds a name and a default value of a single option to be passed
    to the MSTM model in an input file. It is primarily intended to be used
    with the module level dictionary 'opt_dict' which initializes one
    ModelOption for every possible valid option.

    Attributes:
        name (str): the name of the MSTM model option
        default_value: the default value of the model option as listed in MSTM
            documentation

    """
    def __init__(self, name, default_value):
        """Initializes ModelOption with a name and default value"""
        self.name = name
        self.default_value = default_value


opts = [ModelOption('number_spheres', ''),
        ModelOption('sphere_position_file', 'at_bottom'),
        ModelOption('length_scale_factor', 1),
        ModelOption('real_ref_index_scale_factor', 1),
        ModelOption('imag_ref_index_scale_factor', 1),
        ModelOption('real_chiral_factor', 0),
        ModelOption('imag_chiral_factor', 0),
        ModelOption('medium_real_ref_index', 1),
        ModelOption('medium_imag_ref_index', 0),
        ModelOption('medium_real_chiral_factor', 0),
        ModelOption('medium_imag_chiral_factor', 0),
        ModelOption('target_euler_angles_deg', [0, 0, 0]),
        ModelOption('mie_epsilon', 1e-6),
        ModelOption('translation_epsilon', 1e-8),
        ModelOption('solution_epsilon', 1e-8),
        ModelOption('iterations_per_correction', 20),
        ModelOption('max_number_iterations', 2000),
        ModelOption('near_field_translation_distance', 1e6),
        ModelOption('store_translation_matrix', 0),
        ModelOption('fixed_or_random_orientation', 0),
        ModelOption('gaussian_beam_constant', 0),
        ModelOption('gaussian_beam_focal_point', [0, 0, 0]),
        ModelOption('output_file', 'test.dat'),
        ModelOption('run_print_file', 'run_print.dat'),
        ModelOption('write_sphere_data', 1),
        ModelOption('incident_or_target_frame', 0),
        ModelOption('min_scattering_angle_deg', 0),
        ModelOption('max_scattering_angle_deg', 180),
        ModelOption('min_scattering_plane_angle_deg', 0),
        ModelOption('max_scattering_plane_angle_deg', 0),
        ModelOption('delta_scattering_angle_deg', ''),
        ModelOption('number_scattering_angles', ''),
        ModelOption('scattering_angle_file', ''),
        ModelOption('incident_azimuth_angle_deg', 0),
        ModelOption('incident_polar_angle_deg', 0),
        ModelOption('calculate_scattering_coefficients', 1),
        ModelOption('scattering_coefficient_file', 'amn-temp.dat'),
        ModelOption('track_iterations', 1),
        ModelOption('azimuth_average_scattering_matrix', 0),
        ModelOption('calculate_near_field', 0),
        ModelOption('near_field_plane_coord', 1),
        ModelOption('near_field_plane_position', 0),
        ModelOption('near_field_plane_vertices', ''),
        ModelOption('spacial_step_size', 0.1),
        ModelOption('polarization_angle_deg', 0),
        ModelOption('near_field_output_file', 'nf-temp.dat'),
        ModelOption('near_field_output_data', 1),
        ModelOption('plane_wave_epsilon', 1e-3),
        ModelOption('calculate_t_matrix', 1),
        ModelOption('t_matrix_file', 'tmatrix-temp.dat'),
        ModelOption('t_matrix_convergence_epsilon', 1e-6),
        ModelOption('sm_number_processors', 10)]
opt_dict = {o.name: o for o in opts}


class ModelOptionValue(object):
    """A class to associate a ModelOption object with a value set at runtime

    This class holds a ModelOption object and an associated value for a single
    option to be passed to the MSTM model in an input file.

    Attributes:
        option (ModelOption): the ModelOption object to which we want to assign
            a value
        value: the value of that model option

    """
    def __init__(self, option, value=None):
        """Initializes ModelOptionValue with an option and optional value"""
        self.option = option
        if value is None:
            self.value = option.default_value
        else:
            self.value = value

    @classmethod
    def mov_from_name(cls, name, value=None):
        if name in opt_dict:
            return cls(opt_dict[name], value)
        else:
            raise AttributeError('{} is not a valid '
                                 'ModelOption name'.format(name))

    def formatted_option(self):
        return self.option.name + '\n' + str(self.value)


class RunType(Enum):
    FIXED = 0
    RANDOM = 1


class ModelRun(object):
    def __init__(self, name=None, fixed_or_random=RunType.FIXED,
                 option_val_list=None):
        self.name = name
        self.fixed_or_random = fixed_or_random
        if option_val_list is None:
            self.option_val_list = []
            if fixed_or_random == RunType.FIXED:
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('number_spheres'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('sphere_position_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('length_scale_factor'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'real_ref_index_scale_factor'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'imag_ref_index_scale_factor'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('medium_real_ref_index'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('medium_imag_ref_index'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'fixed_or_random_orientation', 0))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('output_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('run_print_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'scattering_coefficient_file'))
            elif fixed_or_random == RunType.RANDOM:
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('number_spheres'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('sphere_position_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('length_scale_factor'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'real_ref_index_scale_factor'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'imag_ref_index_scale_factor'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('medium_real_ref_index'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('medium_imag_ref_index'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'fixed_or_random_orientation', 1))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('output_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('run_print_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name(
                        'scattering_coefficient_file'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('calculate_t_matrix'))
                self.option_val_list.append(
                    ModelOptionValue.mov_from_name('t_matrix_file'))
            else:
                sys.exit('fixed_or_random is an invalid value')
        else:
            self.option_val_list = option_val_list

    def set_option(self, option_name, value=None):
        current_opts = {opt_val.option.name:
                        [opt_val.value, self.option_val_list.index(opt_val)]
                        for opt_val in self.option_val_list}
        if option_name in current_opts:
            self.option_val_list[current_opts[option_name][1]].value = value
        else:
            self.option_val_list.append(
                ModelOptionValue.mov_from_name(option_name, value))

    def formatted_options(self):
        options = ['']
        for option_val in self.option_val_list:
            if option_val is not self.option_val_list[len(
                    self.option_val_list) - 1]:
                options.append(option_val.formatted_option())
                options.append('\n')
            else:
                options.append(option_val.formatted_option())

        return ''.join(options)


class Grain(object):
    def __init__(self, host_sphere=None, rim=None, inclusions=None,
                 host_radius=None, rim_thickness=None, num_inclusions=None):
        self.host_sphere = host_sphere
        self.rim = rim
        self.inclusions = inclusions
        self.host_radius = host_radius
        self.rim_thickness = rim_thickness
        self.num_inclusions = num_inclusions

    @classmethod
    def new_grain(cls, x, y, z, r_host, dr_rim, r_incl, n_incl):
        chost_sphere = Sphere(x, y, z, r_host)
        crim = Sphere(x, y, z, r_host + dr_rim)
        cnum_inclusions = 0
        cinclusions = []
        while cnum_inclusions < n_incl:
            theta = 2 * math.pi * random.random()
            phi = math.acos(2 * random.random()-1)

            arandom = random.random()
            r = (((dr_rim - (2 * r_incl) - 0.0004) * arandom) +
                 ((r_host + r_incl) + 0.0002))
            this_x = (r * math.cos(theta) * math.sin(phi)) + x
            this_y = (r * math.sin(theta) * math.sin(phi)) + y
            this_z = (r * math.cos(phi)) + z

            j = 0
            while j < cnum_inclusions:
                p1 = (this_x, this_y, this_z)
                p2 = (cinclusions[j].x, cinclusions[j].y,
                      cinclusions[j].z)
                distance = math.sqrt(sum([(a - b) ** 2 for a, b in
                                          zip(p1, p2)]))
                if distance <= (2.01 * r_incl):
                    break
                j += 1
            else:
                cinclusions.append(Sphere(this_x, this_y, this_z, r_incl))
                cnum_inclusions += 1

        return cls(host_sphere=chost_sphere, rim=crim,
                   inclusions=cinclusions, host_radius=r_host,
                   rim_thickness=dr_rim, num_inclusions=cnum_inclusions)

    # move grain
    def move_to(self, x, y, z):
        old_x = self.host_sphere.x
        old_y = self.host_sphere.y
        old_z = self.host_sphere.z
        dx = x - old_x
        dy = y - old_y
        dz = z - old_z
        self.host_sphere.x += dx
        self.host_sphere.y += dy
        self.host_sphere.z += dz
        self.rim.x += dx
        self.rim.y += dy
        self.rim.z += dz
        for inclusion in self.inclusions:
            inclusion.x += dx
            inclusion.y += dy
            inclusion.z += dz

    def set_grain_oc(self, host_n, host_k, rim_n, rim_k, incl_n, incl_k):
        self.host_sphere.n = host_n
        self.host_sphere.k = host_k
        self.rim.n = rim_n
        self.rim.k = rim_k
        for inclusion in self.inclusions:
            inclusion.n = incl_n
            inclusion.k = incl_k

    def get_rxyznk(self):
        mylist = [self.host_sphere.get_tabbed_rxyznk(), '\n',
                  self.rim.get_tabbed_rxyznk(), '\n']

        for c, inclusion in enumerate(self.inclusions):
            mylist.append(inclusion.get_tabbed_rxyznk())
            if c < len(self.inclusions) - 1:
                mylist.append('\n')

        return ''.join(mylist)


@attr.s
class Cluster(object):
    grainlist = attr.ib(default=None)

    def get_bounding_sphere(self):
        # max distance between grain centers + 2*rim_r
        # centered on midpoint between grain centers if r1=r2
        max_distance = 0
        center = []
        for grain_a in self.grainlist:
            for grain_b in self.grainlist:
                if grain_a.rim is not None and grain_b.rim is not None:
                    distance = (grain_a.rim.distance_from(
                                another_sphere=grain_b.rim) + grain_a.rim.r +
                                grain_b.rim.r)
                    if distance > max_distance:
                        max_distance = distance
                        center.append((grain_a.rim.x + grain_b.rim.x) / 2)
                        center.append((grain_a.rim.y + grain_b.rim.y) / 2)
                        center.append((grain_a.rim.z + grain_b.rim.z) / 2)
                elif grain_a.rim is not None and grain_b.rim is None:
                    distance = (grain_a.rim.distance_from(
                                another_sphere=grain_b.host) + grain_a.rim.r +
                                grain_b.host.r)
                    if distance > max_distance:
                        max_distance = distance
                        center.append((grain_a.rim.x + grain_b.host.x) / 2)
                        center.append((grain_a.rim.y + grain_b.host.y) / 2)
                        center.append((grain_a.rim.z + grain_b.host.z) / 2)
                elif grain_a.rim is None and grain_b.rim is not None:
                    distance = (grain_a.host.distance_from(
                                another_sphere=grain_b.rim) + grain_a.host.r +
                                grain_b.rim.r)
                    if distance > max_distance:
                        max_distance = distance
                        center.append((grain_a.host.x + grain_b.rim.x) / 2)
                        center.append((grain_a.host.y + grain_b.rim.y) / 2)
                        center.append((grain_a.host.z + grain_b.rim.z) / 2)
                elif grain_a.rim is None and grain_b.rim is None:
                    distance = (grain_a.host.distance_from(
                                another_sphere=grain_b.host) + grain_a.host.r +
                                grain_b.host.r)
                    if distance > max_distance:
                        max_distance = distance
                        center.append((grain_a.host.x + grain_b.host.x) / 2)
                        center.append((grain_a.host.y + grain_b.host.y) / 2)
                        center.append((grain_a.host.z + grain_b.host.z) / 2)

        return Sphere(center[0], center[1], center[2], max_distance/2)

    def get_packing_fraction(self):
        vol_of_grains = 0
        bounding_sphere = self.get_bounding_sphere()
        vol_of_bounding_sphere = (4/3) * math.pi * bounding_sphere.r ** 3
        for grain in self.grainlist:
            if grain.rim is not None:
                vol_of_grains += (4/3) * math.pi * grain.rim.r ** 3
            elif grain.host is not None:
                vol_of_grains += (4/3) * math.pi * grain.host.r ** 3
            else:
                # Should not get here, but may need to handle this later
                pass

        return vol_of_grains/vol_of_bounding_sphere

    # geometric cross section
    def get_geom_xsection(self):
        bounding_sphere = self.get_bounding_sphere()
        return math.pi * bounding_sphere.r ** 2


class Pack(object):
    """A class containing the data generated by a PackLSD run.

    This class holds the metadata and x,y,z coordinates generated by a PackLSD
    run.

    Attributes:


    """
    def __init__(self, dims=None, num_of_particles=None, dispersity=None,
                 sphere_radius=None, sphere_coords=None):
        if sphere_coords is None:
            sphere_coords = [[]]
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
            kwargs = {'newline': ''}

        try:
            with open(a_filename, raccess, **kwargs) as infile:
                rawfile = infile.read().splitlines()
            if not rawfile:
                sys.exit('No data in file ' + a_filename)
        except IOError as e:
            sys.exit('I/O error: file {}: {}'.format(a_filename, e))

        lines = []
        for line in rawfile:
            lines.append(line.split())

        coords = [[float(item) for item in row] for row in lines[6:]]

        return(cls(int(lines[0][0]), int(lines[1][0]), int(lines[1][1]),
                   float(lines[3][0])/2, coords))

    def center_pack(self):
        xtotal = 0
        ytotal = 0
        ztotal = 0
        for line in self.sphere_coords:
            xtotal += line[0]
            ytotal += line[1]
            ztotal += line[2]
        length = len(self.sphere_coords)
        avgs = [xtotal/length, ytotal/length, ztotal/length]
        self.sphere_coords = [[element - avg for element, avg in
                               zip(line, avgs)] for line in self.sphere_coords]

    def rescale_pack(self, new_sphere_radius):
        multiplier = new_sphere_radius / self.sphere_radius
        self.sphere_coords = [[coord * multiplier for coord in row] for row in
                              self.sphere_coords]
        self.sphere_radius = new_sphere_radius


class ModelRunTypeError(Exception):
    def __init__(self, message):
        super(ModelRunTypeError, self).__init__(message)


class RunOutput(object):
    def __init__(self, input_dir_name, theta_i=30, theta_e=0):
        def h(a_mu, a_gamma):
            top = 1 + (2 * a_mu)
            bottom = 1 + (2 * a_gamma * a_mu)
            return top/bottom

        input_dir_name = input_dir_name.strip()
        if not input_dir_name.endswith('/'):
            input_dir_name.append('/')
        self.input_dir_name = input_dir_name
        file_list = glob.glob('{}dat/*.dat'.format(input_dir_name))
        self.wl_dict = dict()
        for file in file_list:
            wl = int(re.search(r'(\d{3,4})nm\.dat', file).group(1))
            self.wl_dict[wl] = SingleRunOutput.from_dat_file(file)
            if self.wl_dict[wl].unpol_qext is None:
                self.wl_dict[wl].incomplete = True
            else:
                self.wl_dict[wl].incomplete = False
        self.wl_list = [*self.wl_dict]
        self.wl_list.sort()

        stats_file = glob.glob('{}*stats.txt'.format(input_dir_name))
        if stats_file:
            try:
                with open(stats_file[0], 'r') as s:
                    for line in s:
                        if re.match(r'runname:', line) is not None:
                            self.runname = re.search(r'runname:\s(.+)\n',
                                                     line).group(1)
                            continue
                        if re.match('cluster packing', line) is not None:
                            self.pack_frac = float(re.search(r':\s(.+)\n',
                                                             line).group(1))
                            continue
                        if re.match('cluster x-section', line) is not None:
                            self.x_sect = float(re.search(r':\s(.+)\n', line
                                                          ).group(1))
                            continue
                        if re.match('cluster bounding', line) is not None:
                            self.bound_r = float(re.search(r':\s(.+)\n', line)
                                                 .group(1))
                            continue
            except IOError as e:
                sys.exit('I/O error: file {}: {}'.format(stats_file, e))
        else:
            raise ModelRunTypeError('This model run does not have a stats file'
                                    '.')

        for wl in self.wl_list:
            mu0 = math.cos(math.radians(theta_i))
            mu = math.cos(math.radians(theta_e))
            self.wl_dict[wl].k = (2 * math.pi) / wl
            self.wl_dict[wl].dcsca = \
                self.wl_dict[wl].s_matrix_dict['150.00']['11']
            self.wl_dict[wl].csca = self.wl_dict[wl].unpol_qsca * self.x_sect
            self.wl_dict[wl].p_g = 4 * math.pi * (1 / self.wl_dict[wl].csca) \
                * self.wl_dict[wl].dcsca
            self.wl_dict[wl].ssa = (self.wl_dict[wl].unpol_qsca /
                                    self.wl_dict[wl].unpol_qext)
            gamma = math.sqrt(1 - self.wl_dict[wl].ssa)
            hmu0 = h(mu0, gamma)
            hmu = h(mu, gamma)
            self.wl_dict[wl].hapke_refl = (self.wl_dict[wl].ssa /
                                           (4 * math.pi)) * (mu0 / (mu0 + mu))\
                * ((1 + 0) * self.wl_dict[wl].p_g + (hmu0 * hmu)-1)

    def write_json(self, filename):
        run_dict = self.__dict__
        for key, val in run_dict.items():
            if isinstance(val, dict):
                for inner_key, inner_val in val.items():
                    if isinstance(inner_val, SingleRunOutput):
                        val[inner_key] = inner_val.__dict__
        try:
            with open(filename, 'w') as o:
                json.dump(run_dict, o, indent=1)
        except IOError as e:
            sys.exit('I/O error: file {}: {}'.format(filename, e))

    def write_pickle(self, filename):
        try:
            with open(filename, 'wb') as p:
                pickle.dump(self, p)
        except IOError as e:
            sys.exit('I/O error: file {}: {}'.format(filename, e))


class SingleRunOutput(object):
    def __init__(self):
        self.run_number = None
        self.n_spheres = None
        self.vol_size_param = None
        self.position_file = None
        self.output_file = None
        self.length_scale_factor = None
        self.real_ref_index_scale_factor = None
        self.imag_ref_index_scale_factor = None
        self.real_chiral_factor = None
        self.imag_chiral_factor = None
        self.incident_or_target_frame = None
        self.theta_min = None
        self.theta_max = None
        self.azimuth_averaged = None
        self.num_scattering_angles = None
        self.eps_mie = None
        self.eps_soln = None
        self.max_iterations = None
        self.real_medium_ref_index = None
        self.imag_medium_ref_index = None
        self.target_euler_rotation_x = None
        self.target_euler_rotation_y = None
        self.target_euler_rotation_z = None
        self.far_field_kr = None
        self.iter_per_corr = None
        self.incident_azimuth_angle = None
        self.incident_polar_angle = None
        self.common_exp_eps = None
        self.sc_file = None
        self.sphere_data_dict = None
        self.unpol_qext = None
        self.unpol_qabs = None
        self.unpol_qsca = None
        self.asym = None
        self.parpol_qext = None
        self.parpol_qabs = None
        self.parpol_qsca = None
        self.perpol_qext = None
        self.perpol_qabs = None
        self.perpol_qsca = None
        self.s_matrix_dict = None

    @classmethod
    def from_dat_file(cls, filename):
        thisrun = SingleRunOutput()
        thisrun.azimuth_averaged = False
        try:
            with open(filename, 'r') as f:
                for line in f:
                    if (re.search('input parameters for run number', line)
                            is not None):
                        thisrun.run_number = int(next(f).strip())
                        continue
                    if re.search('number of spheres', line) is not None:
                        thisrun.n_spheres, thisrun.vol_size_param = \
                            next(f).split()
                        thisrun.n_spheres = int(thisrun.n_spheres)
                        thisrun.vol_size_param = float(thisrun.vol_size_param)
                        continue
                    if re.search('position file', line) is not None:
                        thisrun.position_file = next(f).strip()
                        continue
                    if re.search('output file', line) is not None:
                        thisrun.output_file = next(f).strip()
                        continue
                    if re.search('length, ref\. indx', line) is not None:
                        thisrun.length_scale_factor, \
                            thisrun.real_ref_index_scale_factor, \
                            thisrun.imag_ref_index_scale_factor = \
                            [float(i) for i in next(f).split()]
                        continue
                    if re.search('chiral factors:', line) is not None:
                        thisrun.real_chiral_factor, \
                            thisrun.imag_chiral_factor = \
                            [float(i) for i in next(f).split()]
                        continue
                    if re.search('based on incident frame', line) is not None:
                        thisrun.incident_or_target_frame = 'incident'
                        continue
                    elif re.search('based on target frame', line) is not None:
                        thisrun.incident_or_target_frame = 'target'
                        continue
                    if re.search('thetamin, thetamax', line) is not None:
                        thisrun.theta_min, thisrun.theta_max = \
                            [float(i) for i in next(f).split()]
                        continue
                    if (re.search('scattering matrix is averaged', line) is
                            not None):
                        thisrun.azimuth_averaged = True
                        continue
                    if re.search('number scattering angles', line) is not None:
                        thisrun.num_scattering_angles = int(next(f).strip())
                        continue
                    if re.search('epsmie, epssoln', line) is not None:
                        temp = next(f).strip()
                        captures = re.match(
                            r'(\d*\.\d+E[+|\-]\d{2})\s*(\d*\.\d+E[+|\-]'
                            r'\d{2})\s*(\d+)', temp)
                        thisrun.eps_mie = float(captures.group(1))
                        thisrun.eps_soln = float(captures.group(2))
                        thisrun.max_iterations = int(captures.group(3))
                        continue
                    if re.search('medium refractive index', line) is not None:
                        thisrun.real_medium_ref_index, \
                            thisrun.imag_medium_ref_index = [float(i) for i in
                                                             next(f).split()]
                        continue
                    if re.search('target euler', line) is not None:
                        thisrun.target_euler_rotation_x, \
                            thisrun.target_euler_rotation_y, \
                            thisrun.target_euler_rotation_z = [float(i) for i
                                                               in
                                                               next(f).split()]
                        continue
                    if re.search('far field kr', line) is not None:
                        thisrun.far_field_kr, thisrun.iter_per_corr = \
                            next(f).split()
                        thisrun.far_field_kr = float(thisrun.far_field_kr)
                        thisrun.iter_per_corr = int(thisrun.iter_per_corr)
                        continue
                    if re.search('incident azimuth', line) is not None:
                        thisrun.incident_azimuth_angle, \
                            thisrun.incident_polar_angle = [float(i) for i in
                                                            next(f).split()]
                        continue
                    if re.search('common expansion', line) is not None:
                        thisrun.common_exp_eps = float(next(f).strip())
                        continue
                    if (re.search('scattering coefficients calculated', line)
                            is not None):
                        thisrun.sc_file = next(f).strip()
                        continue
                    if (re.search('calculation results for run', line) is
                            not None):
                        thisrun.run_number = int(next(f).strip())
                        continue
                    if re.search(r'\s*sphere\s*host', line) is not None:
                        headers = line.split()
                        line = next(f)
                        thisrun.sphere_data_dict = dict()
                        while re.match(r'^\s*\d+', line):
                            thisrow = line.split()
                            thisdict = {a: float(b) for a, b in
                                        zip(headers[1:], thisrow[1:])}
                            thisdict['host'] = int(thisdict['host'])
                            thisrun.sphere_data_dict[int(thisrow[0])] = \
                                thisdict
                            line = next(f)
                    if re.search('unpolarized total ext', line) is not None:
                        thisrun.unpol_qext, thisrun.unpol_qabs, \
                            thisrun.unpol_qsca, thisrun.asym = \
                            [float(i) for i in next(f).split()]
                        continue
                    if re.search('parallel total ext', line) is not None:
                        thisrun.parpol_qext, thisrun.parpol_qabs, \
                            thisrun.parpol_qsca = [float(i) for i in
                                                   next(f).split()]
                        continue
                    if re.search('perpendicular total ext', line) is not None:
                        thisrun.perpol_qext, thisrun.perpol_qabs, \
                            thisrun.perpol_qsca = [float(i) for i in
                                                   next(f).split()]
                        continue
                    if (re.search('scattering matrix elements', line) is not
                            None):
                        line = next(f)
                        headers = line.split()
                        thisrun.s_matrix_dict = dict()
                        line = next(f)
                        while re.match(r'^\s*\d+\.\d{2}', line):
                            thisrow = line.split()
                            thisdict = {a: float(b) for a, b in
                                        zip(headers[1:], thisrow[1:])}
                            thisrun.s_matrix_dict[thisrow[0]] = thisdict
                            try:
                                line = next(f)
                            except StopIteration:
                                break
                        continue
        except IOError as e:
            sys.exit('I/O error: file {}: {}'.format(filename, e))

        return thisrun
