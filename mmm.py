# -*- coding: utf-8 -*-
""" Classes for MSTM Model Manager's data model

This module provides the classes necessary for the MVC "Model" portion of the
MSTM Model Manager program.

Classes:
    Sphere
    ModelOptions

"""

import math
import attr
import random
import sys
from enum import Enum


class Sphere:
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
        """Return the geometric cross section of this sphere."""
        return math.pi * self.r ** 2


class ModelOption:
    def __init__(self, name, default_value, group):
        self.name = name
        self.default_value = default_value
        self.group = group


opts = [ModelOption('number_spheres', 'target', ''),
        ModelOption('sphere_position_file', 'target', 'at_bottom'),
        ModelOption('length_scale_factor', 1, 'target'),
        ModelOption('real_ref_index_scale_factor', 1, 'target'),
        ModelOption('imag_ref_index_scale_factor', 1, 'target'),
        ModelOption('real_chiral_factor', 0, 'target'),
        ModelOption('imag_chiral_factor', 0, 'target'),
        ModelOption('medium_real_ref_index', 1, 'target'),
        ModelOption('medium_imag_ref_index', 0, 'target'),
        ModelOption('medium_real_chiral_factor', 0, 'target'),
        ModelOption('medium_imag_chiral_factor', 0, 'target'),
        ModelOption('target_euler_angles_deg', [0, 0, 0], 'target'),
        ModelOption('mie_epsilon', 1e-6, 'solution'),
        ModelOption('translation_epsilon', 1e-8, 'solution'),
        ModelOption('solution_epsilon', 1e-8, 'solution'),
        ModelOption('iterations_per_correction', 20, 'solution'),
        ModelOption('max_number_iterations', 2000, 'solution'),
        ModelOption('near_field_translation_distance', 1e6, 'solution'),
        ModelOption('store_translation_matrix', 0, 'solution'),
        ModelOption('fixed_or_random_orientation', 0, 'misc'),
        ModelOption('gaussian_beam_constant', 0, 'misc'),
        ModelOption('gaussian_beam_focal_point', [0, 0, 0], 'misc'),
        ModelOption('output_file', 'test.dat', 'misc'),
        ModelOption('run_print_file', 'run_print.dat', 'misc'),
        ModelOption('write_sphere_data', 1, 'misc'),
        ModelOption('incident_or_target_frame', 0, 'smatrix'),
        ModelOption('min_scattering_angle_deg', 0, 'smatrix'),
        ModelOption('max_scattering_angle_deg', 180, 'smatrix'),
        ModelOption('min_scattering_plane_angle_deg', 0, 'smatrix'),
        ModelOption('max_scattering_plane_angle_deg', 0, 'smatrix'),
        ModelOption('delta_scattering_angle_deg', '', 'smatrix'),
        ModelOption('number_scattering_angles', '', 'smatrix'),
        ModelOption('scattering_angle_file', '', 'smatrix'),
        ModelOption('incident_azimuth_angle_deg', 0, 'fixed'),
        ModelOption('incident_polar_angle_deg', 0, 'fixed'),
        ModelOption('calculate_scattering_coefficients', 1, 'fixed'),
        ModelOption('scattering_coefficient_file', 'amn-temp.dat', 'fixed'),
        ModelOption('track_iterations', 1, 'fixed'),
        ModelOption('azimuth_average_scattering_matrix', 0, 'fixed'),
        ModelOption('calculate_near_field', 0, 'fixed'),
        ModelOption('near_field_plane_coord', 1, 'fixed'),
        ModelOption('near_field_plane_position', 0, 'fixed'),
        ModelOption('near_field_plane_vertices', '', 'fixed'),
        ModelOption('spacial_step_size', 0.1, 'fixed'),
        ModelOption('polarization_angle_deg', 0, 'fixed'),
        ModelOption('near_field_output_file', 'nf-temp.dat', 'fixed'),
        ModelOption('near_field_output_data', 1, 'fixed'),
        ModelOption('plane_wave_epsilon', 1e-3, 'fixed'),
        ModelOption('calculate_t_matrix', 1, 'random'),
        ModelOption('t_matrix_file', 'tmatrix-temp.dat', 'random'),
        ModelOption('t_matrix_convergence_epsilon', 1e-6, 'random'),
        ModelOption('sm_number_processors', 10, 'random')]
opt_dict = {o.name: o for o in opts}


class ModelOptionValue:
    def __init__(self, option=None, value=None, print_me=True):
        self.option = option
        if value is None:
            self.value = option.default_value
        else:
            self.value = value
        self.print_me = print_me

    @classmethod
    def mov_from_name(cls, name, value=None, print_me=True):
        if name in opt_dict:
            return cls(opt_dict[name], value, print_me)
        else:
            raise AttributeError('{} is not a valid '
                                 'ModelOption name'.format(name))


class RunType(Enum):
    FIXED = 0
    RANDOM = 1


class ModelRun:
    def __init__(self, name=None, fixed_or_random=RunType.FIXED,
                 option_list=None):
        self.name = name
        self.fixed_or_random = fixed_or_random
        if option_list is None:
            self.option_list = {}
            if fixed_or_random == RunType.FIXED:
                option_list.append(
                    ModelOptionValue.mov_from_name('number_spheres'))
                option_list.append(
                    ModelOptionValue.mov_from_name('sphere_position_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name('length_scale_factor'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'real_ref_index_scale_factor'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'imag_ref_index_scale_factor'))
                option_list.append(
                    ModelOptionValue.mov_from_name('medium_real_ref_index'))
                option_list.append(
                    ModelOptionValue.mov_from_name('medium_imag_ref_index'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'fixed_or_random_orientation', 0))
                option_list.append(
                    ModelOptionValue.mov_from_name('output_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name('run_print_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'scattering_coefficient_file'))
            elif fixed_or_random == RunType.RANDOM:
                option_list.append(
                    ModelOptionValue.mov_from_name('number_spheres'))
                option_list.append(
                    ModelOptionValue.mov_from_name('sphere_position_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name('length_scale_factor'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'real_ref_index_scale_factor'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'imag_ref_index_scale_factor'))
                option_list.append(
                    ModelOptionValue.mov_from_name('medium_real_ref_index'))
                option_list.append(
                    ModelOptionValue.mov_from_name('medium_imag_ref_index'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'fixed_or_random_orientation', 1))
                option_list.append(
                    ModelOptionValue.mov_from_name('output_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name('run_print_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name(
                        'scattering_coefficient_file'))
                option_list.append(
                    ModelOptionValue.mov_from_name('calculate_t_matrix'))
                option_list.append(
                    ModelOptionValue.mov_from_name('t_matrix_file'))
            else:
                sys.exit('fixed_or_random is an invalid value')
        else:
            self.option_list = option_list


class Grain:
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
            r = ((dr_rim - (2 * r_incl) - 0.0004) * arandom) + \
                ((r_host + r_incl) + 0.0002)
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

    def print_rxyznk(self, tofile=sys.stdout):
        print(str(self.host_sphere.r) + '\t' + str(self.host_sphere.x) + '\t'
              + str(self.host_sphere.y) + '\t' + str(self.host_sphere.z) +
              '\t' + str(self.host_sphere.n) + '\t' +
              str(self.host_sphere.k), file=tofile)

        print(str(self.rim.r) + '\t' + str(self.rim.x) + '\t' +
              str(self.rim.y) + '\t' + str(self.rim.z) + '\t' +
              str(self.rim.n) + '\t' + str(self.rim.k), file=tofile)

        for inclusion in self.inclusions:
            print(str(inclusion.r) + '\t' + str(inclusion.x) + '\t' +
                  str(inclusion.y) + '\t' + str(inclusion.z) + '\t' +
                  str(inclusion.n) + '\t' + str(inclusion.k), file=tofile)


@attr.s
class Cluster:
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
                        center[0] = (grain_a.rim.x + grain_b.rim.x) / 2
                        center[1] = (grain_a.rim.y + grain_b.rim.y) / 2
                        center[2] = (grain_a.rim.z + grain_b.rim.z) / 2
                elif grain_a.rim is not None and grain_b.rim is None:
                    distance = (grain_a.rim.distance_from(
                                another_sphere=grain_b.host) + grain_a.rim.r +
                                grain_b.host.r)
                    if distance > max_distance:
                        max_distance = distance
                        center[0] = (grain_a.rim.x + grain_b.host.x) / 2
                        center[1] = (grain_a.rim.y + grain_b.host.y) / 2
                        center[2] = (grain_a.rim.z + grain_b.host.z) / 2
                elif grain_a.rim is None and grain_b.rim is not None:
                    distance = (grain_a.host.distance_from(
                                another_sphere=grain_b.rim) + grain_a.host.r +
                                grain_b.rim.r)
                    if distance > max_distance:
                        max_distance = distance
                        center[0] = (grain_a.host.x + grain_b.rim.x) / 2
                        center[1] = (grain_a.host.y + grain_b.rim.y) / 2
                        center[2] = (grain_a.host.z + grain_b.rim.z) / 2
                elif grain_a.rim is None and grain_b.rim is None:
                    distance = (grain_a.host.distance_from(
                                another_sphere=grain_b.host) + grain_a.host.r +
                                grain_b.host.r)
                    if distance > max_distance:
                        max_distance = distance
                        center[0] = (grain_a.host.x + grain_b.host.x) / 2
                        center[1] = (grain_a.host.y + grain_b.host.y) / 2
                        center[2] = (grain_a.host.z + grain_b.host.z) / 2

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


class Pack:
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
        print('nsr:' + str(new_sphere_radius))
        print('ssr:' + str(self.sphere_radius))
        print('multiplier:' + str(multiplier))
        print('sphere_coords:' + str(self.sphere_coords))
        self.sphere_coords = [[coord * multiplier for coord in row] for row in
                              self.sphere_coords]
        print('new_shere_coords:' + str(self.sphere_coords))
        self.sphere_radius = new_sphere_radius
