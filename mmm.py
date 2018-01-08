# -*- coding: utf-8 -*-
""" Classes for MSTM Model Manager's data model

This module provides the classes necessaryfor the MVC "Model" portion of the
MSTM Model Manager program.

Classes:
    Sphere
    ModelOptions

"""

from math import sqrt, pi
import attr
import numpy
import random

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

    def distance_from(self, x1, y1, z1):
        """Return the distance of the center of this sphere from the provided
        coordinates.

        Args:
            x1 -- the x coord of the point from which distance is calculated
            y1 -- the y coord of the point from which distance is calculated
            z1 -- the z coord of the point from which distance is calculated

        Returns:
            A float representing the distance from the provided coords to the
            center of this sphere.

        """
        return sqrt(
            (self.x - x1) ** 2 + (self.y - y1) ** 2 + (self.z - z1) ** 2)

    def geom_x_sec(self):
        """Return the geometric cross section of this sphere."""
        return pi * self.r ** 2

@attr.s
class ModelOptions:
    number_spheres = attr.ib(default="")
    sphere_position_file = attr.ib(default="at_bottom")
    length_scale_factor = attr.ib(default=1)
    real_ref_index_scale_factor = attr.ib(default="1")
    imag_ref_index_scale_factor = attr.ib(default="1")
    real_chiral_factor = attr.ib(default=0)
    imag_chiral_factor = attr.ib(default=0)
    medium_real_ref_index = attr.ib(default=1)
    medium_imag_ref_index = attr.ib(default=0)
    medium_real_chiral_factor = attr.ib(default=0)
    medium_imag_chiral_factor = attr.ib(default=0)
    target_euler_angles_deg = attr.ib(default="0 0 0")
    mie_epsilon = attr.ib(default=1e-6)
    translation_epsilon = attr.ib(default=1e-8)
    solution_epsilon = attr.ib(default=1e-8)
    iterations_per_correction = attr.ib(default=20)
    max_number_iterations = attr.ib(default=2000)
    near_field_translation_distance = attr.ib(default=1e6)
    store_translation_matrix = attr.ib(default=0)
    fixed_or_random_orientation = attr.ib(default=0)
    gaussian_beam_constant = attr.ib(default=0)
    gaussian_beam_focal_point = attr.ib(default="0 0 0")
    output_file = attr.ib(default="test.dat")
    run_print_file = attr.ib(default="run_print.dat")
    write_sphere_data = attr.ib(default=1)
    incident_or_target_frame = attr.ib(default=0)
    min_scattering_angel_deg = attr.ib(default=0)
    max_scattering_angle_deg = attr.ib(default=180)
    min_scattering_plane_angle_deg = attr.ib(default=0)
    max_scattering_plane_angle_deg = attr.ib(default=0)
    delta_scattering_angle_deg = attr.ib(default="")
    number_scattering_angles = attr.ib(default="")
    scattering_angle_file = attr.ib(default="")
    incident_azimuth_angle_deg = attr.ib(default=0)
    incident_polar_angle_deg = attr.ib(default=0)
    calculate_scattering_coefficients = attr.ib(default=1)
    scattering_coefficient_file = attr.ib(default="amn-temp.dat")
    track_iterations = attr.ib(default=1)
    azimuth_average_scattering_matrix = attr.ib(default=0)
    calculate_near_field = attr.ib(default=0)
    near_field_plane_coord = attr.ib(default=1)
    near_field_plane_position = attr.ib(default=0)
    near_field_plane_vertices = attr.ib(default="")
    spacial_step_size = attr.ib(default=0.1)
    polarization_angle_deg = attr.ib(default=0)
    near_field_output_file = attr.ib(default="nf-temp.dat")
    near_field_output_data = attr.ib(default=1)
    plane_wave_epsilon = attr.ib(default=1e-3)
    calculate_t_matrix = attr.ib(default=1)
    t_matrix_file = attr.ib(default="tmatrix-temp.dat")
    t_matrix_convergence_epsilon = attr.ib(default=1e-6)
    sm_number_processors = attr.ib(default=10)

    def is_default(self,name):
        return (getattr(self, name) ==
            getattr(attr.fields(type(self)),name).default)

    def get_formatted(self, setname):
        output = ""
        required_keys = [ "number_spheres", "sphere_position_file",
            "length_scale_factor", "real_ref_index_scale_factor",
            "imag_ref_index_scale_factor", "medium_real_ref_index",
            "medium_imag_ref_index", "fixed_or_random_orientation",
            "output_file", "run_print_file", "scattering_coefficient_file" ]

        if setname == "defaults":
            for key in self.__dict__:
                if self.is_default(key):
                    if len(output) != 0:
                        output = output + "\n"
                    output = output + key + "\n" + str(self.__dict__[key])
            return output
        elif setname == "required":
            #
            # if fixed_or_random_orientation == 1
            # t_matrix_file
            for key in required_keys:
                if len(output) != 0:
                    output = output + "\n"
                output = output + key + "\n" + str(self.__dict__[key])
            if self.__dict__['fixed_or_random_orientation'] == 1:
                output = ("\n" + 't_matrix_file' + "\n" +
                    self.__dict__['t_matrix_file'])
            return output
        elif setname == "non-defaults":
            for key in self.__dict__:
                if not self.is_default(key):
                    if len(output) != 0:
                        output = output + "\n"
                    output = output + key + "\n" + str(self.__dict__[key])
            return output
        elif setname == "required-and-non-defaults":
            for key in self.__dict__:
                if self.is_default(key) or key in required_keys:
                    if len(output) != 0:
                        output = output + "\n"
                    output = output + key + "\n" + str(self.__dict__[key])
            return output
        elif setname == "all":
            for key in self.__dict__:
                if len(output) != 0:
                    output = output + "\n"
                output = output + key + "\n" + str(self.__dict__[key])
            return output
        else:
            return False

@attr.s
class Grain:
    host_sphere = attr.ib(default=None)
    rim = attr.ib(default=None)
    inclusions = attr.ib(default=None)
    host_radius = attr.ib(default=None)
    rim_thickness = attr.ib(default=None)
    num_inclusions = attr.ib(default=None)

    @classmethod
    def new_grain(cls,x,y,z,r_host,dr_rim,r_incl,n_incl):
        chost_sphere = mmm.Sphere(x,y,z,r_host)
        crim = mmm.Sphere(x,y,z,r_host+dr_rim)
        cnum_inclusions = 0
        cinclusions = []
        while cnum_inclusions < n_incl:
            theta = 2*math.pi*random.random()
            phi = math.acos(2*random.random()-1)

            random = random.random()
            r = ((dr_rim - r_incl - 0.0002) * random) + ((r_host + r_incl) +
                    0.0001)
            this_x = r * math.cos(theta) * math.sin(phi)
            this_y = r * math.sin(theta) * math.sin(phi)
            this_z = r * math.cos(phi)

            j = 0
            while j < cnum_inclusions:
                p1 = (this_x, this_y, this_z)
                p2 = (cinclusions[j].x, cinclusions[j].y,
                        cinclusions[j].z)
                distance = math.sqrt(sum([(a - b) ** 2 for a, b in zip(p1,
                    p2)]))
                if distance <= (2 * r_incl):
                    break
                j += 1
            else:
                cinclusions.append(mmm.Sphere(r_incl,this_x,this_y,this_z))
                i += 1
        return cls(host_sphere=chost_sphere, rim=crim,
                inclusions=cinclusions, host_radius = r_host, rim_thickness =
                dr_rim, num_inclusions = cnum_inclusions)


    # move grain
    def move_to(self, x, y, z):
        old_x = host_sphere.x
        old_y = host_sphere.y
        old_z = host_sphere.z
        dx = x - old_x
        dy = y - old_y
        dz = z - old_z
        self.host_sphere.x += dx
        self.host_sphere.y += dy
        self.host_sphere.z += dz
        self.rim.x += dx
        self.rim.y += dy
        self.rim.z += dz
        for inclusion in inclusions:
            inclusion.x += dx
            inclusion.y += dy
            inclusion.z += dz
    # weight percent inclusions?


@attr.s
class Cluster:
    grainlist = attr.ib(default=None)

    # bounding sphere
        # max distance between grain centers + 2*rim_r
        # centered on midpoint between grain centers if r1=r2
    # packing fraction
        # volume of spheres/volume of bounding sphere
    # geometric cross section
        # pi*r_bounding^2
