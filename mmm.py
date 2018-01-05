# -*- coding: utf-8 -*-
""" Classes for MSTM Model Manager's data model

This module provides the classes necessaryfor the MVC "Model" portion of the
MSTM Model Manager program.

Classes:
    Sphere
    SphereGroup

"""

from math import sqrt, pi

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

class SphereGroup:
    def __init__(self, sphere_list=list(), group_name=None, group_n=None,
            group_k=None, group_real_chiral=None,group_imag_chiral=None,
            tmatrix_file=None):
        self.sphere_list = sphere_list
        self.short_name = group_name
        self.group_n = group_n
        self.group_k = group_k
        self.group_real_chiral = group_real_chiral
        self.group_imag_chiral = group_imag_chiral
        self.tmatrix_file = tmatrix_file

        if ((self.group_n is not None)
            or (self.group_k is not None)
            or (self.group_real_chiral is not None)
            or (self.group_imag_chiral is not None)):
            self.apply_oc(group_n, group_k, group_real_chiral,
                    group_imag_chiral)

    def add_sphere(self, a_sphere):
        self.sphere_list.append(a_sphere)

    def apply_oc(self, n=None, k=None, real_chiral=None, imag_chiral=None):
        if n is not None:
            for i in self.sphere_list:
                i.n = n
        if k is not None:
            for i in self.sphere_list:
                i.k = k
        if real_chiral is not None:
            for i in self.sphere_list:
                i.real_chiral = real_chiral
        if imag_chiral is not None:
            for i in self.sphere_list:
                i.imag_chiral = imag_chiral
