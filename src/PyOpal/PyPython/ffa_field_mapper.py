# Routines to generate field maps in Cartesian and cylindrical coordinate
# systems
#
# Copyright (c) 2023, Chris Rogers, STFC Rutherford Appleton Laboratory, Didcot, UK
#
# This file is part of OPAL.
#
# OPAL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# You should have received a copy of the GNU General Public License
# along with OPAL.  If not, see <https://www.gnu.org/licenses/>.
#

"""
Module to hold FFAFieldMapper
"""


import os
import math

import numpy
import matplotlib
import matplotlib.pyplot

import pyopal.objects.field


class FFAFieldMapper():
    """
    Class to make field maps, intended for FFAs/ring geometries
    """
    def __init__(self):
        # for cylindrical field map
        self.r_points = []
        self.phi_points = []

        # for cartesian field map
        self.x_points = []
        self.y_points = []

        self.verbose = 0
        self.cmap = "PiYG"

        self.radial_contours = []
        self.spiral_contours = []

        # for derivative calculations
        self.delta_x = 0.001
        self.delta_t = 0.001
        self.bmax = 1.0

        self.plot_dir = os.getcwd()

        # if empty, accept all tracks; else only those in track_id_list
        self.valid_track_id_list = []
        self.track_orbit_dict = {}

    @classmethod
    def binner(cls, a_list):
        """Build a set of bins based on a list of grid coordinates"""
        return [a_list[0]+(a_list[1]-a_list[0])*(i-0.5) \
                                                  for i in range(len(a_list)+1)]

    def load_tracks(self, track_orbit):
        """
        Load tracks from track_orbit file
        - track_orbit: string file name
        """
        self.track_orbit_dict = {}
        with open(track_orbit, encoding="utf8") as fin:
            for i, line in enumerate(fin.readlines()):
                words = line.split()
                if i < 2: # header lines
                    continue
                track_id = words[0]
                if track_id not in self.track_orbit_dict:
                    if self.valid_track_id_list and \
                       track_id not in self.valid_track_id_list:
                        continue
                    self.track_orbit_dict[track_id] = {
                        "x":[], "px":[], "y":[], "py":[], "r":[], "phi":[]
                    }
                pos_r = (float(words[1])**2+float(words[3])**2)**0.5
                phi = math.atan2(float(words[3]), float(words[1]))
                self.track_orbit_dict[track_id]["x"].append(float(words[1]))
                self.track_orbit_dict[track_id]["px"].append(float(words[2]))
                self.track_orbit_dict[track_id]["y"].append(float(words[3]))
                self.track_orbit_dict[track_id]["py"].append(float(words[4]))
                self.track_orbit_dict[track_id]["r"].append(pos_r)
                self.track_orbit_dict[track_id]["phi"].append(math.degrees(phi))

    def gen_cmap(self, bz_grid):
        """
        Generate the colour mapping for field bz
        - bz_grid: set of bz values. Colour mapping will proceed symmetrically
                   0 in the positive and negative direction. Maximum value is
                   the maximum absolute value in the bz_grid, up to self.bmax.
                   If the grid contains a value greater than self.bmax, then
                   self.bmax becomes the maximum value. To disable the bmax
                   behaviour, set self.bmax to None
        Returns a tuple of (min_bz in the grid,
                            max_bz in the grid,
                            max value of bz to be used in the colormap)
        """
        min_bz = min(bz_grid)
        max_bz = max(bz_grid)
        if self.bmax is None:
            cmax = max(abs(min_bz), abs(max_bz))
        else:
            cmax = self.bmax
        return min_bz, max_bz, cmax

    def field_map_cylindrical(self, axes = None):
        """
        Plot a field map in cylindrical coordinates.

        - axes: matplotlib Axes object or None. If defined, plot the field map
        on axes as a hist2d; if None make a new figure/axes and plot it there.

        Returns the figure, either a new figure or the parent figure to which
        axes belongs.
        """
        r_grid = []
        phi_grid = []
        bz_grid = []

        for radius in self.r_points:
            for phi in self.phi_points:
                r_grid.append(radius)
                phi_grid.append(phi)
                point = (radius*math.cos(math.radians(phi)),
                         radius*math.sin(math.radians(phi)),
                         0,
                         0)
                value = pyopal.objects.field.get_field_value(*point)
                bz_grid.append(value[3])
                if self.verbose > 0:
                    print("Field value at r, phi", radius, round(phi, 2),
                          "point", point,
                          "is B:", value[1:4],
                          "E:", value[4:])
        r_bins = self.binner(self.r_points)
        phi_bins = self.binner(self.phi_points)
        if not axes:
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)

        min_by, max_by, cmax = self.gen_cmap(bz_grid)
        axes.hist2d(phi_grid, r_grid, bins=[phi_bins, r_bins], weights=bz_grid,
                    cmin=min_by, cmax=max_by, cmap=self.cmap, vmin=-cmax, vmax=cmax)
        axes.set_xlabel("$\\phi$ [deg]")
        axes.set_ylabel("r [m]")
        axes.set_title("$B_z$ [T]")
        for contour in self.radial_contours:
            self.draw_cylindrical_radial_contour(axes, contour)
        for contour in self.spiral_contours:
            self.draw_cylindrical_spiral_contour(axes, contour)
        fig_fname = os.path.join(self.plot_dir, "scaling_ffa_map_cyl.png")
        figure.savefig(fig_fname)
        print("Generated cylindrical field map in", fig_fname)
        return figure

    @classmethod
    def draw_cylindrical_radial_contour(cls, axes, contour):
        """
        Draw a purely radial contour on axes
        - axes: matplotlib Axes object to draw on.
        - contour: dictionary (see default_contour for definitions)
        """
        print("Plotting cylindrical radial contour", contour)
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        axes.plot(xlim,
                  [contour["radius"]]*2,
                  linestyle=contour["linestyle"],
                  color=contour["colour"])
        axes.text(xlim[-1],
                  contour["radius"],
                  contour["label"],
                  horizontalalignment="right",
                  va="top",
                  color=contour["colour"])
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

    def draw_cylindrical_spiral_contour(self, axes, contour):
        """
        Draw a radially spiralling contour on axes
        - axes: matplotlib Axes object to draw on.
        - contour: dictionary (see default_contour for definitions)
        """
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        tan_d = math.tan(math.radians(contour["spiral_angle"]))
        phi_points = [math.radians(contour["phi0"]) + \
                         tan_d*math.log(r/contour["r0"]) for r in self.r_points]
        phi_points = [math.degrees(phi) for phi in phi_points]
        axes.plot(phi_points, self.r_points,
                  linestyle=contour["linestyle"],
                  color=contour["colour"])
        axes.text(phi_points[-1],
                  self.r_points[-1],
                  contour["label"],
                  va="top",
                  rotation="vertical",
                  color=contour["colour"])
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

    def field_map_cartesian(self, axes = None):
        """
        Plot a field map in cartesian coordinates.

        - axes: matplotlib Axes object or None. If defined, plot the field map
        on axes as a hist2d; if None make a new figure/axes and plot it there.

        Returns the figure, either a new figure or the parent figure to which
        axes belongs.
        """
        x_grid = []
        y_grid = []
        bz_grid = []
        for pos_x in self.x_points:
            for pos_y in self.y_points:
                x_grid.append(pos_x)
                y_grid.append(pos_y)
                point = (pos_x, pos_y, 0, 0)
                value = pyopal.objects.field.get_field_value(*point)
                bz_grid.append(value[3])
                if self.verbose > 0:
                    print("Field value at point", point,
                          "is B:", value[1:4], "E:", value[4:])
        x_bins = self.binner(self.x_points)
        y_bins = self.binner(self.y_points)
        if not axes:
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)
        min_by, max_by, cmax = self.gen_cmap(bz_grid)
        hist = axes.hist2d(x_grid, y_grid, bins=[x_bins, y_bins], weights=bz_grid,
                    cmin=min_by, cmax=max_by, cmap=self.cmap, vmin=-cmax, vmax=cmax)
        axes.set_xlabel("x [m]")
        axes.set_ylabel("y [m]")
        axes.set_title("$B_{z}$ [T]")
        figure.colorbar(hist[3])
        fig_fname = os.path.join(self.plot_dir, "scaling_ffa_map_cart.png")
        figure.savefig(fig_fname)
        print("Generated cartesian field map in", fig_fname)
        return figure

    def plot_tracks_cartesian(self, axes):
        """
        Plot tracks in cartesian coordinates
        - axes: Axes object to draw on
        """
        for pid, track in self.track_orbit_dict.items():
            pos_x = track["x"]
            pos_y = track["y"]
            axes.plot(pos_x, pos_y)

    def plot_tracks_cylindrical(self, axes):
        """
        Plot tracks in cylindrical coordinates
        - axes: Axes object to draw on
        If phi changes by more than 180 degrees in a single step, the track is
        assumed to have gone through an entire turn (360 <-> 0) and a new plot
        is started
        """
        for pid, track in self.track_orbit_dict.items():
            # find the list of indices where the track loops
            delta_list = [i+1 for i, phi1 in enumerate(track["phi"][1:]) \
                                    if abs(track["phi"][i] - phi1) > 180]
            phi_list_of_lists = numpy.split(track["phi"], delta_list)
            r_list_of_lists = numpy.split(track["r"], delta_list)
            for i in range(len(phi_list_of_lists)):
                axes.plot(phi_list_of_lists[i], r_list_of_lists[i])


    def oned_field_map(self, radius, axes = None):
        """
        Make a one dimensional field map along a line of constant radius
        - radius: line to draw the fields on
        - axes: Axes object to draw on
        """
        bz_points = []

        for phi in self.phi_points:
            point = (radius*math.cos(math.radians(phi)),
                     radius*math.sin(math.radians(phi)),
                     0,
                     0)
            value = pyopal.objects.field.get_field_value(*point)
            bz_points.append(value[3])

        if not axes:
            figure = matplotlib.pyplot.figure()
            axes = figure.add_subplot(1, 1, 1)

        axes.plot(self.phi_points, bz_points)
        #for contour in self.spiral_contours:
        #    self.draw_azimuthal_contour(radius, axes, contour)

        axes.set_xlabel("$\\phi$ [deg]")
        axes.set_ylabel("$B_z$ [T]")
        return axes.figure, bz_points

    @classmethod
    def draw_azimuthal_contour(cls, radius, axes, contour):
        """
        Draw an azimuthal contour (including spiral angle if required)
        - radius: r0 for spiral angle calculations
        - axes: Axes object to draw on
        - contour: contour object
        """
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        tan_d = math.tan(contour["spiral_angle"])
        phi = contour["phi0"] + tan_d*math.log(radius/contour["r0"])
        phi = math.degrees(phi)
        axes.plot([phi, phi], ylim,
                  linestyle=contour["linestyle"],
                  color=contour["colour"])
        axes.text(phi, ylim[-1],
                  contour["label"], va="top",
                  rotation="vertical",
                  color=contour["colour"])
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)


    def get_derivative(self, var1, var2, pos_x, pos_y, pos_z, time):
        """
        Calculate a numerical derivative of the field in cartesian coordinates
        - var1: variable in the numerator, should be one of self.position_variables
        - var2: variable in the denominator (range 1 to 4)
        """
        pos_vec = [pos_x, pos_y, pos_z, time]
        var2 = self.position_variables.index(var2)
        pos_vec[var2] += self.delta_x
        field_plus = pyopal.objects.field.get_field_value(*pos_vec)[var1]
        pos_vec[var2] -= 2*self.delta_x
        field_minus = pyopal.objects.field.get_field_value(*pos_vec)[var1]
        derivative = (field_plus-field_minus)/2.0/self.delta_x
        return derivative

    def get_derivative_int(self, var1, var2, pos_x, pos_y, pos_z, time):
        """
        Calculate a numerical derivative of the field in cartesian coordinates
        - var1: should be an integer range 1 to 6 inclusive for bx, by, bz, ex, ey, ez
        - var2: should be an integer range 0 to 3 inclusive for x, y, z, t
        - x, y, z, t: position at which the derivative should be calculated
        """
        raise RuntimeError("Not implemented")

    def get_div_b(self, pos_x, pos_y, pos_z, time):
        """
        Calculate Div B
        """
        div_b = self.get_derivative("bx", "x", pos_x, pos_y, pos_z, time) + \
                self.get_derivative("by", "y", pos_x, pos_y, pos_z, time) + \
                self.get_derivative("bz", "z", pos_x, pos_y, pos_z, time)
        return div_b

    def get_curl_b(self, pos_x, pos_y, pos_z, time):
        """
        Calculate Curl B
        """
        curl_b = [
            self.get_derivative("by", "z", pos_x, pos_y, pos_z, time) - \
            self.get_derivative("bz", "y", pos_x, pos_y, pos_z, time),
            self.get_derivative("bx", "z", pos_x, pos_y, pos_z, time) - \
            self.get_derivative("bz", "x", pos_x, pos_y, pos_z, time),
            self.get_derivative("bx", "y", pos_x, pos_y, pos_z, time) - \
            self.get_derivative("by", "x", pos_x, pos_y, pos_z, time)
        ]
        return curl_b

    default_radial_contour = {"radius":0.0,
                              "linestyle":"-",
                              "colour":"grey",
                              "oned_plot":False}
    default_spiral_contour = {
        "phi0":0.0, # the contour will pass through a point having cylindrical coordinates phi0, r0
        "r0":0, # the contour will pass through a point having cylindrical coordinates phi0, r0
        "spiral_angle":0, # spiral angle parameter in degrees
        "linestyle":"-", # matplotlib line style
        "colour":"grey", # matplotlib line colour
    }
    position_variables = ["x", "y", "z", "t"]
    field_variables = ["out_of_bounds", "bx", "by", "bz", "ex", "ey", "ez"]
