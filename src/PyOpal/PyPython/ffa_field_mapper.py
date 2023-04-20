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


import os
import math

import matplotlib
import matplotlib.pyplot

import pyopal.objects.field


class FFAFieldMapper(object):
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

        self.plot_dir = os.getcwd()        

    @classmethod
    def binner(self, a_list):
        return [a_list[0]+(a_list[1]-a_list[0])*(i-0.5) \
                                                  for i in range(len(a_list)+1)]

    def load_tracks(self, track_orbit):
        with open(track_orbit) as fin:
            for line in range(2):
                print(fin.readline())

    def gen_cmap(self, bz_grid):
        min_bz = min(bz_grid)
        max_bz = max(bz_grid)
        cmax = max(abs(min_bz), abs(max_bz))
        return min_bz, max_bz, cmax

    def field_map_cylindrical(self, axes = None):
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

    def draw_cylindrical_radial_contour(self, axes, contour):
        print("Plotting cylindrical radial contour", contour)
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        axes.plot(xlim, [contour["radius"]]*2, linestyle=contour["linestyle"], color=contour["colour"])
        axes.text(xlim[-1], contour["radius"], contour["label"], horizontalalignment = "right", va="top", color=contour["colour"])
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

    def draw_cylindrical_spiral_contour(self, axes, contour):
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        tan_d = math.tan(contour["spiral_angle"])
        phi_points = [contour["phi0"] + tan_d*math.log(r/contour["r0"]) for r in self.r_points]
        phi_points = [math.degrees(phi) for phi in phi_points]
        axes.plot(phi_points, self.r_points, linestyle=contour["linestyle"], color=contour["colour"])
        axes.text(phi_points[-1], self.r_points[-1], contour["label"], va = "top", rotation="vertical", color=contour["colour"])
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

    def field_map_cartesian(self, axes = None):
        x_grid = []
        y_grid = []
        bz_grid = []
        for x in self.x_points:
            for y in self.y_points:
                x_grid.append(x)
                y_grid.append(y)
                point = (x, y, 0, 0)
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

    def oned_field_map(self, radius, axes = None):
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

    def draw_azimuthal_contour(self, radius, axes, contour):
        xlim = axes.get_xlim()
        ylim = axes.get_ylim()
        tan_d = math.tan(contour["spiral_angle"])
        phi = contour["phi0"] + tan_d*math.log(radius/contour["r0"])
        phi = math.degrees(phi)
        axes.plot([phi, phi], ylim, linestyle=contour["linestyle"], color=contour["colour"])
        axes.text(phi, ylim[-1], contour["label"], va = "top", rotation="vertical", color=contour["colour"])
        axes.set_xlim(xlim)
        axes.set_ylim(ylim)

    default_radial_contour = {"radius":0.0, "linestyle":"-", "colour":"grey", "oned_plot":False}
    default_spiral_contour = {"phi0":0.0, "r0":0, "tan_delta":0, "line_style":"-", "colour":"grey", "oned_plot":False}
