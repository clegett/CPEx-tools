#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import csv
import sys


def hapke_refl_plot(a_run_output, absolute=False, output_dir=None,
                    include_csv=False):
    refl_dict = dict()
    for wl in a_run_output.wl_list:
        if not a_run_output.wl_dict[wl].incomplete:
            refl_dict[wl] = a_run_output.wl_dict[wl].hapke_refl

    wls = list()
    refls = list()
    fig = plt.figure()
    for wl, refl in refl_dict.items():
        wls.append(float(wl))
        refls.append(float(refl))
    plt.plot(wls, refls, marker='s', markersize=5)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Reflectance')
    if absolute:
        plt.ylim(ymin=0)
    plt.grid(True)
    plt.title('{} Reflectance'.format(a_run_output.runname))
    plt.xticks(np.arange(700, 1800, 100))
    plt.show()
    if output_dir:
        outfile = '{}{}_refl'.format(output_dir, a_run_output.runname)
        if absolute:
            outfile.join('_abs')
        fig.savefig('{}.png'.format(outfile))
        if include_csv:
            try:
                with open('{}.csv'.format(outfile), 'w') as o:
                    csvwriter = csv.writer(o, delimiter=',')
                    csvwriter.writerow(['Wavelength', 'Reflectance'])
                    csvwriter.writerows([[x, y] for x, y in zip(wls, refls)])
            except IOError as e:
                sys.exit('I/O error: file {}{}_refl.csv: {}'
                         ''.format(output_dir, a_run_output.runname, e))


def s11_plots(a_run_output, output_dir=None):
    s11_dict = dict()
    for wl in a_run_output.wl_list:
        inner_s11_dict = dict()
        if not a_run_output.wl_dict[wl].incomplete:
            for angle, angle_dict in a_run_output.wl_dict[wl].s_matrix_dict.items():
                inner_s11_dict[angle] = a_run_output.wl_dict[wl].s_matrix_dict[angle]['11']
            s11_dict[wl] = inner_s11_dict

    for wl, inner_s11_dict in s11_dict.items():
        angles = list()
        s11s = list()
        for angle, s11 in inner_s11_dict.items():
            angles.append(float(angle))
            s11s.append(float(s11))

        angles, s11s = (list(t) for t in zip(*sorted(zip(angles, s11s))))

        angles = (np.pi / 180) * np.asarray(angles)
        s11s = np.asarray(s11s)
        fig = plt.figure()
        ax = plt.subplots(1, 1, subplot_kw=dict(polar=True))[1]
        plot_logpolar(ax, angles, s11s)
        ax.set_title('{}: {} nm'.format(a_run_output.runname, wl))
        if output_dir:
            outfile = '{}{}_{}_s11'.format(output_dir, a_run_output.runname, wl)
            fig.savefig('{}.png'.format(outfile))
        plt.show()


def s11_plot_all(a_run_output, output_dir=None):
    colors = {1700: '#FF0000', 1600: '#FA7500', 1500: '#F6E700', 1400: '#8EF100',
              1300: '#1CED00', 1200: '#00E951', 1100: '#00E4BB', 1000: '#009EE0',
              900: '#0034DB', 800: '#3100D7', 700: '#9300D2'}
    s11_dict = dict()
    for wl in a_run_output.wl_list:
        inner_s11_dict = dict()
        if not a_run_output.wl_dict[wl].incomplete:
            for angle, angle_dict in a_run_output.wl_dict[
                    wl].s_matrix_dict.items():
                inner_s11_dict[angle] = \
                    a_run_output.wl_dict[wl].s_matrix_dict[angle]['11']
            s11_dict[wl] = inner_s11_dict

    fig = plt.figure()
    ax = fig.add_subplot(111, polar=True)
    print(ax)

    max_max10 = None
    min_min10 = None
    for wl, inner_s11_dict in s11_dict.items():
        angles = list()
        s11s = list()
        for angle, s11 in inner_s11_dict.items():
            angles.append(float(angle))
            s11s.append(float(s11))
            angles.append(360 - float(angle))
            s11s.append(float(s11))

        angles, s11s = (list(t) for t in zip(*sorted(zip(angles, s11s))))

        angles = (np.pi / 180) * np.asarray(angles)
        s11s = np.asarray(s11s)
        this_min10 = np.log10(np.min(s11s))
        this_max10 = np.log10(np.max(s11s))
        if min_min10 is None:
            min_min10 = this_min10
        elif this_min10 < min_min10:
            min_min10 = this_min10
        if max_max10 is None:
            max_max10 = this_max10
        elif this_max10 > max_max10:
            max_max10 = this_max10
        plot_logpolar(ax, angles, s11s, min10=min_min10, max10=max_max10,
                      color=colors[wl])
    ax.set_title('S11: {}'.format(a_run_output.runname))
    ax.set_xticks(np.pi / 180. * np.linspace(0, 360, 12, endpoint=False))
    if output_dir:
        outfile = '{}{}_all_s11'.format(output_dir, a_run_output.runname)
        fig.savefig('{}.png'.format(outfile))
    plt.show()


def plot_logpolar(ax, theta, r_, bullseye=None, min10=None, max10=None,
                  **kwargs):
    if min10 is None:
        min10 = np.log10(np.min(r_))
    if max10 is None:
        max10 = np.log10(np.max(r_))
    if bullseye is None:
        bullseye = min10 - np.log10(0.5 * np.min(r_))
    r = np.log10(r_) - min10 + bullseye
    ax.plot(theta, r, **kwargs)
    l = np.arange(np.floor(min10), max10)
    ax.set_rticks(l - min10 + bullseye)
    ax.set_yticklabels(["1e%d" % x for x in l])
    if ax.get_rmax() > (max10 - min10 + bullseye):
        ax.set_rlim(0, ax.get_rmax())
    else:
        ax.set_rlim(0, max10 - min10 + bullseye)
    #ax.viewLim.x1 = np.deg2rad(180)
    #ax.viewLim.x0 = np.deg2rad(0)
    return ax

def plot_all(a_run_output, an_output_dir):
    hapke_refl_plot(a_run_output, absolute=False, output_dir=an_output_dir,
                    include_csv=True)
    hapke_refl_plot(a_run_output, output_dir=an_output_dir, absolute=True)
    s11_plots(a_run_output, an_output_dir)
    s11_plot_all(a_run_output, an_output_dir)
