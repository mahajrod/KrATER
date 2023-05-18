#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'
import math
import sys
import argparse

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from distinctipy import distinctipy

try:
    argrelextrema_imported = True
    from scipy.signal import argrelextrema
except:
    print("WARNING!!! Impossible to import argrelextrema function form scipy.signal. Related functionality is disabled.")
    argrelextrema_imported = False

from KRATER.Routines import FileRoutines


def rgb_tuple_to_hex(rgb_tuple):
    color_code = "#"
    for i in [0, 1, 2]:
        color_code += "{:02X}".format(int(255 * rgb_tuple[i]))

    return color_code


parser = argparse.ArgumentParser()

parser.add_argument("-i", "--input", action="store", dest="input", default=sys.stdin,
                    help="Input 3-column (coverage, gc and count) tab-separated file. Default: stdin")
parser.add_argument("-k", "--kmer_length", action="store", dest="kmer_length", required=True, type=int,
                    help="Length of k-mer. Required")
parser.add_argument("-l", "--lambda", action="store", dest="lambd", type=float, required=True,
                    help="Haploid coverage of the genome. Affects the right border (4 * p * l) of the linear plot. Required")
parser.add_argument("-g", "--gap_coverage", action="store", dest="gap_coverage", type=float,
                    help="Coverage at the gap between peak of errors and unique peak on the kmer histogram."
                         "If not set sizes of GC fractions will not be estimated. Default: not set.")
parser.add_argument("-p", "--ploidy", action="store", dest="ploidy", type=int, default=2,
                    help="Ploidy of the genome. Affects the right border (4 * p * l) of the linear plot. Default: 2")
parser.add_argument("-r", "--border", action="store", dest="border", type=int,
                    help="Right border of the linear plot. If set -p/--ploidy and -l/--lambda options will be ignored. "
                         " Default: not set")
parser.add_argument("-b", "--bin_width", action="store", dest="bin_width", type=int, default=1,
                    help="Width of the bins for linear plot. Default: 1")
parser.add_argument("-m", "--max_ploidy_line", action="store", dest="max_ploidy_line", type=int, default=4,
                    help="Maximum ploidy to indicate by vertical dashed line on histogram."
                         " Ploidies of 4 and higher share same color (black). Default: 4")
parser.add_argument("-o", "--output_prefix", action="store", dest="output_prefix", required=True,
                    help="Prefix of output files")

parser.add_argument("-e", "--extension_list", action="store", dest="extension_list", default=["svg", "png"],
                    type=lambda s: s.split(","),
                    help="Comma-separated list of graphical formats supported by matplotlib. Default: 'png,svg'")

args = parser.parse_args()

output_dir_path_dict = {subdir: Path(args.output_prefix).parent / subdir for subdir in ["tab"] + args.extension_list}
for subdir in output_dir_path_dict:
    FileRoutines.safe_mkdir(output_dir_path_dict[subdir])
# ----- Parameters related to th detection of peaks ----
order = 3
mode = "wrap"
check_peaks_coef = 10
# ------------------------------------------------------

close_right_border = args.boreder if args.border else int(2 * args.ploidy * args.lambd)
far_right_border = args.boreder if args.border else int(10 * args.ploidy * args.lambd)

close_x_bins = list(np.arange(1, close_right_border, args.bin_width))
close_x_bins.append(close_x_bins[-1] + args.bin_width)

far_x_bins = list(np.arange(1, far_right_border, args.bin_width))
far_x_bins.append(far_x_bins[-1] + args.bin_width)

log_x_bins = [2 ** power for power in range(0, 12)]

y_bins = np.arange(0, args.kmer_length + 1, 1)

raw_gc_coverage_count = pd.read_csv(args.input, sep="\t", header=None, index_col=None,
                                    names=["coverage", "gc", "counts"])

# ------------------------GC heatmap------------------------------
fig, ax_array = plt.subplots(3, 2, sharey=True, figsize=(12, 8), dpi=300)

number_of_rows = 3
number_of_columns = 2
min_coverage = min(raw_gc_coverage_count["coverage"])
max_coverage = max(raw_gc_coverage_count["coverage"])

for row, x_bins, x_scale in zip(range(0, number_of_rows),
                                [close_x_bins, far_x_bins, log_x_bins],
                                ["linear", "linear", "log"]):
    for column, color_scale in zip(range(0, number_of_columns),
                                   ["linear", "log"]):
        hist, xedges, y_edges, image = ax_array[row][column].hist2d(raw_gc_coverage_count["coverage"],
                                                                    raw_gc_coverage_count["gc"],
                                                                    weights=raw_gc_coverage_count["counts"], bins=[x_bins, y_bins],
                                                                    norm=mpl.colors.LogNorm() if color_scale == "log" else None)
        ax_array[row][column].set_xlim(xmax=x_bins[-1], xmin=x_bins[0])
        fig.colorbar(image, label="Number of {0}-mers".format(args.kmer_length) if column == number_of_columns-1 else None)

        if column == 0:
            ax_array[row][column].set_ylabel("GC counts")
        if row == (number_of_rows - 1):
            ax_array[row][column].set_xlabel("Coverage")
            ax_array[row][column].set_xscale("log", base=2)
        if row == 0:
            if column == 0:
                ax_array[row][column].set_title("Linear color scale")
            if column == 1:
                ax_array[row][column].set_title("Logarithmic color scale")
        ax_array[row][column].set_xlim(xmin=min_coverage)
        ax_array[row][column].yaxis.set_major_locator(MaxNLocator(integer=True))
plt.suptitle("GC counts vs coverage(frequency) for {0}-mers".format(args.kmer_length), fontweight="bold")
plt.subplots_adjust(hspace=0.2, wspace=0.01, left=0.05, right=0.99)

for ext in args.extension_list:
    fig.savefig(output_dir_path_dict[ext] / "{0}.heatmap.{1}".format(args.output_prefix, ext))

for row in range(0, number_of_rows):
    for column in range(0, number_of_columns):
        ax_array[row][column].grid(visible=True,)
for ext in args.extension_list:
    fig.savefig(output_dir_path_dict[ext] / "{0}.heatmap.with_grid.{1}".format(args.output_prefix, ext))
# ------------------------------------------------------------------

# -------------------GC specific kmer histograms--------------------
colors = distinctipy.get_colors(args.kmer_length + 1)
color_list = list(map(rgb_tuple_to_hex, colors))

ploidy_line_color_list = ["red", "green", "blue", "black"] + (["black"] * (args.max_ploidy_line - 4))

n = int(math.sqrt(args.kmer_length + 2))
m = n
if (n^2) < (args.kmer_length + 2):
    m = n + 1
    if (m * n) < (args.kmer_length + 2):
        n += 1

fig, ax_array = plt.subplots(n + 1, m, sharey=False, sharex=True, figsize=(4 * (n+1), 4 * m), dpi=300)  # Y axis is shared manually as there is a bug with ticks afrter dsharing axis - ticks continue be linked
"""
for row in range(0, n-1): # disconnect Y axis of last subplot in the first row from other subplots
    for column in range(0, m):
        if (row == 0) and (column == m-1):
            ax_array[0][m-1].get_shared_y_axes().remove(ax_array[row][column])
            ax_array[0][m-1].get_shared_x_axes().remove(ax_array[row][column]) 
"""

for gc in range(0, args.kmer_length + 1):
    ax_array[0][0].plot(raw_gc_coverage_count["coverage"][raw_gc_coverage_count["gc"] == gc],
                        raw_gc_coverage_count["counts"][raw_gc_coverage_count["gc"] == gc], color=color_list[gc])

for ploidy, line_color in zip(range(1, args.max_ploidy_line + 1), ploidy_line_color_list):
    ax_array[0][0].axvline(ploidy * args.lambd, color=line_color, linestyle="dashed", zorder=1)
ax_array[0][0].set_xlim(xmin=min_coverage)
ax_array[0][0].set_ylabel("counts", fontsize=20)
ax_array[0][0].grid(visible=True)

last_column_with_plots_in_first_row = 0

gc_maximum_df = []
all_gc_maximum_df = {}
# ----- Maximum counts for
if argrelextrema_imported:
    for gc in range(0, args.kmer_length + 1):
        gc_df = raw_gc_coverage_count[raw_gc_coverage_count["gc"] == gc]
        local_maximums_idx_list = argrelextrema(np.array(gc_df["counts"]), np.greater, order=order, mode=mode)[0]
        local_maximums_idx_list = local_maximums_idx_list if local_maximums_idx_list[0] != 0 else local_maximums_idx_list[1:]
        local_minimums_idx_list = argrelextrema(np.array(gc_df["counts"]), np.less, order=order, mode=mode)[0]
        local_minimums_idx_list = local_minimums_idx_list if local_minimums_idx_list[0] != 0 else local_minimums_idx_list[1:]

        local_maximums_counts = gc_df["counts"].iloc[local_maximums_idx_list]
        local_maximums_coverage = gc_df["coverage"].iloc[local_maximums_idx_list]
        all_gc_maximum_df[gc] = pd.concat([local_maximums_counts, local_maximums_coverage], axis=1)
        all_gc_maximum_df[gc]["gc"] = gc
        all_gc_maximum_df[gc] = all_gc_maximum_df[gc][["coverage", "gc", "counts"]]
        gc_maximum_df.append([local_maximums_coverage.iloc[0], gc, local_maximums_counts.iloc[0]])

    gc_maximum_df = pd.DataFrame(gc_maximum_df, columns=["coverage", "gc", "counts"])
    last_column_with_plots_in_first_row += 1
    ax_array[0][1].plot(gc_maximum_df["coverage"],
                        gc_maximum_df["counts"], )
    ax_array[0][1].scatter(gc_maximum_df["coverage"],
                           gc_maximum_df["counts"],
                           color=color_list)
    for ploidy, line_color in zip(range(1, args.max_ploidy_line + 1), ploidy_line_color_list):
        ax_array[0][1].axvline(ploidy * args.lambd, color=line_color, linestyle="dashed", zorder=0)
    ax_array[0][1].grid(visible=True)
    ax_array[0][1].sharey(ax_array[0][0])
    gc_maximum_df.to_csv(output_dir_path_dict["tab"] / "{0}.gc_first_maximums.tab".format(args.output_prefix))
    all_gc_maximum_df = pd.concat(all_gc_maximum_df.values(), axis=0)
    all_gc_maximum_df.to_csv(output_dir_path_dict["tab"] / "{0}.gc_all_maximums.tab".format(args.output_prefix))
    print(gc_maximum_df)

raw_coverage_count = raw_gc_coverage_count[["coverage", "counts"]].groupby(["coverage"]).sum().reset_index(drop=False)

for gc in range(0, args.kmer_length + 1):
    ax_array[0][m-1].plot(raw_gc_coverage_count["coverage"][raw_gc_coverage_count["gc"] == gc],
                          raw_gc_coverage_count["counts"][raw_gc_coverage_count["gc"] == gc], color=color_list[gc])
    ax_array[0][m-1].plot(raw_gc_coverage_count["coverage"][raw_gc_coverage_count["gc"] == gc],
                          raw_gc_coverage_count["counts"][raw_gc_coverage_count["gc"] == gc], color=color_list[gc], label="All")
    for ploidy, line_color in zip(range(1, args.max_ploidy_line + 1), ploidy_line_color_list):
        ax_array[0][m-1].axvline(ploidy * args.lambd, color=line_color, linestyle="dashed")
ax_array[0][m-1].set_xlim(xmin=min_coverage)
ax_array[0][m-1].set_ylim(ymin=1, ymax=1.1 * max(raw_coverage_count["counts"]))
ax_array[0][m-1].set_ylabel("counts", fontsize=20)
ax_array[0][m-1].yaxis.set_major_locator(MaxNLocator())
ax_array[0][m-1].yaxis.set_tick_params(which='both', labelleft=True)
ax_array[0][m-1].grid(visible=True, axis='both')
ax_array[0][m-1].plot(raw_coverage_count["coverage"],
                      raw_coverage_count["counts"], color="black")


for column in range(last_column_with_plots_in_first_row + 1, m - 1):
    ax_array[0][column].set_axis_off()

last_column_with_plots_in_last_row = args.kmer_length % m

genome_part_dict = {}

for gc in range(0, m * n):
    row = (gc // m) + 1
    column = gc % m

    if gc > args.kmer_length:
        ax_array[row][column].set_axis_off()
    else:
        ax_array[row][column].sharey(ax_array[0][0])
        ax_array[row][column].plot(raw_gc_coverage_count["coverage"][raw_gc_coverage_count["gc"] == gc],
                                   raw_gc_coverage_count["counts"][raw_gc_coverage_count["gc"] == gc],
                                   label="GC: %0i" % gc, color=color_list[gc])
        ax_array[row][column].legend()
        ax_array[row][column].grid(visible=True,)
        if (args.gap_coverage is not None) and (args.gap_coverage > 0):
            genome_part_dict[gc] = np.sum(np.multiply(raw_gc_coverage_count["coverage"][(raw_gc_coverage_count["gc"] == gc) & (raw_gc_coverage_count["coverage"] >= args.gap_coverage)],
                                                      raw_gc_coverage_count["counts"][(raw_gc_coverage_count["gc"] == gc) & (raw_gc_coverage_count["coverage"] >= args.gap_coverage)])) / 2 / args.lambd

        for ploidy, line_color in zip(range(1, args.max_ploidy_line + 1), ploidy_line_color_list):
            ax_array[row][column].axvline(ploidy * args.lambd, color=line_color, linestyle="dashed")
        if column == 0:
            ax_array[row][column].set_ylabel("counts", fontsize=20)
        if (row == n) and (column <= last_column_with_plots_in_last_row):
            ax_array[row][column].set_xlabel("coverage", fontsize=20)
        if (row == n - 1) and (column > last_column_with_plots_in_last_row):
            ax_array[row][column].set_xlabel("coverage", fontsize=20)
            ax_array[row][column].xaxis.set_tick_params(which='both', labelbottom=True)
    ax_array[row][column].set_xlim(xmin=min_coverage)
    plt.suptitle("GC-specific counts vs coverage(frequency) for {0}-mers".format(args.kmer_length),
                 fontsize=32, fontweight="bold")

ax_array[0][0].set_ylim(ymin=1)
plt.subplots_adjust(hspace=0.2, wspace=0.2, left=0.05, right=0.95, bottom=0.05, top=0.95)

for xmax in close_right_border, far_right_border:
    ax_array[0][0].set_xlim(xmin=min_coverage, xmax=xmax)
    ax_array[0][m-1].set_xlim(xmin=min_coverage, xmax=xmax)
    for ext in args.extension_list:
        fig.savefig(output_dir_path_dict[ext] / "{0}.gc_histogramms.max{1}.{2}".format(args.output_prefix, xmax, ext))

ax_array[0][0].set_xscale("log", base=10)
ax_array[0][0].set_yscale("log", base=10)
ax_array[0][m-1].set_xscale("log", base=10)
ax_array[0][m-1].set_yscale("log", base=10)

ax_array[0][0].set_xlim(xmin=min_coverage, xmax=max_coverage)
ax_array[0][m-1].set_xlim(xmin=min_coverage, xmax=max_coverage)
for ext in args.extension_list:
    fig.savefig(output_dir_path_dict[ext] / "{0}.gc_histogramms.log.{1}".format(args.output_prefix, ext))

for power in 3, 4, 5, 6:
    ax_array[0][0].set_xlim(xmin=min_coverage, xmax=10**power)
    ax_array[0][m-1].set_xlim(xmin=min_coverage, xmax=10**power)
    for ext in args.extension_list:
        fig.savefig(output_dir_path_dict[ext] / "{0}.gc_histogramms.max{1}.log.{2}".format(args.output_prefix, 10 ** power, ext))

