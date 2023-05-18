#!/usr/bin/env python
__author__ = 'Sergei F. Kliver'

from pathlib import Path
import numpy as np

try:
    argrelextrema_imported = True
    from scipy.signal import argrelextrema
except:
    print("WARNING!!! Impossible to import argrelextrema function form scipy.signal. Related functionality is disabled.")
    argrelextrema_imported = False

import matplotlib.pyplot as plt

from matplotlib.transforms import Bbox, TransformedBbox, blended_transform_factory
from mpl_toolkits.axes_grid1.inset_locator import BboxPatch, BboxConnector, BboxConnectorPatch

from KRATER.Tools.Abstract import Tool

from KRATER.Tools.Kmers import GenomeScope2


class JellyfishRoutines(Tool):
    """
    Several subcommands are not implemented: query, qhisto, qdump, qmerge, cite
    Not all options were implemented for count subcommand

    """
    def __init__(self, path="", max_threads=4):
        Tool.__init__(self, "jellyfish", path=path, max_threads=max_threads)

    def draw_kmer_distribution(self, histo_file_list, kmer_length, output_prefix, label_list=None,
                               output_formats=("svg", "png"), logbase=10, non_log_low_limit=5, non_log_high_limit=100,
                               order=3, mode="wrap", check_peaks_coef=10, draw_separated_pictures=False,
                               use_second_peak_for_genome_size_estimation=False, dont_show_genome_size_on_plot=False,
                               genomescope2=True, ploidy=2, initial_haploid_coverage=None,
                               genomescope_cmd="genomescope.R", show_confidence_interval=True):
        data_list = []
        parameters_list = []
        out_dir_path = Path(output_prefix).parent
        genomescope2_dir_path = out_dir_path / "genomescope2"
        stat_fd = open("%s.histo.stats" % output_prefix, "w")
        max_selected_counts = 0
        max_counts = 0

        default_labels = map(lambda ind: "S%i" % ind, range(1, len(histo_file_list) + 1))

        histo_file_listtt = [histo_file_list] if isinstance(histo_file_list, str) else histo_file_list
        label_listtt = [label_list] if isinstance(label_list, str) else label_list

        for histo_file, label in zip(histo_file_listtt, label_listtt if label_listtt else default_labels):
            bins, counts = np.loadtxt(histo_file, unpack=True)
            data_list.append((bins, counts))
            genomescope2_sample_dir_path = genomescope2_dir_path / label
            genomescope2_sample_dir_path.mkdir(parents=True, exist_ok=True)
            genomescope2_sample_stats_file = genomescope2_sample_dir_path / (label + "_stats.tsv")

            if genomescope2:
                GenomeScope2.get_genome_size(histo_file, kmer_length, genomescope2_sample_dir_path, label,
                                             ploidy=ploidy, initial_haploid_coverage=initial_haploid_coverage,
                                             draw_fitted_hist=True, testing=True, max_kmer_coverage=100000000,
                                             cmd=genomescope_cmd)
                with open(genomescope2_sample_stats_file, "r") as genomescope2_stat_fd:
                    genomescope2_stat_fd.readline()
                    genomescope2_haplome_coverage, genomescope2_haplome_coverage_min,  \
                    genomescope2_haplome_coverage_max, genomescope2_haplome_coverage_half_conf_len = \
                        list(map(float, genomescope2_stat_fd.readline().strip().split("\t")[1:]))
                    genomescope2_genomesize, genomescope2_genomesize_min, \
                    genomescope2_genomesize_max, genomescope2_genomesize_half_conf_len = \
                        list(map(float, genomescope2_stat_fd.readline().strip().split("\t")[1:]))

            if argrelextrema_imported:

                maximums_to_show, minimums_to_show, \
                    unique_peak_borders, number_of_distinct_kmers, \
                    number_of_distinct_kmers_with_errors,\
                    total_number_of_kmers, total_number_of_kmers_with_errors, \
                    estimated_genome_size, max_estimated_genome_size, min_estimated_genome_size, \
                    estimated_genome_size_half_conf_len = self.extract_parameters_from_histo(counts, bins,
                                                                                             "{}.{}".format(output_prefix, label),
                                                                                             order=order,
                                                                                             mode=mode,
                                                                                             check_peaks_coef=check_peaks_coef,
                                                                                             use_second_peak_for_genome_size_estimation=use_second_peak_for_genome_size_estimation)
                parameters_list.append((maximums_to_show,
                                        minimums_to_show,
                                        unique_peak_borders,
                                        number_of_distinct_kmers,
                                        number_of_distinct_kmers_with_errors,
                                        total_number_of_kmers,
                                        total_number_of_kmers_with_errors,
                                        estimated_genome_size,
                                        max_estimated_genome_size,
                                        min_estimated_genome_size,
                                        estimated_genome_size_half_conf_len))

                unique_peak_width = unique_peak_borders[1] - unique_peak_borders[0] + 1

                unique_peak_borders_mean_multiplicity = self.mean_from_bins(bins[unique_peak_borders[0]: unique_peak_borders[1]+1],
                                                                            counts[unique_peak_borders[0]: unique_peak_borders[1]+1])

                fraction_of_distinct_kmers_with_errors = float(number_of_distinct_kmers_with_errors)/float(number_of_distinct_kmers)
                fraction_of_kmers_with_errors = float(total_number_of_kmers_with_errors)/float(total_number_of_kmers)
            else:
                maximums_to_show, minimums_to_show, \
                    unique_peak_borders, number_of_distinct_kmers, \
                    number_of_distinct_kmers_with_errors,\
                    total_number_of_kmers, total_number_of_kmers_with_errors, \
                    estimated_genome_size, max_estimated_genome_size, min_estimated_genome_size, \
                    estimated_genome_size_half_conf_len, \
                    unique_peak_width, \
                    unique_peak_borders_mean_multiplicity, \
                    fraction_of_distinct_kmers_with_errors, \
                    fraction_of_kmers_with_errors = [None] * 15

            general_stats = "Sample\t%s\n" % label
            general_stats += "Number of distinct kmers\t%s\n" % (str(number_of_distinct_kmers) if number_of_distinct_kmers is not None else "NA")
            general_stats += "Number of distinct kmers with errors\t%s\n" % (str(number_of_distinct_kmers_with_errors) if number_of_distinct_kmers_with_errors is not None else "NA")
            general_stats += "Fraction of distinct kmers with errors\t%s\n" % (("%.3f" % np.around(fraction_of_distinct_kmers_with_errors, decimals=3)) if fraction_of_distinct_kmers_with_errors is not None else "NA")
            general_stats += "Total number of kmers\t%s\n" % (str(total_number_of_kmers) if total_number_of_kmers is not None else "NA")
            general_stats += "Total number of kmers with errors\t%s\n" % (str(total_number_of_kmers_with_errors) if total_number_of_kmers_with_errors is not None else "NA")
            general_stats += "Fraction of kmers with errors\t%s\n" % (("%.3f" % np.around(fraction_of_kmers_with_errors, decimals=3)) if total_number_of_kmers is not None else "NA")
            general_stats += "Kmer multiplicity at first minimum\t%s\n" % (str(minimums_to_show[0][0]) if minimums_to_show else "NA")
            general_stats += "Kmer multiplicity at first maximum\t%s\n" % (str(maximums_to_show[0][0]) if maximums_to_show else "NA")
            general_stats += "Width of first peak\t%s\n" % (str(unique_peak_width) if unique_peak_width is not None else "NA")
            general_stats += "Mean kmer multiplicity in first peak\t%s\n" % (("%.2f" % np.around(unique_peak_borders_mean_multiplicity, decimals=2)) if unique_peak_borders_mean_multiplicity is not None else "NA")
            general_stats += "Estimated genome size (naive), bp\t%s\n" % ("NA" if ((estimated_genome_size is None) or (estimated_genome_size_half_conf_len is None)) else ("%i ± %i" % (estimated_genome_size, estimated_genome_size_half_conf_len)))

            if genomescope2:
                general_stats += "Estimated genome size (genomescope2), bp\t%s\n" % ("NA" if ((genomescope2_genomesize is None) or (genomescope2_genomesize_half_conf_len is None)) else ("%i ± %i" % (int(genomescope2_genomesize), int(genomescope2_genomesize_half_conf_len))))
                general_stats += "Estimated haplome coverage (genomescope2)\t%s\n" % ("NA" if ((genomescope2_haplome_coverage is None) or (genomescope2_haplome_coverage_half_conf_len is None)) else ("%.2f ± %.2f" % (genomescope2_haplome_coverage, genomescope2_haplome_coverage_half_conf_len)))

            stat_fd.write(general_stats)
            print(general_stats)
            print("\n")

            max_bin = max(bins)

            selected_counts = counts[non_log_low_limit-1:non_log_high_limit]
            selected_bins = bins[non_log_low_limit-1:non_log_high_limit]

            current_max_selected_counts = max(selected_counts)
            current_max_counts = max(counts)

            max_selected_counts = current_max_selected_counts if current_max_selected_counts > max_selected_counts else max_selected_counts
            max_counts = current_max_counts if current_max_counts > max_counts else max_counts

            if draw_separated_pictures:
                figure = plt.figure(1, figsize=(5, 5), dpi=300)
                subplot = plt.subplot(1, 1, 1)
                plt.suptitle("Distribution of %i-mers" % kmer_length, fontweight='bold')
                plt.plot(bins, counts)
                plt.xlim(xmin=1, xmax=max_bin)
                plt.ylim(ymin=1, ymax=max_counts)
                plt.xlabel("Multiplicity")
                plt.ylabel("Number of distinct %s-mers" % kmer_length)
                subplot.set_yscale('log', base=logbase)
                subplot.set_xscale('log', base=logbase)

                figure = plt.figure(2, figsize=(5, 5), dpi=300)
                subplot = plt.subplot(1, 1, 1)
                plt.suptitle("Distribution of %s-mers" % kmer_length, fontweight='bold')
                plt.plot(selected_bins, selected_counts)

                plt.xlabel("Multiplicity")
                plt.ylabel("Number of distinct %s-mers" % kmer_length)
                plt.xlim(xmin=non_log_low_limit, xmax=non_log_high_limit)
                plt.ylim(ymin=0, ymax=max_selected_counts)

            final_estimated_genome_size = genomescope2_genomesize if genomescope2 else estimated_genome_size
            final_estimated_genome_size_half_conf_len = genomescope2_genomesize_half_conf_len if genomescope2 else estimated_genome_size_half_conf_len

            if final_estimated_genome_size is None:
                legend = "NA"
            else:
                size_in_terabases = float(final_estimated_genome_size) / float(10 ** 12)
                size_half_conf_len_in_terabases = float(final_estimated_genome_size_half_conf_len) / float(10 ** 12)
                size_in_gigabases = float(final_estimated_genome_size) / float(10 ** 9)
                size_half_conf_len_in_gigabases = float(final_estimated_genome_size_half_conf_len) / float(10 ** 9)
                size_in_megabases = float(final_estimated_genome_size) / float(10 ** 6)
                size_half_conf_len_in_megabases = float(final_estimated_genome_size_half_conf_len) / float(10 ** 6)
                size_in_kilobases = float(final_estimated_genome_size) / float(10 ** 3)
                size_half_conf_len_in_kilobases = float(final_estimated_genome_size_half_conf_len) / float(10 ** 3)

                if size_in_terabases > 1:
                    legend = "%s: %.2f" % (label, size_in_terabases)
                    if show_confidence_interval:
                        legend += "±%.2f" % size_half_conf_len_in_terabases
                    legend += " Tbp"
                elif size_in_gigabases > 1:
                    legend = "%s: %.2f" % (label, size_in_gigabases)
                    if show_confidence_interval:
                        legend += "±%.2f" % size_half_conf_len_in_gigabases
                    legend += " Gbp"
                elif size_in_megabases > 1:
                    legend = "%s: %.2f" % (label, size_in_megabases)
                    if show_confidence_interval:
                        legend += "±%.2f" % size_half_conf_len_in_megabases
                    legend += " Mbp"
                else:
                    legend = "%s: %.2f" % (label, size_in_kilobases)
                    if show_confidence_interval:
                        legend += "±%.2f" % size_half_conf_len_in_kilobases
                    legend += " Kbp"

            for index in range(3, 5):
                figure = plt.figure(index, figsize=(5, 10))
                subplot_list = []
                for i, b, c in zip([1, 2], [bins, selected_bins], [counts, selected_counts]):
                    subplot_list.append(plt.subplot(2, 1, i))
                    plt.suptitle("Distribution of %s-mers" % kmer_length, fontweight='bold', fontsize=13)
                    plt.plot(b, c, label=legend) if (i == 1) and (not dont_show_genome_size_on_plot) else plt.plot(b, c)

                    if i == 1:
                        plt.legend(loc="best")

                    if index == 4:
                        for minimum in minimums_to_show:
                            plt.plot([minimum[0], minimum[0]], [0, minimum[1]], 'r--', lw=1)
                        for maximum in maximums_to_show:
                            plt.plot([maximum[0], maximum[0]], [0, maximum[1]], 'g--', lw=1)

                    plt.ylabel("Number of distinct %s-mers" % kmer_length, fontsize=13)

                    if i == 1:
                        subplot_list[0].set_yscale('log', base=logbase)
                        subplot_list[0].set_xscale('log', base=logbase)
                        plt.xlim(xmin=1, xmax=max_bin)
                        plt.ylim(ymin=1, ymax=max_counts)
                    elif i == 2:
                        plt.ylim(ymin=0, ymax=max_selected_counts)
                        plt.xlim(xmin=non_log_low_limit, xmax=non_log_high_limit)
                        plt.xlabel("Multiplicity", fontsize=15)

                self.zoom_effect(subplot_list[0], subplot_list[1], non_log_low_limit, non_log_high_limit)
                plt.subplots_adjust(hspace=0.12, wspace=0.05, top=0.95, bottom=0.05, left=0.16, right=0.95)

        if draw_separated_pictures:
            for index, name_template in ((1, "%s.logscale.%s"), (2, "%s.no_logscale.%s")):
                plt.figure(index)
                for extension in output_formats:
                    plt.savefig(name_template % (output_prefix, extension))
                plt.close()

        for index in range(3, 5):
            plt.figure(index)
            for extension in output_formats:
                plt.savefig("%s.%s" % ("%s.peaks_and_gaps" % output_prefix if index == 4 else output_prefix, extension))

            plt.close()

    @staticmethod
    def find_peak_indexes_from_histo(counts, order=3, mode="wrap"):
        """
        order:
            How many points on each side to use for the comparison to consider comparator(n, n+x) to be True.
        mode:
            How the edges of the vector are treated. 'wrap' (wrap around) or 'clip' (treat overflow as the same
            as the last
        """
        local_maximums_idx = argrelextrema(counts, np.greater, order=order, mode=mode)[0]
        local_minimums_idx = argrelextrema(counts, np.less, order=order, mode=mode)[0]

        return local_minimums_idx, local_maximums_idx

    def extract_parameters_from_histo(self, counts, bins, output_prefix, order=3, mode="wrap", check_peaks_coef=10,
                                      use_second_peak_for_genome_size_estimation=False):
        """
        check_peaks_coef:
            histogram is checked for presence of additional peaks in range [first_unique_peak, check_peaks_coef*first_unique_peak]
        """

        local_maximums_idx = argrelextrema(counts, np.greater, order=order, mode=mode)[0]
        local_minimums_idx = argrelextrema(counts, np.less, order=order, mode=mode)[0]
        if output_prefix is not None:
            with open("%s.local_maximums" % output_prefix, "w") as out_fd:
                out_fd.write("#multiplicity\tnumber_of_kmers\n")
                for idx in local_maximums_idx:
                    out_fd.write("%i\t%i\n" % (bins[idx], counts[idx]))

            with open("%s.local_minimums" % output_prefix, "w") as out_fd:
                out_fd.write("#multiplicity\tnumber_of_kmers\n")
                for idx in local_minimums_idx:
                    out_fd.write("%i\t%i\n" % (bins[idx], counts[idx]))

        first_unique_peak_idx_idx = 0 if local_maximums_idx[0] != 0 else 1
        # second_unique_peak_idx_idx = 1 if local_maximums_idx[1] != 0 else 2
        first_unique_peak_coverage = bins[local_maximums_idx[first_unique_peak_idx_idx]]
        second_unique_peak_coverage = bins[local_maximums_idx[first_unique_peak_idx_idx+1]]

        max_checked_coverage = check_peaks_coef * first_unique_peak_coverage
        peaks_in_checked_area_idx = [local_maximums_idx[first_unique_peak_idx_idx]]
        minimums_in_checked_area_idx = []

        for i in range(first_unique_peak_idx_idx+1, len(local_maximums_idx)):
            if bins[local_maximums_idx[i]] <= max_checked_coverage:
                peaks_in_checked_area_idx.append(local_maximums_idx[i])
            else:
                break

        for minimum_index in local_minimums_idx:
            if bins[minimum_index] <= max_checked_coverage:
                minimums_in_checked_area_idx.append(minimum_index)

        if len(peaks_in_checked_area_idx) > 1:
            print("WARNING! Additional k-mer peaks were detected with multiplicity (%i, %i]" % (first_unique_peak_coverage,
                                                                                                max_checked_coverage))

        nearest_value_to_first_min_idx = self.find_nearest_scalar(counts[local_maximums_idx[first_unique_peak_idx_idx]:],
                                                                          counts[local_minimums_idx[0]]) + local_maximums_idx[first_unique_peak_idx_idx]

        number_of_distinct_kmers = sum(counts)
        number_of_distinct_kmers_with_errors = sum(counts[0:local_minimums_idx[0]])
        total_number_of_kmers = sum(np.multiply(counts, bins))
        total_number_of_kmers_with_errors = sum(np.multiply(counts[0:local_minimums_idx[0]],
                                                            bins[0:local_minimums_idx[0]]))

        maximums_to_show = [(bins[i], counts[i]) for i in peaks_in_checked_area_idx]
        minimums_to_show = [(bins[i], counts[i]) for i in minimums_in_checked_area_idx]

        estimated_genome_size = 0

        if minimums_to_show:
            for i in range(int(minimums_to_show[0][0]), len(counts)):
                estimated_genome_size += counts[i] * bins[i]
        else:
            print("WARNING: Minimum between error peak and unique peak was not found. "
                  "Estimation of genome size is impossible ")
            estimated_genome_size = None

        genome_coverage_peak = second_unique_peak_coverage if use_second_peak_for_genome_size_estimation else first_unique_peak_coverage
        naive_estimated_genome_size = estimated_genome_size/genome_coverage_peak if estimated_genome_size else None

        max_naive_estimated_genome_size = estimated_genome_size / (genome_coverage_peak - 1) if estimated_genome_size else None
        min_naive_estimated_genome_size = estimated_genome_size / (genome_coverage_peak + 1) if estimated_genome_size else None
        if estimated_genome_size:
            naive_estimated_genome_size_half_conf_len = max(
                                                            [max_naive_estimated_genome_size - naive_estimated_genome_size,
                                                             naive_estimated_genome_size - min_naive_estimated_genome_size])
        else:
            naive_estimated_genome_size_half_conf_len = None
        return maximums_to_show, minimums_to_show, \
               (local_minimums_idx[0], nearest_value_to_first_min_idx), \
               number_of_distinct_kmers, number_of_distinct_kmers_with_errors, \
               total_number_of_kmers, total_number_of_kmers_with_errors, \
               naive_estimated_genome_size, max_naive_estimated_genome_size, min_naive_estimated_genome_size, \
               naive_estimated_genome_size_half_conf_len

    @staticmethod
    def connect_bbox(bbox1, bbox2,
                     loc1a, loc2a, loc1b, loc2b,
                     prop_lines, prop_patches=None):
        if prop_patches is None:
            prop_patches = prop_lines.copy()
            prop_patches["alpha"] = prop_patches.get("alpha", 1) * 0.2

        c1 = BboxConnector(bbox1, bbox2, loc1=loc1a, loc2=loc2a, **prop_lines)
        c1.set_clip_on(False)
        c2 = BboxConnector(bbox1, bbox2, loc1=loc1b, loc2=loc2b, **prop_lines)
        c2.set_clip_on(False)

        bbox_patch1 = BboxPatch(bbox1, **prop_patches)
        bbox_patch2 = BboxPatch(bbox2, **prop_patches)

        p = BboxConnectorPatch(bbox1, bbox2,
                               loc1a=loc1a, loc2a=loc2a, loc1b=loc1b, loc2b=loc2b,
                               **prop_patches)
        p.set_clip_on(False)

        return c1, c2, bbox_patch1, bbox_patch2, p

    def zoom_effect(self, ax1, ax2, xmin, xmax, alpha=0.22, color="gray", **kwargs):
        """
        ax1 : the main axes
        ax2 : the zoomed axes
        (xmin,xmax) : the limits of the colored area in both plot axes.

        connect ax1 & ax2. The x-range of (xmin, xmax) in both axes will
        be marked.  The keywords parameters will be used ti create
        patches.

        """

        trans1 = blended_transform_factory(ax1.transData, ax1.transAxes)
        trans2 = blended_transform_factory(ax2.transData, ax2.transAxes)

        bbox = Bbox.from_extents(xmin, 0, xmax, 1)

        mybbox1 = TransformedBbox(bbox, trans1)
        mybbox2 = TransformedBbox(bbox, trans2)

        prop_patches = kwargs.copy()
        prop_patches["ec"] = "none"
        prop_patches["alpha"] = alpha
        prop_patches["color"] = color

        c1, c2, bbox_patch1, bbox_patch2, p = self.connect_bbox(mybbox1, mybbox2,
                                                                loc1a=3, loc2a=2, loc1b=4, loc2b=1,
                                                                prop_lines=kwargs, prop_patches=prop_patches)

        ax1.add_patch(bbox_patch1)
        ax2.add_patch(bbox_patch2)
        ax2.add_patch(c1)
        ax2.add_patch(c2)
        ax2.add_patch(p)

        return c1, c2, bbox_patch1, bbox_patch2, p

    @staticmethod
    def find_nearest_scalar(array, value, mode="index"):
        # works for 1d array
        idx = (np.abs(array - value)).argmin()
        return idx if mode == "index" else array.flat[idx]

    @staticmethod
    def find_nearest_vector(array, value, mode="index"):
        idx = np.array([np.linalg.norm(x + y) for (x, y) in array - value]).argmin()
        return idx if mode == "index" else array[idx]

    @staticmethod
    def mean_from_bins(bins, counts):
        return sum(np.multiply(bins, counts)) / sum(counts)

    def variance_from_bins(self, bins, counts, mean=None):
        mean_value = mean if mean else self.mean_from_bins(bins, counts)
        deviation = bins - mean_value
        variance = np.sum(np.multiply(np.power(deviation, 2), counts)) / sum(counts)
        return variance

    def std_from_bins(self, bins, counts, mean=None):
        return np.sqrt(self.variance_from_bins(bins, counts, mean=mean))