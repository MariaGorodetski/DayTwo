"""
Spyder Editor

This is a temporary script file.
"""
from __future__ import division
from Bio import SeqIO
import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import math
from collections import defaultdict
import argparse


class FastaAnalyze:

    def __init__(self, input, output):
        #self._input_file = "C:\\Users\\mgorodetski\\Desktop\\daytwoHomework\\FD1179_1.fastq\\FD1179_1.fastq"
        #self._output_file = "C:\Users\mgorodetski\Desktop\daytwoHomework\output"
        self._input_file = input
        self._output_file = output
        self._reads_quality_score_ave = {}
        self._seq_per_loc = defaultdict(list)
        self._count_reads = 0
        self._count_basepairs = 0
        self._cg_count_per_location_dict = defaultdict(int)
        self._total_count_per_location_dict = defaultdict(int)
        self._bp_quality_score_ave = {}

    @staticmethod
    def _convert():
        convert = lambda x: ord(x) - 33
        return convert

    def _get_record(self, fastq_seq):
        fq_id = fastq_seq.next().strip("\n")
        fq_seq = fastq_seq.next().strip("\n")
        fq_sep = fastq_seq.next().strip("\n")
        fq_scores = map(self._convert(), list(fastq_seq.next().strip("\n")))
        return fq_id, fq_seq, fq_scores

    def run_q_1_3(self):
        print "Start parsing fastq file...(q_1_3)"
        count_basepairs = 0
        count_reads = 0
        fastq_seq = open(self._input_file)

        score_dict = {}
        while True:
            try:
                fq_id, fq_seq, fq_scores = self._get_record(fastq_seq)
                if count_basepairs >= 1000000:
                    break
                count_reads += 1
                count_basepairs += len(fq_scores)
                self._reads_quality_score_ave[fq_id] = np.mean(fq_scores)
                score_dict[count_reads] = fq_scores
            except StopIteration:
                break

        df = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in score_dict.iteritems()]))

        self._bp_quality_score_ave = df.mean(1)  # average by location
        self._build_plot_score_per_read()
        self._build_plot_score_per_location()
        self._calc_and_print_error_probabilities()

    @staticmethod
    def _calc_probability(quality):
        return math.pow(10, (-quality / 10))

    def _calc_and_print_error_probabilities(self):
        first_error_prob = self._calc_probability(self._bp_quality_score_ave[0])
        tenth_error_prob = self._calc_probability(self._bp_quality_score_ave[9])
        print "First location error probability:", first_error_prob
        print "Tenth location error probability:", tenth_error_prob

    def _build_plot_score_per_read(self):
        # plot #1 #
        value_set = self._reads_quality_score_ave.values()
        plt.title("Average Quality Score")
        plt.xlabel("Scores")
        plt.ylabel("Number of Reads")
        plt.hist(value_set)
        plt.savefig(self._output_file + os.sep + "average_scores_q1.png")
        plt.close()

    def _build_plot_score_per_location(self):
        # plot #2 #
        value_set = self._bp_quality_score_ave
        plt.title("Average Quality Score per Basepair")
        plt.xlabel("Locations")
        plt.ylabel("Average Score")
        plt.plot(range(len(value_set)), value_set, 'bo')
        plt.savefig(self._output_file + os.sep + "average_scores_bp_q1.png")
        plt.close()

    def run_q_1_4(self):

        print "Start parsing fastq file...(q_1_4)"
        count_basepairs = 0
        count_reads = 0

        fastq_seq = SeqIO.parse(self._input_file, 'fastq')

        for record in fastq_seq:

            # fq_id, fq_seq = record.id, record.seq
            if count_reads >= 1000000:
                break
            count_reads += 1
            self._update_gc_location_dict(record.seq)

        percentage_rate = [(self._cg_count_per_location_dict[i] / self._total_count_per_location_dict[i]) for i in
                           range(len(self._total_count_per_location_dict))]
        self._build_plot_gc_per_location(percentage_rate)

    def _update_gc_location_dict(self, seq_list):
        for i in range(len(seq_list)):
            obj = seq_list[i]
            if obj.lower() == 'g' or obj.lower() == 'c':
                self._cg_count_per_location_dict[i] += 1
            self._total_count_per_location_dict[i] += 1

    def _build_plot_gc_per_location(self, percentage_rate):
        # plot #3 #
        value_set = percentage_rate
        plt.title("GC rich per Basepair")
        plt.xlabel("Locations")
        plt.ylabel("GC rich percentage")
        plt.plot(range(len(value_set)), value_set, 'bo')
        plt.savefig(self._output_file + os.sep + "gc_content_bp_q1.png")
        plt.close()


if __name__ == "__main__":
    print "Start running :)"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    fa = FastaAnalyze(args.input, args.output)
    fa.run_q_1_3()
    fa.run_q_1_4()

    print "Bey bey"