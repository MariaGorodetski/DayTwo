from __future__ import division
import os
import matplotlib.pyplot as plt
import re
import argparse

class SamAnalyze:

    def __init__(self, input_file, output_dir):
        #self._input = "C:\\Users\\mgorodetski\\Desktop\\daytwoHomework\\sam\\r2g1.sam"
        #self._input = "C:\\Users\\mgorodetski\\Desktop\\daytwoHomework\\sam\\r2gs12.sam"
        #self._input = "C:\\Users\\mgorodetski\\Desktop\\daytwoHomework\\sam\\r2gs12_all.sam"
        #self._output = "C:\Users\mgorodetski\Desktop\daytwoHomework\output"
        self._input = input_file
        self._output = output_dir
        self._sam_alignment_dict = {}
        self._reads_per_3k_position_dict = {}

    def _build_location_dict_for_one_ref(self, key=None, ref_length=0):
        self._reads_per_3k_position_dict = {}
        sam_f = open(self._input, 'r')
        count = 0
        for line in sam_f.readlines():
            if re.match("^@(HD|SQ|RG|PG)",line) is not None:
                if "LN:" in line:
                    if key is None:
                        key, ref_length = re.match(".*SN:(\\S+)\\s+LN:(\\d+)", line).groups()
                    self._reads_per_3k_position_dict = {i: 0 for i in range(int(int(ref_length) / 3000) + 1)}
                continue

            line.split()
            ref = line.split()[2]
            if key is not None and not key == ref:
                continue
            count += 1
            position = int(line.split()[3])
            position_range = int(position/3000)

            self._reads_per_3k_position_dict[position_range] += 1

        x = self._reads_per_3k_position_dict.keys()
        label = [i*3000 for i in x]
        y = self._reads_per_3k_position_dict.values()
        plt.bar(x, y)
        plt.xticks(x, label, fontsize=7, rotation=40)
        input_file_name = self._input.split(os.sep)[-1].replace(".sam", "")
        plt.savefig(self._output + os.sep + "reads_per_3k_"+str(self._fix_name_string(key)+"_"+input_file_name)+".png")
        plt.close()
        print count
        sam_f.close()

    def _fix_name_string(self, key):

        return re.match("^(\w+)", key).group(1)

    def _build_location_dict_for_several_ref(self):
        sam_f = open(self._input, 'r')

        length_per_length_dict = {}
        _reads_per_3k_position_dict = {}
        for line in sam_f.readlines():
            if re.match("^@(HD|SQ|RG|PG)", line) is not None:
                if "@SQ" in line:
                    ref_name, ref_length = re.match(".*SN:(\S+)\s+LN:(\\d+)", line).groups()
                    length_per_length_dict[ref_name] = ref_length
                    print ref_name, ref_length

                continue
            else:
                break
        sam_f.close()

        for key in length_per_length_dict.keys():
            self._build_location_dict_for_one_ref( key, length_per_length_dict[key])

    def run_q2_1(self):
        self._build_location_dict_for_one_ref()

    def run_q2_2(self):
        self._build_location_dict_for_several_ref()


if __name__ == "__main__":
    print "Start running :)"
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = parser.parse_args()
    sa = SamAnalyze(args.input, args.output)
    sa.run_q2_1()
    sa.run_q2_2()

    print "Bey bey"
