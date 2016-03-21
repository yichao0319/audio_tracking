#!/usr/bin/env python

import os, sys, math, random, re, math
import list_data

input_dir = "./processed/"
output_dir = "./processed/"


num_train_traces = 3
sample           = 5
num_neighbors    = 5

peak_list = [];
nonpeak_list = [];

for ei in xrange(1,num_train_traces+1):
  filename = "rx.3.%d.samp%d.neighbor%d" % (ei, sample, num_neighbors)

  f = open(input_dir + filename + ".txt", 'r')
  for line in f:
    line = line.strip()
    # print "\n" + line
    m = re.search('(\d )(.*)', line)
    label = int(m.group(1))
    # remains = m.group(2)
    # print "label: " + str(label)
    # print "remain: " + remains

    if label == 0:
      nonpeak_list.append(line)
    elif label == 1:
      peak_list.append(line)

    # m = re.search('(\d+):(-*\d+\.*\d*e*-*\d*)( *)(.*)', remains)
    # ret = m.groups()
    # while len(ret) > 0:
    #   remains = ret[3]
    #   if remains == "":
    #     break

    #   m = re.search('(\d+):(-*\d+\.*\d*e*\d*-*\d*)( *)(.*)', remains)
    #   ret = m.groups()
  f.close()


num_peaks    = len(peak_list)
num_nonpeaks = len(nonpeak_list)
print "  num peak samples: %d" % (num_peaks)
print "  num non-peak samples: %d" % (num_nonpeaks)


train_filename = "%srx.3.train_traces.%d.txt" % (output_dir, num_train_traces)
list_data.store_data(train_filename, peak_list+nonpeak_list)

os.system("python svm_easy.py " + train_filename)

