#!/usr/bin/env python

import os, sys, math, random, re, math
import list_data

input_dir = "./processed/"
output_dir = "./processed/"


# train_peak_idx = sorted(random.sample(xrange(20), 5))
# others = set(xrange(20)) - set(train_peak_idx)

# print train_peak_idx
# print others
# exit()

fold          = 5  ## x-fold cross-validation
sample        = 5
num_neighbors = 10

peak_list = [];
nonpeak_list = [];

for ei in xrange(1,24):
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


num_test_peaks  = int(math.ceil(num_peaks / fold))
num_train_peaks = int(num_peaks - num_test_peaks)

# num_test_nonpeaks  = int(math.ceil(num_nonpeaks / fold))
# num_train_nonpeaks = int(num_nonpeaks - num_test_nonpeaks)
num_train_nonpeaks = 10 * num_train_peaks
num_test_nonpeaks  = num_nonpeaks - num_train_nonpeaks


for fi in xrange(0,fold):

  train_peak_idx = sorted(random.sample(xrange(num_peaks), num_train_peaks))
  train_nonpeak_idx = sorted(random.sample(xrange(num_nonpeaks), num_train_nonpeaks))

  test_peak_idx = list(set(xrange(num_peaks)) - set(train_peak_idx))
  test_nonpeak_idx = list(set(xrange(num_nonpeaks)) - set(train_nonpeak_idx))

  train_peak_list = [peak_list[i] for i in train_peak_idx]
  train_nonpeak_list = [nonpeak_list[i] for i in train_nonpeak_idx]
  train_list = train_peak_list + train_nonpeak_list
  print "    fold%d: train len=%d" % (fi, len(train_list))

  test_peak_list = [peak_list[i] for i in test_peak_idx]
  test_nonpeak_list = [nonpeak_list[i] for i in test_nonpeak_idx]
  test_list = test_peak_list + test_nonpeak_list
  print "    fold%d: test len=%d" % (fi, len(test_list))


  list_data.store_data("%srx.3.all.test%d.txt" % (output_dir, fi), test_list)
  list_data.store_data("%srx.3.all.train%d.txt" % (output_dir, fi), train_list)


  # os.system("python svm_easy.py %srx.3.all.train%d.txt %srx.3.all.test%d.txt" % (output_dir, fi, output_dir, fi))

#
