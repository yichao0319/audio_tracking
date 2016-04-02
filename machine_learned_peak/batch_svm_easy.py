#!/usr/bin/env python

import os, sys, math, random, re, math
import list_data


for ei in xrange(0, 5):
  train_file = "./processed/rx.3.all.train%d.txt" % (ei)
  test_file  = "./processed/rx.3.all.test%d.txt" % (ei)

  for svm_type in [0, 3]:
    for kernel_type in [0, 1, 2, 3]:
      for w0 in [1]:
        for w1 in [1, 10]:
          print "#######################################"
          print ">>> python svm_easy.py %s %s %d %d %d %d" % (train_file, test_file, svm_type, kernel_type, w0, w1)
          print "#######################################"
          os.system("python svm_easy.py %s %s %d %d %d %d" % (train_file, test_file, svm_type, kernel_type, w0, w1))
