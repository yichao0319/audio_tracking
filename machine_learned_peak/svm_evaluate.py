#!/usr/bin/env python

import os, sys, math, random, re
import list_data

test_dir = "./processed/"
pred_dir = "./svm_models/"

filename = sys.argv[1]


f = open(test_dir + filename, 'r')
test_label = []
for line in f:
  line = line.strip()
  # print "\n" + line
  m = re.search('(\d )(.*)', line)
  label = int(m.group(1))
  test_label.append(label)
f.close()

f = open(pred_dir + filename + ".predict", 'r')
pred_label = []
for line in f:
  line = line.strip()
  label = int(line)
  pred_label.append(label)
f.close()


tp = 0
fp = 0
tn = 0
fn = 0
for x in xrange(0,len(pred_label)):
  if test_label[x] == 1 and pred_label[x] == 1:
    tp += 1
  elif test_label[x] == 1 and pred_label[x] == 0:
    fn += 1
  elif test_label[x] == 0 and pred_label[x] == 1:
    fp += 1
  elif test_label[x] == 0 and pred_label[x] == 0:
    tn += 1

prec = float(tp) / (tp + fp)
recall = float(tp) / (tp + fn)
f1 = 2*prec*recall / (prec + recall)

print "tp=%d, tn=%d, fp=%d, fn=%d" % (tp, tn, fp, fn)
print "precision=%.4f, recall=%.4f" % (prec, recall)
print "f1=%.4f" % (f1)
