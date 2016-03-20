#!/usr/bin/python
# -*- coding: utf-8 -*-

#import json, urllib, urllib2
import sys, os, time
import gzip


def load_data(filename):
  list_data = []
  try:
    if filename.endswith('.gz'):
      f = gzip.open(filename, 'rb')
    else:
      f = open(filename, 'r')

    print "  load data: " + filename
    tmp = f.read()
    list_data = tmp.split('\n')
    # print "    "+"\n    ".join(list_data)
    f.close()

  except IOError:
    # print IOError
    print "  file not exist: " + filename
  except:
    print "  cannot load " + filename
    pass

  return list_data


def store_data(filename, list_data):
  try:
    print '  store data: ' + filename
    list_data = sorted(list_data)
    f = open(filename, 'w+')
    for item in list_data:
      if item == "": continue
      if item is None: continue
      f.write(item.encode('utf-8').strip() + "\n")
    f.close()
  except Exception as e:
    print "  cannot store it! type: %s, msg: %s" % (type(e), e)
    pass

