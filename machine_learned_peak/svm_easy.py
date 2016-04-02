#!/usr/bin/env python

import sys
import os
from subprocess import *

if len(sys.argv) <= 1:
	print('Usage: {0} training_file [testing_file]'.format(sys.argv[0]))
	raise SystemExit

# svm, grid, and gnuplot executable files

output_dir = "svm_models/"

is_win32 = (sys.platform == 'win32')
if not is_win32:
	svmscale_exe = "/Users/yichao/tools/libsvm-3.21/svm-scale"
	svmtrain_exe = "/Users/yichao/tools/libsvm-3.21/svm-train"
	svmpredict_exe = "/Users/yichao/tools/libsvm-3.21/svm-predict"
	grid_py = "/Users/yichao/tools/libsvm-3.21/tools/grid.py"
	gnuplot_exe = "/usr/bin/gnuplot"
else:
        # example for windows
	svmscale_exe = r"..\windows\svm-scale.exe"
	svmtrain_exe = r"..\windows\svm-train.exe"
	svmpredict_exe = r"..\windows\svm-predict.exe"
	gnuplot_exe = r"c:\tmp\gnuplot\binary\pgnuplot.exe"
	grid_py = r".\grid.py"

print svmscale_exe
assert os.path.exists(svmscale_exe),"svm-scale executable not found"
assert os.path.exists(svmtrain_exe),"svm-train executable not found"
assert os.path.exists(svmpredict_exe),"svm-predict executable not found"
# assert os.path.exists(gnuplot_exe),"gnuplot executable not found"
assert os.path.exists(grid_py),"grid.py not found"

if len(sys.argv) > 3:
	svm_type = sys.argv[3]
else:
	svm_type = 0

if len(sys.argv) > 4:
	kernel_type = sys.argv[4]
else:
	kernel_type = 2

if len(sys.argv) > 5:
	w0 = sys.argv[5]
else:
	w0 = 1

if len(sys.argv) > 6:
	w1 = sys.argv[6]
else:
	w1 = 1

train_pathname = sys.argv[1]
assert os.path.exists(train_pathname),"training file not found"
file_name = os.path.split(train_pathname)[1]
scaled_file = '{0}{1}.svm{2}.kernel{3}.w{4}{5}.scale'.format(output_dir, file_name, svm_type, kernel_type, w0, w1)
model_file = '{0}{1}.svm{2}.kernel{3}.w{4}{5}.model'.format(output_dir, file_name, svm_type, kernel_type, w0, w1)
range_file = '{0}{1}.svm{2}.kernel{3}.w{4}{5}.range'.format(output_dir, file_name, svm_type, kernel_type, w0, w1)

if len(sys.argv) > 2:
	test_pathname = sys.argv[2]
	test_file_name = os.path.split(test_pathname)[1]
	assert os.path.exists(test_pathname),"testing file not found"
	scaled_test_file = '{0}{1}.svm{2}.kernel{3}.w{4}{5}.scale'.format(output_dir, test_file_name, svm_type, kernel_type, w0, w1)
	predict_test_file = '{0}{1}.svm{2}.kernel{3}.w{4}{5}.predict'.format(output_dir, test_file_name, svm_type, kernel_type, w0, w1)



cmd = '{0} -s "{1}" "{2}" > "{3}"'.format(svmscale_exe, range_file, train_pathname, scaled_file)
print('Scaling training data...')
Popen(cmd, shell = True, stdout = PIPE).communicate()

# cmd = '{0} -svmtrain "{1}" -gnuplot "{2}" "{3}"'.format(grid_py, svmtrain_exe, gnuplot_exe, scaled_file)
cmd = '{0} -svmtrain "{1}" -gnuplot "null" -out "null" "{2}"'.format(grid_py, svmtrain_exe, scaled_file)
print('Cross validation...')
f = Popen(cmd, shell = True, stdout = PIPE).stdout

line = ''
while True:
	last_line = line
	line = f.readline()
	if not line: break
c,g,rate = map(float,last_line.split())

print('Best c={0}, g={1} CV rate={2}'.format(c,g,rate))

cmd = '{0} -c {1} -g {2} -s {3} -k {4} -w0 {5} -w1 {6} "{7}" "{8}"'.format(svmtrain_exe, c, g, svm_type, kernel_type, w0, w1, scaled_file, model_file)
print('Training...')
Popen(cmd, shell = True, stdout = PIPE).communicate()

print('Output model: {0}'.format(model_file))
if len(sys.argv) > 2:
	cmd = '{0} -r "{1}" "{2}" > "{3}"'.format(svmscale_exe, range_file, test_pathname, scaled_test_file)
	print('Scaling testing data...')
	Popen(cmd, shell = True, stdout = PIPE).communicate()

	cmd = '{0} "{1}" "{2}" "{3}"'.format(svmpredict_exe, scaled_test_file, model_file, predict_test_file)
	print('Testing...')
	Popen(cmd, shell = True).communicate()

	print('Output prediction: {0}'.format(predict_test_file))
