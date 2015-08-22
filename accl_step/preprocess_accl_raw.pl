#!/bin/perl

##########################################
## Author: Yi-Chao Chen
##
## - input:
##
## - output:
##   1. locationHeadingTimestamp_since1970
##   2. locationHeadingX
##   3. locationHeadingY
##   4. locationHeadingZ
##   5. locationTrueHeading
##   6. locationMagneticHeading
##   7. locationHeadingAccuracy
##   8. accelerometerTimestamp_sinceReboot
##   9. accelerometerAccelerationX
##   10. accelerometerAccelerationY
##   11. accelerometerAccelerationZ
##   12. gyroTimestamp_sinceReboot
##   13. gyroRotationX
##   14. gyroRotationY
##   15. gyroRotationZ
##
## - e.g.
##   
##
##########################################

use strict;
use POSIX;
use List::Util qw(first max maxstr min minstr reduce shuffle sum);
# use lib "/u/yichao/utils/perl";
# use lib "../utils";

#############
# Debug
#############
my $DEBUG0 = 0;
my $DEBUG1 = 1;
my $DEBUG2 = 1; ## print progress
my $DEBUG3 = 1; ## print output


#############
# Constants
#############


#############
# Variables
#############
my $input_dir  = "./accl_raw";
my $output_dir = "./accl_mat";

my @filenames = ("0820.exp1.accl.walk.50m", "0820.exp2.accl.walk.55m", "0820.exp3.accl.walk.45m", "0820.exp4.accl.walk.50m", "0820.exp5.accl.walk.square");

#############
# check input
#############
if(@ARGV != 0) {
    print "wrong number of input: ".@ARGV."\n";
    exit;
}
# $ARGV[0];


#############
# Main starts
#############

foreach my $filename (@filenames) {
    print "-----------------\n$filename\n" if($DEBUG3);
    print "> not exist: $input_dir/$filename.csv\n" if ! -e "$input_dir/$filename.csv";

    open FH, "$input_dir/$filename.csv" or die $!;
    open FH_OUT, "> $output_dir/$filename.txt" or die $!;
    <FH>;   ## ignore the firt line

    ###########
    ## 1. loggingTime
    ## 2. loggingSample
    ## 3. identifierForVendor
    ## 4. locationHeadingTimestamp_since1970
    ## 5. locationHeadingX
    ## 6. locationHeadingY
    ## 7. locationHeadingZ
    ## 8. locationTrueHeading
    ## 9. locationMagneticHeading
    ## 10. locationHeadingAccuracy
    ## 11. accelerometerTimestamp_sinceReboot
    ## 12. accelerometerAccelerationX
    ## 13. accelerometerAccelerationY
    ## 14. accelerometerAccelerationZ
    ## 15. gyroTimestamp_sinceReboot
    ## 16. gyroRotationX
    ## 17. gyroRotationY
    ## 18. gyroRotationZ
    ###########
    while(<FH>) {
        chomp;
        my @data = split(/,/, $_);
        if(@data != 18) {
            print "number of columns = ".@data." != 18\n";
            exit;
        }

        # print FH_OUT join(",", @data[3..8, 11..13, 15..17])."\n";
        print FH_OUT join(",", @data[3..17])."\n";
    }
    close FH;
    close FH_OUT;
}