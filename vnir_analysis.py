# -*- coding: utf-8 -*-
import sys
import csv
import argparse

# http://sebastianraschka.com/Articles/2014_pca_step_by_step.html

parser = argparse.ArgumentParser(description='Reads in a file to an array')
parser.add_argument("-d", "--directory", dest="myDir", default="",
                    help="Open specified directory")
parser.add_argument("-e", "--extension", dest="myExt", default="txt",
                    help="The extension to match")
args = parser.parse_args()

print(args.myDir + " " + args.myExt)
