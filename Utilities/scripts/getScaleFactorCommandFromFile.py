#!/usr/bin/env python
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str)
parser.add_argument("--gitdiff", action='store_true')
parser.add_argument("--githash", action='store_true')
args = parser.parse_args()

f = np.load(args.input_file, allow_pickle=True)
metaInfo = f["metaData"]
if metaInfo is None:
    print("ERROR! Can't find metaInfo in file")
    exit(1)
print(metaInfo[0])

if args.githash:
    print(metaInfo[1])
if args.gitdiff:
    print(metaInfo[2])

