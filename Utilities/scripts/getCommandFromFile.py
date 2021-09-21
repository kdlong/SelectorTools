#!/usr/bin/env python
import ROOT
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("input_file", type=str)
parser.add_argument("--gitdiff", action='store_true')
parser.add_argument("--githash", action='store_true')
args = parser.parse_args()

rtfile = ROOT.TFile(args.input_file)
metaInfo = rtfile.Get("MetaInfo")
if not metaInfo:
    print("ERROR! Can't find metaInfo in file")
    exit(1)
print(metaInfo.Get("command").GetTitle())

if args.githash:
    print(metaInfo.Get("githash").GetTitle())
if args.gitdiff:
    print(metaInfo.Get("gitdiff").GetTitle())
