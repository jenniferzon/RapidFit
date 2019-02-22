#source /cvmfs/sft.cern.ch/lcg/views/LCG_93/x86_64-slc6-gcc62-opt/setup.sh
from ROOT import TFile, TTree, TString
from ROOT import gSystem, gROOT
from ROOT import *
import sys
import re, array
from array import array

gROOT.ProcessLine("struct weights_t {vector<double> *pvect;}")
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json

# Make it work for Python 2+3 and with Unicode
import io
try:
    to_unicode = unicode
except NameError:
    to_unicode = str

# Read JSON file
with open('/home/s1668588/Bs2JpsiPhi-Run2/ANALYSIS/analysis/fitinputs/fitinputs_st015_trigCat.json') as data_file:
    data_loaded = json.load(data_file)
    data_file.close() # to enable writing to later on

rootfile = TFile.Open("testAngAcc.root", "read")
roottree = rootfile.Get("tree")
tree_orig = weights_t()
roottree.SetBranchAddress("weights", AddressOf(tree_orig,'pvect'))
roottree.GetEntry(0)

## For what datasets create tables
inputs = [ "2015","2016"]

def createAngAcc(year, parameters, trigger):
    print "Creating AngAccWeights_{}_data_015_{}.root".format(year, trigger)
    newfile = TFile.Open("AngAccWeights_{}_data_015_{}.root".format(year, trigger), "recreate");
    newtree = roottree.CloneTree(0);
    tree = tree_orig
    newtree.SetBranchAddress("weights", AddressOf(tree,'pvect'))
    tree.pvect = tree_orig.pvect
    newtree.GetEntry(0)
    for i in range(10):
      tree.pvect[i] = parameters[i]
    newtree.Fill()
    newtree.AutoSave()
    newfile.Write()
    newfile.Save()
    newfile.Close()

for inp in inputs:
    ntuple = data_loaded[inp]["Ntuple"][0]["Name"].replace("_","\_")
    biased  = [float(en["Value"]) for en in data_loaded[inp]["AngularParameterBiased"]]
    unbiased  = [float(en["Value"]) for en in data_loaded[inp]["AngularParameterUnbiased"]]
    createAngAcc(inp, biased, "biased")
    createAngAcc(inp, unbiased, "unbiased")
