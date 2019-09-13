import numpy as np
import pandas as pd
import argparse
import sys
import csv
import os
import scipy.stats as ss
# from collections import OrderedDict


def createNumpyArray(dataFile, useAbs=False):
	df = pd.read_csv(dataFile).rename(columns={"Unnamed: 0":"gene"})
	GenesData = df["gene"].tolist()
	TFsData = df.columns.tolist()
	TFsData.remove("gene")
	dataData = df.values[:, 1:].T.tolist()
	if useAbs:
		dataData = np.abs(dataData).tolist()
	return(dataData,GenesData,TFsData)

def createSysDict(namesFile):
	with open(namesFile, mode='r') as infile:
		reader = csv.reader(infile)
		mydict = {rows[1]:rows[0] for rows in reader}
	return mydict

def parseData(locIn,locOut,colToKeep=1,charsToDrop=7):
	DataData = [];
	GenesData = [];
	TFsData = [];

	count = 0
	for filename in os.listdir(locIn):
		print(filename)
		with open(locIn+"/"+filename) as inf:
			reader = csv.reader(inf, delimiter="\t")
			data = list(zip(*reader))
			genes = data[0]
			values = data[colToKeep]
			tfName = filename[:-charsToDrop]
		DataData.append(values[1:])
		GenesData.append(genes[1:])
		TFsData.append(tfName)
		count = count+1

	saveData(locOut, DataData, GenesData, TFsData)


def saveData(loc, DataData, GenesData, TFsData):
	with open(loc+"/"+'Data.csv', 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(DataData)
	with open(loc+"/"+'GeneNames.csv', 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(GenesData)
	with open(loc+"/"+'TFNames.csv', 'w') as writeFile:
		writer = csv.writer(writeFile)
		for val in TFsData:
			writer.writerow([val])

def computeTopNEdges(locIn,locOut,numEdges):
	dataFile = np.loadtxt(locIn)
	dataFileAbs = np.absolute(dataFile)
	rankedFile = rankEdges(dataFileAbs)
	mask = (rankedFile <= numEdges)*1
	thresholdedEdges = dataFile*mask

	with open(locOut+"/"+'Data.csv', 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(thresholdedEdges)

def rankEdges(data):
	m, n = np.shape(data)
	rank = (m*n+1 - ss.rankdata(data)).reshape((m,n))
	return rank


def main(argv):
	# computeTopNEdges('Input Directory Path','Output Directory Path',Number of Edges)
	# parseData('Input Directory Path','Output Directory Path',Column to Keep,Number of Characters to Drop at end of file name):
	computeTopNEdges('/scratch/mblab/pateln/ThresholdAnalysis_Python/ToPost/RawData/NP_ZEV_15_45_90/netprophet2_network.adjmtr','/scratch/mblab/pateln/ThresholdAnalysis_Python/ToPost/Data/NetProphet_Z_15_45_90_150K/',150000)

if __name__ == "__main__":
	main(sys.argv)
