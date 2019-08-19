import numpy as np
import pandas as pd
import argparse
import sys
import csv
import os
import scipy.stats as ss
# from collections import OrderedDict


def createNumpyArray(dataFile):
	df = pd.read_csv(dataFile).rename(columns={"Unnamed: 0":"gene"})
	GenesData = df["gene"].tolist()
	TFsData = df.columns.tolist()
	TFsData.remove("gene")
	dataData = df.values[:, 1:].T.tolist()
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

	# otherLoc = '/Users/patelni/WolframWorkspaces/Base/ThresholdAnalysis/Python Code/Data/Yeast/Kemmeren'
	saveData(locOut, DataData, GenesData, TFsData)
	
def parseENCODEshRNA(locIn,locOut):
	parseKemmeren(locIn,locOut)

def parseENCODEChIP(locIn,locOut):
	pValsData = [];
	LFCsData = [];
	GenesData = [];
	TFsData = [];

	count = 0
	for filename in os.listdir(locIn):
		print(filename)
		with open(locIn+"/"+filename) as inf:
			reader = csv.reader(inf, delimiter="\t")
			data = list(zip(*reader))
			genes = data[0]
			pvals = data[3]
			lfcs = data[3]
			tfName = filename[:-16]

		uniqueGenes = list(set(genes[1:]))
		genesNP = np.array(genes)
		pvalsNP = np.array(pvals)
		newPValsList = []

		for gene in uniqueGenes:
			geneRowNums = np.where(genesNP == gene)[0]
			genePvals = pvalsNP[geneRowNums]
			newPval = 0
			# print(genePvals)
			for pval in genePvals:
				if(float(pval) > 1): #checks to see if the -log10(pval) is greater than 1 (or that the pvalue is less than 0.1)
					newPval+=float(pval)

			newPval = 10**(-newPval)
			newPValsList.append(newPval)

		pValsData.append(newPValsList)
		LFCsData.append(newPValsList)
		GenesData.append(uniqueGenes)
		TFsData.append(tfName)
		count = count+1

	# otherLoc = '/Users/patelni/WolframWorkspaces/Base/ThresholdAnalysis/Python Code/Data/Yeast/Kemmeren'
	saveData(locOut, pValsData, LFCsData, GenesData, TFsData)


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

def splitFiles(loc):
	DataFile = open(loc+'/Data.csv', 'r')
	GenesFile = open(loc+'/GeneNames.csv', 'r')
	TFsFile = open(loc+'/TFNames.csv', 'r')

	DataReader = csv.reader(DataFile)
	GenesReader = csv.reader(GenesFile)
	TFsReader = csv.reader(TFsFile)

	# ordered_fieldnames = OrderedDict([('field1',None),('field2',None)])

	for row in TFsReader:
		TF = [ str(elem) for elem in row ][0]
		with open(loc+"/Split Files/"+TF + '.txt', 'w') as writeFile:
			writer = csv.writer(writeFile, delimiter = '\t')
			writer.writerow(['#genes','#PVal']) 
			genes = next(GenesReader)
			values = next(DataReader)
			for i in range(len(genes)):
				writer.writerow([genes[i],values[i]]) 

def checkInsufficientTFs(loc):
	DataData,GenesData,TFsData = createNumpyArray(loc)
	insuffTFs = np.zeros(len(TFsData))

	for i in range(len(TFsData)):
		TF = TFsData[i]
		Genes = GenesData[i]
		Data = DataData[i]

		if TF in list(Genes) and Data[list(Genes).index(TF)] > np.log2(0.5):
			insuffTFs[i] = 1

	return insuffTFs

def removeInsufficientTFs(loc):
	insuffTFs = checkInsufficientTFs(loc)

	DataData = []
	GenesData = []
	TFsData = []

	DataFile = open(loc+'/Data.csv', 'r')
	GenesFile = open(loc+'/GeneNames.csv', 'r')
	TFsFile = open(loc+'/TFNames.csv', 'r')

	DataReader = csv.reader(DataFile)
	GenesReader = csv.reader(GenesFile)
	TFsReader = csv.reader(TFsFile)

	for i in range(len(insuffTFs)):
		DataRow = next(DataReader)
		GenesRow = next(GenesReader)
		TFsRow = next(TFsReader)
		if insuffTFs[i] == 0:
			DataData.append( [ float(elem) for elem in DataRow ] )
			GenesData.append( [ str(elem) for elem in GenesRow ] )
			TFsData.append( [ str(elem) for elem in TFsRow ][0] )

	otherLoc = '/Users/patelni/WolframWorkspaces/Base/ThresholdAnalysis/Python Code/Data/Yeast/Kemmeren'
	saveData(otherLoc, DataData, GenesData, TFsData)

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
