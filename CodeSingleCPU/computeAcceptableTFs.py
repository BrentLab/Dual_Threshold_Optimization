import numpy as np
import pandas as pd
import argparse
import random
import sys
from loadData import *
from statistics import *
from thresholdSearch import *
from compileResults import *
from scipy.stats import hypergeom,spearmanr
from datetime import datetime

def parse_args(argv):
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-d","--data_file")
  parser.add_argument("-r","--rand_file")
  parser.add_argument("-o","--output_dir")
  parser.add_argument("-t","--threshold",default = 0.01)
  parsed = parser.parse_args(argv[1:])
  return parsed

def computeCutoff(pVals,threshold):
	if len(pVals) < 10:
		return 0.00001
	position = int(round(len(pVals)*threshold))
	sortedPVals = sorted(pVals)
	cutoff = sortedPVals[position]
	return cutoff

def computeEdges(acceptableTFs):
	edges = []
	for TFData in acceptableTFs:
		TF = TFData[0]
		genesList = TFData[-1]
		genes = genesList.split('\'')[1::2];
		for gene in genes:
			edges.append([TF,gene])

	return edges


def computeTFs(data_file, rand_file, output_dir, threshold=0.01):
	acceptableTFs = []
	cutoffs = []
	randData = []
	RandFile = open(rand_file, 'r')
	randFileReader = csv.reader(RandFile)
	for row in randFileReader:
		randData.append(row[0:2])

	df = pd.DataFrame(randData,columns=['TF','PVal'])

	# Reading in the output file
	outputFile = open(data_file, 'r')
	outputFileReader = csv.reader(outputFile)
	acceptableTFs.append(outputFileReader.next()) #skips header row
	# print(outputFileReader.next())
	for row in outputFileReader:
		TF = row[0]
		HyperPVal = float(row[6])
		FDR = float(row[5])
		locs = df.loc[df['TF'] == TF]
		pValsList = locs['PVal'].values.tolist()
		pVals = [float(num) for num in pValsList]
		cutoff = computeCutoff(pVals,threshold)
		cutoffs.append([TF,cutoff])
		if(HyperPVal < cutoff and FDR <= 0.2):
			acceptableTFs.append(row)

	edges = computeEdges(acceptableTFs)
	uniqueGenes = list(set([row[1] for row in edges]))
	uniqueGenes = [[row] for row in uniqueGenes]

	with open(output_dir+"acceptableTFs.csv", 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(acceptableTFs)

	with open(output_dir+"TFcutoffs.csv", 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(cutoffs)

	with open(output_dir+"edges.csv", 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(edges)

	with open(output_dir+"targets.csv", 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(uniqueGenes)

	with open(output_dir+'summary.txt','w')  as writeFile:
		writeFile.write("Number of acceptableTFs: "+str(len(acceptableTFs)-1))
		writeFile.write("\nNumber of Edges: "+str(len(edges)))
		writeFile.write("\nNumber of Unique Targets: "+str(len(uniqueGenes)))


def main(argv):
	global sysDict,parsed
	parsed = parse_args(argv)
	computeTFs(parsed.data_file,parsed.rand_file,parsed.output_dir,parsed.threshold)


if __name__ == "__main__":
	main(sys.argv)