import numpy as np
import pandas as pd
import argparse
import random
import itertools
import sys
import csv
import time
import multiprocessing as mp
import distutils.core
from scipy.stats import rankdata
from loadData import createSysDict
from thresholdSearch import sysToCommon
from statistics import *

np.set_printoptions(threshold=sys.maxsize)


EPSILON = 10**(-5)

def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-d", "--de_file")
	parser.add_argument("-b", "--bin_file")
	parser.add_argument("-g", "--geneNames_file", default=None)
	parser.add_argument("-n", "--TF_num", type=int)
	parser.add_argument("-o", "--opt_crit", default="pval")
	parser.add_argument("-w", "--rank_width", type=float, default=1.01)
	parser.add_argument("-l", "--tf_list", default=None)
	parser.add_argument("-u", "--genes_universe", default=None)
	parser.add_argument("-r", "--random_iter", default=None, type=int)
	parser.add_argument("--DE_decreasing", default=True)
	parser.add_argument("--Bin_decreasing", default=True)
	parser.add_argument("--DE_pval_lower_bound", type=float, default=1)
	parser.add_argument("--Bin_pval_lower_bound", type=float, default=0.1)
	parser.add_argument("--find_tf_specificity", default="False")
	parser.add_argument("--tuple_index_false_pairing", type = int, default=0)
	parser.add_argument("--output_dir")
	parsed = parser.parse_args(argv[1:])
	return parsed


def str2Bool(boolAsString):
	return bool(distutils.util.strtobool(boolAsString))


def getTargetedTF(DEFile, BinFile, TFNum):
	global deTF, bindingTF
	with open(DEFile, "r") as f:
		DEHeader = f.readline().strip().split(",")[1:]
	with open(BinFile, "r") as f:
		BinHeader = f.readline().strip().split(",")[1:]
	TFIntersection = set(DEHeader) & set(BinHeader)
	if parsed.tf_list is not None:
		TFIntersection &= set(list(np.loadtxt(parsed.tf_list, dtype=str)))

	if(str2Bool(parsed.find_tf_specificity) == False):
		targetedTF = sorted(TFIntersection)[TFNum]
		deTF = targetedTF
		bindingTF = targetedTF
		DEIdx = DEHeader.index(targetedTF)
		BinIdx = BinHeader.index(targetedTF)

	else:
		targetedTF = sorted(TFIntersection)[TFNum]
		tupleList = list(range(0,len(TFIntersection)))
		pairingList1 = [(TFNum,val) for val in tupleList]
		pairingList2 = [(val,TFNum) for val in tupleList]
		set1 = set(pairingList1)
		set2 = set(pairingList2)
		false_pairing_set = set1.union(set2)- set1.intersection(set2)
		false_pairing_list = list(false_pairing_set)
		print(len(false_pairing_list))
		print(false_pairing_list)
		curr_tuple = false_pairing_list[parsed.tuple_index_false_pairing]
		deTF = sorted(TFIntersection)[curr_tuple[0]]
		DEIdx = DEHeader.index(deTF)
		bindingTF = sorted(TFIntersection)[curr_tuple[1]]
		BinIdx = BinHeader.index(bindingTF)
		print(str(deTF) + ":" +  str(bindingTF))

	return (targetedTF, DEIdx, BinIdx)


def getTargetedTFData(dataFile, targetedIdx, targetedTF, useAbs=False):
	data = np.loadtxt(dataFile, usecols=[0, (targetedIdx+1)], 
					dtype=str, delimiter=",")
	genes = data[1:, 0]
	values = data[1:, 1].astype(float)
	if useAbs:
		values = np.abs(values)
	return (values.tolist(), genes.tolist(), targetedTF)


def alignToUniverse(data, GenesUniverse, decreasing=True):
	values, genes, TF = data
	nUnivGenes = len(GenesUniverse)
	if decreasing:
		values2 = np.ones(nUnivGenes) * (min(values) - EPSILON)
	else:
		values2 = np.ones(nUnivGenes) * (max(values) + EPSILON)
	for i,gene in enumerate(GenesUniverse):
		if gene in genes:
			j = genes.index(gene)
			values2[i] = values[j]
	return (list(values2), GenesUniverse, TF)
# def alignToUniverse(data, GenesUniverse):
# 	values, genes, TF = data
# 	idx = [i for i,x in enumerate(genes) if x in GenesUniverse]
# 	return (list(np.array(values)[idx]), list(np.array(genes)[idx]), TF)


def runDualThresholds(DEData, BinData, GenesUniverse):
	TF = DEData[2]
	if parsed.random_iter:
		DEData = randomizeData(DEData, randSeed=parsed.random_iter)
	else:
		DEData = sortData(DEData, str2Bool(parsed.DE_decreasing))
	BinData = sortData(BinData, str2Bool(parsed.Bin_decreasing))

	optimizedResults = optimizeThresholds(DEData, BinData, GenesUniverse)

	#TODO: Add more items to the optimized results when tf_specificity is true
	if(str2Bool(parsed.find_tf_specificity)):
		optimizeResults = optimizedResults[:-1]
		print(deTF + ": check")
		print(bindingTF + ": check")
		print("Optimized results+ = %s" % optimizedResults)
	else:
		print("Optimized results = %s" % optimizedResults[:-1])

	if parsed.random_iter:
		fileName = parsed.output_dir + "/" + TF + "_" + str(parsed.random_iter) + ".csv"
	elif(str2Bool(parsed.find_tf_specificity) == True):
		fileName = parsed.output_dir + "/" + TF + "_" + str(parsed.tuple_index_false_pairing) + ".csv"
	else:
		fileName = parsed.output_dir + "/" + TF + ".csv"
	with open(fileName,'w') as resultFile:
		wr = csv.writer(resultFile)
		wr.writerow(optimizedResults)
	# elif str2Bool(parsed.find_tf_specificity:

def sortData(Data, decreasing=True):
	values, genes, TF = Data
	sorted_idx = np.argsort(values)[::-1] if decreasing else np.argsort(values)
	return (np.array(values)[sorted_idx], np.array(genes)[sorted_idx], TF)


def randomizeData(Data, randSeed=0):
	values, genes, TF = Data
	randSeed += sum([ord(x) for x in TF])
	print("Random seed = %d" % randSeed)
	random.seed(randSeed)
	combined = list(zip(values, genes))
	random.shuffle(combined)
	values, genes = zip(*combined)
	return (np.array(values), np.array(genes), TF)


def generateRanks(values, scaler, threshold, decreasing):
	valuesLength = len(values)
	rankList = []
	currRank = 1
	if decreasing==True:
		newRankedData = -values
		while(currRank < valuesLength):
			if values[int(currRank)] > threshold and values[int(currRank)]!=0:
				rankList.append(int(currRank))
			currRank = round(currRank*scaler + 1)
	else:
		newRankedData = values
		while(currRank < valuesLength):
			if values[int(currRank)] < threshold:
				rankList.append(int(currRank))
			currRank = round(currRank*scaler + 1)
	finalRankList = removeRedundantRanks(rankList, newRankedData)
	return finalRankList


def removeRedundantRanks(rankList, rankedData):
	scipyRanks = rankdata(rankedData, method='max')
	newRankList = []
	for rank in rankList:
		newRankList.append(scipyRanks[rank-1])
	ranks = sorted(set(newRankList))
	return ranks


def optimizeThresholds(DEData, BinData, GenesUniverse):
	global DEGenes, BinGenes
	DEValues, DEGenes, TF = DEData
	BinValues, BinGenes, _ = BinData
	
	if parsed.geneNames_file is None or parsed.geneNames_file == "":
		TFCommon = TF
	else:
		TFCommon = sysDict[TF]

	## TODO: expose response and binding threshold 
	if str2Bool(parsed.DE_decreasing):
		DE_threshold = 0
	else: 
		DE_threshold = parsed.DE_pval_lower_bound
	if str2Bool(parsed.Bin_decreasing):
		Bin_threshold = 0
	else:
		Bin_threshold = parsed.Bin_pval_lower_bound
	DErankList = generateRanks(DEValues, parsed.rank_width, DE_threshold, 
								str2Bool(parsed.DE_decreasing))
	binRankList = generateRanks(BinValues, parsed.rank_width, Bin_threshold, 
								str2Bool(parsed.Bin_decreasing))
	# print('DE list = %s' % DErankList)
	# print('bin list = %s' % binRankList)
	sys.stdout.flush()

	bestIntersection = 0
	boundSubGenes = []
	DESubGenes = []

	bestDEThresh = DErankList[0] if len(DErankList) > 0 else 0
	bestBinThresh = binRankList[0] if len(binRankList) > 0 else 0

	comboRankList = list(itertools.product(binRankList, DErankList))
	pool = mp.Pool(8)
	stats = pool.map(calculateStat, comboRankList)
	pool.close()
	pool.join()

	# updateBest = False
	# optCrit = parsed.opt_crit
	# if len(GenesIntersection) > 1:
	# 	if optCrit=="fdr":
	# 		if FDRBound < bestFDR:
	# 			updateBest = True
	# 		if FDRBound == bestFDR and hyperPval < bestPVal:
	# 			updateBest = True
	# 	if optCrit=="pval":
	# 		if hyperPval < bestPVal:
	# 			updateBest = True
	# 	elif optCrit=="fe":
	# 		if foldEnrichment > bestFE:
	# 			updateBest = True
	# 	elif optCrit=="rr":
	# 		if responseRate > bestRR:
	# 			updateBest = True
	# 	elif optCrit=="rrsk":
	# 		if responseRate > bestRRsk:
	# 			updateBest = True
	# 	elif optCrit=="sc":
	# 		if abs(rho) > abs(bestSC):
	# 			updateBest = True
	# 	elif optCrit=="js":
	# 		if jaccardSim > bestJS:
	# 			updateBest = True

	idx = np.nanargmin(np.array(stats)[:,0])
	bestStat, bestBinThresh, bestDEThresh = stats[idx]

	boundSubGenes = BinGenes[:bestBinThresh]
	DESubGenes = DEGenes[:bestDEThresh]
	GenesIntersection = list(set(DESubGenes) & set(boundSubGenes))
		
	if parsed.geneNames_file is None or parsed.geneNames_file == "":
		intersectionCommon = sorted(GenesIntersection)
	else:
		intersectionCommon = sorted(sysToCommon(GenesIntersection,sysDict))
	
	FDRBound = computeFDRLowerBound(GenesIntersection,GenesUniverse,
									DESubGenes,boundSubGenes)
	hyperPval = computeHyperPVal(GenesIntersection,GenesUniverse,
									DESubGenes,boundSubGenes)
	responseRate = computeResponseRate(GenesIntersection,GenesUniverse,
									DESubGenes,boundSubGenes)
	foldEnrichment = computeFoldEnrichment(GenesIntersection,GenesUniverse,
									DESubGenes,boundSubGenes)
	relativeRisk = computeRelativeRisk(GenesIntersection,GenesUniverse,
									DESubGenes,boundSubGenes)
	jaccardSim = computeJaccardSimilarity(DESubGenes,boundSubGenes)
	out = [TF, TFCommon, len(boundSubGenes), len(DESubGenes), len(GenesIntersection),
			FDRBound, hyperPval, responseRate, relativeRisk, 
			foldEnrichment, jaccardSim, intersectionCommon]
	return out


def calculateStat(threshTuple):
	binThresh, DEThresh = threshTuple
	boundSubGenes = BinGenes[:binThresh]
	DESubGenes = DEGenes[:DEThresh]
	GenesIntersection = list(set(DESubGenes) & set(boundSubGenes))
	if len(GenesIntersection) > 1:
		if parsed.opt_crit.lower() == "fdr":
			stat = computeFDRLowerBound(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
		elif parsed.opt_crit.lower() == "pval":
			stat = computeHyperPVal(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
		elif parsed.opt_crit.lower() == "fe":
			stat = computeFoldEnrichment(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
		elif parsed.opt_crit.lower() == "rr":
			stat = computeRelativeRisk(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
		elif parsed.opt_crit.lower() == "js":
			stat = computeJaccardSimilarity(DESubGenes,boundSubGenes)
	else:
		stat = np.nan
	return (stat, binThresh, DEThresh)


#use slurm id instead of original id for TFs 
def main(argv):
	tStart = time.time()
	global parsed, sysDict, GenesUniverse
	parsed = parse_args(argv)

	## Load data for the targeted TF
	targetedTF, targetedDEIdx, targetedBinIdx = \
		getTargetedTF(parsed.de_file, parsed.bin_file, parsed.TF_num)
	sys.stdout.write("%s\n" % targetedTF)
	sys.stdout.flush()
	DEData = getTargetedTFData(parsed.de_file, targetedDEIdx, 
								targetedTF, useAbs=True)
	BinData = getTargetedTFData(parsed.bin_file, targetedBinIdx, 
								targetedTF, useAbs=True)
	if parsed.geneNames_file is None or parsed.geneNames_file == "":
		sysDict = {}
	else:
		sysDict = createSysDict(parsed.geneNames_file)

	## Map datasets to common gene universe
	if parsed.genes_universe is None or parsed.genes_universe == "":
		GenesUniverse = sorted(set(DEData[1]) & set(BinData[1]))
	else:
		GenesUniverse = sorted(list(np.loadtxt(parsed.genes_universe, dtype=str)))
	# DEData = alignToUniverse(DEData, GenesUniverse)
	# BinData = alignToUniverse(BinData, GenesUniverse)
	DEData = alignToUniverse(DEData, GenesUniverse, parsed.DE_decreasing)
	BinData = alignToUniverse(BinData, GenesUniverse, parsed.Bin_decreasing)

	## Run dual threshold optimization
	runDualThresholds(DEData, BinData, GenesUniverse)
	print("Elapsed: %.5f" % (time.time() - tStart))


if __name__ == "__main__":
	main(sys.argv)
