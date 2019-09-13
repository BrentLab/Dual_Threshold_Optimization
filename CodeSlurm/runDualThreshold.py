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
np.set_printoptions(threshold=np.nan)


def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-d","--de_file")
	parser.add_argument("-b","--bin_file")
	parser.add_argument("-g","--geneNames_file")
	parser.add_argument("-n","--TF_num", type=int)
	parser.add_argument("-o","--opt_crit", default="pval")
	parser.add_argument("-w","--rank_width", type=float, default=1.01)
	parser.add_argument("-l","--tf_list", default=None)
	parser.add_argument("-u","--genes_universe", default=None)
	parser.add_argument("-r","--random_iter", default=None, type=int)
	parser.add_argument("-a","--organism", default="yeast")
	parser.add_argument("-j","--DE_decreasing", default=True)
	parser.add_argument("-k","--Bin_decreasing", default=True)
	parsed = parser.parse_args(argv[1:])
	return parsed


def str2Bool(boolAsString):
	return bool(distutils.util.strtobool(boolAsString))


def getTargetedTF(DEFile, BinFile, TFNum):
	with open(DEFile, "r") as f:
		DEHeader = f.readline().strip().split(",")[1:]
	with open(BinFile, "r") as f:
		BinHeader = f.readline().strip().split(",")[1:]
	TFIntersection = set(DEHeader) & set(BinHeader)
	if parsed.tf_list is not None:
		TFIntersection &= set(list(np.loadtxt(parsed.tf_list, dtype=str)))
	targetedTF = sorted(TFIntersection)[TFNum]
	DEIdx = DEHeader.index(targetedTF)
	BinIdx = BinHeader.index(targetedTF)
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
	##TODO: double check the validity of setting zeros or ones
	if decreasing:
		values2 = np.ones(nUnivGenes) * max(values)
	else:
		values2 = np.ones(nUnivGenes) * min(values)
	for i,gene in enumerate(GenesUniverse):
		if gene in genes:
			j = genes.index(gene)
			values2[i] = values[j]
	return (list(values2), GenesUniverse, TF)


def runDualThresholds(DEData, BinData):
	TF = DEData[-1]
	if parsed.random_iter:
		DEData = randomizeData(DEData, randSeed=parsed.random_iter)
	else:
		DEData = sortData(DEData, str2Bool(parsed.DE_decreasing))
	BinData = sortData(BinData, str2Bool(parsed.Bin_decreasing))

	GenesUniverse = sorted(DEData[1])
	optimizedResults = optimizeThresholds(DEData, BinData, GenesUniverse)
	print(optimizedResults[:-1])

	if parsed.random_iter >= 0:
		fileName = TF + "_" + str(parsed.random_iter) + ".csv"
	else:
		fileName = TF + ".csv"
	with open(fileName,'w') as resultFile:
		wr = csv.writer(resultFile)
		print(optimizedResults)
		wr.writerow(optimizedResults)


def sortData(Data, decreasing=True):
	values, genes, TF = Data
	sorted_idx = np.argsort(values)[::-1] if decreasing else np.argsort(values)
	return (np.array(values)[sorted_idx], np.array(genes)[sorted_idx], TF)


def randomizeData(Data, randSeed=0):
	values, genes, TF = Data
	combined = list(zip(values, genes))
	randSeed += sum([ord(x) for x in TF])
	print(randSeed)
	random.seed(randSeed)
	random.shuffle(combined)
	values, genes = zip(*combined)
	return (np.array(values), np.array(genes), TF)


# def generateRanks(length,scaler,Data,threshold,decreasing):
def generateRanks(values, scaler, threshold, decreasing):
	valuesLength = len(values)
	rankList = []
	currRank = 1
	if decreasing==True:
		newRankedData = [-1*num for num in values]
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
	if parsed.organism == "yeast":
		TFCommon = sysDict[TF]
	else:
		TFCommon = TF

	if str2Bool(parsed.DE_decreasing):
		DE_threshold = 0
	else: 
		if parsed.organism == "yeast":
			DE_threshold = 1
		else:
			DE_threshold = 0.1
	if str2Bool(parsed.Bin_decreasing):
		Bin_threshold = 0
	else:
		Bin_threshold = 0.1
	DErankList = generateRanks(DEValues, parsed.rank_width, DE_threshold, 
								str2Bool(parsed.DE_decreasing))
	binRankList = generateRanks(BinValues, parsed.rank_width, Bin_threshold, 
								str2Bool(parsed.Bin_decreasing))
	print('DE list = %s' % DErankList)
	print('bin list = %s' % binRankList)
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

	##TODO: update argmin/max search
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
	if parsed.organism == "yeast":
		intersectionCommon = sorted(sysToCommon(GenesIntersection,sysDict))
	else:
		intersectionCommon = sorted(GenesIntersection)

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
	sysDict = createSysDict(parsed.geneNames_file)

	## Map datasets to common gene universe
	if parsed.genes_universe is None:
		GenesUniverse = sorted(set(DEData[1]) & set(BinData[1]))
	else:
		GenesUniverse = sorted(list(np.loadtxt(parsed.genes_universe, dtype=str)))
	DEData = alignToUniverse(DEData, GenesUniverse, parsed.DE_decreasing)
	BinData = alignToUniverse(BinData, GenesUniverse, parsed.Bin_decreasing)

	## Run dual threshold optimization
	runDualThresholds(DEData, BinData)
	print("Elapsed: %.5f" % (time.time() - tStart))


if __name__ == "__main__":
	main(sys.argv)