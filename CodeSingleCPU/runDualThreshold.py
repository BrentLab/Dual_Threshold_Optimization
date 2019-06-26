import numpy as np
import pandas as pd
import argparse
import random
import sys
sys.path.append("/scratch/mblab/pateln/ThresholdAnalysis_Python/Code")
from loadData import *
from statistics import *
from thresholdSearch import *
from scipy.stats import hypergeom,spearmanr
from datetime import datetime
from computeAcceptableTFs import *

np.set_printoptions(threshold=np.nan)

def str2Bool(boolAsString):
	return bool(distutils.util.strtobool(boolAsString))

def parse_args(argv):
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-d","--de_dir")
  parser.add_argument("-b","--bin_dir")
  parser.add_argument("-r","--random",default = False)
  parser.add_argument("-g","--geneNames_file",default = "/scratch/mblab/pateln/ThresholdAnalysis_Python/ExtraFiles/YeastCommonAndSystematicGeneNames.csv")
  parser.add_argument("-s","--incl_spearman",default = False)
  parser.add_argument("-o","--opt_crit",default = "pval")
  parser.add_argument("-w","--rank_width",default = 1.01)
  parser.add_argument("-i","--iter_num",default = "")
  parser.add_argument("-l","--tf_list",default = "")
  parser.add_argument("-u","--genes_universe",default = "")
  parser.add_argument("-a","--organism",default = "yeast")
  parser.add_argument("-j","--DE_decreasing",default = True)
  parser.add_argument("-k","--Bin_decreasing",default = True)
  parsed = parser.parse_args(argv[1:])
  return parsed

# from https://stackoverflow.com/questions/3160699/python-progress-bar
# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 10 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rProgress: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()

def generateRanks(length,power,Data,threshold,decreasing):
	rankList = []
	currRank = 1
	if decreasing==True:
		newRankedData = [-1*num for num in Data]
		while(currRank < length):
			if Data[int(currRank)] > threshold and Data[int(currRank)]!=0:
				rankList.append(int(currRank))
			currRank = round(currRank*power + 1)
	else:
		newRankedData = Data
		while(currRank < length):
			if Data[int(currRank)] < threshold:
				rankList.append(int(currRank))
			currRank = round(currRank*power + 1)
	
	finalRankList = removeRedundantRanks(rankList,newRankedData,decreasing)

	return finalRankList

def removeRedundantRanks(rankList,rankedData,decreasing):
	scipyRanks = rankdata(rankedData, method='max')
	newRankList = []
	for rank in rankList:
		newRankList.append(scipyRanks[rank-1])

	ranks = sorted(list(set(newRankList)))
	return ranks


def computeUniverse(DEGenes,BinGenes,givenUniverse = None):
	# Universe computed as intersection of Binding/DE genes
	# return list(set(DEGenes) | set(BinGenes)) # union
	if givenUniverse == None:
		GenesUniverse = list(set(DEGenes) & set(BinGenes)) # intersection
	else:
		file = open(givenUniverse, 'r')
		GenesUniverse = file.readlines()
		for i in range(len(GenesUniverse)):
			GenesUniverse[i] = str.strip(GenesUniverse[i])
	return GenesUniverse

def alignToUniverse(DEData,BinData,GenesUniverse,binIndex,DEIndex,decreasing=True):
	BinDataData,BinGenesData,BinTFsData = BinData
	DEDataData,DEGenesData,DETFsData = DEData

	tempGenes = BinGenesData[binIndex]
	tempData = BinDataData[binIndex]
	for gene in GenesUniverse:
		if gene not in tempGenes:
			tempGenes.append(gene)
			if decreasing == True:
				tempData.append(0)
			else:
				tempData.append(1)
	BinGenesData[binIndex] = tempGenes
	BinDataData[binIndex] = tempData
	newBinData = BinDataData,BinGenesData,BinTFsData
	newDEData = DEDataData,DEGenesData,DETFsData
	return newDEData,newBinData


def getSubset(Genes,Data,otherSet):
	inds,temp = overlap(Genes,otherSet)
	newGenes = Genes[inds]
	newData = [Data[i] for i in inds]
	return(newData,newGenes)

def overlap(a, b):
	# return the indices in a that overlap with b, also returns 
	# the corresponding index in b only works if both a and b are unique! 
	# This is not very efficient but it works
	bool_a = np.in1d(a,b)
	ind_a = np.arange(len(a))
	ind_a = ind_a[bool_a]

	ind_b = np.array([np.argwhere(b == a[x]) for x in ind_a]).flatten()
	return ind_a,ind_b

def sortData(Data,Genes,decreasing=True):
	if decreasing==True:
		DataAbs = [abs(number) for number in Data]
		Genes = np.array(Genes)[np.argsort(DataAbs)[::-1]]
		Data = np.array(Data)[np.argsort(DataAbs)[::-1]]
	else:
		Genes = np.array(Genes)[np.argsort(Data)]
		Data = sorted(Data)
	return(Data,Genes)

def sortDataRandom(Data,Genes):
	c = list(zip(Data,Genes))
	random.shuffle(c)
	data = zip(*c)
	Data = np.array(data[0])
	Genes = np.array(data[1])
	return (Data,Genes)

def sysToCommon(genesList,sysDict):
	intersectionCommon = list(map(sysDict.get, genesList))
	if None in intersectionCommon:
		intersectionCommon = []
		for gene in genesList:
			if sysDict.get(gene) != None:
				intersectionCommon.append(sysDict[gene])
			else:
				intersectionCommon.append(gene)
		return intersectionCommon
	else:	
		return intersectionCommon

def runDualThresholds(DEData,BinData,TFIntersection,TFNum,GenesUniverse_Given,rand,iterNum=0):
	global sysDict,parsed

	if GenesUniverse_Given!=None:
		GenesUniverse = GenesUniverse_Given

	TF = TFIntersection[int(TFNum)]
	# print(TF)
	sys.stdout.flush()
	BinDataData,BinGenesData,BinTFsData = BinData
	DEDataData,DEGenesData,DETFsData = DEData
	
	BinTFIndex = BinTFsData.index(TF)
	DETFIndex = BinTFsData.index(TF)
	DEData,BinData = alignToUniverse(DEData,BinData,GenesUniverse,BinTFIndex,DETFIndex,str2Bool(parsed.Bin_decreasing))
	BinDataData,BinGenesData,BinTFsData = BinData
	DEDataData,DEGenesData,DETFsData = DEData

	if rand:
		DEValues,DEGenes = sortDataRandom(DEDataData[DETFsData.index(TF)], DEGenesData[DETFsData.index(TF)])
	else:
		DEValues,DEGenes = sortData(DEDataData[DETFsData.index(TF)], DEGenesData[DETFsData.index(TF)],str2Bool(parsed.DE_decreasing))
	
	BinValues,BinGenes = sortData(BinDataData[BinTFsData.index(TF)], BinGenesData[BinTFsData.index(TF)],str2Bool(parsed.Bin_decreasing))

	if GenesUniverse_Given==None:
		GenesUniverse = computeUniverse(DEGenes,BinGenes)

	DEValues,DEGenes = getSubset(DEGenes,DEValues,GenesUniverse)
	BinValues,BinGenes = getSubset(BinGenes,BinValues,GenesUniverse)

	optimizedResults = optimizedThresholds(TF,BinGenes,BinValues,DEGenes,DEValues,GenesUniverse)

	if rand==True:
		fileName = "../Results/output_rand.csv"
		with open(fileName,'a') as resultFile:
			wr = csv.writer(resultFile)
			wr.writerow([TF,optimizedResults[6]])
	else:
		fileName = "../Results/output.csv"
		with open(fileName,'a') as resultFile:
			wr = csv.writer(resultFile)
			wr.writerow(optimizedResults)

def optimizedThresholds(TF,BinGenes,BinValues,DEGenes,DEValues,GenesUniverse):
	global sysDict
	if parsed.organism == "yeast":
		TFCommon = sysDict[TF]
	else:
		TFCommon = TF

	if str2Bool(parsed.DE_decreasing) == True:
		DErankList = generateRanks(len(DEValues),float(parsed.rank_width),DEValues,0,str2Bool(parsed.DE_decreasing))
	else:
		if parsed.organism == "yeast":
			DErankList = generateRanks(len(DEValues),float(parsed.rank_width),DEValues,1,str2Bool(parsed.DE_decreasing))
		else:
			DErankList = generateRanks(len(DEValues),float(parsed.rank_width),DEValues,0.1,str2Bool(parsed.DE_decreasing))
	
	if str2Bool(parsed.Bin_decreasing) == True:
		binRankList = generateRanks(len(BinValues),float(parsed.rank_width),BinValues,0,str2Bool(parsed.Bin_decreasing))
	else:
		binRankList = generateRanks(len(BinValues),float(parsed.rank_width),BinValues,0.1,str2Bool(parsed.Bin_decreasing))

	bestFDR = 1
	bestPVal = 1
	bestFE = 0
	bestRR = 0
	bestRRsk = 0
	bestSC = 0
	bestSPVal = 1
	bestJS = 0
	bestIntersection = 0
	boundSubGenes = []
	DESubGenes = []

	if len(DErankList) > 0:
		bestDEThresh = DErankList[0]
	else:
		bestDEThresh = 0
	if len(binRankList) > 0:
		bestBinThresh = binRankList[0]
	else:
		bestBinThresh = 0

	# print('DE list = ',DErankList)
	# print('bin list = ',binRankList)
	sys.stdout.flush()

	for i in range(len(binRankList)):
		binThresh = binRankList[i]
		boundSubGenes = BinGenes[:binThresh]
		boundSubData = BinValues[:binThresh]
		for j in range(len(DErankList)):
			DEThresh = DErankList[j]
			DESubGenes = DEGenes[:DEThresh]
			DESubData = DEValues[:DEThresh]
			GenesIntersection = list(set(DESubGenes) & set(boundSubGenes))
			if len(GenesIntersection) > 0:
				FDRBound = computeFDRLowerBound(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
				hyperPval = computeHyperPVal(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
				responseRate = computeResponseRate(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
				foldEnrichment = computeFoldEnrichment(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
				relativeRisk = computeRelativeRisk(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes)
				# rho,spearmanPval = computeSpearman(boundSubPVals,boundSubGenes,DESubPVals,DESubGenes,GenesIntersection)
				rho,spearmanPval = [0,1]
				jaccardSim = computeJaccardSimilarity(DESubGenes,boundSubGenes)


				optCrit = parsed.opt_crit
				if optCrit=="fdr":
					if FDRBound < bestFDR and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
					if FDRBound == bestFDR and hyperPval < bestPVal and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
					# if FDRBound == bestFDR and len(GenesIntersection) > bestIntersection and len(GenesIntersection) > 1:
					# 	bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
					# 	bestBinThresh = binThresh
					# 	bestDEThresh = DEThresh
				if optCrit=="pval":
					if hyperPval < bestPVal and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
				elif optCrit=="fe":
					if foldEnrichment > bestFE and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
				elif optCrit=="rr":
					if responseRate > bestRR and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
				elif optCrit=="rrsk":
					if responseRate > bestRRsk and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
				elif optCrit=="sc":
					if abs(rho) > abs(bestSC) and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh
				elif optCrit=="js":
					if jaccardSim > bestJS and len(GenesIntersection) > 1:
						bestFDR,bestIntersection,bestPVal,bestFE,bestRR,bestRRsk,bestSC,bestSPVal,bestJS = [FDRBound,len(GenesIntersection),hyperPval,foldEnrichment,responseRate,relativeRisk,rho,spearmanPval,jaccardSim]
						bestBinThresh = binThresh
						bestDEThresh = DEThresh

	boundSubGenes = BinGenes[:bestBinThresh]
	boundSubData = BinValues[:bestBinThresh]
	DESubGenes = DEGenes[:bestDEThresh]
	DESubData = DEValues[:bestDEThresh]
	GenesIntersection = list(set(DESubGenes) & set(boundSubGenes))
	if parsed.organism == "yeast":
		intersectionCommon = sorted(sysToCommon(GenesIntersection,sysDict))
	else:
		intersectionCommon = sorted(GenesIntersection)

	# print(DEValues[bestDEThresh])
	# print(BinValues[bestBinThresh])
	# d = [TF,TFCommon,len(boundSubGenes),len(DESubGenes),len(GenesIntersection),bestPVal,bestRR,bestRRsk,bestFE,bestSC,bestSPVal,bestJS,intersectionCommon]
	d = [TF,TFCommon,len(boundSubGenes),len(DESubGenes),len(GenesIntersection),bestFDR,bestPVal,bestRR,bestRRsk,bestFE,bestJS,intersectionCommon]
	return d

def main(argv):
	global sysDict,parsed
	parsed = parse_args(argv)

	DEData = createNumpyArray(parsed.de_dir)
	BinData = createNumpyArray(parsed.bin_dir)
	sysDict = createSysDict(parsed.geneNames_file)
	TFIntersection = sorted(list(set(DEData[2]) & set(BinData[2])))


	if parsed.tf_list == "":
		TFIntersection = sorted(list(set(DEData[2]) & set(BinData[2])))
	else:
		file = open(parsed.tf_list, 'r')
		TFIntersection = file.readlines()
		for i in range(len(TFIntersection)):
			TFIntersection[i] = str.strip(TFIntersection[i])

	if parsed.genes_universe != "":
		file = open(parsed.genes_universe, 'r')
		GenesUniverse = file.readlines()
		for i in range(len(GenesUniverse)):
			GenesUniverse[i] = str.strip(GenesUniverse[i])
	else:
		GenesUniverse = None

	# Perform standard analysis
	print('Performing DTO analysis')
	fileName = "../Results/output.csv"
	with open(fileName,'a') as resultFile:
		wr = csv.writer(resultFile)
		wr.writerow(['TF','TFCommon','Bound Size','DE Size','Intersection','FDR Lower Bound','HypergeometricPVal','Response Rate','Relative Risk','Fold Enrichment','Jaccard Similarity','Genes'])
	for TFNum in range(len(TFIntersection)):
		update_progress(round((1+TFNum*1.0)/len(TFIntersection),4))
		runDualThresholds(DEData,BinData,TFIntersection,TFNum,GenesUniverse,False)

	if str2Bool(parsed.random)==True:
		# Perform randomized analysis
		print('Performing randomized DTO analysis')
		fileName = "../Results/output_rand.csv"
		numRands = 1000
		for TFNum in range(len(TFIntersection)):
			for iterNum in range(numRands):
				update_progress(round((1+TFNum*numRands+iterNum*1.0)/len(TFIntersection*numRands),4))
				runDualThresholds(DEData,BinData,TFIntersection,TFNum,GenesUniverse,True,iterNum)

		# Compute acceptable TFs
		computeTFs('../Results/output.csv','../Results/output_rand.csv','../Results/')


	# if str2Bool(parsed.random)==True:
	# 	fileName = "../Results/output_rand.csv"
	# else:
	# 	fileName = "../Results/output.csv"
	# 	with open(fileName,'a') as resultFile:
	# 		wr = csv.writer(resultFile)
	# 		wr.writerow("'TF','TFCommon','Bound Size','DE Size','Intersection','FDR Lower Bound','HypergeometricPVal','Response Rate','Relative Risk','Fold Enrichment','Jaccard Similarity','Genes'\n")

	# if parsed.tf_list == "":
	# 	TFIntersection = sorted(list(set(DEData[2]) & set(BinData[2])))
	# else:
	# 	file = open(parsed.tf_list, 'r')
	# 	TFIntersection = file.readlines()
	# 	for i in range(len(TFIntersection)):
	# 		TFIntersection[i] = str.strip(TFIntersection[i])

	# if parsed.genes_universe != "":
	# 	file = open(parsed.genes_universe, 'r')
	# 	GenesUniverse = file.readlines()
	# 	for i in range(len(GenesUniverse)):
	# 		GenesUniverse[i] = str.strip(GenesUniverse[i])
	# else:
	# 	GenesUniverse = None

	# if str2Bool(parsed.random)==True:
	# 	for TFNum in range(len(TFIntersection)):
	# 		for iterNum in range(10):
	# 			runDualThresholds(DEData,BinData,TFIntersection,TFNum,GenesUniverse,iterNum)

	# else:
	# 	for TFNum in range(len(TFIntersection)):
	# 		runDualThresholds(DEData,BinData,TFIntersection,TFNum,GenesUniverse)

if __name__ == "__main__":
	main(sys.argv)