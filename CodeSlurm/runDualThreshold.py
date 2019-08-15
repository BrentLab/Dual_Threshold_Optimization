import numpy as np
import pandas as pd
import argparse
import random
import sys
sys.path.append("/scratch/mblab/yiming.kang/dual_threshold_optimization/CodeSlurm/")
from loadData import *
from statistics import *
from thresholdSearch import *
from scipy.stats import hypergeom,spearmanr
from datetime import datetime

np.set_printoptions(threshold=np.nan)


def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-d","--de_file")
	parser.add_argument("-b","--bin_file")
	parser.add_argument("-r","--random",default = False)
	parser.add_argument("-g","--geneNames_file",default = "/scratch/mblab/pateln/ThresholdAnalysis_Python/ExtraFiles/YeastCommonAndSystematicGeneNames.csv")
	parser.add_argument("-s","--incl_spearman",default = False)
	parser.add_argument("-n","--TF_num")
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

def runDualThresholds(DEData,BinData,TFIntersection,TFNum,GenesUniverse_Given):
	global sysDict,parsed

	if GenesUniverse_Given!=None:
		GenesUniverse = GenesUniverse_Given

	TF = TFIntersection[int(TFNum)]
	print(TF)
	sys.stdout.flush()
	BinDataData,BinGenesData,BinTFsData = BinData
	DEDataData,DEGenesData,DETFsData = DEData
	
	BinTFIndex = BinTFsData.index(TF)
	DETFIndex = BinTFsData.index(TF)
	DEData,BinData = alignToUniverse(DEData,BinData,GenesUniverse,BinTFIndex,DETFIndex,str2Bool(parsed.Bin_decreasing))
	BinDataData,BinGenesData,BinTFsData = BinData
	DEDataData,DEGenesData,DETFsData = DEData

	if parsed.random:
		DEValues,DEGenes = sortDataRandom(DEDataData[DETFsData.index(TF)], DEGenesData)
	else:
		DEValues,DEGenes = sortData(DEDataData[DETFsData.index(TF)], DEGenesData, str2Bool(parsed.DE_decreasing))
	
	BinValues,BinGenes = sortData(BinDataData[BinTFsData.index(TF)], BinGenesData, str2Bool(parsed.Bin_decreasing))

	if GenesUniverse_Given==None:
		GenesUniverse = computeUniverse(DEGenes,BinGenes)

	DEValues,DEGenes = getSubset(DEGenes,DEValues,GenesUniverse)
	BinValues,BinGenes = getSubset(BinGenes,BinValues,GenesUniverse)

	optimizedResults = optimizedThresholds(TF,BinGenes,BinValues,DEGenes,DEValues,GenesUniverse)

	if parsed.random:
		fileName = TF + "_" + parsed.iter_num + ".csv"
	else:
		fileName = TF + ".csv"

	with open(fileName,'w') as resultFile:
		wr = csv.writer(resultFile)
		print(optimizedResults)
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

	print('DE list = ',DErankList)
	print('bin list = ',binRankList)
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

	DEData = createNumpyArray(parsed.de_file)
	BinData = createNumpyArray(parsed.bin_file)
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

	# DEData,BinData = alignToUniverse(DEData,BinData,GenesUniverse)

	runDualThresholds(DEData,BinData,TFIntersection,parsed.TF_num,GenesUniverse)

if __name__ == "__main__":
	main(sys.argv)