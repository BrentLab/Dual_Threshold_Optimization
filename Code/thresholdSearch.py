import numpy as np
import pandas as pd
import argparse
import sys
import random
import os
import glob
from scipy.stats import hypergeom,spearmanr,rankdata
from datetime import datetime
import distutils.core
from loadData import *


# np.set_printoptions(threshold=np.nan)
np.set_printoptions(threshold=sys.maxsize)

def str2Bool(boolAsString):
	return bool(distutils.util.strtobool(boolAsString))

def parse_args(argv):
	parser = argparse.ArgumentParser(description="DTO for identifying convergent evidence of TF binding and TF-perturbation response.")
	parser.add_argument("-d", "--de_file", required=True,
						help="CSV file of transcriptional response levels of TF-perturbation. The response levels are based on differential expression or DE analysis on TF-perturbed expression profile vs. pre-perturbation profile. Those levels can be represented as absolute of log fold-changes (LFCs) or p-values.")
	parser.add_argument("-b","--bin_file", required=True,
						help="CSV file of TF binding strengths. The binding strengths can be represented as occupancy levels of the peaks (e.g. heights of ChIP peaks) or statistical significances (e.g. p-values).")
	parser.add_argument("--DE_decreasing", required=True,
						help="Use 'True' if the response levels should be raneked in descending order. Use 'False' if ranked in ascending order. For example: If using LFCs, set 'True' as higher absolute LFC represents stronger response. If using p-values, set 'False' as lower p-value represents stronger response.")
	parser.add_argument("--Bin_decreasing", required=True,
						help="Use 'True' if the binding strengths should be ranked in descending order. Use 'False' if ranked in ascending order. For example: If using occupancy levels, set 'True' as higher occupancy level represents stronger TF binding. If using p-values, set 'False' as lower p-value represents stronger TF binding.")
	parser.add_argument("-f","--sbatch_loc", required=True,
						help="Directory for output data.")
	parser.add_argument("-u","--genes_universe", default="",
						help="Gene universe list file, which contains all possible genes to be considered as the universe for calculating overlap statistics.")
	parser.add_argument("-g","--geneNames_file", default="",
						help="CSV file of gene name conversion, if the systematic gene names in the data files are preferred to be converted into common gene names.")
	parser.add_argument("-r","--random", default="False",
						help="Use 'False' (default) if running DTO on authentic (non-randomized) input data. Use 'True' if running DTO on randomized input data. The DTO results on randomized data is used to create an empirical null distribution.")
	parser.add_argument("-n","--rand_type", default="split",
						help="Use 'split' (default) to create a null distribution for each TF. Use 'global' to create a single null distribution for all TFs.")
	parser.add_argument("-w","--rank_width", default="1.01",
						help="Rate of stride for defining the search space.")
	parser.add_argument("-o","--opt_crit", default="pval",
						help="Optimization criterion (or objective): hypergeometric p-value (default).")
	parser.add_argument("--DE_pval_lower_bound", default=1,
						help="Lower bound of response p-values. Any genes with p-values greater than the lower bound will not be used for optimization.")
	parser.add_argument("--Bin_pval_lower_bound", default=0.1,
						help="Lower bound of binding p-values. Any genes with p-values greater than the lower bound will not be used for optimization.")
	parser.add_argument("--run_local", action='store_true', default=False,
						help="Flag for running DTO in serial fashion on local machine.")
	parser.add_argument("--find_tf_specificity", default="False", help = "find where the hypergeometric p-val of TFs of interest lie in the distribution of false TF pairings.")
	parsed = parser.parse_args(argv[1:])
	return parsed


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

	tempGenes = BinGenesData
	tempData = BinDataData[binIndex]
	for gene in GenesUniverse:
		if gene not in tempGenes:
			tempGenes.append(gene)
			if decreasing == True:
				tempData.append(0)
			else:
				tempData.append(1)
	BinGenesData = tempGenes
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

def prepDualThresholds(TFIntersection, run_local=False):
	execution = "bash" if run_local else "sbatch"
	codeDir = os.getcwd()

	if not os.path.exists(parsed.sbatch_loc):
		os.makedirs(parsed.sbatch_loc)

	if(str2Bool(parsed.random) == False):
		createSbatchFile(len(TFIntersection),codeDir)
		os.chdir(parsed.sbatch_loc)
		if not os.path.exists('log'):
			os.makedirs('log')
		for filename in glob.glob('*.sbatch'):
			os.system(execution+" "+filename)
	else:
		if(parsed.rand_type == "global"):
			numIterations = (1000/len(TFIntersection)) + 2
			for iterNum in range(numIterations):
				createSbatchFile(len(TFIntersection),codeDir,iterNum,numIterations)
			os.chdir(parsed.sbatch_loc)
			if not os.path.exists('log'):
				os.makedirs('log')
			for filename in glob.glob('*.sbatch'):
				os.system(execution+" "+filename)
		else:
			for TFNum in range(len(TFIntersection)):
				TF = TFIntersection[TFNum]
				if not os.path.exists(parsed.sbatch_loc+TF):
					os.makedirs(parsed.sbatch_loc+TF)
				createSbatchFile(1000,codeDir,1,1,TFNum,TF)
				os.chdir(parsed.sbatch_loc+'/'+TF)
				if not os.path.exists('log'):
					os.makedirs('log')
				for filename in glob.glob('*.sbatch'):
					os.system(execution+" "+filename)
				os.chdir(codeDir)

	os.chdir(codeDir)


def findTfSpecificity(TFIntersection, run_local = False):
	execution = "bash" if run_local else "sbatch"
	codeDir = os.getcwd()

	if not os.path.exists(parsed.sbatch_loc):
		os.makedirs(parsed.sbatch_loc)

	directory_DTORun = parsed.sbatch_loc + "/find_tf_specificity/"
	if not os.path.exists(directory_DTORun):
		os.makedirs(directory_DTORun)

	for TFNum in range(len(TFIntersection)):
		TF = TFIntersection[TFNum]
		if not os.path.exists(directory_DTORun +TF):
			os.makedirs(directory_DTORun +TF)
		createSbatchFile(len(TFIntersection),codeDir, "", 1, TFNum, TF, True)		
		os.chdir(directory_DTORun+'/'+TF)
		if not os.path.exists('log'):
			os.makedirs('log')
		for filename in glob.glob('*.sbatch'):
			os.system(execution+" "+filename)
		os.chdir(codeDir)

def createSbatchFile(numTFs,codeDir,iterNum="",numIters=1,TFNum=1,TF="",TF_Specificity = False):
	global parsed

	if parsed.genes_universe == "":
		universe = '""'
	else:
		universe = parsed.genes_universe
	if parsed.geneNames_file == "":
		geneNames = '""'
	else:
		geneNames = parsed.geneNames_file

	jobName = "runDTO_%A_%a"
	if(str2Bool(parsed.random) == False and TF_Specificity == False):
		f = open(parsed.sbatch_loc+"/runAnalysis.sbatch", 'w')
	elif(str2Bool(parsed.random) == False and TF_Specificity == True):
		f = open(parsed.sbatch_loc+ "/find_tf_specificity/" + TF + "/" +  "runAnalysis.sbatch", 'w')
	else:
		if(parsed.rand_type == "global"):
			f = open(parsed.sbatch_loc+"/runAnalysis_"+str(iterNum)+".sbatch", 'w')
		else:
			f = open(parsed.sbatch_loc+"/"+TF+"/runAnalysis_"+str(iterNum)+".sbatch", 'w')
	f.write("#!/bin/bash\n")
	f.write("#SBATCH -D "+codeDir+"\n")
	f.write("#SBATCH --mem=2G\n")
	if(str2Bool(parsed.random) == False and TF_Specificity == False):
		f.write("#SBATCH -J "+jobName+"\n")
		f.write("#SBATCH -o "+parsed.sbatch_loc +"/log/"+jobName+".out\n")
		f.write("#SBATCH -e "+parsed.sbatch_loc +"/log/"+jobName+".err\n")
		f.write("#SBATCH --array=0-"+str(numTFs-1)+"%50\n")
		f.write("ID=${SLURM_ARRAY_TASK_ID}\n")

	elif(str2Bool(parsed.random) == False and TF_Specificity == True):
		f.write("#SBATCH -J "+jobName+"_"+str(iterNum)+"\n")
		f.write("#SBATCH -o "+ parsed.sbatch_loc+"/find_tf_specificity/"+ TF + "/log/"+jobName+"_"+str(iterNum)+".out\n")
		f.write("#SBATCH -e "+parsed.sbatch_loc+"/find_tf_specificity/"+ TF + "/log/"+jobName+"_"+str(iterNum)+".err\n")
		
		f.write("#SBATCH --array=0-"+str(numTFs*2-3)+"%200\n")
		f.write("ID=${SLURM_ARRAY_TASK_ID}\n")

	else:
		f.write("#SBATCH -J "+jobName+"_"+str(iterNum)+"\n")
		f.write("#SBATCH -o "+parsed.sbatch_loc+"/"+TF+"/log/"+jobName+"_"+str(iterNum)+".out\n")
		f.write("#SBATCH -e "+parsed.sbatch_loc+"/"+TF+"/log/"+jobName+"_"+str(iterNum)+".err\n")
		if(parsed.rand_type == "global"):
			f.write("#SBATCH --array=0-"+str(numTFs-1)+"%"+str(100/numIters)+"\n")
			f.write("ID=${SLURM_ARRAY_TASK_ID}\n")
		else:
			f.write("#SBATCH --array=0-9\n")

	if(str2Bool(parsed.random) == True and parsed.rand_type != "global"):
		f.write("START=$(( SLURM_ARRAY_TASK_ID * 100 ))\n")
		f.write("STOP=$(( START + 99 ))\n")
		f.write("[ \"$STOP\" -eq 999 ] && STOP=1000\n\n")
		f.write("for ID in $( seq $START $STOP ); do\n")
	if(str2Bool(parsed.random) == False and TF_Specificity == False):
		f.write(
			"python runDualThreshold.py --de_file " + parsed.de_file + " --bin_file " + parsed.bin_file + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing) + " --TF_num " + "${ID}" + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --geneNames_file " + geneNames + " --DE_pval_lower_bound " + str(parsed.DE_pval_lower_bound)  + " --Bin_pval_lower_bound " + str(parsed.Bin_pval_lower_bound) + " --output_dir " + parsed.sbatch_loc + "\n")
	
	elif(str2Bool(parsed.random) == False and TF_Specificity == True):
		f.write(
			"python runDualThreshold.py --de_file " + parsed.de_file + " --bin_file " + parsed.bin_file + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing) + " --TF_num " + str(TFNum) + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --geneNames_file " + geneNames + " --DE_pval_lower_bound " + str(parsed.DE_pval_lower_bound)  + " --Bin_pval_lower_bound " + str(parsed.Bin_pval_lower_bound) + " --output_dir " + parsed.sbatch_loc + "/find_tf_specificity/" + TF + " --find_tf_specificity " + str(TF_Specificity) + " --tuple_index_false_pairing " + "${ID}" + "\n")

	else:
		if(parsed.rand_type == "global"):
			f.write(
				"python " + codeDir + "/runDualThreshold.py --de_file " + parsed.de_file + " --bin_file " + parsed.bin_file + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing)+ " --TF_num " + "${ID}" + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --geneNames_file " + geneNames + " --DE_pval_lower_bound " + str(parsed.DE_pval_lower_bound)  + " --Bin_pval_lower_bound " + str(parsed.Bin_pval_lower_bound) + " --random_iter " + str(iterNum) + " --output_dir " + parsed.sbatch_loc + "\n")
			f.write(
				"\tpython " + codeDir + "/runDualThreshold.py --de_file " + parsed.de_file + " --bin_file " + parsed.bin_file + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing) + " --TF_num " + str(TFNum) + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --geneNames_file " + geneNames + " --DE_pval_lower_bound " + str(parsed.DE_pval_lower_bound)  + " --Bin_pval_lower_bound " + str(parsed.Bin_pval_lower_bound) + " --random_iter " + "${ID}" + " --output_dir " + parsed.sbatch_loc + "/" + TF + "\n")
			f.write("done\n")
	f.close()
	

def main(argv):
	global sysDict,parsed
	parsed = parse_args(argv)
	DEData = createNumpyArray(parsed.de_file)
	BinData = createNumpyArray(parsed.bin_file)
	TFIntersection = sorted(list(set(DEData[2]) & set(BinData[2])))
	prepDualThresholds(TFIntersection, parsed.run_local)

	if(str2Bool(parsed.random) == False and str2Bool(parsed.find_tf_specificity) == True):
		findTfSpecificity(TFIntersection, parsed.run_local)


if __name__ == "__main__":
	main(sys.argv)
