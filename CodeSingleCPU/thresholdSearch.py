from loadData import *
from statistics import *
import numpy as np
import pandas as pd
import argparse
import sys
import random
import os
from scipy.stats import hypergeom,spearmanr,rankdata
from datetime import datetime
import distutils.core

np.set_printoptions(threshold=np.nan)

def str2Bool(boolAsString):
	return bool(distutils.util.strtobool(boolAsString))

def parse_args(argv):
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-d","--de_dir")
  parser.add_argument("-b","--bin_dir")
  parser.add_argument("-r","--random")
  parser.add_argument("-g","--geneNames_file")
  parser.add_argument("-f","--sbatch_loc",default = ".")
  parser.add_argument("-w","--rank_width",default = "1.01")
  parser.add_argument("-o","--opt_crit",default = "pval")
  parser.add_argument("-u","--genes_universe",default = "")
  parser.add_argument("-a","--organism",default = "yeast")
  parser.add_argument("-n","--rand_type",default = "split")
  parser.add_argument("-j","--DE_decreasing")
  parser.add_argument("-k","--Bin_decreasing")

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

def prepDualThresholds(TFIntersection):
	codeDir = os.getcwd()
	# for i in range(5):
	if(str2Bool(parsed.random) == False):
		createSbatchFile(len(TFIntersection),codeDir)
		os.chdir(parsed.sbatch_loc)
		for filename in os.listdir('.'):
			os.system("sbatch "+filename)
	else:
		if(parsed.rand_type == "global"):
			numIterations = (1000/len(TFIntersection)) + 2
			for iterNum in range(numIterations):
				createSbatchFile(len(TFIntersection),codeDir,iterNum,numIterations)
			os.chdir(parsed.sbatch_loc)
			for filename in os.listdir('.'):
				os.system("sbatch "+filename)
		else:
			for TFNum in range(len(TFIntersection)):
				TF = TFIntersection[TFNum]
				os.makedirs(parsed.sbatch_loc+TF)
				createSbatchFile(1000,codeDir,1,1,TFNum,TF)
				os.chdir(parsed.sbatch_loc+'/'+TF)
				for filename in os.listdir('.'):
					os.system("sbatch "+filename)
				os.chdir(codeDir)

	os.chdir(codeDir)

def createSbatchFile(numTFs,codeDir,iterNum="",numIters=1,TFNum=1,TF=""):
	global parsed

	if parsed.genes_universe == "":
		universe = '""'
	else:
		universe = parsed.genes_universe

	jobName = "dualThresholds_%A_%a"
	if(str2Bool(parsed.random) == False):
		f = open(parsed.sbatch_loc+"/runAnalysis.sbatch", 'w')
	else:
		if(parsed.rand_type == "global"):
			f = open(parsed.sbatch_loc+"/runAnalysis_"+str(iterNum)+".sbatch", 'w')
		else:
			f = open(parsed.sbatch_loc+'/'+TF+"/runAnalysis_"+str(iterNum)+".sbatch", 'w')
	f.write("#!/bin/sh\n")
	f.write("#SBATCH -D ./\n")
	f.write("#SBATCH --mem=2G\n")
	if(str2Bool(parsed.random) == False):
		f.write("#SBATCH -J "+jobName+"\n")
		f.write("#SBATCH -o "+jobName+".out\n")
		f.write("#SBATCH -e "+jobName+".err\n")
		f.write("#SBATCH --array=0-"+str(numTFs)+"%50\n")
		f.write("ID=${SLURM_ARRAY_TASK_ID}\n")
	else:
		f.write("#SBATCH -J "+jobName+"_"+str(iterNum)+"\n")
		f.write("#SBATCH -o "+jobName+"_"+str(iterNum)+".out\n")
		f.write("#SBATCH -e "+jobName+"_"+str(iterNum)+".err\n")
		if(parsed.rand_type == "global"):
			f.write("#SBATCH --array=0-"+str(numTFs)+"%"+str(100/numIters)+"\n")
			f.write("ID=${SLURM_ARRAY_TASK_ID}\n")
		else:
			f.write("#SBATCH --array=0-9\n")
		
	f.write("module load numpy\n")
	f.write("module load pandas\n\n")

	if(str2Bool(parsed.random) == True and parsed.rand_type != "global"):
		f.write("START=$(( SLURM_ARRAY_TASK_ID * 100 ))\n")
		f.write("STOP=$(( START + 99 ))\n")
		f.write("[ \"$STOP\" -eq 999 ] && STOP=1000\n\n")
		f.write("for ID in $( seq $START $STOP ); do\n")

	if(str2Bool(parsed.random) == False):
		f.write("python " + codeDir + "/runDualThreshold.py --de_dir " + parsed.de_dir + " --bin_dir " + parsed.bin_dir + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing) + " --TF_num " + "${ID}" + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --organism " + parsed.organism + "\n")
	else:
		if(parsed.rand_type == "global"):
			f.write("python " + codeDir + "/runDualThreshold.py --de_dir " + parsed.de_dir + " --bin_dir " + parsed.bin_dir + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing)+ " --TF_num " + "${ID}" + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --organism " + parsed.organism + " --iter_num " + str(iterNum) + " --random True" + "\n")
		else:
			f.write("\tpython " + codeDir + "/runDualThreshold.py --de_dir " + parsed.de_dir + " --bin_dir " + parsed.bin_dir + " --DE_decreasing " + str(parsed.DE_decreasing) + " --Bin_decreasing " + str(parsed.Bin_decreasing) + " --TF_num " + str(TFNum) + " --rank_width " + parsed.rank_width + " --opt_crit " + parsed.opt_crit + " --genes_universe " + universe + " --organism " + parsed.organism + " --iter_num " + "${ID}" + " --random True" + "\n")
			f.write("done\n")
	f.close()

	

def main(argv):
	global sysDict,parsed
	parsed = parse_args(argv)

	DEData = createNumpyArray(parsed.de_dir)
	BinData = createNumpyArray(parsed.bin_dir)
	sysDict = createSysDict(parsed.geneNames_file)

	TFIntersection = sorted(list(set(DEData[2]) & set(BinData[2])))

	prepDualThresholds(TFIntersection)


if __name__ == "__main__":
	main(sys.argv)