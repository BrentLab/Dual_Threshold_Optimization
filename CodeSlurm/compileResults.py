import os
import sys
import numpy as np
import pandas as pd
import glob
import argparse
import shutil
from loadData import *
from statistics import *
from thresholdSearch import *
import ast

csv.field_size_limit(sys.maxsize)

def parse_args(argv):
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-i", "--input_dir")
  parser.add_argument("-o", "--output_file")
  parser.add_argument("-r", "--rand_type", default="global")
  parsed = parser.parse_args(argv[1:])
  return parsed

def compileTargetsFixed(loc):
	resultsFile = open(loc, 'r')
	resultsReader = csv.reader(resultsFile)
	targetsData = []
	next(resultsReader)
	for row in resultsReader:
		TF = row[0]
		targets = row[-1:][0]
		# if row[5] != "NA" and float(row[5]) < 0.2 and float(row[6]) < 0.01:
		if row[5] != "NA" and float(row[5]) < 0.2:
			targetsData.append([TF,targets])

	return targetsData

def compileFromCSV(loc,threshold):
	os.chdir(loc)
	targetsData = []
	print("TFs with >= 1000 targets: " + loc)
	for file in glob.glob("*.csv"):
		if file != "results.csv":
			csvfile = open(file)
			resultsReader = csv.reader(csvfile)
			for row in resultsReader:
				TF = row[0]
				TFCommon = row[1]
				pVal = row[6]
				targets = row[-1:][0]
				numTargets = int(row[4])
				if row[5] != "NA" and float(row[5]) < 0.2 and float(row[6]) < threshold:
					if numTargets >= 1000:
						print(TFCommon)
					targetsData.append([TF,targets,pVal])
	return targetsData

# def compileFromCSV2(loc):
# 	os.chdir(loc)
# 	targetsData = []
# 	print("Set 2 TFs with >= 1000 targets")
# 	for file in glob.glob("*.csv"):
# 		if file != "results.csv":
# 			csvfile = open(file)
# 			resultsReader = csv.reader(csvfile)
# 			for row in resultsReader:
# 				TF = row[0]
# 				TFCommon = row[1]
# 				targets = row[-1:][0]
# 				if row[5] != "NA" and float(row[5]) < 0.2 and float(row[6]) < 0.00000610:
# 					if len(targets) >= 1000:
# 						print(TFCommon)
# 					targetsData.append([TF,targets])
# 	return targetsData

def computeEdgeUnion(list1,list2,name1,name2):
	network = list1.copy()
	for edge in list2:
		if(not edgeInNetwork(edge,network)):
			network.append(edge)

	newNetwork = computePosession(list1,list2,network,name1,name2)
	return newNetwork

def computePosession(list1,list2,network,name1,name2):
	newNetwork = []
	for edge in network:
		newEdge = edge
		if (edge in list1) and (edge in list2):
			# newEdge.append(name1 + " and " + name2)
			newEdge.append(1)
			newEdge.append(1)
		elif (edge in list1) and not (edge in list2):
			# newEdge.append(name1)
			newEdge.append(1)
			newEdge.append(0)
		elif (edge in list2) and not (edge in list1):
			# newEdge.append(name2)
			newEdge.append(0)
			newEdge.append(1)
		newNetwork.append(newEdge)

	return newNetwork

def computeEdgeIntersection(list1,list2,name1,name2):
	network = []
	for edge in list1:
		if(edgeInNetwork(edge,list2)):
			network.append(edge)

	return network

def edgeInNetwork(newEdge,network):
	for edge in network:
		if edge[0] == newEdge[0] and edge[1] == newEdge[1]:
			return True
	return False

def computeNumLargerPvals(TFIntersection,data1,data2):
	TFList1 = []
	TFList2 = []
	for row1 in data1:
		TF = row1[0]
		if row1[0] in TFIntersection:
			for row2 in data2:
				if TF==row2[0]:
					if float(row1[2])<float(row2[2]):
						TFList1.append(TF)
					else:
						TFList2.append(TF)

	return TFList1,TFList2

def compileTargetIntersection(loc1,loc2,name1="Harbison v Kemmeren",name2="Harbison v ZEV",thresh1 = 1,thresh2 = 1):
	# loc1Data = compileTargetsFixed(loc1)
	loc1Data = compileFromCSV(loc1,thresh1)
	loc2Data = compileFromCSV(loc2,thresh2)

	loc1TFs = [row[0] for row in loc1Data]
	loc1Genes = [row[1] for row in loc1Data]
	loc1Targets = []
	loc1Edges = []
	for i in range(len(loc1Data)):
		TF = loc1Data[i][0]
		sublist = loc1Genes[i]
		temp = ast.literal_eval(sublist)
		temp = [n.strip() for n in temp]
		for item in temp:
			loc1Targets.append(item)
			loc1Edges.append([TF,item])
	loc1UniqueTargets = set(loc1Targets)

	loc2TFs = [row[0] for row in loc2Data]
	loc2Genes = [row[1] for row in loc2Data]
	loc2Targets = []
	loc2Edges = []
	for i in range(len(loc2Data)):
		TF = loc2Data[i][0]
		sublist = loc2Genes[i]
		temp = ast.literal_eval(sublist)
		temp = [n.strip() for n in temp]
		for item in temp:
			loc2Targets.append(item)
			loc2Edges.append([TF,item])
	loc2UniqueTargets = set(loc2Targets)

	TFUnion = list(set(loc1TFs) | set(loc2TFs))
	TFIntersection = list(set(loc1TFs) & set(loc2TFs))
	combinedUniqueTargetsUnion = list(set(loc1UniqueTargets) | set(loc2UniqueTargets))
	combinedUniqueTargetsIntersection = list(set(loc1UniqueTargets) & set(loc2UniqueTargets))
	combinedUniqueEdgesUnion = computeEdgeUnion(loc1Edges,loc2Edges,name1,name2)
	combinedUniqueEdgesIntersection = computeEdgeIntersection(loc1Edges,loc2Edges,name1,name2)

	set1TFsBetterPVals,set2TFsBetterPVals = computeNumLargerPvals(TFIntersection,loc1Data,loc2Data)

	with open("/Users/patelni/Desktop/Harbison v Kemmeren and ZEV15.csv", "w") as f:
		writer = csv.writer(f)
		writer.writerow(["Systematic", "Common", name1, name2])
		writer.writerows(combinedUniqueEdgesUnion)

	print("Significant Results (FDR < 0.2 & HPval < 0.01)")
	print("----------------------------------------------")
	print("Dataset 1:")
	print("Number of TFs: ",len(loc1Data))
	print("Number of total targets: ",len(loc1Targets))
	print("Number of unique targets: ",len(loc1UniqueTargets))
	print("Number of Edges: ",len(loc1Edges))
	print()
	print("Dataset 2:")
	print("Number of TFs: ",len(loc2Data))
	print("Number of total targets: ",len(loc2Targets))
	print("Number of unique targets: ",len(loc2UniqueTargets))
	print("Number of Edges: ",len(loc2Edges))
	print()
	print("Combined:")
	print("Number of TFs (Union): ",len(TFUnion))
	print("Number of TFs (Intersection): ",len(TFIntersection))
	print("Number of Unique Targets (Union): ",len(combinedUniqueTargetsUnion))
	print("Number of Unique Targets (Intersection): ",len(combinedUniqueTargetsIntersection))
	print("Number of Unique Edges (Union)",len(combinedUniqueEdgesUnion))
	print("Number of Unique Edges (Intersection)",len(combinedUniqueEdgesIntersection))
	print()
	print("Number of Acceptable TFs with smaller pVals in set 1: ",len(set1TFsBetterPVals))
	print("Number of Acceptable TFs with smaller pVals in set 2: ",len(set2TFsBetterPVals))


def main(argv):
	global parsed
	parsed = parse_args(argv)

	if parsed.rand_type == "global":
		header = ['TF','TFCommon','Bound Size','DE Size','Intersection',
				'FDR Lower Bound','HypergeometricPVal','Response Rate',
				'Relative Risk','Fold Enrichment','Jaccard Similarity','Genes']
		os.system("echo %s >> %s" % (",".join(header), parsed.output_file))
		for csvFile in glob.glob("%s/*.csv" % parsed.input_dir):
			os.system("cat %s >> %s" %(csvFile, parsed.output_file))

	else:
		TFDirs = glob.glob("%s/*/" % parsed.input_dir)
		for i, TFDir in enumerate(TFDirs):
			sys.stdout.write("%d " % (i+1))
			for csvFile in glob.glob("%s/*.csv" % TFDir):
				os.system("cut -d',' -f1,7 %s >> %s" % (csvFile, parsed.output_file))


if __name__ == "__main__":
	main(sys.argv)