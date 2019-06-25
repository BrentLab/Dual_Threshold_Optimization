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
  parser.add_argument("-a","--edges_one")
  parser.add_argument("-b","--edges_two")
  parser.add_argument("-o","--output_dir")
  parsed = parser.parse_args(argv[1:])
  return parsed

def computeEdgeUnion(list1,list2):
	network = list1[:]
	for edge in list2:
		if(not edgeInNetwork(edge,network)):
			network.append(edge)

	newNetwork = computePosession(list1,list2,network)
	return newNetwork

def edgeInNetwork(newEdge,network):
	for edge in network:
		# if edge[0] == newEdge[0]:
		if edge[0] == newEdge[0] and edge[1] == newEdge[1]:
			return True
	return False

def computePosession(list1,list2,network):
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

def main(argv):
	global sysDict,parsed
	parsed = parse_args(argv)

	edgesTotal = []
	edges1 = []
	edges2 = []

	edgesOneFile = open(parsed.edges_one, 'r')
	edgesOneFileReader = csv.reader(edgesOneFile)
	edgesOneFileReader.next()
	for row in edgesOneFileReader:
		edges1.append(row)

	edgesTwoFile = open(parsed.edges_two, 'r')
	edgesTwoFileReader = csv.reader(edgesTwoFile)
	edgesTwoFileReader.next()
	for row in edgesTwoFileReader:
		edges2.append(row)

	edgesTotal = computeEdgeUnion(edges1,edges2)
	edgesTotal.insert(0,['TF','Gene','Analysis 1','Analysis 2'])


	with open(parsed.output_dir+"edgesBinary.csv", 'w') as writeFile:
		writer = csv.writer(writeFile)
		writer.writerows(edgesTotal)


if __name__ == "__main__":
	main(sys.argv)