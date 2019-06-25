from thresholdSearch import *
import numpy as np
import pandas as pd
import sys
from scipy.stats import hypergeom,spearmanr

def computeHyperPVal(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes):
	pval = hypergeom.sf(len(GenesIntersection)-1, len(GenesUniverse), len(DESubGenes), len(boundSubGenes))
	return pval

def computeFoldEnrichment(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes):
	if len(boundSubGenes) > 0 and len(DESubGenes) > 0:
		FE = len(GenesIntersection)/((len(DESubGenes)*len(boundSubGenes)*1.0)/len(GenesUniverse))
	else:
		FE = "NA"
	return FE

def computeResponseRate(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes):
	if len(boundSubGenes) > 0:
		return float(len(GenesIntersection))/len(boundSubGenes)
	else:
		return "NA"

def computeRelativeRisk(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes):
	if len(GenesIntersection) > 0 and len(boundSubGenes) > 0 and len(DESubGenes) > len(GenesIntersection) and len(GenesUniverse) > len(boundSubGenes):
		responseRate = float(len(GenesIntersection))/len(boundSubGenes)
		unboundResponsive = float(len(DESubGenes)-len(GenesIntersection))/(len(GenesUniverse)-len(boundSubGenes))
		return np.log2(responseRate/unboundResponsive)
	else:
		return "NA"

def computeSpearman(boundSubPVals,boundSubGenes,DESubPVals,DESubGenes,GenesIntersection):
	DEPVals,DEGenes = getSubset(DESubGenes,DESubPVals,GenesIntersection)
	BinPVals,BinGenes = getSubset(boundSubGenes,boundSubPVals,GenesIntersection)
	DEPVals,DEGenes = sortDataGenes(DEPVals,DEGenes)
	BinPVals,BinGenes = sortDataGenes(BinPVals,BinGenes)
	# pass in pvalues sorted alphabetically by gene name
	rho,pval = spearmanr(BinPVals,DEPVals)
	return(rho,pval)

# from https://dataconomy.com/2015/04/implementing-the-five-most-popular-similarity-measures-in-python/
def computeJaccardSimilarity(DEGenes,BoundGenes):
	intersection_cardinality = len(set.intersection(*[set(DEGenes), set(BoundGenes)]))
	union_cardinality = len(set.union(*[set(DEGenes), set(BoundGenes)]))
	if union_cardinality!=0:
		return intersection_cardinality/float(union_cardinality)
	else:
		return "NA"

def computeFDRLowerBound(GenesIntersection,GenesUniverse,DESubGenes,boundSubGenes):
	sensitivity = 0.8
	boundNum = len(boundSubGenes) - len(GenesIntersection)/sensitivity
	DENum = len(DESubGenes) - len(GenesIntersection)/sensitivity
	boundNum = max(0,boundNum)
	DENum = max(0,DENum)
	numerator = boundNum * DENum

	denominator = len(GenesUniverse) * len(GenesIntersection)

	return(numerator/denominator)