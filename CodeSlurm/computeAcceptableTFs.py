import sys
import pandas as pd
import argparse

def parse_args(argv):
  parser = argparse.ArgumentParser(description="")
  parser.add_argument("-d", "--data_file")
  parser.add_argument("-r", "--rand_file")
  parser.add_argument("-o", "--output_dir")
  parser.add_argument("-t", "--threshold", default = 0.01)
  parsed = parser.parse_args(argv[1:])
  return parsed

def computeCutoff(pVals, threshold):
	if len(pVals) < 10:
		return 0.00001
	position = int(round(len(pVals) * threshold))
	return sorted(pVals)[position]

def computeEdges(df):
	edges_df = pd.DataFrame(columns=["TF","Gene"])
	for _, row in df.iterrows():
		TF = row["TF"]
		genes = row["Genes"].split('\'')[1::2]
		for gene in genes:
			edges_df = edges_df.append(pd.Series({"TF": TF, "Gene": gene}),
										ignore_index=True)
	return edges_df


def main(argv):
	parsed = parse_args(argv)

	acceptableTFsIdx = []
	cutoffs_df = pd.DataFrame(columns=["TF", "HypergeometricPValCutoff"])
	data_df = pd.read_csv(parsed.data_file)
	rand_df = pd.read_csv(parsed.rand_file, names=["TF", "PVal"])
	for i, row in data_df.iterrows():
		TF = row["TF"]
		randPVals = rand_df.loc[rand_df["TF"] == TF, "PVal"].values
		cutoff = computeCutoff(randPVals, parsed.threshold)
		cutoffs_df = cutoffs_df.append(pd.Series({"TF": TF, 
									"HypergeometricPValCutoff": cutoff}), 
									ignore_index=True)
		print("%.2e, %.2e; %.2e, .2" % (row["HypergeometricPVal"], cutoff, row["FDR Lower Bound"]))
		if i > 4:
			sys.exit()
		if(row["HypergeometricPVal"] < cutoff and row["FDR Lower Bound"] <= 0.2):
			acceptableTFsIdx.append(i)

	acceptableTFs_df = data_df.iloc[acceptableTFsIdx]
	edges_df = computeEdges(acceptableTFs_df)
	targets_df = pd.DataFrame({"Network_targets": sorted(pd.unique(edges_df["Gene"]))})

	if not os.path.exists(parsed.output_dir):
		os.makedirs(parsed.output_dir)
	acceptableTFs_df.to_csv(parsed.output_dir + "/acceptableTFs.csv", index=False)
	cutoffs_df.to_csv(parsed.output_dir + "/TFcutoffs.csv", index=False)
	edges_df.to_csv(parsed.output_dir + "/edges.csv", index=False)
	targets_df.to_csv(parsed.output_dir + "/targets.csv", index=False)

	with open(parsed.output_dir+'/summary.txt','w') as f:
		f.write("Number of acceptableTFs: %d" % acceptableTFs_df.shape[0])
		f.write("\nNumber of Edges: %d" % edges_df.shape[0])
		f.write("\nNumber of Unique Targets: %s" % targets_df.shape[0])

if __name__ == "__main__":
	main(sys.argv)