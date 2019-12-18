import sys
import numpy as np
import pandas as pd
import scipy.stats as ss
import argparse


def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--input_file", 
						help="TF x gene matrix in csv format. First row and column represent genes and TFs.")
	parser.add_argument("-o","--output_file", 
						help="Gene x TF matrix in csv format with low rank edges removed.")
	parser.add_argument("-t","--rank_threshold", type=int, default=None)
	parsed = parser.parse_args(argv[1:])
	return parsed


def main(argv):
	parsed = parse_args(argv)
	df = pd.read_csv(parsed.input_file, index_col=0)
	m,n = df.shape

	df = np.abs(df)
	ranked = ss.rankdata(df, method="max").reshape((m,n))
	ranked = np.max(ranked) + 1 - ranked
	if parsed.rank_threshold:
		df[ranked > parsed.rank_threshold] = 0
	print("Found %d top edges" % np.count_nonzero(df.values))
	
	df.T.to_csv(parsed.output_file)


if __name__ == "__main__":
	main(sys.argv)