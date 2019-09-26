import sys
import numpy as np
import pandas as pd
import scipy.stats as ss
import argparse


def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--input_file")
	parser.add_argument("-o","--output_file")
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
		ranked[ranked > parsed.rank_threshold] = m*n
	ranked = ranked / float(m*n)
	print("Found %d top edges" % np.sum(ranked < 1))
	
	df = pd.DataFrame(data=ranked, index=df.index, columns=df.columns).T
	df.to_csv(parsed.output_file)


if __name__ == "__main__":
	main(sys.argv)