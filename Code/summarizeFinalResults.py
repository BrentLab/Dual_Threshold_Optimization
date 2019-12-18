import sys
import os
import argparse


def parse_args(argv):
	parser = argparse.ArgumentParser(description="")
	parser.add_argument("-i","--input_dirpath", required=True,
						help="Directory of DTO results including both authentic model and randomized models.")
	parser.add_argument("-l","--tfs_filepath", default=None,
						help="File of TF list, for which will be analyzed, e.g. if a particular family of TFs are of interest.")
	parser.add_argument("-e","--user_email", default=None,
						help="Email address for SLURM notification")
	parser.add_argument("--run_local", action='store_true', default=False,
						help="Flag for running DTO in serial fashion on local machine.")
	parsed = parser.parse_args(argv[1:])
	return parsed
	

def generateSbatch(input_dirpath, tfs_filepath):
	input_name = os.path.basename(input_dirpath.strip("/"))
	sbatch_filepath = "tmp/summarize_{}.sh".format(input_name)
	if tfs_filepath is None:
		tf_list = ""
	else:
		tf_list = " -l {}".format(tfs_filepath)
	with open(sbatch_filepath, "w") as f:
		f.write("#!/bin/bash\n#SBATCH --mem=4G\n#SBATCH -J sumDTO\n#SBATCH -o tmp/log_{0}.out\n#SBATCH -e tmp/log_{0}.err\n".format(input_name))
		f.write("python -u compileResults.py -i {0}/authentic_model/ -o {0}/authentic_model.csv\npython -u compileResults.py -i {0}/random_models/ -o {0}/random_models.csv -r split\npython computeAcceptableTFs.py -d {0}/authentic_model.csv -r {0}/random_models.csv -o {0}/summary/{1}\ncat {0}/summary/summary.txt\n".format(input_dirpath, tf_list))
	return sbatch_filepath


def submitSbatch(sbatch_filepath, user_email, run_local=False):
	if run_local:
		cli_job = "bash {}".format(sbatch_filepath)
	else:
		if user_email is None:
			cli_job = "sbatch {}".format(sbatch_filepath)
		else:
			cli_job	= "sbatch --mail-type=END,FAIL --mail-user={0} {1}".format(user_email, sbatch_filepath)
	os.system(cli_job)


def main(argv):
	parsed = parse_args(argv)
	sbatch_filepath = generateSbatch(parsed.input_dirpath, parsed.tfs_filepath)
	submitSbatch(sbatch_filepath, parsed.user_email, parsed.run_local)


if __name__ == "__main__":
	main(sys.argv)