import hmm_module
import numpy as np
import pandas as pd
import argparse
import pdb, json

def get_dist_set(pos_list):
    dist_set = set()
    for positions in pos_list:
        dist = [positions[i+1]-positions[i] for i in range(0,len(positions)-1)]
        dist_set.update(dist)
    return dist_set


def posterior_decoding_fb(hmm, file_name, outfile, thresh, chromID, arcID): # for real data
	with open(file_name) as fin:
		header = fin.readline().strip().split()
	names = header[3:]

	# run posterior decoding for each individual
	for i in range(0, len(names)):
		print('analyzing %s' % names[i])
		hap = pd.read_csv(file_name, header=0, delimiter = "\t", usecols=[0, 1, 2, (i+3)], na_values=".")
		obs = np.zeros(hap.shape[0], dtype=int)
		for pos in range(0, hap.shape[0]):
			obs[pos] = int(hap.iloc[pos, 3] == 1 and hap.iloc[pos, 1] <= thresh and hap.iloc[pos, 2] > 0)
		
		path, probs = hmm.forward_backward_scaled(obs, np.array(hap.iloc[:,0]))
		starts, ends, l_probs = hmm.get_starts_ends(path, np.array(hap.iloc[:,0]), probs)
		for idx,start in enumerate(starts):
			oline = [chromID, starts[idx], ends[idx], l_probs[idx], names[i], arcID]
			outfile.write("\t".join([str(x) for x in oline]) + "\n")
		
	# close file
	outfile.close()


def init_hmm(r, t, m, file_name, thresh, e_11, e_01):
	# transition and emission probabilities 
	# 0: modern, 1: archaic
	a_01 = r*(t-1)*m
	a_10 = r*(t-1)*(1-m)
	prior = np.array([1-m,m])
	transition = np.array([[1-a_01, a_01],
						   [a_10, 1-a_10]])
	emission = np.array([[1-e_01, e_01],
						 [1-e_11, e_11]])
	
	# for precomputing transition matrices
	positions = pd.read_csv(file_name, header=0, delimiter = "\t", usecols=[0])
	l_positions = list(positions.iloc[:,0])
	dist_set = get_dist_set([l_positions])

	return hmm_module.HMM(prior, transition, emission, dist_set, thresh)

if __name__=='__main__':
	#params
	r = 1e-8 # recomb rate per bp per gen
	t = 1900 # admixture time in generations
	m = 0.05 # admixture proportion
	info_thresh = 0.05 
	hmm_thresh = 0.9 # call tract if posterior probability is higher than this value
	
	# Racimo et al. MBE 2017
#	e_11 = .050702 # from simulation
#	e_01 = 7.405e-4 # from simulation

	# SEGUIN-ORLANDO et al. Science 2014
	e_11 = 0.0155
	e_01 = 6.67e-9

	#read command line arguments
	parser = argparse.ArgumentParser(description="infer admixure tracts")
	parser.add_argument("-i", type=str, dest="i",help="absolute path of the input file")
	parser.add_argument("-chr", type=str, dest="chr",help="chromosome ID")
	parser.add_argument("-arc", type=str, dest="arc",help="archaic population ID")
	parser.add_argument("-o", type=str, dest="o",help="absolute path of the output path")
	parser.add_argument("-r", type=float, dest="r",help="recombination rate (per bp per generation)", default=r)
	parser.add_argument("-t", type=float, dest="t",help="admixture time (in generations)", default=t)
	parser.add_argument("-m", type=float, dest="m",help="admixture proportion", default=m)
	parser.add_argument("-it", type=float, dest="it",help="threshold at which to call a site informative", default=info_thresh)
	parser.add_argument("-ht", type=float, dest="ht",help="confidence level at which to call tracts", default=hmm_thresh)
	parser.add_argument("-json", type=str, dest="js",help="a json file with trained HMM parameters", default=None)
	args=parser.parse_args()	
	r = args.r
	t = args.t
	m = args.m
	info_thresh = args.it
	hmm_thresh = args.ht
	infile_path = args.i
	chromID = args.chr
	arcID = args.arc
	outfile_path = args.o
	jsonfile = args.js

	if jsonfile:
		with open(jsonfile) as f:
			data = json.load(f)
		m = data["starting_probabilities"][1]

	print('initializing HMM')
	hmm = init_hmm(r, t, m, infile_path, hmm_thresh, e_11, e_01)
	
	print('posterior decoding')
	outfile = open(outfile_path,'w')
	posterior_decoding_fb(hmm, infile_path, outfile, info_thresh, chromID, arcID)


