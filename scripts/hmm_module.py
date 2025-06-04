import numpy as np
import pdb
#from scipy.optimize import fmin_l_bfgs_b as minimize
#import time
#import copy

class HMM:

    def __init__(self, prior, transition, emission, dist_set,thresh): 
        self.prior = prior
        self.transition = transition
        self.emission = emission
        self.trans_dict = self.compute_trans_mats(dist_set,self.transition)
        self.thresh = thresh
        
    def __del__(self):
        del self.prior
        del self.transition
        del self.emission
        del self.trans_dict

    def compute_trans_mats(self,dist_set,A):
        trans_dict = dict()
        for dist in dist_set:
            trans_dict[dist] = np.linalg.matrix_power(A,dist)
        return trans_dict

    def forward_backward_scaled(self, observed_states, positions):
        # Initialize
        num_hidden_states = self.transition.shape[0]
        num_sites = len(positions)
        alpha_table = np.zeros([num_hidden_states,num_sites])
        beta_table = np.zeros([num_hidden_states,num_sites])
        gamma_table = np.zeros([num_hidden_states,num_sites])
        path = np.zeros(num_sites,dtype=int)
        probs = np.zeros(num_sites)
        alpha_scales = np.zeros(num_sites)
        beta_scales = np.zeros(num_sites)

        # Alpha pass
        alpha_table[:,0] = self.prior * self.emission[:,observed_states[0]]
        alpha_scales[0] = max(alpha_table[:,0])
        alpha_table[:,0] = alpha_table[:,0]/alpha_scales[0]
        for t in range(1,num_sites): 
            trans_mat = self.trans_dict[positions[t]-positions[t-1]]
            for s in range(0,num_hidden_states):
                alpha_table[s,t] = sum(alpha_table[:,t-1] * trans_mat[:,s]) * self.emission[s,observed_states[t]]
            alpha_scales[t] = max(alpha_table[:,t])
            alpha_table[:,t] = alpha_table[:,t]/alpha_scales[t]
        
        # Beta pass
        beta_table[:,-1] = 1
        beta_scales[-1] = 1
        for t in range(num_sites-2,-1,-1): 
            trans_mat = self.trans_dict[positions[t+1]-positions[t]]
            for s in range(0,num_hidden_states):
                beta_table[s,t] = sum(trans_mat[s,:] * self.emission[:,observed_states[t+1]] * beta_table[:,t+1])
            beta_scales[t] = max(beta_table[:,t])
            beta_table[:,t] = beta_table[:,t]/beta_scales[t]

        # Compute gammas
        for t in range(0,num_sites):
            denom = sum(alpha_table[:,t] * beta_table[:,t])
            gamma_table[:,t] = (alpha_table[:,t] * beta_table[:,t])/denom

        # Reconstruct best path from gamma table
        path[np.where(gamma_table[1,:]>self.thresh)[0]] = 1
        for i in range(0,len(probs)):
            probs[i] = gamma_table[path[i],i]

        return path, probs
    
    def stitch_tracts(self, path, length_thresh, dist_thresh):
        num_sites = len(path)
        # compute starts and ends
        starts, ends = np.array([]),np.array([])
        prev,curr = 0,0
        for i in range(0,num_sites):
            curr = path[i]
            if curr==1 and prev==0:
                starts = np.append(starts,i)
            elif (curr==0 and prev==1) or (i==num_sites-1 and curr==1):
                ends = np.append(ends,i-1)
            prev = curr

        # stitch short tracts
        stitched_starts, stitched_ends = np.array([]), np.array([])
        stitched_starts = np.append(stitched_starts,starts[0])
        for i in range(1,len(starts)):
            prev_end = ends[i-1]
            curr_start = starts[i]
            if curr_start-prev_end > dist_thresh:
                stitched_ends = np.append(stitched_ends,prev_end)
                stitched_starts = np.append(stitched_starts, curr_start)
        stitched_ends = np.append(stitched_ends, ends[i]) # add the last end

       # remove short tracts
        lengths = stitched_ends-stitched_starts
        bad_ind = np.where(lengths<length_thresh)[0]
        stitched_starts = np.delete(stitched_starts, bad_ind)
        stitched_ends = np.delete(stitched_ends, bad_ind)

        # stitched path
        stitched_path = np.zeros(num_sites,dtype=int)
        for i in range(0,len(stitched_starts)):
            stitched_path[stitched_starts[i]:stitched_ends[i]+1] = 1
        
        return stitched_path

    def get_starts_ends(self, path, pos, probs):
        num_sites = len(path)
        starts, ends, l_probs = np.array([], dtype=int), np.array([], dtype=int), np.array([], dtype=float)
        tmp_l_probs = np.array([], dtype=float)

        prev,curr = 0,0
        for i in range(0,num_sites):
            curr = path[i]
            
            if curr == 1:
                tmp_l_probs = np.append(tmp_l_probs, probs[i])

            if (prev == 0 and curr == 1):
                starts = np.append(starts, pos[i])
            elif (prev == 1 and curr == 0):
                ends = np.append(ends, pos[i-1])
                l_probs = np.append(l_probs, np.median(tmp_l_probs))
                tmp_l_probs = np.array([], dtype=int)

            # test if reach to the last variant
            if (i==num_sites-1 and curr==1):
                ends = np.append(ends, pos[i])
                l_probs = np.append(l_probs, np.median(tmp_l_probs))

            prev = curr

        return starts, ends, l_probs

