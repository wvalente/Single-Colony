from __future__ import print_function

# -*- coding: utf-8 -*-
"""
Created on Wed Feb 25 01:57:55 2015

@author: BJ Valente
"""

"""
HW5 - implementing a HMM

starting details:


"""
#!/usr/bin/env python

#################################PART 1############
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import copy

import math
import collections 
import itertools
    
A = SeqIO.read("CP003331.gbk", "genbank") #read in the Vibrio cholerae Chr. 2
seqA = Seq(str(A.seq), IUPAC.unambiguous_dna)

"""
Create a constructor class 'HMM' that contains method for setting probabilities
as well as performing the log-transform on those same probabilities. As in Durbin,
it's important to be using additive probabilities, rather than multiplicative, 
so the log-transform makes things additive. 

Probabilities are stored in dictionaries for fast, easy access later on. Better
mutability as well. 
"""
class HMM:
    def __init__(self, init_prob, trans_prob, emit_prob):
        self._name = 'Vibrio'        
        self.init_prob = {} #dictionary -  e.g. k = "1" , v = 0.996
        self.trans_prob = {} #dictionary - e.g. k = "(1,1)", v = 0.999  
        self.emit_prob = {} #dictionary -  e.g. k = "(1,'A'), v = 0.291
    
    def set_init(self, prob):
        (state, initial_prob) = prob #takes in an ordered pair      
        self.init_prob[state] = initial_prob
        assert sum(self.init_prob.values()) <= 1
    
    def set_trans(self, prob):
        ((state1,state2), transition_prob) = prob #takes in an ordered pair      
        self.trans_prob[(state1,state2)] = transition_prob
        assert sum(self.init_prob.values()) <= 1    
    
    def set_emit(self,prob):
        ((state,observation), emission_prob) = prob #takes in an ordered pair      
        self.emit_prob[(state,observation)] = emission_prob
        #assert sum(self.emit_prob.values()) <= 1
    
    #make probabilities additive, also "copy" to make efficient r/w characteristics
    #for the dictionaries that are transformed
    def log_transform(self,probability):
        log_prob = copy.copy(probability) #makes shallow copy
        try: 
            neg_inf = float("-inf") 
        except ValueError: 
             # On Python 2.5 or older that was handled in C code, 
             # and failed on Windows XP 32bit 
            neg_inf = - 1E400 
        for key in log_prob: 
            prob = log_prob[key] 
            if prob > 0: 
                log_prob[key] = math.log(log_prob[key]) 
            else: 
                log_prob[key] = neg_inf  
        return log_prob
    
    def __repr__(self):
        return 'HMM({!r:})'.format(self._name)
        
    """
    Employs the Viterbi algorithm on an HMM to readjust the transition probabilities
    without adjusting initial or emission probabilities
    Based on the explanations in class and Durbin et al. (esp. pgs 55-59)
    names of objects typically those used in Durbin pg. 56
    """    
    def ViterbiIterative(self, repeats = None): #modified itertools recipe   
        if repeats == 1:
            return self.Viterbi(seqA)
        if repeats != None:
            for i in range(1,repeats+1): 
                print('\n',"Iteration",i,":",'\n')
                self.Viterbi(seqA)
    
    def Viterbi(self, seqA):  
        #find improved parameter estimates for the transition probabilities
        #find the highest probability underlying state sequence
         
        V_probs = {} # viterbi probabilities
        pred_state_seq = collections.OrderedDict() #predicted state for a given sequence index i 
        possible_state_probs = {}        
        states = ['2', '1']
        log_prob_init = self.log_transform(self.init_prob)
        log_prob_emit = self.log_transform(self.emit_prob)
        log_prob_trans = self.log_transform(self.trans_prob)
        
        """
        Recursive portion: Vl(i) = El(Xi) + maxk(Vk(i-1) + Akl) (Durbin pg 77)
        Vl(i) == prob of the most prob path ending in state l w/ obs i
        El(Xi) == emission prob in state l of the observation Xi
        Vk(i-1) == prob of the most prob path ending in state k with obs i - 1
        Akl == trans prob of state k --> state l
        Note: sequence extends from 1 --> L in Durbin
        """
    #cur_state replaces l, prev_state replaces k
        for i in range(0, len(seqA)):        
            for cur_state in states: 
                ElXi = log_prob_emit[(cur_state, seqA[i])]                        
                max_prob = 0
                if i == 0:     #initialize! for i = 0 essentially
                    max_prob = log_prob_init[cur_state]
                else:
                    for prev_state in states:
                        Akl = log_prob_trans[(prev_state, cur_state)]
                        Vk_prev = V_probs[(prev_state, i - 1)]
                        cur_prob = Akl + Vk_prev                    
                        possible_state_probs[prev_state] = cur_prob                  
                    max_prob = max(possible_state_probs.values())
                V_probs[(cur_state, i)] = ElXi + max_prob  
                 
                if i > 0:
                    for state in possible_state_probs:
                        if possible_state_probs[state] == max_prob:
                            pred_state_seq[(i-1, cur_state)] = state
                            break
                        
        """
        Termination portion: P(x,most_probable_path) = maxk(Vk(L) + Ak0)
        """
        all_probs = {} # dict of all probs for all states as defined in v
        for state in states:
            all_probs[state] = V_probs[(state, len(seqA)-1)]
        
        state_path_prob = max(all_probs.values())
    
        last_state = ""
        for state in all_probs:
            if all_probs[state] == state_path_prob:
                last_state = state
    
        assert last_state != "", "Couldn't find the last state!"
    
        """
        Traceback: (i = L...1): most_prob_path(i-1) = pointer(most_prob_path(i))
        a pointer points "backwards" for observation == i on the path
        """
        traceback = []
        seq_nts = list(range(1, len(seqA)))
        seq_nts.reverse()  # L...1
        state = last_state
        traceback.append(state)        
        
        for i in seq_nts:
            state = pred_state_seq[(i-1), state]
            traceback.append(state)
        
        traceback.reverse() #flip it back to the right way, in place
        
        state_appearances = collections.Counter()
        segments = collections.Counter()
        state_2 = []
        for s in range(len(traceback)-1):
            state = traceback[s]
            next_state = traceback[s+1]    
            state_appearances[state] +=1
            if state == '1' and next_state == '2':
                state_2.append(s)
            if state == '2' and next_state == '1':
                state_2.append(s)
            if s == (len(traceback)-1):
                if state != next_state:
                    segments[state] +=1
                    segments[next_state] +=1
                else:
                    segments[state] +=1
            elif state != next_state:
                segments[state] +=1 
        """
        Recalculate 'improved' transition probabilities given results in segments
        and in state_appearances
        """
        new_trans11 = (state_appearances['1'] - segments[1])/state_appearances['1']
        new_trans12 = segments['1']/state_appearances['1']
        new_trans21 = segments['2']/state_appearances['2']
        new_trans22 = (state_appearances['2'] - segments[2])/state_appearances['2']
        
        self.set_trans((('1','1'),new_trans11))
        self.set_trans((('1','2'),new_trans12))
        self.set_trans((('2','1'),new_trans21))
        self.set_trans((('2','2'),new_trans22))
        """
        Output: # of states of each of the two types, number of segments
        new trans_prob to be used in the next iteration
        """
        print("State Histogram:")        
        print('1',"=",state_appearances['1'])
        print('2',"=",state_appearances['2'],'\n')
        print("Segment Histogram:")
        print('1','=', segments['1'])
        print('2','=', segments['2'],'\n')
        print("Transition Probabilities:")
        print("1,1", "=", new_trans11)
        print("1,2", "=", new_trans12)
        print("2,1", "=", new_trans21)
        print("2,2", "=", new_trans22,'\n')
        print("Segment List:")        
        args = [iter(state_2)]*2
        start_end = itertools.zip_longest(*args,fillvalue=None)
        for (s,e) in start_end:
            print(s,e)
        
        return self
        
    
VibrioHMM = HMM({},{},{})

VibrioHMM.set_init(("1",0.996))
VibrioHMM.set_init(("2",0.004)) 
for trans in [(('1','1'),0.999),(('1','2'),0.001), (('2','1'),0.01), (('2','2'),0.99)]:
    VibrioHMM.set_trans(trans) 

for emission in [(('1','A'),0.291),(('1','T'),0.291), (('1','C'),0.209), (('1','G'),0.209),\
                (('2','A'),0.169),(('2','T'),0.169), (('2','C'),0.331), (('2','G'),0.331)]:
    VibrioHMM.set_emit(emission)
            
            

with open('CP003331.fna', 'r') as fna_file:
    print(fna_file.readline())
    
VibrioHMM.ViterbiIterative(10)





