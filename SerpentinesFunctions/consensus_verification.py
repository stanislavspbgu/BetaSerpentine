from decimal import Decimal

import sys
sys.setrecursionlimit(5000)
import math


def fisher_test(NStrands, NArcs, expNStrands, expNArcs):
    a = (math.factorial(NStrands + NArcs) * math.factorial(expNArcs + expNStrands) * math.factorial(NStrands + expNStrands) * math.factorial(NArcs + expNArcs))
    b = math.factorial(NArcs + NStrands + expNArcs + expNStrands) * math.factorial(NArcs) * math.factorial(NStrands) * math.factorial(expNArcs) * math.factorial(expNStrands)
    return float(Decimal(a) / Decimal(b))

def verify_consensus(NStrands, NArcs, sequence):
    random_frequences = {'A' : 0.619, 'C' : 0.624, 'D' : 0.421, 'E' : 0.453, 'F' : 0.591,
                         'G' : 0.576, 'H' : 0.599, 'I' : 0.623, 'K' : 0.442, 'L' : 0.631,
                         'M' : 0.617, 'N' : 0.635, 'P' : 0.035, 'Q' : 0.635, 'R' : 0.431,
                         'S' : 0.613, 'T' : 0.594, 'V' : 0.659, 'W' : 0.616, 'Y' : 0.612}
    output_sequence = ''
    output_p_value = ''
    output_structure = ''
    consensus = {'strands' : [], 'arcs' : []}
    for i in range(len(sequence)):
        output_sequence += ' ' + sequence[i] + ' '
        if NArcs[i] + NStrands[i] > 0 and sequence[i] in random_frequences: ####modificattion
            NSerp = NArcs[i] + NStrands[i]
            Sfreq = random_frequences[sequence[i]]
            expNStrand = int(Sfreq * NSerp)
            expNArc = int((1-Sfreq) * NSerp)
            p = fisher_test(NStrands[i], NArcs[i], expNStrand, expNArc)
            if p < 0.001:
                output_p_value += (' * ')
            else: output_p_value += ('   ')

            if NStrands[i] >= NArcs[i]:
                output_structure += (' S ')
                if p < 0.001: consensus['strands'].append(i + 1)
            else:
                output_structure += (' A ')
                if p < 0.001: consensus['arcs'].append(i + 1)
        else:
            output_structure += (' - ')
            output_p_value += ('   ')
    #print consensus
    return [output_p_value, output_structure, output_sequence, consensus]

def CompareToConsensus(consensus, serpentine_structure):

    CStrands = consensus['strands']
    #print 'CStrands ', CStrands
    CArcs = consensus['arcs']

    Strands = []
    Arcs = []
    for i in range(len(serpentine_structure)):
        if i % 2 == 0:
            Strands += serpentine_structure[i]
        else: Arcs += serpentine_structure[i]

    Intersection = list(set(CArcs) & set(Arcs)) #list(set(a) & set(b))
    Index = float(len(Intersection))/len(Arcs)
    return Index

