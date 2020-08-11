import os

A = [0,0,0,1]
C = [0,0,1,0]
G = [0,1,0,0]
T = [1,0,0,0]
def toOneHot(seq,seqID):
    '''translate seq to ont-hot'''
    seq = seq.strip()
    list_hot = []
    for coden in seq:
        if coden == 'A':
            list_hot.append(A)
        elif coden == 'C':
            list_hot.append(C)
        elif coden == 'G':
            list_hot.append(G)
        elif coden == 'T':
            list_hot.append(T)
    return list_hot