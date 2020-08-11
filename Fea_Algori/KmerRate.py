import numpy as np
from Cap_Tools.GetFa import GetFasta

_DNA = ['A', 'C', 'G', 'T']
_Kmer_list_ = []

_8mer_list_ = []
for dna1 in _DNA:
    for dna2 in _DNA:
        for dna3 in _DNA:
            for dna4 in _DNA:
                for dna5 in _DNA:
                    for dna6 in _DNA:
                        for dna7 in _DNA:
                            for dna8 in _DNA:
                                _8mer_list_.append(dna1 + dna2 + dna3 + dna4 + dna5 + dna6 + dna7 + dna8)

def getKmerRate(seq, seqid):

    print (_8mer_list_)
    seq = seq.strip()
    print (seq)
    _num_list_  = np.zeros(65536)
    print (len(_num_list_))
    for index in range(len(seq)):
        subseq = seq[index:index+8]
        if (len(subseq)) == 8:
            for listseq in _8mer_list_:
                if listseq == subseq:
                    _num_list_[_8mer_list_.index(listseq)] += 1
    _freq_list_ = []
    for chips in _num_list_:
        temp = chips/(len(seq)-7)
        _freq_list_.append(temp)

    # _freq_list1_ = [str(i) for i in _freq_list_]
    # feature = '\t'.join(_freq_list1_)
    # print (feature)
    # return feature
    print (_freq_list_)
    return _freq_list_