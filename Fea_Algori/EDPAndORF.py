import numpy as np
from Fea_Extract.GetFa import Codon2AA2

## AA list
_AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

## 2mer list
_DNA = ['A', 'C', 'G', 'T']
_Kmer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        _Kmer_list.append(dna1+dna2)

## 3mer list
_3mer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        for dna3 in _DNA:
            _3mer_list.append(dna1+dna2+dna3)

## IUPAC code
_IUPAC = {'A':'A','C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'M':'AC', 'K':'GT', 'S':'CG', 'W':'AT', 'H':'ACT', \
          'B':'CGT', 'V':'ACG', 'D':'AGT', 'N':'ACGT'}

def IUPAC_2mer(seq):
    '''Return a list of all possible 2mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            kmer_list.append(dna1+dna2)
    return kmer_list

def IUPAC_3mer(seq):
    '''Return a list of all possible 3mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            for dna3 in _IUPAC[seq[2]]:
                if Codon2AA2(dna1+dna2+dna3) != "J":
                    kmer_list.append(dna1+dna2+dna3)
    return kmer_list

def GetKmerEDP(seq):
    """Get feature of 2-mer EDP"""
    Kmer = {}
    for aa in _Kmer_list:
        Kmer[aa] = 1e-9

    sum_Kmer = 1e-9 * 16

    if(len(seq) > 3):
        for i in range(0,len(seq)-1) :
            ## IUPAC kmer
            if seq[ i:(i+2) ] not in _Kmer_list:
                tmp_kmer_list = IUPAC_2mer(seq[i:(i+2)])
                for tmp_kmer in tmp_kmer_list:
                    Kmer[ tmp_kmer ] += 1.0 / len(tmp_kmer_list)
            else:
                Kmer[ seq[ i:(i+2)] ] += 1.0
            sum_Kmer += 1.0

        H = 0.0
        for (k,v) in Kmer.items():
            Kmer[k] /= sum_Kmer
            Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
            H += Kmer[k]

        EDP = {}
        for (k,v) in Kmer.items():
            EDP[k] = Kmer[k] / H

        outline = ''
        for (k,v) in EDP.items():
            outline += str(v) + "\t"

        return outline

    else:
        return GetKmerEDP_Default()

def GetKmerEDP_Default():
    '''kmer entropy density'''
    Kmer = {}
    for aa in _Kmer_list:
        Kmer[aa] = 1e-9

    sum_Kmer = 1e-9 * 16

    H = 0.0
    for (k,v) in Kmer.items():
        Kmer[k] /= sum_Kmer
        Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
        H += Kmer[k]

    EDP = {}
    for (k,v) in Kmer.items():
        EDP[k] = Kmer[k] / H

    outline = ''
    for (k,v) in EDP.items():
        outline += str(v) + "\t"

    return outline

def GetEDP(seq, transcript_len):
    '''get features including: ORF length, ORF ratio, ORF EDP of codon'''

    Codon = {}
    for aa in _AA_list:
        Codon[aa] = 1e-9

    sum_codon = 1e-9 * 20

    if (len(seq) > 3):
        num = len(seq) // 3
        for i in range(0, num):
            if Codon2AA2(seq[i * 3:(i + 1) * 3]) == "J":
                continue
            ## IUPAC codon
            elif Codon2AA2(seq[i * 3:(i + 1) * 3]) == "Z":
                tmp_kmer_list = IUPAC_3mer(seq[i * 3:(i + 1) * 3])
                for tmp_kmer in tmp_kmer_list:
                    Codon[Codon2AA2(tmp_kmer)] += 1.0 / len(tmp_kmer_list)
                sum_codon += 1.0
            else:
                Codon[Codon2AA2(seq[i * 3:(i + 1) * 3])] += 1.0
                sum_codon += 1.0

        H = 0.0
        for (k, v) in Codon.items():
            Codon[k] /= sum_codon
            Codon[k] = -Codon[k] * np.log2(Codon[k])
            H += Codon[k]

        EDP = {}
        for (k, v) in Codon.items():
            EDP[k] = Codon[k] / H

        outline = ''
        for (k, v) in EDP.items():
            outline += str(v) + "\t"

        return outline

def GetEDP_noORF():
    '''return the default values'''

    Codon = {}
    for aa in _AA_list:
        Codon[aa] = 1e-9

    sum_codon = 1e-9 * 20

    H = 0.0
    for (k,v) in Codon.items():
        Codon[k] /= sum_codon
        Codon[k] = -Codon[k] * np.log2(Codon[k])
        H += Codon[k]

    EDP = {}
    for (k,v) in Codon.items():
        EDP[k] = Codon[k] / H

    outline = ''
    for i in range(20):
        outline += str(0) + "\t"

    return outline
