import os
import sys
import Fea_Extract.toKmerRate as kr
import Fea_Extract.GetFa as gf

# def GetInputEntry(inputfile):
# #     """Get entry from input file"""
# #     SeqID,SeqList = gf.GetFasta(inputfile)
# #     return (SeqID,SeqList)

def GetFileName(inputfile):
    """Get name of inputfile as target variable quantity"""
    filename = inputfile.split('/')
    return filename[-1]

def CreateFeaFile(inputfile):
    """Create feature file of inputfile"""

    SeqID, SeqList = gf.GetFasta(inputfile)
    features = []
    for entryid,entry in zip(SeqID,SeqList):
        feature = kr.ExtractFea(entry,entryid)
        features.append(feature)
    return features