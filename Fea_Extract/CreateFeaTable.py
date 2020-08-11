import os
import sys
import Fea_Extract.FeaExtract as fe
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
        feature = fe.ExtractFea(entry,entryid)
        feature_str = entryid + '\t' + feature
        features.append(feature_str.split('\t'))
    return features


