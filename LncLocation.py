# -*- coding:utf-8 -*-

"""
Author : Feng Shiyao
Usage: Please refer to the README file
Intro : A local application for predicting the subcellular location of lncRNAs
"""

import os
import sys
import numpy as np
from sklearn.externals import joblib
from GetFeature import SeqFeature
from sklearn import preprocessing

class_dict = {'0': 'cytoplasm', '1': 'exosome', '2': 'nucleus', '3': 'ribosome'}


def detecting_python_version_AND_output_help_info():
    if sys.version[0] == '2':
        print("""\nVersionError: The current python version: '%s',
                    You must use python version 3 or later to run this script!\n""" % ((sys.version).split('\n')[0]))
        exit(0)
    else:
        pass

    try:
        if sys.argv[1] == "--help":
            print_help_info()
        else:
            pass
    except:
        print_help_info()


def print_help_info():
    print("""
    Usage: python LncLocation.py input_file output_profix RNAfold_path

    Input_file must be fasta format

    Users can use the output_profix key value to specify the output folder of the predicted results file

    A demo program is ready in the package and the user can test it with instructions after installing RNAfold:
        python LncLocation.py ./test.fa result
    The test results will be saved in the result folder

    Other details, please read the README file!
""")
    exit(0)


def GetFileName(inputfile):
    """Get name of inputfile as target variable quantity"""
    filename = inputfile.split('/')
    return filename[-1]


def GetTempfile(inputfile):
    f1 = open('./script/paticipate', 'r')
    temp_seq = f1.read()
    f1.close()

    with open(inputfile, 'r') as f:
        data0 = f.read()
        data = data0.strip()
        f.close()
        # print(data)

    tmp_path = os.path.join(os.path.abspath('.'), 'script')
    if not os.path.exists(tmp_path):
        os.mkdir(tmp_path)
    temp_file = os.path.join(tmp_path, "temp_" + GetFileName(inputfile))
    with open(temp_file, 'w') as tmp_data:
        tmp_data.write(data + '\n' + temp_seq)
        tmp_data.close()
    return temp_file, GetFileName(inputfile)


def CreatePredictFile(inputfile, output_folder, profix):
    """Create feature file of inputfile"""
    tmpdir = os.path.join(os.path.abspath('.'), output_folder)
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    result_file = os.path.join(tmpdir, "Predict_" + profix)
    with open(result_file, 'w') as fout:
        minmax = preprocessing.MinMaxScaler()
        seqID, feature = SeqFeature(inputfile)
        features = minmax.fit_transform(feature.values)
        clf = joblib.load("./svm_predictor.m")
        result = clf.predict_proba(features)
        count = 0
        for entryid, res in zip(seqID, result):
            count += 1
            if count == len(seqID):
                break
            fout.write(entryid + ':' + '\n\t'
                       + 'cytoplasm:' + str(res[0]) + '\t'
                       + 'exosome:' + str(res[1]) + '\t'
                       + 'nucleus:' + str(res[2]) + '\t'
                       + 'ribosome:' + str(res[3]) + '\n\t'
                       + "predicted results:" + class_dict[str(np.argmax(res))] + '\n')
            print(entryid + ':' + class_dict[str(np.argmax(res))] + '\n')
        fout.close()
        print("Performing completed! Please check output folder!")


def GetPredictFile(inputfile, output_profix):
    temp_file, profix = GetTempfile(inputfile)
    CreatePredictFile(inputfile=temp_file, output_folder=output_profix, profix=profix)


detecting_python_version_AND_output_help_info()
inputfile = sys.argv[1]
output_profix = sys.argv[2]

# inputfile = "./test.fa"
# output_profix = "result"

if __name__ == '__main__':
    GetPredictFile(inputfile, output_profix)