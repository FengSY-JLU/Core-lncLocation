# -*- coding:utf-8 -*-

import re
import os
import sys
import multiprocessing
import Fea_Extract.GetFa as gf


# 计算
def Calculate_fre(seq):
    eightNC_num = len(seq) - 8 + 1
    eighttuple_matrix_dict = dict()
    eighttuple_matrix_dict = eighttuple_matrix_dict.fromkeys(eighttuple_matrix, 0)

    for each_seq_th in range(eightNC_num):
        octNC = seq[each_seq_th:each_seq_th + 8]
        try:
            eighttuple_matrix_dict[octNC] += 1
        except:
            pass
    feature = []
    numth = 0
    for eighttuple_th in eighttuple_matrix:
        numth += 1
        fre = eighttuple_matrix_dict[eighttuple_th] / float(eightNC_num)
        if numth != 4107:
            feature.append('%.6f' % fre)
        else:
            feature.append('%.6f' % fre)
    return feature


##获得排序联体特征
sorteighttuple = open('./script/sort4107.txt')
eighttuple_matrix = []

for eachsort in sorteighttuple:
    eachsort = eachsort.split('	')
    eigthtuple = eachsort[1]
    eighttuple_matrix.append(eigthtuple)

sorteighttuple.close()

def ExtractFea(seq, seqid):
    return Calculate_fre(seq)