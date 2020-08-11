import Fea_Extract.GetFa as gf
import pandas as pd
import FinderF as ff
import Fea_Extract.CreateFeaTable as cft
import Autoencoder as ae
import Fea_Extract.CreateTupleTable as ctt

#
# # os.system("RNAfold.exe  <test.fa   >result.txt")
#
# r_script = '''
# library("LncFinder")
# library("seqinr")
# Seqs <- read.fasta(file = "E:/workspace/FengSY/PyCharmProjects/Core_lnclocation/test.fa", seqtype = "DNA",as.string = TRUE)
# RNAfold.path <- '"E:/Program Files (x86)/ViennaRNA Package/RNAfold.exe"'
# SS.seq_1 <- run_RNAfold(Seqs[1:2], RNAfold.path = RNAfold.path, parallel.cores = 2)
# my_features <- extract_features(SS.seq_1, label = NULL,
# 	                            SS.features = TRUE,
# 			                    format = "SS",
#                 	            frequencies.file = "human",
#                                 parallel.cores = 2)
# write.csv(my_features, file = "finder_feature.csv")
# '''

# r_script = '''
# library("LncFinder")
# Seqs <- "E:/workspace/FengSY/PyCharmProjects/Core_lnclocation/test_data/test.fa"
# RNAfold.path <- '"E:/Program Files (x86)/ViennaRNA Package/RNAfold.exe"'
# SS.seq_1 <- run_RNAfold(Seqs[1:2], RNAfold.path = RNAfold.path, parallel.cores = 2)
# '''

# finder = importr('LncFinder')
# result = r["extract_features"](temp_data)
# print (result)

# robjects.r(r_script)
#
#

# for entryid, entry in zip(SeqID, SeqList):
#     feature = fe.ExtractFea(entry, entryid)
#
# print (feature)
#

test_data = "./test.fa"
SeqID,SeqList = gf.GetFasta(test_data)

def SeqFeature(seq_data):

    # finder
    df1 = ff.ExtractFinderF(input_file=seq_data)
    df1.columns = ['0', '1', '2', '3', '4',
                   '5', '6', '7', '8', '9',
                   '10', '11', '12', '13', '14',
                   '15', '16', '17', '18', '19',]
    # print(df1)

    # Bio Feature
    feature_bio = cft.CreateFeaFile(seq_data)
    df2 = pd.DataFrame(feature_bio)
    df2 = df2.drop([21,40],axis=1)
    seqID = df2[0]
    # df2.rename(columns={'0' : 'SeqID'}, inplace=True)
    # print(df2)
    # df2.to_csv("./test.csv",index=0)

    df1_temp = df1.drop('0',axis=1)
    df = pd.concat([df2, df1_temp], axis=1, ignore_index=True)
    # print(df)
    # df.to_csv("./test.csv",index=0)

    # advanced feature
    encoder = ae.getTemFeature(df)
    feature_temp = pd.concat([df, encoder], axis=1, ignore_index=True)
    # feature_temp.to_csv("./test11.csv", index=0)
    # print(feature_temp)

    # RFE feature
    sortRFE = open('./script/sortRFE.txt')
    RFE_matrix = []

    for eachsort in sortRFE:
        eachsort = eachsort.split('\n')
        RFE = eachsort[0]
        RFE_matrix.append(RFE)

    sortRFE.close()

    # 字符转换为整型
    new_RFE_matrix = []
    for n in RFE_matrix:
      new_RFE_matrix.append(int(n))

    # print(new_RFE_matrix)

    feature_RFE = feature_temp.iloc[:,new_RFE_matrix]
    # print (feature_RFE)

    # Tuple feature
    feature_tuple = ctt.CreateFeaFile(seq_data)
    df3 = pd.DataFrame(feature_tuple)
    # print(df3)

    # merge all features
    feature_all = pd.concat([df3, feature_RFE], axis=1, ignore_index=True)
    # print (feature_all.shape)

    return seqID, feature_all