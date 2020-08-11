import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import pandas as pd

def CreateFinder(input_file):

    r_script = '''
    library("LncFinder")
    library("seqinr")
    Seqs <- read.fasta(file = "{inputfile}", seqtype = "DNA",as.string = TRUE)
    SS.seq_1 <- run_RNAfold(Seqs, RNAfold.path = "RNAfold", parallel.cores = 2)
    my_features <- extract_features(SS.seq_1, label = NULL,
                                    SS.features = TRUE, 
                                    format = "SS",
                                    frequencies.file = "human",
                                    parallel.cores = 2)
    write.csv(my_features, file = "finder_feature.csv")
    '''.format(inputfile = input_file)
    return r_script

def ExtractFinderF(input_file):
    r_script = CreateFinder(input_file)
    robjects.r(r_script)
    df = pd.read_csv("./finder_feature.csv", header=0)
    return df