from Fea_Algori import FickettFea as ft, Hexamer as ha, EDPAndORF as fe, UTRCovAndGCCont as ug
import os

def ExtractFea(seq, seqid, log_hexammer=None):
    if log_hexammer == None:
        dir = os.path.dirname(os.path.abspath(__file__))
        log_hexammer = os.path.join(dir,"../Fea_Algori/human.hexamer.logscore")
    logscore_dict = ha.ReadLogScore(log_hexammer)
    return GetFeature(seq, seqid, logscore_dict)

def GetFeature(seq, seqid, logscore_dict, ORFFea=True, EdpFea=True,
               FicFea=True, EDPkmerFea=True, UTRkmerFea=True, HexFea=True):

    # search ORF
    seq = seq.strip()
    ORF, UTR5, UTR3 = ug.GetORF_UTR(seq)
    transcript_len = len(seq)

    # ORF features
    if ORFFea:
        ORF_fea = str( len(ORF) ) + "\t" + str( len(ORF) * 1.0 / transcript_len )
    else:
        ORF_fea = ""

    # The EDP of ORF
    if EdpFea :
        if len(ORF) < 6:
            EDP_fea = fe.GetEDP_noORF()
        else:
            EDP_fea = fe.GetEDP(ORF, transcript_len)
    else:
        EDP_fea = ""

    if EDPkmerFea :
        Kmer_EDP_fea = fe.GetKmerEDP(ORF)
    else:
        Kmer_EDP_fea = ""

    # UTR feature
    if UTRkmerFea :
        UTR5_fea = str(len(UTR5) * 1.0 / len(seq)) + "\t" + str(ug.GetGC_Content(UTR5))
        UTR3_fea = str(len(UTR3) * 1.0 / len(seq)) + "\t" + str(ug.GetGC_Content(UTR3))
    else:
        UTR5_fea = ""
        UTR3_fea = ""

    # Hexamer feature
    ave_hexamer_score = ha.HexamerScore2(ORF, logscore_dict)
    if HexFea :
        hex_score = str(ave_hexamer_score)
    else:
        hex_score = ""

    # Fickett feature
    if FicFea :
        A_pos_fea = ft.GetBasePositionScore(seq, 'A')
        C_pos_fea = ft.GetBasePositionScore(seq, 'C')
        G_pos_fea = ft.GetBasePositionScore(seq, 'G')
        T_pos_fea = ft.GetBasePositionScore(seq, 'T')
        base_ratio = ft.GetBaseRatio(seq)

        fickett_fea = "\t".join([str(A_pos_fea), str(C_pos_fea), str(G_pos_fea), str(T_pos_fea), str(base_ratio[0]), str(base_ratio[1]), str(base_ratio[2]), str(base_ratio[3])])
    else:
        fickett_fea = ""

    feature = "\t".join([EDP_fea, ORF_fea, Kmer_EDP_fea, UTR5_fea, UTR3_fea, hex_score, fickett_fea])

    return feature

