def GetBasePositionScore(seq, base):
    '''compute base(A, C, G, T) position score based
       on Fickett's methods
    '''
    base_num = [0, 0, 0]
    tmp = len(seq) // 3

    for i in range(0, tmp):
        for j in range(0, 3):
            if seq[j + 3*i] == base:
                base_num[j] += 1

    base_pos_score = max(base_num) * 1.0 / (min(base_num) + 1)

    return base_pos_score

def GetBaseRatio(seq):
    '''calculate the A, C, G, T percentage of ORF'''
    A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9
    for i in range(0, len(seq)):
        if seq[i] == 'A':
            A += 1
        elif seq[i] == 'C':
            C += 1
        elif seq[i] == 'G':
            G += 1
        elif seq[i] == 'T':
            T += 1

    A = A * 1.0 / (len(seq) + 4e-9)
    C = C * 1.0 / (len(seq) + 4e-9)
    G = G * 1.0 / (len(seq) + 4e-9)
    T = T * 1.0 / (len(seq) + 4e-9)

    return [A, C, G, T]
