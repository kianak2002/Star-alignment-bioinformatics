def global_align(x, y, s_match, s_mismatch, s_gap):
    A = []
    for i in range(len(y) + 1):
        A.append([0] * (len(x) + 1))
    for i in range(len(y) + 1):
        A[i][0] = s_gap * i
    for i in range(len(x) + 1):
        A[0][i] = s_gap * i
    for i in range(1, len(y) + 1):
        for j in range(1, len(x) + 1):
            A[i][j] = max(
                A[i][j - 1] + s_gap,
                A[i - 1][j] + s_gap,
                A[i - 1][j - 1] + (s_match if (y[i - 1] == x[j - 1] and y[i - 1] != '-') else 0) + (
                    s_mismatch if (y[i - 1] != x[j - 1] and y[i - 1] != '-' and x[j - 1] != '-') else 0) + (
                    s_gap if (y[i - 1] == '-' or x[j - 1] == '-') else 0)
            )
    align_X = ""
    align_Y = ""
    i = len(x)
    j = len(y)
    while i > 0 or j > 0:
        current_score = A[j][i]
        if i > 0 and j > 0 and (
                ((x[i - 1] == y[j - 1] and y[j - 1] != '-') and current_score == A[j - 1][i - 1] + s_match) or
                ((y[j - 1] != x[i - 1] and y[j - 1] != '-' and x[i - 1] != '-') and current_score == A[j - 1][
                    i - 1] + s_mismatch) or
                ((y[j - 1] == '-' or x[i - 1] == '-') and current_score == A[j - 1][i - 1] + s_gap)
        ):
            align_X = x[i - 1] + align_X
            align_Y = y[j - 1] + align_Y
            i = i - 1
            j = j - 1
        elif i > 0 and (current_score == A[j][i - 1] + s_gap):

            align_X = x[i - 1] + align_X

            align_Y = "-" + align_Y

            i = i - 1

        else:

            align_X = "-" + align_X

            align_Y = y[j - 1] + align_Y

            j = j - 1

    return (align_X, align_Y, A[len(y)][len(x)])


def create_distance_matrix(seqs):
    distance_matrix = []
    for i in range(len(seqs)):
        temp = []
        for j in range(len(seqs)):
            if i != j:
                temp.append(global_align(seqs[i], seqs[j], 3, -1, -2))
        distance_matrix.append(temp)
    # print(distance_matrix)
    return distance_matrix


def choose_center(distance_matrix):
    seqs_score = []
    for i in range(len(distance_matrix)):
        sum = 0
        for j in range(len(distance_matrix[i])):
            # print(distance_matrix[i][j])
            sum += distance_matrix[i][j][2]
        seqs_score.append(sum)
    # print(seqs_score)
    max_seq = max(seqs_score)
    return seqs_score.index(max_seq)


def multi_align(seqs, center_index):
    for freq in range(len(seqs) - 1):
        for i in range(len(seqs)):
            if i != center:
                seq_new_1, seq_new_2, be_dard_nakhor = global_align(seqs[center], seqs[i], 3, -1, -2)
                seqs[center] = seq_new_1
                seqs[i] = seq_new_2
    return seqs


def score_MSA(seqs):
    sum_msa = 0
    for i in range(len(seqs)):
        for j in range(i+1, len(seqs)):
            sum_msa += score_two_seq(seqs[i], seqs[j])
    return sum_msa



def score_two_seq(seq1, seq2):
    score = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i] != '-':
            score += 3
        elif seq1[i] == seq2[i] == '-':
            score += 0
        elif (seq1[i] == '-' and seq2[i] != '-') or (seq1[i] != '-' and seq2[i] == '-'):
            score += -2
        else:
            score += -1
    return score


if __name__ == '__main__':
    n = input()  # how many sequences
    sequences = []  # an array for all sequences
    for i in range(int(n)):
        sequences.append(input())
    distance_matrix = create_distance_matrix(sequences)
    center = choose_center(distance_matrix)
    center_seq = sequences[center]
    sequences_new = multi_align(sequences, center)
    score = score_MSA(sequences_new)
    print(score)
    for sequence in sequences_new:
        print(sequence)
    # print(score_two_seq('WRYIAMRE-QYES--', '--YI-MQEVQQE--R'))
##### 13, 8, 16, 0, 7, 4