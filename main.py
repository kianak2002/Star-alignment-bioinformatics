import numpy as np


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
                # for m in range(len(temp)):
                #    temp[m][2] = distance(temp[m][0], temp[m][1])

        distance_matrix.append(temp)
    # print(distance_matrix)
    return distance_matrix


def distance(seq1, seq2):
    match = 0
    mismatch = 0
    for i in range(len(seq1)):
        if seq1[i] == seq2[i] != '-':
            match += 0
        if seq1[i] != '-' and seq2[i] != '-' and seq2[i] != seq1[i]:
            mismatch += 0
    return mismatch / match


'''
it has the distance matrix
chooses the one with the max score as the center
'''


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


'''
align all sequences one by one using their alignment with the center sequence
'''


def multi_align(seqs, center_index):
    list_first_aligns = []
    for k in range(len(seqs)):
        if k != center_index:
            align1, align2, be_dard_nakhor = global_align(seqs[center_index], seqs[k], 3, -1, -2)
            list_first_aligns.append([align1, align2])  # list of all aligns with the center sequence
    # print(list_first_aligns)

    center_seq = list_first_aligns[0][0]
    align = []
    align.append(list_first_aligns[0][0])
    align.append(list_first_aligns[0][1])
    # print(align)
    for i in range(len(seqs) - 2):
        # print(center_seq, 'before')
        center_seq, next_align = change_aligned_sequences(list_first_aligns[i + 1], center_seq)
        # print(center_seq, ' after')
        for j in range(len(align) - 1):
            not_useful, changed_align = change_aligned_sequences([align[0], align[j + 1]], center_seq)
            align[j + 1] = changed_align
        align[0] = center_seq
        align.append(next_align)
        # print(align)
        # print(align, 'align', i)
    return align

    # for freq in range(len(seqs) - 1):
    #     for i in range(len(seqs)):
    #         if i != center_index:
    #             seq_new_1, seq_new_2, be_dard_nakhor = global_align(seqs[center_index], seqs[i], 3, -1, -2)
    #             seqs[center_index] = seq_new_1
    #             seqs[i] = seq_new_2
    # return seqs


'''
it gets the aligned sequence with center and the new center
returns the new center and new sequence
'''


def change_aligned_sequences(seq_aligned, seq_center):
    seq_new = []
    center_new = []
    j = 0
    z = 0
    diff_center = abs(len(seq_aligned[0]) - len(seq_center))
    count1 = len(seq_center)
    count2 = len(seq_aligned[0])

    for i in range(200000):
        if j >= len(seq_aligned[0]):
            center_new.append(seq_center[z])
            seq_new.append(seq_center[z])
            count1 -= 1
        elif z >= len(seq_center):
            center_new.append(seq_aligned[0][j])
            seq_new.append(seq_aligned[1][j])
            count2 -= 1
        elif seq_center[z] == '-' != seq_aligned[0][j]:
            seq_new.append('-')
            center_new.append('-')
            count1 -= 1
        elif seq_aligned[0][j] == '-' != seq_center[z]:
            center_new.append('-')
            seq_new.append(seq_aligned[1][j])
            z -= 1
            j += 1
            count2 -= 1
        else:
            seq_new.append(seq_aligned[1][j])
            center_new.append(seq_center[z])
            j += 1
            count1 -= 1
            count2 -= 1
        z += 1

        if count1 == count2 == 0:
            break

    center_new = ''.join([str(elem) for elem in center_new])
    seq_new = ''.join([str(elem) for elem in seq_new])

    # print(center_new, seq_new, 'change_align')
    return center_new, seq_new


'''
calculate the score for each pair and add them together
'''


def score_MSA(seqs):
    sum_msa = 0
    for i in range(len(seqs)):
        for j in range(i + 1, len(seqs)):
            # print(score_two_seq(seqs[i], seqs[j]), seqs[i], seqs[j])
            sum_msa += score_two_seq(seqs[i], seqs[j])
    return sum_msa
    # for i in range(len(seqs)):
    #     for j in range(len(seqs[0])):
    #         sum +=


'''
the score for each two characters
'''


def score_block(first, second):
    score = 0
    if first == second == '-':
        score += 0
    elif first == second != '-':
        score += 3
    elif (first == '-' and second != '-') or (first != '-' and second == '-'):
        score += -2
    else:
        score += -1
    return score


'''
gets two sequences
calculate their score 
'''


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


def sort_aligns(center_index, last_list):
    new_list = []
    for i in range(len(last_list)):
        if i == center_index:
            new_list.append(last_list[0])
        elif i < center_index:
            new_list.append(last_list[i + 1])
        else:
            new_list.append(last_list[i])
    return new_list


def block_multi_align(blocks):
    distance_matrix_b = create_distance_matrix(blocks)
    center_block = choose_center(distance_matrix_b)
    block_new = multi_align(blocks, center_block)
    score = score_MSA(block_new)
    block_new = sort_aligns(center_block, block_new)
    return block_new, score


def choose_block(seqs, score_msa):
    last_score = score_msa
    for i in range(len(seqs[0]), 1, -1):
        for k in range(len(seqs[0]) - i + 1):
            block = []
            for j in range(len(seqs)):
                block.append(seqs[j][k:k + i])
            if check_block_improvement(block):
                new_blocks = remove_gaps(block)
                aligned_block, score = block_multi_align(new_blocks)
                # print(block, 'previous chosen block')
                # print(aligned_block, ' changed block')
                # if block != aligned_block or block[0][0] == 'E':
                    # print("score jadide:", score, "score ghadimie:", score_MSA(block))
                # print(aligned_block ,'khar\n', block)
                # print(score, score_MSA(block), "hehhhh")
                if score > score_MSA(block):
                    for n in range(len(seqs)):
                        new_seq = seqs[n][:k - 1] + aligned_block[n] + seqs[n][k + i + 1:]
                        # seqs[n][k:k+i] = aligned_block[n]
                        seqs[n] = new_seq
        last_score = score_MSA(seqs)
    return seqs, last_score


'''
checks if the chosen block is an improving block 
'''


def check_block_improvement(block):
    # block_np = np.array(block)
    # print(block_np, "bbbbbbb")
    # block_T = block_np.transpose()
    # print(block_T, "ahhhh")
    # for i in range(len(block_T)):
    #     if max(block_T[i]) == min(block_T[i]):  # to check if whole row is equal
    #         return False
    for i in range(len(block[0])):
        check = False
        for j in range(len(block)-1):
            if block[j][i] != block[j+1][i]:
                check = True
        if not check:
            return False
    return True


def remove_gaps(block):
    new_block = []
    for i in range(len(block)):
        temp = []
        for j in range(len(block[i])):
            if block[i][j] != '-':
                temp.append(block[i][j])
        temp = ''.join([str(elem) for elem in temp])
        new_block.append(temp)

    return new_block


if __name__ == '__main__':
    n = input()  # how many sequences
    sequences = []  # an array for all sequences
    for i in range(int(n)):
        sequences.append(input())
    distance_matrix = create_distance_matrix(sequences)
    center = choose_center(distance_matrix)
    # center_seq = sequences[center]
    sequences_new = multi_align(sequences, center)
    score = score_MSA(sequences_new)
    # print(score)
    # print(center_seq)
    sequences_new = sort_aligns(center, sequences_new)
    # best = block_replace(sequences, score)
    # sequences_new = block_multi_align(sequences)
    # for sequence in sequences_new:
    #     print(sequence)

    seqs_new, score_new = choose_block(sequences_new, score)
    print(score_new)
    for sequence in seqs_new:
        print(sequence)
    # print(global_align('FFVSANPW', 'NAHTAFL', 3, -1, -2))
    # print(global_align('EKKITGYTT', 'NAHTAFL', 3, -1, -2))
    # print(global_align('LNELQQ', 'NAHTAFL', 3, -1, -2))
    # print(global_align('EKKITGYTT', 'FFVSANPW', 3, -1, -2))
    # print(global_align('EKKITGYTT', 'LNELQQ', 3, -1, -2))
    # print(global_align('EKKITGYTT', 'NAHTAFL', 3, -1, -2))
    # print(global_align('FFVSANPW', 'EKKITGYTT', 3, -1, -2))
    # print(global_align('FFVSANPW', 'LNELQQ', 3, -1, -2))
    # print(global_align('FFVSANPW', 'NAHTAFL', 3, -1, -2))
    # print(global_align('LNELQQ', 'EKKITGYTT', 3, -1, -2))
    # print(global_align('LNELQQ', 'FFVSANPW', 3, -1, -2))
    # print(global_align('LNELQQ', 'NAHTAFL' , 3, -1, -2))

    # choose_block(sequences_new, [])
    # change_aligned_sequences(['-MNAHTA-FLL', 'NMFFVSANPW'], 'M-NAHT-AFL')
    # print(score_two_seq('WRYIAMRE-QYES--', '--YI-MQEVQQE--R'))
##### 13, 8, 16, 0, 7, 4
