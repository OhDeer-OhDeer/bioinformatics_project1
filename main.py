match_premium = 1
mismatch_penalty = -1
indel_penalty = -1

def calc_table(table, i, j):
    vert_score = table[i - 1][j] + indel_penalty
    horiz_score = table[i][j-1] + indel_penalty
    diag_score = table[i-1][j-1] + mismatch_penalty
    if seq1[i] == seq2[j]:
        diag_score = table[i-1][j-1] + match_premium
    return max(vert_score, horiz_score, diag_score)
    return max(vert_score, horiz_score)

def print_table(table, seq1, seq2):
    # print table
    # top labels
    print("           ", end='')  # align top labels to table
    for i in range(0, len(seq2)):
        print(seq2[i], end='        ')  # top labels
    print()  # newline

    # side labels and content
    for i in range(0, len(seq1)):
        print(seq1[i], end='  ')  # side labels

        # horiz arrows
        for j in range(0, len(seq2)):
            print("  ---> ", end='')
            if table[i][j] >= 0:
                print(" ", end='')

            print(str(table[i][j]), end='')  # current place
        print()  # newline

        # vertical and diagonal arrows
        print("   ", end='')
        for j in range(0, len(seq2)):
            if (i < len(seq1) - 1 #and j < len(seq2) - 1
                    and seq1[i + 1] == seq2[j]):
                print(" \\", end='')        # diagonal arrow 1
            else:
                print("  ", end='')
            print("      |", end='')    # vertical arrow 1
        print() # newline

        print("   ", end='')
        for j in range(0, len(seq2)):
            if (i < len(seq1) - 1 #and j < len(seq2)
                    and seq1[i + 1] == seq2[j]):
                print("    \\", end='')        # diagonal arrow 2
            else:
                print("     ", end='')
            print("   |", end='')    # vertical arrow 2
        print()  # newline

        print("   ", end='')
        for j in range(0, len(seq2)):
            print("       \\/", end='') # vertical/diagonal arrow 3
        print() # newline

# Align given sequences
if __name__ == '__main__':
    # leading zero represents empty first row/col of table
    seq1 = "0" + "TACGGGTAT" # v
    seq2 = "0" + "GGACGTACG" # w

    table = [[0 for w_base in seq2] for v_base in seq1]

    for j in range(0, len(seq2)):
        table[0][j] = j * indel_penalty
    for i in range(0, len(seq1)):
        table[i][0] = i * indel_penalty

    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            table[i][j] = calc_table(table, i, j)

    print_table(table, seq1, seq2)
