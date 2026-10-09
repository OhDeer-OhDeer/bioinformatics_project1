# ===== CLASSES =====
# Node with value and direction tracing which neighbor it got that value from.
# Takes two-item list of value and direction for initialization. This makes it
# easy to pipe a result from calc_table directly into creating a new node.
class Node:
    def __init__(self, val_and_dir):
        self.value = val_and_dir[0]
        self.direction = val_and_dir[1]



# ===== GLOBALS =====
# Scores for match, mismatch, and gap.
match_premium = 2
mismatch_penalty = -1
indel_penalty = -2
# Sequences to align. Leading zero represents empty first row/col of table.
seq1 = "0" + "HEAGAWGHEE" # v
seq2 = "0" + "LLLVLLLAAA" # w
"""
Sample values for testing:
Query:
    HEAGAWGHEE
Testing vals:
    PAWHEAE
    HEAGAWGHEE  <- highest-ranked: identical
    HEAGAWGHDE
    LLLVLLLAAA  <- lowest-ranked: almost nothing in common
    HEAGTWGHEE
"""
# Table containing value and direction info for each base pair alignment.
table = [[Node([0, "UNINITIALIZED"]) for w_base in seq2] for v_base in seq1]



# ===== FUNCTIONS =====
# Calculates value of node at position [i][j] in table.
def calc_table(i, j):
    # Calculate which neighbor has highest score
    vert_score = table[i - 1][j].value + indel_penalty
    horiz_score = table[i][j-1].value + indel_penalty
    diag_score = table[i-1][j-1].value + mismatch_penalty
    # Recalc diagonal if there's a match
    if seq1[i] == seq2[j]:
        diag_score = table[i-1][j-1].value + match_premium
    # Determine maximum and return.
    max_score = max(vert_score, horiz_score, diag_score, 0)
    if max_score == 0:
        return [0, "ZERO"]
    elif diag_score == max_score:
        return [diag_score, "DIAGONAL"]
    elif vert_score == max_score:
        return [vert_score, "VERTICAL"]
    else:
        return [horiz_score, "HORIZONTAL"]


# Recursively traces back optimal alignment.
# Return values: [0] -> seq1 local alignment, [1] -> seq2 local alignment,
#                [2] -> [start index in seq1, start index in seq2]
def backtrace(i, j):
    curr_node = table[i][j]
    # Recurse, then add appropriate bases to return value
    if curr_node.direction == "VERTICAL":
        next_recur = backtrace(i-1, j)
        return [next_recur[0] + seq1[i], next_recur[1] + " ", next_recur[2]]
    elif curr_node.direction == "HORIZONTAL":
        next_recur = backtrace(i, j - 1)
        return [next_recur[0] + " ", next_recur[1] + seq2[j], next_recur[2]]
    elif curr_node.direction == "DIAGONAL":
        next_recur = backtrace(i-1, j - 1)
        return [next_recur[0] + seq1[i], next_recur[1] + seq2[j],
                next_recur[2]]
    # Base case: reached start of local alignment
    else:
        return ["", "", [i, j]]


# Align given sequences
if __name__ == '__main__':
    # 1. CALCULATE TABLE
    # Populate edges table with zeros (redundant?)
    for j in range(0, len(seq2)):
        table[0][j] = Node([0, ""])
    for i in range(0, len(seq1)):
        table[i][0] = Node([0, ""])

    # Calculate table values
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            table[i][j] = Node(calc_table(i, j))

    # 2. BACKTRACE
    # Find the largest table value
    max_score = -1
    max_row = -1
    max_col = -1
    for i in range(0, len(seq1)):
        for j in range(0, len(seq2)):
            if table[i][j].value > max_score:
                max_score = table[i][j].value
                max_row = i
                max_col = j

    # Recurse back from the largest value to calculate local alignment.
    # See backtrace() function for how local_aln is formatted.
    local_aln = backtrace(max_row, max_col)

    # 3. PRINT
    # Recover start/end points of local alignment
    # Start vals are minimum 1 because index zero is row/col of 0's
    seq1_start = max(local_aln[2][0], 1)
    seq2_start = max(local_aln[2][1], 1)
    seq1_end = max_row + 1
    seq2_end = max_col + 1

    # Formatting to align printed sequences
    seq1_padding = " " * max(seq2_start - seq1_start, 0)
    seq2_padding = " " * max(seq1_start - seq2_start, 0)

    # Print full alignment
    print("SCORE = " + str(max_score))
    print("Sequence 1: " + seq1_padding + seq1[1:seq1_start] + " --> "
          + local_aln[0] + " <-- " + seq1[seq1_end:])
    print("Sequence 2: " + seq2_padding + seq2[1:seq2_start] + " --> "
          + local_aln[1] + " <-- " + seq2[seq2_end:])

