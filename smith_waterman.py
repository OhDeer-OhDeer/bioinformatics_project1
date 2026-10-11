import globals

# ===== CLASSES =====


# ===== CLASSES =====
# Node with value and direction tracing which neighbor it got that value from.
# Takes two-item list of value and direction for initialization. This makes it
# easy to pipe a result from calc_table directly into creating a new node.
class Node:
    def __init__(self, val_and_dir):
        self.value = val_and_dir[0]
        self.direction = val_and_dir[1]

# Alignment information to be returned by align(s1, s2) function.
class Alignment:
    def __init__(self, sequences, alignments, start_indices, end_indices, score):
        self.seq1 = sequences[0]
        self.seq2 = sequences[1]
        self.seq1_aln = alignments[0]
        self.seq2_aln = alignments[1]
        self.seq1_start = start_indices[0]
        self.seq2_start = start_indices[1]
        self.seq1_end = end_indices[0]
        self.seq2_end = end_indices[1]
        self.score = score



# ===== GLOBALS =====
# Scores for match, mismatch, and gap.
match_premium = 2
mismatch_penalty = -1
indel_penalty = -2
# Sequences to align. Expects leading zero to represent empty first row/column
seq1 = ""   # Assigned only in align()
seq2 = ""   # Assigned only in align()
# Table in which nodes are stored.
table = []  # Assigned in calc_node() as called by calc_table()



# ===== FUNCTIONS, IN ORDER OF BEING CALLED BY align(s1, s2) OR __main__ =====
# Calculates the value of each node in the table.
def calc_whole_table():
    # initialize table
    global table
    table = [[Node([0, "BLANK"]) for w_base in seq2] for v_base in seq1]
    # Calculate table values
    for i in range(1, len(seq1)):
        for j in range(1, len(seq2)):
            table[i][j] = Node(calc_node(i, j))


# Calculates value of node at position [i][j] in table.
def calc_node(i, j):
    # Calculate which neighbor has highest score
    vert_score = table[i - 1][j].value + indel_penalty
    horiz_score = table[i][j-1].value + indel_penalty
    diag_score = table[i-1][j-1].value + mismatch_penalty

    # Recalc diagonal if there's a match
    if seq1[i] == seq2[j]:
        diag_score = table[i-1][j-1].value + match_premium

    # Determine maximum and return
    max_score = max(vert_score, horiz_score, diag_score, 0)
    if max_score == 0:
        return [0, "ZERO"]
    elif diag_score == max_score:
        return [diag_score, "DIAGONAL"]
    elif vert_score == max_score:
        return [vert_score, "VERTICAL"]
    else:
        return [horiz_score, "HORIZONTAL"]


# Returns information about optimal local alignment.
def find_local_aln():
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
    # local_aln[0] = aligned subsequences, local_aln[1] = start indices
    local_aln = backtrace(max_row, max_col)

    alignment = Alignment(
        [seq1, seq2],   # unaligned full sequences
        local_aln[0],   # aligned subsequences
        local_aln[1],   # start indices
        [max_row, max_col], # end indices
        max_score # score
    )
    return alignment


# Recursively traces back optimal alignment.
# Return values: [0] -> [seq1 local alignment, seq2 local alignment]
#                [1] -> [Start index in seq1, start index in seq2]
def backtrace(i, j):
    curr_node = table[i][j]
    # Recurse, then add appropriate bases to the return value
    if curr_node.direction == "VERTICAL":
        next_recur = backtrace(i-1, j)
        return [ [next_recur[0][0] + seq1[i], next_recur[0][1] + " "],
                next_recur[1] ]
    elif curr_node.direction == "HORIZONTAL":
        next_recur = backtrace(i, j - 1)
        return [ [next_recur[0][0] + " ", next_recur[0][1] + seq2[j]],
                next_recur[1] ]
    elif curr_node.direction == "DIAGONAL":
        next_recur = backtrace(i-1, j - 1)
        return [ [next_recur[0][0] + seq1[i], next_recur[0][1] + seq2[j]],
                next_recur[1] ]
    # Base case: reached start of local alignment
    else:
        #       [0]       [1]
        return [["", ""], [i, j]]


# Prints information about the alignment `a`.
def print_alignment(a):
    # Formatting to align printed sequences
    seq1_padding = " " * max(a.seq2_start - a.seq1_start, 0)
    seq2_padding = " " * max(a.seq1_start - a.seq2_start, 0)

    # Print full alignment
    print("SCORE = " + str(a.score))
    print("Sequence 1: " + seq1_padding + a.seq1[1:a.seq1_start + 1] + " --> "
          + a.seq1_aln + " <-- " + a.seq1[a.seq1_end + 1:])
    print("Sequence 2: " + seq2_padding + a.seq2[1:a.seq2_start + 1] + " --> "
          + a.seq2_aln + " <-- " + a.seq2[a.seq2_end + 1:])



# ===== PROGRAM RUNNERS =====
# Returns an Alignment object for two parameterized sequences.
# This is how this file is meant to be called by external programs.
def align(s1, s2):
    global seq1, seq2, table
    seq1 = "0" + s1
    seq2 = "0" + s2
    calc_whole_table()
    return find_local_aln()

# Print optimal local alignment and score of two manually defined sequences.
# Allows this file to be run on its own for testing purposes.
if __name__ == '__main__':
    alignment = align("HEAGAWGHEE", "PAWHEAE")
    print_alignment(alignment)



