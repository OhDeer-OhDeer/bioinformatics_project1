import globals
import smith_waterman
from timeit import default_timer as timer


# Run full smith-waterman search on all database sequences
if __name__ == '__main__':
    start = timer()

    aln_list = []
    for s in globals.seq_db:
        # Append        optimal alignment    of the query   against sequence s
        aln_list.append(smith_waterman.align(globals.query, s))
    aln_list.sort(key = lambda x: x.score, reverse=True)

    for a in aln_list:
        print("\nalignment of " + a.seq1[1:] + " and " + a.seq2[1:])
        smith_waterman.print_alignment(a)

    end = timer()
    print("time is " + str(end - start))