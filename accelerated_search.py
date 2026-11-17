import globals
import smith_waterman
import kmer
from timeit import default_timer as timer

# retrieve candidates with k-mers and run smith waterman on candidates
def run(k, N):

    start = timer()

    index = kmer.build_kmer_index(globals.seq_db, k)

    for key in index:
        print(key + ": " + str(index[key]))

    # candidate_ids = kmer.retrieve_candidates(globals.query, index, k, N)

    # aln_list = []
    # for s in candidate_ids:
    #     aln_list.append(smith_waterman.align(globals.query, s))
    
    # aln_list.sort(key = lambda x: x.score, reverse=True)

    # for a in aln_list:
    #     print("\nalignment of " + a.seq1[1:] + " and " + a.seq2[1:])
    #     smith_waterman.print_alignment(a)

    end = timer()
    print("time is " + str(end - start))

    return end - start
