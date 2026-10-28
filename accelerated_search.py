import globals
import smith_waterman
import kmer
import argparse
from timeit import default_timer as timer

# retreive candidates with k-mers and run smith waterman on candidates
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--k", type=int, default=3,)
    parser.add_argument("--N", type=int, default=10)

    args = parser.parse_args()

    k = args.k
    N = args.N

    start = timer()

    index = kmer.build_kmer_index(globals.seq_db, k)

    candidate_ids = kmer.retrieve_candidates(globals.query, index, k, N)

    aln_list = []
    for s in candidate_ids:
        aln_list.append(smith_waterman.align(globals.query, s))
    
    aln_list.sort(key = lambda x: x.score, reverse=True)

    for a in aln_list:
        print("\nalignment of " + a.seq1[1:] + " and " + a.seq2[1:])
        smith_waterman.print_alignment(a)

    end = timer()
    print("time is " + str(end - start))


if __name__ == '__main__':
    main()