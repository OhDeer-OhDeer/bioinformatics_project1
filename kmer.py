import globals

# break every database protein into overlapping k-mers and construct mapping
# key = kmer, value = set of proteins that contain that kmer
def build_kmer_index(seq_db, k):
    d = {"": []}
    for seq in seq_db:
        # For each kmer in seq, add seq to that kmer's dict entry
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer not in d:
                d[kmer] = []
            d[kmer].append(seq)
    return d

# Break the query into k-mers, count shared k-mers for each database protein,
# and retain the top N candidate proteins
def retrieve_candidates(query, index, k, N):
    # For each unique kmer in query, add seq to list
    query_kmer_list = []
    for i in range(0, len(query) - k):
        kmer = query[i:i + k]
        # Unsure if we should double-count k-mers that show up multiple times;
        # for now, it doesn't.
        if kmer not in query_kmer_list:
            query_kmer_list.append(kmer)


     # Count shared kmers for each database protein
    scores = {}

    for kmer in query_kmer_list:

        if kmer in index:
            for seq in index[kmer]:

                if seq not in scores:
                    scores[seq] = 0

                scores[seq] += 1

    # produce sorted list of tuples: (sequence, score) where score is num of shared kmers
    sorted_candidates = sorted(
        scores.items(), # get dict items as tuples
        key=lambda x: x[1], # sort by score (second item)
        reverse=True # highest to lowest
    )

    # select top N
    candidates = []
    for seq, score in sorted_candidates[:N]:
        candidates.append(seq)

    return candidates