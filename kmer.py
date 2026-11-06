import globals

# break every database protein into overlapping k-mers and construct mapping
# key = kmer, value = set of proteins that contain that kmer
def build_kmer_index(seq_db, k):
    d = {"": []}
    for seq in seq_db:
        # For each kmer in seq, add seq to that kmer's dict entry
        for i in range(0, len(seq) - k):
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
    # TODO: count shared k-mers for each database protein
    for kmer in index:
        pass
    # TODO: retain the top N candidate proteins

    return ""