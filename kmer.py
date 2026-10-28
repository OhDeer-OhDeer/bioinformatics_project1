import globals

# break every database protein into overlapping k-mers and construct mapping
# key = kmer, value = set of proteins that contain that kmer
def build_kmer_index(seq_db, k):
    return ""

# Break the query into k-mers, count shared k-mers for each database protein, and retain the top N candidate proteins
def retrieve_candidates(query, index, k, N):
    return ""