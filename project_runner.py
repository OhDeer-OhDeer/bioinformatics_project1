# where we perform the experiments and make plots
import exhaustive_search
import accelerated_search
import globals


print("Running exhaustive search")
exhaust_run = exhaustive_search.run()


print("")
print("Running accelerated search on k = 2, 3, 4")

runtimes = []

for k in range(2, 6):
    # set N, subject to change
    run = accelerated_search.run(k, 10)
    runtimes.append(run)



#TODO: find optimal k or set a k to explore effect of varying N
opt_k = 0

print("")
print("Accelerated search with varying N")
max_N = len(globals.seq_db)

for N in range(1, max_N+1):
        accelerated_search.run(opt_k, N)


#TODO: required figures
# Search runtime versus k.
# Number of candidate sequences versus k.
# Retrieval performance versus k.
# Retrieval performance versus candidate-set size N.
# At least one example of a final Smith-Waterman alignment for a retrieved protein