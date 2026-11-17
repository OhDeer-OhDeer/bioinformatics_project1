# bioinformatics_project1

The aim of this project is to code a protein-search system and rank the results with the Smith-Waterman Algorithm.

Two methods of search:
- Exhaustive search: search entire protein database to find query
- Accelerated search: break the query into k-mers, count shared k-mers for each database protein, and retain
the top N candidate proteins. Run Smith-Waterman only on the candidates

## Project structure
```
bioinformatics_project1/
│
├── accelerated_search.py  # run accelerated search
├── exhaustive_search.py   # run exhaustive search
├── globals.py             # contains protein database
├── kmer.py                # methods to break query into k-mers and count shared k-mers
├── project_runner.py      # file to run all the experiments
└── smith_waterman.py      # Smith-Waterman algorithm code
```

## How to run this project
1. Ensure `globals.py` contains the dataset you want to search
2. run `project_runner.py`. It will produce figures and printed results to the console
