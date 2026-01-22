## To-do List

- remove all need for **RCall**
- iterated filtering
- weighted particle filter
- multiparameter particle filter? should be trivial from existing
- potentially move `accumvars` to the rprocess plugin (closer to where it is used)
- panel pfilter? may be trivial from existing (would need independent resampling)
  - is a panel pomp simply a vector of pomps with a parameter embedding?
  - or a hierarchical structure (better)
- `DoubleFloats.Double64` for 128-bit likelihood computations
- embeddings
- explore autodiff capabilities
- explore parallelization capabilities
