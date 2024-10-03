## To-do List

- plug-ins for *rprocess*
- accumulator variables
- iterated filtering
- weighted particle filter
- multiparameter particle filter? should be trivial from existing
- panel pfilter? may be trivial from existing (would need independent resampling)
- `DoubleFloats.Double64` for 128-bit likelihood computations
- documentation
- explore autodiff capabilities
- explore parallelization capabilities
- `logdmeasure` ~~instead of~~ `dmeasure` ~~with~~ `give_log=true`
- ~~speed tests for simple cases~~
- ~~separate simulations into new~~ `SimPompObject` ~~type~~
- ~~implement~~ `pfilter` ~~with new~~ `PfilterdPompObject` ~~type~~
- ~~separate parameters from~~ `PompObject`
- ~~rearrange so~~ `val_array` ~~arrays are sims x params x times~~
- ~~remove most sanity checks from workhorse functions: beef up documentation accordingly~~
- `melt` ~~methods will need to be adjusted~~
