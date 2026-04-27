import DataFrames: DataFrame

"""
    parus_data

*Parus major* (Great Tit) census (all individuals), Wytham Wood, Oxfordshire.
Global Population Dynamics Database dataset #10163.
(NERC Centre for Population Biology, Imperial College (2010)
The Global Population Dynamics Database Version 2.
http://www.sw.ic.ac.uk/cpb/cpb/gpdd.html).

Original source:
McCleery, R. & Perrins, C. (1991)
Effects of predation on the numbers of Great Tits, *Parus major*.
In: Bird Population Studies,
edited by Perrins, C.M., Lebreton, J.-D. & Hirons, G.J.M.
Oxford. Univ. Press. pp. 129--147.
"""
const parus_data = DataFrame(
    year = [1960, 1961, 1962, 1963, 1964, 1965, 1966, 1967, 1968, 1969,
            1970, 1971, 1972, 1973, 1974, 1975, 1976, 1977, 1978, 1979,
            1980, 1981, 1982, 1983, 1984, 1985, 1986],
    pop = [148, 258, 185, 170, 267, 239, 196, 132, 167, 186,
           128, 227, 174, 177, 137, 172, 119, 226, 166, 161,
           199, 306, 206, 350, 214, 175, 211]
)
