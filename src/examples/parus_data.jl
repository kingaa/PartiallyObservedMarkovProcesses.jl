import DataFrames: DataFrame
import CSV: File

const parus_data = DataFrame(
    File(
        IOBuffer("""
## Parus major (Great Tit) census (all individuals)
## Wytham Wood, Oxfordshire
## Global Population Dynamics Database dataset #10163.
## (NERC Centre for Population Biology, Imperial College (2010)
## The Global Population Dynamics Database Version 2.
## http://www.sw.ic.ac.uk/cpb/cpb/gpdd.html).
##
## Original source:
## McCleery, R. & Perrins, C. (1991)
## Effects of predation on the numbers of Great Tits, Parus major.
## In: Bird Population Studies,
## edited by Perrins, C.M., Lebreton, J.-D. & Hirons, G.J.M.
## Oxford. Univ. Press. pp. 129--147.
##
year;pop
1960;148
1961;258
1962;185
1963;170
1964;267
1965;239
1966;196
1967;132
1968;167
1969;186
1970;128
1971;227
1972;174
1973;177
1974;137
1975;172
1976;119
1977;226
1978;166
1979;161
1980;199
1981;306
1982;206
1983;350
1984;214
1985;175
1986;211
"""),
        comment="#",delim=";"
    )
)
