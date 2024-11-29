using POMP
using RCall

P = rmca()

R"""
library(tidyverse)
$(melt(P)) |>
   ggplot(aes(x=N,y=P))+
   geom_path()+
   scale_x_log10()+
   scale_y_log10()+
   theme_bw()
"""

R"""ggsave(filename="rmca-01.png",width=7,height=7)"""

R"""
$(melt(P)) |>
   select(time,N,P) |>
   pivot_longer(-time) |>
   ggplot(aes(x=time,y=value))+
   geom_line()+
   facet_grid(name~.,scales="free_y")+
   theme_bw()
"""

R"""ggsave(filename="rmca-02.png",width=7,height=4)"""

pfilter(P,Np=1000).logLik
