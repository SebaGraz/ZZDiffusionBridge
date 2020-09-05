library(ggplot2)
library(tidyverse)
theme_set(theme_light())

setwd("~/Sync/DOCUMENTS/onderzoek/phd-studenten/sebastiano/comparison_zz_mala")

library(readr)
benchmark <- read_csv("Sync/DOCUMENTS/onderzoek/phd-studenten/sebastiano/comparison_zz_mala/benchmark2.csv") %>% 
  mutate(yr = y/runtime) %>% mutate(sampler=fct_relevel(sampler, "MALA", "BPS","ZZ_var", "ZZ")) #%>% mutate(alpha=as.character(alpha))
head(benchmark)

#%>%filter(alpha==0.1) 



p <- benchmark %>% #filter(nbatches==50) %>% 
   ggplot(aes(x=sampler, y =yr)) +
   geom_bar(stat='identity')  + facet_grid(alpha ~ stat) + ylab("ESS per second")+
   theme(strip.text.y = element_text(angle = -90)) + coord_flip()

p_free <- benchmark %>% filter(stat %in% c('ess_min','ess_x_T2')) %>% 
  ggplot(aes(x=sampler, y =yr)) +
  geom_bar(stat='identity')  + facet_grid(stat ~ alpha) + ylab("ESS per second")+
  theme(strip.text.y = element_text(angle = -90)) + coord_flip()

p_free

p1 <- benchmark %>% filter(alpha==0.1) %>% 
  ggplot(aes(x=sampler, y =yr)) +
  geom_bar(stat='identity')  + facet_grid(. ~ stat, scales='free') + ylab("ESS per second")+
  theme(strip.text.y = element_text(angle = -90)) + coord_flip()
p1
 
p2 <- benchmark %>% filter(alpha==0.3) %>% 
  ggplot(aes(x=sampler, y =yr)) +
  geom_bar(stat='identity')  + facet_grid(. ~ stat, scales='free') + ylab("ESS per second")+
  theme(strip.text.y = element_text(angle = -90)) + coord_flip()
p2


p3 <- benchmark %>% filter(alpha==0.5) %>% 
  ggplot(aes(x=sampler, y =yr)) +
  geom_bar(stat='identity')  + facet_grid(. ~ stat, scales='free') + ylab("ESS per second")+
  theme(strip.text.y = element_text(angle = -90)) + coord_flip()
p3



pdf("comp_zz_mala_bps.pdf",width=7,height=5)
show(p)
dev.off()

pdf("comp_zz_mala_bps_free.pdf",width=7,height=5)
show(p_free)
dev.off()

pdf("comp_zz_mala_bps_free_alphasmall.pdf",width=7,height=5)
show(p1)
dev.off()

pdf("comp_zz_mala_bps_free_alphamedium.pdf",width=7,height=5)
show(p2)
dev.off()

pdf("comp_zz_mala_bps_free_alphalarge.pdf",width=7,height=5)
show(p3)
dev.off()