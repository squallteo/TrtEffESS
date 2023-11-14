rm(list=ls())
source("logOR.R")
source("logOR2.R")


ttt <- inner_join(out, out2, by = c("Correction", "sizeIU")) %>%
  mutate(diff=nSubj - nSubj2)


png("diff.png", width = 800, height = 400)
ggplot(ttt, aes(x=sizeIU, y=diff, group=Correction, color=Correction)) + geom_line(size = 1.2) + 
  xlab("Size of Information Unit") + ylab("Difference: ESS in Total Number of Subjects") + theme_bw() +
  theme(axis.text = element_text(size = 15),
        axis.title = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, vjust = -8, size = 15),
        legend.text = element_text(size = 15)
  )
dev.off()