library(ggplot2)
library(gridExtra)
library(dplyr)

summary_data = read.delim('bgi10_stats.txt', na=0)
call_data = read.delim('bgi10_call.txt')

all_data = cbind(summary_data, call_data)

all_data[7] = NULL
all_data[7] = NULL

all_data$tw_per_site = all_data$theta_w / all_data$call
all_data$pi_per_site = all_data$pi / all_data$call

write.csv(all_data, 'bgi10_summary_stats.csv', row.names=FALSE)

all_data = subset(all_data, category != 'noncoding_noUCNEs'
  & category != 'cds_frameshift' & category != 'cds_non_frameshift'
  & category != 'noncoding')

all_data$category = factor(all_data$category,
  levels=rev(c('CDS', '0fold', 'nonsense', 'UCNE', '4fold', 'ALL', 'introns', 'intergenic', 'AR')))

# reset 0fold 4fold nonsense indel tajd to 0

# 0fold
all_data[2,]$tajD = NA
all_data[3,]$tajD = NA
all_data[4,]$tajD = NA
all_data[2,]$pi_per_site = NA
all_data[3,]$pi_per_site = NA
all_data[4,]$pi_per_site = NA
all_data[2,]$tw_per_site = NA
all_data[3,]$tw_per_site = NA
all_data[4,]$tw_per_site = NA

# 4fold
all_data[6,]$tajD = NA
all_data[7,]$tajD = NA
all_data[8,]$tajD = NA
all_data[6,]$pi_per_site = NA
all_data[7,]$pi_per_site = NA
all_data[8,]$pi_per_site = NA
all_data[6,]$tw_per_site = NA
all_data[7,]$tw_per_site = NA
all_data[8,]$tw_per_site = NA

# theta_w
tw = ggplot(all_data, aes(x=category, y=tw_per_site, colour=variation))+
  geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab(expression(theta[w])) +
  theme(legend.title=element_blank(), legend.position='none',
  axis.text.x=element_text(angle=45, hjust=1))

# pi
pi = ggplot(all_data, aes(x=category, y=pi_per_site, colour=variation))+
  geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab(expression(pi)) +
  theme(legend.title=element_blank(), legend.position='none',
  axis.text.x=element_text(angle=45, hjust=1))

# tajD
tajd = ggplot(all_data, aes(x=category, y=tajD, colour=variation))+
  geom_point(stat='identity', position = position_dodge(width=0.9), size = 2) +
  theme_bw() +
  xlab('')  + ylab("Tajima's D") +
  #ylim(-1.3, -0.2) +
  theme(legend.title=element_blank(), legend.position=c(0.2, 0.3),
  axis.text.x=element_text(angle=45, hjust=1))

# pdf(file='gt_summary_stats.pdf', width=9, height=3)
#
# grid.arrange(tajd, tw, pi, nrow=1)
#
# dev.off()

pdf(file='gt_tajd.pdf', width=3, height=3)

tajd

dev.off()

png(file='gt_summary_stats.png', width=9, height=3, units='in', res=360)

grid.arrange(tajd, tw, pi, nrow=1)

dev.off()