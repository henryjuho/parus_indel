# Title     : TODO
# Objective : TODO
# Created by: henryjuho
# Created on: 28/11/2017

library(ggplot2)

fiftykb = read.delim('distance_bin_beds/bin_summaries.txt')
fiftykb$window_size = 50000

fivekb = read.delim('distance_bin_beds_5kb/bin_summaries_5kb.txt')
fivekb$window_size = 5000

onekb = read.delim('distance_bin_beds_1kb/bin_summaries_1kb.txt')
onekb$window_size = 1000

all_data = rbind(fiftykb, fivekb, onekb)

ins_data = subset(all_data, type=='ins', select=c(bin, count, window_size))
colnames(ins_data) = c('bin', 'ins_count', 'window_size')

del_data = subset(all_data, type=='del', select=c(bin, count, window_size))
colnames(del_data) = c('bin', 'del_count', 'window_size')

dist_rdi = data.frame(cbind(ins_data, del_data$del_count))
colnames(dist_rdi) = c('bin', 'ins_count', 'window_size', 'del_count')
dist_rdi = subset(dist_rdi, ins_count > 5 & del_count > 5)

dist_rdi$rdi = as.numeric(as.character(dist_rdi$del_count)) / as.numeric(as.character(dist_rdi$ins_count))
dist_rdi$distance = as.numeric(as.character(dist_rdi$bin)) * dist_rdi$window_size

fifty_test_d = subset(dist_rdi, window_size==50000)
fifty_test = cor.test(fifty_test_d$rdi, fifty_test_d$distance, method='spearman', exact=NULL)
fifty_test

five_test_d = subset(dist_rdi, window_size==5000)
five_test = cor.test(five_test_d$rdi, five_test_d$distance, method='spearman', exact=NULL)
five_test

one_test_d = subset(dist_rdi, window_size==1000)
one_test = cor.test(one_test_d$rdi, one_test_d$distance, method='spearman', exact=NULL)
one_test

rdi_plot = ggplot(dist_rdi, aes(x=distance, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      geom_smooth(method='lm') +
      ylab('rDI')  + xlab('Distance from gene') +
      facet_wrap(~window_size, scale='free_y')

png(file='distance_rdi_counts.png', width=9, height=3, units='in', res=360)

rdi_plot

dev.off()