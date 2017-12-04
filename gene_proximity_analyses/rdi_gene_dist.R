# Title     : TODO
# Objective : TODO
# Created by: henryjuho
# Created on: 28/11/2017

library(ggplot2)

fiftykb = read.delim('distance_bin_beds/bin_summaries.txt')
fiftykb$window_size = 50000
fiftykb$window_id = '50kb'

fivekb = read.delim('distance_bin_beds_5kb/bin_summaries_5kb.txt')
fivekb$window_size = 5000
fivekb$window_id = '5kb'

onekb = read.delim('distance_bin_beds_1kb/bin_summaries_1kb.txt')
onekb$window_size = 1000
onekb$window_id = '1kb'

fivekb_nc = read.delim('distance_bin_beds_5kb_noncoding/bin_summaries_5kb_nc.txt')
fivekb_nc$window_size = 5000
fivekb_nc$window_id = '5kb nc'

fivekb_nc_ucne = read.delim('distance_bin_beds_5kb_noncoding_UCNE/bin_summaries_5kb_nc_UCNE.txt')
fivekb_nc_ucne$window_size = 5000
fivekb_nc_ucne$window_id = '5kb nc UCNE'

all_data = rbind(fiftykb, fivekb, onekb, fivekb_nc, fivekb_nc_ucne)

all_data

ins_data = subset(all_data, type=='ins', select=c(bin, count, window_size, window_id))
colnames(ins_data) = c('bin', 'ins_count', 'window_size', 'window_id')

del_data = subset(all_data, type=='del', select=c(bin, count, window_size, window_id))
colnames(del_data) = c('bin', 'del_count', 'window_size', 'window_id')

dist_rdi = data.frame(cbind(ins_data, del_data$del_count))
colnames(dist_rdi) = c('bin', 'ins_count', 'window_size', 'window_id', 'del_count')
dist_rdi = subset(dist_rdi, ins_count > 10 & del_count > 10)

dist_rdi$rdi = as.numeric(as.character(dist_rdi$del_count)) / as.numeric(as.character(dist_rdi$ins_count))
dist_rdi$distance = as.numeric(as.character(dist_rdi$bin)) * dist_rdi$window_size

fifty_test_d = subset(dist_rdi, window_size==50000)
fifty_test = cor.test(fifty_test_d$rdi, fifty_test_d$distance, method='spearman', exact=NULL)
fifty_test

five_test_d = subset(dist_rdi, window_size==5000 & window_id=='5kb')
five_test = cor.test(five_test_d$rdi, five_test_d$distance, method='spearman', exact=NULL)


one_test_d = subset(dist_rdi, window_size==1000)
one_test = cor.test(one_test_d$rdi, one_test_d$distance, method='spearman', exact=NULL)
one_test

five_test_d_nc = subset(dist_rdi, window_size==5000 & window_id=='5kb nc')
five_test_nc = cor.test(five_test_d_nc$rdi, five_test_d_nc$distance, method='spearman', exact=NULL)

five_test_d_nc_ucne = subset(dist_rdi, window_size==5000 & window_id=='5kb nc UCNE')
five_test_nc_ucne = cor.test(five_test_d_nc_ucne$rdi, five_test_d_nc_ucne$distance, method='spearman', exact=NULL)

five_test
five_test_nc
five_test_nc_ucne

rdi_plot = ggplot(dist_rdi, aes(x=distance, y=rdi)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      geom_smooth(method='lm') +
      ylab('rDI')  + xlab('Distance from gene') +
      facet_wrap(~window_id, scale='free_y', nrow=2)

png(file='distance_rdi_counts.png', width=6, height=6, units='in', res=360)

rdi_plot

dev.off()