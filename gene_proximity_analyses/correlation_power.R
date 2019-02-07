library(ggplot2)
library(gridExtra)
library(gtools)
library(viridis)

# getting theta anavar data estimates into shape

theta_data = read.csv('gt_nc_v_ar_2kb_wind_cds.results.csv')
window_size = 2

theta_data$bin = factor(theta_data$bin, levels=mixedsort(as.character(subset(theta_data, sel_type=='sel' & var_type=='ins')$bin)))
theta_data$distance = as.numeric(theta_data$bin) * window_size
theta_data = subset(theta_data, sel_type=='sel', select=c('theta', 'var_type', 'distance'))

theta_cors = NULL
str(theta_data)

# getting pi estimates into shape

sum_file = read.delim('gt_prox_summary_stats.txt')
call_file = read.delim('gt_prox_call.txt')

pi_data = cbind(sum_file, call_file)
pi_data[7] = NULL
pi_data[7] = NULL

pi_data$dist = as.numeric(sapply(strsplit(as.character(pi_data$category), 'bin'), "[[", 2)) * 2

pi_data$pi_per_site = pi_data$pi / pi_data$callable

pi_data = subset(pi_data, variation == 'DEL' | variation == 'INS', select=c('pi_per_site', 'variation', 'dist'))
colnames(pi_data) <- c('theta', 'var_type', 'distance')
pi_data$var_type = as.factor(tolower(pi_data$var_type))
str(pi_data)

for(data in c('theta', 'pi')){

    if(data == 'theta'){target_set=theta_data}else{target_set=pi_data}

    for(variant in c('ins', 'del')){

        var_data = subset(target_set, var_type==variant)

        for(i in 1:length(var_data$theta)-1){

            trimmed_set = var_data[i:length(var_data$theta),]

            test = cor.test(as.numeric(trimmed_set$theta), trimmed_set$distance, method='spearman', exact=NULL)
            p = test$p.value
            rho = test$estimate
            if(p > 0.05){
                sig = 'p > 0.05'
            }else{
                sig = 'p <= 0.05'}

            out_row = c(i * window_size, p, rho, variant, sig, data)
            theta_cors = rbind(theta_cors, out_row)
        }
    }
}
theta_cors = as.data.frame(theta_cors)
colnames(theta_cors) = c('min_dist', 'p', 'rho', 'var_type', 'p_cat', 'param')

pal = viridis(n=3)[1:3]

theta_plot = ggplot(theta_cors, aes(x=as.numeric(as.character(min_dist)), y=as.numeric(as.character(rho)), colour=toupper(var_type), shape=p_cat)) +
      geom_point(stat='identity', size = 2) +
      theme_bw() +
      #xlim(0,10000)+
      ylab(expression(rho))  + xlab('Distance of first bin from exon (kb)') +
      facet_wrap(~param, nrow=1, labeller=label_parsed) +
      theme(legend.title = element_blank(), legend.position = c(0.15, 0.15),
            legend.background = element_blank(), legend.box = "horizontal") +
      scale_colour_manual(values=pal)


pdf(file='theta_cor_power.pdf', width=6, height=3)

theta_plot

dev.off()
#
# #sum stats plot
# sum_file = read.delim('../gene_proximity_analyses/gt_prox_summary_stats.txt')
# call_file = read.delim('../gene_proximity_analyses/gt_prox_call.txt')
#
# all_data = cbind(sum_file, call_file)
# all_data[7] = NULL
# all_data[7] = NULL
#
# all_data$dist = as.numeric(sapply(strsplit(as.character(all_data$category), 'bin'), "[[", 2)) * 2
#
# all_data$pi_per_site = all_data$pi / all_data$callable
#
# ins_theta = subset(all_data, variation == 'INS', select=c(category, pi_per_site, dist))
# colnames(ins_theta) = c('bin', 'ins_theta', 'distance')
#
# cor.test(as.numeric(ins_theta$ins_theta), ins_theta$distance, method='spearman', exact=NULL)
#
# del_theta = subset(all_data, variation == 'DEL', select=c(category, pi_per_site, dist))
# colnames(del_theta) = c('bin', 'del_theta', 'distance')
#
# cor.test(as.numeric(del_theta$del_theta), del_theta$distance, method='spearman', exact=NULL)
#
# rdi_data = as.data.frame(cbind(ins_theta, del_theta$del_theta))
# rdi_data$rdi = rdi_data$del_theta / rdi_data$ins_theta
#
# cor.test(as.numeric(rdi_data$rdi), as.numeric(rdi_data$bin), method='spearman', exact=NULL)
#
# theta_plot = ggplot(subset(all_data, variation!='SNP' & variation!='INDEL'), aes(x=dist, y=pi_per_site, colour=variation)) +
#       geom_point(stat='identity', size = 2) +
#       theme_bw() +
#       geom_smooth(method='lm', se=FALSE)+
#       #xlim(0,10000)+
#       ylab(expression(pi))  + xlab('Distance from exons (kb)') +
#       theme(legend.title = element_blank(), legend.position = c(0.15, 0.8), legend.background = element_blank()) +
#       scale_colour_manual(values=pal)
#
# # gamma_plot = ggplot(subset(all_data, variation!='SNP' & variation!='INDEL'), aes(x=dist, y=tajD, colour=variation)) +
# #       geom_point(stat='identity', size = 2) +
# #       theme_bw() +
# #       #xlim(0,10000)+
# #       ylab("Tajima's D")  + xlab('Distance from CDS (kb)') +
# #       theme(legend.title = element_blank(), legend.position = 'none', legend.background = element_blank())
# #
# # rdi_plot = ggplot(rdi_data, aes(x=distance, y=rdi)) +
# #       geom_point(stat='identity', size = 2) +
# #       theme_bw() +
# #       #xlim(0,10000)+
# #       geom_smooth(method='lm') +
# #       ylab('rDI')  + xlab('Distance from CDS (kb)')
#
# pdf('linked_sel_pi.pdf', width=3, height=3)
#
# theta_plot
#
# dev.off()
#
# pdf('dist_gamma_plot.pdf', width=3, height=3)
#
# gamma_plot
#
# dev.off()
#
#
