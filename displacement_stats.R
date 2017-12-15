# script for loading in processed displacement maps for each subject and computing statistical maps
# resampled to 1mm isotropic for the group analysis
# first step: one-tailed t-statistics
# second step: visualization of voxel-wise plots of displacement, and histogram
# third step: consider nonparametric testing alternatives

rm(list=ls())

#install.packages('fslr')
#install.packages('ff')
#install.packages('ffbase')
#install.packages('ggplot2')

library(fslr)
library(ff)
library(ffbase)
library(ggplot2)

# TODO
df <- read.csv('~/khangrp/users/ppark/OHBM_dbm/subjects_temp.tsv', header = FALSE, sep = ',' )
names(df)[1] <- 'subjid'

# read in nlin_displacement files that are warped into standard space
# apply bet mask
atlas_path <- '/eq-nas/jlau/EpilepsyDatabase/standard/projects/atlases/MNI152_1mm/t1/t1.nii.gz'
atlas_nii <- readNIfTI2(atlas_path)

atlas_mask_path <- '/eq-nas/jlau/EpilepsyDatabase/standard/projects/atlases/MNI152_1mm/t1/brainmask.nii.gz'
mask_nii <- readNIfTI2(atlas_mask_path)

df$mni_disp_path <- sprintf('/home/jlau/khangrp/users/ppark/OHBM_dbm/sub-%s/subj_MNI_displacement_nlin_1mm.nii.gz', df$subjid)
#df$mni_ljac_path <- sprintf('/eq-nas/jlau/EpilepsyDatabase/Projects/test_geo_dbm/%s/Preop/Processed/MetricsToMNI/subj_MNI_log_jacobian_nlin_1mm.nii.gz', df$subjid)

nSubj <- length(df[,1]);

vol_dim <- dim(atlas_nii@.Data)

ff_disp <- ff(0, dim=c(vol_dim,nSubj) ) # ff makes it more memory efficient (probably could be optimized further)

for (i in 1:nSubj)
{
  current_disp_nii <- readNIfTI2(df$mni_disp_path[i])
  #current_ljac_nii <- readNIfTI2(df$mni_ljac_path[i])
  print(as.character(df$subjid[i]))
  ff_disp[,,,i] <- current_disp_nii@.Data
  #ff_ljac[,,,i] <- current_ljac_nii@.Data
}

#########################################################################################
# compute mean/stdev for displacement and log jacobian
#########################################################################################
  mean_disp_nii <- atlas_nii
  mean_disp_vol <- ffapply(X=ff_disp,MARGIN=c(1,2,3),AFUN=mean,RETURN=TRUE)
  mean_disp_nii@.Data <- mean_disp_vol[]
  #orthographic(mean_nii)
  writeNIfTI2(mean_disp_nii, '~/khangrp/users/ppark/OHBM_dbm/stats/OHBM_mean.nii')
  
  stdev_disp_nii <- atlas_nii
  stdev_disp_vol <- ffapply(X=ff_disp,MARGIN=c(1,2,3),AFUN=sd,RETURN=TRUE)
  stdev_disp_nii@.Data <- stdev_disp_vol[]
  #orthographic(mean_nii)
  writeNIfTI2(stdev_disp_nii, '~/khangrp/users/ppark/OHBM_dbm/stats/OHBM_sd.nii')
  
  # FIGURE: Distribution of Mean Displacements (NOTE: above zero! since at least some displacement due to registration error)
  # TODO: overlay each subject's displacement -- model after Google Drive geom_histogram
  hist(mean_disp_vol[mask_nii@.Data > 0.05],breaks=200,main="Histogram of Mean Displacements within Brain Mask")
  abline(v=mean(mean_disp_vol[mask_nii@.Data > 0.05]))
  
#########################################################################################
# one-sided, unpaired wilcoxon (nonparametric) test looking for difference from one (1.0)
# null hypothesis: no difference from one (mu = 1.0 mm)
#########################################################################################
  wilcox_disp_p_vol <- array(0, vol_dim) # set to zeros by default
  wilcox_disp_p_vol <- ffapply( AFUN=f <- function(x) {
                                            as.numeric(wilcox.test(x, mu = 1.0, alternative="greater", paired=FALSE)$p.value)
                                          }, X = ff_disp, MARGIN = c(1,2,3), RETURN = TRUE)
  wilcox_disp_nii <- atlas_nii
  wilcox_disp_nii@.Data <- wilcox_disp_p_vol[]
  wilcox_disp_fdr <- p.adjust(wilcox_disp_p_vol[], "fdr")
  wilcox_disp_fdr_vol <- array(wilcox_disp_fdr,vol_dim) # reshape into 3D volume
  wilcox_disp_fdr_bin_vol <- array(0, vol_dim)
  wilcox_disp_fdr_bin_vol[ mask_nii@.Data > 0.5 ] <- wilcox_disp_fdr_vol[ mask_nii@.Data > 0.5 ] < 0.05
  orthographic(wilcox_disp_fdr_bin_vol)
  (sum(wilcox_disp_fdr_bin_vol) / sum( mask_nii@.Data ) )*100
  # statistical significance does not mean clinical significance

### create fdr masked displacement map
fdr_masked_nii <- atlas_nii;
fdr_masked_vol <- array(0,c(182,218,182))
fdr_masked_vol[ mask_nii@.Data > 0.5 ] <- mean_disp_vol[ mask_nii@.Data > 0.5 ] * wilcox_disp_fdr_bin_vol[ mask_nii@.Data > 0.5 ]
fdr_masked_nii@.Data <- fdr_masked_vol[]
writeNIfTI2(fdr_masked_nii, '~/khangrp/users/ppark/OHBM_dbm/stats/OHBM_wilcoxon_masked.nii')

#########################################################################################
# compute mean for displacement in an ROI
#########################################################################################
# from Harvard-Oxford ROIs thresholded at 25
ROI_HO_lobes_path <- '~/khangrp/projects/unsorted/test_geo_dbm/HarvardOxford-combined-lobar-maxprob-thr25-1mm.nii.gz'
ROI_HO_lobes_nii <- readNIfTI2(ROI_HO_lobes_path)
# orthographic(ROI_HO_lobes_nii)
num_labels <- range(ROI_HO_lobes_nii@.Data)[2] # 0...6
# csv file is here: /cluster/software/vasst-dev/pipeline/cfg/labels
df_lobes <- read.table('~/khangrp/projects/unsorted/test_geo_dbm/HarvardOxford-combined-lobar.csv', header=TRUE, sep=',')

# initialize csv files
df_lobes$mean <- rep(0,num_labels)
df_lobes$stdev <- rep(0,num_labels)
df_lobes$lowCI <- rep(0,num_labels)
df_lobes$highCI <- rep(0,num_labels)
df_lobes$t_stat <- rep(0,num_labels)
df_lobes$t_pval <- rep(0,num_labels)
df_lobes$w_vval <- rep(0,num_labels)
df_lobes$w_pval <- rep(0,num_labels)

ff_label_mean <- ff(0, dim=c(num_labels,nSubj) ) # ff makes it more memory efficient (probably could be optimized further)

for (i in 1:nSubj) { # compute mean for each subject and each label

  temp_vol <- ff_disp[,,,i]
  for (j in 1:num_labels) {
    #hist(temp_vol[ ROI_HO_lobes_nii@.Data == j ] # other details about data if necessary
    ff_label_mean[j,i] <- mean(temp_vol[ ROI_HO_lobes_nii@.Data == j ])
  }
}
# image(ff_label_mean[]) # show heat map of different labels
for (j in 1:num_labels) {
  print(j)
  t <- t.test(ff_label_mean[j,], mu = 1, alternative = 'greater' , paired=FALSE)
    df_lobes$mean[j] <- mean(ff_label_mean[j,])
    df_lobes$stdev[j] <- sd(ff_label_mean[j,])
    df_lobes$lowCI[j] <- t$conf.int[1]
    df_lobes$highCI[j] <- t$conf.int[2]
    df_lobes$t_stat[j] <- t$statistic
    df_lobes$t_pval[j] <- t$p.value
  w <- wilcox.test(ff_label_mean[j,], mu = 1, alternative = 'greater', paired=FALSE)
    df_lobes$w_vval[j] <- w$statistic
    df_lobes$w_pval[j] <- w$p.value
}
# export table to file
write.table(df_lobes, '/home/jlau/khangrp/users/ppark/OHBM_dbm/stats/HarvardOxford-combined-lobar-stats.csv', sep = ",", quote = FALSE, row.names = FALSE)

#############################################################################################################

# perform again for full set of Harvard-Oxford labels (113 total)
ROI_HO_path <- '~/khangrp/projects/unsorted/test_geo_dbm/HarvardOxford-combined-maxprob-thr25-1mm.nii.gz'
ROI_HO_nii <- readNIfTI2(ROI_HO_path)
num_labels <- range(ROI_HO_nii@.Data)[2] # 0...113

df_HO <- read.table('~/khangrp/projects/unsorted/test_geo_dbm//HarvardOxford-combined-NetworkLabelList_fullname.csv', header=TRUE, sep=',')

# initialize csv files
df_HO$mean <- rep(0,num_labels)
df_HO$stdev <- rep(0,num_labels)
df_HO$lowCI <- rep(0,num_labels)
df_HO$highCI <- rep(0,num_labels)
df_HO$t_stat <- rep(0,num_labels)
df_HO$t_pval <- rep(0,num_labels)
df_HO$w_vval <- rep(0,num_labels)
df_HO$w_pval <- rep(0,num_labels)

ff_label_mean <- ff(0, dim=c(num_labels,nSubj) ) # ff makes it more memory efficient (probably could be optimized further)

for (i in 1:nSubj) { # compute mean for each subject and each label

  temp_vol <- ff_disp[,,,i]
  for (j in 1:num_labels) {
    #hist(temp_vol[ ROI_HO_lobes_nii@.Data == j ] # other details about data if necessary
    ff_label_mean[j,i] <- mean(temp_vol[ ROI_HO_nii@.Data == j ])
  }
}
# image(ff_label_mean[]) # show heat map of different labels
for (j in 1:num_labels) {
  print(j)
  t <- t.test(ff_label_mean[j,], mu = 1, alternative="greater", paired=FALSE)
    df_HO$mean[j] <- mean(ff_label_mean[j,])
    df_HO$stdev[j] <- sd(ff_label_mean[j,])
    df_HO$lowCI[j] <- t$conf.int[1]
    df_HO$highCI[j] <- t$conf.int[2]
    df_HO$t_stat[j] <- t$statistic
    df_HO$t_pval[j] <- t$p.value
  w <- wilcox.test(ff_label_mean[j,], mu = 1, alternative="greater", paired=FALSE)
    df_HO$w_vval[j] <- w$statistic
    df_HO$w_pval[j] <- w$p.value
}
# export table to file
write.table(df_HO, '/home/jlau/khangrp/users/ppark/OHBM_dbm/stats/HarvardOxford-combined-HO-stats.csv', sep = ",", quote = FALSE, row.names = FALSE)
