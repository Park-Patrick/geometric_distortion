#clear
rm(list=ls())

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
  print('Usage: park_displacement_stats.R comp_type')
  print('i.e. park_displacement_stats.R A_vs_B')
  stop()
}
print(args[1])

comp_type <- args[1]

# non cmd line use
# comp_type <- 'comp_type_here'

#import libraries
library(fslr)
library(ff)
library(ffbase)
library(ggplot2)

# Reading in comma seperated file, assign to variable df (data frame)
# read in subject list
base_dir = '/home/ROBARTS/ppark/Documents/geo_nlin_only'
subjects_tsv <- '~/Documents/geo_nlin_only/subjects_temp.tsv'

############################
# project folder structure #
############################
# displacement files: base_dir/comp_type/sub-num/abc.nii.gz
# all other required files: base_dir/periphery/ (shared across all comparison types)
# files created by script (stats): base_dir/comp_type_stats/ <-- these folders must be created prior to running

df <- read.csv(subjects_tsv, header = FALSE, sep = ',' )
# print(df)

# Sets the first column name. 
names(df)[1] <- 'subjid'

# read in nlin_displacement files that are warped into standard space
# apply bet mask
atlas_path <- sprintf('%s/periphery/t1.nii.gz',base_dir)
atlas_nii <- readNIfTI2(atlas_path)


atlas_mask_path <- sprintf('%s/periphery/brainmask.nii.gz',base_dir)
mask_nii <- readNIfTI2(atlas_mask_path)

# Dynamically create file paths based on subject number from df
df$mni_disp_path <- sprintf('%s/%s/sub-%s/subj_atlas_displacement_nlin_1mm.nii.gz', base_dir, comp_type,df$subjid)

# Get number of rows in df
nSubj <- length(df[,1]);

# Calling feature of nifti file using @ method. Get 3D dimensions
vol_dim <- dim(atlas_nii@.Data)

ff_disp <- ff(0, dim=c(vol_dim, nSubj) )

for (i in 1:nSubj)
{
  current_disp_nii <- readNIfTI2(df$mni_disp_path[i])
  #current_ljac_nii <- readNIfTI2(df$mni_ljac_path[i])
  print(as.character(df$subjid[i]))
  ff_disp[,,,i] <- current_disp_nii@.Data
  #ff_ljac[,,,i] <- current_ljac_nii@.Data
}
print(0)
#########################################################################################
# compute mean/stdev for displacement and log jacobian
#########################################################################################
mean_disp_nii <- atlas_nii
mean_disp_vol <- ffapply(X=ff_disp,MARGIN=c(1,2,3),AFUN=mean,RETURN=TRUE)
mean_disp_nii@.Data <- mean_disp_vol[]
#orthographic(mean_nii)
mean_disp_nii_path <- sprintf('%s/%s_stats/OHBM_mean.nii', base_dir, comp_type)
writeNIfTI2(mean_disp_nii, mean_disp_nii_path)

stdev_disp_nii <- atlas_nii
stdev_disp_vol <- ffapply(X=ff_disp,MARGIN=c(1,2,3),AFUN=sd,RETURN=TRUE)
stdev_disp_nii@.Data <- stdev_disp_vol[]
#orthographic(mean_nii)
stdev_disp_nii_path <- sprintf('%s/%s_stats/OHBM_sd.nii', base_dir, comp_type)
writeNIfTI2(stdev_disp_nii, stdev_disp_nii_path)
  
# FIGURE: Distribution of Mean Displacements (NOTE: above zero! since at least some displacement due to registration error)
# TODO: overlay each subject's displacement -- model after Google Drive geom_histogram
# hist(mean_disp_vol[mask_nii@.Data > 0.05],breaks=200,main="Histogram of Mean Displacements within Brain Mask")
# abline(v=mean(mean_disp_vol[mask_nii@.Data > 0.05]))

#########################################################################################
# one-sided, unpaired wilcoxon (nonparametric) test looking for difference from one (1.0)
# null hypothesis: no difference from one (mu = 1.0 mm)
#########################################################################################
wilcox_disp_p_vol <- array(0, vol_dim) # set to zeros by default
print(1)
wilcox_disp_p_vol <- ffapply( AFUN=f <- function(x) {
                                        as.numeric(wilcox.test(x, mu = 0.7, alternative="greater", paired=FALSE)$p.value)
                                      }, X = ff_disp, MARGIN = c(1,2,3), RETURN = TRUE)
wilcox_disp_nii <- atlas_nii
wilcox_disp_nii@.Data <- wilcox_disp_p_vol[]
wilcox_disp_fdr <- p.adjust(wilcox_disp_p_vol[], "fdr")
wilcox_disp_fdr_vol <- array(wilcox_disp_fdr,vol_dim) # reshape into 3D volume
wilcox_disp_fdr_bin_vol <- array(0, vol_dim)
wilcox_disp_fdr_bin_vol[ mask_nii@.Data > 0.5 ] <- wilcox_disp_fdr_vol[ mask_nii@.Data > 0.5 ] < 0.05
# orthographic(wilcox_disp_fdr_bin_vol)
(sum(wilcox_disp_fdr_bin_vol) / sum( mask_nii@.Data ) )*100
# statistical significance does not mean clinical significance
print(2)

### create fdr masked displacement map
fdr_masked_nii <- atlas_nii
fdr_masked_vol <- array(0,c(182,218,182))
fdr_masked_vol[ mask_nii@.Data > 0.5 ] <- mean_disp_vol[ mask_nii@.Data > 0.5 ] * wilcox_disp_fdr_bin_vol[ mask_nii@.Data > 0.5 ]
fdr_masked_nii@.Data <- fdr_masked_vol[]
fdr_masked_nii_path <- sprintf('%s/%s_stats/OHBM_wilcoxon_masked.nii', base_dir, comp_type)
writeNIfTI2(fdr_masked_nii, fdr_masked_nii_path)

print(3)

#########################################################################################
# compute mean for displacement in an ROI
#########################################################################################
# from Harvard-Oxford ROIs thresholded at 25
ROI_HO_lobes_path <- sprintf('%s/periphery/HarvardOxford-combined-lobar-maxprob-thr25-1mm.nii.gz', base_dir)
ROI_HO_lobes_nii <- readNIfTI2(ROI_HO_lobes_path)
# orthographic(ROI_HO_lobes_nii)
num_labels <- range(ROI_HO_lobes_nii@.Data)[2] # 0...6
# csv file is here: /cluster/software/vasst-dev/pipeline/cfg/labels

df_lobes_path <- sprintf('%s/periphery/HarvardOxford-combined-lobar.csv', base_dir)
df_lobes <- read.table(df_lobes_path, header=TRUE, sep=',')

# initialize csv files
df_lobes$mean <- rep(0,num_labels)
df_lobes$stdev <- rep(0,num_labels)
df_lobes$lowCI <- rep(0,num_labels)
df_lobes$highCI <- rep(0,num_labels)
df_lobes$t_stat <- rep(0,num_labels)
df_lobes$t_pval <- rep(0,num_labels)
df_lobes$w_vval <- rep(0,num_labels)
df_lobes$w_pval <- rep(0,num_labels)
df_lobes$t_pval_adj <- rep(0,num_labels)
df_lobes$w_pval_adj <- rep(0, num_labels)

ff_label_mean <- ff(0, dim=c(num_labels,nSubj) ) # ff makes it more memory efficient (probably could be optimized further)
print(4)


for (i in 1:nSubj) { # compute mean for each subject and each label

temp_vol <- ff_disp[,,,i]
for (j in 1:num_labels) {
  #hist(temp_vol[ ROI_HO_lobes_nii@.Data == j ] # other details about data if necessary
  ff_label_mean[j,i] <- mean(temp_vol[ ROI_HO_lobes_nii@.Data == j ])
}
}
# image(ff_label_mean[]) # show heat map of different labels
for (j in 1:num_labels) {
# print(j)
t <- t.test(ff_label_mean[j,], mu = 0.7, alternative = 'greater' , paired=FALSE)
  df_lobes$mean[j] <- mean(ff_label_mean[j,])
  df_lobes$stdev[j] <- sd(ff_label_mean[j,])
  df_lobes$lowCI[j] <- t$conf.int[1]
  df_lobes$highCI[j] <- t$conf.int[2]
  df_lobes$t_stat[j] <- t$statistic
  df_lobes$t_pval[j] <- t$p.value
w <- wilcox.test(ff_label_mean[j,], mu = 0.7, alternative = 'greater', paired=FALSE)
  df_lobes$w_vval[j] <- w$statistic
  df_lobes$w_pval[j] <- w$p.value

  }

df_lobes$t_pval_adj <- p.adjust(df_lobes$t_pval, "fdr")
df_lobes$w_pval_adj <- p.adjust(df_lobes$w_pval, "fdr")
# export table to file
df_lobes_tbl_path <- sprintf('%s/%s_stats/HarvardOxford-combined-lobar-stats.csv', base_dir, comp_type)
write.table(df_lobes, df_lobes_tbl_path, sep = ",", quote = FALSE, row.names = FALSE)

###### 
print(5)

# perform again for full set of Harvard-Oxford labels (113 total)
ROI_HO_path <- sprintf('%s/periphery/HarvardOxford-combined-maxprob-thr25-1mm.nii.gz', base_dir)
ROI_HO_nii <- readNIfTI2(ROI_HO_path)
num_labels <- range(ROI_HO_nii@.Data)[2] # 0...113

df_HO_path <- sprintf('%s/periphery/HarvardOxford-combined-NetworkLabelList_fullname.csv', base_dir)
df_HO <- read.table(df_HO_path, header=TRUE, sep=',')

# initialize csv files
df_HO$mean <- rep(0,num_labels)
df_HO$stdev <- rep(0,num_labels)
df_HO$lowCI <- rep(0,num_labels)
df_HO$highCI <- rep(0,num_labels)
df_HO$t_stat <- rep(0,num_labels)
df_HO$t_pval <- rep(0,num_labels)
df_HO$w_vval <- rep(0,num_labels)
df_HO$w_pval <- rep(0,num_labels)
df_HO$t_pval_adj <- rep(0,num_labels)
df_HO$w_pval_adj <- rep(0,num_labels)

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
  # print(j)
  t <- t.test(ff_label_mean[j,], mu = 0.7, alternative="greater", paired=FALSE)
    df_HO$mean[j] <- mean(ff_label_mean[j,])
    df_HO$stdev[j] <- sd(ff_label_mean[j,])
    df_HO$lowCI[j] <- t$conf.int[1]
    df_HO$highCI[j] <- t$conf.int[2]
    df_HO$t_stat[j] <- t$statistic
    df_HO$t_pval[j] <- t$p.value
  w <- wilcox.test(ff_label_mean[j,], mu = 0.7, alternative="greater", paired=FALSE)
    df_HO$w_vval[j] <- w$statistic
    df_HO$w_pval[j] <- w$p.value
}

df_HO$t_pval_adj <- p.adjust(df_HO$t_pval, "fdr")
df_HO$w_pval_adj <- p.adjust(df_HO$w_pval, "fdr")


# export table to file
df_HO_tbl_path <- sprintf('%s/%s_stats/HarvardOxford-combined-HO-stats.csv', base_dir, comp_type)
write.table(df_HO, df_HO_tbl_path, sep = ",", quote = FALSE, row.names = FALSE)

print('stage-6')
#########

#############################################################################################################

# perform again for full set of Cerebellum labels
# /usr/share/fsl/5.0/data/atlases/Cerebellum/Cerebellum-MNIfnirt-maxprob-thr25-1mm.nii.gz
# xml label file: /usr/share/fsl/5.0/data/atlases/Cerebellum_MNIfnirt.xml
ROI_cerebellum_path <- sprintf('%s/periphery/Cerebellum-MNIfnirt-maxprob-thr25-1mm.nii.gz', base_dir)
ROI_cerebellum_nii <- readNIfTI2(ROI_cerebellum_path)
num_labels <- range(ROI_cerebellum_nii@.Data)[2] # 0...27

df_cerebellum_path <- sprintf('%s/periphery/FSLCerebellum-combined-NetworkLabelList_fullname.csv', base_dir)
df_cerebellum <- read.table(df_cerebellum_path, header=TRUE, sep=',')

# initialize csv files
df_cerebellum$mean <- rep(0,num_labels)
df_cerebellum$stdev <- rep(0,num_labels)
df_cerebellum$lowCI <- rep(0,num_labels)
df_cerebellum$highCI <- rep(0,num_labels)
df_cerebellum$t_stat <- rep(0,num_labels)
df_cerebellum$t_pval <- rep(0,num_labels)
df_cerebellum$w_vval <- rep(0,num_labels)
df_cerebellum$w_pval <- rep(0,num_labels)
df_cerebellum$t_pval_adj <- rep(0,num_labels)
df_cerebellum$w_pval_adj <- rep(0, num_labels)

ff_label_mean <- ff(0, dim=c(num_labels,nSubj) ) # ff makes it more memory efficient (probably could be optimized further)

for (i in 1:nSubj) { # compute mean for each subject and each label
  
  temp_vol <- ff_disp[,,,i]
  for (j in 1:num_labels) {
    #hist(temp_vol[ ROI_cerebellum_lobes_nii@.Data == j ] # other details about data if necessary
    ff_label_mean[j,i] <- mean(temp_vol[ ROI_cerebellum_nii@.Data == j ])
  }
}
# image(ff_label_mean[]) # show heat map of different labels
for (j in 1:num_labels) {
  # print(j)
  t <- t.test(ff_label_mean[j,], mu = 0.7, alternative="greater", paired=FALSE)
  df_cerebellum$mean[j] <- mean(ff_label_mean[j,])
  df_cerebellum$stdev[j] <- sd(ff_label_mean[j,])
  df_cerebellum$lowCI[j] <- t$conf.int[1]
  df_cerebellum$highCI[j] <- t$conf.int[2]
  df_cerebellum$t_stat[j] <- t$statistic
  df_cerebellum$t_pval[j] <- t$p.value
  w <- wilcox.test(ff_label_mean[j,], mu = 0.7, alternative="greater", paired=FALSE)
  df_cerebellum$w_vval[j] <- w$statistic
  df_cerebellum$w_pval[j] <- w$p.value
}
df_cerebellum$t_pval_adj <- p.adjust(df_cerebellum$t_pval, "fdr")
df_cerebellum$w_pval_adj <- p.adjust(df_cerebellum$w_pval, "fdr")
# export table to file
df_cerebellum_tbl_path <- sprintf('%s/%s_stats/FSL-combined-cerebellum-stats.csv', base_dir, comp_type)
write.table(df_cerebellum, df_cerebellum_tbl_path, sep = ",", quote = FALSE, row.names = FALSE)

################
print(7)

##################################################################################################################
# # Now for the whole cerebellum label


# ROI_whole_cerebellum_path <- sprintf('%s/periphery/Cerebellum-MNIfnirt-maxprob-thr25-1mm_bin.nii.gz', base_dir)
# ROI_whole_cerebellum_nii <- readNIfTI2(ROI_whole_cerebellum_path)
# num_labels <- range(ROI_whole_cerebellum_nii@.Data)[2] # 1 whole whole_cerebellum label

# label_mean <- rep(0, c(nSubj) ) # ff makes it more memory efficient (probably could be optimized further)

# for (i in 1:nSubj) { # compute mean for each subject and each label
  
#   temp_vol <- ff_disp[,,,i]
#   label_mean[i] <- mean(temp_vol[ ROI_whole_cerebellum_nii@.Data == 1 ])
# }

# t.test(label_mean, mu = 0.7, alternative="greater", paired=FALSE)
# mean(label_mean)
# t$conf.int[1]
# t$conf.int[2]
# t$statistic
# t$p.value
# w <- wilcox.test(label_mean, mu = 0.7, alternative="greater", paired=FALSE)
# w$statistic
# w$p.value
# # Not significantly different across the whole cerebellum
# # 

# df_cerebellum[df_cerebellum$w_pval < 0.05,c("name","side","mean","stdev")]
# df_HO[grep("*occipital*", df_H  # print(j)O$name),c("name","side","mean","stdev")]
