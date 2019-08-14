# Read in masks and convert from 3D to 1D vector
mask_dir <- '/Users/zoe/Documents/Bruckert_fMRI/FMRI/LI_Zoe/masks/'

left_mask <- read.csv(paste0(mask_dir, 'Frontal_left00000'), sep='')
right_mask <- read.csv(paste0(mask_dir, 'Frontal_right00000'), sep='')

left_mask <- as.vector(t(left_mask))
right_mask <- as.vector(t(right_mask))

# Loop through subjects to read in data, convert to vector and mask
data_dir <- '/Users/zoe/Documents/Bruckert_fMRI/FMRI/LI_Zoe/pe1_data/'

file_list <- list.files(data_dir, pattern='*00000')

leftdata <- matrix(data=NA, nrow=length(file_list), ncol=length(which(left_mask==1)))
rightdata <- matrix(data=NA, nrow=length(file_list), ncol=length(which(right_mask==1)))

for (i in 1:length(file_list)){
  
  mydata <- read.csv(paste0(data_dir, file_list[i]), sep=' ')
  mydata <- mydata[ , -92] # For some reason a row of NAs has appeared at the end of the data since I added the percentage signal change thing. I'm ignoring it.
  mydata <- as.vector(t(mydata))
  
  leftdata[i, ] <- mydata[which(left_mask==1)]
  rightdata[i, ] <- mydata[which(right_mask==1)]

}





