# --------------
#
# R script for extracting tails data from xls files
#
#
# 2016-09-08 GJT
#
# ---------------
#
# Based on previous scripts 
#
#
# Data file in (XLS format):
#
# - First column is a column of measured lengths
#
# - Each column has a single lane of intensity data, obtained from applying poly(A) gel
#   analysis protocol, with header (numeric) time at which taken
#
# - First row under column titles has a flag (0 or 1) depending on whether lane is to be
#   included in simulation or not: note that this is used to flag lanes which are 
#   of poor quality (flag 0) or good quality (flag 1)
#
# - Each row after this has length (in col 1) and measured intensities per lane
#
# - Note that multiple gel scans can be included in this dataset in the following way
#
#   Length        [Gel1Time1] [Gel1Time2] [Gel1Time3] [Gel2Time1] [Gel2Time2] ....
#   Flag               1           1           1           0           1      ....
#   [Gel1Length1]      X           X           X
#        ...           X           X           X
#   [Gel1LengthX]      X           X           X
#   [Gel2Length1]                                          X           X      ....
#        ...                                               X           X      ....
#   [Gel2LengthY]                                          X           X      ....
#   [Gel3Length1]
#        ...
#  
#   i.e. multiple gels can have their data added to spreadsheet (without overlap)
#   with their corresponding lengths in column 1
#   note that cells MUST be empty if not part of that gel
#
#
# Algorithm
#
# 1. Read whole XLS file containing all lanes 
#
# 2. Then per lane:
#
#   i.   Extract non-empty rows of spreadsheet
#   ii.  Fit smoothing spline through intensity 
#   iii. Evaluate spline at lengths 1-600
#
# 3. Sort lanes in time order
#
# 4. Build mRNA list variable with following components:
#    mRNA:       intensity data (intensities>0, otherwise NA)
#    times:      sorted times from spreadsheet
#    time_mask:  list of flags per time 
# 
# 5. Output RData file with same stem filename as input XLS file
#
#
#

library(gdata) # Note this uses perl to extract data from xls file
library(stats)

input_files <- dir(pattern="xls")

for(i in input_files){

   print(paste("Input file",i))
   
   input <- read.xls(xls=i,verbose=T,check.names=F,stringsAsFactors=F) # verbose=T gives debug information
   
   n_lanes <- dim(input)[2]-1
   
   mRNA <- list()
   
   mRNA$mRNA      <- array(NA,dim=c(600,n_lanes))
   mRNA$times     <- as.numeric(names(input)[-1])   # fill in times from col headers
   mRNA$time_mask <- input[2,-1]                    # fill in mask from relevant row
   
   lengths  <- input[-c(1,2), 1]
   mRNAdata <- input[-c(1,2),-1]
   
   for(j in 2:(n_lanes+1)){
   
      print(paste("Column",j,"at time",mRNA$times[j-1]))
	  
	  lanedata    <- mRNAdata[,j-1]     # extract relevant column
	  
	  na.mask     <- !is.na(lanedata)   # non-missing-data mask in lane
	  
	  lanelengths <- lengths[na.mask]   # lengths for this lane
     
      lanedata    <- lanedata[na.mask]
	  
	  # fit and evaluate spline at 1:600
	  
	  xyfunction  <- smooth.spline(x=lanelengths,y=lanedata,all.knots=T) 
	  mRNA$mRNA[,j-1] <- predict(xyfunction,1:600)$y
	  
#      windows()
      x11()
	  
	  plot(lanelengths,lanedata,type="p",col=2)
	  lines(1:600,mRNA$mRNA[,j-1])
	  
	  # remove data with lane intensity <0
	  
	  zero_mask <- mRNA$mRNA[,j-1]<0
	  mRNA$mRNA[zero_mask,j-1] <- NA
	     
   }

   time_order <- order(mRNA$times) # sort in time order
   
   mRNA$mRNA      <- mRNA$mRNA[,time_order]
   mRNA$times     <- mRNA$times[time_order]
   mRNA$time_mask <- mRNA$time_mask[time_order]
   
   out_file <- sub("xls$","RData",i)
   
   save(mRNA,file=out_file)

}


