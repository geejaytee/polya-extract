# polya-extract

Poly(A)-tail gel scanned data to tail distribution statistics

This repository holds code for automated processing of poly(A)-tail gel scanned data to tail distribution statistics for use with poly(A)-tail distributional modelling and fitting

This code converts previously scanned data (in XLS format) into RData files. For format of input XLS file see R script header.

Algorithm is as follows:

1. Read whole XLS file containing all lanes 
2. Then per lane:
   - Extract non-empty rows of spreadsheet
   - Fit smoothing spline through intensity 
   - Evaluate spline at lengths 1-600
3. Sort lanes in time order
4. Build mRNA list variable with following components:
   - mRNA:       intensity data (intensities>0, otherwise NA)
   - times:      sorted times from spreadsheet 
   - time_mask:  list of flags per time 
 
5. Output RData file with same stem filename as input XLS file

These output files are used by the modelling and parameter fitting codes for obtaining parameters for the poly(A) model.

A test file (test.xls) in the correct format is supplied for code test purposes.

Notes:
- code runs in base R, but uses gdata package to extract Excel data. gdata itself requires a working Perl installation. For Windows, Strawberry Perl works perfectly; for Unix/Linux versions there is a native Perl that should be found if installed correctly when gdata is added to R.
