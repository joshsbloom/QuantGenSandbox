# You have two strains of yeast, the lab reference strain BY and a vineyard isolate RM 
# Is there a difference in how the yeast grow in YNB + glucose media?

# 1)-------------------------------------------------
# 120uL cultures strains were grown overnight to saturation in a 96-well plate called 'Plate080525_YNB_GLU'
# where BY is in well A01 and RM is in A02
#
# Use the Biomek robot to inccoluate the two yeast into another 96-well plate called 'Plate080625perm_YNB_GLU' to obtain replicate growth curves and address the question above.
# the 96-well plate format has rows A-H and columns 1-12, all numbers less than 10 need to have a preceding 0, for example B03, C09, etc for the Biomek to understand
#
# Make a file for the Biomek robot to innoculate 1 uL of saturated culture into 30 different random wells in a target plate for each of the two strains. 
# If you grow up this plate you will have 30 replicate growth curves for each strain.
# The edges of the plate should contain media blanks. This will evaluate whether your media was contaminated. 
# They are also more greatly affected by evaporation and yield unreliable growth curves.
# The Biomek requires a comma-delimited csv file as input with the following column heading, it should also contain DOS new line characters:
# Source.Plate, Source.Wells, Strain, Target.Plate, Target.Wells, Volume

# hints use an R data frame, set a seed (set it to 10),  use paste0(), rep(), toupper(), and sprintf() to construct the well names, and sample() to randomize them
# use the write.table() function, use the appropriate arguments to write.table() for col.names, row.names, sep, quote, and eol 
# for extra credit write the whole code in 4 lines 


# 2) --------------------------------------------------
# You've completed the growth curve experiment. Is there a difference in how the yeast grow?
# read in the resulting data 'growthCurves_01.csv' into R, use the file 'randomize_01.csv' to map the wells to the different strains
# 1. Use the growth.gcFitSpline() function from the QurvE R package to fit a curve to data for each well. Use functions in the lubridate R package, to convert time into hours with decimal units (e.g. 1.52 hours etc)
# Extract the doubling time during log phase growth. Visualize some of the fits.
# 2. Visualize the doubling times as a histogram with different colors for the two strains. Use ggplot2 to make this visualization.
# 3. Use R's t.test() function to test whether the distributions whether the means of the two groups come from the same distribution. Should you use a paired t-test?
# 4. You show your results to Josh and he's concerned about outliers, your data not being normally distributed, and not having equal variances. He suggest you perform a permutation test.
# What assumptions does the permutation test make?
# Use R's t.test() function to generate a null distribution of t statistics permuting the assignment of strain to doubling time 1000 times. Calculate an approximate p-value given this empirical null
# distribution of test statistics and your observed test statistic from your experiment. 
# 5. Given your results Josh asks you do to a power calculation. He wants to know how small of a fitness difference in doubling time we are powered to see between BY and other strains.
# He wants to know what the power would be given 3,6,12,24,48,96,and 384 replicate growth curves for each strain. Test increments of 3 minutes from 60 minutes to 240 minutes.
# Assume BY has a true mean doubling time of 90 minutes with an SD of 18 minutes (compare that to what you observed experimentally), and that the other strains also have the same SD.
# Calculate power as the fraction of 1000 tests with p<.05 for each combination of parameters. Use ggplot to plot the results. On the x-axis show the difference in doubling times 
# between BY and the hypothetical strain that you are powered to see in minutes. IF you want to have power to detect a ten minute difference in doubling time between BY and another strain
# how many replicates of each strain should you use.
# For visualization, hint, convert sample size to a factor to aid visualization.

