Rmimic
==============

**<p align="center" style="padding-left: 2em; padding-right: 2em;">This package simply provides wrappers around well-established and validated packages (e.g. stats, lawstat,   
    psychometric, MBESS, lme4, psych, epitools) to make their use and output more user friendly.</p>**

This R package mimics the functionalities of popular commercial statistics software packages in that it will compute
multiple-related tests. However, rather than outputting all results blindly, the package uses the specified/
automatically determined test correction/adjustments and outputs the results in APA format. Additionally, 
when appropriate, functions will automatically break down the data to consider all potential comparisons
(A vs B, A vs C, B vs C). 

This package deviates from the typical R ethos such that "super functions" will output
results to the console window in a style that roughly mimics the output of popular statistics software
packages. Test results are also available as environmental variables and console outputs can be suppressed
using verbose calls to the function.

Several functions are also available through the R Studio Addins dropdown menu to provide a graphic user interface similar to popular commercial statistics packages. You may have to restart R Studio after installing the package for the functions to show up. Functions accessed through the dropdown menu will pop up a graphical user interface for selecting variables and parameters, when the results of the function printed to the console window along with the syntax.
<p align="center"><img src="/screencaps/screencap_PopGUI.png?raw=true" width="400" alt="screencap PopGUI"><img src="/screencaps/screencap_PopGUI2.png?raw=true" width="400" alt="screencap PopGUI2"></p>

Installation
------------
To use this package, from R run the following commands:
```r
    tryCatch(library(devtools), error=function(e){install.packages("devtools"); library(devtools)})
    devtools::install_github("mattpontifex/Rmimic", force=TRUE); library(Rmimic)  
```

**Common Issues**
* You may be prompted to update packages, you can select 'None' (3) to just install Rmimic. 
* You may get a backports error installing devtools, this may require you to restart R and then the code below.
After completing that step, you can repeat the installation instructions above.
```r
    install.packages("backports") 
```
* You may get an error indicating that R cannot remove prior installation of a package. 
You can try installing the package again to see if that fixes it. Most times you will need to run the code below, 
then manually delete the folder and then run the install packages command.
```r
    find.package('farver') 
```



Super Function List
------------
These functions mimic the overarching outputs of popular commercial statistics software packages that lump together
several related or inherently sequential tests.


* **RmimicAnova**: Function that computes a SPSS style univariate ANOVA with effect size and confidence
intervals using the ezANOVA function. Main effects and interactions are automatically decomposed using the
specified post-hoc corrections. 
```r
    anovaresult <- Rmimic::RmimicAnova(data = PlantGrowth, dependentvariable='weight',
                            subjectid=NULL, between='group', within=NULL, sphericity='Greenhouse-Geisser', 
                            posthoc='False Discovery Rate Control', verbose=TRUE)  
```
<p align="center"><img src="/screencaps/screencap_RmimicAnova.png?raw=true" width="600" alt="screencap RmimicAnova"></p>


* **RmimicLMcontrast**: Compute SPSS style results for regression analysis with effect size and confidence intervals. This function takes stats::lm fits for a base model and the model of interest and calculates statistics for the model of interest relative to the base model. This function is also able to take stats::glm binomial family model fits for logistic regression.
```r
    basefit <- lm(mpg ~ am + wt, data = mtcars)
    fit <- lm(mpg ~ am + wt + qsec, data = mtcars)
    regresult <- Rmimic::RmimicLMcontrast(basefit, fit, 
                            confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE)
```
<p align="center"><img src="/screencaps/screencap_RmimicLMcontrast1.png?raw=true" width="600" alt="screencap RmimicLMcontrast1"><img src="/screencaps/screencap_RmimicLMcontrast2.png?raw=true" width="600" alt="screencap RmimicLMcontrast2"></p>

* **RmimicTtest**: Function that computes SPSS style t-tests with effect size and confidence intervals. 
Optional parameters are also provided to compute non-parametric t-tests with appropriate non-parametric
effect size estimates. For parametric test the function automatically determines if the variances are equal
using levene's test and outputs the correct statistcs. The function can handle factors with more than 2 levels
and will perform t-tests for each comparison with post-hoc comparison corrections.
```r
    ttestresult <- Rmimic::RmimicTtest(PlantGrowth, dependentvariable='weight',  
                            subjectid=NULL, between='group', within=NULL,  
                            nonparametric=FALSE, posthoc='Holm-Bonferroni', verbose=TRUE)   
```
<p align="center"><img src="/screencaps/screencap_RmimicTtest.png?raw=true" width="600" alt="screencap RmimicTtest"></p>

* **RmimicChisquare**: Function that computes SPSS style results for Chi-square analysis with odds ratios
and confidence intervals. For samples less than 1000, Fishers exact test statistic is used if possible.
The function can handle outcomes with more than 2 levels and will perform comparisons for each pair of outcomes.
```r
    tempdata <- data.frame("Age"="8","Pet"="Dog","Freq"=282)
    tempdata <- rbind(tempdata, data.frame("Age"="30","Pet"="Dog","Freq"=199))
    tempdata <- rbind(tempdata, data.frame("Age"="8","Pet"="Cat","Freq"=170))
    tempdata <- rbind(tempdata, data.frame("Age"="30","Pet"="Cat","Freq"=240))
    chisquareresult <- RmimicChisquare( 
                            x='Age', y='Pet',  data=tempdata, 
                            posthoc='False Discovery Rate Control',
                            confidenceinterval=0.95, studywiseAlpha=0.05, 
                            planned=FALSE, verbose=TRUE)

    chisquareresult <- Rmimic::RmimicChisquare(x='Sex', y='Survived', data=Titanic,
                            posthoc='False Discovery Rate Control', planned=FALSE, 
                            confidenceinterval=0.95, studywiseAlpha=0.05, verbose=TRUE)  
```
<p align="center"><img src="/screencaps/screencap_RmimicChisquare1.png?raw=true" width="600" alt="screencap RmimicChisquare1"></p>
<p align="center"><img src="/screencaps/screencap_RmimicChisquare2.png?raw=true" width="600" alt="screencap RmimicChisquare2"></p>

* **correlation**: Function that computes SPSS style correlations or partial correlations, with optional
parameters for the approach (pearson (default), spearman, or kendall).
```r
    tempdata <- data.frame("X" = runif(100), "Y" = runif(100), "Z" = runif(100))  
    corresult <- Rmimic::correlation(variables=c('X', 'Y', 'Z'), partial=FALSE, 
                            data=tempdata, method='pearson', listwise=TRUE, studywiseAlpha=0.05, 
                            confidenceinterval=0.95, verbose=TRUE)  
```
<p align="center"><img src="/screencaps/screencap_Correlation.png?raw=true" width="600" alt="screencap Correlation"></p>

* **descriptives**: Function that computes SPSS style descriptive statistics and frequencies.
```r
    tempdata <- data.frame("Group" = sample(1:2,100, replace=TRUE), "X" = runif(100), "Y" = runif(100))
    desc <- Rmimic::descriptives(variables=c('X','Y'), groupvariable=c("Group"), data=tempdata, verbose=TRUE) 
``` 
<p align="center"><img src="/screencaps/screencap_descriptives.png?raw=true" width="600" alt="screencap descriptives"></p>


* **RmimicMLAnova**: Function that computes a SPSS style univariate ANOVA with effect size and confidence intervals using a multi-level model from the lme4 function. Main effects and interactions are automatically decomposed using the
specified post-hoc corrections. 
```r
    workingdatabase <- Rmimic::alertness
    workingdatabase <- workingdatabase[which(workingdatabase$Condition == 'Condition2'),]
    anovaresult <- Rmimic::RmimicMLAnova(data = workingdatabase,dependentvariable = "Alertness",
                              subjectid = "PartID", between = "Group", within = c("Time"), 
                              randomintercept = c("PartID"))
```
<p align="center"><img src="/screencaps/screencap_RmimicMLAnova1.png?raw=true" width="600" alt="screencap RmimicMLAnova1"></p>
<p align="center"><img src="/screencaps/screencap_RmimicMLAnova2.png?raw=true" width="600" alt="screencap RmimicMLAnova2"></p>
<p align="center"><img src="/screencaps/screencap_RmimicMLAnova3.png?raw=true" width="600" alt="screencap RmimicMLAnova3"></p>
<p align="center"><img src="/screencaps/screencap_RmimicMLAnova4.png?raw=true" width="600" alt="screencap RmimicMLAnova4"></p>



## Main Function List
* **clusterthreshold1d**: Function that calculates contiguous clusters of locations in a 1D array that are
above or below some threshold and of some minimum cluster size (i.e., a cluster of 30 points all below 0.05).
```r
    tempdata <- c(0.2, 0.3, 0.05, 0.04, 0.06, 0.08, 0.009, 0.05, 0.02, 0.03, 0.08, 0.1, 0.4)  
    Rmimic::clusterthreshold1d(tempdata, crit = 0.05, clustersize = 3, direction = 'LessThan')
    # returns: c(0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0)
```

* **computechange**: Function to compute the change or difference between within subjects conditions. A new 
variable is returned to the data containing the change calculation. Optional parameters are included for
computing PercentChange, specifying the factor level to use as the baseline, and for other factors to control
for.
```r
    mockdatabase <- data.frame("ID" = rep_len(1:20,60), 
                            "Time" = c(rep_len("Time1",20),rep_len("Time2",20),rep_len("Time3",20)),
                            "X" = runif(60))
    mockdatabase <- Rmimic::computechange(mockdatabase, dependentvariable='X', subjectid='ID', within='Time')
```

* **identifyoutliers**: Function to identify outliers based upon the interquartile range (as SPSS does for
boxplots) and replace those values with NA.
```r
    tempdata <- runif(100,1,10); tempdata[6] <- 1000; tempdata[10] <- 1000  
    tempdata <- Rmimic::identifyoutliers(tempdata, iqrlimit = 3, verbose=TRUE)  
```

* **multipleimputation**: Function that uses the mice package to replace missing data points.
```r
    tempdata <- data.frame("X" = runif(100), "Y" = runif(100))  
    tempdata[6,'X'] <- NA; tempdata[10,'Y'] <- NA;   
    tempdata <- Rmimic::multipleimputation(tempdata, imputations=10)  
```

* **ezANOVA2text**: Function to output ezANOVA results in APA style format with effect sizes and confidence
intervals.
```r
    result <- ez::ezANOVA(data=elashoff,dv=Alertness,wid=PartID,
                            between=Group,within=.(Drug,Dose),type=3,detailed=TRUE,return_aov=TRUE)
    result <- Rmimic::ezANOVA2text(result, numparticipants=16, feffect="Generalized Eta Squared", 
                            sphericity="Greenhouse-Geisser", confidenceinterval=0.95, studywiseAlpha=0.05)
```

* **lmer2text**: Function to output lmerTest::lmer results in APA style format with effect sizes and confidence
intervals. The function will also work with lme4:lmer models, however these models do not report probability values.
```r
    fit <- lmerTest::lmer(Alertness ~ Group*Drug*Dose + (1 | PartID), data=Rmimic::elashoff)
    result <- Rmimic::lmer2text(fit, df="Kenward-Roger", numparticipants=16, numfactors=4)
```

* **ttest2text**: Function that takes a t-test result from the stats package and outputs the t-test result for
use in an APA style manuscript (i.e., t(18) = 2.3, p = 0.031) with proper rounding. Supports independent and
paired t-tests for both parametric and nonparametric data.
```r
    comparison1 <- c(2, 3, 4, 6, 8, 9, 10, 11, 13, 15); comparison2 <- c( 5, 7, 9, 10, 13, 15, 16, 17, 18, 20)  
    ttestresult <- stats::t.test(x= comparison1, y=comparison2, paired=FALSE, var.equal=TRUE)  
    ttestresult$effectsize <- ttestresult$statistic * sqrt((1/length(comparison1)) + (1/length(comparison2)))  
    tempout <- Rmimic::ttest2text(ttestresult, verbose=TRUE)  
    #returns: t(18) = 2.3, p = 0.031, ds = -1.04.
```

## Accessory Function List
* **cell2span**: Function that takes a vector of data table cells and creates a character string of a particular
length. The text empty will result in an open span in the output.
* **decimalplaces**: Function to obtain the number of decimal places of precision in a vector.
* **determineallpossiblecombinations**: Function to determine all possible combinations of an input array. For
instance, an array containing A, B, and C could be assessed looking at A, B, C, A:B, A:C, B:C, or A:B:C.
* **fuzzyP**: Function to round P values for reporting. Because reporting p = 0.912 to three digits of precision
is a bit silly.
* **posthoc2text**: Function to output ANOVA posthoc results in APA style format with effect sizes and confidence
intervals.
* **table2console**: Function that takes a data table and prints a nicely formatted table to the console.
* **typewriter**: Function to control the text outputted to the console to create a consistent indent.
