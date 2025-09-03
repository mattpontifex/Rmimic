Rmimic
==============

**<p align="center" style="padding-left: 2em; padding-right: 2em;">This package simply provides wrappers around well-established and validated packages (e.g. stats, lawstat, psychometric, MBESS, lme4, psych, epitools) to make their use and output more user friendly.</p>**
------------

Why should I use this package?
------------
Because using R for statistics should not be hard. This package:
* Provides a standard operating approach for running statistics in R with results that generally align with 
those that would be obtained through SPSS, SAS, and/or Systat.
* Provides automatic posthoc decomposition of significant effects.
* Provides statistical results, effect sizes, effect size confidence intervals, and interpretations in APA format.
* Provides a graphical user interface for point and click analysis.


What is this?
------------
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

* **lmerPosthoc**: Function that performs posthoc decomposition of univariate ANOVA with effect size and confidence intervals
 using a multi-level model from the lme4 function. Interactions are decomposed multiple ways (A holding B, B holding A) and
 superseeding interactions suppress lower level effects tests (no posthoc test of A:B if A:B:C is significant). Tests of A:B
 can still be obtained using the planned parameter if desired. Optional ability to perform bootstrapping of the model is also
 included. Main effects and interactions are automatically decomposed using the specified post-hoc corrections. 
```r
    workingdatabase <- Rmimic::alertness
    workingdatabase <- workingdatabase[which(workingdatabase$Condition == 'Condition2'),]
    fit <- lmerTest::lmer(Alertness ~ Group*Time + (1 | PartID), data = workingdatabase)
    results <- Rmimic::lmerPosthoc(fit, dependentvariable="Alertness", subjectid='PartID',
                between=c('Group'), within=c('Time'), covariates=NULL,
                planned=NULL, posthoccorrection="False Discovery Rate Control", 
                bootstrap = list('repetitions'=100, 'subsample'=0.96, 'method'='default'),
                progressbar=TRUE)
    Rmimic::lmerEffectsSummarize(results, tag='', show='html', outputfile="test.html")
```
<p align="center"><img src="/screencaps/screencap_lmerPosthoc1.png?raw=true" width="700" alt="screencap lmerPosthoc1"></p>
<p align="center"><img src="/screencaps/screencap_lmerPosthoc2.png?raw=true" width="700" alt="screencap lmerPosthoc2"></p>
<p align="center"><img src="/screencaps/screencap_lmerPosthoc3.png?raw=true" width="700" alt="screencap lmerPosthoc3"></p>
<p align="center"><img src="/screencaps/screencap_lmerPosthoc4.png?raw=true" width="700" alt="screencap lmerPosthoc4"></p>


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
    chisquareresult <- Rmimic::RmimicChisquare(x='Age', y='Pet',  data=tempdata, 
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



* **mediate2text**: Function that outputs mediation results in a more intelligible format. 
```r
    workingdatabase <- Rmimic::gradwakefulness
    fitM <- stats::lm(CaffeineConsumption ~ HoursAwake, data=workingdatabase)
    fitY <- stats::lm(Wakefulness ~ HoursAwake + CaffeineConsumption, data=workingdatabase)
    fitMed <- mediation::mediate(fitM, fitY, treat='HoursAwake', mediator='CaffeineConsumption',
                              boot=FALSE, sims=1000, conf.level=0.95)
    res <- Rmimic::mediate2text(fitMed, studywiseAlpha=0.05)
```
<p align="center"><img src="/screencaps/screencap_mediate2text.png?raw=true" width="600" alt="screencap mediate2text"></p>


