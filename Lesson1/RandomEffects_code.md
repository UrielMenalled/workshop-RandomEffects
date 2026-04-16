Random effects demo
================

# Loading

``` r
library(lmerTest)
```

    ## Loading required package: lme4

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
library(performance)
library(glmmTMB)
library(DHARMa)
```

    ## This is DHARMa 0.4.7. For overview type '?DHARMa'. For recent changes, type news(package = 'DHARMa')

``` r
library(ADER)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.1.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ tidyr::expand() masks Matrix::expand()
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ✖ tidyr::pack()   masks Matrix::pack()
    ## ✖ tidyr::unpack() masks Matrix::unpack()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
options(contrasts = c("contr.sum", "contr.poly")) #setting contrasts
```

Today we will be using two data frames:

1.  A subset of data published in Menalled et al. (2023):
    <https://doi.org/10.1038/s41598-023-43987-x>

    This experiment looks at changes in weed communities across 10 cover
    crop treatments (5 summer cover crop treatments, 5 winter cover crop
    treatments). Today, we will just work with winter cover crops
    (i.e. tilled control, canola, cereal rye (CR), hairy vetch (HV), and
    HV x CR mix).

2.  Richness of species in 45 sampling stations along the coastline of
    The Netherlands, measured by two researches of the the RIKZ
    (Rijksinstituut voor Kust en Zee), i.e., the Dutch Institute for
    Coastal and Marine Management.

# Cover Crop example

``` r
#Each year is stored separately, I'll load them in, calculate total weed biomass, and combined the data frames
WintData1_tmp <- read.csv("https://ecommons.cornell.edu/server/api/core/bitstreams/adb24a0c-0b2c-4bd6-90ad-e203c270bba3/content")
WintData2_tmp <- read.csv("https://ecommons.cornell.edu/server/api/core/bitstreams/f05bf7b9-3c15-48ec-9f28-0c978cf0560b/content")

str(WintData1_tmp) #loaded in nicely
str(WintData2_tmp) #loaded in nicely

#We don't want to do weed community analysis today, just simple models of total weed biomass
WintData1 <-
  WintData1_tmp %>% 
  rowwise() %>%
  mutate(WeedBiomass = sum(c_across(-c(Trial, Year, Site, Block, Plot, CoverCrop, CoverCropBiomass, canola, 
                                       cereal.rye, hairy.vetch)), na.rm = FALSE)) %>% 
  ungroup() %>% 
  select(Year, Site, Block, Plot, CoverCrop,WeedBiomass)

WintData2 <-
  WintData2_tmp %>% 
  rowwise() %>%
  mutate(WeedBiomass = sum(c_across(-c(Trial, Year, Site, Block, Plot, CoverCrop, CoverCropBiomass, canola, 
                                       cereal.rye, hairy.vetch)), na.rm = FALSE)) %>%
  ungroup() %>% 
  select(Year, Site, Block, Plot, CoverCrop,WeedBiomass)

#Combine data frames
WintData_tmp <- bind_rows(WintData1,WintData2)

WintData <- WintData_tmp %>% mutate(across(-WeedBiomass, as.factor))
str(WintData) #looks good!
```

\#Modeling

I want to see how weed biomass is impacted by cover crop treatment,
basically WeedBiomass~CoverCrop. However, my experiment has clear
nesting! I must account for non-independence through a model with more
fixed and/or random effects. Looking at my experimental design (picture
below) we can see that there are **three** *total* levels.

1.  Site and Year are at an equal level because both sites are in both
    years and visa versa. there are 2 levels of each site and year

2.  Block is within Site and Year. I know this because we used different
    fields in each site year, so block 1 in Musgrave during year 1 **is
    not** the same as block 1 in Musgrave during year 2. Because each
    block is within a unique site year, there are 2 sites x 2 years x 4
    blocks = 16 groups.

3.  Cover crop treatments are within blocks, this is the level at which
    I made my observations.

As nice it would be to account for Site and Year through a random
effect, we don’t have enough levels of either group. Remember we need at
least three groups per variable for it to be a random effect. Thus,
these variables will be fixed effects. On the other hand, site.year has
4 groups (iffy, but let’s just try) and block has 16 levels across
site-years (nice, that’s a solid amount of groups for a RE). As for
CoverCrop, this is the level of observation, it does not need to be a
random effect.

<figure>
<img src="Example_MapDesign.png" width="452"
alt="Here is a nice picture of the map and experimental design of the data we are working with today" />
<figcaption aria-hidden="true">Here is a nice picture of the map and
experimental design of the data we are working with today</figcaption>
</figure>

``` r
#create Site.Year variable
WintData$Site.Year <- paste(WintData$Site,WintData$Year,sep = ".")
WintData$Site.Year <- as.factor(WintData$Site.Year)


#let's get a sense of what we are working with...
#We'll probably be working with a Tweedie glmer because 
plot(WeedBiomass~1, data = WintData) #a cluster of points with near zero weed biomass, probably 
```

![](RandomEffects_code_files/figure-gfm/modeling-1.png)<!-- -->

``` r
plot(WeedBiomass~CoverCrop*Site*Year, data = WintData) #non-normally distributed around treatments
```

![](RandomEffects_code_files/figure-gfm/modeling-2.png)<!-- -->![](RandomEffects_code_files/figure-gfm/modeling-3.png)<!-- -->![](RandomEffects_code_files/figure-gfm/modeling-4.png)<!-- -->

``` r
#let's get modeling

#simple lmer first
WeedBio_Mod1 <- lmer(WeedBiomass~Site*Year + CoverCrop + (1|Site.Year/Block), data = WintData)
```

    ## boundary (singular) fit: see help('isSingular')

    ## Warning in as_lmerModLT(model, devfun): Model may not have converged with 1
    ## eigenvalue close to zero: -1.3e-09

``` r
#singularity issues from RE of Block:Site.Year = 0, could be okay if anova DF are good
summary(WeedBio_Mod1) #RE captured correctly because number of groups in RE follows design
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: WeedBiomass ~ Site * Year + CoverCrop + (1 | Site.Year/Block)
    ##    Data: WintData
    ## 
    ## REML criterion at convergence: 728.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.8458 -0.4970  0.0584  0.3143  5.5951 
    ## 
    ## Random effects:
    ##  Groups          Name        Variance Std.Dev.
    ##  Block:Site.Year (Intercept)    0.00   0.000  
    ##  Site.Year       (Intercept)   71.32   8.445  
    ##  Residual                    1095.66  33.101  
    ## Number of obs: 79, groups:  Block:Site.Year, 16; Site.Year, 4
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)   29.088      5.632  71.000   5.165 2.11e-06 ***
    ## Site1         12.471      5.632  71.000   2.214 0.030007 *  
    ## Year1          2.263      5.632  71.000   0.402 0.688998    
    ## CoverCrop1   -20.480      7.414  71.000  -2.762 0.007306 ** 
    ## CoverCrop2   -25.782      7.414  71.000  -3.477 0.000869 ***
    ## CoverCrop3    -1.049      7.414  71.000  -0.141 0.887940    
    ## CoverCrop4   -27.448      7.414  71.000  -3.702 0.000419 ***
    ## Site1:Year1    3.147      5.632  71.000   0.559 0.578100    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) Site1  Year1  CvrCr1 CvrCr2 CvrCr3 CvrCr4
    ## Site1        0.006                                          
    ## Year1        0.006  0.006                                   
    ## CoverCrop1  -0.005 -0.005 -0.005                            
    ## CoverCrop2  -0.005 -0.005 -0.005 -0.246                     
    ## CoverCrop3  -0.005 -0.005 -0.005 -0.246 -0.246              
    ## CoverCrop4  -0.005 -0.005 -0.005 -0.246 -0.246 -0.246       
    ## Site1:Year1  0.006  0.006  0.006 -0.005 -0.005 -0.005 -0.005
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
plot(WeedBio_Mod1) #terrible, as expected
```

![](RandomEffects_code_files/figure-gfm/modeling-5.png)<!-- -->

``` r
check_model(WeedBio_Mod1) #I'm mainly interested in pp_check (looks bad) and random effects check (looks okay)
```

![](RandomEffects_code_files/figure-gfm/modeling-6.png)<!-- -->

``` r
anova(WeedBio_Mod1)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##           Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## Site        5373  5373.0     1    71  4.9039   0.03001 *  
    ## Year         177   176.9     1    71  0.1615   0.68900    
    ## CoverCrop 112957 28239.2     4    71 25.7737 3.271e-13 ***
    ## Site:Year    342   342.0     1    71  0.3122   0.57810    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
#now glmer
WeedBio_Mod2 <- glmmTMB(WeedBiomass~Site*Year + CoverCrop + (1|Site.Year/Block),
                        family = tweedie(link = "log"),
                        data = WintData)
summary(WeedBio_Mod2) #RE captured correctly because number of groups in RE follows design
```

    ##  Family: tweedie  ( log )
    ## Formula:          WeedBiomass ~ Site * Year + CoverCrop + (1 | Site.Year/Block)
    ## Data: WintData
    ## 
    ##       AIC       BIC    logLik -2*log(L)  df.resid 
    ##     503.3     531.7    -239.6     479.3        67 
    ## 
    ## Random effects:
    ## 
    ## Conditional model:
    ##  Groups          Name        Variance  Std.Dev. 
    ##  Block:Site.Year (Intercept) 2.786e-02 1.669e-01
    ##  Site.Year       (Intercept) 5.935e-10 2.436e-05
    ## Number of obs: 79, groups:  Block:Site.Year, 16; Site.Year, 4
    ## 
    ## Dispersion parameter for tweedie family (): 2.23 
    ## 
    ## Conditional model:
    ##             Estimate Std. Error z value Pr(>|z|)    
    ## (Intercept)  2.06149    0.16613  12.409  < 2e-16 ***
    ## Site1        0.77195    0.16085   4.799 1.59e-06 ***
    ## Year1       -0.03369    0.15743  -0.214   0.8306    
    ## CoverCrop1   0.03726    0.29381   0.127   0.8991    
    ## CoverCrop2  -1.46188    0.37533  -3.895 9.82e-05 ***
    ## CoverCrop3   0.72355    0.28169   2.569   0.0102 *  
    ## CoverCrop4  -1.91508    0.31750  -6.032 1.62e-09 ***
    ## Site1:Year1  0.09991    0.15392   0.649   0.5163    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
simulateResiduals(WeedBio_Mod2, plot = TRUE) #much better, good enough to publish
```

![](RandomEffects_code_files/figure-gfm/modeling-7.png)<!-- -->

    ## Object of Class DHARMa with simulated residuals based on 250 simulations with refit = FALSE . See ?DHARMa::simulateResiduals for help. 
    ##  
    ## Scaled residual values: 0.752 0.788 0.356 0.636 0.232 0.416 0.3 0.952 0.208 0.58 0.124 1 0.084 0.396 0.116 0.084 0.188 0.01347419 0.716 0.412 ...

``` r
check_model(WeedBio_Mod2) #pp_check still shows some error around 0. All RE look good. Good enough!
```

    ## `check_outliers()` does not yet support models of class `glmmTMB`.

![](RandomEffects_code_files/figure-gfm/modeling-8.png)<!-- -->

``` r
car::Anova(WeedBio_Mod2, type = 3)
```

    ## Analysis of Deviance Table (Type III Wald chisquare tests)
    ## 
    ## Response: WeedBiomass
    ##                Chisq Df Pr(>Chisq)    
    ## (Intercept) 153.9873  1  < 2.2e-16 ***
    ## Site         23.0327  1  1.593e-06 ***
    ## Year          0.0458  1     0.8306    
    ## CoverCrop   130.3030  4  < 2.2e-16 ***
    ## Site:Year     0.4213  1     0.5163    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

# Beach example

The RIKZ data set is a great way to look at denominator degrees of
freedom. In this experiment, researchers looked at species richness
(Richness) across the altitude (NAP) and sun exposure (Exposure).
Exposure is the same with each beach. Because you have repeated
measurements at the beach level, non-independence at this level should
be accounted for, let’s use a random effect. In this case, observations
at this level (Exposure) should have denominator degrees of freedom that
less than the beach level (i.e., less than 9).

<img src="RIKZ_Design.png" width="466" />

``` r
data("RIKZ")
RIKZ_Wide <- RIKZ
str(RIKZ_Wide) #C2 - I5 are species
```

    ## 'data.frame':    45 obs. of  89 variables:
    ##  $ Sample       : int  1 2 3 4 5 6 7 8 9 10 ...
    ##  $ C1           : num  4 0 0 0 1 0 0 0 0 0 ...
    ##  $ P1           : num  0 0 0 0 0 0 0 0 0 2 ...
    ##  $ P2           : num  0 1 3 0 0 0 0 0 0 0 ...
    ##  $ P3           : num  0 0 0 0 0 0 1 0 21 2 ...
    ##  $ P4           : num  0 0 0 0 0 0 1 0 11 0 ...
    ##  $ P5           : num  0 0 0 0 0 1 0 0 3 0 ...
    ##  $ P6           : num  1 0 0 0 0 0 0 0 0 0 ...
    ##  $ P7           : num  0 0 0 0 0 0 0 0 34 7 ...
    ##  $ P8           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ P9           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ P10          : num  0 2 1 0 1 0 0 0 0 2 ...
    ##  $ P11          : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ P12          : num  0 0 0 0 0 1 3 0 14 0 ...
    ##  $ P13          : num  0 0 0 0 0 0 0 5 16 0 ...
    ##  $ P14          : num  0 0 0 0 0 889 21 10 49 0 ...
    ##  $ P15          : num  4 0 0 0 0 0 0 0 1 0 ...
    ##  $ P16          : num  32 22 16 79 15 0 3 7 1 0 ...
    ##  $ P17          : num  0 0 0 0 0 0 0 0 0 16 ...
    ##  $ P18          : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ P19          : num  0 0 1 0 0 0 0 0 1 0 ...
    ##  $ P20          : num  0 0 0 0 0 0 0 0 0 2 ...
    ##  $ P21          : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ P22          : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ P24          : num  0 0 0 92 0 0 0 0 0 0 ...
    ##  $ P25          : num  50 12 7 1 0 0 0 0 0 0 ...
    ##  $ N1           : num  0 0 0 0 0 0 0 1 0 0 ...
    ##  $ CR1          : num  2 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR2          : num  0 0 0 2 0 0 0 0 0 0 ...
    ##  $ CR3          : num  0 1 2 0 25 1 0 0 0 0 ...
    ##  $ CR4          : num  0 8 25 0 13 0 0 0 0 0 ...
    ##  $ CR5          : num  0 0 0 1 1 0 175 163 0 0 ...
    ##  $ CR6          : num  0 0 0 1 0 0 0 0 0 0 ...
    ##  $ CR7          : num  0 0 1 0 0 0 0 0 0 0 ...
    ##  $ CR8          : num  0 0 0 0 0 0 0 0 0 1 ...
    ##  $ CR9          : num  0 0 0 1 4 0 0 0 0 0 ...
    ##  $ CR10         : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR11         : num  0 0 0 0 0 0 0 0 2 0 ...
    ##  $ CR12         : num  0 0 0 0 0 25 35 4 17 0 ...
    ##  $ CR13         : num  0 0 0 0 0 0 0 0 1 0 ...
    ##  $ CR14         : num  30 3 2 4 2 0 0 1 1 1 ...
    ##  $ CR15         : num  9 0 0 15 4 0 1 1 0 0 ...
    ##  $ CR16         : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR17         : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR18         : num  0 0 2 0 0 0 0 0 0 1 ...
    ##  $ CR19         : num  8 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR20         : num  0 0 2 0 0 0 0 0 0 0 ...
    ##  $ CR21         : num  0 1 7 0 1 0 0 0 0 0 ...
    ##  $ CR22         : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR23         : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR24         : num  0 0 0 0 0 0 0 0 0 1 ...
    ##  $ CR25         : num  0 1 0 0 0 0 0 0 0 0 ...
    ##  $ CR26         : num  2 0 0 1 0 0 0 0 0 0 ...
    ##  $ CR27         : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ CR28         : num  0 0 1 0 0 0 0 0 0 0 ...
    ##  $ M1           : num  0 0 0 0 0 0 0 0 2 1 ...
    ##  $ M2           : num  0 0 0 0 0 0 0 0 2 0 ...
    ##  $ M3           : num  0 0 0 0 0 0 0 0 0 1 ...
    ##  $ M4           : num  0 0 0 0 0 0 0 0 0 1 ...
    ##  $ M5           : num  0 0 0 0 0 0 0 0 0 1 ...
    ##  $ M6           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ M7           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ M8           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ M9           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ M10          : num  0 0 0 0 0 0 0 0 16 6 ...
    ##  $ M11          : num  0 0 0 0 0 0 1 0 18 0 ...
    ##  $ M12          : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ M14          : num  0 1 0 0 0 0 0 0 0 1 ...
    ##  $ M15          : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ M16          : num  0 0 0 0 0 0 0 0 0 2 ...
    ##  $ M17          : num  0 0 0 0 0 0 0 0 1 0 ...
    ##  $ I1           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ I2           : num  0 0 0 0 0 1 0 0 0 0 ...
    ##  $ I3           : num  1 0 0 2 0 24 0 0 0 0 ...
    ##  $ I4           : num  0 0 0 0 0 2 0 0 0 0 ...
    ##  $ I5           : num  0 0 0 0 0 0 0 0 0 0 ...
    ##  $ week         : int  1 1 1 1 1 1 1 1 1 1 ...
    ##  $ angle1       : int  32 62 65 55 23 129 126 52 26 143 ...
    ##  $ angle2       : int  96 96 96 96 96 89 89 89 89 89 ...
    ##  $ exposure     : int  10 10 10 10 10 8 8 8 8 8 ...
    ##  $ salinity     : num  29.4 29.4 29.4 29.4 29.4 29.6 29.6 29.6 29.6 29.6 ...
    ##  $ temperature  : num  17.5 17.5 17.5 17.5 17.5 20.8 20.8 20.8 20.8 20.8 ...
    ##  $ NAP          : num  0.045 -1.036 -1.336 0.616 -0.684 ...
    ##  $ penetrability: num  254 227 237 249 252 ...
    ##  $ grainsize    : num  222 200 194 221 202 ...
    ##  $ humus        : num  0.05 0.3 0.1 0.15 0.05 0.1 0.1 0.1 0.15 0 ...
    ##  $ chalk        : num  2.05 2.5 3.45 1.6 2.45 2.5 1.85 1.7 2.3 2.6 ...
    ##  $ sorting1     : num  69.8 59 59.2 67.8 57.8 ...
    ##  $ Beach        : int  1 1 1 1 1 2 2 2 2 2 ...

``` r
RIKZ_Long <-
  RIKZ_Wide %>% 
  mutate(across(C1:I5, ~ if_else(. > 0, 1, 0)),
         Sample = rep(1:5,9),
         Richness = rowSums(across(C1:I5)),
         Exposure = exposure) %>% 
  select(Beach, Sample, Exposure, NAP, Richness) 

RIKZ_Mod1 <- lmer(Richness ~ NAP + Exposure +  (1|Beach), data = RIKZ_Long)
summary(RIKZ_Mod1)
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Richness ~ NAP + Exposure + (1 | Beach)
    ##    Data: RIKZ_Long
    ## 
    ## REML criterion at convergence: 225.5
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3462 -0.5025 -0.2306  0.2150  4.2746 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Beach    (Intercept) 0.3378   0.5812  
    ##  Residual             9.3735   3.0616  
    ## Number of obs: 45, groups:  Beach, 9
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)  37.2974     5.5547  6.7064   6.715 0.000330 ***
    ## NAP          -2.6941     0.4700 41.5377  -5.732 9.98e-07 ***
    ## Exposure     -3.0005     0.5417  6.7290  -5.539 0.000996 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) NAP   
    ## NAP       0.018       
    ## Exposure -0.996 -0.047

``` r
anova(RIKZ_Mod1)
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## NAP      307.99  307.99     1 41.538  32.857  9.98e-07 ***
    ## Exposure 287.54  287.54     1  6.729  30.676 0.0009955 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

When you have observations at multiple levels , like what we see in the
RIKZ data, singularity in the random effects can cause issues with
denominator degrees of freedom calculation. This is a *huge deal* as
ANOVA is using degrees of freedom to calculate p-values.

``` r
#randomly shuffling beach around so that it accounts for no variance
RIKZ_Long <-
  RIKZ_Long %>% 
  mutate(Beach2 = rep(sample(1:5), 9))
RIKZ_Mod2 <- lmer(Richness ~ NAP + Exposure + (1|Beach2),
                  data = RIKZ_Long)
```

    ## boundary (singular) fit: see help('isSingular')

``` r
summary(RIKZ_Mod2) #singularity warning shows that beach accounts for no variance
```

    ## Linear mixed model fit by REML. t-tests use Satterthwaite's method [
    ## lmerModLmerTest]
    ## Formula: Richness ~ NAP + Exposure + (1 | Beach2)
    ##    Data: RIKZ_Long
    ## 
    ## REML criterion at convergence: 225.6
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -1.3869 -0.5507 -0.2733  0.2471  4.2901 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev.
    ##  Beach2   (Intercept) 0.000    0.000   
    ##  Residual             9.649    3.106   
    ## Number of obs: 45, groups:  Beach2, 5
    ## 
    ## Fixed effects:
    ##             Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)  37.2909     5.1878 42.0000   7.188 7.83e-09 ***
    ## NAP          -2.7252     0.4716 42.0000  -5.779 8.26e-07 ***
    ## Exposure     -2.9988     0.5060 42.0000  -5.926 5.07e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##          (Intr) NAP   
    ## NAP       0.019       
    ## Exposure -0.996 -0.051
    ## optimizer (nloptwrap) convergence code: 0 (OK)
    ## boundary (singular) fit: see help('isSingular')

``` r
anova(RIKZ_Mod2) #df for exposure is WRONGE, notice that Satterthwaite's method is used
```

    ## Type III Analysis of Variance Table with Satterthwaite's method
    ##          Sum Sq Mean Sq NumDF DenDF F value    Pr(>F)    
    ## NAP      322.24  322.24     1    42  33.396 8.257e-07 ***
    ## Exposure 338.86  338.86     1    42  35.118 5.074e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
anova(RIKZ_Mod2, ddf = "Kenward-Roger") #Kenward-Roger denominator degrees of freedom will help solve this issue, often for real data—when there is some very small amount of variance—the impact of a Kenward-Roger is stronger. Unfortunately, there isn't a way to do a Kenward-Rodger in car::Anova(), at least to my knowledge...
```

    ## Type III Analysis of Variance Table with Kenward-Roger's method
    ##          Sum Sq Mean Sq NumDF  DenDF F value    Pr(>F)    
    ## NAP      313.96  313.96     1 39.683  32.538 1.262e-06 ***
    ## Exposure 338.83  338.83     1 38.005  35.116 7.191e-07 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
