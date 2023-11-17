## ----global-options, include=FALSE--------------------------------------------
knitr::opts_chunk$set(comment=NULL)
knitr::opts_chunk$set(fig.width=6, fig.height=4)
suppressPackageStartupMessages(library(TestDimorph))

## ----Table 2--------------------------------------------------------------------------------------
Table.02=function () 
{
library(TestDimorph)
options(width=100) # This option just for output from Rmarkdown
NHANES_univariate<<-extract_sum(NHANES_1999,test='uni',run=FALSE) # BMXWT (Body mass)
univariate(NHANES_univariate,es_anova = "eta2",pairwise = TRUE)
}

Table.02()

## ----Table 3--------------------------------------------------------------------------------------
Table.03=function()
{
library(TestDimorph)
NHANES_multivariate<<-extract_sum(NHANES_1999,test='multi',run=FALSE)
multivariate(NHANES_multivariate)
}

Table.03()

## ----Table 4--------------------------------------------------------------------------------------
Table.04=function()
{
library(TestDimorph)
print(univariate(NHANES_univariate, type_anova='III'))
t_greene(NHANES_univariate,plot = TRUE,padjust ="fdr")
}

Table.04()

## ----Table 5--------------------------------------------------------------------------------------
Table.05=function()
{
library(TestDimorph)
to_van_Vark=extract_sum(Howells,test='van',run=F)
van_vark(to_van_Vark)
}

Table.05()

## ----Table 6--------------------------------------------------------------------------------------
Table.06=function () 
{
# Comparisons of femur head diameter in four populations
library(TestDimorph)
df <- data.frame(
  Pop = c("Turkish", "Bulgarian", "Greek", "Portuguese"),
  m = c(150.00, 82.00, 36.00, 34.00),
  f = c(150.00, 58.00, 34.00, 24.00),
  M.mu = c(49.39, 48.33, 46.99, 45.20),
  F.mu = c(42.91, 42.89, 42.44, 40.90),
  M.sdev = c(3.01, 2.53, 2.47, 2.00),
  F.sdev = c(2.90, 2.84, 2.26, 2.90)
)
print(aov_ss(x = df, CI=0.95),digits=6)
}

Table.06()

## ----Table 7,fig.width=6, fig.height=3------------------------------------------------------------
Table.07=function (i.which=13) 
{
library(TestDimorph)
print(MI_index(Cremains_measurements[i.which,],B=1000,rand=F,verbose=F,plot=T))
print(MI_index(Cremains_measurements[i.which,],index_type='NI',
   B=1000,rand=F,plot=T,verbose=F))
print(D_index(Cremains_measurements[i.which,],B=1000,rand=F,verbose=F,plot=T))
print(Hedges_g(Cremains_measurements[i.which,],B=1000,rand=F,verbose=F))
}

Table.07()

