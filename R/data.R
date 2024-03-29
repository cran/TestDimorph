#' List format for the baboon.parms_df for multivariate analysis
#'
#' @format A list of 5 matrices (R.res, M.mu, F.mu, M.sdev, and F.sdev) and
#' two vectors (m and f)
#' \describe{ \item{R.res}{pooled within group
#' correlation matrix}
#' \item{M.mu}{Means of LDL and apo B in different sub-species for males}
#' \item{F.mu}{Means of LDL and apo B in different sub-species for females}
#' \item{m}{Male sample sizes} \item{f}{Female sample sizes}
#' \item{M.sdev}{Standard deviations for males} \item{F.sdev}{Standard
#' deviations for females} }
#' @seealso \link{baboon.parms_df}
"baboon.parms_list"


#' data frame format for the baboon.parms_df for multivariate analysis
#'
#'
#' @description  A dataset containing summary statistics for low density lipoprotein (LDL) and
#' apolipoprotein B (apo B) levels in 604 baboons measured on two different diets:
#' a basal diet and a high cholesterol, saturated fat diet. The baboons were
#' classified into one of two subspecies and a hybrid of the two subspecies
#' (Papio hamadryas anubis, P.h. cynocephalus, or hybrid). Each animal was
#' measured on each of the two diets.
#'
#' @format A data frame with 12 rows and 8 variables \describe{
#' \item{Trait}{Apolipoprotein B and LDL on two diets} \item{Sub}{Sub-species or hybrid}
#' \item{M.mu}{Means of LDL and apo B in different sub-species for males}
#' \item{F.mu}{Means of LDL and apo B in different sub-species for females}
#' \item{m}{Male sample sizes} \item{f}{Female sample sizes}
#' \item{M.sdev}{Standard deviations for males} \item{F.sdev}{Standard
#' deviations for females}
#'
#'  }
#'
#' @note The baboon data collection were supported by NIH grant HL28972 and
#' NIH contract HV53030 to the Southwest Foundation for Biomedical Research
#' (Now: Texas Biomedical Research Institute), and funds from the Southwest
#' Foundation for Biomedical Research
#' @references
#'
#' Konigsberg LW (1991). An historical note on the t-test for differences in
#' sexual dimorphism between populations. American journal of physical
#' anthropology, 84(1), 93–96.
#'
#'
"baboon.parms_df"

#' Pooled within group correlation matrix for baboon data
#' @format A 4*4 numerical matrix
#' @seealso \link{baboon.parms_list}
"baboon.parms_R"

#' The Howells' craniometric data
#'
#' A subset of a dataset that consists of 82 craniometric measurements taken
#' from approximately two thousands and half human crania from 28
#' geographically diverse populations. The full data set can be found in
#' \url{https://rdrr.io/github/geanes/bioanth/man/howell.html}
#'
#' @format A data frame with 441 rows and 10 variables:
#'  \describe{
#' \item{Sex}{'M' for male and 'F' for female}
#' \item{Pop}{Populations' names}
#' \item{GOL}{Glabello occipital length}
#' \item{NOL}{Nasio occipital length}
#' \item{BNL}{Bastion nasion length}
#' \item{BBH}{Basion bregma height}
#' \item{XCB}{Maximum cranial breadth}
#' \item{XFB}{Maximum frontal breadth}
#' \item{ZYB}{Bizygomatic breadth}
#' \item{AUB}{Biauricular breadth}
#'   }
#' @references
#'
#' Howells WW. (1989). Skull Shapes and the Map. Craniometric Analyses in the
#' Dispersion of Modern Homo. Papers of the Peabody Museum of Archaeology and
#' Ethnology, vol. 79, pp. 189. Cambridge, Mass.: Peabody Museum.
#'
#' Howells WW. (1995). Who's Who in Skulls. Ethnic Identification of Crania from
#' Measurements. Papers of the Peabody Museum of Archaeology and Ethnology,
#' vol. 82, pp. 108. Cambridge, Mass.: Peabody Museum.
#'
#' Howells, W. W. (1973). Cranial Variation in Man: A Study by Multivariate
#' Analysis of Patterns of Difference Among Recent Human Populations (Vol. 67).
#' Cambridge, MA: Peabody Museum of Archaeology and Ethnology.
#'
#' Howells, W. W. (1996). Notes and Comments: Howells' craniometric data on the
#' internet. American Journal of Physical Anthropology, 101(3), 441-442
"Howells"

#' Summary of the Howells' craniometric data
#'
#' Summary statistics of the Howells' data subset.
#'
#' @format A data frame with 32 rows and 8 variables:
#' \describe{
#' \item{Trait}{Measured feature} \item{Pop}{Population name}
#' \item{M.mu}{Means of males}
#' \item{F.mu}{Means of females}
#' \item{m}{Male sample sizes} \item{f}{Female sample sizes}
#' \item{M.sdev}{Standard deviations for males} \item{F.sdev}{Standard
#' deviations for females}
#' }
#' @references
#'
#' \link{Howells}
#'
"Howells_summary"
#' List format of \link{Howells_summary} for multivariate analysis
#'
#' @format A list of 5 matrices (R.res, M.mu, F.mu, M.sdev, and F.sdev) and
#' two vectors (m and f) with structure similar to \link{baboon.parms_list}
#'
"Howells_summary_list"

#' Pooled within group correlation matrix for Howells' data
#' @format A 8*8 numerical matrix
"Howells_R"
#' Pooled within-group variance-covariance matrix for Howells' data
#' @format A 8*8 numerical matrix
#' @seealso \link{Howells}
"Howells_V"

#' Measurements from calcined postcranial materials.
#'
#' Part of Table 3 from Cavazzuti et al. (2019).
#'
#' @format A data frame with 22 rows and 8 variables:
#' \describe{
#' \item{Trait}{Measured feature}
#' \item{M.mu}{Means of males}
#' \item{F.mu}{Means of females}
#' \item{m}{Male sample sizes} \item{f}{Female sample sizes}
#' \item{M.sdev}{Standard deviations for males} \item{F.sdev}{Standard
#' deviations for females}\item{D}{published value for Chakraborty and Majumder's
#' (1982) measure of sexual dimorphism.}
#' }
#' @references
#'
#' Cavazzuti, Claudio, et al. (2019) "Towards a new osteometric method for sexing
#' ancient cremated human remains. Analysis of Late Bronze Age and Iron Age samples
#' from Italy with gendered grave goods." PloS one 14.1: e0209423.
#'
#' Chakraborty, R., & Majumder, P. P. (1982). On Bennett's measure of sex
#' dimorphism. American journal of physical anthropology, 59(3), 295-298.
"Cremains_measurements"

#' Australia
#'
#' Raw data from Joseph Birdsell's 1938 survey. The data is from two regions
#' (B1 and B19), see Gilligan and Bulbeck (2007) for a map of the regions.
#' Data downloaded from Dr. Peter Brown's website:
#' \url{https://www.peterbrown-palaeoanthropology.net/resource.html}
#'
#'
#' @format A data frame with 94 rows and 9 variables:
#'  \describe{
#' \item{Pop}{(Region) ("B1" = Southwest Australia, "B19" = Northeast Australia),
#'  see Gilligan and Bulbeck (2007)}
#' \item{Sex}{Sex coded as "F" or "M"}
#' \item{Weight.kg}{body weight in kilograms}
#' \item{Stature.mm}{Standing height in millimeters}
#' \item{Hum.Lgth}{Humeral length in millimeters}
#' \item{Rad.Lgth}{Radius length in millimeters}
#' \item{Fem.Lgth}{Femoral length in millimeters}
#' \item{Tib.Lgth}{Tibial length in millimeters}
#' \item{Bi.illiac}{Bi-illiac breadth in millimeters}
#'   }
#' @references
#'
#' Gilligan, I., & Bulbeck, D. (2007). Environment and morphology in Australian
#' Aborigines: A re-analysis of the Birdsell database. American Journal of Physical
#' Anthropology, 134(1), 75-91.
#'
"Australia"

#' NHANES 1999
#'
#' @description Raw data from 1999-2000 NHANES (National Health and Nutrition
#' Examination Survey). Centers for Disease Control and Prevention (CDC). National
#' Center for Health Statistics (NCHS). National Health and Nutrition Examination
#' Survey Data. Hyattsville, MD: U.S. Department of Health and Human Services, Centers
#' for Disease Control and Prevention, 2020, \url{https://www.cdc.gov/nchs/nhanes/index.htm}
#'
#'
#' @format A data frame with 1430 rows and 5 variables:
#'  \describe{
#' \item{Sex}{(RIAGENDR) Sex coded as "F" or "M"}
#' \item{Pop}{(RIDRETH1) Self-reported race, coded as "Black" = Non-Hispanic Black,
#' "Mex.Am" = Mexican American, or "White" = Non-Hispanic White}
#' \item{BMXWT}{Body weight in kilograms}
#' \item{BMXHT}{Standing height in centimeters}
#' \item{BMXARML}{Upper arm length in centimeters}
#'   }
#' @note
#'
#' This is not the complete dataset. It is selected so that age in years is
#' greater than or equal to 20 and less than or equal to 40
#'
"NHANES_1999"

#' Heuristic data
#' @description Heuristic data from Fidler and Thompson (2001)
#' @format A data frame with 24 rows and 3 variables:
#'  \describe{
#' \item{Sex}{'M' for male and 'F' for female}
#' \item{Pop}{Populations' names}
#' \item{x}{Dependent variable}
#'   }
#' @references
#'
#' Fidler, Fiona, and Bruce Thompson. "Computing correct confidence intervals
#' for ANOVA fixed-and random-effects effect sizes." Educational and
#' Psychological Measurement 61.4 (2001): 575-604.
#'

"FT"

#' Hypothetical set of unbalanced data
#' @description Example data set from Shaw and Mitchell-Olds (1993)
#' @format A data frame with 11 rows and 3 variables:
#'  \describe{
#' \item{Sex}{'M' for male and 'F' for female}
#' \item{Pop}{Populations' names}
#' \item{x}{Dependent variable}
#'   }
#' @references
#'
#' Shaw, Ruth G., and Thomas Mitchell-Olds. "ANOVA for unbalanced data: an overview.
#' " Ecology 74.6 (1993): 1638-1645.
#'

"SMO"
