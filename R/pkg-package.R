#' hdm: High-Dimensional Metrics
#'
#' This package implements methods for estimation and inference in a high-dimensional setting.
#'
#' \tabular{ll}{ Package: \tab hdm\cr Type: \tab Package\cr Version: \tab
#' 0.1\cr Date: \tab 2015-05-25\cr License: \tab GPL-3\cr } This package provides efficient estimators 
#' and  uniformly valid confidence intervals for various low-dimensional causal/structural parameters  
#' appearing in high-dimensional approximately sparse models. The package
#' includes functions for fitting heteroskedastic robust Lasso regressions with non-Gaussian erros and 
#' for instrumental variable (IV) and treatment effect estimation in a
#' high-dimensional setting. Moreover, the methods enable valid post-selection
#' inference. Moreover, a theoretically grounded, data-driven choice of the penalty level is provided.
#'
#' @name hdm-package
#' @aliases hdm-package hdm
#' @docType package
#' @author Victor Chernozhukov, Christian Hansen, Martin Spindler\cr
#'
#' Maintainer: Martin Spindler <spindler@@mea.mpisoc.mpg.de>
#' @references A. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369-2429.
#' @references A. Belloni, V. Chernozhukov and C. Hansen (2013). Inference for
#' high-dimensional sparse econometric models. In Advances in Economics and
#' Econometrics: 10th World Congress, Vol. 3: Econometrics, Cambirdge
#' University Press: Cambridge, 245-295.
#' @references A. Belloni, V. Chernozhukov, C. Hansen (2014). Inference on
#' treatment effects after selection among high-dimensional controls. The
#' Review of Economic Studies 81(2), 608-650.
#' @importFrom stats binomial
#' @importFrom stats coef
#' @importFrom stats confint 
#' @importFrom stats cor
#' @importFrom stats glm
#' @importFrom stats lm
#' @importFrom stats model.frame
#' @importFrom stats model.matrix
#' @importFrom stats model.response
#' @importFrom stats pnorm
#' @importFrom stats predict
#' @importFrom stats printCoefmat
#' @importFrom stats qnorm
#' @importFrom stats qt
#' @importFrom stats quantile
#' @importFrom stats rexp
#' @importFrom stats rnorm
#' @importFrom stats var
#' @importFrom stats vcov
#' @importFrom stats as.formula
#' @importFrom stats contrasts
#' @importFrom stats formula
#' @importFrom stats terms
#' @importFrom stats update
#' @importFrom Formula model.part
#' @importFrom methods is
NULL

#' Growth data set
#'
#' Data set of growth compiled by Barro Lee.
#'
#' The data set contains growth data of Barro-Lee. The Barro Lee data consists
#' of a panel of 138 countries for the period 1960 to 1985. The dependent
#' variable is national growth rates in GDP per capita for the periods
#' 1965-1975 and 1975-1985. The growth rate in GDP over a period from \eqn{t_1} to \eqn{t_2}
#' is commonly defined as \eqn{\log(GDP_{t_1}/GDP_{t_2})}. The number of covariates is p=62.
#' The number of complete observations is 90.
#'
#' @source The full data set and further details can be found at
#' \url{http://www2.nber.org/pub/barro.lee/} and,
#' \url{https://www.bristol.ac.uk//Depts//Economics//Growth//barlee.htm}.
#'
#' @name Growth Data
#' @aliases Growth Example GrowthData GDP
#' @docType data
#' @format Dataframe with the following variables:
#' \describe{
#' \item{outcome}{dependent variable: national growth rates in GDP per capita
#' for the periods 1965-1975 and 1975-1985}
#' \item{x}{covariates which might influence growth}}
#' @references R.J. Barro, J.W. Lee (1994). Data set for a panel of 139
#' countries. NBER.
#'
#' @references R.J. Barro, X. Sala-i-Martin (1995). Economic Growth. McGrwa-Hill, New York.
#' @keywords datasets Grwoth GDP
#' @examples
#' data(GrwothData)
NULL


#' Pension 401(k) data set
#'
#' Data set on financial wealth and 401(k) plan participation
#'
#' The sample is drawn from the 1991 Survey of Income and Program Participation
#' (SIPP) and consists of 9,915 observations. The observational units are
#' household reference persons aged 25-64 and spouse if present. Households are
#' included in the sample if at least one person is employed and no one is
#' self-employed. The data set was analysed in Chernozhukov and Hansen (2004)
#' and Belloni et al. (2014) where further details can be found. They examine
#' the effects of 401(k) plans on wealth using data from the Survey of Income
#' and Program Participation using 401(k) eligibility as an instrument for
#' 401(k) participation.
#'
#' @name pension
#' @aliases 401(k) plans wealth data pension
#' @docType data
#' @format Dataframe with the following variables (amongst others): \describe{
#' \item{p401}{participation in 401(k)}
#' \item{e401}{eligibility for 401(k)}
#' \item{a401}{401(k) assets}
#' \item{tw}{total wealth (in US $)}
#' \item{tfa}{financial assets (in US $)}
#' \item{net_tfa}{net financial assets (in US $)}
#' \item{nifa}{non-401k financial assets (in US $)}
#' \item{net_nifa}{net non-401k financial assets}
#' \item{net_n401}{net non-401(k) assets (in US $)}
#' \item{ira}{individual retirement account (IRA)}
#' \item{inc}{income (in US $)}
#' \item{age}{age}
#' \item{fsize}{family size}
#' \item{marr}{married}
#' \item{pira}{participation in IRA} \item{db}{defined benefit pension}
#' \item{hown}{home owner} \item{educ}{education (in years)} \item{male}{male}
#' \item{twoearn}{two earners} %\item{i1-i7}{} %\item{a1-a5}{}
#' \item{nohs, hs, smcol, col}{dummies for education: no high-school, high-school, some college, college}
#' \item{hmort}{home mortage (in US $)}
#' \item{hequity}{home equity (in US $)} \item{hval}{home value (in US $)}}
#' @references V. Chernohukov, C. Hansen (2004). The impact of 401(k)
#' participation on the wealth distribution: An instrumental quantile
#' regression analysis. The Review of Economic and Statistics 86 (3), 735--751.
#'
#' @references A. Belloni, V. Chernozhukov, I. Fernandez-Val, and C. Hansen (2014). Program
#' evaluation with high-dimensional data. Working Paper.
#' @keywords datasets 401(k) pension
#' @examples
#' data(pension)
NULL

#' AJR data set
#'
#' Dataset on settler mortality.
#'
#' Data set was analysed in Acemoglu et al. (2001). A detailed description of the data can be found at \url{https://economics.mit.edu/faculty/acemoglu/data/ajr2001} 
#'
#' @name AJR
#' @docType data
#' @format \describe{ 
#' \item{Mort}{Settler mortality}
#' \item{logMort}{logarithm of Mort}
#' \item{Latitude}{Latitude}
#' \item{Latitude2}{Latitude^2}
#' \item{Africa}{Africa}
#' \item{Asia}{Asia}
#' \item{Namer}{North America}
#' \item{Samer}{South America}
#' \item{Neo}{Neo-Europes}
#' \item{GDP}{GDP}
#' \item{Exprop}{Average protection against expropriation risk} 
#' }
#' @references D. Acemoglu, S. Johnson, J. A. Robinson  (2001). Colonial origins of comparative development: an empirical investigation.
#' American Economic Review, 91, 1369--1401.
#' @keywords datasets
#' @examples
#' data(AJR)
NULL

#' Eminent Domain data set
#'
#' Dataset on judicial eminent domain decisions.
#'
#' Data set was analyzed in Belloni et al. (2012).  They estimate the effect of judicial eminent domain decisions on economic outcomes with instrumental variables (IV) in a setting high a large set of potential IVs. 
#' A detailed decription of the data can be found at 
#' \url{https://www.econometricsociety.org/publications/econometrica/2012/11/01/sparse-models-and-methods-optimal-instruments-application} 
#' The data set contains four "sub-data sets" which differ mainly in the dependent variables: repeat-sales FHFA/OFHEO house price index for metro (FHFA) and non-metro (NM) area, the Case-Shiller home price index (CS), 
#' and state-level GDP from the Bureau of Economic Analysis - all transformed with the logarithm. The structure of each subdata set is given above.
#' In the data set the following variables and name conventions are used:
#' "numpanelskx_..." is the number of panels with at least k members with the characteristic following the "_". 
#' The probability controls (names start with "F_prob_") follow a similar naming convention and give the probability of observing a panel with characteristic given following second "_" given the characteristics of the pool of judges available to be assigned to the case. 
#' 
#' Characteristics in the data for the control variables or instruments:
#' \describe{
#' \item{noreligion}{judge reports no religious affiliation}
#' \item{jd_public}{judge's law degree is from a public university}
#' \item{dem}{judge reports being a democrat}
#' \item{female}{judge is female}
#' \item{nonwhite}{judge is nonwhite (and not black)}
#' \item{black}{judge is black}
#' \item{jewish}{judge is  Jewish}
#' \item{catholic}{judge is Catholic}
#' \item{mainline}{baseline religion}
#' \item{protestant}{belongs to a protestant church}
#' \item{evangelical}{belongs to an evangelical church}
#' \item{instate_ba}{judge's undergraduate degree was obtained within state}
#' \item{ba_public}{judge's undergraduate degree was obtained at a public university}
#' \item{elev}{judge was elevated from a district court}
#' \item{year}{year dummy (reference category is one year before the earliest year in the data set (excluded))}
#' \item{circuit}{dummy for the circuit level (reference category excluded)}
#' \item{missing_cy_12}{a dummy for whether there were no cases in that circuit-year}
#' \item{numcasecat_12}{the number of takings appellate decisions}
#'}
#' @name EminentDomain
#' @docType data
#' @format \describe{ 
#' \item{y}{economic outcome variable}
#' \item{x}{set of exogenous variables}
#' \item{d}{eminent domain decisions}
#' \item{z}{set of potential instruments}
#' }
#' @references D. Belloni, D. Chen, V. Chernozhukov and C. Hansen (2012).
#' Sparse models and methods for optimal instruments with an application to
#' eminent domain. \emph{Econometrica} 80 (6), 2369--2429.
#' @keywords datasets
#' @examples
#' data(EminentDomain)
NULL

#' BLP data set
#'
#' Automobile data set from the US.
#'
#' Data set was analysed in Berry, Levinsohn and Pakes (1995). The data stem from annual issues of the Automotive News Market Data Book. 
#' The data set inlcudes information on all models marketed during the the period beginning 1971 and ending in 1990 cotaining 2217 model/years from 997 distinct models.
#' A detailed description is given in BLP (1995, 868--871). The internal function \code{constructIV} constructs instrumental variables along the lines described and used in BLP (1995).
#'
#' @name BLP
#' @docType data
#' @format \describe{ 
#' \item{model.name}{model name}
#' \item{model.id}{model id}
#' \item{firm.id}{firm id}
#' \item{cdid}{cdid}
#' \item{id}{id}
#' \item{price}{log price}
#' \item{mpg}{miles per gallon}
#' \item{mpd}{miles per dollar}
#' \item{hpwt}{horse power per weight}
#' \item{air}{air conditioning (binary variable)}
#' \item{space}{size of the car}
#' \item{share}{market share}
#' \item{outshr}{share s0}
#' \item{y}{outcome variable defined as log(share) - log(outshr) }
#' \item{trend}{time trend}
#' }
#' @references S. Berry, J. Levinsohn, A. Pakes (1995). Automobile Prices in Market EquilibriumD. Econometrica, 63(4), 841--890.
#' @keywords datasets
#' @examples
#' data(BLP)
NULL

#' cps2012 data set
#'
#' Census data from the US for the year 2012.
#'
#' The CPS is a monthly U.S. household survey conducted jointly by the U.S. Census Bureau and the Bureau of Labor Statistics. The data comprise the year 2012. 
#' This data set was used in Mulligan and Rubinstein (2008). 
#' The sample comprises white non-hipanic, ages 25-54, working full time full year (35+ hours per week at least 50 weeks), exclude living in group quarters, 
#'  self-employed, military, agricultural, and private household sector, allocated earning, inconsistent report on earnings and employment, missing data.
#' @name cps2012
#' @docType data
#' @format \describe{ 
#' \item{lnw}{log of hourly wage (annual earnings / annual hours)}
#' \item{female}{female indicator}
#' \item{married status}{ six indicators: widowed, divorced, separated, nevermarried, and married (omitted)}
#' \item{education attainment}{six indicators: hsd08, hsd911, hsg, cg, ad, and sc (omitted)}
#' \item{region indicators}{four indicators: mw, so, we, and ne (omitted)}
#' \item{potential experience}{(max[0, age - years of education - 7]): exp1, exp2 (divided by 100), exp3 (divided by 1000), exp4 (divided by 10000)}
#' \item{weight}{March Supplement sampling weight}
#' \item{year}{CPS year}
#' }
#' @references  C. B. Mulligan and Y. Rubinstein (2008). Selection, investment, and women's relative wages over time. The Quarterly Journal of Economics, 1061--1110.
#' @keywords datasets
#' @examples
#' data(BLP)
NULL