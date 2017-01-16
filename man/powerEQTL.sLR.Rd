\name{powerEQTL.SLR}
\alias{powerEQTL.SLR}

\title{
Power Calculation for EQTL Analysis Based on Simple Linear Regression
}
\description{
Power calculation for eQTL analysis that tests if a SNP is associated to a gene probe by using simple linear regression.
}
\usage{
powerEQTL.SLR(
  MAF,
  typeI = 0.05,
  nTests = 2e+05,
  slope = 0.13,
  myntotal = 200,
  mystddev = 0.13,
  verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{MAF}{
Minor allele frequency.
}
  \item{typeI}{
Type I error rate for testing if a SNP is associated to a gene probe.
}
  \item{nTests}{
integer. Number of tests in eQTL analysis.
}
  \item{slope}{
Slope \eqn{\beta_1}{beta_1} of the simple linear regression
\deqn{y_i = \beta_0+\beta_1 x_i + \epsilon_i,}{y_i = beta_0 + beta_1 * x_i + epsilon_i,}
where \eqn{y_i}{y_i} is the gene expression level of the \eqn{i}{i}-th subject,
\eqn{x_i}{x_i} is the genotype of the \eqn{i}{i}-th subject, and
\eqn{\epsilon_i}{epsilon_i} is the random error term. Additive coding for genotype is used.
}
  \item{myntotal}{
integer. Number of subjects.
}
  \item{mystddev}{
Standard deviation of the random error term \eqn{\epsilon}{epsilon} in simple linear regression.
}
  \item{verbose}{
logic. indicating if intermediate results should be output.
}
}
\details{
To test if a SNP is associated with a gene probe, we use the simple linear regression
\deqn{y_i = \beta_0+\beta_1 x_i + \epsilon_i,}{y_i = beta_0 + beta_1 * x_i + epsilon_i,}
where \eqn{y_i}{y_i} is the gene expression level of the \eqn{i}{i}-th subject,
\eqn{x_i}{x_i} is the genotype of the \eqn{i}{i}-th subject, and
\eqn{\epsilon_i} is the random error term. Additive coding for genotype is used. To test if the SNP is associated with the gene probe, we test the null hypothesis \eqn{H_0: \beta_1=0}{H_0: beta_1 = 0}.

Denote \eqn{p}{p} as the minor allele frequency (MAF) of the SNP. Under Hardy-Weinberg equilibrium, we can calculate the variance of genotype of the SNP: 
\eqn{\sigma^2_x=2 p (1-p)}{sigma^2_x = 2 * p * (1 - p)}, 
where \eqn{\sigma^2_x}{sigma^2_x} is the variance of the
predictor (i.e. the SNP) \eqn{x_i}{x_i}.

We then can use Dupont and Plummer's (1998) power/sample size calculation formula to calculate the minimum detectable slope, adjusting for multiple testing.
}
\value{
power of the test after Bonferroni correction for multiple testing.
}
\references{
Dupont, W.D. and Plummer, W.D.. Power and Sample Size Calculations for Studies Involving Linear Regression. Controlled Clinical Trials. 1998;19:589-601.
}
\author{
Xianjun Dong <XDONG@rics.bwh.harvard.edu>,
Tzuu-Wang Chang <Chang.Tzuu-Wang@mgh.harvard.edu>,
Scott T. Weiss <restw@channing.harvard.edu>,
Weiliang Qiu <stwxq@channing.harvard.edu>
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ssEQTL.SLR}, \link{minSlopeEQTL.SLR}}
}
\examples{

powerEQTL.SLR(
  MAF = 0.1,
  typeI = 0.05,
  nTests = 2e+05,
  slope = 0.13,
  myntotal = 176,
  mystddev = 0.13,
  verbose = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

