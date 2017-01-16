\name{ssEQTL.ANOVA}
\alias{ssEQTL.ANOVA}
\title{
Sample Size Calculation for EQTL Analysis Based on Un-Balanced One-Way ANOVA
}
\description{
Sample size calculation for eQTL analysis that tests if a SNP is associated to a gene probe by using un-balanced one-way ANOVA.
}
\usage{
ssEQTL.ANOVA(
  MAF,
  typeI = 0.05,
  nTests = 2e+05,
  mypower = 0.8,
  mystddev = 0.13,
  deltaVec = c(0.13, 0.13))
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
  \item{mypower}{
Desired power for the eQTL analysis.
}
  \item{mystddev}{
Standard deviation of gene expression levels in one group of subjects.
Assume all 3 groups of subjects (mutation homozygote, heterozygote, wild-type homozygote) have the same standard deviation of gene expression levels.
}
  \item{deltaVec}{
A vector having 2 elements. The first element is equal to
\eqn{\mu_2-\mu_1}{mu_2 - mu_1} and the second elementis equalt to
\eqn{\mu_3-\mu_2}{mu_3 - mu_2}, where \eqn{\mu_1}{mu_1} is
the mean gene expression level for the mutation homozygotes,
\eqn{\mu_2}{mu_2} is the mean gene expression level for the heterozygotes,
and \eqn{\mu_3}{mu_3} is the mean gene expression level for the 
wild-type gene expression level.
}
}
\details{
The assumption of the ANOVA approach is that the association of a SNP to a gene probe is tested by using un-balanced one-way ANOVA (e.g. Lonsdale et al. 2013).
According to SAS online document
\url{https://support.sas.com/documentation/cdl/en/statug/63033/HTML/
default/viewer.htm#statug_power_a0000000982.htm}, the power calculation formula is
\deqn{power=Pr\left(\left.F\geq F_{1-\alpha}\left(k-1, N-k\right)\right| F\sim
F_{k-1, N-k, \lambda}\right),}{power = Pr(F >= F(1 - alpha, k - 1, N - k)| F ~ F(k - 1, N - k, lambda)),}
where \eqn{k=3}{k = 3} is the number of groups of subjects, \eqn{N}{N} is the total number
of subjects, \eqn{F_{1-\alpha}\left(k-1, N-k\right)}{F_{1 - alpha}(k - 1, N - k)} is the 
\eqn{100(1-\alpha)}{100 * (1 - alpha)}-th percentile of central F distribution with degrees of freedoms \eqn{k-1}{k - 1} and \eqn{N-k}{N - k}, and \eqn{F_{k-1, N-k, \lambda}}{F_{k - 1, N - k, lambda}} is the non-central F distribution
with degrees of freedoms \eqn{k-1}{k - 1} and \eqn{N-k}{N - k} and non-central parameter (ncp)
\eqn{\lambda}{lambda}. The ncp \eqn{\lambda}{lambda} is equal to
\deqn{\lambda=\frac{N}{\sigma^2}\sum_{i=1}^{k} w_i \left(\mu_i-\mu\right)^2,}{
lambda = N * sum(wi * (mu_i - mu)^2, i = 1,.., k)/sigma^2,
} where \eqn{\mu_i}{mu_i} is the mean gene expression level
for the \eqn{i}{i}-th group of subjects, \eqn{w_i}{w_i} is the weight for the \eqn{i}{i}-th group of subjects, \eqn{\sigma^2}{sigma^2} is the variance of the random errors in ANOVA (assuming each group has equal variance), and
\eqn{\mu}{mu} is the weighted mean gene expression level
\deqn{\mu=\sum_{i=1}^{k}w_i \mu_i.}{mu = sum(w_i * mu_i, i = 1, ..., k).}
The weights \eqn{w_i}{w_i} are the sample proportions for the 3 groups of subjects. Hence, 
\eqn{\sum_{i=1}^{3}w_i = 1}{sum(w_i, i = 1, 2, 3) = 1}.

}
\value{
sample size required for the eQTL analysis to achieve the desired power.
}
\references{
Lonsdale J and Thomas J, et al. The Genotype-Tissue Expression (GTEx) project. Nature Genetics, 45:580-585, 2013.
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
\code{\link{minEffectEQTL.ANOVA}, \link{powerEQTL.ANOVA}, \link{powerEQTL.ANOVA2}, \link{ssEQTL.ANOVA2}}
}
\examples{
ssEQTL.ANOVA(MAF = 0.1,
       typeI = 0.05,
       nTests = 200000,
       mypower = 0.8,
       mystddev = 0.13,
       deltaVec = c(0.13, 0.13))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

