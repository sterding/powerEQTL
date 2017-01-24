\name{minMAFeQTL.ANOVA}
\alias{minMAFeQTL.ANOVA}
\title{
Calculation of Minimum Detectable Minor Allele Frequency for EQTL Analysis Based on Un-Balanced One-Way ANOVA
}
\description{
Calculation of minimum detectable minor allele frequency (MAF) for eQTL analysis that tests if a SNP is associated to a gene probe by using un-balanced one-way ANOVA.
}
\usage{
minMAFeQTL.ANOVA(effsize, 
                 typeI = 0.05,
                 nTests = 200000,
                 myntotal = 200,
                 mypower = 0.8,
                 verbose = TRUE)
}
\arguments{
  \item{effsize}{
Effect size \eqn{\delta/\sigma}{delta / sigma}.
}
  \item{typeI}{
Type I error rate for testing if a SNP is associated to a gene probe.
}
  \item{nTests}{
integer. Number of tests in eQTL analysis.
}
  \item{myntotal}{
integer. Number of subjects.
}
  \item{mypower}{
Desired power for the eQTL analysis.
}
  \item{verbose}{
logic. indicating if intermediate results should be output.
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

We assume that \eqn{\mu_2-\mu_1=\mu_3-\mu_2=\delta}{mu_2 - mu_1 = mu_3 - mu_2 = delta}, 
where \eqn{\mu_1}{mu_1},
\eqn{\mu_2}{mu_2}, and \eqn{\mu_3}{mu_3} are the
mean gene expression level for mutation homozygotes, heterozygotes,
and wild-type homozygotes, respectively.

Denote \eqn{p}{p} as the minor allele frequency (MAF) of a SNP. Under Hardy-Weinberg equilibrium, we have genotype frequencies: \eqn{p_2=p^2}{p_2 = p^2}, \eqn{p_1=2 p q}{p_1 = 2 * p * q},
and \eqn{p_0=q^2}{p_0 = q^2}, where \eqn{p_2}{p_2}, \eqn{p_1}{p_1}, and \eqn{p_0}{p_0} are genotype
for mutation homozygotes, heterozygotes, and wild-type homozygotes, respectively,
\eqn{q=1-p}{q = 1 - p}. Then ncp can be simplified as
\deqn{ncp=2pq N\left(\frac{\delta}{\sigma}\right)^2,}{ncp = 2 * p * q * N * (delta/sigma)^2.}



}
\value{
minimum detectable MAF.
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
\seealso{
\code{\link{powerEQTL.ANOVA}, \link{powerEQTL.ANOVA2}, \link{ssEQTL.ANOVA}, \link{ssEQTL.ANOVA2}}
}
\examples{

minMAFeQTL.ANOVA(effsize = 1, 
                 typeI = 0.05,
                 nTests = 200000,
                 myntotal = 234,
                 mypower = 0.8,
                 verbose = TRUE)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ method }

