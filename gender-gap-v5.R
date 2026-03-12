####################

# Homework: light code modifications + comments

####################

# This empirical example uses the data from the CPS to illustrate the use of distribution regression and counterfactual distributions
# in the estimation of the gender gap

# Authors: V. Chernozhukov and I. Fernandez-Val

# Data source: IPUMS-CPS through Paul Schrimpf
# URL for data: https://cps.ipums.org/cps/

# Description of the data: the sample selection and variable contruction follow

# Mulligan, Casey B., and Yona Rubinstein. "Selection, investment, and women's relative wages over time." The Quarterly Journal of Economics (2008): 1061-1110.

# Sample selection: white non-hipanic, ages 25-54, working full time full year (35+ hours per week at least 50 weeks), exclude living in group quarters, 
# self-employed, military, agricultural, and private household sector, allocated earning, inconsistent report on earnings and employment, missing data

# The variables used in the analysis include:

# - lnw: log of hourly wage (annual earnings / annual hours)
# - female: female indicator
# - 6 married status indicators: widowed, divorced, separated, nevermarried, and married (omitted)
# - 6 education attainment indicators: hsd08, hsd911, hsg, cg, ad, and sc (omitted)
# - 4 region indicators: mw, so, we, and ne (omitted)
# - Quartic in potential experience (max[0, age - years of education - 7]): exp1, exp2 (divided by 100), exp3 (divided by 1000), exp4 (divided by 10000)
# - weight: March Supplement sampling weight
# - year: CPS year

# Updated on 04/12/2021

# This versions uses simulataneous confidence bands for all the distribution functions

# We use the package discreteQ available at https://github.com/bmelly/discreteQ

# clear environment
rm(list = ls());

# set path
filepath = "~/Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology/Documents/2 Class/14.382/14.382_pset3"
setwd(filepath)

# required packages and data
options(warn=-1); # sets warnings off;
library(discreteQ);
library(xtable);

load('cps2015.rdata');
attach(data);
 
########## Descriptive statistics
# set up the table
vars <- c("lnw","female","married","widowed","separated","divorced","nevermarried","lhs","hsg","sc","cg","ad","ne","mw","so","we","exp1");
options(digits=2);
dstats <- cbind(sapply(data[,vars], weighted.mean, weight), apply(data[female==1,vars], 2, weighted.mean, weight[female==1]), apply(data[female==0,vars], 2, weighted.mean, weight[female==0]));
# print summary stats table
print(xtable(dstats),file = "gender_gap_summary_stats.tex", include.rownames = FALSE);



# Parameters for estimation;
nind  <- 39; # Number of distribution regressions to estimate conditional and counterfactual distribution;
taus <- c(3:(nind-2))/(nind+1); # quantile indexes of interest, compute as [0,1] equally spaced gridpoints
alpha <- 0.10; # significance level for confidence bands

######### Estimation

# formula specification
form  <- lnw ~ widowed + divorced + separated + nevermarried +(lhs+hsg+cg+ad)*(exp1+exp2+exp3+exp4) + mw + so + we;
reg <- model.matrix(form, data=data); # matrix of regressors for DR
ys <- quantile(lnw, c(.02,c(1:nind)/(nind+1), .98)); # thresholds for DR estimation as the quantiles of the wage distribution
R <- 200; #number of bootstrap repetitions
my_cl <- parallel::makePSOCKcluster(24, setup_timeout = 0.5); # number of proccessors for parallel computing
set.seed(14382);

# Estimation (take a while to run);
dr.fit <- discreteQ(y=lnw, d=1-female, x=reg, w=weight, decomposition = TRUE, q.range=range(taus),
                 method="logit", bsrep=R, alpha=alpha, ys=ys, cl= my_cl);
# Define the treatment (d) as male to obtain positive values for the gender gap


########## Figure 1: wage distributions
data |>
  filter(female == 1) |> 
  ggplot(aes(x=lnw)) +
  geom_histogram(fill = "purple", color = "white") +
  theme_minimal() +
  labs(
    x = "Log of hourly wage",
    y = "Count",
    title = "Women"
  )
ggsave(
  "women_hist.png",
  width = 3.2,
  height = 2.2,
  units = "in",
  dpi = 300
)

data |>
  filter(female == 0) |> 
  ggplot(aes(x=lnw)) +
  geom_histogram(fill = "green", color = "white") +
  theme_minimal() +
  labs(
    x = "Log of hourly wage",
    y = "Count",
    title = "Men"
  )
ggsave(
  "women_hist.png",
  width = 3.2,
  height = 2.2,
  units = "in",
  dpi = 3
)

############## Figure 2: counterfactual distributions
png("gender_gap_dist.png", width = 8, height = 6, units = "in", res = 300)

# set up plot with the x axis being the wage range
plot(range(ys), c(0,1), type="n",col=1, xlab="Log of hourly wage", lwd=2,
     ylab="Probability", 
     main="Observed and Counterfactual Distributions (90% CI) ",
     sub=" ")

# plot lower and upper bounds for the female distribution (F0)
polygon(
  c(ys,rev(ys)),
  c(dr.fit$ub.F0(ys), rev(dr.fit$lb.F0(ys))), 
  density=100, 
  order=F, 
  col='light blue', 
  lty = 1, 
  lwd = 1)
lines(ys, dr.fit$F0(ys), lwd=1, col = 'dark green'  );

# plot lower and upper bounds for the counterfactual distribution (Fc)
polygon(c(ys,rev(ys)),c(dr.fit$ub.Fc(ys),rev(dr.fit$lb.Fc(ys))), density=100, border=F, 
        col='light grey', lty = 1, lwd = 1);
lines(ys, dr.fit$Fc(ys), lwd=1, col = 'dark grey'  );

# plot lower and upper bounds for the male distribution (F1)
polygon(c(ys,rev(ys)),c(dr.fit$ub.F1(ys),rev(dr.fit$lb.F1(ys))), density=100, border=F, 
        col='light green', lty = 1, lwd = 1);
lines(ys, dr.fit$F1(ys), lwd=1, col = 4  )

legend(quantile(ys, .01), 1, c(' ', ' ',' ' ), col = c('light blue','light green','light grey'), lwd = c(4,4,4), horiz = F, bty = 'n');
legend(quantile(ys, .01), 1, c('Observed women distribution','Observed men distribution', 'Counterfactual distribution'), col = c(4,'dark green', 'dark grey'), lwd = c(1,1,1), horiz = F, bty = 'n');


dev.off();

############## Figure 3: Gender gap by quantile
png("gender_gap_quantiles.png", width = 8, height = 6, units = "in", res = 300)

par(mfrow=c(1,1));

# set up plot
plot(range(taus), range(ys),  type="n",col=1, xlab="Quantile index", lwd=2,
     ylab="Log of hourly wage", 
     main="Observed and Counterfactual Quantiles (90% CI) ",
     sub=" ");

# plot upper bounds by quantile
polygon(c(taus,rev(taus)),c(dr.fit$ub.Q0(taus),rev(dr.fit$lb.Q0(taus))), density=100, border=F, 
        col='light blue', lty = 1, lwd = 1);
lines(taus, dr.fit$Q0(taus), lwd=1, col = 4  );

# plot center by quantile
polygon(c(taus,rev(taus)),c(dr.fit$ub.Qc(taus),rev(dr.fit$lb.Qc(taus))), density=100, border=F, 
        col='light grey', lty = 1, lwd = 1);
lines(taus, dr.fit$Qc(taus), lwd=1, col = 'dark grey'  );

# plot lower bounds by quantile
polygon(c(taus,rev(taus)),c(dr.fit$ub.Q1(taus),rev(dr.fit$lb.Q1(taus))), density=100, border=F, 
        col='light green', lty = 1, lwd = 1);
lines(taus, dr.fit$Q1(taus), lwd=1, col = 'dark green'  );


legend(min(taus), max(ys), c(' ', ' ',' '), col = c('light blue','light green','light grey'), lwd = c(4,4,4), horiz = F, bty = 'n');
legend(min(taus), max(ys), c('Observed women quantiles', 'Observed men quantiles', 'Counterfactual quantiles'), col = c(4,'dark green','dark grey'), lwd = c(1,1,1), horiz = F, bty = 'n');


dev.off();


############## Figure 4: Decomposing the gender gap
png("gender_gap_decomposition.png", width = 8, height = 6, units = "in", res = 300)

par(mfrow=c(3,1));

plot(dr.fit, which="observed", support="continuous",main="Gender wage gap (90% CI)", ylim = c(-.15,.5));
abline(h=0, col ="light grey")
plot(dr.fit, which="unexplained", support="continuous",main="Discrimination (90% CI)", ylim = c(-.15,.5));
abline(h=0, col ="light grey")
plot(dr.fit, which="composition", support="continuous",main="Composition (90% CI)", ylim = c(-.15,.5));
abline(h=0, col ="light grey")

dev.off();





