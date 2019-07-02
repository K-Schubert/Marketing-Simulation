func2 <- function(N, d_hyp, d_true, ss, dd){

  #Distribution of d_true
  x1 = rnorm(N[dd], d_true, 1) # expectation = d_true, sd = 1, N samples
  meanx1 = mean(x1)
  
  #Distribution of Control condition  
  x0 = rnorm(N[dd], 0, 1) # expectation = 0, sd = 1, N samples
  meanx0 = mean(x0)
  
  #calculate Standard Error as in R t.test function
  v = ((N[dd]-1)*var(x1)+(N[dd]-1)*var(x0))/(2*N[dd]-2) #not SE but avg deviation between 2 groups (pooled var)
  
  SE =  sqrt(v*(1/N[dd]+1/N[dd]))
  lowdiff = mean(x1 - x0) - 1.96*SE #different way of estimating effect size: more conservative (lower bound of 95% CI) => IGNORE
  
  #Calculate t-values
  Tval[dd,ss,times] = (meanx1-meanx0)/(sqrt((var(x1) + var(x0))/2)*sqrt(2/N[dd]))
  p[dd,ss,times] = 2*pt(-Tval[dd,ss,times], (2*N[dd]-2)) # p-value
  d_emp[dd,ss,times] = 2*Tval[dd,ss,times]/sqrt(N[dd]) #d-emp
  
  TvalLow[dd,ss,times] = (lowdiff)/(sqrt((var(x1) + var(x0))/2)* sqrt(2/N[dd])) # lower bounds (effect based on lower bound) => IGNORE
  d_empLow[dd,ss,times] = 2*TvalLow[dd,ss,times]/sqrt(N[dd]) # lower bounds => IGNORE
  
  L0x[dd,ss,times] = dt(Tval[dd,ss,times], (2*N[dd]-2), 0) # likelihood of 0 given data, density at observed t-value given df and non centrality param (0)
  
  if(!is.na(d_hyp[dd])){
    
    ncp = d_hyp[dd]*sqrt(N[dd])/sqrt(2) #non-centrality param (of t distribution)
    L1x[dd,ss,times] = dt(Tval[dd,ss,times], (2*N[dd]-2), ncp) #likelihood of obs t-value given hypothesied effect           size for transformed distr given ncp
    Lhypx[dd,ss,times] = dt(ncp, (2*N[dd]-2), ncp) #density at median of same distribution (.39)
    LR[dd,ss,times] = L1x[dd,ss,times]/L0x[dd,ss,times] # likelihood ratio (obs t-value and null hyp) (LRT>4 => obs effect size is not from the H0 distr)
    Lhyp[dd,ss,times] = Lhypx[dd,ss,times]/L1x[dd,ss,times] #likelihood ration of the median of hypothesied distr to obs effect size => how far obs effect size from the H0 effect size
  }
}

func1 <- function(d_hyp, dd, alpha, pow, N){
  
  if (!is.na(d_hyp[dd])){
    # N[dd] = round(((qnorm(1-alpha)+qnorm(pow))/(d_hyp[dd])*(sqrt(2)))^2) }
    N[dd] = round(2*(qnorm(1-alpha/2)+qnorm(pow))^2/d_hyp[dd]^2)
  }
  
  # if d_hyp IS NA => if no d_hyp => take random sample size (no hyp in marketing research) to simulate usual way of       consumer researcher
  if (is.na(d_hyp[dd])){
    N[dd] = sample(seq(10,100),1)
  }
}

Lhyp[dd, ss, ] <- replicate(2, func2(N=392, d_hyp=0.2, d_true=0.2, ss=1, dd=1))

replicate(Nsamples, func1(d_hyp=0.2, dd=1, alpha, pow, N))

