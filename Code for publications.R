## ----load libraries-------------------------------------------------
library(tidyverse)
library(pomp)
library(panelPomp)
library(readxl)
library(doParallel)
library(doRNG)
library(readxl)
## ----define rprocess-------------------------------------------------
rproc2 <- Csnippet("
                    double beta,  iota2, foi,  dw, births, S_births, V_births, betaSab, foiSab;
                    double rate[9], trans[38];
             
                    // transmission rate
                    beta = R0*(gamma+mu)*(sigma+mu)/sigma;
                    // force of infection
                    iota2 = rbinom(1, iota/10); #iota is constrained between 0-0.1
                    foi = beta*pow(I+iota2,alpha)/pop;
                    // white noise (extrademographic stochasticity)
                    dw = rgammawn(sigmaSE,dt);
                    
                    
                    // transmission rate
                    betaSab = SabR0*(gamma+mu)*(sigma+mu)/sigma;
                    // force of infection
                    foiSab = betaSab*ISab/pop;
                    
                    rate[0] = foi;  //infection rate (stochastic)
                    rate[1] = mu;             // natural S death
                    rate[2] = foiSab; // second VDPV2 infection rate
                    rate[3] = sigma;        // rate of ending of latent stage
                    rate[4] = mu;             // natural E death
                    rate[5] = gamma;        // recovery
                    rate[6] = mu;             // natural I death
                    rate[7] = zeta; // erlang
                    rate[8] = mu;
                   // Rprintf(rate[7]);
                   // Rprintf(zeta);
                    
                    // Poisson births
                    births = rpois(mu*pop*dt);
                    S_births = rbinom(births,IPV_effxcov);// IPV efficacy and RI IPV cov
                    V_births = births - S_births;
                    
                    // transitions between classes
                    
                    // vaccination
                    if (t==vacc_time1) {
                    double child_vaccS, child_vaccV,child_revS,child_revV,child_RS,child_RV;
                    child_vaccS = rbinom(S,cov_sia * eff_opv2);
                    child_vaccV = rbinom(V,cov_sia * eff_opv2);
                    child_revS = rbinom(child_vaccS,kappa);
                    child_revV = rbinom(child_vaccV,kappa);
                    child_RS = child_vaccS-child_revS;
                    child_RV = child_vaccV-child_revV;
                    reulermultinom(3,S-child_vaccS,&rate[0],dt,&trans[0]);
                    reulermultinom(3,V-child_vaccV,&rate[0],dt,&trans[7]);
                    reulermultinom(1,R+child_RS+child_RV,&rate[6],dt,&trans[10]);
                    reulermultinom(2,ESab + child_revS +child_revV,&rate[3],dt,&trans[23]);
                    reulermultinom(2,E,&rate[3],dt,&trans[3]);
                    S += S_births   - trans[0] - trans[1] - trans[2]-child_vaccS;
                    V += V_births - trans[7] - trans[8]-child_vaccV;
                    R += trans[5] - trans[10] + trans[25] + child_RS + child_RV;
                    ESab += trans[2] - trans[23] - trans[24]+ child_revS+child_revV;
                    E += trans[0] + trans[7] - trans[3] - trans[4];
                    } 
                    else if (t==vacc_time2) {
                    double child_vaccS, child_vaccV,child_revS,child_revV,child_RS,child_RV;
                    child_vaccS = rbinom(S,cov_sia2 * eff_opv2);
                    child_vaccV = rbinom(V,cov_sia2 * eff_opv2);
                    child_revS = rbinom(child_vaccS,kappa);
                    child_revV = rbinom(child_vaccV,kappa);
                    child_RS = child_vaccS-child_revS;
                    child_RV = child_vaccV-child_revV;
                    reulermultinom(3,S-child_vaccS,&rate[0],dt,&trans[0]);
                    reulermultinom(3,V-child_vaccV,&rate[0],dt,&trans[7]);
                    reulermultinom(1,R+child_RS+child_RV,&rate[6],dt,&trans[10]);
                    reulermultinom(2,ESab + child_revS +child_revV,&rate[3],dt,&trans[23]);
                    reulermultinom(2,E,&rate[3],dt,&trans[3]);
                    S += S_births   - trans[0] - trans[1] - trans[2]-child_vaccS;
                    V += V_births - trans[7] - trans[8]-child_vaccV;
                    R += trans[5] - trans[10] + trans[25] + child_RS + child_RV;
                    ESab += trans[2] - trans[23] - trans[24]+ child_revS+child_revV;
                    E += trans[0] + trans[7] - trans[3] - trans[4];
                    } 
                    else {
                    reulermultinom(3,S,&rate[0],dt,&trans[0]);
                    reulermultinom(3,V,&rate[0],dt,&trans[7]);
                    reulermultinom(1,R,&rate[6],dt,&trans[10]);
                    reulermultinom(2,ESab,&rate[3],dt,&trans[23]);
                    reulermultinom(2,E,&rate[3],dt,&trans[3]);
                    S += S_births   - trans[0] - trans[1] - trans[2];
                    V += V_births - trans[7] - trans[8] - trans[9];
                    R += trans[5] - trans[10] + trans[25]; //cehck 25
                    ESab += trans[2] + trans[9] - trans[23] - trans[24];
                    E += trans[0] + trans[7] - trans[3] - trans[4];
                    }
                    
                    reulermultinom(2,I,&rate[5],dt,&trans[5]);
                    
                    reulermultinom(2,H1,&rate[7],dt,&trans[11]);
                    reulermultinom(2,H2,&rate[7],dt,&trans[13]);
                    reulermultinom(2,H3,&rate[7],dt,&trans[15]);
                    reulermultinom(2,H4,&rate[7],dt,&trans[17]);
                    reulermultinom(2,H5,&rate[7],dt,&trans[19]);
                    reulermultinom(2,H6,&rate[7],dt,&trans[21]);
                    
                    reulermultinom(2,ISab,&rate[5],dt,&trans[25]);
                    //
                    
                                      I += trans[3] - trans[5] - trans[6];
                     
                    
                    H1 += trans[0] - trans[11] - trans[12];           // true incidence
                    H2 += trans[11] - trans[13] - trans[14]; 
                    H3 += trans[13] - trans[15] - trans[16]; 
                    H4 += trans[15] - trans[17] - trans[18]; 
                    H5 += trans[17] - trans[19] - trans[20];
                    H6 += trans[19] - trans[21] - trans[22];
                    C += trans[21];
                    
                    
                    ISab += trans[23] - trans[25] - trans[26];
                    W += (dw - dt)/sigmaSE;  // standardized i.i.d. white noise
                    ")


## ----define rinit-------------------------------------------------
rinit2 <- Csnippet("
 double m = pop;
  S = nearbyint(m*(1-IPV_effxcov));
  E = 0;
  I = 0;
  R = 0;
  V = nearbyint(m*IPV_effxcov);
  W = 0;
  C = 0;
  H1 = 0;
  H2 = 0;
  H3 = 0;
  H4 = 0;
  H5 = 0;
  H6 = 0;
  ESab = 0;
  ISab = 0;
                  ")


## ----define ddmeasure-------------------------------------------------
dmeas2 <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1.0e-18;
  if (cases > 0.0) {
      lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)-pnorm(cases-0.5,m,sqrt(v)+tol,1,0)+tol;
  } else {
    lik = pnorm(cases+0.5,m,sqrt(v)+tol,1,0)+tol;
  }
  if (give_log) lik = log(lik);
")
## ----define rmeasure-------------------------------------------------
rmeas2 <- Csnippet("
  double m = rho*C;
  double v = m*(1.0-rho+psi*psi*m);
  double tol = 1.0e-18;
  cases = rnorm(m,sqrt(v)+tol);
  if (cases > 0.0) {
    cases = nearbyint(cases);
  } else {
    cases = 0.0;
  }
")
####---build pomp framework-------------------------------------------------
#Here data_full is the observation dataset, with columns for:
#sub-national location (Province),
#time in days from start of outbreak (time)
#per-day incidence of reported cases (cases)

max_time <- max(data_full$time) #maximum time for model - use time to last observation
pomp(data = data.frame(time=1:max_time, cases = NA),
     t0 = 0,
     time="time",
     rprocess=euler(rproc2,delta.t=1),
     rinit=rinit2,
     dmeasure=dmeas2,
     rmeasure=rmeas2,
     partrans=parameter_trans(
       log=c("mu","sigma","gamma","R0", "sigmaSE", 
             "SabR0", "vacc_time1", "vacc_time2", "psi"),
       logit=c("rho", "kappa", "eff_opv2", "cov_sia", "cov_sia2", "alpha",
               "iota"),
     ),
     accumvars=c("C"),
     statenames=c("S","E","I","R","V","C","H1","H2","H3","H4","H5","H6","W", "ISab", "ESab"),
     paramnames=c("R0","mu","sigma","gamma","alpha", "iota", 
                  "rho","psi","zeta", "sigmaSE", "pop", "SabR0", "vacc_time1", "vacc_time2",
                  "cov_sia", "cov_sia2", "eff_opv2", "kappa", 
                  "IPV_effxcov"
     )
) -> m1

#### ----panelpomp-construction-----------------------------------------------

####Initialize list of pomps for panelpomp
provs <- unique(c(data_full$Province)) #unique subnational units from data
U <- length(provs)
poList <- setNames(vector(mode="list",length=U),
                   nm=paste0("unit",1:U))
for (i.u in seq_len(U)) {
  data_u <- data.frame(subset(data_full,(Province==provs[i.u]),
                              select = c("time", "cases")))
  con_u <- pomp(data_u,
                times="time",
                t0= 0
  )
  poList[[i.u]] <- m1
  poList[[i.u]]@data <- con_u@data
}

## ----starting parameters-----------------------------------------------

params_1 <-c(
  "mu" = 0,
  "sigma"=1/4,
  "gamma"=1/14, # estimated to be shorter in Taj / Congo
  "psi" = 0,
  "rho" = 1/2000,
  "zeta" = 0.329,
  "alpha" = 1,
  "sigmaSE" = 0,
  "SabR0" = 0.9,
  "cov_sia" = 1,
  "cov_sia2" = 1,
  "eff_opv2" = 0.6,
  "kappa" = 1,
  "iota" = 0.05,
  "IPV_effxcov" = 0.6*0.5,
  "R0" = 1.5
)

#specific params matrix

j <- c("population size" = 0, #enter values
       "vacc_time1" = 0,  #enter values 
       "vacc_time2" = 0) #enter values

params_sp_mat <- rbind(
  matrix(j,
         nrow=length(j),
         ncol=length(poList),
         dimnames=list(names(j),names(poList))))

#create pomp
panelPomp(object=poList,
          shared=params_1,
          specific=params_sp_mat) -> pan_pol_vac
####################################


