## Bidirectional vs. Unidirectional IV-CLPM, Given Unidirectional Causal Effects

## Aim: Compare the LRT statistics in bidirectional and unidirectional models

# clear workspace
rm(list=ls(all=TRUE))

# SET THE WORKING DIRECTORY
# getwd()
resfile <- 'resfile_bidir_v_unidir_IV-CLPM.dta'  # output file for the results

# PART I. MODEL FITTING  ---------------------------------------------------------

# a function for power calc based on the central and noncentral chisq distributions
getchipow <- function(alpha,df,Tval) {
  ca <- qchisq(alpha,df,ncp=0,lower.tail=F)
  # critical value given alpha
  power <- pchisq(ca,df,ncp=Tval,lower.tail=F)
  power
}

# Libraries
library(MASS)
Sys.setenv(OMP_NUM_THREADS = parallel::detectCores())
library(OpenMx)

# Sample size
N <- 1000

# Alpha
alpha <- 0.05

# Number of occasions in simulated time-series
T_ <- 150

# Number of variables (except the IVs)
nv <- T_*2     # X and Y * T_

# Check ID?  
doID <- T

# IVx --> X_t
px <- .05

# IVy --> Y_t
py <- .05

# r(IVx,IVy)
rIV <- 0.25      # Correlation of the 2 IVs

# X_t --> Y_{t+1}
bxy <- .4     # 1st-order Causal path coef.

# Y_t --> X_{t+1}
byx <- 0     # 1st-order Causal path coef. 

# Residual Cross-sectional Correlation of X and Y
r_xy <- .3   

# AR1
bx  <- .8      # AR1 of X 
by  <- .8      # AR1 of Y 


## CLPM Panel Time Lag
startT <- 100              # Wave 1 across all models = T100
Tdis <- seq(1:50)          # Wave 2 changed from T101 to T150
# 50 models will be fitted to each time-series generated using the parameters specified above


## Number of cells in the factorial design
nrows <- length(Tdis)
nrows


# A matrix to store the parameters of interest
keep <- matrix(0,nrows,42)
# 42 parameters
colnames(keep) <-  c(
  ## Data-generating parameters
  'bxy','byx','bxx','byy','rxy','Tdis',
  ## LRT Statistics
  ## Unidirectional Effect
  'NCP_distal_xy','NCP_prox_xy1','NCP_prox_xy2',
  'NCP_prox12_xy','NCP_dist_prox1_xy','NCP_dist_prox2_xy','NCP_all_xy',
  ## LRTs for the bidirectional effects
  'NCP_distal_bd','NCP_prox1_bd','NCP_prox2_bd',
  'NCP_prox12_bd','NCP_dist_prox1_bd','NCP_dist_prox2_bd','NCP_all_bd',
  ## R-sq - to check that the simulated variables have sensible relationship
  'IVx->X1','IVx->X2','IVy->Y1','IVy->Y2',
  'X1->Y1','X2->Y2','Y1->X1','Y2->X2','X1->Y2','Y1->X2',
  # ## Recovered model estimates
  'b_distal_xy','b_prox1_xy','b_prox2_xy','b_distal_yx','b_prox1_yx','b_prox2_yx',
  ## Standardized path coefficients
  'std_distal_xy','std_prox1_xy','std_prox2_xy','std_distal_yx','std_prox1_yx','std_prox2_yx'
)


# --------------------------------------------------------------------------------------------


# start loop ---------------------------------------------------------------------------------

# start time
time_start <- Sys.time()

irow <- 0

for (idis in Tdis) {       # Panel Lag
  
  irow  <-  irow+1
  
  # Data-generating parameters
  keep[irow,1:5] <-  c(bxy, byx, bx, by, r_xy)
  
  ## Data simulation ----------------------------------------------------
  
  vx <- 1               # Residual variance of X
  vy <-  1              # Residual Variance of Y
  vxy <- r_xy           # Correlation of the Residuals
  
  # Beta Matrix (except the column for the IV)
  # X   YX
  # XY  Y
  BE_ <-   matrix(0,nv,nv)        # Full matrix initiated as zeroes
  BE_xy <- matrix(0,T_,T_)        # Lower Off-diagonal block: X-->Y
  BE_yx <- matrix(0,T_,T_)        # Upper Off-diagonal block: Y-->X
  BE_x <-  matrix(0,T_,T_)        # Top diagonal block for AR of X
  BE_y <-  matrix(0,T_,T_)        # Bottom diagonal block for AR of Y
  
  #
  for (i in 2:T_) { BE_x[i,i-1] <- bx }    # First order AR of X
  for (i in 2:T_) { BE_y[i,i-1] <- by }    # First order AR of Y
  for (i in 2:T_) { BE_xy[i,i-1] <- bxy }  # First order causal effect of X on Y
  for (i in 2:T_) { BE_yx[i,i-1] <- byx }  # First order causal effect of Y on X
  
  # assemble
  BE_[1:T_, 1:T_] <- BE_x                      # ARx  [1:150, 1:150]
  BE_[(1+T_):(2*T_), 1:T_] <- BE_xy            # X->Y [151:300, 1:150]
  BE_[1:T_, (1+T_):(2*T_)] <- BE_yx            # Y->X [1:150, 151:300]
  BE_[(1+T_):(2*T_), (1+T_):(2*T_)] <- BE_y    # ARy  [151:300, 151:300]
  
  # check
  # BE_[1:5,1:5]
  # BE_[151:155,1:5]
  # BE_[1:5,151:155]
  # BE_[151:155,151:155]
  
  ## Additional columns for IVx-->X, IVy-->Y
  Gamma <- matrix(0,nv+2,2)                     
  Gamma[1:(nv+2),1] <- c(0,0,rep(px,T_),rep(0,T_))  # IVx -> X, rows 3:152 are X1:X150                    
  Gamma[1:(nv+2),2] <- c(0,0,rep(0,T_),rep(py,T_))  # IVy -> Y, rows 153:302 are Y1:Y150
  
  # check
  # Gamma[1:6,]
  # Gamma[150:155,]
  
  # Bind the gamma columns to the beta matrix
  BE_ <- cbind(Gamma[3:(nv+2),], BE_)  
  BE_ <- rbind(rep(0,(nv+2)), rep(0,(nv+2)), BE_)      # Row 1 and 2 are for the IVs
  
  # check
  # BE_[1:7,1:7]
  
  # Psi Matrix - construct similarly as with Beta
  PS_ <- matrix(0,nv,nv)
  PS_x <- matrix(0,T_,T_)
  PS_y <- matrix(0,T_,T_)
  PS_xy <- matrix(0,T_,T_)
  
  #
  PS_x[1,1] <- vx    
  PS_y[1,1] <- vy
  PS_xy[1,1] <- vxy
  for (i in 2:T_) { PS_x[i,i] <- vx }     # Residual Variance of X
  for (i in 2:T_) { PS_y[i,i] <- vy }     # Residual Variance of Y
  for (i in 2:T_) { PS_xy[i,i] <- vxy }   # Residual Covariance of X and Y
  # assemble
  PS_[1:T_, 1:T_] <- PS_x
  PS_[(1+T_):(2*T_), 1:T_] <- PS_xy
  PS_[1:T_,(1+T_):(2*T_)] <- PS_xy
  PS_[(1+T_):(2*T_), (1+T_):(2*T_)] <- PS_y
  # Add additional rows and columns the IVs 
  PS_ <- cbind(c(rep(0,(nv))), c(rep(0,(nv))), PS_)
  PS_ <- rbind(c(rep(0,(nv+2))), c(rep(0,(nv+2))), PS_)
  
  # variance of the IVs
  PS_[1,1] <- PS_[2,2] <- 1 
  # covariance of the IVs
  PS_[2,1] <- PS_[1,2] <- rIV
  
  # check
  # PS_[1:6,1:6]
  # PS_[150:155, 1:6]
  
  # Matrices for the algebra
  I_ <- diag(nv+2)
  IBE_ <- solve(I_-BE_)
  
  # Expected Covarianc Matrix
  S_ <- IBE_%*%PS_%*%t(IBE_)
  rownames(S_) <- colnames(S_) <- c('IVx', 'IVy', paste('X',1:T_, sep=""),paste('Y',(1:T_), sep=""))
  
  # check
  # round(cov2cor(S_[1:6, 1:6]),3)
  # round(cov2cor(S_[153:157, 153:157]),3)
  # round(cov2cor(S_[c(1,2,102:105,252:255),c(1,2,102:105,252:255)]),3)
  
  ## Model Fitting -----------------------------------------------------
  
  # The two panels of CLPM
  Ts <- c(startT, startT+idis)
  keep[irow,6] <- idis
  
  # Variables for CLPM:
  isel <- c(1, 2,                        # IVx, IVy
            Ts[1]+2 , Ts[1]+T_+2,      # X and Y at Wave 1
            Ts[2]+2, Ts[2]+T_+2)       # X and Y at Wave 2
  Smod <- S_[isel, isel]
  
  # check covariances 
  if(bxy==max(bxys) & byx==max(byxs) & bx==max(bxs) & by==max(bys) & r_xy==min(rxys)) {
    # covariances at the closest wave2
    if (idis == min(Tdis)) {
      MinDis <- paste("Shortest panel lag: Time lag =",idis)
      MinCov <- round(Smod,3)
      MinCor <- round(cov2cor(Smod),3)
    }
    # covariances at the farthest wave
    if (idis == max(Tdis)) {
      MaxDis <- paste("Longest panel lag: Time lag =",idis)
      MaxCov <- round(Smod,3)
      MaxCor <- round(cov2cor(Smod),3)
    }
  }
  
  #
  # Number of observed variables (ny) = Number of underlying latent scores (ne)
  ny <- 6
  ne <- 6  
  #
  # Means
  mu_mod <- rep(0,ny)
  # Data generation
  dat_mod <- mvrnorm(N, mu_mod, Sigma=Smod*(N/(N-1)), emp=T)  
  dat_mod <- as.data.frame(dat_mod) 
  vnames <- colnames(dat_mod) <-  c('IVx', 'IVy', paste(c('X','Y'),Ts[1],sep=""), paste(c('X','Y'),Ts[2],sep=""))
  vnames1 <- vnames
  
  # Matrices for mxModel
  
  mxfreeBE <- matrix(c(
    # px  py  x1  y1  x2  y2   #
    F,  F,  F,  F,  F,  F,  # px
    F,  F,  F,  F,  F,  F,  # py
    T,  F,  F,  T,  F,  F,  # x1
    F,  T,  T,  F,  F,  F,  # y1
    T,  F,  T,  T,  F,  T,  # x2
    F,  T,  T,  T,  T,  F   # y2
  ),ne,ne, byrow=T)
  #
  mxlabelsBE <- matrix(c(
    NA,   NA,  NA,     NA,     NA,   NA,
    NA,   NA,  NA,     NA,     NA,   NA,
    'px1',NA,  NA,     'byx1', NA,   NA,
    NA,  'py1','bxy1', NA,     NA,   NA,
    'px2',NA,  'bxx',  'byx12',NA,   'byx2', 
    NA,  'py2','bxy12','byy',  'bxy2',NA
  ),ne,ne, byrow=T)   
  #
  mxvaluesBE <- matrix(c(
    # px  py   x1   y1   x2   y2   #
    0.00,0.00,0.00,0.00,0.00,0.00, # px
    0.00,0.00,0.00,0.00,0.00,0.00, # py
    0.05,0.00,0.00,0.09,0.00,0.00, # x1
    0.00,0.05,0.09,0.00,0.00,0.00, # y1
    0.05,0.00,0.60,0.04,0.00,0.06, # x2
    0.00,0.05,0.04,0.60,0.06,0.00  # y2
  ),ne,ne, byrow=T)
  # starting values
  #
  mxfreePS <- matrix(c(
    # px  py   x1   y1   x2   y2   #
    T,T,F,F,F,F,  # px
    T,T,F,F,F,F,  # py
    F,F,T,T,F,F,  # x1
    F,F,T,T,F,F,  # y1
    F,F,F,F,T,T,  # x2
    F,F,F,F,T,T   # y2
  ), ne,ne, byrow=T)
  mxlabelsPS <- matrix(c(
    'Vpx','Cxyp',NA,NA,NA,NA,
    'Cxyp','Vpy',NA,NA,NA,NA,
    NA,NA,'Vx1','Cxy1',NA,NA,
    NA,NA,'Cxy1','Vy1',NA,NA,
    NA,NA,NA,NA,'Vx2','Cxy2',
    NA,NA,NA,NA,'Cxy2','Vy2'
  ), ne,ne, byrow=T)
  mxvaluesPS <- matrix(c(
    1.00, 0.0, 0.00, 0.00, 0.00, 0.00,
    0.0, 1.00, 0.00, 0.00, 0.00, 0.00,
    0.00, 0.00, 0.90, 0.30, 0.00, 0.00,
    0.00, 0.00, 0.30, 0.90, 0.00, 0.00,
    0.00, 0.00, 0.00, 0.00, 0.90, 0.30,
    0.00, 0.00, 0.00, 0.00, 0.30, 0.90
  ), ne,ne, byrow=T)
  
  # Lambda matrix connecting observed to latent score 
  # (loadings fixed at 1, assuming no measurement error)
  mxvaluesLY <- matrix(0,ny,ne)
  mxvaluesLY[1:ny,1:ny] <- diag(ny)
  mxfreeLY <- matrix(F,ny,ne)
  mxlabelsLY <- matrix(NA,ny,ne)
  
  #
  # Regression models... effect sizes.
  #
  dat1 <- mvrnorm(N,rep(0,6),Sigma=Smod,emp=T)
  as.data.frame(dat1) -> dat1
  colnames(dat1)  <-  c('IVx','IVy','X1','Y1','X2','Y2')
  
  r2s <- matrix(0,10,1)
  r1 <- summary(lm(X1~IVx,data=dat1))$r.squared   # IVx --> X1
  r2 <- summary(lm(X2~IVx,data=dat1))$r.squared   # IVx --> X2
  r3 <- summary(lm(Y1~IVy,data=dat1))$r.squared   # IVy --> Y1
  r4 <- summary(lm(Y2~IVy,data=dat1))$r.squared   # IVy --> Y2
  r5 <- summary(lm(Y1~X1,data=dat1))$r.squared   # X1 --> Y1
  r6 <- summary(lm(Y2~X2,data=dat1))$r.squared   # X2 --> Y2
  r7 <- summary(lm(X1~Y1,data=dat1))$r.squared   # Y1 --> X1
  r8 <- summary(lm(X2~Y2,data=dat1))$r.squared   # Y2 --> X2
  r9 <- summary(lm(Y2~X1,data=dat1))$r.squared   # X1 --> Y2
  r10 <- summary(lm(X2~Y1,data=dat1))$r.squared   # Y1 --> X2
  r2s[,1] <- c(r1,r2,r3,r4,r5,r6,r7,r8,r9,r10)
  rownames(r2s) <-  c(  'IVx->X1','IVx->X2','IVy->Y1','IVy->Y2',
                        'X1->Y1','X2->Y2','Y1->X1','Y2->X2','X1->Y2','Y1->X2')
  keep[irow,21:30] <- r2s
  #
  # 
  # Means
  mu_mod <- matrix(0,1,ny)
  colnames(mu_mod) <- vnames1
  rownames(Smod) <- colnames(Smod) <- vnames1
  
  # Full Model
  MRmodel1 <- mxModel('MR',
                      mxMatrix(type='Full',nrow=ne,ncol=ne,
                               free=mxfreeBE,value=mxvaluesBE,labels=mxlabelsBE,dimnames=list(vnames1,vnames1),name='BE'),
                      # Beta
                      mxMatrix(type='Iden',nrow=ne,ncol=ne,name='I'),
                      #
                      mxMatrix(type='Full',nrow=ny,ncol=ne,
                               free=mxfreeLY,value=mxvaluesLY,labels=mxlabelsLY,dimnames=list(vnames1,vnames1),name='LY'),
                      # Lambda
                      mxMatrix(type='Symm',nrow=ne,ncol=ne,
                               free=mxfreePS,value=mxvaluesPS,labels=mxlabelsPS,dimnames=list(vnames1,vnames1),name='PS'),
                      # Psi
                      mxMatrix(type='Full',nrow=1,ncol=ny,
                               free=T,value=0,
                               labels=c('meIVx','meIVy','meX1','meY1','meX2','meY2'),name='mu'),
                      # Means
                      mxAlgebra(expression=solve(I-BE),name='iBE'),
                      mxAlgebra(expression=LY%*%iBE,  name='LYB'),
                      #
                      mxAlgebra(expression=LYB%*%(PS)%*%t(LYB),name='Sigma')     
  )
  #
  expF <- mxExpectationNormal(covariance="MR.Sigma", means="MR.mu", dimnames=vnames1)
  #
  #
  fitFunction <- mxFitFunctionML()
  dat_  <-  mxData( observed=Smod, type="cov", means=mu_mod, numObs=N)
  Model_MR <-  mxModel("MR1", MRmodel1,expF, fitFunction, dat_)
  #
  MRModel1out <- mxTryHard(Model_MR, 50)
  #
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel1out, details=T)
    test$non_identified_parameters
  }
  #
  # summary(MRModel1out)
  
  ## LRT ---------------------------------------------------------------
  
  ### Distal X>Y -----------------------------------------------------------
  Model_MR2  <-  omxSetParameters(MRModel1out, labels=c('bxy12'), values=c(0), free=c(F))
  MRModel2out <- mxTryHard(Model_MR2,50)
  mxCompare(MRModel1out, MRModel2out)
  T1  <- as.numeric(mxCompare(MRModel1out, MRModel2out)[2,7:8])
  if (T1[1]<.000001)  {T1[1] = 0}   # ncp = 0
  keep[irow,7] <- T1[1]
  getchipow(alpha,T1[2],T1[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel2out, details=T)
    test$non_identified_parameters
  }
  
  ### Proximal X>Y @T1 ------------------------------------------------------
  Model_MR3  <-  omxSetParameters(MRModel1out, labels=c('bxy1'), values=c(0), free=c(F))
  MRModel3out <- mxTryHard(Model_MR3,50)
  mxCompare(MRModel1out, MRModel3out)
  T2  <- as.numeric(mxCompare(MRModel1out, MRModel3out)[2,7:8])
  if (T2[1]<.000001)  {T2[1] = 0}   # ncp = 0
  keep[irow,8] <- T2[1]
  getchipow(alpha,T2[2],T2[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel3out, details=T)
    test$non_identified_parameters
  }
  
  ### Proximal X>Y @T2 ------------------------------------------------------
  Model_MR4  <-  omxSetParameters(MRModel1out, labels=c('bxy2'), values=c(0), free=c(F))
  MRModel4out <- mxTryHard(Model_MR4,50)
  mxCompare(MRModel1out, MRModel4out)
  T3  <- as.numeric(mxCompare(MRModel1out, MRModel4out)[2,7:8])
  if (T3[1]<.000001)  {T3[1] = 0}   # ncp = 0
  keep[irow,9] <- T3[1]
  getchipow(alpha,T3[2],T3[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel4out, details=T)
    test$non_identified_parameters
  }
  
  
  
  ### All Three X>Y -----------------------------------------------------------
  Model_MR5  <-  omxSetParameters(MRModel1out, labels=c('bxy1','bxy2','bxy12'), values=c(0), free=c(F))
  MRModel5out <- mxTryHard(Model_MR5, 50)
  mxCompare(MRModel1out, MRModel5out)
  T4  <- as.numeric(mxCompare(MRModel1out, MRModel5out)[2,7:8])
  if (T4[1]<.000001)  {T4[1] = 0}   # ncp = 0
  keep[irow,13] <- T4[1]
  getchipow(alpha,T4[2],T4[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel5out, details=T)
    test$non_identified_parameters
  }
  
  
  ### Distal Bidir -----------------------------------------------------------
  Model_MR6  <-  omxSetParameters(MRModel1out, labels=c('bxy12','byx12'), values=c(0), free=c(F))
  MRModel6out <- mxTryHard(Model_MR6,50)
  mxCompare(MRModel1out, MRModel6out)
  T5  <- as.numeric(mxCompare(MRModel1out, MRModel6out)[2,7:8])
  if (T5[1]<.000001)  {T5[1] = 0}   # ncp = 0
  keep[irow,14] <- T5[1]
  getchipow(alpha,T5[2],T5[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel6out, details=T)
    test$non_identified_parameters
  }
  
  ### Proximal Bidir @T1 ------------------------------------------------------
  Model_MR7  <-  omxSetParameters(MRModel1out, labels=c('bxy1','byx1'), values=c(0), free=c(F))
  MRModel7out <- mxTryHard(Model_MR7,50)
  mxCompare(MRModel1out, MRModel7out)
  T6  <- as.numeric(mxCompare(MRModel1out, MRModel7out)[2,7:8])
  if (T6[1]<.000001)  {T6[1] = 0}   # ncp = 0
  keep[irow,15] <- T6[1]
  getchipow(alpha,T6[2],T6[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel7out, details=T)
    test$non_identified_parameters
  }
  
  ### Proximal Bidir @T2 ------------------------------------------------------
  Model_MR8  <-  omxSetParameters(MRModel1out, labels=c('bxy2','byx2'), values=c(0), free=c(F))
  MRModel8out <- mxTryHard(Model_MR8,50)
  mxCompare(MRModel1out, MRModel8out)
  T7  <- as.numeric(mxCompare(MRModel1out, MRModel8out)[2,7:8])
  if (T7[1]<.000001)  {T7[1] = 0}   # ncp = 0
  keep[irow,16] <- T7[1]
  getchipow(alpha,T7[2],T7[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel8out, details=T)
    test$non_identified_parameters
  }
  
  
  ### Proximal Bidir @T1+T2 ---------------------------------------------------
  Model_MR9  <-  omxSetParameters(MRModel1out, labels=c('bxy1','bxy2','byx1','byx2'), values=c(0), free=c(F))
  MRModel9out <- mxTryHard(Model_MR9,50)
  mxCompare(MRModel1out, MRModel9out)
  T8  <- as.numeric(mxCompare(MRModel1out, MRModel9out)[2,7:8])
  if (T8[1]<.000001)  {T8[1] = 0}   # ncp = 0
  keep[irow,17] <- T8[1]
  getchipow(alpha,T8[2],T8[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel9out, details=T)
    test$non_identified_parameters
  }
  
  ### Bidir Distal + Proximal @T1 ---------------------------------------------
  Model_MR10  <-  omxSetParameters(MRModel1out, labels=c('bxy1','bxy12','byx1','byx12'), values=c(0), free=c(F))
  MRModel10out <- mxTryHard(Model_MR10,50)
  mxCompare(MRModel1out, MRModel10out)
  T9  <- as.numeric(mxCompare(MRModel1out, MRModel10out)[2,7:8])
  if (T9[1]<.000001)  {T9[1] = 0}   # ncp = 0
  keep[irow,18] <- T9[1]
  getchipow(alpha,T9[2],T9[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel10out, details=T)
    test$non_identified_parameters
  }
  
  ### Bidir Distal + Proximal @T2 ----------------------------------------------
  Model_MR11  <-  omxSetParameters(MRModel1out, labels=c('bxy2','bxy12','byx2','byx12'), values=c(0), free=c(F))
  MRModel11out <- mxTryHard(Model_MR11,50)
  mxCompare(MRModel1out, MRModel11out)
  T10  <- as.numeric(mxCompare(MRModel1out, MRModel11out)[2,7:8])
  if (T10[1]<.000001)  {T10[1] = 0}   # ncp = 0
  keep[irow,19] <- T10[1]
  getchipow(alpha,T10[2],T10[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel11out, details=T)
    test$non_identified_parameters
  }
  
  
  ### Bidir All Three -----------------------------------------------------------
  Model_MR12  <-  omxSetParameters(MRModel1out, labels=c('bxy1','bxy2','bxy12','byx1','byx2','byx12'), values=c(0), free=c(F))
  MRModel12out <- mxTryHard(Model_MR12, 50)
  mxCompare(MRModel1out, MRModel12out)
  T11  <- as.numeric(mxCompare(MRModel1out, MRModel12out)[2,7:8])
  if (T11[1]<.000001)  {T11[1] = 0}   # ncp = 0
  keep[irow,20] <- T11[1]
  getchipow(alpha,T11[2],T11[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel12out, details=T)
    test$non_identified_parameters
  }
  
  
  ## Some more unidirectional tests
  ### Proximal X>Y @T1+T2 ---------------------------------------------------
  Model_MR13  <-  omxSetParameters(MRModel1out, labels=c('bxy1','bxy2'), values=c(0), free=c(F))
  MRModel13out <- mxTryHard(Model_MR13,50)
  mxCompare(MRModel1out, MRModel5out)
  T12  <- as.numeric(mxCompare(MRModel1out, MRModel13out)[2,7:8])
  if (T12[1]<.000001)  {T12[1] = 0}   # ncp = 0
  keep[irow,10] <- T12[1]
  getchipow(alpha,T12[2],T12[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel13out, details=T)
    test$non_identified_parameters
  }
  
  ### Distal + Proximal X>Y @T1 ---------------------------------------------
  Model_MR14  <-  omxSetParameters(MRModel1out, labels=c('bxy1','bxy12'), values=c(0), free=c(F))
  MRModel14out <- mxTryHard(Model_MR14,50)
  mxCompare(MRModel1out, MRModel14out)
  T13  <- as.numeric(mxCompare(MRModel1out, MRModel14out)[2,7:8])
  if (T13[1]<.000001)  {T13[1] = 0}   # ncp = 0
  keep[irow,11] <- T13[1]
  getchipow(alpha,T13[2],T13[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel14out, details=T)
    test$non_identified_parameters
  }
  
  ### Distal + Proximal X>Y @T2 ----------------------------------------------
  Model_MR15  <-  omxSetParameters(MRModel1out, labels=c('bxy2','bxy12'), values=c(0), free=c(F))
  MRModel15out <- mxTryHard(Model_MR15,50)
  mxCompare(MRModel1out, MRModel15out)
  T14  <- as.numeric(mxCompare(MRModel1out, MRModel15out)[2,7:8])
  if (T14[1]<.000001)  {T14[1] = 0}   # ncp = 0
  keep[irow,12] <- T14[1]
  getchipow(alpha,T14[2],T14[1])
  #
  if (doID) {
    test <- mxCheckIdentification(MRModel15out, details=T)
    test$non_identified_parameters
  }
  
  # 
  ### Effect Estimates --------------------------------------------------
  
  # 
  bxy1_hat  <- summary(MRModel1out)$par[5,5]
  bxy12_hat <- summary(MRModel1out)$par[7,5]
  byx1_hat  <- summary(MRModel1out)$par[8,5]
  byx12_hat <- summary(MRModel1out)$par[9,5]
  bxy2_hat  <- summary(MRModel1out)$par[11,5]
  byx2_hat  <- summary(MRModel1out)$par[12,5]
  # 
  keep[irow, 31:36]  <-  c(bxy12_hat, bxy1_hat, bxy2_hat, 
                           byx12_hat, byx1_hat, byx2_hat)
  
  
  ### Std Path Coefficients ----------------------------------------------------
  
  # Refit with Std data
  Std_Model <- mxModel("MRstd", MRmodel1,expF, fitFunction,
                       mxData( observed=cov2cor(Smod), type="cov", means=mu_mod, numObs=N))
  #
  StdModelout <- mxTryHard(Std_Model, 50)
  
  bxy1_std  <- summary(StdModelout)$par[5,5]
  bxy12_std <- summary(StdModelout)$par[7,5]
  byx1_std  <- summary(StdModelout)$par[8,5]
  byx12_std <- summary(StdModelout)$par[9,5]
  bxy2_std  <- summary(StdModelout)$par[11,5]
  byx2_std  <- summary(StdModelout)$par[12,5]
  # 
  keep[irow, 37:42]  <-  c(bxy12_std, bxy1_std, bxy2_std, 
                           byx12_std, byx1_std, byx2_std)
  
  #
  #
  
  
  print(c(' irow ............',irow,nrows))
  print(c(' irow ............',irow,nrows))
  print(c(' irow ............',irow,nrows))
  
} 


# end loop ------------------------------------------------------------------------

## Time cost
time_end <- Sys.time()
print(time_start)
print(time_end)

print(paste("Time cost for full simulation =", (time_end - time_start)))


## Check covariances

## At shortest time lag
MinDis
MinCov
MinCor

write.csv(MinCor, "cor_at_min_lag.csv")

## At longest time lag
MaxDis
MaxCov
MaxCor

write.csv(MaxCor, "cor_at_max_lag.csv")


## Save the results

write.table(keep, file=resfile) 


