## Correlation of the residuals IV-CLPM with Varying Time intervals

## Simulate time-series with varying resid corr for 150 time points.
## Fit a series of IV-CLPM with 2 panels but increasing panel lag
## T1 in the models is set at T=100 of the time-series.
## T2 in the models increases sequentially from T101 to T150 of the time-series.

# clear workspace
rm(list=ls(all=TRUE))

# SET THE WORKING DIRECTORY
getwd()
resfile <- 'resfile_rescor_IV-CLPM.dta'  # output file for the results

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
px <- .08

# IVy --> Y_t
py <- .08

# r(IVx,IVy)
rIV <- 0.25      # Correlation of the 2 IVs

# X_t --> Y_{t+1}
bxy <- .2     # 1st-order Causal path coef.

# Y_t --> X_{t+1}
byx <- .2     # 1st-order Causal path coef. 

# AR1
bx  <- .7      # AR1 of X 
by  <- .7      # AR1 of Y 

# Varying parameters -----------------------------------------------------------

# Residual Cross-sectional Correlation of X and Y
rxys <- c(.1,.3)   

## CLPM Panel Time Lag
startT <- 100              # Wave 1 across all models = T100
Tdis <- seq(1:50)          # Wave 2 changed from T101 to T150
# 50 models will be fitted to each time-series generated using the parameters specified above


## Number of cells in the factorial design
nrows <- length(rxys) * length(Tdis)
nrows


# A matrix to store the parameters of interest
keep <- matrix(0,nrows,7)
# 7 parameters
colnames(keep) <-  c(
  ## Data-generating parameters
  'rxy','Tdis',
  # ## Recovered model estimates
  'b_distal_xy','b_prox1_xy','b_prox2_xy',
  'r_xy1','r_xy2'
)


# --------------------------------------------------------------------------------------------


# start loop ---------------------------------------------------------------------------------

# start time
time_start <- Sys.time()

irow <- 0

for (r_xy in rxys) {         # Residual Correlation
  for (idis in Tdis) {       # Panel Lag
    
    irow  <-  irow+1
    
    # Data-generating parameters
    keep[irow,1:2] <-  c(r_xy,idis)
    
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
    
    # Variables for CLPM:
    isel <- c(1, 2,                      # IVx, IVy
              Ts[1]+2 , Ts[1]+T_+2,      # X and Y at Wave 1
              Ts[2]+2, Ts[2]+T_+2)       # X and Y at Wave 2
    Smod <- S_[isel, isel]
    
    ## check covariances 
    # round(cov2cor(Smod),3)
    
    
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
    
    ## Effect Estimates --------------------------------------------------
    
    # 
    bxy1_hat  <- summary(MRModel1out)$par[5,5]
    bxy12_hat <- summary(MRModel1out)$par[7,5]
    bxy2_hat  <- summary(MRModel1out)$par[11,5]
    # 
    keep[irow, 3:5]  <-  c(bxy12_hat, bxy1_hat, bxy2_hat)
    
    # Correlation Path between the Residuals
    
    res_cor <- cov2cor(MRModel1out@output$matrices$MR.PS)
    
    r_xy1   <- as.numeric(res_cor[4,3]) 
    r_xy2   <- as.numeric(res_cor[6,5]) 
    
    keep[irow, 6:7]  <-  c(r_xy1, r_xy2)
    
    
    
    print(c(' irow ............',irow,nrows))
    print(c(' irow ............',irow,nrows))
    print(c(' irow ............',irow,nrows))
    
  } 
}



# end loop ------------------------------------------------------------------------

## Time cost
time_end <- Sys.time()
print(time_start)
print(time_end)

print(paste("Time cost for full simulation =", (time_end - time_start)))


## Save the results

write.table(keep, file=resfile) 




# PART II. RESULTS ---------------------------------------------------------------

keep <- read.table("resfile_rescor_IV-CLPM.dta")
colnames(keep) <-  c(
  ## Data-generating parameters
  'rxy','Tdis',
  ## Recovered parameter estimates
  'b_distal_xy','b_prox1_xy','b_prox2_xy',
  'r_xy1','r_xy2'
)

psych::describe(keep)


# Plot --------------------------------------------------------------------------

library(tidyverse)
library(RColorBrewer)
library(kableExtra)

keep <- as_tibble(keep)


# facet panel labels
rxy_labs = c("r_exy = 0.1","r_exy = 0.3")
names(rxy_labs) = c("0.1","0.3")

keep %>% 
  pivot_longer(cols = c('b_prox1_xy','b_distal_xy','b_prox2_xy','r_xy1','r_xy2'), 
               names_to = "parameter",
               values_to = "estimate") %>% 
  mutate(`Parameter Type` = as_factor(ifelse(parameter %in% c("b_prox1_xy","b_distal_xy","b_prox2_xy"),
                                             "Regression", "Correlation")),
         Parameter = as_factor(parameter),
         Parameter = fct_recode(Parameter, "Proximal_W1" = "b_prox1_xy",
                                "Distal" = "b_distal_xy",
                                "Proximal_W2" = "b_prox2_xy",
                                "r_exy1" = "r_xy1",
                                "r_exy2" = "r_xy2"),
         Parameter = fct_relevel(Parameter, "Proximal_W1")) %>% 
  ggplot(aes(Tdis, estimate, color = Parameter)) +
  geom_line(aes(linetype = `Parameter Type`), stat = "summary") +
  geom_point(aes(shape = `Parameter Type`), size = 1, stat = "summary") +
  facet_wrap(~ rxy, labeller = labeller(rxy = rxy_labs)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = brewer.pal(n=7,name="Dark2")[c(1:3,6:7)]) +
  coord_cartesian(xlim = c(0,40), ylim = c(min(keep$r_xy2), max(keep$b_prox1_xy))) +
  scale_y_continuous(breaks = seq(-.6,.6,0.3)) +
  labs(title = "IV-CLPM", y = "Estimate", x = "Time Interval between Study Waves") +
  theme_bw(14) + 
  theme(legend.position = "right")

