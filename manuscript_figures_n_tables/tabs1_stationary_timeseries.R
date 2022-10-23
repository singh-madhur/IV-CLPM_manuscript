## Simulated Stationary Time Series

## Aim: Table to show stationary covariances through the time period used for model fitting.

## Simulate time-series with varying causal, AR, and resid corr, for 150 time points.
## T1 in the models is set at T=100 of the time-series.
## T2 in the models increases sequentially from T101 to T150 of the time-series.

## Clear the workspace
rm(list=ls(all=TRUE))


# Library
library(xlsx)



## Data parameters --------------------------------------------------------------
# 

N <- 1000 # sample size
px <- .08 # IVx --> Xt
py <- .08 # IVy --> Yt
bx <- .7 # ARx
by <- .7 # ARy
bxy <- .2 # X_t --> Y_{t+1}
byx <- .2 # Y_t --> X_{t+1}
r_xy <- .1   # residual correlation of X and Y
rIV <- .25

#
T_ <- 150   # total time points
nv <- T_*2  # total observed variables (except IV)


## Data simulation ----------------------------------------------------

vx <- 1               # Residual variance of X
vy <-  1              # Residual Variance of Y
vxy <- r_xy           # Residual non-causal covariance

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
BE_[1:5,1:5]
BE_[151:155,1:5]
BE_[1:5,151:155]
BE_[151:155,151:155]

## Additional columns for IVx-->X, IVy-->Y
Gamma <- matrix(0,nv+2,2)                     
Gamma[1:(nv+2),1] <- c(0,0,rep(px,T_),rep(0,T_))  # IVx -> X, rows 3:152 are X1:X150                    
Gamma[1:(nv+2),2] <- c(0,0,rep(0,T_),rep(py,T_))  # IVy -> Y, rows 153:302 are Y1:Y150

# check
Gamma[1:6,]
Gamma[150:155,]

# Bind the gamma columns to the beta matrix
BE_ <- cbind(Gamma[3:(nv+2),], BE_)  
BE_ <- rbind(rep(0,(nv+2)), rep(0,(nv+2)), BE_)      # Row 1 and 2 are for the IVs

# check
BE_[1:7,1:7]

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
PS_[1:6,1:6]
PS_[150:155, 1:6]

# Matrices for the algebra
I_ <- diag(nv+2)
IBE_ <- solve(I_-BE_)

# Expected Covarianc Matrix
S_ <- IBE_%*%PS_%*%t(IBE_)
rownames(S_) <- colnames(S_) <- c('IVx', 'IVy', paste('X',1:T_, sep=""),paste('Y',(1:T_), sep=""))


## Covariance -------------------------------------------------------------------

## Select the time points for cov table
Ts <- c(100,101,125,150)  # Wave 1, Nearest Wave 2, Midway, Farthest Wave 2

isel <- cbind(1,2)                      # IVx, IVy
for (ii in 1:length(Ts)) {
  isel <- cbind(isel, Ts[ii]+2)         # X at the selected Ts
}
for (jj in 1:length(Ts)) {
  isel <- cbind(isel, Ts[jj]+T_+2)      # Y at the selected Ts
}

isel

# Extract the selected variables from the full covariance matrix
Smod <- S_[isel, isel]
round(Smod, 4)

datcor <- round(cov2cor(Smod),4)
datcor[upper.tri(datcor, diag = F)] <- NA 

datcor

# Save -------------------------------------------------------------------------

resfile <- 'stationaryCov.xlsx'

xlsx::write.xlsx(datcor, file = resfile, showNA = FALSE)
