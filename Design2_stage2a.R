arrID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
tic <- proc.time()

library(OmicsPLS)
library(PO2PLS)
library(magrittr)
library(tidyverse)
library(pryr)
library(glmnet)
library(caret)
library(MASS)

# 1. load X and C --------------------------------------------------------------
source('00_function.R')
load('X_all_raw_new.RData') # X is generated from PO2PLS once and fixed
load('C_raw.RData') # C_ij iid b(1,0.8)

r=8; rx=2; ry=1
Nh=1000
a=t(rep(2,r)); b=t(rep(0,r))

alpha.G <- 0.3

for(N in c(200,2000)){
  for(dim in 1:2){
    p <- c(100, 2000)[dim]
    q <- c(15, 25)[dim]
    for(noi in 1:2){
      alpha.F <- c(.1,.6)[noi]
      
# 2.  generate YZ by X ---------------------------------------------------------
      
      # training
      X.train <- X_all_raw_new[1:N,1:p] %>% scale(scale = F)
      # test
      X.test <- X_all_raw_new[(2000+1):(2000+Nh), 1:p] %>% scale(scale = F)
      X.all <- rbind(X.train, X.test)
      
      #### Generate Y ####
      W.all <- prcomp(X.all)$rotation[,1:(1+r-1)]
      Tt.all <- X.all %*% W.all
      var.X <- cov(X.all)
      Sig.Tt <- t(W.all) %*% var.X %*% W.all
      C <- C_raw[1:q, 1:r]
      
      sig2F <- alpha.F/(1-alpha.F) * (C %*% Sig.Tt %*% t(C)) %>% diag() %>% mean()
      Ff <- matrix(rnorm((N+Nh)*q, mean=0, sd=sqrt(sig2F)), nrow=(N+Nh), ncol=q)
      Y.all <- Tt.all %*% t(C) + Ff
      Y.train <- Y.all[1:N, ]
      Y.test <- Y.all[(N+1):(N+Nh), ]
      
      Tt.train <- Tt.all[1:N, ] %>% as.data.frame
      Tt.test <- Tt.all[(N+1):(N+Nh), ] %>% as.data.frame
      
      #### Generate Z ####
      
      # generate l
      l <- prcomp(Tt.all %*% t(C))$rotation[,1]
      
      # variance of G
      sig2G <- alpha.G/(1-alpha.G) * (t(l) %*% C %*% diag(diag(Sig.Tt)) %*% t(C) %*% l)
      
      # generate Z
      Gg <- matrix(rnorm((N+Nh)*1, mean=0, sd=sqrt(sig2G)), nrow=(N+Nh), ncol=1)
      Z.all <- Tt.all %*% t(C) %*% l + Gg
      
      Z.train <- Z.all[1:N, ]
      Z.test <- Z.all[(N+1):(N+Nh), ]
      
      #### APPROACH 1: PRS ####
      cat("Starting PRS method...\n")
      l.min <- sapply(1:8, function(ii){
        cv.glmnet(X.train, Y.train[,ii], alpha = 1)$lambda.min
      })
      
      cat("Best Lambda Lasso =", "\n")
      print(summary(l.min))
      l.min <- quantile(l.min, 0.25) # use the same lambda in each column of Y
      
      ## define function of fitting Lasso to Y(each column) and X
      fit.LASSO <- function(Y.col){
        glmnet(X.train, Y.col, alpha=1, lambda=l.min)
      }
      
      models.LASSO <- lapply(1:ncol(Y.train), function(i){
        fit.LASSO(Y.train[,i])
      })
      
      Y.train.pred <- sapply(models.LASSO, 
                             function(model){
                               predict(model, newx=X.train)
                             }
      )
      
      Y.test.pred <- sapply(models.LASSO, 
                            function(model){
                              predict(model, newx=X.test)
                            }
      )
      
      ## Stage 2: Ridge Z~Y
      fit.Ridge <- cv.glmnet(Y.train, Z.train, alpha = 0)
      cat("\nBest Lambda Ridge =", fit.Ridge$lambda.min, "\n")
      
      Z.train.pred <- predict(fit.Ridge, newx = Y.train.pred, s = "lambda.min")
      Z.test.pred <- predict(fit.Ridge, newx = Y.test.pred, s = "lambda.min")
      
      sst.Z <- sum((Z.test - mean(Z.test))^2)
      sse.Z <- sum((Z.test - Z.test.pred)^2)
      rsq.Z.prs <- 1 - sse.Z / sst.Z
      
      rmse.Z.prs.train <- sqrt(mean((Z.train - Z.train.pred)^2))
      rmse.Z.prs.test <- sqrt(mean((Z.test - Z.test.pred)^2))
      
      #### APPROACH 2: O2PLS ####
      
      ## follow the previous procedure
      fit.o2m <- o2m(X.train, Y.train, n=r,nx=rx,ny=ry)
      Tt.train.estim <- X.train %*% fit.o2m$W. %>% as.data.frame
      Tt.test.estim <- X.test %*% fit.o2m$W. %>% as.data.frame
      
      # Z~T linear regression
      Tt.o2mfit <- lm(Z.train ~ . -1, data = Tt.train.estim)
      Z.train.pred <- predict(Tt.o2mfit, newdata = Tt.train.estim)
      Z.test.pred <- predict(Tt.o2mfit, newdata = Tt.test.estim)
      
      sst.Z <- sum((Z.test - mean(Z.test))^2)
      sse.Z <- sum((Z.test - Z.test.pred)^2)
      rsq.Z.o2m <- 1 - sse.Z / sst.Z
      
      rmse.Z.o2m.train <- sqrt(mean((Z.train - Z.train.pred)^2))
      rmse.Z.o2m.test <- sqrt(mean((Z.test - Z.test.pred)^2))
      
      #### APPROACH 3: PO2PLS ####
      fit.po2m <- PO2PLS(X.train, Y.train, r=r,rx=rx,ry=ry,steps=1e3)
      Tt.train.estim <- X.train %*% fit.po2m$par$W %>% data.frame
      Tt.test.estim <- X.test %*% fit.po2m$par$W %>% data.frame
      
      # Z~T linear regression
      Tt.po2mfit <- lm(Z.train ~ . -1, data = Tt.train.estim)
      Z.train.pred <- predict(Tt.po2mfit, newdata = Tt.train.estim)
      Z.test.pred <- predict(Tt.po2mfit, newdata = Tt.test.estim)
      
      sst.Z <- sum((Z.test - mean(Z.test))^2)
      sse.Z <- sum((Z.test - Z.test.pred)^2)
      rsq.Z.po2m <- 1 - sse.Z / sst.Z
      rsq.Z.po2m
      
      rmse.Z.po2m.train <- sqrt(mean((Z.train - Z.train.pred)^2))
      rmse.Z.po2m.test <- sqrt(mean((Z.test - Z.test.pred)^2))
      
      rsq.Z <- tibble(PRS = rsq.Z.prs, O2M = rsq.Z.o2m, PO2M = rsq.Z.po2m)
      rmse.Z <- tibble(
        set = c('train', 'train', 'train', 'test', 'test', 'test'),
        model = c('PRS', 'O2M', 'PO2M', 'PRS', 'O2M', 'PO2M'),
        value = c(rmse.Z.prs.train, rmse.Z.o2m.train, rmse.Z.po2m.train, rmse.Z.prs.test, rmse.Z.o2m.test, rmse.Z.po2m.test)
      )
      save(rsq.Z, file=paste0("outp_N_",N,"_dim_",dim,"_Noi_",noi,"_ID_",arrID, ".RData"))
      save(rmse.Z, file=paste0("rmse_N_",N,"_dim_",dim,"_Noi_",noi,"_ID_",arrID, ".RData"))
    }
  }
}

proc.time() - tic
gc()
mem_used()
