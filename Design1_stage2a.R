arrID <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))
tic <- proc.time()

library(OmicsPLS)
library(magrittr)
library(tidyverse)
library(pryr)
library(glmnet)
library(caret)
library(PO2PLS)

source('00_function.R')

Nh=1000 # test set
alpha.z <- 0.3

for(N in c(200,2000)){
  for(dim in 1:2){
    p <- c(100, 2000)[dim]
    q <- c(15, 25)[dim]
    for(noi in 1:2){ # noise part
      alpha.x <- c(.1,.95)[noi]
      alpha.y <- c(.1,.6)[noi]
      
      for (alpha.tu in c(0.1,0.4)){ # heterogeneity
        alpha <- c(alpha.x, alpha.y, alpha.tu)
        
# 1. Generate XYZ --------------------------------------------------------------
        r=5; rx=5; ry=5 # joint & specific components
        a=t(rep(2,r))
        
        X <- matrix(0, nrow=N, ncol=p); Y <- matrix(0, nrow=N, ncol=q)
        params.true <- generate_params(X, Y, r, rx, ry, alpha=alpha, type='unit')
        sig2G <- alpha.z/(1-alpha.z) * (a%*%params.true$SigT%*%t(a))
        sig2G <- sig2G %>% as.numeric()
        params.true$sig2G <- sig2G
        params.true$a <- a
        
        dat <- generate_data_inclZ(N+Nh, params=params.true, distr=rnorm)
        dat$Z <- dat$Z %>% as.matrix()
        
        ## training data
        X.train <- dat$X[1:N, ] %>% scale(scale = F)
        Y.train <- dat$Y[1:N, ] %>% scale(scale = F)
        Z.train <- dat$Z[1:N, ] %>% scale(scale = F)
        ## test data
        X.test <- dat$X[(N+1):(N+Nh), ] %>% scale(scale = F)
        Y.test <- dat$Y[(N+1):(N+Nh), ] %>% scale(scale = F)
        Z.test <- dat$Z[(N+1):(N+Nh), ] %>% scale(scale = F)
       
# 2. 3 models ------------------------------------------------------------------
        
        ###  -------------------------------  ###
        #### APPROACH 1: PRS (Ridge & LASSO) ####
        ###  -------------------------------  ###
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
        
        ## Ridge Z~Y
        fit.Ridge <- cv.glmnet(Y.test.pred, Z.test, alpha = 0)
        cat("\nBest Lambda Ridge =", fit.Ridge$lambda.min, "\n")
        
        Z.train.pred <- predict(fit.Ridge, newx = Y.train.pred, s = "lambda.min")
        Z.test.pred <- predict(fit.Ridge, newx = Y.test.pred, s = "lambda.min")
        
        sst.Z <- sum((Z.test - mean(Z.test))^2)
        sse.Z <- sum((Z.test - Z.test.pred)^2)
        rsq.Z.prs <- 1 - sse.Z / sst.Z
        
        rmse.Z.prs.train <- sqrt(mean((Z.train - Z.train.pred)^2))
        rmse.Z.prs.test <- sqrt(mean((Z.test - Z.test.pred)^2))
        
        ###  -----------------  ###
        #### APPROACH 2: O2PLS ####
        ###  -----------------  ###
        
        ## choose r before fitting model
        choose_r <- function(r_try = 1:(q-rx)){
          r_hats <- sapply(r_try, function(ii){
            fit.o2m <- o2m(X.train, Y.train, n=ii,nx=0,ny=0)
            T.estim <- X.train %*% fit.o2m$W.
            X.test.A <- X.test[1:(Nh/2),]
            X.test.B <- X.test[(Nh/2+1):Nh,]
            Z.test.A <- Z.test[1:(Nh/2)]
            Z.test.B <- Z.test[-(1:(Nh/2))]
            T.test.A.estim <- X.test.A %*% fit.o2m$W.
            T.test.B.estim <- X.test.B %*% fit.o2m$W.
            
            # Z~T regression
            T.dat <- data.frame(V = T.test.A.estim)
            T.o2mfit <- lm(Z.test.A ~ ., data = T.dat) # estimation of a and a0
            Z.test.A.pred <- predict(T.o2mfit, newdata = data.frame(V = T.test.A.estim))
            Z.test.B.pred <- predict(T.o2mfit, newdata = data.frame(V = T.test.B.estim))
            
            sst.Z <- sum((Z.test.B - mean(Z.test.B))^2)
            sse.Z <- sum((Z.test.B - Z.test.B.pred)^2)
            rsq.Z.o2m <- 1 - sse.Z / sst.Z
            return(rsq.Z.o2m)
          })
          r_try[which.max(r_hats)]
        }
        r <- choose_r(1:(q-rx))
        
        cat("Starting O2PLS method...\n")
        fit.o2m <- o2m(X.train, Y.train, n=r,nx=rx,ny=ry)
        T.estim <- X.train %*% fit.o2m$W.
        T.test.estim <- X.test %*% fit.o2m$W.
        
        # Z~T regression
        T.dat <- data.frame(V = T.test.estim)
        T.o2mfit <- lm(Z.test ~ ., data = T.dat) # estimation of a and a0
        Z.test.pred <- predict(T.o2mfit, newdata = data.frame(V = T.test.estim))
        
        sst.Z <- sum((Z.test - mean(Z.test))^2)
        sse.Z <- sum((Z.test - Z.test.pred)^2)
        rsq.Z.o2m <- 1 - sse.Z / sst.Z
        
        rmse.Z.o2m.train <- sqrt(mean((Z.train - Z.train.pred)^2))
        rmse.Z.o2m.test <- sqrt(mean((Z.test - Z.test.pred)^2))
        
        ###  ------------------  ###
        #### APPROACH 3: PO2PLS ####
        ###  ------------------  ###
        
        cat("Starting PO2PLS method...\n")
        fit.po2m <- PO2PLS(X.train, Y.train, r=r,rx=rx,ry=ry,steps=1e3) # r is same as APPROACH 2
        T.estim <- X.train %*% fit.po2m$par$W
        T.test.estim <- X.test %*% fit.po2m$par$W
        
        # Z~T regression
        T.dat <- data.frame(V = T.test.estim)
        T.po2mfit <- lm(Z.test ~ ., data = T.dat) # estimation of a and a0
        Z.test.pred <- predict(T.po2mfit, newdata = data.frame(V = T.test.estim))
        
        sst.Z <- sum((Z.test - mean(Z.test))^2)
        sse.Z <- sum((Z.test - Z.test.pred)^2)
        rsq.Z.po2m <- 1 - sse.Z / sst.Z
        
        rmse.Z.po2m.train <- sqrt(mean((Z.train - Z.train.pred)^2))
        rmse.Z.po2m.test <- sqrt(mean((Z.test - Z.test.pred)^2))

        # save RData
        rsq.Z <- tibble(PRS = rsq.Z.prs, O2M = rsq.Z.o2m, PO2M = rsq.Z.po2m)
        save(rsq.Z, file=paste0("outp_N_",N,"_dim_",dim,"_Noi_",noi,"_alphatu_",alpha.tu, "_ID_",arrID, ".RData"))
        
        rmse.Z <- tibble(
          set = c('train', 'train', 'train', 'test', 'test', 'test'),
          model = c('PRS', 'O2M', 'PO2M', 'PRS', 'O2M', 'PO2M'),
          value = c(rmse.Z.prs.train, rmse.Z.o2m.train, rmse.Z.po2m.train, rmse.Z.prs.test, rmse.Z.o2m.test, rmse.Z.po2m.test)
        )
        save(rmse.Z, file=paste0("rmse_N_",N,"_dim_",dim,"_Noi_",noi,"_alphatu_",alpha.tu, "_ID_",arrID, ".RData"))
      }
    }
  }
}

proc.time() - tic
gc()
mem_used()
