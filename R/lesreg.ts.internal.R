se = 0.1

RQObjectiveFunction <- function(beta,
                                out,
                                model,
                                data,
                                obs,
                                pval,
                                emp_qnt,
                                varPredict = F){

  #This is optimization function, that uses underlying cpp files to calculate VaR.
  #It comes in two main modes:
  # - out = 1, where it returns RQ (cost function) value
  # - out = 2, where it returns calculated VaR and Hit values

  if (varPredict) {
    VaR = rep(0, obs + 1)
  } else {
    VaR = rep(0, obs)
    Hit = VaR
  }

  VaR[1] <- 1*emp_qnt

  VaR <- if (model == 1) {
    caviar_SAV(beta, data, VaR[1], VaR, obs, varPredict)
  } else if (model == 2) {
    caviar_AS(beta, data, VaR[1], VaR, obs, varPredict)
  } else if (model == 3) {
    -caviar_GARCH(beta, data, abs(VaR[1]), VaR, obs, varPredict)
  } else if (model == 4) {
    caviar_SAV12(beta, data, VaR[1], VaR, obs, varPredict)
  }else if (model == 5){
    caviar_AS12(beta, data, VaR[1], VaR, obs, varPredict)
  }

  if (!varPredict){
    Hit <- pval - (data < VaR)
    if (out == 1) {
      RQ = t(Hit) %*% (data - VaR)
      if (is.infinite(RQ))
        RQ = 1e+100
      return(RQ)
    } else if (out == 2) {
      return(cbind(VaR, Hit))
    }
  } else if (varPredict) {
    return(VaR)
  }

}

caviarOptim <- function(data,
                        model = 1,
                        pval = 0.01,
                        REP = 5,
                        MAXITER = 500,
                        predict = F,
                        initial = NULL){
  #This function is modified from the previous caviarOptim function, the definition of VaR have been changed into the tau
  #level quantile, instead of the inverse of the quantile
  #
  # Input:
  #
  # - data         - numeric         - data for the CAViaR model to calculate risk
  #                  (convertible to numerical)
  # - model        - numeric         - model type (1 - SAV, 2 - AS, 3 - GARCH, 4 - SAV12, 5-AS12) (defaults to 1)
  # - pval         - numeric         - level of value at risk probability (defaults to 0.01)
  # - REP          - numeric         - number of hessian optimization procedure repeats
  # - MAXITER      - numeric         - maximum number of optimization procedure iterations
  # - predict      - logical         - whether to return only predicted value for next period or whole list
  # - initial      - numeric         - initial estimate of quantile parameter
  #
  # Output:
  # - outL        - list:
  #   - bestVals     - numeric         - best cost function values
  #   - bestPar      - numeric         - best parameters
  #   - VaR          - numeric         - calculated Value at Risk
  #   - bestRQ       - numeric         - best cost function among all
  #   - VarPredict   - numeric         - predicted value at risk for the next period (tail(VaR,1))

  if (length(data) < 2) {
    stop("data should be series of at least 300 length")
  } else if (length(data) < 300) {
    warning(paste0("for the best empirical quantile estimation, it is suggested to provide data of at least 300 length. Only ",
                   length(data), " samples provided"))
  }

  tryCatch(data <- as.numeric(data),
           error = function(e){
             print("Can't convert to numeric! Provide data convertible to numerical!")
             stop()
           })

  if(is.null(initial)){
    if (!model %in% c(1,2,3,4,5)) {
      stop('Wrong MODEL selected')
    } else if (model == 1 || model == 3) {
      initialTargetVectors <- matrix(runif(10000*3),ncol=3)
      nInitialCond = 10
    } else if (model == 2) {
      initialTargetVectors <- matrix(runif(10000*4),ncol=4)
      nInitialCond = 15
    } else if (model == 4) {
      initialTargetVectors <- matrix(runif(10000*4,-2,2),ncol=4)
      nInitialCond = 15
    } else if (model == 5){
      initialTargetVectors <- matrix(runif(10000*6,-2,2),ncol=6)
      nInitialCond = 20
    }
    obs <- length(data)

    VaR <- rep(0, obs)
    Hit <- VaR

    emp_qnt <- if (obs < 300){
      quantile(data[1:obs], pval)
    } else {
      quantile(data[1:300], pval)
    }

    RQfval <- apply(initialTargetVectors, 1, RQObjectiveFunction, 1, model, data, obs, pval, emp_qnt)

    BestInitialCond <- initialTargetVectors[order(RQfval), ][1:nInitialCond, ]}else{
      obs <- length(data)
      emp_qnt <- if (obs < 300){
        quantile(data[1:obs], pval)
      } else {
        quantile(data[1:300], pval)
      }
      BestInitialCond<- matrix(initial, nrow = 1)
    }
  para_nums <- 3*(model %in% c(1,3)) + 4*(model%in%c(2,4)) + 6*(model == 5)
  RQoptim <- cbind(NA, BestInitialCond, NA)

  met <- 'L-BFGS-B'
  met_hes <- "Nelder-Mead"
  con <- list(maxit = MAXITER)
  if(model == 1){
    low = c(-Inf, -1+ 1e-3, -Inf)
    up = c(Inf,1-1e-3,Inf)
  }else if(model == 2){
    low = c(0, -1+ 1e-3, 0, 0)
    up = c(Inf,1-1e-3,Inf,Inf)
  }else if (model == 3){
    low = c(1e-3, 1e-3,0)
    up = c(Inf, 1-1e-3, Inf)
  } else if (model == 4){
    low = c(-Inf,-1+ 1e-3,-Inf,-Inf)
    up = c(Inf, 1 - 1e-3, Inf, Inf)
  } else if (model == 5){
    low = c(-Inf,-1+ 1e-3,-Inf,-Inf,-Inf,-Inf)
    up = c(Inf, 1 - 1e-3, Inf, Inf, Inf, Inf)
  }

  con <- list(maxit = MAXITER)

  for (i in 1:nrow(BestInitialCond)) {
    # initial optimization procedure
    vOptim <- optim(BestInitialCond[i,],
                    RQObjectiveFunction,
                    out = 1,
                    model = model,
                    data = data,
                    obs = obs,
                    pval = pval,
                    emp_qnt = emp_qnt,
                    method = met,
                    lower = low,
                    upper = up,
                    control = con)

    RQoptim[i, 1] <- vOptim$value
    RQoptim[i, 2:(para_nums+1)] <- vOptim$par
    if(model == 3){
      beta_tmp<- vOptim$par+ stats::runif(length(vOptim$par), min = -se, max = se)
      beta_tmp[2] <- vOptim$par[2]
      beta_tmp[which(beta_tmp<=0)] <- vOptim$par[which(beta_tmp<=0)]
    }else{
      beta_tmp <- vOptim$par+ stats::rnorm(length(vOptim$par), mean = 0, sd = se)
      beta_tmp[2] <- vOptim$par[2]}


    # repeating optim on hessian optimiziation to get more accurate results
    for(j in 1:REP){
      # once optimizing previous best param (using BFGS)
      vOptim <- optim(beta_tmp,
                      RQObjectiveFunction,
                      out = 1,
                      model = model,
                      data = data,
                      obs = obs,
                      pval = pval,
                      emp_qnt = emp_qnt,
                      method = met_hes,
                      control = con)
      # secondly optimizing previous best param with afore method
      if(model == 3){
        beta_tmp<- vOptim$par+ stats::runif(length(vOptim$par), min = -se, max = se)
        beta_tmp[2] <- vOptim$par[2]
        beta_tmp[which(beta_tmp<=0)] <- vOptim$par[which(beta_tmp<=0)]
      }else{
        beta_tmp <- vOptim$par+ stats::rnorm(length(vOptim$par), mean = 0, sd = se)
        beta_tmp[2] <- vOptim$par[2]}
      vOptim <- optim(beta_tmp,
                      RQObjectiveFunction,
                      out = 1,
                      model = model,
                      data = data,
                      obs = obs,
                      pval = pval,
                      emp_qnt = emp_qnt,
                      method = met,
                      lower = low,
                      upper = up,
                      control = con)
    }
    #checking whether hessian optimization converged
    if((RQoptim[i, 1] - vOptim$value)>0){
      RQoptim[i, 1] <- vOptim$value
      RQoptim[i,  2:(para_nums+1)] <- vOptim$par
    } else {
      RQoptim[i,ncol(RQoptim)] <- j
    }
  }

  bestPar <- RQoptim[order(RQoptim[,1]), ,drop = F][1, 2:(para_nums+1)]
  VaR<- RQObjectiveFunction(beta = bestPar, out = 2, model = model, data = data,obs = obs, pval = pval,
                            emp_qnt = emp_qnt, varPredict = T)
  VarPredict <- VaR %>% tail(1)
  VaR <- VaR[-length(VaR)]

  if (predict)
    return(VarPredict)
  else {
    outL <- list(bestVals   = RQoptim,
                 bestPar    = bestPar,
                 VaR        = VaR,
                 bestRQ     = RQoptim[order(RQoptim[,1]), ,drop = F][1, 1],
                 VarPredict = VarPredict
    )
    return(outL)
  }

}

ESObjectiveFunction = function(theta,out,model,data,pseudo_y,obs,emp_es,weights,esPred = F){
  #This is optimization function, that uses underlying cpp files to calculate ES.
  #It comes in two main modes:
  # - out = 1, where it returns SE (cost function) value
  # - out = 2, where it returns calculated ES

  ES = rep(0, obs + esPred)
  ES[1] = emp_es

  ES <- if (model == 1) {
    SAV_ES(theta, data, ES[1], ES, obs, esPred)
  } else if (model == 2) {
    AS_ES(theta, data, ES[1], ES, obs, esPred)
  } else if (model == 3) {
    GARCH_ES(theta, data, ES[1], ES, obs, esPred)
  } else if (model == 4){
    SAV12_ES(theta, data, ES[1], ES, obs, esPred)
  } else if (model == 5){
    AS12_ES(theta, data, ES[1], ES, obs, esPred)
  }

  if(!esPred){
    if(out==1){
      return(mean(weights*(pseudo_y - ES)^2 ) )
    }else if(out == 2){
      return(ES)
    }
  }else if(esPred){return(ES)}
}

truncated_value = function(dat,method,tau){
  if(method == 1){
    val = mean((dat^2)[which(dat<=0)])
  }else if(method == 2){
    df = stats::density(dat, bw = "SJ")
    lower_int_boundary = min(df$x)
    pdf = stats::approxfun(df, yleft = 0, yright = 0)

    truncated_point = 0
    div = 1000

    h = (truncated_point - lower_int_boundary) / (div - 1)
    b_approx = seq(lower_int_boundary, truncated_point + h, h)

    y0 = pdf(b_approx)
    y2 = b_approx^2 * y0

    m2 = sum(y2[-1] + y2[-div]) / 2 * h
    val = m2/tau
  }
  return(val)
}

ESEst_Inf_fair = function(qmodel, model, pval, data, weights, predict, REP, MAXITER,covMethod, truncationMethod,order,ini_spec = F){
  tlist_1 = c()
  tlist_2 = c()
  obs = length(data)
  if(ini_spec){
    initials = get_sieve_initial(data = data,qmodel = qmodel, emodel = model, order = order,tau = pval)
    caviar_obj = caviarOptim(data = data, model = qmodel, pval = pval, REP = REP, MAXITER = MAXITER,initial = initials$beta_0)
    beta_hat = caviar_obj$bestPar
    #if(qmodel == 1){beta_hat = c(-beta_hat[1], beta_hat[2],-beta_hat[3])}else if (qmodel == 2){beta_hat = c(-beta_hat[1], beta_hat[2],-beta_hat[3],-beta_hat[4])}
    q_hat = caviar_obj$VaR
    pseudo_y = q_hat +  (data<=q_hat)*(data - q_hat)/pval
    emp_es = (obs <= 300)*mean(data[which(data<=q_hat)] ) + (obs>=300)*mean(data[1:300][which(data[1:300]<=q_hat[1:300])])

    theta_tmp = initials$theta_0
    ES_tmp = get_ES_value(data = data, theta = theta_tmp, emodel = model)
    theta_tmp[1] = theta_tmp[1] - max(ES_tmp)*(max(ES_tmp)>0)

    BestInitialCond <- matrix(theta_tmp,nrow = 1)

    ESoptim <- cbind(NA, BestInitialCond, NA)}else{
      caviar_obj = caviarOptim(data = data, model = qmodel, pval = pval, REP = REP, MAXITER = MAXITER,initial = NULL)
      q_hat = caviar_obj$VaR
      beta_hat = caviar_obj$bestPar
      pseudo_y = q_hat +  (data<=q_hat)*(data - q_hat)/pval
      emp_es = (obs <= 300)*mean(data[which(data<=q_hat)] ) + (obs>=300)*mean(data[1:300][which(data[1:300]<=q_hat[1:300])])
      if(model %in% c(1,2,3)){
        initialTargetVectors <- matrix(runif(10000*3,-1,1)*(model!=3) + runif(10000*3,0,1)*(model==3)
                                       ,ncol=4*(model==2) + 3*(model %in% c(1,3)))} else if (model == 4){
                                         initialTargetVectors <- matrix(runif(10000*4,-2,2),ncol=4)
                                       }else if (model == 5){
                                         initialTargetVectors <- matrix(runif(10000*6,-2,2),ncol=6)
                                       }
      nInitialCond = 10
      ESval <- apply(initialTargetVectors, 1, ESObjectiveFunction, 1, model, data, pseudo_y,obs, emp_es, weights, F)
      BestInitialCond <- initialTargetVectors[order(ESval), ][1:nInitialCond, ]
      ESoptim <- cbind(NA, BestInitialCond, NA)
    }

  met <- 'L-BFGS-B'
  #met <- 'BFGS'
  #met_hes <- 'SANN'
  met_hes <- "Nelder-Mead"
  con <- list(maxit = MAXITER)
  if(model == 1){
    low = c(-Inf, 1e-3, -Inf)
    up = c(-1e-3,1-1e-3,-1e-3)
  }else if(model == 2){
    low = c(-Inf, -1+1e-3, -Inf, -Inf)
    up = c(0,1-1e-3,0,0)
  }else if (model == 3){
    low = c(1e-3,1e-3,0)
    up = c(Inf, 1-1e-3, Inf)
  }else if (model == 4){
    low = c(-Inf,-1+1e-3,-Inf,-Inf)
    up = c(Inf, 1 - 1e-3, Inf, Inf)
  }else if (model == 5){
    low = c(-Inf,-1+1e-3,-Inf,-Inf,-Inf,-Inf)
    up = c(Inf, 1 - 1e-3, Inf, Inf, Inf, Inf)
  }

  for (i in 1:nrow(BestInitialCond)) {
    # initial optimization procedure
    t1 = Sys.time()
    vOptim <- optim(BestInitialCond[i,],
                    ESObjectiveFunction,
                    out = 1,
                    model = model,
                    data = data,
                    pseudo_y = pseudo_y,
                    obs = obs,
                    emp_es = emp_es,
                    weights = weights,
                    method = met,
                    control = con,
                    lower = low,
                    upper = up)
    t2 = Sys.time()
    tlist_1 = c(tlist_1,t2 -t1)
    ESoptim[i, 1] <- vOptim$value
    ESoptim[i, 2:(ncol(BestInitialCond)+1)] <- vOptim$par



    # repeating optim on hessian optimiziation to get more accurate results
    for(j in 1:REP){
      # once optimizing previous best param (using BFGS)
      t1 = Sys.time()
      if(model == 3){
        theta_tmp<- vOptim$par+ stats::runif(length(vOptim$par), min = -se, max = se)
        theta_tmp[2] <- vOptim$par[2]
        theta_tmp[which(theta_tmp<=0)] <- vOptim$par[which(theta_tmp<=0)]
      }else{
        theta_tmp = vOptim$par + stats::rnorm(length(vOptim$par), mean = 0, sd = se)
        theta_tmp[2] = vOptim$par[2]}
      vOptim <- optim(theta_tmp,
                      ESObjectiveFunction,
                      out = 1,
                      model = model,
                      data = data,
                      pseudo_y = pseudo_y,
                      obs = obs,
                      emp_es = emp_es,
                      weights = weights,
                      method = met_hes,
                      control = con)
      t2 = Sys.time()
      tlist_2 = c(tlist_2,t2 - t1)
      # secondly optimizing previous best param with afore method
      t1 = Sys.time()
      if(model == 3){
        theta_tmp<- vOptim$par+ stats::runif(length(vOptim$par), min = -se, max = se)
        theta_tmp[2] <- vOptim$par[2]
        theta_tmp[which(theta_tmp<=0)] <- vOptim$par[which(theta_tmp<=0)]
      }else{
        theta_tmp = vOptim$par + stats::rnorm(length(vOptim$par), mean = 0, sd = se)
        theta_tmp[2] = vOptim$par[2]}
      vOptim <- optim(theta_tmp,
                      ESObjectiveFunction,
                      out = 1,
                      model = model,
                      data = data,
                      pseudo_y = pseudo_y,
                      obs = obs,
                      emp_es = emp_es,
                      weights = weights,
                      method = met,
                      control = con,
                      lower = low,
                      upper = up)
      t2 = Sys.time()
      tlist_2 = c(tlist_2,t2 - t1)
    }
    #checking whether hessian optimization converged
    if(ESoptim[i, 1] - vOptim$value>=0){
      ESoptim[i, 1] <- vOptim$value
      ESoptim[i, 2:(ncol(BestInitialCond)+1)] <- vOptim$par
    } else {
      ESoptim[i,ncol(ESoptim)] <- j
    }
  }

  bestPar <- ESoptim[order(ESoptim[,1]), ,drop = F][1, 2:(ncol(BestInitialCond)+1)]

  ES = ESObjectiveFunction(bestPar, 2, model, data, pseudo_y, obs, emp_es, weights = rep(1,length(data)), predict)
  if(predict == T){ES0 = ES[-length(ES)]}else{ES0 = ES}
  # VC Estimation
  residuals = pseudo_y - ES0
  denom = q_hat - ES0
  winsorize = quantile(abs(denom),0.005)
  denom[which(abs(denom)<=winsorize)] = winsorize
  eta_hat = (data - q_hat)/denom
  nuisance_truncated = truncated_value(eta_hat,truncationMethod,pval)
  if(covMethod == 1){
    conditional_variance = residuals^2
  }else if (covMethod == 2){
    conditional_variance = (q_hat - ES)^2*(1/pval*nuisance_truncated - 1)
  }
  if(nuisance_truncated>=10){
    warning(paste0('Yield extreme truncated quantity estimator! Approach ', covMethod,' can not fit, switch to default method.'))
    conditional_variance = residuals^2
  }
  if(min(conditional_variance)<=0){
    warning(paste0('Yield negative variance estimator! Approach ', covMethod,' can not fit, switch to default method.'))
    conditional_variance = residuals^2
  }
  derivative1 = rep(0,obs)
  derivative2 = rep(0,obs)
  derivative3 = rep(0,obs)
  derivative4 = rep(0,obs)
  derivative5 = rep(0,obs)
  derivative6 = rep(0,obs)

  Sigma = matrix(0, nrow = length(bestPar), ncol = length(bestPar))
  Omega = Sigma
  gradient = matrix(0, nrow = obs, ncol = length(bestPar))

  for( t in 2:obs){
    if(model == 1){
      #SAV
      derivative1[t] = 1 + bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + abs(data[t-1])

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t])

      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])

    }
    else if(model == 2){
      # AS
      derivative1[t] = 1+ bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + data[t-1]*(data[t-1]>0)
      derivative4[t] = bestPar[2]*derivative4[t-1] - data[t-1]*(data[t-1]<0)

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t], derivative4[t])
      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])
    }
    else if (model == 3){
      #GARCH
      derivative1[t] = (1 + 2*bestPar[2]*ES[t-1]*derivative1[t-1]) / (2*ES[t])
      derivative2[t] = (ES[t-1]^2 + 2*bestPar[2]*ES[t-1]*derivative2[t-1]) / (2*ES[t])
      derivative3[t] = (2*bestPar[2]*ES[t-1]*derivative3[t-1] + data[t-1]^2) / (2*ES[t])

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t])

      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])
    }
    else if (model == 4){
      derivative1[t] = 1 + bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + abs(data[t-1])
      derivative4[t] = if(t==2){0}else{bestPar[2]*derivative4[t-1] + abs(data[t-2])}

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t],derivative4[t])
      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])
    }
    else if (model == 5){
      derivative1[t] = 1+ bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + data[t-1]*(data[t-1]>0)
      derivative4[t] = bestPar[2]*derivative4[t-1] - data[t-1]*(data[t-1]<0)
      derivative5[t] = if(t==2){0}else{bestPar[2]*derivative5[t-1] + data[t-2]*(data[t-2]>0)}
      derivative6[t] = if(t==2){0}else{bestPar[2]*derivative6[t-1] - data[t-2]*(data[t-2]<0)}

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t], derivative4[t],derivative5[t], derivative6[t])
      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])
    }
  }

  Sigma = Sigma/obs
  Omega = Omega/obs
  VC = solve(Sigma)%*%Omega%*%solve(Sigma)

  ci_low = bestPar - 1.96*sqrt(diag(VC)/obs)
  ci_up = bestPar + 1.96*sqrt(diag(VC)/obs)

  res = cbind(ci_low, ci_up, bestPar)
  colnames(res) = c('ci_low', 'ci_up', 'c_est')
  truncated = nuisance_truncated
  time = list(t1= tlist_1,t2 = tlist_2)
  if(!predict){
    return(list(res = res, ES = ES,q_hat = q_hat,residuals = residuals, VC = VC,gradient = gradient,pseudo_y = pseudo_y,emp_es = emp_es,time = time,beta_hat = beta_hat,truncated = truncated,denom = denom))
  }else{
    return(list(res = res, ES = ES,q_hat = q_hat,residuals = residuals, VC = VC,gradient = gradient,pseudo_y = pseudo_y,emp_es = emp_es,time = time,beta_hat = beta_hat,truncated = truncated,denom = denom))
  }
}

ESEst_Inf_wls_fair = function(pseudo_y,emp_es, model, pval, data, weights, predict, REP, MAXITER,warmstart,covMethod, truncationMethod,q_hat){
  obs = length(data)

  if(!is.null(warmstart)){
    ES_tmp = get_ES_value_joint(data = data, theta = warmstart, emodel = model,qt0 =0)
    warmstart[1] = warmstart[1] - max(ES_tmp)*(max(ES_tmp)>0)

    initialTargetVectors <- matrix(warmstart, ncol = length(warmstart))
    ESval <- apply(initialTargetVectors, 1, ESObjectiveFunction, 1, model, data, pseudo_y,obs, emp_es, weights, F)
    BestInitialCond <- matrix(initialTargetVectors,nrow = 1)
    ESoptim <- cbind(NA, BestInitialCond, NA)}else{
      emp_es = (obs <= 300)*mean(data[which(data<=q_hat)] ) + (obs>=300)*mean(data[1:300][which(data[1:300]<=q_hat[1:300])])
      if(model %in% c(1,2,3)){
        initialTargetVectors <- matrix(runif(10000*3,-1,1)*(model!=3) + runif(10000*3,0,1)*(model==3)
                                       ,ncol=4*(model==2) + 3*(model %in% c(1,3)))} else if (model == 4){
                                         initialTargetVectors <- matrix(runif(10000*4,-2,2),ncol=4)
                                       }else if (model == 5){
                                         initialTargetVectors <- matrix(runif(10000*6,-2,2),ncol=6)
                                       }
      nInitialCond = 10
      ESval <- apply(initialTargetVectors, 1, ESObjectiveFunction, 1, model, data, pseudo_y,obs, emp_es, weights, F)
      BestInitialCond <- initialTargetVectors[order(ESval), ][1:nInitialCond, ]
      ESoptim <- cbind(NA, BestInitialCond, NA)
    }

  met <- 'L-BFGS-B'
  met_hes <- "Nelder-Mead"
  con <- list(maxit = MAXITER)
  if(model == 1){
    low = c(-Inf, 1e-3, -Inf)
    up = c(-1e-3,1-1e-3,-1e-3)
  }else if(model == 2){
    low = c(-Inf, -1+1e-3, -Inf, -Inf)
    up = c(0,1-1e-3,0,0)
  }else if (model == 3){
    low = c(1e-3,-1+1e-3,0)
    up = c(Inf, 1-1e-3, Inf)
  }else if (model == 4){
    low = c(-Inf,-1+1e-3,-Inf,-Inf)
    up = c(Inf, 1 - 1e-3, Inf, Inf)
  }else if (model == 5){
    low = c(-Inf,-1+1e-3,-Inf,-Inf,-Inf,-Inf)
    up = c(Inf, 1 - 1e-3, Inf, Inf, Inf, Inf)
  }

  for (i in 1:nrow(BestInitialCond)) {
    # initial optimization procedure
    vOptim <- optim(BestInitialCond[i,],
                    ESObjectiveFunction,
                    out = 1,
                    model = model,
                    data = data,
                    pseudo_y = pseudo_y,
                    obs = obs,
                    emp_es = emp_es,
                    weights = weights,
                    method = met,
                    control = con,
                    lower = low,
                    upper = up)

    ESoptim[i, 1] <- vOptim$value
    ESoptim[i, 2:(ncol(initialTargetVectors)+1)] <- vOptim$par



    # repeating optim on hessian optimiziation to get more accurate results
    for(j in 1:REP){
      # once optimizing previous best param (using BFGS)
      if(model == 3){
        theta_tmp<- vOptim$par+ stats::runif(length(vOptim$par), min = -se, max = se)
        theta_tmp[2] <- vOptim$par[2]
        theta_tmp[which(theta_tmp<=0)] <- vOptim$par[which(theta_tmp<=0)]
      }else{
        theta_tmp = vOptim$par + stats::rnorm(length(vOptim$par), mean = 0, sd = se)
        theta_tmp[2] = vOptim$par[2]}
      vOptim <- optim(theta_tmp,
                      ESObjectiveFunction,
                      out = 1,
                      model = model,
                      data = data,
                      pseudo_y = pseudo_y,
                      obs = obs,
                      emp_es = emp_es,
                      weights = weights,
                      method = met_hes,
                      control = con
      )
      # secondly optimizing previous best param with afore method
      if(model == 3){
        theta_tmp<- vOptim$par+ stats::runif(length(vOptim$par), min = -se, max = se)
        theta_tmp[2] <- vOptim$par[2]
        theta_tmp[which(theta_tmp<=0)] <- vOptim$par[which(theta_tmp<=0)]
      }else{
        theta_tmp = vOptim$par + stats::rnorm(length(vOptim$par), mean = 0, sd = se)
        theta_tmp[2] = vOptim$par[2]}
      vOptim <- optim(theta_tmp,
                      ESObjectiveFunction,
                      out = 1,
                      model = model,
                      data = data,
                      pseudo_y = pseudo_y,
                      obs = obs,
                      emp_es = emp_es,
                      weights = weights,
                      method = met,
                      control = con,
                      lower = low,
                      upper = up)
    }
    #checking whether hessian optimization converged
    if(ESoptim[i, 1] - vOptim$value>=0){
      ESoptim[i, 1] <- vOptim$value
      ESoptim[i, 2:(ncol(initialTargetVectors)+1)] <- vOptim$par
    } else {
      ESoptim[i,ncol(ESoptim)] <- j
    }
  }

  bestPar <- ESoptim[order(ESoptim[,1]), ,drop = F][1, 2:(ncol(BestInitialCond)+1)]

  ES = ESObjectiveFunction(bestPar, 2, model, data, pseudo_y, obs, emp_es, weights, predict)
  if(predict == T){ES0 = ES[-length(ES)]}else{ES0 = ES}
  # VC Estimation
  residuals = pseudo_y - ES0
  denom_wls = q_hat - ES
  winsorize_wls = quantile(denom_wls, 0.005)
  denom_wls[which(abs(denom_wls)<=winsorize_wls)] = winsorize_wls
  eta_hat = (data - q_hat)/denom_wls
  nuisance_truncated = truncated_value(eta_hat,truncationMethod,pval)
  if(covMethod == 1){
    conditional_variance = residuals^2
  }else if (covMethod == 2){
    conditional_variance = (q_hat - ES)^2*(1/pval*nuisance_truncated - 1)
  }
  if(nuisance_truncated>=10){
    warning(paste0('Yield extreme truncated quantity estimator! Approach ', covMethod,' can not fit, switch to default method.'))
    conditional_variance = residuals^2
  }
  if(min(conditional_variance)<=0){
    warning(paste0('Yield negative variance estimator! Approach ', covMethod,' can not fit, switch to default method.'))
    conditional_variance = residuals^2
  }
  derivative1 = rep(0,obs)
  derivative2 = rep(0,obs)
  derivative3 = rep(0,obs)
  derivative4 = rep(0,obs)
  derivative5 = rep(0,obs)
  derivative6 = rep(0,obs)

  Sigma = matrix(0, nrow = length(bestPar), ncol = length(bestPar))
  Omega = Sigma
  gradient = matrix(0, nrow = obs, ncol = length(bestPar))

  for( t in 2:obs){
    if(model == 1){
      #SAV
      derivative1[t] = 1 + bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + abs(data[t-1])

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t])

      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])*weights[t]
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])*(weights[t]^2)

    }
    else if(model == 2){
      # AS
      derivative1[t] = 1+ bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + data[t-1]*(data[t-1]>0)
      derivative4[t] = bestPar[2]*derivative4[t-1] - data[t-1]*(data[t-1]<0)

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t], derivative4[t])
      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])*weights[t]
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])*(weights[t]^2)
    }
    else if (model == 3){
      #GARCH
      derivative1[t] = (1 + 2*bestPar[2]*ES[t-1]*derivative1[t-1]) / (2*ES[t])
      derivative2[t] = (ES[t-1]^2 + 2*bestPar[2]*ES[t-1]*derivative2[t-1]) / (2*ES[t])
      derivative3[t] = (2*bestPar[2]*ES[t-1]*derivative3[t-1] + data[t-1]^2) / (2*ES[t])

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t])

      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])*weights[t]
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])*(weights[t]^2)
    }
    else if (model == 4){
      derivative1[t] = 1 + bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + abs(data[t-1])
      derivative4[t] = if(t==2){0}else{bestPar[2]*derivative4[t-1] + abs(data[t-2])}

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t],derivative4[t])
      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])*weights[t]
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])*(weights[t])^2
    }
    else if (model == 5){
      derivative1[t] = 1+ bestPar[2]*derivative1[t-1]
      derivative2[t] = ES[t-1] + bestPar[2]*derivative2[t-1]
      derivative3[t] = bestPar[2]*derivative3[t-1] + data[t-1]*(data[t-1]>0)
      derivative4[t] = bestPar[2]*derivative4[t-1] - data[t-1]*(data[t-1]<0)
      derivative5[t] = if(t==2){0}else{bestPar[2]*derivative5[t-1] + data[t-2]*(data[t-2]>0)}
      derivative6[t] = if(t==2){0}else{bestPar[2]*derivative6[t-1] - data[t-2]*(data[t-2]<0)}

      gradient[t,] = c(derivative1[t], derivative2[t], derivative3[t], derivative4[t],derivative5[t], derivative6[t])
      Sigma = Sigma + gradient[t,]%*%t(gradient[t,])*weights[t]
      Omega = Omega + conditional_variance[t]*gradient[t,]%*%t(gradient[t,])*(weights[t])^2
    }
  }

  Sigma = Sigma/obs
  Omega = Omega/obs
  VC = (solve(Sigma)%*%Omega%*%solve(Sigma))

  ci_low = bestPar - 1.96*sqrt(diag(VC)/obs)
  ci_up = bestPar + 1.96*sqrt(diag(VC)/obs)

  res = cbind(ci_low, ci_up, bestPar)
  colnames(res) = c('ci_low', 'ci_up', 'c_est')
  truncated = nuisance_truncated
  if(!predict){
    return(list(res = res, ES = ES,residuals = residuals, VC = VC,gradient = gradient,Sigma = Sigma,Omega = Omega,truncated = truncated))
  }else{
    return(list(res = res, ES = ES,residuals = residuals, VC = VC,gradient = gradient,Sigma = Sigma,Omega = Omega,truncated = truncated))
  }
}

get_ES_value_joint = function(data,theta,emodel,qt0){
  obs = length(data)
  emp_es = qt0 - 0.5 #just to make sure strictly smaller

  ES = rep(0,obs)
  ES[1] = emp_es

  ES <- if (emodel == 1) {
    SAV_ES(theta, data, ES[1], ES, obs, F)
  } else if (emodel == 2) {
    AS_ES(theta, data, ES[1], ES, obs, F)
  } else if (emodel == 3) {
    GARCH_ES(theta, data, ES[1], ES, obs, F)
  } else if (emodel == 4){
    SAV12_ES(theta, data, ES[1], ES, obs, F)
  } else if (emodel == 5){
    AS12_ES(theta, data, ES[1], ES, obs, F)
  }
  return(ES)
}
