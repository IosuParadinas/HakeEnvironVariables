# Script Info -------------------------------------------------------------


# Script Info -------------------------------------------------------------

# Function to perform forward model selection based on AIC for scam models with concave terms (IM-20-MISSION)

# Function to perform forward model selection based on AIC for scam models with concave terms (IM-20-MISSION)

# The function arguments are: 
#   basef   : base model formula, usually  y ~ 1
#   vars    : list of potential variables to be included
#   dat     : data set
#   vmax    : maximum number of variables to be included. If NULL, there is no upper limit
#   k_cross : number of cross validation sets
#   p_train : percent of data used in training datasets
#   k       : knots
#   sp      : either the fixed smoothing parameter or NULL (if sp are estimated) 
#   t.test.alpha : model selection is performed based on a paired t.test performed across the obtained pseudoR2, MAE and RMSE at eack k_cross subset. 
#   N_test_positive : required minimim number of t.test below t.test.alpha to accept the new model (1, 2 or 3)
# Iosu Paradinas (based on code from Leire Ibaibarriaga)
# Last modifications 17/3/2023


CrossVal.scam <- function(basef, vars,  dat, vmax=NULL, k=10, sp=NULL, set_cross =TRUE,
                          k_cross = 10, p_train = 0.8,t.test.alpha=0.05,N_test_positive=2,
                          base_n_scams = 0){
  require("scam")
  # initialise the vector of selected variables
  svars <- NULL 
  # initialise the list of selected models at each step
  smod <- list() 
  smod[[1]] = basef
  
  if(set_cross == T){
  k_cross = k_cross
  random_sample <- createDataPartition(1:nrow(dat),p = p_train, list = FALSE,times = k_cross)
  }else{
    k_cross = length(partitions)
  }
  # maximum number of variables to be included
  vmax <- min(vmax, length(vars)) 
  # base model 
  for(crossV in 1:k_cross){
    # creating training data as 80% of the dataset
    training_dataset  <- dat[random_sample[,crossV], ]
    testing_dataset <- dat[-random_sample[,crossV], ]
    
    m0 <- scam(formula(basef), family=binomial(link="logit"), data=training_dataset)
    #m0 <- scam(formula(basef), family=poisson(link="log"), data=training_dataset)
    
    predictions0 <- predict(m0, newdata = testing_dataset,type="response",)
    
    # computing model performance metrics
    perf0 = data.frame( R2 = R2(predictions0, testing_dataset$y),
                        RMSE = RMSE(predictions0, testing_dataset$y),
                        MAE = MAE(predictions0, testing_dataset$y))
    if(crossV==1){performance0 = perf0}else{performance0 = rbind(performance0,perf0)}
    
  } # final model 0
  
  #mean_perf0 = apply(performance0,2,mean)
  
  # loop: add one variable at each step 
  for (i in 1:vmax){
    
    print(paste("Fitting models with", i, "terms"))
    print(svars)
    perf.list <- list() # list of models with i terms
    for (j in 1:length(vars)){
      if (! vars[j] %in% svars){
        tempvars <- c(svars, vars[j])
        formu <- paste(c(basef, paste0("s(",tempvars,",m=2,bs='cv',k=",k,")")), collapse="+")
        cat("\t", formu, "\n")
        if (is.null(sp)){
          
          for(crossV in 1:k_cross){
            training_dataset  <- dat[random_sample[,crossV], ]
            testing_dataset <- dat[-random_sample[,crossV], ]
            
            m = try(scam(formula(formu), family=binomial(link="logit"), data=training_dataset,sp=NULL), silent=T)
            #m = try(scam(formula(formu), family=poisson(link="log"), data=training_dataset,sp=NULL), silent=T)
            
            # print message if model with sp NULL is not working 
            if (any(class(m)=="try-error")){
              print(paste("Model selection did not work for species SP", i))
              
              mysp <- 0.00001
                            #fit the model with assigned sp
              m <- scam(formula(formu), family=binomial(link="logit"), data=training_dataset,sp=rep(mysp,length(tempvars)+base_n_scams)) # by default aic.tol=
              #m <- scam(formula(formu), family=poisson(link="log"), data=training_dataset,sp=rep(mysp,length(tempvars))) # by default aic.tol=
            }
            
            predictions <- predict(m, newdata = testing_dataset,type="response")
            
            # computing model performance metrics
            perf = data.frame( R2 = R2(predictions, testing_dataset$y),
                               RMSE = RMSE(predictions, testing_dataset$y),
                               MAE = MAE(predictions, testing_dataset$y))
            if(crossV==1){performance = perf}else{performance = rbind(performance,perf)}
            
          }
        }else{
          for(crossV in 1:k_cross){
            
            training_dataset  <- dat[random_sample[,crossV], ]
            testing_dataset <- dat[-random_sample[,crossV], ]
            
            m <- scam(formula(formu), family=binomial(link="logit"), data=training_dataset, sp=rep(sp,length(tempvars)+base_n_scams))
            #m <- scam(formula(formu), family=poisson(link="log"), data=training_dataset, sp=rep(sp,length(tempvars)))
            predictions <- predict(m, newdata = testing_dataset,type="response")
            
            # computing model performance metrics
            perf = data.frame( R2 = R2(predictions, testing_dataset$y),
                               RMSE = RMSE(predictions, testing_dataset$y),
                               MAE = MAE(predictions, testing_dataset$y))
            if(crossV==1){performance = perf}else{performance = rbind(performance,perf)}
            
          }
        }
        perf.list[[vars[j]]] <- performance
      }
    }
    
    ##### mirar de todos lasasas
    ea = t(sapply(perf.list, function(x){apply(x,2,mean)})) 
    metrics = which.max(ea[,1]);which.min(ea[,2]);which.min(ea[,3]) 
    predictive_best = round(mean(metrics))
    
    perf.list[[predictive_best]]
    
    ### maybe t-test for comparing the two groups?
    
    R2_test = t.test(performance0$R2,perf.list[[predictive_best]]$R2 , paired = T,alternative = "less")
    mae_test = t.test(performance0$MAE,perf.list[[predictive_best]]$MAE , paired = T,alternative = "greater")
    rmse_test = t.test(performance0$RMSE,perf.list[[predictive_best]]$RMSE , paired = T,alternative = "greater")
    
    if(sum(c(R2_test$p.value,mae_test$p.value,rmse_test$p.value)<t.test.alpha)>=N_test_positive){
      #aic.min <- aic.vec[which.min(aic.vec)]
      #if ( (aic0-unname(aic.min)) > aic.tol){
      svars <- c(svars, rownames(ea)[predictive_best]) # selected vars with best performance
      #svars <- c(svars, names(aic.min)) # selected vars with smallest AIC
      
      smod[[i+1]] <- paste(c(smod[[i]], paste0("s(",row.names(ea)[predictive_best],",m=2,bs='cv',k=",k,")")), collapse="+") # selected model
      performance0 <- perf.list[[predictive_best]]
      # smod[[i+1]] <- m.list[[names(aic.min)]] # selected model
      # aic0 <- unname(aic.min)
    }else{
      print(paste("The more complex models did not improve the current one"))
      break
    }
    
  }
  
  final_model <- scam(formula(smod[[length(smod)]]), family=binomial(link="logit"), data=dat)
  #final_model <- scam(formula(smod[[length(smod)]]), family=poisson(link="log"), data=dat)
  
  names(smod) <- c("null", svars)
  out <- list(smod=smod, svars=svars,final_model=final_model)
  return(out)
  
}

