# Labor Exchange Simulation
lex <- function( N = 12,            # N individuals
                 k = 1,             # lag length
                 groups = 5,        # N kin groups
                 tsteps = 10,       # iterations
                 alpha = -2,        # log odds base rate
                 betaK =  2,        # kin effect on dyads
                 betaS =  0,        # skill/knowledge effect on dyads
                 betaH =  1.2,      # household size effect
                 betaW = -1.5,      # wealth effect
                 betaL,             # lag effects
                 betaV =  0,        # village effects 
                 betaR =  0,        # reputation effects
                 return_data = F,   # return xlist, d, and v
                 choose_egos = F,   # manually select egos (T/F)
                 egos,              # vector of ego ids
                 Vweights = c(0.5,0.5),
                 burn_in = 10,
                 add_defectors = F, 
                 N_defectors = 1, 
                 def_strength = -10,
                 def_rate = 1,
                 seed  = 777 )                   

{  
    require(tidyverse)
    # functions
    logit <- function(p) log(p/(1-p))
    inv_logit <- function(x) 1/(1+exp(-x))
    
    # DATA STRUCTURE
    # create
    d <- data.frame( t(combn(1:N, 2)) )
    colnames(d) <- c('i','j')
    N_dyads <- nrow(d)
    d$did <- 1:N_dyads
    d$exchange <- NA
    d$pool <- NA
    
    # create attributes
    set.seed(seed)
    v <- data.frame(
        id = 1:N,
        hhsize = scale( rpois(N,2)+1 ),
        kgroup = sample(LETTERS[1:groups], N, replace = T), 
        wealth = scale( rlnorm(N,0,3) ), 
        skills = scale( rlnorm(N,0,1)+1 ), 
        reputa = scale( rbinom(N, 5, prob=0.5)), 
        villag = sample(LETTERS[c(3,7)], N, replace = T, prob = Vweights), 
        defect = 0
    )
    
    # defection 
    if( add_defectors == T ) { v[ sample(1:N, size = N_defectors, replace = F), 'defect'] <- def_strength }
    
    # kin 
    d1 <- merge(d,  v, by.x = 'j', by.y = 'id')
    d  <- merge(d1, v, by.x = 'i', by.y = 'id')
    d$kin   <- as.integer(d$kgroup.x == d$kgroup.y)        # same kin group
    d$res   <- as.integer(d$villag.x == d$villag.y)        # village co residence
    d$skill.diff <- d$skills.x + d$skills.y
    
    # defective dyads 
    d$defect <- d$defect.x + d$defect.y
    
    # rename 
    d <- d %>% select(i, j, did, exchange, pool, kin, res, skill.diff, defect,
                      hhsize.i = hhsize.y, hhsize.j = hhsize.x, 
                      kgroup.i = kgroup.y, kgroup.j = kgroup.x, 
                      wealth.i = wealth.y, wealth.j = wealth.x, 
                      skills.i = skills.y, skills.j = skills.x, 
                      reputa.i = reputa.y, reputa.j = reputa.x, 
                      villag.i = villag.y, villag.j = villag.x, 
                      defect.i = defect.y, defect.j = defect.x
                      )
    
    # defective dyads 
    #d$defect <- d$defect.i + defect.j
    
    # set parameters
    alpha <- alpha
    betaH <- betaH
    betaL <- betaL
    betaK <- betaK
    betaS <- betaS
    betaR <- betaR
    betaV <- betaV
    betaW <- betaW
    
    # create a list
    xlist <- list()
    
    # INITIALIZATION
    # sample ego
    for(b in 1:burn_in) {
        
        if ( choose_egos == T ) { ego <- sample(egos, 1) }
        else if ( choose_egos == F) { ego <- sample(1:N, 1) }
        
        # tag eligible
        for(i in 1:nrow(d)) {
            x <- as.numeric(ego %in% d[i, 1:2])
            d[i, 'pool'] <- x 
        }
        
        # initialize exchanges 
        d$exchange <- rbinom(N_dyads, 1, d$pool * inv_logit(alpha + 
                                                                betaK*d$kin + 
                                                                betaV*d$res +
                                                                betaS*d$skill.diff + 
                                                                betaH*v$hhsize[ego] + 
                                                                betaW*v$wealth[ego] + 
                                                                betaR*v$reputa[ego] ))
        # track time, ego, probability
        d$t <- b
        d$ego <- ego
        d$prob <- d$pool * inv_logit(alpha + 
                                         betaK*d$kin + 
                                         betaV*d$res +
                                         betaS*d$skill.diff + 
                                         betaH*v$hhsize[ego] + 
                                         betaW*v$wealth[ego] + 
                                         betaR*v$reputa[ego] )
        d$phase <- 'burn'
        d$lag <- NA
        
        xlist[[b]] <- d 
    }
    
    # GENERATE EXCHANGE PROBABILITIES
    # run loop
    start <- burn_in+1
    steps <- burn_in+tsteps
    
    
    for(t in start:steps) {
        
        if ( choose_egos == T ) { ego <- sample(egos, 1) }
        else if ( choose_egos == F) { ego <- sample(1:N, 1) }
        
        # reset
        d$exchange <- NA
        d$pool <- NA
        
        for(i in 1:nrow(d)) {
            x <- as.numeric(ego %in% d[i, 1:2])
            d[i, 'pool'] <- x 
        }
        
        d$exchange <- rbinom(n = N_dyads, 
                             size = 1, 
                             prob = d$pool * inv_logit( alpha + 
                                                            def_rate*d$defect + 
                                                            betaK*d$kin + 
                                                            betaV*d$res +
                                                            betaS*d$skill.diff + 
                                                            betaH*v$hhsize[ego] + 
                                                            betaW*v$wealth[ego] + 
                                                            betaR*v$reputa[ego] +
                                                            betaL*apply(matrix(unlist(lapply(xlist[(t-1):(t-k)], '[[', 'exchange')), ncol=length((t-1):(t-k))), 1, sum)  ))
        # track time, ego, and probability
        d$t <- t
        d$ego <- ego
        d$prob <- d$pool * inv_logit( alpha + 
                                          def_rate*d$defect + 
                                          betaK*d$kin + 
                                          betaV*d$res +
                                          betaS*d$skill.diff + 
                                          betaH*v$hhsize[ego] + 
                                          betaW*v$wealth[ego] + 
                                          betaR*v$reputa[ego] +
                                          betaL*apply(matrix(unlist(lapply(xlist[(t-1):(t-k)], '[[', 'exchange')), ncol=length((t-1):(t-k))), 1, sum) )
        d$phase <- 'sim'
        d$lag <- apply(matrix(unlist(lapply(xlist[(t-1):(t-k)], '[[', 'exchange')), ncol=length((t-1):(t-k))), 1, sum) 
        
        xlist[[t]] <- d 
    }
    
    if( return_data == F) { return(bind_rows(xlist)) } else { return(list(bind_rows(xlist), v, d)) }
    
}

# Example simulation: kinship and lag effects only 
SKL <- lex( N = 25,
            k =  4,
            groups  =   6,
            tsteps  =  50, 
            burn_in =  10, 
            alpha   =  -2, 
            betaK   =   2,
            betaL   =   2,
            betaR   =   1.2,
            betaS   =  -0.5,
            betaV   =   0,
            betaH   =   0.5, 
            betaW   =   0, 
            seed    =  27, 
            return_data = T,
            add_defectors = F, 
            N_defectors = 5)


SKL[[1]] %>% filter(pool==1 & phase=='sim') 


#### Three Scenarios: Kinship + Village; HHsize + Reputation; Defection 

# In scenario A we include strong effects of kinship and village.
simA <- lex(N = 25, 
            k = 1,
            groups =   6, 
            tsteps = 400, 
            alpha  =  -2, 
            betaK  =   2, 
            betaV  =   2, 
            betaS  =   0, 
            betaH  =   0, 
            betaR  =   0, 
            betaW  =   0, 
            betaL  =   0, 
            add_defectors = F, 
            return_data = T, 
            seed = 27,
            Vweights = c(0.4,0.6)) 

# In scenario B, we leave out the dyadic effects and only include individual effects. In this case, I am going to focus on reputation and household size, to keep these simple. This will serve as a sort NULL model. 
simB <- lex(N = 25, 
            k = 1,
            groups =   6, 
            tsteps = 400, 
            alpha  =  -2, 
            betaK  =   0, 
            betaV  =   0, 
            betaS  =   0, 
            betaH  =   2, 
            betaR  =   2, 
            betaW  =   0, 
            betaL  =   0, 
            add_defectors = F, 
            return_data = T, 
            seed = 27,
            Vweights = c(0.4,0.6)) 


# And the final C scenario will include the same dyadic effects as simA, but will include 3 pure defectors with a def rate of 1. 
simC <- lex(N = 25, 
            groups =   6, 
            tsteps = 400, 
            alpha  =  -2, 
            betaK  =   2, 
            betaV  =   2, 
            betaS  =   0, 
            betaH  =   0, 
            betaR  =   0, 
            betaW  =   0, 
            betaL  =   0, 
            add_defectors = T,
            N_defectors = 3, 
            def_strength = -6, 
            def_rate = 1,
            return_data = T, 
            Vweights = c(0.4,0.6)) 



