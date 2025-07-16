simprob <- function(ndose, targetE, targetT, u1, u2, randomtype) {
  if (OBD == 0) {
    obd <-sample(ndose, 1)
  } else {
    obd <- OBD
  }
  mtd <- obd # initialization
  prob_plateau <- 0.5
  obd.temp <- 0
  jj <- kk <- rep(0, ndose)
  uu <- runif(1)
  if (uu < prob_plateau && obd < ndose) {
    # generate plateau cases
    while(obd.temp != obd || kk[obd] > (targetT + 0.05) || jj[obd] < (targetE + 0.05)) {
      mtd <- (obd - 1) + sample(ndose - obd + 1, 1)
      D <- 1:ndose
      if (mtd < ndose) {
        temp <- sort(runif(length(D) - mtd, targetT, 1))
        bornesup <- max(targetT + 0.10, temp[length(D) - mtd])
      }
      if (mtd == ndose) {
        bornesup <- targetT + (1 - targetT) * rbeta(1, 0.5, 1)
      }
      mtd.temp <- 0
      while (mtd.temp != mtd) {
        kk <- sort(runif(ndose, 0, bornesup))
        mtd.temp <- which.min(abs(kk - targetT))
      }
      bornesup <- runif(1, targetE + 0.10, 1.0)
      #med=obd
      #bornesup<-max(runif(med,targetE+0.10,0.9))
      jj[1:obd] <- sort(runif(obd, 0, bornesup))
      jj[obd:ndose] <- jj[obd]
      utility <- u1 * jj + (1 - kk) * u2
      obd.temp <- which.max(utility[1:mtd])
    }
  }
  if (uu >= prob_plateau && obd < ndose) {
    # generate non-plateau cases 
    while (obd.temp != obd || kk[obd] > (targetT + 0.05) || jj[obd] < (targetE + 0.05)) {
      mtd <- (obd - 1) + sample(ndose - obd + 1, 1)
      D <- 1:ndose
      if (mtd < ndose) {
        temp <- sort(runif(length(D) - mtd, targetT, 1))
        bornesup <- max(targetT + 0.10,temp[length(D) - mtd])
      }
      if(mtd == ndose) {
        bornesup <- targetT + (1 - targetT) * rbeta(1, 0.5, 1)
      }
      mtd.temp <- 0
      while(mtd.temp != mtd) {
        kk <- sort(runif(ndose, 0, bornesup))
        mtd.temp <- which.min(abs(kk - targetT))
      }
      med <- (obd - 1) + sample(ndose - obd + 1, 1)
      bornesup <- runif(1, targetE + 0.10, 1.0)
      #bornesup<-max(runif(med,targetE+0.10,0.9))
      jj[med] <- max(runif(ndose, 0, bornesup))
      if(med > 1) {
        jj[1:(med - 1)] <- sort(runif(med - 1, 0, jj[med]))
      }
      if(med < ndose) {
        jj[(med + 1):ndose] <- sort(runif(ndose - med, 0, jj[med]), decreasing = TRUE)
      }
      utility <- u1 * jj + (1 - kk) * u2
      obd.temp <- which.max(utility[1:mtd])
    }
  }
  if (obd == ndose) {
    while (obd.temp != obd || kk[obd] > (targetT + 0.05) || jj[obd] < (targetE + 0.05)) {
      mtd <- (obd - 1) + sample(ndose - obd + 1, 1)
      D <- 1:ndose
      if (mtd < ndose) {
        temp <- sort(runif(length(D) - mtd, targetT, 1))
        bornesup <- max(targetT + 0.10, temp[length(D) - mtd])
      }
      if (mtd == ndose) {
        bornesup <- targetT + (1 - targetT) * rbeta(1, 0.5, 1)
      }
      mtd.temp <- 0
      while(mtd.temp != mtd) {
        kk <- sort(runif(ndose, 0, bornesup))
        mtd.temp <- which.min(abs(kk - targetT))
      }
      bornesup <- runif(1, targetE + 0.10, 1.0)
      #med=ndose
      #bornesup<-max(runif(med,targetE+0.10,0.9))
      jj <- sort(runif(ndose, 0, bornesup))
      utility <- u1 * jj + (1 - kk) * u2
      obd.temp <- which.max(utility[1:mtd])
    }
  }
  return(list(pE = jj, pT = kk, obd = obd, mtd = mtd))
}

peestimate <- function(yE, n) {
  ndose <- length(yE)
  lik <- rep(0, ndose)
  pe <- (yE + 0.05) / (n + 0.1)
  p.e <- matrix(NA, ndose, ndose)

  for (i in 1:ndose) {
    if (i == 1) {
      x <- seq(ndose, 1, by = -1)
    } else {
      x <- c(1:(i - 1), seq(ndose, i))
    }

    p.e[i, ] <- ufit(pe, lmode = i, x = x, w = n + 0.5)[[2]]
    lik[i] <- prod(dbinom(yE, n, p.e[i, ]))
  }

  lik <- lik / sum(lik)
  pe <- t(p.e) %*% lik + 0.01 * seq(1, ndose)

  return(pe)
}

pava <- function(x, wt = rep(1, length(x))) {
  n <- length(x)

  if (n <= 1) {
    return(x)
  }

  if (any(is.na(x)) || any(is.na(wt))) {
    stop("Missing values in 'x' or 'wt' not allowed")
  }

  lvlsets <- 1:n

  repeat {
    viol <- as.vector(diff(x)) < 0

    if (!any(viol)) {
      break
    }

    i <- min((1:(n - 1))[viol])
    lvl1 <- lvlsets[i]
    lvl2 <- lvlsets[i + 1]
    ilvl <- (lvlsets == lvl1 | lvlsets == lvl2)

    x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])
    lvlsets[ilvl] <- lvl1
  }

  x
}

get.boundary <- function(target, targetE, ncohort, cohortsize = 3,
                         p.saf = NA, p.tox = NA, cutoff.eli = 0.95, 
                         cutoff.eli.E = 0.90) {

  # if the user does not provide p.saf and p.tox, use the default values
  if (is.na(p.saf)) p.saf <- 0.6 * target
  if (is.na(p.tox)) p.tox <- 1.4 * target

  ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
  npts <- ncohort * cohortsize
  ntrt <- NULL
  b.e <- NULL
  b.d <- NULL
  elim <- NULL
  elimE <- NULL

  for (n in (1:ncohort) * cohortsize) {
    error.min <- 3
    for (m1 in 0:(n - 1)) {
      for (m2 in (m1 + 1):n) {
        error1 <- pbinom(m1, n, target) + 1 - pbinom(m2 - 1, n, target)
        error2 <- 1 - pbinom(m1, n, p.saf)
        error3 <- pbinom(m2 - 1, n, p.tox)
        error <- error1 + error2 + error3

        if (error < error.min) {
          error.min <- error
          cutoff1 <- m1
          cutoff2 <- m2
        }
      }
    }

    ntrt <- c(ntrt, n)
    b.e <- c(b.e, cutoff1)
    b.d <- c(b.d, cutoff2)

    elimineed <- 0  # indicating whether elimination is needed
    elimineedE <- 0
    if (n < 3) {
      elim <- c(elim, NA)
      elimE <- c(elimE, NA)  # require treating at least 3 patients before eliminating a dose
    } else {
      for (ntox in 3:n) {  # determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        if (1 - pbeta(target, ntox + 1, n - ntox + 1) > cutoff.eli) {
          elimineed <- 1
          break
        }
      }

      if (elimineed == 1) {
        elim <- c(elim, ntox)
      } else {
        elim <- c(elim, NA)  # set the elimination boundary large such that no elimination will actually occurs
      }

      for (neff in n:0) {
        if (pbeta(targetE, neff + 1, n - neff + 1) > cutoff.eli.E) {
          elimineedE <- 1
          break
        }
      }

      if (elimineedE == 1) {
        elimE <- c(elimE, neff)
      } else {
        elimE <- c(elimE, NA)
      }
    }
  }

  for (i in 1:length(b.d)) {
    if (!is.na(elim[i]) && (b.d[i] > elim[i])) {
      b.d[i] <- elim[i]
    }
  }

  boundaries <- rbind(ntrt, elim, b.d, b.e, elimE)
  rownames(boundaries) <- c("Number of patients treated", "Eliminate if # of DLT >=",
                            "Deescalate if # of DLT >=", "Escalate if # of DLT <=", "Eliminate if # of Eff <=")
  colnames(boundaries) <- rep("", ncohort)

  return(boundaries)
}

f1 <- function(x, bn.m1, bn.m2, rho) {
  ff <- dnorm(x, bn.m1, 1) * (1 - pnorm(0, bn.m2 + rho * (x - bn.m1), sqrt(1 - rho^2)))
  return(ff)
}

ji3plus3 <- function(x, y, n, pT, eps1, eps2, pE) {
  if ((x / n < pT - eps1) && (y / n <= pE)) {
    dec <- "E"
  } else if ((x / n < pT - eps1) && (y / n > pE)) {
    dec <- "S"
  } else if ((x / n >= pT - eps1 && x / n <= pT + eps2) && (y / n <= pE)) {
    dec <- "E"
  } else if ((x / n >= pT - eps1 && x / n <= pT + eps2) && (y / n > pE)) {
    dec <- "S"
  } else if ((x / n > pT + eps2) && ((x - 1) / n < pT - eps1)) {
    dec <- "S"
  } else {
    dec <- "D"
  }
  return(dec)
}

betavar <- function(a, b) {
  resp <- a * b / ((a + b) ^ 2 * (a + b + 1))
  return(resp)
}

	## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
stein.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95) {
# if the user does not provide p.saf and p.tox, use the default values
		if(is.na(p.saf)) p.saf=0.75*target;
		if(is.na(p.tox)) p.tox=1.25*target;

### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
		npts = ncohort*cohortsize;
		ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
		for(n in (1:ncohort)*cohortsize)
		{
			error.min=3;
			for(m1 in 0:(n-1))
			{
				for(m2 in (m1+1):n)
				{
 
						error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
						error2 = 1-pbinom(m1, n, p.saf);
						error3 = pbinom(m2-1, n, p.tox);
					
					error=error1+error2+error3;
					if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
				}
			}
			ntrt = c(ntrt, n);
			b.e = c(b.e, cutoff1);
			b.d = c(b.d, cutoff2);
			
			elimineed=0; # indicating whether elimination is needed
			if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
			else
			{
				for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
				{
					if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
				}
				if(elimineed==1) { elim = c(elim, ntox); }
				else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
			}
		}
		for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
		boundaries = rbind(ntrt, elim, b.d, b.e);
		rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
								 "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
		colnames(boundaries) = rep("", ncohort);
		
		return(boundaries);
	}

tepi.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95) {
    
    # if the user does not provide p.saf and p.tox, use the default values
    if(is.na(p.saf)) p.saf=0.6*target;
    if(is.na(p.tox)) p.tox=1.4*target;
    
    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      error.min=3;
      for(m1 in 0:(n-1))
      {
        for(m2 in (m1+1):n)
        {
          
          error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
          error2 = 1-pbinom(m1, n, p.saf);
          error3 = pbinom(m2-1, n, p.tox);
          
          error=error1+error2+error3;
          if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
        }
      }
      ntrt = c(ntrt, n);
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);
      
      elimineed=0; # indicating whether elimination is needed
      if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
    colnames(boundaries) = rep("", ncohort);
    
    return(boundaries);
  }

  ## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
utpi.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95)
  {
    
    # if the user does not provide p.saf and p.tox, use the default values
    if(is.na(p.saf)) p.saf=0.6*target;
    if(is.na(p.tox)) p.tox=1.4*target;
    
    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      error.min=3;
      for(m1 in 0:(n-1))
      {
        for(m2 in (m1+1):n)
        {
          
          error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
          error2 = 1-pbinom(m1, n, p.saf);
          error3 = pbinom(m2-1, n, p.tox);
          
          error=error1+error2+error3;
          if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
        }
      }
      ntrt = c(ntrt, n);
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);
      
      elimineed=0; # indicating whether elimination is needed
      if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
    colnames(boundaries) = rep("", ncohort);
    
    return(boundaries);
  }
  ## to make get.oc as self-contained function, we copied functions get.boundary() and select.mtd() here.
boinet.boundary <- function(target, ncohort, cohortsize=3, p.saf=NA, p.tox=NA,  cutoff.eli=0.95)
  {
    
    # if the user does not provide p.saf and p.tox, use the default values
    if(is.na(p.saf)) p.saf=0.6*target;
    if(is.na(p.tox)) p.tox=1.4*target;
    
    ### numerical search for the boundaries that minimize decision errors of dose escalation/deescalation
    npts = ncohort*cohortsize;
    ntrt=NULL; b.e=NULL; b.d=NULL; elim=NULL;
    for(n in (1:ncohort)*cohortsize)
    {
      error.min=3;
      for(m1 in 0:(n-1))
      {
        for(m2 in (m1+1):n)
        {
          
          error1 = pbinom(m1, n, target)+1-pbinom(m2-1, n, target);
          error2 = 1-pbinom(m1, n, p.saf);
          error3 = pbinom(m2-1, n, p.tox);
          
          error=error1+error2+error3;
          if(error<error.min) {error.min=error; cutoff1=m1; cutoff2=m2;}
        }
      }
      ntrt = c(ntrt, n);
      b.e = c(b.e, cutoff1);
      b.d = c(b.d, cutoff2);
      
      elimineed=0; # indicating whether elimination is needed
      if(n<3) { elim = c(elim, NA); }  # require treating at least 3 patients before eliminating a dose
      else
      {
        for(ntox in 3:n) #determine elimination boundary, prior beta(1,1) is used in beta-binomial model
        {
          if(1-pbeta(target, ntox+1, n-ntox+1)>cutoff.eli) {elimineed=1; break;}
        }
        if(elimineed==1) { elim = c(elim, ntox); }
        else { elim = c(elim, NA); } # set the elimination boundary large such that no elimination will actually occurs
      }
    }
    for(i in 1:length(b.d)) { if(!is.na(elim[i]) && (b.d[i]>elim[i])) b.d[i]=elim[i]; }
    boundaries = rbind(ntrt, elim, b.d, b.e);
    rownames(boundaries) = c("Number of patients treated", "Eliminate if # of DLT >=", 
                             "Deescalate if # of DLT >=",  "Escalate if # of DLT <=");
    colnames(boundaries) = rep("", ncohort);
    
    return(boundaries);
  }

pos<-function(pE,yE,yT,n,u1,u2,u){
      f<-dbeta(pE,yE+1,n-yE+1)
      f<-f*pbeta((u1*pE+u2-u)/u2,yT+1,n-yT+1)
      return(f)
  }


upmgen <- function(tbound,ebound,parabeta,n,ntox,neff){
  upm <- matrix(0,length(tbound)-1,length(ebound)-1)
  for (i in 1:(length(tbound)-1)){
    for (j in 1:(length(ebound)-1)){
      upm[i,j]=( pbeta(tbound[i+1],parabeta[1]+ntox,parabeta[2]+n-ntox) - 
                   pbeta(tbound[i],parabeta[1]+ntox,parabeta[2]+n-ntox) ) * 
        ( pbeta(ebound[j+1],parabeta[1]+neff,parabeta[2]+n-neff) - 
            pbeta(ebound[j],parabeta[1]+neff,parabeta[2]+n-neff) ) / 
        ((tbound[i+1]-tbound[i])*(ebound[j+1]-ebound[j]))
    }
  }
  return(upm)
}

maptoDEC <- function(index) {
  if (index %in% c("UU", "UL")) {
    DEC <- "D"
  } else if (index %in% c("EL", "EU", "LU")) {
    DEC <- "S"
  } else {
    DEC <- "E"
  }
  return(DEC)
}

upmgen <- function(tbound,ebound,parabeta,n,ntox,neff){
  upm <- matrix(0,length(tbound)-1,length(ebound)-1)
  for (i in 1:(length(tbound)-1)){
    for (j in 1:(length(ebound)-1)){
      upm[i,j]=( pbeta(tbound[i+1],parabeta[1]+ntox,parabeta[2]+n-ntox) - 
                   pbeta(tbound[i],parabeta[1]+ntox,parabeta[2]+n-ntox) ) * 
        ( pbeta(ebound[j+1],parabeta[1]+neff,parabeta[2]+n-neff) - 
            pbeta(ebound[j],parabeta[1]+neff,parabeta[2]+n-neff) ) / 
        ((tbound[i+1]-tbound[i])*(ebound[j+1]-ebound[j]))
    }
  }
  return(upm)
}

maptoDEC <- function(index) {
  if (index %in% c("UU", "UL")) {
    DEC <- "D"
  } else if (index %in% c("EL", "EU", "LU")) {
    DEC <- "S"
  } else {
    DEC <- "E"
  }
  return(DEC)
}