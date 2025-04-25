# solver for generally non-linear root functions
deSolverNonLinRootFinder <- function(z,
                                     z.new,
                                     k1,
                                     k4,
                                     t,
                                     h,
                                     parms,
                                     root.func,
                                     init.subdiv = 8L,
                                     last.root = NULL,
                                     last.root.offset = 1.0e-4,
                                     tol_uniroot = 1e-8) {
  # sets up interpolation polynomial (hermite, check Hairer 1993, equation 6.7)
  ff3 <- h * (k1 + k4) - 2.0 * (z.new - z)
  ff2 <- 3.0 * (z.new - z) - h * (k4 + 2.0 * k1)
  ff1 <- k1 * h
  r.grid <- seq(from = 0,
                to = 1,
                length.out = init.subdiv)
  rf.left <- root.func(t, z, parms)
  left.modified <- -1L
  if (!is.null(last.root)) {
    if (last.root$which.root > 0 &&
        last.root$root.type == 0) {
      # if there is a root in previous last step and it is due to a non-linear root function
      left.modified <- last.root$which.root
      rf.left[last.root$which.root] <- root.func(
        t + last.root.offset * h,
        z + last.root.offset * ff1
        + last.root.offset ^ 2 * ff2
        + last.root.offset ^ 3 * ff3,
        parms
      )[last.root$which.root] # offset the coordinate corresponding to the last root by a small amount in order to not get the initial point as the root
    }
  }
  i.done <- -1L
  #print(rf.left)
  for (i in 2:init.subdiv) {
    # divide into many intervals and check if root exists in the bracketed interval
    r <- r.grid[i]
    rf.right <-
      root.func(t + h * r, z + r * ff1 + r ^ 2 * ff2 + r ^ 3 *
                  ff3, parms)
    rp <- rf.left * rf.right
    #print(c(rf.left,rf.right))
    if (any(rp < 0.0)) { # changed from any(rp <= 0.0) - 19.09.2024
      # if there is at least one coordinate with root being bracketed by the subinterval
      i.done <- i
      break
    }
    
    rf.left <- rf.right
    
  }
  
  if (i.done > 0) {
    # if a root has been bracketed
    rs <-
      which(rp < 0.0, TRUE) # check which coordinates of the vector root function this correspond to
    candidate <- 2.0
    which.root <- -1L
    for (rt in rs) {
      obj <- function(rr) {
        return(root.func(t + h * rr, z + rr * ff1 + rr ^ 2 * ff2 + rr ^ 3 * ff3, parms)[rt]) # select the part of the root function to the coordinates from above
      }
      interval <-
        r.grid[(i.done - 1):i.done] # define the intervals where roots have been bracketed
      if (i.done == 2 &&
          rt == left.modified) {
        interval[1] <-
          last.root.offset # if a root is found in the first interval, offset the time in order to not get the root at initial time point
      }
      
      ret <- uniroot(obj, interval = interval, tol = tol_uniroot)
      #print(ret)
      if (candidate > ret$root) {
        # in order to get out the smallest root among the roots found
        candidate <- ret$root
        which.root <- rt
      }
      
    }
    return(list(
      which.root = which.root,
      root.time = candidate,
      root.type = 0
    ))
  } else {
    return(list(
      which.root = -1L,
      root.time = 2.0,
      root.type = 0
    ))
  }
}


deSolverLinRootFinder <- function(z,
                                  z.new,
                                  k1,
                                  k4,
                                  t,
                                  h,
                                  parms,
                                  lin.root.func, # this should return A + B*y
                                  precision_real_root = 1.0e-13,
                                  last.root = NULL,
                                  last.root.offset = 1.0e-4) {
  
  linRootA <- lin.root.func(t, numeric(length(z)), parms)
  
  # sets up interpolation of A + B*y
  ff3 <- lin.root.func(t, h * (k1 + k4) - 2.0 * (z.new - z), parms) - linRootA #linRootB %*% (h*(k1+k4) - 2.0*(z.new-z))
  ff2 <- lin.root.func(t, 3.0 * (z.new - z) - h * (k4 + 2.0 * k1), parms) - linRootA #linRootB %*% (3.0*(z.new-z) - h*(k4+2.0*k1))
  ff1 <- lin.root.func(t, (k1 * h), parms) - linRootA #linRootB %*% (k1*h)
  ff0 <- lin.root.func(t, z, parms) #linRootA + linRootB %*% z
  #print(cbind(ff0,ff1,ff2,ff3))
  candidate <- 2.0
  which.root <- -1
  for (i in 1:length(ff0)) {
    out <- polyroot(c(ff0[i], ff1[i], ff2[i], ff3[i]))
    #print((out))
    real.roots <-
      Re(out[abs(Im(out)) < precision_real_root]) # check if there exists a real root
    #print(real.roots)
    lower.bound <- 0.0
    if (last.root$root.type == 1 &&
        last.root$which.root == i)
      lower.bound <-
      last.root.offset # if last root is also of linear type and for the same coordinate --> do this in order to not get the same root at initial time again
    if (length(real.roots) > 0) {
      allowed.roots <-
        real.roots[real.roots >= lower.bound & real.roots <= 1.0]
      if (length(allowed.roots) > 0) {
        if (min(allowed.roots) < candidate) {
          candidate <- min(allowed.roots)
          which.root <- i
        }
      }
    }
  }
  #print(candidate)
  #if(candidate<1.0) stop("")
  return(list(
    which.root = which.root,
    root.time = candidate,
    root.type = 1
  ))
}



#lin root occurs when either element of linRootA + linRootB %*% y is equal to zero
#then, lin.root.event.func gets called.

deSolverRoot <- function(y,
                         times,
                         func,
                         parms = NULL,
                         rtol = 1.0e-4,
                         atol = 1.0e-4,
                         root.func = NULL,
                         event.func = NULL,
                         lin.root.func = NULL,
                         lin.root.event.func = NULL,
                         last.root.offset.lin.root.finder = 1.0e-10,
                         last.root.offset.non.lin.root.finder = 1.0e-10,
                         precision_real_root_lin_root_finder = 1.0e-13,
                         num_subdiv_non_lin_root_finder = 8L,
                         tol_uniroot = 1e-8,
                         h.max = 0.3) {
  
  D <- length(y)
  z <- y
  Tmax <- times[length(times)]
  h <- h.max
  
  samples <- matrix(0.0, length(times), D + 1)
  t <- times[1]
  samples[1, 1] <- t
  samples[1, 2:(D + 1)] <- y
  samples.count <- 2L
  
  s.times <- c(times, Tmax + 1000.0) # used for while loop below
  
  k1 <- func(t, z, parms)
  
  PI.alpha <- 0.7 / 3.0
  PI.beta <- 0.4 / 3.0
  
  # if (!is.null(root.func)) {
  #r1 <- root.func(t,z,parms)
  
  # }
  
  # Comment - 14.08.24: Commented out everything related to storing event information, seems like this lead to memory issues when having situations with large dimension and frequent events
  
  # events.block.size <- 1000L
  # event.samples <- matrix(0.0, events.block.size, 3 + D)
  # event.count <- 0L
  
  last.root.obj <- list(
    which.root = -1L,
    root.time = 2.0,
    root.type = -1L
  )
  
  curr.root.obj <- last.root.obj
  
  n.evals <- 1L
  old.err <- 1.0
  
  while (t < Tmax) {
    not.done <- TRUE
    ntrials <- 0
    h <- min(h, Tmax - t)
    while (not.done && ntrials < 20L) {
      # Runge Kutta 3(2)
      ntrials <- ntrials + 1L
      ev1 <- z + 0.5 * h * k1
      k2 <- func(t + 0.5 * h, ev1, parms)
      
      ev2 <- z + 0.75 * h * k2
      k3 <- func(t + 0.75 * h, ev2, parms)
      
      z.new <-
        z + h * ((2.0 / 9.0) * k1 + (1.0 / 3.0) * k2 + (4.0 / 9.0) *
                   k3)
      k4 <- func(t + h, z.new, parms)
      z.low <-
        z + h * ((7.0 / 24.0) * k1 + 0.25 * k2 + (1.0 / 3.0) * k3 + (1.0 / 8.0) *
                   k4)
      err <-
        max(abs(z.new - z.low) / (atol + rtol * pmax(abs(z), abs(z.new))))
      n.evals <- n.evals + 3L
      if (is.finite(err)) {
        if (err < 1.0) {
          not.done <- FALSE
          break
        } else {
          h <- h * max(0.1, 0.9 * err ^ (-1.0 / 3.0))
        }
      } else {
        h <- 0.1 * h
      }
    }
    
    if (not.done) {
      print(k1)
      print(z)
      print(t)
      stop("problem with making progress!")
    }
    #print(paste0("h = ",h))
    
    curr.root.obj <-
      list(
        which.root = -1L,
        root.time = 2.0,
        root.type = -1L
      )
    
    if (!is.null(root.func)) {
      rt.out <-
        deSolverNonLinRootFinder(
          z = z,
          z.new = z.new,
          k1 = k1,
          k4 = k4,
          parms = parms,
          t = t,
          h = h,
          root.func = root.func,
          last.root = last.root.obj,
          last.root.offset = last.root.offset.non.lin.root.finder,
          init.subdiv = num_subdiv_non_lin_root_finder,
          tol_uniroot = tol_uniroot
        )
      
      curr.root.obj <- rt.out
      
    }
    
    if (!is.null(lin.root.func)) {
      
      rt.out <-
        deSolverLinRootFinder(
          z = z,
          z.new = z.new,
          k1 = k1,
          k4 = k4,
          parms = parms,
          t = t,
          h = h,
          lin.root.func = lin.root.func, 
          last.root = last.root.obj,
          last.root.offset = last.root.offset.lin.root.finder,
          precision_real_root = precision_real_root_lin_root_finder
        )
      
      if (rt.out$root.time < curr.root.obj$root.time) {
        curr.root.obj <- rt.out
      }
      
    }
    
    last <- min(1.0, curr.root.obj$root.time)
    #print(t)
    #print(curr.root.obj)
    #print(z.new)
    #if(curr.root.obj$which.root<0 && z.new[5] != model$region.id(z.new[3])){
    #  stop("")
    #}
    
    # sets up interpolation polynomial
    ff3 <- h * (k1 + k4) - 2.0 * (z.new - z)
    ff2 <- 3.0 * (z.new - z) - h * (k4 + 2.0 * k1)
    ff1 <- k1 * h
    
    # collect samples until min(root.time,t+h)
    while (t + h * last  >= s.times[samples.count] - 1.0e-14) {
      r <- (s.times[samples.count] - t) / h
      samples[samples.count, 1] <- s.times[samples.count]
      samples[samples.count, 2:(D + 1)] <-
        z + r * ff1 + r ^ 2 * ff2 + r ^ 3 * ff3
      samples.count <- samples.count + 1L
    }
    
    # do event if required
    if (curr.root.obj$which.root > 0) {
      # interpolated state at event
      z.at.event <- z + last * ff1 + last ^ 2 * ff2 + last ^ 3 * ff3
      #print("z.at.event")
      #print(z.at.event)
      # store info
      # event.count <- event.count + 1L
      
      # if (event.count > nrow(event.samples)) {
      #   event.samples <-
      #     rbind(event.samples, matrix(0.0, events.block.size, 3 + D))
      # }
      # 
      # event.samples[event.count,] <-
      #   c(t + h * last,
      #     curr.root.obj$which.root,
      #     curr.root.obj$root.type,
      #     z.at.event)
      
      
      #r4.before.event <- root.func(t+h*last,z.at.event,parms)
      if (curr.root.obj$root.type == 0) {
        z.new <- event.func(t + h * last, z.at.event, parms)
      } else {
        z.new <-
          lin.root.event.func(t + h * last, z.at.event, parms, lin.root.func)
      }
      
      k4 <- func(t + h * last, z.new, parms)
      n.evals <- n.evals + 1L
      
    }
    
    # prepare for next step
    k1 <- k4
    #if(!is.null(root.func)) r1 <- r4
    t <- t + h * last
    z <- z.new
    last.root.obj <- curr.root.obj
    #print("z last")
    #print(z)
    
    
    # adapt h using PI controller
    if (err > 1.0e-6) {
      h.fac <- 0.95 * (err ^ (-PI.alpha)) * (old.err ^ PI.beta)
      h <- h * min(2.0, max(0.2, h.fac))
    } else {
      h <- h.max # new safeguard to handle linear dynamics
    }
    
    h <- min(h, h.max)
    old.err <- err
    
    
  }
  
  colnames(samples) <- c("time", paste0("z", 1:D))
  # if (!is.null(root.func)) {
  #   colnames(event.samples) <-
  #     c("time", "rootDim", "rootType", paste0("z", 1:D))
  #   return(list(
  #     samples = samples,
  #     event.samples = event.samples[1:event.count,],
  #     n.evals = n.evals
  #   ))
  # } else {
  #   return(list(samples = samples, n.evals = n.evals))
  # }
  return(list(samples = samples, n.evals = n.evals))
}
