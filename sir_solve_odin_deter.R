### SOlving SIR model with Odin and odin.dust


## Method 1
##  Here we supply all the parameters  and initial conditions

setwd("~/codes/dengue_vaccination_singapore/try/scripts")


sir <- odin::odin({
  ## Derivatives
  deriv(S) <- mu*(S+I+R) -beta*S*I/(S+I+R) -mu*S
  deriv(I) <- beta*S*I/(S+I+R)  - gamma*I - mu*I
  deriv(R) <- gamma*I - mu*R
  
  ## Initial conditions
  initial(S) <- 10000
  initial(I) <- 10
  initial(R) <- 0
  
  ## parameters
  beta<- interpolate(tt, beta_t, "spline")
  gamma <- user()
  mu<-user()
  beta_t[]<-user()
  tt[]<-user()
  dim(beta_t)<- user()
  dim(tt)<- user()
  
})

t <- seq(0, 5000, by=1)
tt=t
beta_t=0.3+ 0*sin(tt)

sir_model <- sir$new(tt=tt, beta_t=beta_t, gamma=0.1, mu=0.01)

y <- sir_model$run(t)


## Plotting
sir_col <- c("#8c8cd9", "#cc0044", "#999966")
par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
matplot(y[, 1], y[, -1], xlab = "Time", ylab = "Number of individuals",
        type = "l", col = sir_col, lty = 1)
legend("topright", lwd = 1, col = sir_col, legend = c("S", "I", "R"), bty = "n")


matplot(rowSums(y[,-1]), type="l")


## Method 2
##Here we supply all the parameters/ some part of the parameters as user-defined 


sir <- odin::odin({
  ## Derivatives
  deriv(S) <- -beta*S*I/(S+I+R) 
  deriv(I) <- beta*S*I/(S+I+R)  - gamma*I
  deriv(R) <- gamma*I
  
  ## Initial conditions
  initial(S) <- S0
  initial(I) <- I0
  initial(R) <- R0
  
  ## parameters & Initial conditions
  beta  <- user()
  gamma <- 0.1
  S0 <- 10000
  I0 <- user()
  R0 <- 0
  
})


sir_model <- sir$new(beta=0.25, I0=10 )
t <- seq(0, 100, length.out = 50000)
y <- sir_model$run(t)

plot(t,y[,3])



#### Method 3

### Write the model in a different file like "sirdeter.R" and call odin

sir_model= odin::odin("sirdeter.R")

sir_model1 <- sir_model$new(beta=0.25, I0=10 )
t <- seq(0, 100, length.out = 50000)
y <- sir_model1$run(t)

plot(t,y[,3])




##### SIR using only 2D array in ODIN LINK: https://mrc-ide.github.io/odin/articles/odin.html


sir_2D_array <- odin::odin({
  deriv(S[]) <- -beta[i]*S[i]*I[i]/N 
  deriv(I[]) <- beta[i]*S[i]*I[i]/N - gamma[i]*I[i] 
  deriv(R[]) <- gamma[i]*I[i]
  
  
  initial(S[]) <- S0[i]
  initial(I[]) <- I0[i]
  initial(R[]) <- R0[i]
  
  
  S0[] <- user()
  I0[] <- user()
  R0[] <- user()
  beta[] <- user()
  gamma[] <- user()
  
  N <- sum(S[])+sum(I[])+sum(R[])
  # print("N: {N}")
  print("{S[1]} {S[2]} {S[3]}")
  
  n_dim <- user()
 
  
  dim(S) <- n_dim
  dim(I) <- n_dim
  dim(R) <- n_dim
  dim(S0) <- n_dim
  dim(I0) <- n_dim
  dim(R0) <- n_dim
  
  dim(beta) <- n_dim
  dim(gamma) <- n_dim
  
  # config(base) <- "lv4"
},  debug_enable = TRUE)



pars <- list(beta = c(0.5, 0.6, 0.3, 0.4),
             n_dim=4,
             gamma = c(0.2, 0.2, 0.1, 0.2),
             S0 = c(1000, 1100, 1300, 1200),
             I0 = c(1, 1, 1, 1),
             R0 = c(0,0,0,0)
             )

mod_sir_2D_array <- sir_2D_array$new(user = pars)

t <- seq(0, 3, length.out = 4)
out <- mod_sir_2D_array$run(t)


matplot(out[, 6], type='l')



######################################### 3D array  ###########################################

sir_3D_array <- odin::odin({
  
  foi[1:n_sero]<-beta_m[i]*I_m[i]
  
  S_diag[1:n_age,1:n_sero, 1:n_sero] <- I[i,j]
  
  # I_ij_sumk[1:n_age,1:n_sero] <- sum(I_ij[,,k])
  
  deriv(S[1]) <- lambda_hum*N - sum(foi[])*S[i]/N - death_hum*S[i] - age_rate[i]*S[i] 
  
  deriv(S[2:n_age]) <- - sum(foi[])*S[i]/N - death_hum*S[i] + age_rate[i-1]*S[i-1]- age_rate[i]*S[i] 
  
  deriv(I[1,1:n_sero]) <- foi[j]*S[i]/N - death_hum*I[i,j] - age_rate[i]*I[i,j] 
  deriv(I[2:n_age,1:n_sero]) <- foi[j]*S[i]/N - death_hum*I[i,j] +age_rate[i-1]*I[i-1,j] - age_rate[i]*I[i,j]
  
  deriv(I_ij[1:n_age,1:n_sero,1:n_sero]) <- I_ij[i,j,k] + S_diag[i,j,k]
  
  
  deriv(I_m[1:n_sero])<- sum(I_ij[i,,]) + sum(x[,i])
  
  # S_diag_zero[,,] <- if (j==k) 0 else S_diag_zero[i,j,k]
  # print("{S_diag[1,1,1]} {S_diag[1,1,2]}")
  
  # print("{I_ij[1,2,3]}")
  
  # S_diag[,] <- if (i==j) 0 else S_diag[i,j]
  
  x[1:n_age,1:n_sero] <- sum(I_ij[i,j,])
  
  print("{x[1,1]} {x[1,2]}")
  
  initial(S[]) <- S0[i]
  initial(I[,]) <- I0[i,j]
  initial(I_m[]) <- I_m0[i]
  initial(I_ij[,,]) <- I_ij0[i,j,k]
  
  
 
  
  S0[] <- user()
  I0[,] <- user()
  I_m0[] <- user()
  I_ij0[,,] <- user()
  beta_m[] <- user()
  death_hum <- user()
  age_rate[] <- user()
  lambda_hum <- user()
  
  N <- sum(S[])+sum(I[,])
  # print("x: {x[1]} {x[2]} {x[3]} {x[4]} {x[5]}")
  # print("{foi[1]} {foi[2]} {foi[3]}")
  
  n_age <- user()
  n_sero <- user()
  
  
  dim(S) <- n_age
  dim(I) <- c(n_age, n_sero)
  dim(I_m) <- n_sero
  dim(I_ij) <- c(n_age, n_sero, n_sero)
  dim(foi) <- n_sero
  dim(S0) <- n_age
  dim(I0) <- c(n_age, n_sero)
  dim(I_m0) <- n_sero
  dim(I_ij0) <- c(n_age, n_sero, n_sero)
  dim(S_diag) <- c(n_age, n_sero, n_sero)
  # dim(I_ij_sumk) <- c(n_age,n_sero)
  dim(beta_m) <- n_sero
  dim(age_rate) <- n_age
  dim(x) <- c(n_age, n_sero)
  
  # config(base) <- "lv4"
},  debug_enable = TRUE)

I_ij0=array(, dim=c(3,4,4))
I_ij0[1,,]= matrix(1:16, nrow=4, ncol=4, byrow=TRUE)
I_ij0[2,,]= matrix(17:32, nrow=4, ncol=4, byrow=TRUE)
I_ij0[3,,]= matrix(33:48, nrow=4, ncol=4, byrow=TRUE)


pars <- list(beta_m = c(0.5, 0.6, 0.3, 0.4),
             n_age=3,
             n_sero=4,
             death_hum=1/100,
             lambda_hum=1/10,
             age_rate=c(0.1,0.2,0.3),
             S0 = c(1000, 1100, 1300),
             I0 = matrix(1:12, nrow=3, ncol=4, byrow = TRUE),
             I_m0 = c(10,10,10,10),
             I_ij0 =I_ij0 # array(1:48, dim=c(4, 4, 3))
)

mod_sir_3D_array <- sir_3D_array$new(user = pars)

t <- seq(0, 3, length.out = 4)
out <- mod_sir_3D_array$run(t)


matplot(out[, 6], type='l')


















####################### Solving Generalized Lotka-Volterra Model

gen <- odin::odin({
  deriv(y[]) <- r[i] * y[i] * (1 - sum(ay[i, ]))
  initial(y[1:4]) <- y0[i]
  
  y0[] <- user()
  r[] <- user()
  a[, ] <- user()
  ay[, ] <- a[i, j] * y[j]
  
  dim(r) <- user()
 
  n_spp <- length(r)
  
  # print("dim(r): {dim(r)}")
  
  dim(y) <- n_spp
  dim(y0) <- n_spp
  dim(a) <- c(n_spp, n_spp)
  dim(ay) <- c(n_spp, n_spp)
  
  # config(base) <- "lv4"
},  debug_enable = TRUE)



pars <- list(r = c(1.00, 0.72, 1.53, 1.27),
             a = rbind(c(1.00, 1.09, 1.52, 0.00),
                       c(0.00, 1.00, 0.44, 1.36),
                       c(2.33, 0.00, 1.00, 0.47),
                       c(1.21, 0.51, 0.35, 1.00)),
             y0 = c(0.3013, 0.4586, 0.1307, 0.3557))



mod <- gen$new(user = pars)

t <- seq(0, 2000, length.out = 10001)
y <- mod$run(t)
pairs(y[, -1], panel = lines, col = "#00000055", lwd = 0.2)


















