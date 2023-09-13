### File name should not contain any special character

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