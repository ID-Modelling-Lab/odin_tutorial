# hello
library(odin.dust)
sirodindust <- odin.dust::odin_dust({
  
  
  ## Definition of the time-step and output as "time"
  dt <- user()
  initial(time) <- 0
  update(time) <- (step + 1) * dt
  
  # print("{step}")
  ## Core equations for transitions between compartments:
  update(S) <- S - n_SI
  update(I) <- I + n_SI - n_IR
  update(R) <- R + n_IR
  
  ## Individual probabilities of transition:
  p_SI <- 1 - exp(-beta * I / N * dt) # S to I
  p_IR <- 1 - exp(-gamma * dt) # I to R
  
  ## Draws from binomial distributions for numbers changing between
  ## compartments:
  n_IR <- rbinom(I, p_IR)
  n_SI <- rbinom(S, p_SI)
  
  ## Total population size
  N <- S + I + R
  
  ## Initial states:
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- 0
  
  ## User defined parameters - default in parentheses:
  S_ini <- user(1000)
  I_ini <- user(10)
  beta <- user(0.2)
  gamma <- user(0.1)
},debug_enable=TRUE)

### Initializing the model
sir_model<-sirodindust$new(pars = list(dt = 1,
                                       S_ini = 1000,
                                       I_ini = 10,
                                       beta = 0.2,
                                       gamma = 0.1),
                           time = 1,
                           n_particles = 10L,
                           n_threads = 4L,
                           seed = 10L)

# sir_model$run(10) ### will genrate at time=10, the values of state variables

end_step=10
dt=1
steps <- seq(0, end_step/dt, by = 1)

### sir_odin_dust
n_times <- 200
x <- array(NA, dim = c(sir_model$info()$len, n_particles=10, n_times))

for (t in seq_len(n_times)) {
  x[ , , t] <- sir_model$run(t)
}
time <- x[1, 1, ]
x <- x[-1, , ]

par(mar = c(4.1, 5.1, 0.5, 0.5), las = 1)
cols <- c(S = "#8c8cd9", I = "#cc0044", R = "#999966")
matplot(time, t(x[1, , ]), type = "l",
        xlab = "Time", ylab = "Number of individuals",
        col = cols[["S"]], lty = 1, ylim = range(x))
matlines(time, t(x[2, , ]), col = cols[["I"]], lty = 1)
matlines(time, t(x[3, , ]), col = cols[["R"]], lty = 1)
legend


incidence <- read.table("sir_incidence_data.csv", header = TRUE, sep = ",")
# true_history <- readRDS("sir_true_history.rds")

sir_data <- mcstate::particle_filter_data(data = incidence,
                                          time = "day",
                                          initial_time = 0,
                                          rate = 1/dt)

