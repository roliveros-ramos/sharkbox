
# Data
# Building simulator

model = smc(order|trip ~ b(block_sp) +
              mc(species, by=block_sp) + mc(size, by=species), data=dat0)

summary(model)

plot(model)

simulator = predict(model) # create a function: creates trajectories

# Simulating data
x0 = sim(N=300)
x1 = sim(N=1000)
plot(x0) # method for 1 to 3 classes
plot(x1) # method for 1 to 3 classes


n = 100
unloading_lengths = rdunif(n=n, min=300, max=1500)

unloading_data = sim(N=unloading_lengths) # a list with n simulations

# Sampling design experiments

