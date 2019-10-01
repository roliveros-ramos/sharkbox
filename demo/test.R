if(!require(sharkbox))
  devtools::install_github("roliveros-ramos/sharkbox")
library(sharkbox)
library(fields)

data_file = system.file("data/supersamples_demo.csv", package="sharkbox")
dat = read.csv(data_file)

# creating blocks
mod_block0 = blocks(order|trip ~ species, data=dat)
mod_block1 = blocks(order|trip ~ group, data=dat)
dat$block0 = fitted(mod_block0) # add fitted groups
dat$block1 = fitted(mod_block1) # add fitted groups
# re-name "block0" model levels
levels(dat$block0) = c("sharks", "dorado", "billfishes", "tunas", "swordfishes")
# groups are not equal to blocks based on groups
table(group=dat$group, block=dat$block1)
# species are not equal to blocks based on species
table(species=dat$species, block=dat$block0)
# comparison of two block classifications
table(model0=dat$block0, model1=dat$block1)

# a model
# start with simpler model
model = smc(order|trip ~ b(block0) + mc(species, by=block0) + m(group, by=species) + mc(size, by=group), data=dat)
model1 = smc(order|trip ~ b(block1) + mc(size, species, by=block0) + m(group, by=species) + mc(size, by=group), data=dat)

# another model
model0 = smc(order|trip ~ b(block0) + mc(species), data=dat)

# create simulators
sim = predict(model)
sim0 = predict(model0)

# make predictions
x = sim(300)
y = sim0(300)

# create simulated data
n = 50
unloading_lengths = rdunif(n=n, min=300, max=1500)

unloading_data = sim(N=unloading_lengths) # a list with n simulations

# Sampling
sampling(unloading_data[[1]], skip = 3, start = 0, fraction = 0.8)
