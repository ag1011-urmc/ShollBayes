model{
    for(i in 1:N){
        y[i] ~ dnegbin(lambda[i], phi)
        eta[i] = ifelse(x[i] <= g[Cell_ID[i]], a1[Cell_ID[i]] * (g[Cell_ID[i]] - x[i])^k1 + t[Cell_ID[i]], a2[Cell_ID[i]] * (x[i] - g[Cell_ID[i]])^k2 + t[Cell_ID[i]])
        log(lambda[i]) <- eta[i]
    }
    for(j in 1:N.Cell){
        a1[j] ~ dnorm(a1.image[Image_ID[j]], a1.tau.cell) T(, 0)
        a2[j] ~ dnorm(a2.image[Image_ID[j]], a2.tau.cell) T(, 0)
        g[j] ~ dnorm(g.image[Image_ID[j]], g.tau.cell) T(0, 98)
        t[j] ~ dnorm(t.image[Image_ID[j]], t.tau.cell) T(0,)
        k1[j] ~ dnorm(k1.image[Image_ID[j]], k1.tau.cell) T(1,)
        k2[j] ~ dnorm(k2.image[Image_ID[j]], k2.tau.cell) T(1,)
    }
    for(k in 1:N.Image){
        a1.image[k] ~ dnorm(a1.animal[Animal_ID[k]], a1.tau.image) T(, 0)
        a2.image[k] ~ dnorm(a2.animal[Animal_ID[k]], a2.tau.image) T(, 0)
        g.image[k] ~ dnorm(g.animal[Animal_ID[k]], g.tau.image) T(0, 98)
        t.image[k] ~ dnorm(t.animal[Animal_ID[k]], t.tau.image) T(0,)
        k1.image[k] ~ dnorm(k1.animal[Animal_ID[k]], k1.tau.image) T(1,)
        k2.image[k] ~ dnorm(k2.animal[Animal_ID[k]], k2.tau.image) T(1,)
    }
    for(l in 1:N.Animal){
        a1.animal[l] ~ dnorm(a1.pop, a1.tau.animal) T(, 0)
        a2.animal[l] ~ dnorm(a2.pop, a2.tau.animal) T(, 0)
        g.animal[l] ~ dnorm(g.pop, g.tau.animal) T(0, 98)
        t.animal[l] ~ dnorm(t.pop, t.tau.animal) T(0,)
        k1.animal[k] ~ dnorm(k1.pop, k1.tau.animal) T(1,)
        k2.animal[k] ~ dnorm(k2.pop, k2.tau.animal) T(1,)
    }

    #cell level
    a1.tau.cell <- pow(a1.sd.cell, -2)
    a2.tau.cell <- pow(a2.sd.cell, -2)
    g.tau.cell <- pow(g.sd.cell, -2)
    t.tau.cell <- pow(t.sd.cell, -2)
    k1.tau.cell <- pow(k1.sd.cell, -2)
    k2.tau.cell <- pow(k2.sd.cell, -2)

    a1.sd.cell ~ dt(0, 5000000, 4) T(0, ) # mean(location), precision(tau), df(nu)
    a2.sd.cell ~ dt(0, 5000000, 4) T(0, )
    g.sd.cell ~ dt(0, 0.5, 4) T(0, )
    t.sd.cell ~ dt(0, 50, 4) T(0, )
    k1.sd.cell ~ dt(0, 1, 4) T(0, ) ###########
    k2.sd.cell ~ dt(0, 1, 4) T(0, ) ###########

    #image level
    a1.tau.image <- pow(a1.sd.image, -2)
    a2.tau.image <- pow(a2.sd.image, -2)
    g.tau.image <- pow(g.sd.image, -2)
    t.tau.image <- pow(t.sd.image, -2)
    k1.tau.image <- pow(k1.sd.image, -2)
    k2.tau.image <- pow(k2.sd.image, -2)

    a1.sd.image ~ dt(0, 5000000, 4) T(0, ) # mean, precision, df
    a2.sd.image ~ dt(0, 5000000, 4) T(0, )
    g.sd.image ~ dt(0, 0.5, 4) T(0, )
    t.sd.image ~ dt(0, 50, 4) T(0, )
    k1.sd.image ~ dt(0, 1, 4) T(0, ) ###########
    k2.sd.image ~ dt(0, 1, 4) T(0, ) ###########

    #animal level
    a1.tau.animal <- pow(a1.sd.animal, -2)
    a2.tau.animal <- pow(a2.sd.animal, -2)
    g.tau.animal <- pow(g.sd.animal, -2)
    t.tau.animal <- pow(t.sd.animal, -2)
    k1.tau.animal <- pow(k1.sd.animal, -2)
    k2.tau.animal <- pow(k2.sd.animal, -2)

    a1.sd.animal ~ dt(0, 5000000, 4) T(0, ) # mean, precision, df
    a2.sd.animal ~ dt(0, 5000000, 4) T(0, )
    g.sd.animal ~ dt(0, 0.5, 4) T(0, )
    t.sd.animal ~ dt(0, 50, 4) T(0, )
    k1.sd.animal ~ dt(0, 1, 4) T(0, ) ###########
    k2.sd.animal ~ dt(0, 1, 4) T(0, ) ###########

    #pop level
    a1.pop ~ dnorm(0, 100000) T(, 0)
    a2.pop ~ dnorm(0, 100000) T(, 0)
    g.pop ~ dnorm(0, 0.01) T(0, 98)
    t.pop ~ dnorm(0, 0.25) T(0,)
    # log-gamma prior for \phi
}
