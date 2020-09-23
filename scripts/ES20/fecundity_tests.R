a_f1 <- pl_1_es$strategy$a_f1
a_f2 <- pl_1_es$strategy$a_f2
hmat <- pl_1_es$strategy$hmat


h <- seq(from = 0, to = hmat, length.out = 100)
fec_alloc <- a_f1/(1 + exp(a_f2 * (1 - h/hmat)))


plot(h, fec_alloc)
