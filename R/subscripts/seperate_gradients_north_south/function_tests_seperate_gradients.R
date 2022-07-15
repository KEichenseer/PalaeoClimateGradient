gradient1 <- function(lat,A,DKA,Q,M) A + DKA/((1+(exp(Q*(lat-M)))))
gradient2 <- function(lat,A,DKA,Q,M) A + DKA/(-(1+(exp(Q*(lat-M)))))

gradient_two_A <- function(lat,A1,A2,DKA,Q1,M1,Q2,M2) A1 + (DKA-A1)/((1+(exp(Q1*(lat-M1))))) + (DKA-A2)/(-(1+(exp(Q2*(lat-M2)))))
gradient_one_A <- function(lat,A,DKA,Q1,M1,Q2,M2) A + (DKA-A)/((1+(exp(Q1*(lat-M1))))) + (DKA-A)/(-(1+(exp(Q2*(lat-M2)))))

DKA <- 30
A1 <- 10
Q1 <- 0.2
M1 <- 45

A2 <- 15
Q2 <- 0.2
M2 <- 45
plot(lat,(DKA-A1)/((1+(exp(Q1*(lat-M1))))))
plot(lat,(DKA-A2)/(-(1+(exp(Q2*(lat-M2))))))

lat <- seq(-90,90,1)
grad1 <- gradient1(lat,0,30,.075,45)
#grad2 <- gradient2(lat,0,30,.1,-45)
grad2a <- gradient_two_A(lat,10,-10,30,.1,45,.1,-45)
grad1a <- gradient_one_A(lat,-10,30,.1,45,.1,-35)

plot(lat,grad1,type = "l", ylim = c(-10,60))
#points(lat,grad2,type = "l", col = "red")
points(lat,grad1a,type = "l", lty = 1, lwd = 2, col = rgb(.75,0,0,0.5))
points(lat,grad2a,type = "l", lty = 1, lwd = 2, col = rgb(0,0,0.75,0.5))

abline(v=0)
abline(h=29)
