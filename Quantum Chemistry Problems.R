##QUANTUM CHEMISTTRY PROBLEMS 
### Number 2
L<- 1
x<- seq(from = 0, to = L, by = 0.0025) #mimic continuum
psi_org <- sqrt(1/L)*(0.25 - (x/L - 0.5)^2 )  ## get original psi

plot(x,psi_org)

n<- 1
f_n_x <- sqrt(2/L)*sin(n*pi*x/L)

c_n <- (-2*sqrt(2)/(n^3 * pi^3)) * ((-1)^n - 1) 
f_1 <- f_n_x*c_n 



plot(x,f_1, type = 'l', ylab = 'f1(x)', xlab = 'x')
lines(x,psi_org,col="green")


###Truncate at 4
n<- 3
c_n <- (-2*sqrt(2)/(n^3 * pi^3)) * ((-1)^n - 1) 
f_n_x <- sqrt(2/L)*sin(n*pi*x/L)
f_3 <- f_n_x*c_n
sum_fs <- f_1 + f_3
lines(x, sum_fs, col = "red")

legend(0.7, 0.25, legend=c("f1", "Psi Org", "f1+f3"), col=c("black","green", "red"), lty=1)

## NUMBER 3: Numberical integration 

##############################################
# PART (a) ####################################
##############################################

# check to see if  bessel J works the way we think
integrand<- function (x) {cos(-sin(x))}
I<- integrate(f = integrand, lower  = 0, upper = pi)
# it does so now proceed
delta <- 0.0025
x<- seq(from = -10, to = 10, by = delta)
#x<- seq(from = 0, to = 500, by = delta) ##just checking
f_x<- sqrt(sqrt(1/pi)) *exp(-x*x/2) ## just checking something
z<- pi*x
g_x_nA <-  besselJ(pi*abs(x), 0)*exp(-abs(x)/2) ##just checking the function
g<- function(x) {(besselJ(pi*abs(x), nu = 0)*exp(-abs(x)/2))} ## incase alternative approach will
f <- function(x){sqrt(sqrt(1/pi)) *exp(-x*x/2)}

## want to normalise g in <g|g>
g_norm<- function(x) {(besselJ(pi*abs(x), nu = 0)*exp(-abs(x)/2))^2}
f_norm <- function(x){(sqrt(sqrt(1/pi)) *exp(-x*x/2))^2}

#intergral using package to check approach
g_Area <- integrate(g_norm, -10, 10)
sqrt(g_Area[["value"]])

##Write trapeziod rule/solve with trapezoid rule 
sum_x <- matrix(nrow = 10, ncol = 1)
for (i in 1:10) {
  delta <- 0.0025
  x<- seq(from = -i , to = i, by= delta)
  g_x_nAp <- (besselJ(pi*abs(x), 0)*exp(-abs(x)/2) )
  g_x_nA_norm <- g_x_nAp%*%g_x_nAp
  sum_x[i,] <- (delta)*(g_x_nA_norm)
}

##check what g(0)  is 
plot(sum_x, type = 'l')

#converges rather quickly 
sqrt(sum_x[9,]) - sqrt(g_Area[["value"]]) #difference is acceptable 
A<- sqrt(sum_x[9,])
A<- 1/A ## to get constant
A
g_x_A <- A*g_x_nA

 #g_x_A <- g_Area[["value"]]*g_x_nA
 plot(x, g_x_A, type ='l')
 
 ##############################################
 ###### PART B##############################################
 ##############################################
 
 ##Normalise f(x)
 x<- seq(from = -10 , to = 10, by= delta)
 sum_fx <- matrix(nrow = 10, ncol = 1)
 for (i in 1:10) {
   delta <- 0.0025
   x<- seq(from = -i , to = i, by= delta)
   f_x<- sqrt(sqrt(1/pi)) *exp(-x*x/2)
   f_x_nA_norm <- f_x%*%f_x
   sum_fx[i,] <- (delta)*(f_x_nA_norm)
 }
 plot(sum_fx, type ='l')
 f_x<- sqrt(sqrt(1/pi)) *exp(-x*x/2)
 
 f_norm <- function(x){(sqrt(sqrt(1/pi)) *exp(-x*x/2))^2}
 f_area <- integrate(f_norm, lower = -10, upper = 10)
 f_area[["value"]]
 
 x<- seq(from = -10 , to = 10, by= delta) 
plot(x, g_x_A, ylim=c(-0.5, 1.25), type = 'l', ylab = 'g(x)', xlab = 'x')
lines(x, f_x, col = "red")

## f*g
fg <- f_x*g_x_A
lines(x,fg, col = "green")
legend(4, 1.2, legend=c("g(x)", "f(x)", "g(x)*f(x)"), col=c("black","red", "green"), lty=1)
       
### Want to check why g(0) = 1
besselJ(pi*abs(0), 0)*exp(-abs(0)/2)
besselJ(pi*abs(0), 0)
sqrt(sqrt(1/pi))
########################################################
##### PART C ############################
####################################################################################
sum_fg_x <- matrix(nrow = 20, ncol = 1)
for (i in 1:20) {
  delta <- 0.0025
  x<- seq(from = -i , to = i, by= delta)
  f_x_1<- sqrt(sqrt(1/pi)) *exp(-x*x/2)
  g_x_2 <- A*(besselJ(pi*abs(x), 0)*exp(-abs(x)/2) )
  fg2<- g_x_2*f_x_1
  sum_fg_x[i-1,] <- (delta)*sum(fg2)
}

fg_norm <- sum_fg_x[19,]
plot(sum_fg_x, type = 'l')

g_norm<- function(x) {A*(besselJ(pi*abs(x), nu = 0)*exp(-abs(x)/2))^2}
f_norm <- function(x){(sqrt(sqrt(1/pi)) *exp(-x*x/2))^2}
over_lap_func <- function(x){A*((besselJ(pi*abs(x), nu = 0)*exp(-abs(x)/2))) * ((sqrt(sqrt(1/pi)) *exp(-x*x/2)))}

over_area <- integrate(over_lap_func, -20,20)
over_area[["value"]]
over_area[["value"]] - sum_fg_x[19,]
sum_fg_x[19,]

################PART D#############################
##Will call on variables we already know because will rearrange the <phi|phi> 
##Let us respecify for clarification 
A_A<- (sum_x[9,])
fg_norm <- sum_fg_x[19,]
##get sum of f
x<- seq(from = -10 , to = 10, by= delta)
sum_fx <- matrix(nrow = 10, ncol = 1)
for (i in 1:10) {
  delta <- 0.0025
  x<- seq(from = -i , to = i, by= delta)
  f_x<- sqrt(sqrt(1/pi)) *exp(-x*x/2)
  f_x_nA_norm <- f_x%*%f_x
  sum_fx[i,] <- (delta)*(f_x_nA_norm)
}
plot(sum_fx, type ='l')

f_dot_area  <-sum_fx[9,] ## <f|f>

imag<- complex(real = 0, imaginary = pi/2)
conj_imag <- complex(real = 0, imaginary = -pi/2)

phase_shift <- imag+conj_imag ## is zero so can ignore

phi_area <- 1 + A_A
phi_area<- sqrt(phi_area)
phi_area<- 1/phi_area  

phi_approx <- f_dot_area + A_A
phi_approx<- sqrt(phi_approx)
phi_approx<- 1/phi_approx  
phi_area - phi_approx


#########################################################
phi_x <- f_x+ exp(Im(imag))*g_x_nA
phi_x_mag <- t(Conj(phi_x))%*%phi_x
(phi_x_mag[1:100])
plot(phi_x_mag,ylim=c(-0.000001, 0.000001))
sum_phi_x_mag<- matrix(nrow = 10, ncol = 1)
for (i in 1:10) {
  delta <- 0.0025
  x<- seq(from = -i , to = i, by= delta)
  f_x_1<- sqrt(sqrt(1/pi)) *exp(-x*x/2)
  g_x_2 <- (besselJ(pi*abs(x), 0)*exp(-abs(x)/2) )
  phi_x_loop <- f_x_1+ exp(imag)*g_x_2
  phi_x_mag_loop <- (delta)*Conj(phi_x)%*%phi_x
  sum_phi_x_mag[i,] <- (phi_x_mag_loop)
}
Re(sum_phi_x_mag[9,])
sum_phi_x_mag[9,]

phi_funct <- function(x){(exp(imag)* (besselJ(pi*abs(x), nu = 0)*exp(-abs(x)/2))) * ((sqrt(sqrt(1/pi)) *exp(-x*x/2)))}

integer()

I <- integrate(f = Re(phi_funct), -10,10) + integrate(f = Im(phi_funct), -10,10)

plot(Re(sum_phi_x_mag), type = 'l')



rep(c(1,2,3), each = 2)

##### Number 6 (b)
t <- 0
delt_E <- 3*pi*pi/2
delt_E<-1/delt_E
x<- seq(from  = 0, to = 1, by = 0.005)
psi_2 <- sin(2*x*pi) 
psi_1 <- sin(x*pi)
over_all_psi <- psi_2 + psi_1
norm_psi_t1 <- over_all_psi*over_all_psi
plot(x,norm_psi_t1, xlab = "x", ylab = 'psi(x)*psi(x)')

## t at pi/2
t <- pi/2
num <-complex(real = 0, imaginary = (-t*delt_E))
psi_2 <- sin(2*x*pi)*exp(num)
psi_1 <- sin(x*pi)
over_all_psi <- psi_2 + psi_1
norm_psi_t2 <- Conj(over_all_psi)*over_all_psi
lines(x,norm_psi_t2,col="green")
psi_2

### t = pi
t <- pi
num <-complex(real = 0, imaginary = (-t))
psi_2 <- sin(2*x*pi)*exp(num)
psi_1 <- sin(x*pi)
over_all_psi <- psi_2 + psi_1
norm_psi_t3 <- Conj(over_all_psi)*over_all_psi
lines(x,norm_psi_t3,col="red")


