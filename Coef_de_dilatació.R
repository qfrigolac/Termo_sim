#Termodinàmica i Mecànica Estadística - Treball de Simulació
#Pau Sánchez Pascual i Joaquim Frigola Casals

############################ Diferents n ############################
alp_25 <-na.exclude(scan("Dades/300k_25.txt",sep = ","))
length(alp_25)
hist(alp_25,breaks=10,main="",xlab="Coeficient de dilatació (n=25)",ylab="Freqüència")
mean(alp_25)
t.test(alp_25,mu=1/300)
var(alp_25)


alp_50=na.exclude(scan("Dades/300k_50.txt",sep = ","))
length(alp_50)
hist(alp_50,breaks=10,main="",xlab="Coeficient adiabàtic",ylab="Freqüència")
mean(alp_50)
t.test(alp_50,mu=1/300)
var(alp_50)


alp_100=na.exclude(scan("Dades/300k_100.txt",sep = ","))
length(alp_100)
hist(alp_100,breaks=10,main="",xlab="Coeficient de dilatació (n=100)",ylab="Freqüència")
mean(alp_100)
t.test(alp_100,mu=1/300)
var(alp_100)

#################### Diferents temperatures ########################
alp_100k=na.exclude(scan("Dades/100k.txt",sep = ","))
length(alp_100k)
hist(alp_100k,breaks=10,main="",xlab="Coeficient de dilatació",ylab="Freqüència")
mean(alp_100k)
t.test(alp_100k,mu=1/100)
var(alp_100k)


alp_200k=na.exclude(scan("Dades/200k.txt",sep = ","))
length(alp_200k)
hist(alp_200k,breaks=10,main="",xlab="Coeficient de dilatació",ylab="Freqüència")
mean(alp_200k)
t.test(alp_200k,mu=1/200)
var(alp_200k)


alp_400k=na.exclude(scan("Dades/400k.txt",sep = ","))
length(alp_400k)
hist(alp_400k,breaks=10,main="",xlab="Coeficient de dilatació",ylab="Freqüència")
mean(alp_400k)
t.test(alp_400k,mu=1/400)
var(alp_400k)


alp_500k=na.exclude(scan("Dades/500k.txt",sep = ","))
length(alp_500k)
hist(alp_500k,breaks=10,main="",xlab="Coeficient de dilatació",ylab="Freqüència")
mean(alp_500k)
t.test(alp_500k,mu=1/500)
var(alp_500k)




punts=c(0.01390627,0.008378939,0.00542364,0.004275957,0.003550829)
Temp=c(100,200,300,400,500)
teo = function(x){1/x}

#El següents gràfics el fem amb l'origin
plot(teo,from=50,to=550,type='l',y=0.02,xlab="T (K)", ylab="Coef. adiabatic")
points(Temp,punts)


p=c(2.3897215783519408e-18,
    2.396296202290923e-18,
    2.4162346789302282e-18,
    2.4187215140633072e-18,
    1.6521649789569902e-18,
    1.668678225448954e-18,
    1.654679337630708e-18,
    1.657353800877854e-18,
    1.6653491688766922e-18,
    1.4443212293095254e-18,
    1.4605178130686356e-18,
    1.4656425158267226e-18,
    1.4622204592417303e-18,
    1.4639843349287476e-18,
    1.3434936837586252e-18,
    1.336435994506296e-18,
    1.332649069840184e-18,
    1.3302403495919165e-18,
    1.3316975269588793e-18,
    1.2957718204997624e-18,
    1.2976032587162596e-18,
    1.296944640705904e-18,
    1.2930239876995575e-18,
    1.2941921690085968e-18,
    1.269017059539235e-18,
    1.2708081058778388e-18,
    1.2673976453041416e-18,
    1.2693324440256035e-18,
    1.2703032559604826e-18,
    1.3040075387377173e-18,
    1.2994944338818565e-18,
    1.3045514992992336e-18,
    1.3003393469456358e-18,
    1.2971998399204385e-18)
inv_v=c(1.0,
        1.0,
        1.0,
        1.0,
        0.5,
        0.5,
        0.5,
        0.5,
        0.5,
        0.3333333333333333,
        0.3333333333333333,
        0.3333333333333333,
        0.3333333333333333,
        0.3333333333333333,
        0.25,
        0.25,
        0.25,
        0.25,
        0.25,
        0.2,
        0.2,
        0.2,
        0.2,
        0.2,
        0.16666666666666666,
        0.16666666666666666,
        0.16666666666666666,
        0.16666666666666666,
        0.16666666666666666,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285,
        0.14285714285714285)


reg=lm(p ~ inv_v)
summary(reg)
lin=function(x){1.332e-18*x+1.038e-18}
plot(lin,from=0.1, to=1,type='l',xlab="1/V (m-3)",ylab="p (Pa)")
points(inv_v,p)
