wi=5
# color blind safe palette from https://jfly.uni-koeln.de/color/
set style line 1 pt 7 lw wi lc "black"  
set style line 2 pt 3 lw wi lc "#56b4e9"
set style line 3 pt 5 lw wi lc "#d55e00"
set style line 4 lw wi lc "#0072b2"
set style line 5 lw wi lc "#f0e442"
set style line 6 dt 3 lw wi lc "black"  
set style line 7 dt 3 lw wi lc "#56b4e9"
set style line 8 dt 3 lw wi lc "#d55e00"
set style line 9 dt 3 lw wi lc "#0072b2"
set style line 10 dt 3 lw wi lc "#f0e442"
set style line 101  lw wi lc "#00c000"

set style line 22 pt 3 lw wi ps 2 lc "#56b4e9"

#for filling
set style line 11 lw 1 lc "black"  
set style line 12 lw 1 lc "#56b4e9"
set style line 13 lw 1 lc "#d55e00"
set style line 14 lw 1 lc "#0072b2"
set style line 15 lw 1 lc "#f0e442"
set style line 16 pt 5 lw wi lc "red"

#for colorbar
#set palette defined (0 "white",\
#                     0.2 "red",\
#                     0.4 "orange-red",\
#                     0.6 "orange",\
#                     0.8 "yellow",\
#                     1.0 "light-green",\
#                     1.2 "green")
set palette defined (0 "white",\
                     0.2 "purple",\
                     0.4 "blue",\
                     0.6 "green",\
                     0.8 "yellow",\
                     1.0 "orange",\
                     1.2 "red")

set style fill  transparent solid 0.5 noborder 

set key reverse Left
set tics format "$%g$" front

d2r=pi/180

linear(x,a,b)=a+b*x
invlinear(z,a,b)=(z-a)/b

res(x,a,b,c)=sqrt(a*a/x+b*b+c*c/(x*x))

gauss(x,pos,sigma,h)=h*exp(-(x-pos)**2/(2*sigma**2))

chi2(x,ndf,A)=A * ndf / (2**(ndf/2) * gamma(ndf/2)) * (ndf*x)**(ndf/2 - 1) * exp(-ndf*x/2)

mp=.93827231
mC=mp*12
mup=2.79277
array a[16]=[0.15721e-1,0.38732e-1,0.36808e-1,0.14671e-1,-0.43277e-2, -0.97752e-2,-0.68908e-2,-0.27631e-2,-0.63538e-3,0.71809e-4,0.18441e-13,0.75066e-4,0.51069e-4,0.14308e-4,0.23170e-5,0.68465e-6]


Ep(theta,E)=E/(1+E/mp*(1-cos(theta*d2r)))
Q2(theta,E)=4*E*Ep(theta,E)*sin(theta/2*d2r)**2
tau(Q2)=Q2/(4*mp)
epsilon(Q2,theta)=1/(1+2*(1+tau(Q2)*tan(theta*d2r/2)**2))
#values from Povh and Rith particles and nuclei
#alpha(Q2)=sqrt(Q2)*2.5/197
R=8

rutherford(theta,E)=(1./137.*.197/(2*E*sin(theta/2*d2r)**2))**2/100
mott(theta,E,beta)=rutherford(theta,E)*(1-beta**2*sin(theta/2*d2r)**2)*Ep(theta,E)/E
#F(Q2)=3*1/alpha(Q2)**3*(sin(alpha(Q2)) -alpha(Q2)*cos(alpha(Q2)) )
denom=sum [i=1:16] (-1)**i * a[i]/(i**2*pi**2)
num(Q2) = sum [i=1:16] (-1)**i * a[i] /(i**2 * pi**2 -Q2/197**2*R**2)
F(Q2)= sin(sqrt(Q2)/197*R)/(sqrt(Q2)/197*R)*num(Q2)/denom

rosenbluthC(theta,E,beta)=mott(theta,E,beta)*F(Q2(theta,E))**2
# Arrington https://arxiv.org/pdf/0707.1861.pdf
a1 = 3.439
a2 = -1.602
a3 = 0.068
b1 = 15.055
b2 = 48.061
b3 = 99.304
b4 = 0.012
b5 = 8.650
#GE(Q2) = (1 + a1*Q2 + a2*Q2**2 + a3*Q2**3)/(1+ b1*Q2+b2*Q2**2 + b3*Q2**3 + b4*Q2**4 + b5*Q2**5);
#GM(Q2) = GE(Q2)*mup
GE(Q2) = (1+Q2/0.71)**-2
GM(Q2) = GE(Q2)*mup
rosenbluth(theta,E,beta)=mott(theta,E,beta)*(epsilon(Q2(theta,E),theta)*GE(Q2(theta,E)/1e6)**2+tau(Q2(theta,E))*GM(Q2(theta,E)/1e6)**2)/(epsilon(Q2(theta,E),theta)*(1+tau(Q2(theta,E))))

set xtics nomirror
set bars 3.0