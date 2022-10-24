implicit none 
dimension x(0:100001),u_n(0:100001),u_n1(0:100001)
dimension aa(0:100001),bb(0:100001),cc(0:100001) 
dimension dd(0:100001),ee(0:100001),ff(0:100001)
dimension uu(0:100001),ex(0:100001),errmax(0:100001) 
dimension k1(0:100001),k1x(0:100001),k1xx(0:100001) 
dimension k1xxx(0:100001),k2(0:100001),k2x(0:100001) 
dimension k2xx(0:100001),k2xxx(0:100001),f1(0:100001) 
dimension f1x(0:100001),f1xx(0:100001),f1xxx(0:100001) 
dimension f2(0:100001), xhalf(0:100001),f1t(0:100001)
dimension f2x(0:100001),f2xx(0:100001),f2xxx(0:100001) 
dimension k1half(0:100001),k1xhalf(0:100001),f2t(0:100001) 
dimension k1xxhalf(0:100001),k2xxhalf(0:100001) 
dimension k2half(0:100001),k2xhalf(0:100001)
dimension cf1(0:100001),cf2(0:100001),cf3(0:100001) 
dimension cf4(0:100001),cf5(0:100001),cf6(0:100001)
double precision h,dt,f1,ff1,fff1,D1,D2,DD1,DD2,AA1,AA2,A1,A2,cf5 
double precision k1,k2,rou,p1,p2,q1,q2,pp1,pp2,qq1,qq2,y1,y2,y3 
double precision x_o,x_1,u_n,u_n1,aa,bb,cc,dd,ee,ff,g5,gg5,ggg5 
double precision uu,ex,errmax,err,pi,x,beta_1,beta_2,cf1,cf2,cf3 
double precision cf4,f1xx,f1xxx,f1x,f2,f2x,f2xx,f2xxx,f1t,f2t 
double precision B1,B2,Z2,Z3,xhalf,s1,ss1,s2,ss2,s3,ss3 
double precision BB1,BB2,z ,Z1 ,i1,i2,i3 
double precision k1x,k2x,k1xx,k2xx,k1xxx,k2xxx,k1xhalf,k1half 
double precision k2half,k2xhalf,k1xxhalf,k2xxhalf,cf6 
double precision c1,c2,c3,c4,c5,c6,c7,c8,c9,c10 
double precision cc1,cc2,cc3,cc4,cc5,cc6,cc7,cc8,cc9,cc10 
integer kk,Nt,jj,mm,j,m,NN : 
 
NN=80 
x_o=0D0 
x_1=2D0 
dt=0.0001 
Nt=1D0/dt 
h=(x_1-x_o)/NN 
m=(NN)/2 
pi=3.1415926535897932384626D0 
rou = dt/(h*h) 

cc Grid Points 
do j=0,NN 
x(j)=j*h 
xhalf(j)=(j+0.5D0)*h 
enddo

cc Initial Condition 
do j=0,m-1 
u_n(j)=sin(2D0*pi*x(j)) 
enddo 
u_n(m)=sin(2D0*pi*x(m)) 
do j=m+1,NN 
u_n(j)=sin(4D0*pi*x(j)) 
enddo

cc Time Iteration 
kk=1 62: 
cc Defining of k(x) 
1 do j=0,NN 
k1(j)=2D0*(x(j)*x(j)+x(j)+1D0) 
k1half(j)=2D0*(xhalf(j)*xhalf(j)+xhalf(j)+1D0)
k1x(j)=2D0*(2D0*x(j)+1D0) 
k1xhalf(j)=2D0*(2D0*xhalf(j)+1D0) 
k1xx(j)=4D0 
k1xxhalf(j)=4D0 
k1xxx(j)=0D0 
k2(j)=(x(j)*x(j)+x(j)+1D0) 
k2half(j)=(xhalf(j)*xhalf(j)+xhalf(j)+1D0) 
k2x(j)=(2D0*x(j)+1D0) 
k2xhalf(j)=(2D0*xhalf(j)+1D0) 
k2xx(j)=2D0 
k2xxhalf(j)=2D0 
k2xxx(j)=0D0 
enddo

cc Defining f(x) and its required derivatives 
do j=0,NN 
f1(j)=((8D0*pi*pi*(x(j)*x(j)+x(j)+1D0)-1D0)*sin(2D0*pi*x(j)).
& 4D0*pi*(2D0*x(j)+1D0)*cos(2D0*pi*x(j)))*exp(-(kk-0.5D0)*dt) 
f1t(j)=-((8D0*pi*pi*(x(j)*x(j)+x(j)+1D0)-1D0)*sin(2D0*pi*x(j)).
& 4D0*pi*(2D0*x(j)+1D0)*cos(2D0*pi*x(j)))*exp(-(kk-0.5D0)*dt) 
f1x(j)=((16D0*pi*pi*(2D0*x(j)+1D0))*sin(2D0*pi*x(j))+ 
& (16D0*pi*pi*pi*(x(j)*x(j)+x(j)+1D0)-10D0*pi)*cos(2D0*pi*x(j)))* 
& exp(-(kk-0.5D0)*dt) 
f1xx(j)=((48D0*pi*pi*pi*(2D0*x(j)+1D0))*cos(2D0*pi*x(j)).
& 4D0*pi*pi*((8D0*pi*pi*(x(j)*x(j)+x(j)+1D0)-13D0)*
& sin(2D0*pi*x(j))))*exp(-(kk-0.5D0)*dt)
f2(j)=((16D0*pi*pi*(x(j)*x(j)+x(j)+1D0)-1D0)*sin(4D0*pi*x(j)).
& 4D0*pi*(2D0*x(j)+1)*cos(4D0*pi*x(j)))*exp(-(kk-0.5D0)*dt)
f2t(j)=-((16D0*pi*pi*(x(j)*x(j)+x(j)+1D0)-1D0)*sin(4D0*pi*x(j)).
& 4D0*pi*(2D0*x(j)+1D0)*cos(4D0*pi*x(j)))*exp(-(kk-0.5D0)*dt) 
f2x(j)=((32D0*pi*pi*(2D0*x(j)+1D0))*sin(4D0*pi*x(j))+ 
& (64D0*pi*pi*pi*(x(j)*x(j)+x(j)+1D0)-12D0*pi)*cos(4D0*pi*x(j)))* 
& exp(-(kk-0.5D0)*dt) 
f2xx(j)=(192D0*pi*pi*pi*(2D0*x(j)+1D0)*cos(4D0*pi*x(j)).
& 16D0*pi*pi*(16D0*pi*pi*(x(j)*x(j)+x(j)+1D0)-7D0)*
& sin(4D0*pi*x(j)))*exp(-(kk-0.5D0)*dt)


cc Functions for the RHS in the scheme :
cf1(j)=((-0.5D0+h/12D0*k1x(0)/k1(0)+h*h/48D0*(3D0*k1(0)* 
& k1xx(0)-4D0*k1x(0)*k1x(0))/(k1(0)*k1(0)))*f1(0)+ 
& (-h/6D0+h*h/24D0*k1x(0)/k1(0))*f1x(0)-h*h/24D0*f1xx(0)) 
cf2(j)=(-1D0+h*h/24D0*(k1(j)*k1xx(j)-k1x(j)*k1x(j))/ 
& (k1(j)*k1(j)))*f1(j)+h*h/12D0*k1x(j)/k1(j)*f1x(j).
& h*h/12D0*f1xx(j) 
cf3(j)=(-0.5D0-h/12D0*k1x(m)/k1(m)+h*h/48D0*(3D0*k1(m)*k1xx(m).
& 4D0*k1x(m)*k1x(m))/(k1(m)*k1(m)))*f1(m) + 
& (h/6D0+h*h/24D0*k1x(m)/k1(m))*f1x(m)-h*h/24D0*f1xx(m)
cf4(j)=((-0.5D0+h/12D0*k2x(m)/k2(m)+h*h/48D0*(3D0*k2(m)* 
& k2xx(m)-4D0*k2x(m)*k2x(m))/(k2(m)*k2(m)))*f2(m)+ 
& (-h/6D0+h*h/24D0*k2x(m)/k2(m))*f2x(m)-h*h/24D0*f2xx(m))
cf5(j)=(-1D0+h*h/24D0*(k2(j)*k2xx(j)-k2x(j)*k2x(j))/ 
& (k2(j)*k2(j)))*f2(j)+h*h/12D0*k2x(j)/k2(j)*f2x(j).
& h*h/12D0*f2xx(j) 
cf6(j)=((-0.5D0-h/12D0*k2x(NN)/k2(NN)+h*h/48D0*(3D0*k2(NN)* 
& k2xx(NN)-4D0*k2x(NN)*k2x(NN))/(k2(NN)*k2(NN)))*f2(NN)+ 
& (h/6D0+h*h/24D0*k2x(NN)/k2(NN))*f2x(NN)-h*h/24D0*f2xx(NN)) 147: 148: 
enddo 
c1=5D0/12D0-(h/12D0)*k1x(0)/k1(0)+((h*h)/48D0)* 
& (4D0*k1x(0)*k1x(0)-3D0*k1(0)*k1xx(0))/(k1(0)*k1(0)) 
c2=1D0/12D0 
c3=1D0+(h/k1(0))*((h/24D0)*((2D0*k1x(0)*k1x(0))/ 
& (k1(0))-k1xx(0))+h*h/48D0*(5D0*k1(0)*k1x(0)*k1xx(0).
& 4D0*k1x(0)*k1x(0)*k1x(0))/(k1(0)*k1(0))) 
p1=-h/12D0+h*h/24D0*k1x(0)/k1(0)

cc Triadiagonal System 
aa(0)=-0D0 168: bb(0)=c1+rou*k1half(0)/2D0 
cc(0)=-(c2-rou*k1half(0)/2D0) 
dd(0)=(c1-rou*k1half(0)/2D0)*u_n(0)+ 
& (c2+rou*k1half(0)/2D0)*u_n(1) 
& -dt/h*c3*2D0*pi*exp(-(kk-0.5D0)*dt)*k1(0) 
& -dt*cf1(0)+p1*(-2D0*pi*exp(-(kk-0.5D0)*dt))
do j=1,m-1 
c4=1D0/12D0
c5=5D0/6D0 
c6=1D0/12D0 
B1= c4+(h/24D0)*(k1x(j-1)/k1(j-1)) 
B2= c6-(h/24D0)*(k1x(j+1)/k1(j+1)) 
D1=2D0*k1xhalf(j)*k1xhalf(j)/k1half(j)-k1xxhalf(j) 
D2=2D0*k1xhalf(j-1)*k1xhalf(j-1)/k1half(j-1)-k1xxhalf(j-1) 
aa(j)=-(B1-rou*(k1half(j-1)/2D0)+dt/48D0*D2) 
bb(j)=(c5+rou*(k1half(j-1)/2D0)+rou*(k1half(j)/2D0).
& dt/48D0*(D1+D2)) 
cc(j)=-(B2-rou*(k1half(j)/2D0)+dt/48D0*D2) 
dd(j)=(B1+rou*(k1half(j-1)/2D0)-dt/48D0*D2)*u_n(j-1)+ 
& (c5-rou*(k1half(j-1)/2D0)-rou*(k1half(j)/2D0)+dt/48D0*(D1+D2))* 
& u_n(j)+(B2+rou*(k1half(j)/2D0)-dt/48D0*D1)*u_n(j+1)-dt*cf2(j) 
enddo 
c8=5D0/12D0+(h/12D0)*k1x(m)/k1(m)+((h*h)/48D0)* 
& ((4D0*k1x(m)*k1x(m)-3D0*k1(m)*k1xx(m))/(k1(m)*k1(m))) 
c7=1D0/6D0+h/24D0*k1x(m)/k1(m) 
c9=1D0+(h/k1(m))*((h/24D0)*((2D0*k1x(m)*k1x(m))/ 
& (k1(m))-k1xx(m))-h*h/48D0*(5D0*k1(m)*k1x(m)*k1xx(m).
& 4D0*k1x(m)*k1x(m)*k1x(m))/(k1(m)*k1(m))) 
cc7=5D0/12D0-(h/12D0)*k2x(m)/k2(m)+((h*h)/48D0)* 
& ((4D0*k2x(m)*k2x(m)-3D0*k2(m)*k2xx(m))/(k2(m)*k2(m))) 
cc8=1D0/6D0-h/24D0*k2x(m)/k2(m) 
cc9=1D0+(h/k2(m))*((h/24D0)*(2D0*k2x(m)*k2x(m)/ 
& (k2(m))-k2xx(m))+h*h/48D0*(5D0*k2(m)*k2x(m)*k2xx(m).
& 4D0*k2x(m)*k2x(m)*k2x(m))/(k2(m)*k2(m))) 
z=c9/cc9 
Z1=h*h/16D0+h*h*h/48D0*k1x(m)/k1(m) 
Z2=h*h/16D0-h*h*h/48D0*k2x(m)/k2(m) 
s3=1D0/(k1(m)/k2(m)+h/(4D0*k2(m))*(k1x(m)-k2x(m)*k1(m)/k2(m))) 
s1=2D0-s3*k1(m)/k2(m) 
s2=-(s1+s3) 
ss1=2D0/(k2(m)/k1(m)-1D0/h*(h*h/(2D0*k1(m))*(k2x(m)-k1x(m)* 
& k2(m)/k1(m))-h*k2(m)/k1(m))) 
ss3=2D0-ss1*k2(m)/k1(m) 227: ss2=-(ss1+ss3) 
i1= c7-s1*Z1/(h*h)-ss1*z*Z2/(h*h) 
i2= c8+z*cc7-s2*Z1/(h*h)-ss2*z*Z2/(h*h) 
i3= z*cc8-z*ss3*Z2/(h*h)-s3*Z1/(h*h)
aa(m)=-(i1-k1half(m-1)*rou/2D0) 
bb(m)=i2+k1half(m-1)*rou/2D0+z*k2half(m)*rou/2D0 
cc(m)=-(i3-z*k2half(m)*rou/2D0) 
dd(m)=(i1+k1half(m-1)*rou/2D0)*u_n(m-1)+ 
& (i2-k1half(m-1)*rou/2D0-z*k2half(m)*rou/2D0)*u_n(m) + 
& (i3+z*k2half(m)*rou/2D0)*u_n(m+1)-dt*(cf3(m)+z*cf4(m)) 
& +dt*(s3/(2D0*k2(m))*(f2t(m)-f1t(m))+ 
& z*ss1/(2D0*k1(m))*(f1t(m)-f2t(m))) 
do j=m+1,NN 244: 
cc4=1D0/12D0 
cc5=5D0/6D0 
cc6=1D0/12D0 
BB1= c4+(h/24D0)*(k2x(j-1)/k2(j-1)) 
BB2= c6-(h/24D0)*(k2x(j+1)/k2(j+1)) 
DD1=2D0*k2xhalf(j)*k2xhalf(j)/k2half(j)-k2xxhalf(j) 
DD2=2D0*k2xhalf(j-1)*k2xhalf(j-1)/k2half(j-1)-k2xxhalf(j-1)
aa(j)=-(BB1-rou*(k2half(j-1)/2D0)+dt/48D0*DD2) 
bb(j)=(c5+rou*(k2half(j-1)/2D0)+rou*(k2half(j)/2D0)-dt/48D0* 
& (DD1+DD2))
cc(j)=-(BB2-rou*(k2half(j)/2D0)+dt/48D0*DD1) 
dd(j)=(BB1+rou*(k2half(j-1)/2D0)-dt/48D0*DD2)*u_n(j-1)+ 
& (cc5-rou*(k2half(j-1)/2D0)-rou*(k2half(j)/2D0)+dt/48D0* 264: &(DD1+DD2))*u_n(j)+(BB2+rou*(k2half(j)/2D0)-dt/48D0*DD1)* 
& u_n(j+1)-dt*cf5(j) 266: 267: 268: 
enddo  
  
cc2=5D0/12D0+(h/12D0)*k2x(NN)/k2(NN)+((h*h)/48D0)*  
&  ((4D0*k2x(NN)*k2x(NN)-3D0*k2(NN)*k2xx(NN))/(k2(NN)*k2(NN)))  
cc1=1D0/12D0  
cc3=1D0+(h/k2(NN))*((h/24D0)*((2D0*k2x(NN)*k2x(NN))/  
&  (k2(NN))-k2xx(NN))-h*h/48D0*(5D0*k2(NN)*k2x(NN)*k2xx(NN). 
&  4D0*k2x(NN)*k2x(NN)*k2x(NN))/(k2(NN)*k2(NN)))  
  
  
p2=h/12D0+h*h/24D0*(k2x(NN)/k2(NN))  
    
aa(NN)=-(cc1-rou*k2half(NN-1)/2D0)  
bb(NN)=(cc2+rou*k2half(NN-1)/2D0)  
cc(NN)=-0D0  
dd(NN)=(cc1+rou*k2half(NN-1)/2D0)*u_n(NN-1)+  
&  (cc2-rou*k2half(NN-1)/2D0)*u_n(NN)+  
&  dt/h*cc3*4D0*pi*exp(-(kk-0.5D0)*dt)*k2(NN)-dt*cf6(NN) +  
&  p2*(-4D0)*pi*exp(-(kk-0.5D0)*dt)  

ee(0)=cc(0)/bb(0)  
ff(0)=dd(0)/bb(0)  
  
do jj=1,NN  
ff(jj)=(dd(jj)+aa(jj)*ff(jj-1))/(bb(jj)-aa(jj)*ee(jj-1))  
ee(jj)=cc(jj)/(bb(jj)-aa(jj)*ee(jj-1))  
enddo  
uu(NN+1)=0D0  
  
do jj=0,NN  
mm=NN-jj  
uu(mm)=ff(mm)+ee(mm)*uu(mm+1)  
enddo  
  
cc  Numerical solution at time n+1 level  
  
do j=0,NN  
u_n1(j)=uu(j)  
enddo  
  
cc  Exact solution and error  
do j=0,m  
ex(j)=exp(-dt*kk)*sin(2D0*pi*x(j))  
enddo  
do j=m+1,NN  
ex(j)=exp(-dt*kk)*sin(4D0*pi*x(j))  
enddo  
 
  
do j=1,NN-1  
err=abs(u_n1(j)-ex(j))  
enddo  
OPEN(unit=6,file='variable_version__5_80_4.dat')  
print *, err  
  
kk=kk+1    
if (kk.eq.NT+1) goto 2  
do j=0,NN  
u_n(j)=u_n1(j)  
enddo  
goto 1  
  
2 end  
