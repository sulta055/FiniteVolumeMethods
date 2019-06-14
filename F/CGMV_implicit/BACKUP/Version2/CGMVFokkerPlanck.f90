!Coarse Grain Momentum Volume(CGMV) algorithm for Isotropic Cosmic Ray-Fokker Planck Equation

module FokkerPlanck_mod
use constants_mod
implicit none

integer,parameter::subcycleOption=1 !0: off, 1:on
real*8::fexact(0:np-1),dt2,dt3,dtf
real*8,parameter::cour2=0.5
integer::nt2

contains

subroutine CRinit()

integer::i,j
real*8::pres

open(unit=10,file='fp.txt')
open(unit=13,file='fx.txt')

open(unit=110,file='ng_p.txt')
open(unit=111,file='ng_x.txt')


!----------------------------------------------
!Initialization
!----------------------------------------------
!Momentum space bounds
pmin=1.D-3
pmax=1.D10

print*,'np,pmin,pmax=',np,pmin,pmax

nparticle=0.

!Generate momentum mesh points
do i=-2,np+3
  if(momentumMeshType==1)then
   px(i)=pmin+i*(pmax-pmin)/np
  else if(momentumMeshType==2)then
   px(i)=pmin*(pmax/pmin)**(real(i)/real(np))
  end if
  !print*,'i,px(i)=',i,px(i)
end do


!Momentum mesh spacing
do i=-1,np+2
 dp(i)=0.5*(px(i+1)-px(i-1))
 dp2(i)=px(i+1)-px(i)
 !print*,'i,p,dp=',i,px(i),dp(i)
end do


!Initialize distribution function
do i=0,np-1
 do j=-1,nx
  !if(px(i)>=0.1)then
  ! f(i,j)=px(i)**(-4.5)!*exp(-((xmin+j*dx)/(10.*dx))**2. )
  !else 
  ! f(i,j)=1.d-50
  !end if
  f(i,j)=px(i)**(-6.)

  !Gaussian in momentum space, uniform in physical space
  !f(i,j)=1000.*exp(-( (pmin+i*dp(i)-0.5*(pmax-pmin))/(10*dp(i)) )**2. )

  !Gaussian in physicsl space, uniform in momentum space
  !f(i,j)=1000.*exp(-( (xmin+j*dx-0.5*(xmax-xmin))/(50.*dx) )**2. )
  
  !Power law in moment space, uniform in phsycial space  
  !if(i<52)then
  !  f(i,j)=50.*px(i)
  !else 
  !  f(i,j)=50.*px(51)*(px(i)/px(51))**(-5.)
  !end if

  !f(i,j)=1.D-20


 end do
 f(-1,j)=f(0,j)
 f(np,j)=f(np-1,j)
end do

!----------------------------------------------------------------------
!Piece-wise power law interpolation: f_(p),j= f_i-1,j*(p/p_i-1)^q_i,j
!----------------------------------------------------------------------
do j=0,nx
 do i=1,np-1
   !Intrabin Spatial Index
   q(i,j)=-log(f(i-1,j)/f(i,j))/log(px(i-1)/px(i))
 end do
   q(0,j)=q(1,j)
   q(np,j)=q(np-1,j)
end do 


!do i=1,np-1
! print*,'i,q=',i,q(i,250)
!end do


!Normalize Distribution function such that CR pressure is equal to downstream fluid pressure
do j=0,nx-1
 do i=0,np-1
  Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp2(i)
 end do
end do

pres=(gam-1.)*(u1(nx,3)-0.5*u1(nx,2)*u1(nx,2)/u1(nx,1))

do i=0,np-1
 do j=0,nx-1
   f(i,j)=f(i,j)*(pres/Pc(j))
  end do
end do

!------------------------------------------------
!CR distribution function moments: n_i,j, g_i,j
!------------------------------------------------
do j=0,nx
 do i=1,np
   dw=px(i)/px(i-1)
   n(i,j)=f(i-1,j)*(px(i-1)**3.)/(3.-q(i,j))
   n(i,j)=n(i,j)*(dw**(3.-q(i,j))-1.)

   g(i,j)=f(i-1,j)*(px(i-1)**4.)/(4.-q(i,j))
   g(i,j)=g(i,j)*(dw**(4.-q(i,j))-1.)
 end do
 n(0,j)=n(1,j)
 g(0,j)=g(1,j)
end do 

if(mod(i,tSkip)==0)then

do i=1,np-1
 write(110,*) i,px(i),n(i,shockCell(1)-2),g(i,shockCell(1)-2),q(i,shockCell(1)-2)
end do

do j=0,nx-1
 write(111,*) j,xmin+dx*j,n(10,j),g(10,j),q(10,j)
end do

end if

!Compute CR Pressure
Pc=0.
do j=0,nx-1
  do i=0,np-1
    !Pc(j)=Pc(j)+(4.*cl/3.)*(px(i)**3.)*f(i,j)*dp(i)
    Pc(j)=Pc(j)+(4.*cl/3.)*g(i,j)
  end do
end do
Pc(-1)=Pc(0)
Pc(-2)=Pc(0)
Pc(nx)=Pc(nx-1)
Pc(nx+1)=Pc(nx-1)

!do j=0,nx-1
!   pres=(gam-1.)*(u1(j,3)-0.5*u1(j,2)*u1(j,2)/u1(j,1))
!  print*,'j,Pc,Pg=',j,Pc(j),pres
!end do

print*,'Done initializing CR distribution.'

end subroutine

subroutine CRevolve()
!Evolve Distribution Function
integer::i,j,k
real*8::Fn(-1:np),Fp(-1:np),qtemp
real*8::ntot,gtot,ftot,fluxsum1,fluxsum2

!dudx=-0.5
go to 1111
if(subcycleOption==1)then
!Momentum advection time step
!dt2=dt*2.
dt2=dt  
do j=0,nx-1
 if(abs(dudx(j))>1.d-20)then
  do k=1,np-1
   dt2=min(dt2,cour2*3.*dp(k)/(px(k)*abs(dudx(j))))
  end do
 !else
 ! dt2=min(dt2,dt)
 end if
end do

if(dt2<dt)then
 nt2=floor(dt/dt2)
 dtf=dt-nt2*dt2
 print*,'Momentum advection time step, dt2=',dt2
 print*,'# of momentum sweep sub cycles=',nt2+1
else
  nt2=1
end if
dt3=dt2

else if(subcycleOption==0)then
 nt2=0
 dt2=dt
 dt3=dt
 dtf=0.
end if
1111 continue

if(subcycleOption==1)then
!Momentum advection time step
dt2=dt*2.  
do j=0,nx-1
 if(abs(dudx(j))>1.d-20)then
  do k=0,np-1
   dt2=min(dt2,cour2*3.*dp(k)/(px(k)*abs(dudx(j))))
  end do
 else
  dt2=min(dt2,dt)
 end if
end do

if(dt2<dt)then
 nt2=floor(dt/dt2)
 dtf=dt-nt2*dt2
 print*,'Momentum advection time step, dt2=',dt2
 print*,'# of momentum sweep sub cycles=',nt2+1
else
  nt2=1
end if
dt3=dt2

else if(subcycleOption==0)then
 nt2=0
 dt2=dt
 dt3=dt
end if

!print*,'CHECKPOINT1'
!------------------------------------------
!Momentum Space Evolution
!------------------------------------------ 
do i=1,nt2+1
  do j=0,nx-1
   !---------------------------
   !Upwinded flux
   !---------------------------
   do k=1,np-2
    dw=px(k)/px(k-1)
    if(dudx(j)<=0._8)then
     Fn(k)=f(k-1,j)*dw**(-q(k,j))
     Fn(k)=Fn(k)*(1./3.)*dudx(j)*(px(k)**3.)
    else
     Fn(k)=f(k,j)
     Fn(k)=Fn(k)*(1./3.)*dudx(j)*(px(k)**3.)
    end if
    Fp(k)=px(k)*Fn(k)

    !Momentum Space Diffusion contribution
    if(q(k,j)>0._8)then
      Fn(k)=Fn(k)-q(k,j)*px(k)*D(k)*f(k-1,j)*dw**(-q(k,j))
    else
      Fn(k)=Fn(k)-q(k+1,j)*px(k)*D(k)*f(k,j)
    end if
   end do
   !No-flux momentum boundary conditions
   Fn(0)=0._8
   Fn(np-1)=0._8
   Fp(0)=0._8
   Fp(np-1)=0._8
 
   fluxsum1=0._8
   fluxsum2=0._8
   do k=1,np-1
    fluxsum1=fluxsum1+(Fn(k)-Fn(k-1))
    fluxsum2=fluxsum2+(Fp(k)-Fp(k-1))
   end do
   if(abs(fluxsum1)>1.d-13 .or. abs(fluxsum2)>1.d-13 )then
     !print*,'Flux check: sum_i=0..np-1 (Fn), sum_i=0..np-1 (Fp)=',fluxsum1,fluxsum2
   end if

   do k=1,np-1
     !print*,'k=',k
     n(k,j)=n(k,j)+dt2*(Fn(k)-Fn(k-1))!+dt2*Qi(k,j,1)
     !print*,'n=',n(k,j)
     g(k,j)=g(k,j)+dt2*(Fp(k)-Fp(k-1)) &
            +dt2*(-(1./3.)*dudx(j)+q(k,j)*D(k)/(px(k)*px(k)))*g(k,j)!+dt2*Qi(k,j,2)     
   end do

   !Ghost cell values
   n(0,j)=n(1,j)
   n(np,j)=n(np-1,j)
   g(0,j)=g(1,j)
   g(np,j)=g(np-1,j)

  end do
  if(i==nt2) dt2=dtf
end do
 
!------------------------------------------
!Physical Space Evolution
!------------------------------------------  
do i=1,nt2+1
 do j=1,np-1
 
  ftemp2=0.
 
  !Evolve zeroth moment: n_i,j
  do k=-1,nx
   ftemp2(k)=n(j,k)
  end do  

  call computeSpatialCoefficients(j,1)
  call tridiagSolveSpace(nx,ftemp2)

  ftemp2(-1)=ftemp2(0)  
  ftemp2(nx)=ftemp2(nx-1)

  do k=-1,nx
  n(j,k)=ftemp2(k)   
  end do

  !Evolve first moment: g_i,j
  do k=-1,nx
   ftemp2(k)=g(j,k)
  end do  

  call computeSpatialCoefficients(j,2)
  call tridiagSolveSpace(nx,ftemp2)

  ftemp2(-1)=ftemp2(0)  
  ftemp2(nx)=ftemp2(nx-1)

  do k=-1,nx
  g(j,k)=ftemp2(k)   
  end do

  !Periodic Boundary
  n(j,-1)=n(j,nx-1)
  g(j,-1)=g(j,nx-1)
  n(j,nx)=n(j,0)
  g(j,nx)=g(j,0)



 end do

 if(i==nt2) dt3=dtf
end do

!print*,'CHECKPOINT2'

!----------------------------------------------
!Reconstruct distribution function from moments
!----------------------------------------------
do j=0,nx-1   
 do k=1,np!-1
  dw=px(k)/px(k-1)

  !Update Spectral Index
  !print*,'j,k,q0,g/np=',j,k,q(k,j),g(k,j)/n(k,j)/px(k-1)
  call rootSolveNewton(q(k,j),g(k,j)/n(k,j)/px(k-1),dw,qtemp)  
  
  if(abs(qtemp-3.)<1.d-3 .or. abs(qtemp-4.)<1.d-3)then
   !print*,'Root solver returned bad value...adjusting.'
   qtemp=qtemp+0.01  
  end if
  q(k,j)=qtemp

  !Update Distribution function bin edge values
  f(k-1,j)=(3.-q(k,j))*n(k,j)/(px(k-1)**3.)
  f(k-1,j)=f(k-1,j)/(dw**(3.-q(k,j))-1.)
 end do
 !f(np-1,j)=f(np-2,j)*dw**(-q(k,j)) 
 f(-1,j)=f(0,j)
 f(np,j)=f(np-1,j)
end do

go to 125
ftot=0.
ntot=0._8
gtot=0._8

do j=0,nx-1
do k=1,np-1
  ftot=ftot+f(k,j)
  ntot=ntot+n(k,j)
  gtot=gtot+g(k,j)
end do
end do
print*,'ftot,ntot,gtot=',ftot,ntot,gtot
125 continue

do k=1,np-1
 write(110,*) k,px(k),n(k,shockCell(1)-2),g(k,shockCell(1)-2),q(k,shockCell(1)-2)
end do

do j=0,nx-1
 write(111,*) j,xmin+dx*j,n(5,j),g(5,j),q(5,j)
end do

!-----------------------------------------------
!Update CR pressure
!-----------------------------------------------
Pc=0.
do j=0,nx-1
  do k=0,np-1
    Pc(j)=Pc(j)+(4.*cl/3.)*g(k,j)
  end do
end do
Pc(-1)=Pc(0)
Pc(-2)=Pc(0)
Pc(nx)=Pc(nx-1)
Pc(nx+1)=Pc(nx-1)

end subroutine CRevolve


subroutine computeSpatialCoefficients(pk,option)

!Input variables
integer::pk,xk,option
integer::jj
real*8::Wplus(-1:nx),Wminus(-1:nx)
real*8::B2,C2,B1,C1,pi2,pi2L,w2
real*8::delta(-1:nx),Cp(-1:nx)

!-------------------------------------------------------------
!Weights: w_j, delta_j, W^+_j+1/2, W^-_j+1/2 
!-------------------------------------------------------------
do jj=-1,nx
 B2=0.5*(E(jj)+E(jj+1))
 C2=Kxx(jj,pk)

 if(abs(C2)<1.d-20)then
  if(B2<0._8)then
   delta(jj)=1.
  else
   delta(jj)=0._8
  end if
 else
  w2=dx*B2/C2
  if(abs(w2)<1.d-3)then
    delta(jj)=0.5-(1./12.)*w2   
  else
    delta(jj)=(1./w2)-(1./(exp(w2)-1.))
  end if
 end if

end do

do jj=-1,nx
  B2=0.5*(E(jj)+E(jj+1))
  Cp(jj)=Kxx(jj,pk)/dx
  Wminus(jj)=B2*delta(jj)
  Wplus(jj)=B2*(1.-delta(jj))
end do

!---------------------------------
!No-flux Boundary condition 
!---------------------------------
Wplus(0)=0.
Wminus(0)=0.
Cp(0)=0.
Wplus(-1)=0.
Wminus(-1)=0.
Cp(-1)=0.

!Right boundary open
!Wplus(nx-1)=0.
!Wminus(nx-1)=0.
!Cp(nx-1)=0.



!------------------------------------------------------------
!Co-efficients for tridiagonal system: ~A_j,~B_j,~C_j
!------------------------------------------------------------
do jj=0,nx-1
  Atx(jj)=dt3/dx
  Atx(jj)=Atx(jj)*(Wminus(jj-1)-Cp(jj-1))
  if(methodType==2)then
   Atx(jj)=Atx(jj)*0.5 
  end if
  
  Ctx(jj)=dt3/dx
  Ctx(jj)=-Ctx(jj)*(Wplus(jj)+Cp(jj))
  if(methodType==2)then
   Ctx(jj)=Ctx(jj)*0.5
  end if

   Btx(jj)= Wminus(jj)-Wplus(jj-1)-Cp(jj)-Cp(jj-1)
   Btx(jj)=Btx(jj)*dt3/dx
   Btx(jj)=-Btx(jj)+1.

  if(methodType==2)then 
   Btx(jj)=0.5*(Btx(jj)+1.)
  end if

  !print*,'jj,Atx,Btx,Ctx=',jj,Atx(jj),Btx(jj),Ctx(jj)
end do

end subroutine computeSpatialCoefficients

subroutine tridiagSolveSpace(n,ftemp) !from Numerical Recipes book

!Input variables
integer::n
real*8::ftemp(-1:n)

real*8::bet,gam(0:n-1),r(0:n-1)

integer::kk


do kk=0,n-1
 r(kk)=ftemp(kk)
 if(methodType==2)then
   r(kk)=r(kk)+ftemp(kk)-Atx(kk)*ftemp(kk-1)-Btx(kk)*ftemp(kk)-Ctx(kk)*ftemp(kk+1)
 end if
end do


bet=Btx(0)
ftemp(0)=r(0)/bet

!Upper triangular decomposition and forward substitution
do kk=1,n-1
  gam(kk)=Ctx(kk-1)/bet
  bet=Btx(kk)-Atx(kk)*gam(kk)
  ftemp(kk)=(r(kk)-Atx(kk)*ftemp(kk-1))/bet
end do

!Back substitution
do kk=n-2,0,-1
 ftemp(kk)=ftemp(kk)-gam(kk+1)*ftemp(kk+1)
end do

end subroutine tridiagSolveSpace

!------------------------------------------------------------------!
!------------------------------------------------------------------!
!Momentum Diffusion co-efficient
function D(pk) result(fx)

integer,intent(in)::pk
real*8::fx,Dpp

 Dpp=0._8
 fx=Dpp

end function D

!Spatial Advection co-efficient
function E(xk) result(ex)
integer,intent(in)::xk
real*8::ex

!ex=-1.
ex=-u1(xk,2)/u1(xk,1)

end function E

!Spatial Diffusion co-efficient
function Kxx(xk,pk) result(kx)

integer,intent(in)::xk,pk
real*8::kx

!kx=0._8
kx=0.5/u1(xk,1)
!kx=0.1*(px(pk)**0.51)
end function Kxx

!CR particle number Injection term
function Qi(pk,xk,option) result(sx)
 
integer,intent(in)::pk,xk,option
real*8::sx

real*8,parameter::alpha=2.
real*8::eps=0.001
real*8::mp=0.1

real*8::cs2,rho1,rho2,pinj,us,pres,wx
real*8::ut1,ut2
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Shock Speed
ut1=-2.*sqrt(gam)    !u1(shockCell(1),2)/u1(shockCell(1),1)
ut2=0.               !u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1)
us=(1./3.)*(4.*ut2-ut1)     !(1./3.)*2.*sqrt(gam)
!pre-shock density
rho1=1.!u1(shockCell(1),1)
!Post-shock sound speed
pres=(gam-1.)*(u1(shockCell(1)-1,3) &
-0.5*u1(shockCell(1)-1,2)*u1(shockCell(1)-1,2)/u1(shockCell(1)-1,1))
rho2=u1(shockCell(1)-1,1)
cs2=sqrt(gam*pres/rho2)

!Compute injection momentum
pinj=alpha*cs2
pin=np*log(pinj/pmin)/log(pmax/pmin)

!print*,'pinj,pin=',pinj

!Spatial weight function
wx=(dx*xk-dx*(shockCell(1)-1))/(4.*dx)
wx=exp(-wx**2.)/sqrt(3.141592*(4.*dx)**2.)
!if(xk==shockCell(1)-1)then
!  wx=1.
!else
!  wx=0.
!end if

!Injection source term
if(pinj<= px(pk) .and. pinj>= px(pk-1))then
 sx=wx*eps*rho1*us/mp 
else
 sx=0._8
end if

if(option==2)then
sx=sx*pinj
end if

end function Qi

function fluidS(xk) result(fx)
integer,intent(in)::xk
real*8::fx

real*8,parameter::alpha=2.
real*8::eps=0.005
real*8::mp=0.1

real*8::cs2,rho1,rho2,pinj,us,pres,wx
real*8::ut1,ut2
integer::pin

!------------------------------
!Flux fraction injection model
!------------------------------
!Shock Speed
ut1=-2.*sqrt(gam)    !u1(shockCell(1),2)/u1(shockCell(1),1)
ut2=0.               !u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1)
us=(1./3.)*(4.*ut2-ut1)     !(1./3.)*2.*sqrt(gam)
!pre-shock density
rho1=1.!u1(shockCell(1),1)
!Post-shock sound speed
!pres=55.78043
!pres=(gam-1.)*(u1(shockCell(1)-2,3) &
!-0.5*u1(shockCell(1)-2,2)*u1(shockCell(1)-2,2)/u1(shockCell(1)-2,1))
!rho2=u1(shockCell(1)-2,1)
pres=(gam-1.)*(u1(shockCell(1)-1,3) &
-0.5*u1(shockCell(1)-1,2)*u1(shockCell(1)-1,2)/u1(shockCell(1)-1,1))
rho2=u1(shockCell(1)-1,1)
cs2=sqrt(gam*pres/rho2)

!Compute injection momentum
!pinj=lambda*sqrt(gam*pres/ 3.976468)
pinj=alpha*cs2
pin=np*log(pinj/pmin)/log(pmax/pmin)

!Spatial weight function
wx=(dx*xk-dx*(shockCell(1)-1))/(4.*dx)
wx=exp(-wx**2.)/sqrt(3.141592*(4.*dx)**2.)
!if(xk==shockCell(1)-1)then
!  wx=1.
!else
!  wx=0.
!end if

fx=0.5*eps*wx*(pin**2.)*rho1*us

end function fluidS
!------------------------------------------------------------------!
!------------------------------------------------------------------!

subroutine rootSolveNewton(qi0,gnp,dw,qi)
real*8::qi0,gnp,dw,qi
integer::niter=50
real*8::q0
real*8::itertol=1.D-5
integer::i,iterCount=0

!Intial guess
q0=qi0

do i=1,niter
  qi=q0-psi(q0,gnp,dw)/psip(q0,gnp,dw)
  if(debug==1)then
   print*,'iteration #',i+iterCount*20
   print*,'q0,qi=',q0,qi
  end if 

  if(abs(qi-q0)<itertol)then
    if(debug==1)then
    print*,'Iterations converged.'
    print*,'qi0,qi, #iteration=',qi0,qi,i+iterCount*20
    end if 
    EXIT
  end if
  if(abs(qi-q0)>itertol .and. i==niter-1 .and. iterCount<10)then
    niter=niter+20
    iterCount=iterCount+1
  end if

  if(abs(qi-q0)>itertol .and. i==niter-1 .and. iterCount==10)then
     print*,'Iterations fail to converge...terminating suboutine.'
     qi=qi0
     EXIT
  end if

  q0=qi
end do


end subroutine rootSolveNewton

function psi(qi,gnp,dw) result(fx)
real*8,intent(in)::qi,gnp,dw
real*8::fx,t1,t2,t3,t4

t1=3.-qi
t2=(dw**(4.-qi))-1.
t3=-gnp*(4.-qi)
t4=(dw**(3.-qi))-1.

fx=t1*t2+t3*t4

end function psi

function psip(qi,gnp,dw) result(fx)
real*8,intent(in)::qi,gnp,dw
real*8::fx,t1,t2,t3,t4

t1=1.-(dw**(4.-qi))
t2=-(3.-qi)*log(dw)*(dw**(4.-qi))
t3=gnp*((dw**(3.-qi))-1.)
t4=gnp*log(dw)*(4.-qi)*(dw**(3.-qi))

fx=t1+t2+t3+t4

end function psip


end module FokkerPlanck_mod