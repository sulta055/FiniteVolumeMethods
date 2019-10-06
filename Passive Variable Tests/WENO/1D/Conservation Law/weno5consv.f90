!Passive Variable Transport Equation
!Conservation Law: p_t + f_x=0, f=v p, where v is the fluid veocity 
!Spatial domain: x in [xmin,xmax].
!Fluxes at cell-interfaces are computed using
!the 5th order WENO scheme, with (Lax-Friedrichs) flux splitting.
!TVD 3rd Order and 4th Order Runge-Kutta time-integrators are available for use.
!Reference:Jiang & Wu, J Comp. Phys. 150 (1999)

program scalarHyperbolicWENO_1D
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=100 !number of cells
integer,parameter::nt=200 !total time steps
real,parameter::eps=1.d-6
real,parameter::cfl=0.4 !cfl number
real,parameter::RK_option=2 !1:TVD RK3, 2: RK4

integer::is
real::xmin,xmax,dx,dt,x,t
integer::i,j

real::cr(0:2,0:2)
real::p(-3:nx+3)
real::p0(-3:nx+3),p1(-3:nx+3),p2(-3:nx+3)
real::flux0(0:nx),flux1(0:nx),flux2(0:nx),flux3(0:nx)  !f_i+1/2
!*********************************************************************

xmin=0.
xmax=1.
dx=(xmax-xmin)/real(nx)

!time step for linear advection with unit speed
!dt=cfl*dx
t=0.
open(unit=10,file='output_consv.txt')

!Set initial state
call initialize()
print*,'Done initialization.'

!Set polynomial interpolation coefficients
cr(0,0)=1./3.
cr(0,1)=-7./6.
cr(0,2)=11./6.
cr(1,0)=-1./6.
cr(1,1)=5./6.
cr(1,2)=1./3.
cr(2,0)=1./3.
cr(2,1)=5./6.
cr(2,2)=-1./6.
 
do i=1,nt

  !compute time-step size
  !dt=cfl*dx/abs(v(0))
  !do j=1,nx
  !  dt=min(dt,cfl*dx/abs(v(j)))
  !end do
  dt=cfl*dx/abs(p(0))
  do j=1,nx
    dt=min(dt,cfl*dx/abs(p(j)))
  end do


  print*,'Time Step,dt=',i,dt
  !Compute cell interface left-right states: v_i+1/2^(+),v_i+1/2^(-) 
 
  if(RK_option==1)then
  !Update cell-average value using 3rd order TVD 
  !Runge-Kutta time integration

  !Step 1:

  call compute_flux(flux0,p)
  do j=1,nx
    p0(j)=p(j)-(dt/dx)*(flux0(j)-flux0(j-1))
  end do  

  call bound(p0)
  p=p0
  
  !Step2:
  call compute_flux(flux1,p0)
  do j=1,nx
    p1(j)=p0(j)-(1./4.)*(dt/dx)*(-3.*(flux0(j)-flux0(j-1))+(flux1(j)-flux1(j-1)))
  end do

  call bound(p1)
  p=p1

  !Step3:
  call compute_flux(flux2,p1)
  do j=1,nx 
    p(j)=p1(j)-(1./12.)*(dt/dx)*(-(flux0(j)-flux0(j-1))-(flux1(j)-flux1(j-1))&
        +8.*(flux2(j)-flux2(j-1)))
  end do

  else if(RK_option==2)then
  !Update cell-average value using 4th order 
  !Runge-Kutta time integration

  !Step 1:
  call compute_flux(flux0,p)
  do j=1,nx
    p0(j)=p(j)-(1./2.)*(dt/dx)*(flux0(j)-flux0(j-1))
  end do  

  call bound(p0)
  p=p0

  !Step2:
  call compute_flux(flux1,p0)
  do j=1,nx
    p1(j)=p0(j)-(1./2.)*(dt/dx)*(-(flux0(j)-flux0(j-1))+(flux1(j)-flux1(j-1)))
  end do

  call bound(p1)
  p=p1

  !Step3:
  call compute_flux(flux2,p1)
  do j=1,nx 
    p2(j)=p1(j)-(1./2.)*(dt/dx)*(-(flux1(j)-flux1(j-1))+2.*(flux2(j)-flux2(j-1)))
  end do

  call bound(p2)
  p=p2

  !Step4:
  call compute_flux(flux3,p2)
  do j=1,nx 
    p(j)=p2(j)-(1./6.)*(dt/dx)*((flux0(j)-flux0(j-1))+2.*(flux1(j)-flux1(j-1))&
          -4.*(flux2(j)-flux2(j-1))+(flux3(j)-flux3(j-1)))
  end do


  end if

  call bound(p)

  t=t+dt

  do j=1,nx
    x=xmin+(j-1)*dx
    write(10,*) x,p(j),pexact(x) 
  end do

end do



print*,'Done.'

close(unit=10)


contains
!*********************************************************************
!Subroutine and Function Definitions
!*********************************************************************

subroutine initialize()
integer::i
real::x

!Step function
do i=1,nx
  x=xmin+(i-1)*dx
  if(x>=0.35 .and. x<=0.65)then
    p(i)=1.
  else  
   p(i)=-1.
  end if
end do

call bound(p)

p0=p

end subroutine initialize

subroutine compute_flux(flux,p)
!Local Variables
real::flux(0:nx),p(-3:nx+3)
real::a1,b1,c1,d1
real::a2,b2,c2,d2
integer::i

do i=0,nx

  a1=fplus(p(i-1),i-1)-fplus(p(i-2),i-2)
  b1=fplus(p(i),i)-fplus(p(i-1),i-1)
  c1=fplus(p(i+1),i+1)-fplus(p(i),i)
  d1=fplus(p(i+2),i+2)-fplus(p(i+1),i+1)

  a2=fminus(p(i+3),i+3)-fminus(p(i+2),i+2)
  b2=fminus(p(i+2),i+2)-fminus(p(i+1),i+1)
  c2=fminus(p(i+1),i+1)-fminus(p(i),i)
  d2=fminus(p(i),i)-fminus(p(i-1),i-1)

  flux(i)=(1./12.)*(-f(p(i-1),i-1)+7.*f(p(i),i)+7.*f(p(i+1),i+1)-f(p(i+2),i+2))&
          -phi(a1,b1,c1,d1)+phi(a2,b2,c2,d2)

end do



end subroutine compute_flux


function phi(a,b,c,d) result(phix)
!Inputs
real,intent(in)::a,b,c,d
!Outputs
real::phix
!Local Variables
real::w0,w2,alpha0,alpha1,alpha2,alphasum,IS0,IS1,IS2

IS0=13.*(a-b)**2+3.*(a-3.*b)**2
IS1=13.*(b-c)**2+3.*(b+c)**2
IS2=13.*(c-d)**2+3.*(3.*c-d)**2

alpha0=1./((eps+IS0)**2)
alpha1=6./((eps+IS1)**2)
alpha2=3./((eps+IS2)**2)

alphasum=alpha0+alpha1+alpha2

w0=alpha0/alphasum
w2=alpha2/alphasum

phix=(1./3.)*w0*(a-2.*b+c)+(1./6.)*(w2-0.5)*(b-2.*c+d)

end function phi

function f(px,i) result(fx)

real,intent(in)::px
integer,intent(in)::i
real::fx

!advective flux
!fx=v(i)*p
fx=0.5*px*px

end function f

function fplus(px,i) result(fx)

real,intent(in)::px
integer,intent(in)::i
real::fx,alpha

!alpha:=max(|f'(p)|) 
!alpha=max(abs(v(i)),abs(v(i+1)),abs(v(i-1))) 
alpha=max(abs(p(i)),abs(p(i+1)),abs(p(i-1))) 

fx=0.5*(f(px,i)+alpha*px)

end function fplus

function fminus(px,i) result(fx)

real,intent(in)::px
integer,intent(in)::i
real::fx,alpha

!alpha:=max(|f'(p)|) 
!alpha=max(abs(v(i)),abs(v(i+1)),abs(v(i-1)))  
alpha=max(abs(p(i)),abs(p(i+1)),abs(p(i-1))) 

fx=0.5*(f(px,i)-alpha*px)

end function fminus

function pexact(x) result(px)
!Inputs
real,intent(in)::x
!Outputs
real::px

integer::i

i=1+floor((x-xmin)/dx)

if(x>=0.4+v(i)*t .and. x<0.5+v(i)*t)then
  px=1.
else
  px=0.
end if

end function pexact

function v(i) result(vx)
!Inputs
integer,intent(in)::i
!Outputs
real::vx

real::x

!Uniform advection to the left
!vx=1.

!Non-uniform Advection (Burgers equation)
vx=1.

end function v

subroutine bound(p)
integer::j
real::p(-3:nx+3)

!Outflow boundary condition
do j=0,3
  p(-j)=p(1)
  p(nx+j)=p(nx)
end do

end subroutine bound

end program scalarHyperbolicWENO_1D
