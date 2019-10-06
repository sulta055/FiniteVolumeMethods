!Passive Variable Transport Equation
!Pure Advcetion: p_t +v p_x=0
!Spatial domain: x in [xmin,xmax].
!Cell-interfaces values are approximated using
!the 5th order WENO scheme.
!TVD 3rd Order and 4th Order Runge-Kutta time-integrators are available for use.
!Reference:Jiang & Wu, J Comp. Phys. 150 (1999)

program scalarHyperbolicWENO_1D
implicit none

!*********************************************************************
!Global Variable Definitions
!*********************************************************************
integer,parameter::nx=50 !number of cells
integer,parameter::nt=1 !total time steps
real,parameter::eps=1.d-6
real,parameter::cfl=0.4 !cfl number
real,parameter::RK_option=1 !1:TVD RK3, 2: RK4

integer::is
real::xmin,xmax,dx,dt,x,t
integer::i,j

real::cr(0:2,0:2)
real::p(-3:nx+3)
real::p1(-3:nx+3),p2(-3:nx+3),p3(-3:nx+3)
real::flux0_plus(0:nx),flux1_plus(0:nx),flux2_plus(0:nx),flux3_plus(0:nx)  !f_i+1/2
real::flux0_minus(0:nx),flux1_minus(0:nx),flux2_minus(0:nx),flux3_minus(0:nx)
real::vavg_plus,vavg_minus
real::temp1,temp2,temp3,temp4
real::Lp0,Lp1,Lp2,Lp3
!**********************************************************************************

xmin=0.
xmax=1.
dx=(xmax-xmin)/real(nx)

!time step for linear advection with unit speed
!dt=cfl*dx
t=0.
open(unit=10,file='output_adv.txt')

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
  dt=cfl*dx/abs(v(0))
  do j=1,nx
    dt=min(dt,cfl*dx/abs(v(j)))
  end do

  print*,'Time Step,dt=',i,dt
  !Compute cell interface left-right states: v_i+1/2^(+),v_i+1/2^(-) 
 
  if(RK_option==1)then
  !Update cell-average value using 3rd order TVD 
  !Runge-Kutta time integration

  call compute_flux(flux0_plus,flux0_minus,p)
  do j=1,nx
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
   
    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
          -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    p1(j)=p(j)+dt*Lp0
  end do  

  call bound(p1)

  print*,'j	p_j	 p_j-1/2^+	pj+1/2^-'
  do j=1,nx
   print,j,p(j),
  end do
  
  !Step2:
  call compute_flux(flux1_plus,flux1_minus,p1)
  do j=1,nx
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
       -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    Lp1=-(v(j)/dx)*( (temp1*flux1_plus(j)+temp2*flux1_minus(j))&
       -(temp3*flux1_plus(j-1)+temp4*flux1_minus(j-1)) )

    p2(j)=p1(j)+(1./4.)*dt*(-3.*Lp0+Lp1)
  end do

  call bound(p2)

  !Step3:
  call compute_flux(flux2_plus,flux2_minus,p2)
  do j=1,nx 
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
 
    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
       -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    Lp1=-(v(j)/dx)*( (temp1*flux1_plus(j)+temp2*flux1_minus(j))&
       -(temp3*flux1_plus(j-1)+temp4*flux1_minus(j-1)) )

    Lp2=-(v(j)/dx)*( (temp1*flux2_plus(j)+temp2*flux2_minus(j))&
       -(temp3*flux2_plus(j-1)+temp4*flux2_minus(j-1)) )

    p(j)=p2(j)+(1./12.)*dt*(-Lp0-Lp1+8.*Lp2)
  end do

  else if(RK_option==2)then
  !Update cell-average value using 4th order 
  !Runge-Kutta time integration

  !Step 1:

  call compute_flux(flux0_plus,flux0_minus,p)
  do j=1,nx
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
   
    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
          -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    p1(j)=p(j)+(1./2.)*dt*Lp0
  end do  

  call bound(p1)

  !Step2:
  call compute_flux(flux1_plus,flux1_minus,p1)
  do j=1,nx
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))

    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
       -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    Lp1=-(v(j)/dx)*( (temp1*flux1_plus(j)+temp2*flux1_minus(j))&
       -(temp3*flux1_plus(j-1)+temp4*flux1_minus(j-1)) )

    p2(j)=p1(j)+(1./2.)*dt*(-Lp0+Lp1)
  end do

  call bound(p2)

  !Step3:
  call compute_flux(flux2_plus,flux2_minus,p2)
  do j=1,nx 
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
 
    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
       -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    Lp1=-(v(j)/dx)*( (temp1*flux1_plus(j)+temp2*flux1_minus(j))&
       -(temp3*flux1_plus(j-1)+temp4*flux1_minus(j-1)) )

    Lp2=-(v(j)/dx)*( (temp1*flux2_plus(j)+temp2*flux2_minus(j))&
       -(temp3*flux2_plus(j-1)+temp4*flux2_minus(j-1)) )

    p3(j)=p2(j)+(1./2.)*dt*(-Lp1+2.*Lp2)
  end do

  call bound(p3)

  !Step4:
 call compute_flux(flux3_plus,flux3_minus,p3)
  do j=1,nx 
    vavg_plus=0.5*(v(j)+v(j+1))
    vavg_minus=0.5*(v(j-1)+v(j))
    temp1=0.5*(1+sign(1.0,vavg_plus))
    temp2=0.5*(1-sign(1.0,vavg_plus))
    temp3=0.5*(1+sign(1.0,vavg_minus))
    temp4=0.5*(1-sign(1.0,vavg_minus))
 
    Lp0=-(v(j)/dx)*( (temp1*flux0_plus(j)+temp2*flux0_minus(j))&
       -(temp3*flux0_plus(j-1)+temp4*flux0_minus(j-1)) )

    Lp1=-(v(j)/dx)*( (temp1*flux1_plus(j)+temp2*flux1_minus(j))&
       -(temp3*flux1_plus(j-1)+temp4*flux1_minus(j-1)) )

    Lp2=-(v(j)/dx)*( (temp1*flux2_plus(j)+temp2*flux2_minus(j))&
       -(temp3*flux2_plus(j-1)+temp4*flux2_minus(j-1)) )
    
    Lp3=-(v(j)/dx)*( (temp1*flux3_plus(j)+temp2*flux3_minus(j))&
       -(temp3*flux3_plus(j-1)+temp4*flux3_minus(j-1)) )
    
    p(j)=p3(j)+(1./6.)*dt*(Lp0+2.*Lp1-4.*Lp2+Lp3)
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
integer::j
real::x

!Step function
do j=1,nx
  x=xmin+(j-1)*dx
!  if(x>=0.35 .and. x<=0.65)then
!   p(j)=1.
!  else  
!   p(j)=-1.
!  end if
  if(x<=0.5)then
    p(j)=1.0
  else
    p(j)=0.0
  end if

end do

call bound(p)

end subroutine initialize

subroutine compute_flux(flux_plus,flux_minus,p)
!Local Variables
real::flux_plus(0:nx),flux_minus(0:nx),p(-3:nx+3)
real::a1,b1,c1,d1
real::a2,b2,c2,d2
integer::i

do i=0,nx

  a1=f(p(i-1),i-1)-f(p(i-2),i-2)
  b1=f(p(i),i)-f(p(i-1),i-1)
  c1=f(p(i+1),i+1)-f(p(i),i)
  d1=f(p(i+2),i+2)-f(p(i+1),i+1)

  a2=f(p(i+3),i+3)-f(p(i+2),i+2)
  b2=f(p(i+2),i+2)-f(p(i+1),i+1)
  c2=f(p(i+1),i+1)-f(p(i),i)
  d2=f(p(i),i)-f(p(i-1),i-1)

  flux_plus(i)=(1./12.)*(-f(p(i-1),i-1)+7.*f(p(i),i)+7.*f(p(i+1),i+1)-f(p(i+2),i+2))&
          -phi(a1,b1,c1,d1)
  flux_minus(i)=(1./12.)*(-f(p(i-1),i-1)+7.*f(p(i),i)+7.*f(p(i+1),i+1)-f(p(i+2),i+2))&
          +phi(a2,b2,c2,d2)
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

function f(p,j) result(fx)

real,intent(in)::p
integer,intent(in)::j
real::fx

!advective flux
fx=p

end function f

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

function v(j) result(vx)
!Inputs
integer,intent(in)::j
!Outputs
real::vx
real::x

!Uniform advection to the left
!vx=1.

!Non-uniform advection advection
vx=p(j)


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
