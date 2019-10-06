module constants_mod
implicit none

!-----------------------------------------------------------------------------------------
integer,parameter::nt=4000!000
integer,parameter::nx=400
integer,parameter::ny=400
real*8,parameter::xmin=0.0,xmax=1.0,ymin=0.0,ymax=1.0,cour=0.8

integer,parameter::debug=0 !0:off 1:on
integer,parameter::discontinuityDir=1!1:x, 2:y, 3:oblique, 4:partial grid circle 
integer,parameter::boundaryType=1 !1:outflow 2:periodic in x, outflow in y, 3:reflecting 4:periodic, 5:peridoc in x reflecting in y 6:reflecting on left x, outflow  ontherwise
integer,parameter::initOption=0
integer,parameter::fluxType=3!1:Roe 2:HLL 3:HLLI
integer,parameter::HLLOption=0

integer,parameter::HLLIentropymode=1
integer,parameter::HLLIslowmode=1
integer,parameter::HLLIalfvenmode=1
integer,parameter::HLLIfastmode=0

integer,parameter::gravity=0 !0:off 1:on , uniform external gravitational field

integer,parameter::perturbationType=1!1:Sinusoidal 2:Random
integer,parameter::uniformAdvection_test=0
integer,parameter::KH_test=0 !0:off 1:on
integer,parameter::RT_test=0 !0:off 1: vertical 2:radial
integer,parameter::jet_test=0 !0:off 1:on
integer,parameter::cloudShock=1 !0:off 1:on
integer,parameter::randforce_test=0 !0:off 1:on
real,parameter::eps=0.00001

real*8,parameter::tol=1.d-20
real*8,parameter::min_pres=1.d-5
real*8,parameter::min_dens=1.d-5
real*8,parameter::pi=3.14159265359

real*8::gx,gy !gravitaional field components
real*8::divU(-1:nx+2,-1:ny+2)
real*8::shockCell(-1:nx+2,-1:ny+2)

integer,parameter::logOutput=1

!-------------------------------------------------------------------------------------------
integer,parameter::tracerType=2 !1:Mass Flux 2: Velocity Field
integer::trdebug=0!0:off, 1:on
integer,parameter::init_option=3! 1:match fluid density distribution 2:uniformly distribute tracers in cells [minCell,maxCell] 3:Match fluid distrinution in cells [minCell,maxCell] 
integer,parameter::trboundaryType=2 !1:periodic, 2:outflow, 3:periodic in x & outflow in y, 4:reflecting on left x, outflow  ontherwise
integer,parameter::Ninit=200000 !fixed # of tracers
integer,parameter::trInject=0 !0:off 1:on
integer::N 
integer,parameter::minCellx=100,maxCellx=300,minCelly=170,maxCelly=230 
integer,parameter::interpolation_option=1 ;!nearest_cell , 2:(slope limited)linear interpolation
integer,parameter::rand_displace=0 !offsets initial tracer placement from cell center by small random amount
integer,parameter::density_option=1 !1:regular 2: CIC

real*8::tracer_mass

!-------------------------------------------------------------------------------------------
integer,parameter::tSkip=100
!-------------------------------------------------------------------------------------------


end module constants_mod
