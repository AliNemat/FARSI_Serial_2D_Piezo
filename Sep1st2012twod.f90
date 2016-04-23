
 program June5th2012twod
 use M_General_2D          
 use M_Platform_Constant_2D
 use M_Mesh_2D          

   
    implicit none

parameter iterate=2


dOUBLE PRECISION Tx(2:nx-1,1:ny-1),Ty(1:nx-1,2:ny-1)

dOUBLE PRECISION miu,miudrop,miuair,rodrop,roair
double precision dt,div,divmax,beta,maxp,pdif
double precision xbar,ybar,ubar,vbar,gx,omegaz


double precision period,hait,HH,divmax2,Tmass,gap,eps
double precision xbarold,ybarold,ubarold,vbarold,omegazold,torder,teta,sumphi,sumphi0,yfree
double precision ubargh,vbargh,omegazgh,anac,anacexp,tolen,acgx,acgy
double precision acgxold,acgyold,anacold,tetaold,leg,ks,sumugh,sumvgh
double precision tetaplot,xgage,yfreeg,xbar0,ybar0,xbarold2,ybarold2,sumu,sumv,sumuold,sumvold,Drag,Lift,Dragold,Liftold

double precision,dimension (1:3)                  ::   lzero,lnew,ox,oy
double precision,dimension (1:nx,0:ny)            ::   advectu,abcx
double precision,dimension (0:nx,1:ny)            ::   advectv,abcy
double precision,dimension (1:12,1:4)             ::   pxold2,pyold2,px,py,pxold,pyold 
double precision,dimension (0:nx+1,-1:ny+1)       ::   ugh,uold,uext,fximp,fximpold
double precision,dimension (-1:nx+1,0:ny+1)       ::   vgh,vold,vext,fyimp,fyimpold

double precision,dimension (-1:nx+1,-1:ny+1)      ::   phiold2,phioldc
double precision,dimension (1:nx-1,1:ny-1)        ::   tmp1,tmp2,mom,cmomx,cmomy,cbalast,fbox,fboy,fmoorx,fmoory,vort
double precision,dimension (0:nx,0:ny)            ::   ro,rt,p,miuv,Ib,Ic,Icold,Ibold,ro2,miuw,row,div2
double precision,dimension (0:nx,0:ny,iterate)    ::   presv 

real HV
integer i,j,k,kk,tp,count,it,ADVECT,plot,solver,viscose,tstep,number2,tpstar,TCount
integer isp(1:2),jsp(1:2),ksp(1:2),tetkind,numimp,ibdy,contact,sumnew,sumold
CHARACTER(LEN=6) :: NNumber
CHARACTER(LEN=19):: FileName


!! input parameters !!
call Platform_Constant_Ini() 
include 'par.txt'
!!!!!!!!!!!!!!!!!!!!!!!
Call meshgenerator()

!! intial values !! 
do i=1,2 
Lzero(i)=sqrt(  (ox(i)-px(i,4))**2 + (oy(i)-py(i,4))**2  )   !-0.01
end do 

eps=2*hy(ny/2)  !2*hy(ny/2)
gap=0.3*hy(ny/2)
count=0

!xbar0=( rosolid*( 0.5*(px(1,1)+px(2,2)) )+percent*rocon*( 0.5*(px(3,3)+px(2,2)) ) )/(rosolid+percent*rocon)
!ybar0=( rosolid*( 0.5*(py(1,1)+py(2,2)) )+percent*rocon*( 0.5*(py(3,3)+py(2,2)) ) )/(rosolid+percent*rocon)

xbar=lx/4.0  ; ybar=Ly/2.0

call Inicondition(xbar,ybar,p,Ib,Ic,Px,Py,gap,teta,yfree,dt,percent)
call boundarycond(0,dt,HH,period,Landa,hait,yfree,yfreeg,xgage)
presv=0
contact=1
!! end intial values !!

 
OPEN(25,file='Input_Piezo.plt') 
OPEN(35,file='result.plt')
OPEN(45,file='inside.plt')
OPEN(75,file='tetherforce.plt')
OPEN(105,file='moment.plt')
OPEN(125,file='neweq.plt')
open(135,file='boundary condition.plt')
open(145,file='boundary condition2.plt')



write(25,*) 'variables="t","ybar","Drag","Lift"'
write(35,*) 'variables="x","y","u","v","phi","solid","ro","p"'
write(45,*) 'variables="t","teta","sumphi","yfree","yfreeg","angularacc","accx","accy"'
write(75,*) 'variables="time","ftether1","ftether2"'
write(105,*) 'variables="tp","phit"'
write(125,*) 'variables="x","y","uext","vext","phi","ib"'
write(135,*) 'variables="x","y","uimp","vimp","ib"'
write(145,*) 'variables="time","sum"'

fximp=0 ; fyimp=0
sumu=0 ; sumv=0 ; Drag=0 ; Lift=0 
      !! main loop of the code !!! 
      do tp=1,tstep !time step
      print*,tp
      contact=contact+1
      !! parameters for Make solution second order in time !! 
      uold=u  ; vold=v  
      phiold2=phi ; Icold=Ic ; Ibold=Ib
      
      xbarold=xbar     ; ybarold=ybar   ; pxold=px ; pyold=py    
      ubarold=ubar     ; vbarold=vbar    
      omegazold=omegaz
      
      acgxold=acgx ; acgyold=acgy
      anacold=anac       
      tetaold=teta
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      fximpold=fximp
      fyimpold=fyimp
      
      sumuold=sumu ; sumvold=sumv  
      Dragold=Drag ; Liftold=Lift
      !! Second order time step !! 
      do torder=1,2
        tpstar=tp-1+torder
        ugh=u ; vgh=v  ; omegazgh=omegaz  ;  ubargh=ubar ; vbargh=vbar 
        
    

pxold2=px ;pyold2=py 
numimp=1 
xbarold2=xbar
ybarold2=ybar

divmax2=1

TCount=int(100000)+int(10)*tp+int(1)*torder 
write(unit=NNumber, fmt='(I6)') TCount
FileName='Iteration'//NNumber//'.plt'
open(unit=TCount,file=FileName) 
Write(TCount,*) 'variables="numimp","sumu","sumv","Div"'

 
 !! only large number to go inside the loop !!
!do while (numimp.le.iterate) !! Loop for iterating on solid correction !! 
sumugh=100000000
sumvgh=100000000
do while ( numimp.lt.iterate.AND.(abs(sumu-sumugh)+abs(sumv-sumvgh)).gt.0.001)  
 !! density and viscosity definition !! 
    ro(:,:)=1000
    miuv(:,:)=0.001  
    !call property(Ib,Ic,ro,rodrop,roair,rosolid,rocon,gap)
    !call property(Ib,Ic,miuv,miudrop,miuair,miusolid,miucon,gap)
     
!     do i=0,nx ;do j=0,ny 
!     ro(i,j)=roair+HV( -phi(i,j),eps )*(rodrop-roair)
!     ro(i,j)=ro(i,j)+Ib(i,j)*( rosolid-ro(i,j) )
!     ro(i,j)=ro(i,j)+Ic(i,j)*( rocon )
!     end do ;end do 
!     
!      
!     do i=0,nx ;do j=0,ny 
!     miuv(i,j)=miuair+HV( -phi(i,j),eps )*(miudrop-miuair)
!     miuv(i,j)=miuv(i,j)+Ib(i,j)*( miusolid-miuv(i,j) )
!     miuv(i,j)=miuv(i,j)+Ic(i,j)*( miucon )
!     end do ; end do 
      
!! writing the initial domian parameters!!       
if (tpstar.eq.1.AND.numimp.eq.1) then     
write(35,*) 'zone i=',nx-1,' j=',ny-1
 Do j=1,ny-1 ; Do i=1,nx-1
write(35,122) x(i),y(j),0.5*(u(i,J)+u(i+1,j)),0.5*(v(i,j)+V(i,j+1)),phi(i,j),Ib(i,j),ro(i,j),0.0
 122  format (8(1x,e15.7))
end do ; end do  
end if  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! advection and viscose terms in Navier-Stokes equation!!! 
call advection (advectu,advectv)
call viscosity (ro,miuv,Tx,Ty)
!call tether(px,py,x,y,nx,ny,Ox,Oy,ks,Fmoorx,Fmoory,Lzero,Ic,Lnew,xbar,ybar,ubar,xbar0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fbox=0 ; fboy=0 
do i=1,nx-1 ; do j=1,ny-1 
fbox(i,j)=0 ! fmoorx(i,j)
fboy(i,j)=0 ! fmoory(i,j)
end do ; end do 

sumold=sumnew
!call impactV(fximp,fyimp,u,v,ubar,vbar,xbar,ybar,omegaz,x,y,nx,ny,Ib,tp,numimp,r2,dt)
           



!!!!!!!!!!!!! initial guess for velocity without pressure term!!!!!!!!!!!!!
Do i=2,nx-1 ; Do j=1, ny-1 
U(i,j)=Ugh(i,j)+advectu(i,j)*dt+Tx(i,j)*dt+gx*dt+( fbox(i,j)+fximp(i,j) )*dt/(  0.5*(ro(i,j)+ro(i-1,j))  )/( hx(i)*hy(j) )

end do ; end do
Do j=2,ny-1 ; Do i=1,nx-1
V(i,j)=Vgh(i,j)+advectv(i,j)*dt+Ty(i,j)*dt+gy*dt+( fboy(i,j)+fyimp(i,j) )*dt/(  0.5*(ro(i,j)+ro(i,j-1))  )/( hx(i)*hy(j) ) 
end do ; end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


 
print*,"dt=",dt
print *, "alpha=",anac, "differecne of alpha=",anacexp-anac
anacexp=anac

call property(Ib,Ic,ro,rodrop,roair,rosolid,rocon,gap)
call boundarycond(tpstar,dt,HH,period,Landa,hait,yfree,yfreeg,xgage) !! ? !!


!!!!!!!!!!!!!!!!!!!!!!!!!!solve pressure SOR start !!!!!!!!!!!!!!!!!
!p(:,:)=presv(:,:,numimp)       !!! this is just  for better initial guess 
call poisson (ro,dt,pdif,p,beta)        
!presv(:,:,numimp)=p(:,:)       !!! For saving initial guess for pressure in next time step                                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!solve pressure SOR end !!!!!!!

 !!! start update velocity !!!!!!!!!!!!!!!!!!!
do i=2,nx-1 ; do j=1,ny-1 

U(i,j)=U(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )*dt/( 0.5*(hx(i)+hx(i-1))  )
end do ; end do 

do i=1,nx-1 ; do j=2,ny-1 

V(i,j)=V(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )*dt/( 0.5*(hy(j)+hy(j-1))  )
end do ; end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call boundarycond(tpstar,dt,HH,period,Landa,hait,yfree,yfreeg,xgage)
 







!!!!! computing maximum Divergence in the domain before  Solid correction in "Solid Subroutine !!!
divmax2=0
do i=2,nx-1 ;do j=2,ny-1 
div=(U(i+1,j)-U(i,j))/hx(i)+(v(i,j+1)-v(i,j))/hy(j)

if(abs(div).gt.divmax2)then
divmax2=abs(div)
end if 

end do ; end do 

WRITE(*,*)'DIV before  solid Correct=',DIVMAX2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Drag=0 ; Lift=0 
do i=2,nx-1 ; do j=2,ny-1 
Drag=Drag+ 0.5*(ro(i,j)+ro(i-1,j))*0.5*(Ib(i,j)+Ib(i-1,j))*(x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )* &
&          ( +advectu(i,j)-(  1/(  0.5*(ro(i,j)+ro(i-1,j))  )   )*( p(i,j)-p(i-1,j) )/( 0.5*(hx(i)+hx(i-1))  )+Tx(i,j)+gx )

Lift=Lift+ 0.5*(ro(i,j)+ro(i,j-1))*0.5*(Ib(i,j)+Ib(i,j-1))*(y(j)-y(j-1))*( 0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )* &
&          ( +advectv(i,j)-(  1/(  0.5*(ro(i,j)+ro(i,j-1))  )   )*( p(i,j)-p(i,j-1) )/( 0.5*(hy(j)+hy(j-1))  )+Ty(i,j)+gy )

end do ; end do 
sumugh=sumu
sumvgh=sumv


 !! immersed boundary method for solid correction !!    
call solid(ro,Ib,Ic,xbar,ybar,ubar,vbar,dt,omegaz,Px,Py,ox,oy,lzero,&
  &     rosolid,rocon,gap,gy,Tmass,teta,percent,pxold2,pyold2,tp,numimp,fximp,fyimp,xbarold2,ybarold2,sumu,sumv)
 !! End immersed boundary method for solid correction !! 
!!!!! computing maximum Divergence in the domain after Solid correction in "Solid Subroutine !!!
divmax2=0
do i=2,nx-1 ;do j=2,ny-1 
div=(U(i+1,j)-U(i,j))/hx(i)+(v(i,j+1)-v(i,j))/hy(j)

if(abs(div).gt.divmax2)then
divmax2=abs(div)
end if 

end do ; end do 

WRITE(*,*)'DIVERGENCEiterate=',DIVMAX2 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

acgx=(ubar-ubargh)/dt ; acgy=(vbar-vbargh)/dt 

anac=0.1*(omegaz-omegazgh)/dt + 0.9*anac
!anac=(omegaz-omegazgh)/dt


Write(TCount,*)  numimp,sumu,sumv,divmax2
numimp=numimp+1

end do !! end of iteration !!

if (numimp.lt.5) then 
  close (unit=TCount,STATUS='delete') 
else 
  close (unit=TCount) 
end if 





 contact=ibdy
!if (contact.eq.ibdy) then 

!end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!! Solving Level set equation !!!!!!!!!!!!!! 
!call levelset(phioldc,dt,numimp)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!! correct interface position for applying proper contact angle !!!!!!!!!!!!!!
!!if (contact.eq.ibdy) then 
!call contactangle(uext,vext,dt,Ib,gap)
!call levelset2(uext,vext,0.25*dt,tp,ib)
!contact=0
!!endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!!!!!!!!!!!!! Solving reinitialization equation for keeping Level set function as a distance function !!!!!!!!!!!   
!call reinitialize(eps,dt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!! contact agnle modeling !!!!!!!!!!!!!!!
!


end do   !! torder !!
!!!!!!!!!!!!!!!!!!! averaging to consequence time step for reaching second order in time for velocity field !!!!!!!!!!!

!! fluid values !!
u=0.5d0*(u+uold)                   ; v=0.5d0*(v+vold)  
!! property  of solid and fluid (density and viscosity)               
phi=0.5d0*(phi+phiold2)           
Ic=0.5d0*(Icold+Ic)  ; Ib=0.5d0*(Ibold+Ib)
!!!!
!! solid and tether positions/velocity/accelaration !! 
xbar=0.5d0*(xbarold+xbar)          ; ybar=0.5d0*(ybarold+ybar)  
px=0.5d0*(pxold+px)                ; py=0.5d0*(pyold+py)         
ubar=0.5d0*(ubarold+ubar)          ; vbar=0.5d0*(vbarold+vbar)         
omegaz=0.5d0*(omegazold+omegaz)
 
acgx=0.5*(acgxold+acgx)          ; acgy=0.5*(acgyold+acgy)
anac=0.5*(anacold+anac)   
!! just for plotting !!!       
 teta=0.5*(tetaold+teta)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

fximp=0.5*(fximpold+fximp)
fyimp=0.5*(fyimpold+fyimp)

sumu=0.5*(sumu+sumuold) 
sumv=0.5*(sumv+sumvold)

Drag=0.5*(Drag+Dragold)
Lift=0.5*(Lift+Liftold)




!! Post processing start !! 


!!!!!!!!!!!!!!!!!! computing divergence at the end of marching 1 step in time !!!!!!!!!!!!!!!!!!!!!!!!!!!
divmax2=0
div2=0
do i=2,nx-1 ;do j=2,ny-1 
div2(i,j)=(U(i+1,j)-U(i,j))/hx(i)+(v(i,j+1)-v(i,j))/hy(j)

!if(abs(div2(i,j)).gt.divmax2)then
!divmax2=abs(div2(i,j))
!else
!end if 

end do ; end do
WRITE(*,*)'DIVERGENCE Second order in time =', divmax2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!! writing  1 !!!!!!!!!!!!!!!!!!!!
tetaplot=180/pi*ACOS( ( px(2,2)-px(1,1) )/sqrt( (px(2,2)-px(1,1))**2+( py(2,2)-py(1,1) )**2+0.000001 )  )
!write (25,*) tp*dt,xbar,ybar,tetaplot,omegaz,Tmass
 !write (25,*) tp*dt,xbar,ybar,tetaplot,(sumu-0.0)/dt/(0.5*rodrop*(0.375**2)*2*r2),(sumv-0.0)/dt/(0.5*rodrop*(0.375**2)*2*r2)
! write (25,*) tp*dt,ybar,(sumu)/dt/(0.5*rodrop*(0.375**2)*2*r2),(sumv)/dt/(0.5*rodrop*(0.375**2)*2*r2),&
! &            Drag/(0.5*rodrop*(0.375**2)*2*r2),Lift/(0.5*rodrop*(0.375**2)*2*r2)
 write (25,133) tp*dt,ybar,Drag/(0.5*rodrop*(0.375**2)*2*r2),Lift/(0.5*rodrop*(0.375**2)*2*r2)

133  format (4(1x,e15.7))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i=1,nx-1 ; do j=1,ny-1 
vort(i,j)= 0.25*(( v(i,j)    -v(i-1,j)   )/(x(i)-x(i-1))   &
        &       +( v(i+1,j)  -v(i,j)     )/(x(i+1)-x(i))   &
        &       +( v(i,j+1)  -v(i-1,j+1) )/(x(i)-x(i-1))   &
        &       +( v(i+1,j+1)-v(i,j+1)   )/(x(i+1)-x(i))  )&
        & -0.25*(( u(i,j)    -u(i,j-1)   )/(y(j)-y(j-1))   &
        &       +( u(i+1,j)  -u(i+1,j-1) )/(y(j)-y(j-1))   &
        &       +( u(i,j+1)  -u(i,j)     )/(y(j+1)-y(j))   &
        &       +( u(i+1,j+1)-u(i+1,j)   )/(y(j+1)-y(j))  ) 
        
end do ; end do  

!!!!!!!!!!!!!!! writing  2 !!!!!!!!!!!!!!!!!!!!
count=count+1
if (count.eq.plot.AND.DIVMAX2.lt.500)then 
  count=0
  print*,"data is wrriten"
  write(35,*) 'zone i=',nx-1,' j=',ny-1
  Do j=1,ny-1 ; Do i=1,nx-1
    write(35,122) x(i),y(j),0.5*(u(i,J)+u(i+1,j)),0.5*(v(i,j)+V(i,j+1)),vort(i,j),Ib(i,j),ro(i,j),abs(div2(i,j))
  end do ; end do 

end if 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!! writing  3 !!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!! checking mass conservation in Level Set method !!!!!!!!!!!!!!!!!!!
sumphi=0 
do i=1,nx-1 ; do j=1,ny-1 
if (phi(i,j).gt.0) then 
sumphi=sumphi+1
end if 
end do ; end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
write(45,*) tp*dt,teta*180/pi,sumphi,yfree,yfreeg,anac,acgx,acgy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!!!!!! write for tethers 4 (not now )!!!!!!!!!!!!!!!!!!!!!!!
write (75,*) tp*dt, ( max(Lzero(1),Lnew(1))-Lzero(1) ),( max(Lzero(2),Lnew(2))-Lzero(2) )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
!! Post processing End  !! 


 

end do !time step !! end of the main loop of the code !! 

write(*,*)'end'
read (*,*)

    end program 

!!!!!!!!!!!!!!!!!!!!!!! subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




Subroutine Inicondition(xbar,ybar,p,Ib,Ic,Px,Py,gap,teta,yfree,dt,percent)
use M_General_2D,  only: phi,u,v,x,y,nx,ny,lx,pi
use M_Platform_Constant_2D, only:r2,len,floatx
use M_Mesh_2D   ,  only: hx,hy


implicit none 
integer i,j,tpstar,tp,ipr

double precision p(0:nx,0:ny)
double precision Px(1:12,1:4),Py(1:12,1:4)
double precision Ic(0:nx,0:ny),Ib(0:nx,0:ny)
double precision gap,dt,teta,yfree,percent
double precision xphi(0:2500),yphi(0:2500),phiges,xs(4),ys(4),xbar,ybar



Ib(0:nx,0:ny)=0
Ic(0:nx,0:ny)=0
phi=10000


!! constructing initial Level set values (Level set goes around the obeject !!!!!! 
do i=0,500 
xphi(i)=( floatx-r2/sin(75*pi/180) )*dble(i)/500
yphi(i)=yfree 
end do 

do i=1,500   
xphi(i+500)=( floatx+r2/sin(75*pi/180) ) +( Lx-(floatx+r2/sin(75*pi/180))  )*dble(i)/500
yphi(i+500)=yfree 
end do 


xs(1)= floatx-r2/sin(90*pi/180) 
ys(1)=yfree


xs(2)=px(2,2)-r2*cos(0*pi/180)
ys(2)=py(2,2)-r2*sin(0*pi/180)


xs(3)=xs(2)+2*r2*cos(0*pi/180)
ys(3)=ys(2)+2*r2*sin(0*pi/180)

xs(4)=floatx+r2/sin(90*pi/180) 
ys(4)=yfree


do i=1,500
xphi(i+1000)=xs(1)+(xs(2)-xs(1))*dble(i)/500
yphi(i+1000)=ys(1)+(ys(2)-ys(1))*dble(i)/500
end do 

do i=1,500
xphi(i+1500)=xs(2)+(xs(3)-xs(2))*dble(i)/500
yphi(i+1500)=ys(2)+(ys(3)-ys(2))*dble(i)/500
end do 

do i=1,500
xphi(i+2000)=xs(3)+(xs(4)-xs(3))*dble(i)/500
yphi(i+2000)=ys(3)+(ys(4)-ys(3))*dble(i)/500
end do 



OPEN(85,file='point.plt')

write(85,*) 'variables="xphi","yphi"'

do i=0,2500 
write (85,*) xphi(i),yphi(i)
end do 


do i=1,nx-1 ; do j=1,ny-1 

do ipr=0,2500
phiges=sqrt(   ( x(i)-xphi(ipr) )**2  +  (  y(j)-yphi(ipr)  )**2   )

if (phiges.lt.phi(i,j) ) then 
phi(i,j)=phiges
end if 

end do
phiges=1000

end do ; end do 



do i=1,nx-1 ; do j=1,ny-1 

if (x(i).lt.xs(1)) then
          if ( y(j).gt.yfree ) then 
          phi(i,j)=-phi(i,j)
          end if 
end if           
!! I made it comment for not initial angle 
!if (x(i).ge.xs(1).AND.X(i).lt.xs(2) ) then  
! 
!          if (y(j).gt.(ys(1)-tan(75*pi/180)*(x(i)-xs(1)) )   ) then 
!          phi(i,j)=-phi(i,j)
!          end if
!end if
          
if (x(i).ge.xs(2).AND.X(i).lt.xs(3) ) then  
          
         if (y(j).gt.(ys(2)+tan(0*pi/180)*(x(i)-xs(2)) )   ) then 
         phi(i,j)=-phi(i,j)
         end if 
end if 
!! I made it comment for not initial angle           
!if (x(i).gt.xs(4).AND.x(i).le.xs(3) ) then 
!         if (y(j).gt.(ys(3)-tan(75*pi/180)*(x(i)-xs(3)) )   ) then 
!         phi(i,j)=-phi(i,j)
!         end if         
!end if 
if (x(i).gt.xs(4)) then
          if ( y(j).gt.yfree ) then 
          phi(i,j)=-phi(i,j)
          end if  
end if

!if (x(i).lt.xs(4).AND.x(i).ge.xs(3).AND.y(j).lt.yfree ) then
!         if (y(j).gt.(ys(3)-tan(75*pi/180)*(x(i)-xs(3)) )   ) then 
!         phi(i,j)=-phi(i,j)
!         end if 
!end if  



 
end do ; end do



phi=-phi


!! end of constructing Level set values !! 


!! find to grids occupied by the  solid !! 

!call insidef(px(1,1),px(2,2),py(1,1),py(2,2),Ib,x,y,r2,gap,1.0/1.0*len,nx,ny,dt,teta)
!call insidef(px(3,3),px(2,2),py(3,3),py(2,2),Ic,x,y,r2,gap,percent*len,nx,ny,dt,teta)
call insidef2(xbar,ybar,Ib,x,y,r2,gap,nx,ny)

v=0 ; U=0 ; p=0 !!initial values for velocity and pressure 



return 
end 



subroutine boundarycond(tpstar,dt,HH1,period,Landa,hait,yfree,yfreeg,xgage)
use M_General_2D,  only : phi,u,v,nx,ny,x,y,pi
use M_Mesh_2D,     only : hx,hy       
implicit none 
integer i,j,tpstar,imax,jmax,Cb,igage
double precision dt
double precision period,landa,hait,HH,vala1,vala2,yfree,kwave,wwave,sumvv,wavegen
double precision xgage,yfreeg,HH1
!! boundary condition for Velocity (Including generating wave ) and Level Set function !! 
!! generated wave from Stokes second order solution  but only velocity in x direction is given velocity in y direction will be 
!! satisifed by mass conservation  
!! inflow out flow on left and top !! 
!! Slip boundary condition on the others !! 

wavegen=4

if ( (tpstar*dt).lt.(3*period) ) then 
HH=(tpstar*dt)/(3*period)*HH1
else 
HH=HH1
end if 
!!!!!!!!!!!!!!!!!! finding the wave height at the  wave gage place !!!!!!!!!!!!!!!!!
igage=1 !! for problems with out free surface 
do  i=1,nx-1
 
  if (x(i).le.xgage.AND.x(i+1).ge.xgage) then 
     igage=i
     exit 
  end if 
end do 
yfreeg=0
do j=1,ny-1
  if (0.5*( phi(igage+1,j)+phi(igage,j) ).le.0.AND.0.5*( phi(igage+1,j+1)+phi(igage,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
    yfreeg=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(igage+1,j+1)+phi(igage,j+1) )-0.5*( phi(igage+1,j)+phi(igage,j) ) ) & 
    &    *( -0.5*( phi(igage+1,j)+phi(igage,j) ) ) 
    exit 
  end if 
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


 
!!!!!!!!! indicator for error !!!!!!!
if (yfreeg.eq.0)then 
print*,"errorgage" 
end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!! boundary condition for Level set function only extrapolating !!!!!!!!!!!!!!!!!
do j=1,ny-1 
phi(0,j)=2*phi(1,j)-phi(2,j)
phi(-1,j)=3*phi(1,j)-2*phi(2,j)

phi(nx,j)=2*phi(nx-1,j)-phi(nx-2,j)
phi(nx+1,j)=3*phi(nx-1,j)-2*phi(nx-2,j)
end do 

do i=-1,nx+1 
phi(i,0)=2*phi(i,1)-phi(i,2)
phi(i,-1)=3*phi(i,1)-2*phi(i,2)

phi(i,ny)=2*phi(i,ny-1)-phi(i,ny-2)
phi(i,ny+1)=3*phi(i,ny-1)-2*phi(i,ny-2)
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

kwave=2*pi/landa
wwave=2*pi/period

!!!!!!!!!!!!!!! tracking the free surface place in left hand side of the domain for using in wave generator !!!!!!!!!!!!!
yfree=0
do j=1,ny-1
  if (0.5*( phi(1,j)+phi(0,j) ).le.0.AND.0.5*( phi(1,j+1)+phi(0,j+1) ).ge.0 ) then             !! may be more precsiness can be done !! 
    yfree=y(j)+ ( y(j+1)-y(j) )/ ( 0.5*( phi(1,j+1)+phi(0,j+1) )-0.5*( phi(1,j)+phi(0,j) ) ) *( -0.5*( phi(1,j)+phi(0,j) ) ) 
    exit 
  end if 
end do  
if (yfree.eq.0)then 
print*,"error" 
end if





     
!!!!!! analytical solution for second order wave, the velocity for the air is fake (only for continuty )
sumvv=0
if (wavegen.eq.1 ) then 

  do j=1 ,ny-1 

    if (phi(0,j).lt.0 ) then

    u(1,j)=HH/2* wwave *cos(-wwave*(tpstar*dt))*cosH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave*hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
 
    u(0,j)=HH/2* wwave*cos(-wwave*(tpstar*dt))*cosH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
                        
    else

    u(1,j)=HH/2* wwave *cos(-wwave*(tpstar*dt))*cosH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave*hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
 
    u(0,j)=HH/2* wwave*cos(-wwave*(tpstar*dt))*cosH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
    & 3/4*((pi*HH)**2)/(period*landa)*cos(-2*wwave*(tpstar*dt))*cosH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
    
    end if 
 sumvv=sumvv+u(1,j)*hy(j)*1.0d0 

  !U(1,j)=0      
  !U(0,j)=-U(2,j)
   u(nx,j)=0 
   u(nx+1,j)=-u(nx-1,j) 
  end do 

else if (wavegen.eq.2) then 

   do j=1 ,ny-1 

     if (phi(0,j).gt.0 ) then
     
     u(1,j)=HH*y(j)*sin(kwave*0-wwave*tpstar*dt)      
     u(0,j)=HH*y(j)*sin(kwave*-2*x(1)-wwave*tpstar*dt)
     
     else 
     
     u(1,j)=HH*( 2*yfree-y(j)   )*sin(kwave*0-wwave*tpstar*dt)     
     u(0,j)=HH*( 2*yfree-y(j)   )*sin(kwave*-2*x(1)-wwave*tpstar*dt)
     
     end if 
   sumvv=sumvv+u(1,j)*hy(j)*1.0d0
   
   !U(1,j)=0      
   !U(0,j)=-U(2,j)
   u(nx,j)=0 
   u(nx+1,j)=-u(nx-1,j) 
  end do
else if (wavegen.eq.3) then 
 
  do j=1 ,ny-1 
   U(1,j)=0      
   U(0,j)=-U(2,j)
   u(nx,j)=0 
   u(nx+1,j)=-u(nx-1,j) 
  end do
  
 else 
 do j=1 ,ny-1 
   u(1,j)   =0.15 !0.000375 !u(nx-1,j)      
   u(0,j)   =0.15 !0.000375 !u(nx-2,j)
   u(nx,j)  =u(nx-1,j) 
   u(nx+1,j)=u(nx-2,j) 
  end do
  
  

end if 

!! at first smaller and then use them for boundary of bigger !!!


!!!slip!! 
Do i=0, nx+1
U(i,0)=U(i,1)
U(i,-1)=U(i,2)
U(i,ny)=U(i,ny-1)
U(i,ny+1)=U(i,ny-2)
end do

!no-slip!! cavity 
!Do i=0, nx+1
!U(i,0)=-U(i,1)
!U(i,-1)=-U(i,2)
!U(i,ny)=-U(i,ny-1)  !!1.0 
!U(i,ny+1)=-U(i,ny-2) !!1.0
!end do



if (wavegen.eq.1.OR.wavegen.eq.2) then
   Do i=1, nx-1 
   V(i,1)=0
   V(i,0)=-V(i,2)
   V(i,ny)=sumvv/( real(nx-1) *hx(i) )         !!!!!!!!!!!!!inflow out flow  boundary condition !!!!!!!!!!!!!!
   V(i,ny+1)=V(i,ny)
   end do
else
   Do i=1, nx-1 
   V(i,1)=0
   V(i,0)=-V(i,2)
   V(i,ny)=0
   V(i,ny+1)=-V(i,ny-1)
   end do
end if 



if (wavegen.eq.5) then
 
   do j=0,ny+1
   
   if (phi(0,j).lt.0 ) then
   v(0,j)=HH/2*wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )

   v(-1,j)=HH/2* wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( (y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* ((y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
   else 
   v(0,j)=HH/2*wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )

   v(-1,j)=HH/2* wwave *sin(-wwave*(dble(tpstar)*dt))*sinH(kwave* ( -(y(j)-yfree)+hait) )/sinH(kwave *hait)+ &
   & 3/4*((pi*HH)**2)/(period*landa)*sin(-2*wwave*(tpstar*dt))*sinH( 2*kwave* (-(y(j)-yfree)+hait) )/( (sinH(kwave*hait))**4 )
   end if
   
   V(nx,j)=V(nx-1,j)
   V(nx+1,j)=V(nx-2,j)
   
   end do 
else

!! slip !!
!   do j=0,ny+1 
!   
!   V(0,j)=V(1,j)
!   V(-1,j)=V(2,j)
!   V(nx,j)=V(nx-1,j)
!   V(nx+1,j)=V(nx-2,j)
!   
!   end do

!! No-Slip cavity !!
!do j=0,ny+1 
!   
!   V(0,j)=-V(1,j)
!   V(-1,j)=-V(2,j)
!   V(nx,j)=-V(nx-1,j)
!   V(nx+1,j)=-V(nx-2,j)
!   
! end do

! periodic !!!!!!!!!
!do j=0,ny+1 
!   
!   V(nx,j)=V(1,j)
!   V(nx+1,j)=V(2,j)
!   V(0,j)=V(nx-1,j)
!   V(-1,j)=V(nx-2,j)
!   
!   end do

!! out flow 
do j=0,ny+1 
   
   V(0,j)=V(1,j)
   V(-1,j)=V(2,j)
   V(nx,j)=V(nx-1,j)
   V(nx+1,j)=V(nx-2,j)
   
   end do
end if    
 


return 
end 

subroutine levelset(phioldc,dt,numimp)
use M_General_2D,  only:phi,u,v,nx,ny 
use M_Mesh_2D,     only: hx,hy                  
implicit none 
Double precision phioldc(-1:nx+1,-1:ny+1)
double precision, DIMENSION (0:nx,0:ny)       ::Dpphix,Dmphix,Dpphiy,Dmphiy
double precision, DIMENSION (1:nx-1,1:ny-1)   ::phix,phiy,Lphin
Double precision Lphis,dt
integer i,j,kk,numimp

!if (numimp.eq.2) then 
phioldc(1:nx-1,1:ny-1)=phi(1:nx-1,1:ny-1)
!end if 


!! Second ordeer ENO convective terms for Advection terms of Level set equation !! 
!! Second order in time !! 

do kk=1,2  !!prediction correction method!!


do i=0,nx ; do j=0,ny                                
Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  0.5*(hx(i)+hx(i+1))  )
Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  0.5*(hx(i)+hx(i-1))  )
Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  0.5*(hy(j)+hy(j+1))  )
Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  0.5*(hy(j)+hy(j-1))  )

end do ;end do 

do i=1,nx-1 ; do j=1,ny-1    !!A!!

if (0.5*( u(i,j)+u(i+1,j) ).gt.0.d0) then


  if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
phix(i,j)=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
   else
phix(i,j)=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
   end if 
    
else

   if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
phix(i,j)=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
   else
phix(i,j)=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
   end if 

end if 



if (0.5*( V(i,j)+V(i,j+1) ).gt.0.d0) then


   if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
phiY(i,j)=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
   else
phiY(i,j)=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
   end if 
    
else

  if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
phiY(i,j)=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
   else
phiY(i,j)=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
   end if 


end if 


     
end do ;end do  !!A!!

if (kk.eq.1) then

do i=1,nx-1 ; do j=1,ny-1 
Lphin(i,j)=  -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
         &   -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
           
phi(i,j)=phioldc(i,j)+dt*( -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
                     & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) )
                
end do ;end do 

                  
end if                        

   if (kk.eq.2) then
   
   do i=1,nx-1 ; do j=1,ny-1 
 
Lphis=   -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)&
       & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
            
phi(i,j)=phioldc(i,j)+0.5*dt*(Lphis+Lphin(i,j) )                       
   end do ;end do 

   end if 


end do  !!prediction corection method!! 
 
call boundarycondphi()                    

RETURN
END 


subroutine reinitialize(eps,dt)
use M_General_2D,  only:phi,u,v,nx,ny,pi 
use M_Mesh_2D,     only: hx,hy    
implicit none
Double precision eps
Double precision phixm,phiym,phixp,phiyp,phixr,phiyr,sphi,s,dtau,AA
double precision, Dimension (0:nx,0:ny)       ::  Dpphix,Dmphix,Dpphiy,Dmphiy
double precision, Dimension (1:nx-1,1:ny-1)   ::  phiold,phin,lphin
real sgn
Double precision lphis,bb,avgh,dt,HH,period,Landa,hait
integer i,j,kk,m

! Ref; Journal of computational physics vol. 152 pp-493-516 (1999)

phiold(1:nx-1,1:ny-1)=phi(1:nx-1,1:ny-1)
avgh=1000 !! large number!! 
do i=-1,nx+1 ; do j=-1, ny+1 
bb=min (hx(i),hy(j))
if (bb.lt.avgh) then
avgh=bb
end if 
end do ; end do 

dtau=0.5*avgh  !! fictious time step !


do kk=1,6

do m=1,2   !!prediction correction method for time step the same as Level set for advection terms  !! 


do i=0,nx ; do j=0,ny 
Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  0.5*(hx(i)+hx(i+1))  )
Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  0.5*(hx(i)+hx(i-1))  )
Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  0.5*(hy(j)+hy(j+1))  )
Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  0.5*(hy(j)+hy(j-1))  )

end do ;end do 

do i=1,nx-1 ; do j=1,ny-1 !!A!!  !!A!!


!sphi=phiold(i,j)/(  dsqrt(phiold(i,j)*phiold(i,j)+h*h)  )
!sphi=sgn(phiold)



if (phiold(i,j).gt.eps) then
sphi=1.d0
else if (phiold(i,j).lt.-eps) then
sphi=-1.d0
else
sphi=phiold(i,j)/eps -(1/pi)*dsin(pi*phiold(i,j)/eps) 
end if 

 

   if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
phixm=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
   else
phixm=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
   end if 
    

    if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
phixp=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
   else
phixp=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
   end if 



    if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
phiYm=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
   else
phiYm=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
   end if 
    


   if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
phiYp=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
   else
phiYp=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
   end if
  
  
   

   
        if (sphi*phixp.ge.0.d0.AND.sphi*phixm.ge.0.d0) then
   phixr=phixm
   else if (sphi*phixp.le.0.d0.AND.sphi*phixm.le.0.d0) then 
   phixr=phixp
   else if (sphi*phixp.gt.0.d0.AND.sphi*phixm.lt.0.d0) then
   phixr=0.0
   else if (sphi*phixp.lt.0.d0.AND.sphi*phixm.gt.0.d0) then
   s=sphi*( abs(phixp)-abs(phixm) )/(phixp-phixm)
                      if (s.gt.0.0) then
                      phixr=phixm
                   else
                      phixr=phixp 
                   end if
   end if
   
   
        if (sphi*phiyp.ge.0.d0.AND.sphi*phiym.ge.0.d0) then
   phiyr=phiym
   else if (sphi*phiyp.le.0.d0.AND.sphi*phiym.le.0.d0) then 
   phiyr=phiyp
   else if (sphi*phiyp.gt.0.d0.AND.sphi*phiym.lt.0.d0) then
   phiyr=0.0
   else if (sphi*phiyp.lt.0.d0.AND.sphi*phiym.gt.0.d0) then
   s=sphi*( abs(phiyp)-abs(phiym) )/(phiyp-phiym)
                      if (s.gt.0.0) then
                      phiyr=phiym
                   else
                      phiyr=phiyp 
                  end if
    end if
    
        
   
   
   if (m.eq.1) then
   lphin(i,j)=sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
   phin(i,j)=phi(i,j)
   phi(i,j)=phi(i,j)+dtau*sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
   
   end if
   
   
   if (m.eq.2) then
   lphis   =sphi*(  1.d0-dsqrt(phiyr*phiyr+phixr*phixr)  )
   phi(i,j)=phin(i,j)+0.5*dtau*(  lphis+lphin(i,j)  )                   
   end if 

end do ; end do 

call boundarycondphi() 


end do 



end do !!fictious time step!!

return 
END

subroutine solid(ro,Ib,Ic,xbar,ybar,ubar,vbar,dt,omegaz,Px,Py,ox,oy,lzero, &
&          rosolid,rocon,gap,gy,Tmass,teta,percent,pxold2,pyold2,tp,numimp,uimp,vimp,xbarold2,ybarold2,sumu,sumv)
use M_General_2D          ,  only:u,v,x,y,nx,ny,Ly
use M_Platform_Constant_2D,  only:r2,len
use M_mesh_2D             ,  only:hx,hy 
implicit none
integer nz,tp
double precision,intent (in)                          :: gy,percent,dt,rocon,rosolid,gap
double precision,intent (in)   ,dimension(1:3)        :: Lzero,ox,oy
double precision,intent (in)                          :: uimp(0:nx+1,-1:ny+1),vimp(-1:nx+1,0:ny+1)
double precision,intent (in)   ,dimension(0:nx,0:ny)  :: ro
double precision,intent (out)                         :: Tmass
double precision,intent (inout)                       :: ubar,vbar,teta,omegaz,xbar,ybar,xbarold2,ybarold2,sumu,sumv
double precision,intent (inout),dimension(1:12,1:4 )  :: pxold2,pyold2,px,py
double precision,intent (inout),dimension(0:nx,0:ny)  :: Ib,Ic


!! local variables !!

double precision,dimension (1:12,1:4)        ::pxn,pyn
double precision Ixx,Iyy,Ixy,Izz, Izzb,sumIzb,AAA,BBB,CCC,mom(1:nx-1,1:ny-1),pi
double precision AAA2,BBB2
integer i,j,k,numimp
!!!!!!!!!!!!!!!!!!!!!!!


pi=3.14159265


sumu =0 ;sumv =0 ;   Izz =0 ; sumIzb=0 ;    !!!sumIx=Hx , SumIy=Hy, SumIz=Hz !!
 Izzb=0   


!if (numimp.eq.1) then 
 
!! tracking the two points of the solid !! 
Pxn(1,1)=Pxold2(1,1) +dt*( ubar   -omegaz*(Py(1,1)-ybar) )
Pyn(1,1)=Pyold2(1,1) +dt*( vbar   +omegaz*(Px(1,1)-xbar) )

Pxn(2,2)=Pxold2(2,2) +dt*( ubar   -omegaz*(Py(2,2)-ybar) )
Pyn(2,2)=Pyold2(2,2) +dt*( vbar   +omegaz*(Px(2,2)-xbar) )

Px(1,1)=Pxn(1,1); Px(2,2)=Pxn(2,2)
py(1,1)=pyn(1,1); py(2,2)=pyn(2,2) 

!! for none unifrom solid !! 
px(3,3)= (1-percent)*px(2,2)+percent*px(1,1) 
py(3,3)= (1-percent)*py(2,2)+percent*py(1,1) 

xbar=xbarold2+ubar*dt
ybar=ybarold2+vbar*dt


!! tether position tracking !!
do i=1,2 
 Pxn(i,4)=Pxold2(i,4)+dt*( ubar   -omegaz*(Py(i,4)-ybar) )
 Pyn(i,4)=Pyold2(i,4)+dt*( vbar   +omegaz*(Px(i,4)-xbar) )
 
end do

do i=1,2 
 Px(i,4)=Pxn(i,4)
 py(i,4)=pyn(i,4)  
end do




 
!! solid CG tracking !! 

!xbar=( rosolid*( 0.5*(px(1,1)+px(2,2)) )+percent*rocon*( 0.5*(px(3,3)+px(2,2)) ) )/(rosolid+percent*rocon)
!ybar=( rosolid*( 0.5*(py(1,1)+py(2,2)) )+percent*rocon*( 0.5*(py(3,3)+py(2,2)) ) )/(rosolid+percent*rocon)

!! find the grids occupied by the solid !!! 
!call insidef(px(1,1),px(2,2),py(1,1),py(2,2),Ib,x,y,r2,gap,1.0/1.0*len,nx,ny,dt,teta)
!call insidef(px(3,3),px(2,2),py(3,3),py(2,2),Ic,x,y,r2,gap,percent*len,nx,ny,dt,teta)
call insidef2(xbar,ybar,Ib,x,y,r2,gap,nx,ny)


!call extforce(px,py,ox,oy,sumextu,sumextv,sumextIzb,lzero,anac,acgx,acgy,xbar,ybar,tp,dt)
         

!! calculating the summation of linear and Angular momentum of the Structure !! 
!end if 
 
         

mom=0 ; Tmass=0 

do i=1,nx-1 ; do j=1,ny-1 



!AAA=hx(i)*hy(j)*  0.5d0*( u(i,j)+u(i+1,j) )*(rocon*Ic(i,j)+rosolid*Ib(i,j))  
!BBB=hx(i)*hy(j)*  0.5d0*( v(i,j)+v(i,j+1) )*(rocon*Ic(i,j)+rosolid*Ib(i,j))     
!AAA=hx(i)*hy(j)*  0.5d0*( u(i,j)+u(i+1,j) )*(ro(i,j)*Ib(i,j))  
!BBB=hx(i)*hy(j)*  0.5d0*( v(i,j)+v(i,j+1) )*(ro(i,j)*Ib(i,j))     
AAA=(x(i)-x(i-1))*( 0.5*(y(j)+y(j+1))-0.5*(y(j)+y(j-1)) )*u(i,j)*0.5*(ro(i-1,j)+ro(i,j))*0.5*(Ib(i-1,j)+Ib(i,j)) 
BBB=(y(j)-y(j-1))*( 0.5*(x(i)+x(i+1))-0.5*(x(i)+x(i-1)) )*v(i,j)*0.5*(ro(i,j-1)+ro(i,j))*0.5*(Ib(i,j-1)+Ib(i,j))     

AAA2=0 !hx(i)*hy(j)*  0.5d0*( uimp(i,j)+uimp(i+1,j) )*(rocon*Ic(i,j)+rosolid*Ib(i,j))  
BBB2=0 !hx(i)*hy(j)*  0.5d0*( vimp(i,j)+vimp(i,j+1) )*(rocon*Ic(i,j)+rosolid*Ib(i,j))  

sumu=sumu+AAA-AAA2 
sumv=sumv+BBB-BBB2 

!Tmass=Tmass+ hx(i)*hy(j)*(rocon*Ic(i,j)+rosolid*Ib(i,j)) 
Tmass=Tmass+ hx(i)*hy(j)*(ro(i,j)*Ib(i,j)) 
mom(i,j)=(x(i)-xbar)*(BBB-BBB2) -(y(j)-ybar)*(AAA-AAA2) 

sumIzb=sumIzb+ mom(i,j)  
  
!Izzb=Izzb+ ( (x(i)-xbar)**2+ (y(j)-ybar)**2 )*hx(i)*hy(j)*(rocon*Ic(i,j)+rosolid*Ib(i,j))   
Izzb=Izzb+ ( (x(i)-xbar)**2+ (y(j)-ybar)**2 )*hx(i)*hy(j)**(ro(i,j)*Ib(i,j))   
 


end do ; end do 



!! calculating linear and angular velocity of the solid !! 

omegaz=0 ! (sumIzb)/(Izzb)


ubar  =0 !(sumu)/Tmass 
vbar  =0 !(sumv)/Tmass  ! +dt/Tmass*(-125*(Ybar-Ly/2))

print*,"mass and moment",  Tmass,Izzb



!! Correcting Navier stokes prediction for velocity field !!                                                                
do i=2,nx-1 ; do j=2,ny-1 

u(i,j)=u(i,j)+0.5d0*( Ib(i-1,j)+Ib(i,j) )*( ( ubar  -omegaz*(y(j)-ybar) )-u(i,j) ) !! mistake solved!!
v(i,j)=v(i,j)+0.5d0*( Ib(i,j-1)+Ib(i,j) )*( ( vbar  +omegaz*(x(i)-xbar) )-v(i,j) )

end do ; end do !! be careful if the solid is near boundry then boundry condition of velocity should change !!! 


print*,"omega=",omegaz

return
 End
 

!!!!!!!!!!!!! smoothed heaviside function !!!!!!!!!!!!!
function HV(phi,eps)

double precision phi,pi,eps


pi=3.141592654D0
if (phi.gt.eps) then
HV=1.d0
else if (phi.lt.-eps) then
HV=0.d0
else
HV=0.5d0*(  1.d0+ phi/eps +(1/pi)*dsin(pi*phi/eps)  )
end if 

return
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!subroutine for finding weather a point is inside the solid or not !!!
subroutine insidef(x1,x2,y1,y2,Ib1,x,y,rr,gap,len1,nx,ny,dt,teta)   
 
implicit none
double precision x1,x2,y1,y2,xbar1,ybar1,Ib1(0:nx,0:ny),x(-1:nx+1),y(-1:ny+1),gap,len1
double precision ax,ay,bx,by,rr,d,tt,cx,cy,Ib2(0:nx,0:ny),Ib3(0:nx,0:ny),d1,d2,tt1,tt2,teta,L1,L2,x3,y3,dt,sumIbb
integer i,j,nx,ny

sumIbb=0.d0 

xbar1=0.5d0*(x1+x2) ; ybar1=0.5d0*(y1+y2) 


L1=dsqrt ( (x2-x1)**2 + (y2-y1)**2 )
if (y1.ge.y2) then 
teta=dasin( (x1-x2)/L1 )
else 
teta=dasin( (x2-x1)/L1 )
end if 
x3=xbar1- 0.5d0*L1*dcos(teta)
y3=ybar1+ 0.5d0*L1*dsin(teta)
L2=dsqrt ( (x3-xbar1)**2 + (y3-ybar1)**2 ) 

do i=0,nx ; do j=0,ny 

d1= abs( (x2-x1)*(y1-y(j))-(x1-x(i))*(y2-y1) )/L1
!Ib1(i,j)=0.5d0+0.5d0*(rr-d1)/dsqrt((rr-d1)**2+gap*gap)
Ib1(i,j)=0.5d0+0.5d0*( (rr-d1)**3 + 1.5d0* gap*gap *(rr-d1) )  / ( (rr-d1)*(rr-d1)+gap*gap )**1.5d0  

d2= abs( (x3-xbar1)*(ybar1-y(j))-(xbar1-x(i))*(y3-ybar1) )/L2

!Ib2(i,j)=0.5d0+0.5d0*(0.5d0*len1-d2)/dsqrt((0.5d0*len1-d2)**2+gap*gap) 
Ib2(i,j)=0.5d0+0.5d0*( (0.5d0*len1-d2)**3 + 1.5d0* gap*gap *(0.5d0*len1-d2) )  / ( (0.5d0*len1-d2)*(0.5d0*len1-d2)+gap*gap )**1.5d0  
 
Ib1(i,j)=Ib1(i,j)*Ib2(i,j)


if (Ib1(i,j).le.0.001) then 
Ib1(i,j)=0.d0
else if (Ib1(i,j).ge.0.999) then 
Ib1(i,j)=1.0d0 
end if

end do ; end do  
 

return 

End 

subroutine insidef2(x1,y1,Ib1,x,y,rr,gap,nx,ny)   
 
implicit none
double precision x1,y1,Ib1(0:nx,0:ny),x(-1:nx+1),y(-1:ny+1),gap
double precision rr,d1
integer i,j,nx,ny




do i=0,nx ; do j=0,ny 

d1=sqrt( (x(i)-x1)**2 + (y(j)-y1)**2)
!Ib1(i,j)=0.5d0+0.5d0*(rr-d1)/dsqrt((rr-d1)**2+gap*gap)
Ib1(i,j)=0.5d0+0.5d0*( (rr-d1)**3 + 1.5d0* gap*gap *(rr-d1) )  / ( (rr-d1)*(rr-d1)+gap*gap )**1.5d0  

if (Ib1(i,j).le.0.001) then 
Ib1(i,j)=0.d0
else if (Ib1(i,j).ge.0.999) then 
Ib1(i,j)=1.0d0 
end if

end do ; end do  
 

return 

End 

!!!!!!!!!!!!!!!!!! subroutine for boundary condition of  level set function !!!!!!!!!!!!!!!!!!!!
subroutine boundarycondphi()
use M_General_2D,  only:phi,nx,ny 
implicit none 
integer i,j



do j=1,ny-1 
phi(0,j)=2*phi(1,j)-phi(2,j)
phi(-1,j)=3*phi(1,j)-2*phi(2,j)

phi(nx,j)=2*phi(nx-1,j)-phi(nx-2,j)
phi(nx+1,j)=3*phi(nx-1,j)-2*phi(nx-2,j)
end do 

do i=-1,nx+1 
phi(i,0)=2*phi(i,1)-phi(i,2)
phi(i,-1)=3*phi(i,1)-2*phi(i,2)

phi(i,ny)=2*phi(i,ny-1)-phi(i,ny-2)
phi(i,ny+1)=3*phi(i,ny-1)-2*phi(i,ny-2)
end do 

return 
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!!!!!!!!!!!!!!!! finding tether forces ( not in the code now) !!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!! subroutine for adding tower and nacelle mass (not in the code now ) !!!!!!!!!
subroutine addmass( xbar,ybar,x,y,cbalast,nx,ny)
parameter npoint=8
 implicit none 
 double precision xbar,ybar,length,x(-1:nx+1),y(-1:ny+1)
 double precision, dimension (1:npoint)           :: lminx,lminy,lmin
 double precision, dimension (1:nx-1,1:ny-1) :: cbalast      
 integer nx,ny,i,j,w,ww
 
 !! for adding tower and nacelle weight !! 
 cbalast=0 
 lmin=1000 ; lminx=1000  ; lminy=1000 
 
 do i=1,nx-1 ; do j=1,ny-1 

  length= dsqrt( (x(i)-xbar)**2 + (y(j)-ybar)**2  )

    do  w=1,npoint
      if (length.lt.lmin(w) ) then 
      
           if (w.ne.npoint) then 
             do  ww=npoint-1,w,-1
              lminx(ww+1)=lminx(ww) ; lminy(ww+1)=lminy(ww)  ;lmin(ww+1)=lmin(ww)
             end do 
           end if 
       lminx(w)=i ; lminy(w)=j  ;lmin(w)=length
      
        exit 
      
      end if 
    end do 
 
 

 end do ; end do 
 
 
 do w=1,npoint 
  cbalast(lminx(w),lminy(w))=1.0/real(npoint)  
 end do 
 
 return 
 
 end 
 


subroutine contactangle(uext,vext,dt,Ib,gap)
use M_General_2D,  only:phi,x,y,nx,ny 
implicit none
dOUBLE PRECISION gap,dt
double precision,dimension (0:nx+1,-1:ny+1)       ::   uext
double precision,dimension (-1:nx+1,0:ny+1)       ::   vext
double precision,dimension (0:nx,0:ny)            ::   Ib,divext,divextx,divexty,Iphi

integer i,j

do i=0,nx ; do j=0,ny 

Iphi(i,j)=(0.5d0+0.5d0*( phi(i,j)**3 + 1.5d0* gap*gap *phi(i,j) )  / ( phi(i,j)*phi(i,j)+gap*gap )**1.5d0)
 
end do  ; end do 


divext=0 ; divextx=0 ; divexty=0 ; uext=0 ; vext=0 
do i=1,nx-1 ; do j=1,ny-1 
if ( ib(i,j)    .ne.0.0.AND.ib(i+1,j)  .ne.0.0.AND.ib(i-1,j)  .ne.0.0.AND.ib(i,j-1)  .ne.0.0.AND.ib(i,j+1).ne.0.0.AND.& 
&    ib(i+1,j+1).ne.0.0.AND.ib(i+1,j-1).ne.0.0.AND.ib(i-1,j+1).ne.0.0.And.ib(i-1,j-1).ne.0.0.AND. &
&    ib(i,j)    .ne.1.0.AND.ib(i+1,j)  .ne.1.0.AND.ib(i-1,j)  .ne.1.0.AND.ib(i,j-1)  .ne.1.0.AND.ib(i,j+1).ne.1.0.AND.& 
&    ib(i+1,j+1).ne.1.0.AND.ib(i+1,j-1).ne.1.0.AND.ib(i-1,j+1).ne.1.0.And.ib(i-1,j-1).ne.1.0.AND. &
&    Iphi(i,j).gt.0.001.AND. Iphi(i,j).lt.0.999) then
 
!divextx(i,j)=sqrt ( ( (Ib(i,j)-Ib(i-1,j))/( x(i)-x(i-1) ) )**2 + ( (Ib(i,j)-Ib(i,j-1))/( y(j)-y(j-1) ) )**2 )
divextx(i,j)=sqrt ( ( ( 0.25*(Ib(i,j)+Ib(i-1,j)+Ib(i,j+1)+Ib(i-1,j+1))-0.25*(Ib(i,j-1)+Ib(i-1,j-1)+Ib(i,j)+Ib(i-1,j))   ) &
& /( 0.5*(y(j+1)-y(j-1)) ) )**2 + ( ( Ib(i,j)-Ib(i-1,j) )/( x(i)-x(i-1) ) )**2 )


divexty(i,j)=sqrt ( ( ( 0.25*( Ib(i,j)+Ib(i+1,j)+Ib(i,j-1)+Ib(i+1,j-1) )-0.25*(Ib(i,j)+Ib(i-1,j)+Ib(i,j-1)+Ib(i-1,j-1)) )/( 0.5*(x(i+1)-x(i-1)) ) )**2 + &
( (Ib(i,j)-Ib(i,j-1))/(y(j)-y(j-1)) )**2  ) 
 
!if (divextx(i,j).gt.0.0001.AND.divexty(i,j).gt.0.0001) then 
uext(i,j)=8*( (0.5-0.5*(Ib(i,j)+Ib(i-1,j)))**3 )*(Ib(i,j)-Ib(i-1,j))/( x(i)-x(i-1) )/( divextx(i,j)+0.001**2)
vext(i,j)=8*( (0.5-0.5*(Ib(i,j)+Ib(i,j-1)))**3 )*(Ib(i,j)-Ib(i,j-1))/( y(j)-y(j-1) )/( divexty(i,j)+0.001**2) 
else 
uext(i,j)=0
vext(i,j)=0
end if 

end do ; end do 



 return 
 end 
 
 
 subroutine levelset2(dt,tp,ib)
use M_General_2D,  only:phi,u,v,x,y,nx,ny
use M_Mesh_2D,     only: hx,hy    
implicit none 
Double precision phiold(-1:nx+1,-1:ny+1)
double precision, DIMENSION (0:nx,0:ny)       ::Dpphix,Dmphix,Dpphiy,Dmphiy
double precision, DIMENSION (1:nx-1,1:ny-1)   ::phix,phiy,Lphin
Double precision Lphis,dt,dphi,dphimax,avgh,bb,dtau
double precision,dimension (0:nx,0:ny)       ::   Ib
integer i,j,kk,LL,tp


!! since level set and contact angle equatio have the same shape Level set subroutine is used for contact angle equation 
!! but with fictious time steps !! 
avgh=1000 !! large number!! 
do i=-1,nx+1 ; do j=-1, ny+1 
bb=min (hx(i),hy(j))
if (bb.lt.avgh) then
avgh=bb
end if 
end do ; end do 

dtau=0.5*avgh
dtau=dt
do LL=1,4

phiold(-1:nx+1,-1:ny+1)=phi(-1:nx+1,-1:ny+1)


do kk=1,2  !!prediction correction method!!


do i=0,nx ; do j=0,ny                                
Dpphix(i,j)=( phi(i+1,j)-phi(i,j)   )/(  0.5*(hx(i)+hx(i+1))  )
Dmphix(i,j)=( phi(i,j)  -phi(i-1,j) )/(  0.5*(hx(i)+hx(i-1))  )
Dpphiy(i,j)=( phi(i,j+1)-phi(i,j)   )/(  0.5*(hy(j)+hy(j+1))  )
Dmphiy(i,j)=( phi(i,j)  -phi(i,j-1) )/(  0.5*(hy(j)+hy(j-1))  )

end do ;end do 

do i=1,nx-1 ; do j=1,ny-1    !!A!!

if (0.5*( u(i,j)+u(i+1,j) ).gt.0.d0) then


  if (  abs(  Dmphix(i,j)-Dmphix(i-1,j) ).lt.abs(  Dpphix(i,j)-Dpphix(i-1,j) )   ) then
phix(i,j)=Dmphix(i,j)+  0.5*(  Dmphix(i,j)-Dmphix(i-1,j)  ) 
   else
phix(i,j)=Dmphix(i,j)+  0.5*(  Dpphix(i,j)-Dpphix(i-1,j)  )
   end if 
    
else

   if (  abs(  Dmphix(i+1,j)-Dmphix(i,j) ).lt.abs(  Dpphix(i+1,j)-Dpphix(i,j) )   ) then
phix(i,j)=Dpphix(i,j)-  0.5*(  Dmphix(i+1,j)-Dmphix(i,j)  ) 
   else
phix(i,j)=Dpphix(i,j)-  0.5*(  Dpphix(i+1,j)-Dpphix(i,j)  )
   end if 

end if 



if (0.5*( V(i,j)+V(i,j+1) ).gt.0.d0) then


   if (  abs(  DmphiY(i,j)-DmphiY(i,j-1) ).lt.abs(  DpphiY(i,j)-DpphiY(i,j-1) )   ) then
phiY(i,j)=DmphiY(i,j)+  0.5*(  DmphiY(i,j)-DmphiY(i,j-1)  ) 
   else
phiY(i,j)=DmphiY(i,j)+  0.5*(  DpphiY(i,j)-DpphiY(i,j-1)  )
   end if 
    
else

  if (  abs(  DmphiY(i,j+1)-DmphiY(i,j) ).lt.abs(  DpphiY(i,j+1)-DpphiY(i,j) )   ) then
phiY(i,j)=DpphiY(i,j)-  0.5*(  DmphiY(i,j+1)-DmphiY(i,j)  ) 
   else
phiY(i,j)=DpphiY(i,j)-  0.5*(  DpphiY(i,j+1)-DpphiY(i,j)  )
   end if 


end if 


     
end do ;end do  !!A!!

if (kk.eq.1) then

do i=1,nx-1 ; do j=1,ny-1 
Lphin(i,j)=  -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
         &   -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
           
phi(i,j)=phi(i,j)+dtau*( -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)+&
                     & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) )
                
end do ;end do 

                  
end if                        

   if (kk.eq.2) then
   
   do i=1,nx-1 ; do j=1,ny-1 
 
Lphis=   -0.5*( u(i,j)+u(i+1,j) )* phix(i,j)&
       & -0.5*( V(i,j)+V(i,j+1) )* phiy(i,j) 
            
phi(i,j)=phiold(i,j)+0.5*dtau*(Lphis+Lphin(i,j) )                       
   end do ;end do 

   end if 


end do  !!prediction corection method!! 
 
!call boundarycondphi()
!
!write(125,*) 'zone i=',nx-2,' j=',ny-2
!do j=1,ny-2 ;  Do i=1,nx-2
!
!write(125,560) x(i),y(j),0.5*(u(i,j)+u(i+1,j)),0.5*(v(i,j)+v(i,j+1)),phi(i,j),Ib(i,j)
!end do ; end do
!560 format (6(1x,e15.7)) 


end do 




                   

RETURN
END 



subroutine property(Ib,Ic,romiu,romiudrop,romiuair,romiusolid,romiucon,gap)

use M_General_2D,  only: phi,x,y,nx,ny
implicit none 

integer i,j 
dOUBLE PRECISION gap,romiudrop,romiusolid,romiuair,romiucon
double precision,dimension (0:nx,0:ny)       ::   Ib,Ic,Iphi,romiu


!! this subroutime is used for assiginig density and viscosity for the whole domain !! 



do i=0,nx ; do j=0,ny 

Iphi(i,j)=(0.5d0+0.5d0*( phi(i,j)**3 + 1.5d0* gap*gap *phi(i,j) )  / ( phi(i,j)*phi(i,j)+gap*gap )**1.5d0)

!if (Iphi(i,j).lt.0.001) then 
!romiu(i,j)=romiudrop
!else if (Iphi(i,j).gt.0.999) then
!romiu(i,j)=romiuair
!else  
romiu(i,j)=romiudrop !+Iphi(i,j)*(romiuair-romiudrop)
  
!end if 

!romiu(i,j)=romiu(i,j)+Ib(i,j)*( romiusolid-romiu(i,j) )
!romiu(i,j)=romiu(i,j)+Ic(i,j)*( romiucon )
end do ; end do 




return 

end   

   
 subroutine advection (advectu,advectv)  

use M_General_2D,  only:u,v,nx,ny
use M_Mesh_2D,     only:hx,hy  
 implicit none 
 double precision ux,uy,vx,vy
 double precision,dimension (1:nx,0:ny)       ::   dpux,dmux,dpuy,dmuy,advectu
 double precision,dimension (0:nx,1:ny)        ::   dpvx,dmvx,dpvy,dmvy,advectv

 
 integer i,j,k
 !!!!!!!!!!!! second order ENO method for advection terms !!!!!!!!!!!!
 
 do i=1,nx ; do j=0,ny 
 Dpux(i,j)=( u(i+1,j)-u(i,j)   )/hx(i)
 Dmux(i,j)=( u(i,j)  -u(i-1,j) )/hx(i-1)
 Dpuy(i,j)=( u(i,j+1)-u(i,j)   )/( 0.5d0*( hy(j+1)+hy(j)) )
 Dmuy(i,j)=( u(i,j)  -u(i,j-1) )/( 0.5d0*( hy(j)+hy(j-1)) )
 
 end do ;end do 

 do i=2,nx-1 ; do j=1,ny-1    !!A!!
 
 
 !!U1    
 if (u(i,j).gt.0.d0) then

    if (  abs(  Dmux(i,j)-Dmux(i-1,j) ).lt.abs(  Dpux(i,j)-Dpux(i-1,j) )   ) then
 ux=Dmux(i,j)+  0.5*(  Dmux(i,j)-Dmux(i-1,j)  ) 
    else
 ux=Dmux(i,j)+  0.5*(  Dpux(i,j)-Dpux(i-1,j)  )
    end if 
    
 else

    if (  abs(  Dmux(i+1,j)-Dmux(i,j) ).lt.abs(  Dpux(i+1,j)-Dpux(i,j) )   ) then
 ux=Dpux(i,j)-  0.5*(  Dmux(i+1,j)-Dmux(i,j)  ) 
    else
 ux=Dpux(i,j)-  0.5*(  Dpux(i+1,j)-Dpux(i,j)  )
    end if 

 end if 


!!U2
 if (  (0.25*( V(i,j)+V(i,j+1)+v(i-1,j)+v(i-1,j+1) )).gt.0.d0) then

    if (  abs(  DmuY(i,j)-DmuY(i,j-1) ).lt.abs(  DpuY(i,j)-DpuY(i,j-1) )   ) then
 uY=DmuY(i,j)+  0.5*(  DmuY(i,j)-DmuY(i,j-1)  ) 
    else
 uY=DmuY(i,j)+  0.5*(  DpuY(i,j)-DpuY(i,j-1)  )
    end if 
    
 else

    if (  abs(  DmuY(i,j+1)-DmuY(i,j) ).lt.abs(  DpuY(i,j+1)-DpuY(i,j) )   ) then
 uY=DpuY(i,j)-  0.5*(  DmuY(i,j+1)-DmuY(i,j)  ) 
    else
 uY=DpuY(i,j)-  0.5*(  DpuY(i,j+1)-DpuY(i,j)  )
    end if 

 end if 
 
 
 
 advectu(i,j)=-(  u(i,j)*uX+0.25*( V(i,j)+V(i,j+1)+v(i-1,j)+v(i-1,j+1) )*uY  )

 end do ; end do
  
 do i=0,nx ; do j=1,ny 
 Dpvx(i,j)=( v(i+1,j)-v(i,j)   )/( 0.5d0*( hx(i+1)+hx(i)) )
 Dmvx(i,j)=( v(i,j)  -v(i-1,j) )/( 0.5d0*( hx(i)+hx(i-1)) )
 Dpvy(i,j)=( v(i,j+1)-v(i,j)   )/hy(j)
 Dmvy(i,j)=( v(i,j)  -v(i,j-1) )/hy(j-1)
 
 end do ;end do 

 do i=1,nx-1 ; do j=2,ny-1  !!A!!
 
 
!!V1
 if (  (0.25*( u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1) )).gt.0.d0) then

    if (  abs(  Dmvx(i,j)-Dmvx(i-1,j) ).lt.abs(  Dpvx(i,j)-Dpvx(i-1,j) )   ) then
 vx=Dmvx(i,j)+  0.5*(  Dmvx(i,j)-Dmvx(i-1,j)  ) 
    else
 vx=Dmvx(i,j)+  0.5*(  Dpvx(i,j)-Dpvx(i-1,j)  )
    end if 
    
 else

    if (  abs(  Dmvx(i+1,j)-Dmvx(i,j) ).lt.abs(  Dpvx(i+1,j)-Dpvx(i,j) )   ) then
 vx=Dpvx(i,j)-  0.5*(  Dmvx(i+1,j)-Dmvx(i,j)  ) 
    else
 vx=Dpvx(i,j)-  0.5*(  Dpvx(i+1,j)-Dpvx(i,j)  )
    end if 

 end if 



!!V2
 if ( V(i,j).gt.0.d0) then

    if (  abs(  DmvY(i,j)-DmvY(i,j-1) ).lt.abs(  DpvY(i,j)-DpvY(i,j-1) )   ) then
 vY=DmvY(i,j)+  0.5*(  DmvY(i,j)-DmvY(i,j-1)  ) 
    else
 vY=DmvY(i,j)+  0.5*(  DpvY(i,j)-DpvY(i,j-1)  )
    end if 
    
 else

    if (  abs(  DmvY(i,j+1)-DmvY(i,j) ).lt.abs(  DpvY(i,j+1)-DpvY(i,j) )   ) then
 vY=DpvY(i,j)-  0.5*(  DmvY(i,j+1)-DmvY(i,j)  ) 
    else
 vY=DpvY(i,j)-  0.5*(  DpvY(i,j+1)-DpvY(i,j)  )
    end if 

end if 
 
 
 advectv(i,j)=-(0.25*(  u(i,j)+u(i+1,j)+u(i,j-1)+u(i+1,j-1)  )*vX+ V(i,j)*vY)

 end do ; end do 
 
 
 return 
 
  end 


  



Subroutine viscosity (ro,miuv,Tx,Ty)
use M_General_2D,  only:u,v,nx,ny 
use M_Mesh_2D,     only:hx,hy  
implicit none 
dOUBLE PRECISION  Txxr,Txxl,Tyxd,Tyxu,Tyyu,Tyyd,Txyr,Txyl
double precision,dimension (0:nx,0:ny)            ::   ro,miuv
double precision  Tx(2:nx-1,1:ny-1),Ty(1:nx-1,2:ny-1)
integer i,j

!!!!!!!!!!!! simple centeral difference for viscose terms !!!!!!!!!!!!!!!!!!

Do i=2,nx-1 ; Do j=1, ny-1 

 Txxr=2*miuv(i,j)  *( U(i+1,j)-U(i,j)   )/hx(i)
 
 Txxl=2*miuv(i-1,j)*( U(i,j)  -U(i-1,j) )/hx(i-1)

 Tyxu=0.25*( miuv(i,j)+miuv(i,j+1)+miuv(i-1,j)+miuv(i-1,j+1) )*&
   &(      ( U(i,j+1)-U(i,j)     )/( 0.5*(hy(j)+hy(j+1)) )+ ( V(i,j+1)-V(i-1,j+1) )/( 0.5*(hx(i)+hx(i-1)) )     )
                 
 Tyxd=0.25*( miuv(i,j-1)+miuv(i,j)+miuv(i-1,j-1)+miuv(i-1,j) )*&
   &(      ( U(i,j)-U(i,j-1)     )/( 0.5*(hy(j)+hy(j-1)) )+ ( V(i,j)  -V(i-1,j)   )/( 0.5*(hx(i)+hx(i-1)) )     )

  

 Tx(i,j)= (  (Txxr-Txxl)/( 0.5*( hx(i)+hx(i-1)) ) + (Tyxu-Tyxd)/hy(j)  )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 
   
end do ; end do   




Do j=2,ny-1 ; Do i=1,nx-1 


Tyyu=2*miuv(i,j)  *( V(i,j+1)-V(i,j)   )/hy(j)
 
 Tyyd=2*miuv(i,j-1)*( V(i,j)  -V(i,j-1) )/hy(j-1)
 
 
 Txyr=0.25*( miuv(i,j)+miuv(i,j-1)+miuv(i+1,j)+miuv(i+1,j-1) )*&
  & (     ( V(i+1,j)-V(i,j)     )/( 0.5*(hx(i)+hx(i+1)) )+ ( U(i+1,j)  -U(i+1,j-1) )/( 0.5*(hy(j)+hy(j-1))  )   )
  
 Txyl=0.25*( miuv(i-1,j)+miuv(i-1,j-1)+miuv(i,j)+miuv(i,j-1) )*&
  & (     ( V(i,j)-V(i-1,j)     )/( 0.5*(hx(i)+hx(i-1)) )+ ( U(i,j)    -U(i,j-1)   )/( 0.5*(hy(j)+hy(j-1))  )   )
  
 

 Ty(i,j)= (  (Tyyu-Tyyd)/( 0.5*( hy(j)+hy(j-1)) ) + (Txyr-Txyl)/hx(i) )/(  0.5*(ro(i,j)+ro(i,j-1))  )

 
end do ; end do 
 

return 
end 


!!!!!!!!!!!!!! simple iterative method for solving poisson equation !!!!!!!!!!!!
subroutine poisson (ro,dt,pdif,p,beta)
use M_General_2D,  only:u,v,nx,ny 
use M_Mesh_2D,     only:hx,hy
implicit none 

double precision dt,beta,maxp,pdif
double precision,dimension (0:nx,0:ny)            ::   ro,miuv,p,pold
double precision,dimension (1:nx-1,1:ny-1)        ::   Apx,Amx,Apy,Amy,AP,Q
integer i,j,number2
 


do i=1,nx-1 ; do j=1,ny-1 
                         
 
 if (i==nx-1) then
 Apx(i,j)=0
 else
 
 Apx(i,j)=1/( 0.5d0*(hx(i)+hx(i+1))*hx(i) )/(  0.5*(ro(i,j)+ro(i+1,j))  )
 end if 
 
 if (i==1)then 
 Amx(i,j)=0       
 else 
 
 Amx(i,j)=1/( 0.5d0*(hx(i)+hx(i-1))*hx(i) )/(  0.5*(ro(i,j)+ro(i-1,j))  )
 end if 
 
 if (j==ny-1) then
 Apy(i,j)=0
 else
 
 Apy(i,j)=1/( 0.5d0*(hy(j)+hy(j+1))*hy(j) )/(  0.5*(ro(i,j)+ro(i,j+1))  )

 end if 
 if (j==1) then
 Amy(i,j)=0
 else
 
 Amy(i,j)=1/( 0.5d0*(hy(j)+hy(j-1))*hy(j) )/(  0.5*(ro(i,j)+ro(i,j-1))  )
 
 end if 
 
 
 
 AP(i,j)=-( Apx(i,j)+Amx(i,j)+Apy(i,j)+Amy(i,j) )   
 
 Q(I,J)=(  (U(I+1,J)-U(I,J))/hx(i)+(V(I,J+1)-V(I,J))/hy(j) )/dt
 

 end do ; end do 
 
 
 number2=0
 maxp=1.5
           do while (maxp.gt.pdif.AND.number2.lt.10000)

 maxp=0 
 number2=number2+1
 
 do i=1,nx-1 ; do j=1,ny-1 

 pold(i,j)=p(i,j)
 p(i,j)=beta*( ( Apx(i,j)*P(I+1,j)+Amx(i,j)*P(I-1,j)+Apy(i,j)*P(I,j+1)+Amy(i,j)*P(I,j-1)-Q(I,j)  )/(-AP(i,j))  ) &
         & +(1-beta)*p(i,j)
 if (abs( pold(i,j)-p(i,j) ).gt.maxp) then 
 maxp=abs( pold(i,j)-p(i,j) )
 end if 
 
 end do ; end do
 
            end do  !!while !
                        
 
   print *,"pdif=",maxp,"number it=",number2
   
   
return 
end     
 
 
 
 



Subroutine extforce(px,py,ox,oy,sumextu,sumextv,sumextIzb,lzero,anac,acgx,acgy,xbar,ybar,tp,dt)
implicit none 
parameter tetnum=2
integer tp
dOUBLE PRECISION, intent(in)  ::  Px(1:12,1:4),Py(1:12,1:4),Ox(1:3),Oy(1:3),Lzero(1:3)
dOUBLE PRECISION, intent(in)  ::  acgx,acgy,anac,xbar,ybar,dt
double precision, intent(out) ::  sumextu,sumextv,sumextIzb
   
!!local variables !!
double precision fmoorx(1:tetnum),fmoory(1:tetnum),ks,lnew(1:tetnum),towerh,mtower,rx,ry,fbx,fby,Momen,gy,fthrust
integer i
!include 'fwtpa.txt'



!acgx=0 
!acgy=0 
!anac=0 
 !fthrust=tp*(5*dt)*fthrust        
do i=1,tetnum   
   Lnew(i)=sqrt(  (ox(i)-px(i,4))**2 + (oy(i)-py(i,4))**2  )
   Fmoorx(i)=Ks* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Ox(i)-Px(i,4) )/Lnew(i)
   Fmoory(i)=Ks* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Oy(i)-Py(i,4) )/Lnew(i)
end do 


!! towerh should be distance from CG!!

     rx=towerh*( px(1,1)-px(2,2) )/sqrt( (px(1,1)-px(2,2))**2+(py(1,1)-py(2,2))**2)
     ry=towerh*( py(1,1)-py(2,2) )/sqrt( (px(1,1)-px(2,2))**2+(py(1,1)-py(2,2))**2)
     
     fbx=fthrust-mtower*(acgx-ry*anac)
     fby=mtower*gy-mtower*(acgy+rx*anac)
     
     Momen=-Mtower*towerh**2*anac+mtower*gy*rx-fthrust*ry
     
     sumextu=(Fmoorx(1)+Fmoorx(2)+fbx)*dt
     sumextv=(Fmoory(1)+Fmoory(2)+fby)*dt
     sumextIzb=(Momen+ &
              & (px(1,4)-xbar)*fmoory(1)-(py(1,4)-ybar)*fmoorx(1)+ &
              & (px(2,4)-xbar)*fmoory(2)-(py(2,4)-ybar)*fmoorx(2)  )*dt

!      sumextu=(Fmoorx(1)+Fmoorx(2))*dt
!      sumextv=(Fmoory(1)+Fmoory(2))*dt
!      sumextIzb=( &
!              & (px(1,4)-xbar)*fmoory(1)-(py(1,4)-ybar)*fmoorx(1)+ &
!              & (px(2,4)-xbar)*fmoory(2)-(py(2,4)-ybar)*fmoorx(2)  )*dt
              
 return 
 end 
 
subroutine tether(px,py,x,y,nx,ny,Ox,Oy,ks,Fmoorx,Fmoory,Lzero,Ic,Lnew,xbar,ybar,ubar,xbar0)


implicit none 
parameter npoint=16,tetnum=1
dOUBLE PRECISION x(-1:nx+1),y(-1:ny+1),Px(1:12,1:4),Py(1:12,1:4),Ks 
dOUBLE PRECISION length,lmin(1:npoint,1:3),xbar,ybar,ubar,xbar0,cvis
double precision,dimension (0:nx,0:ny)       :: Ic
double precision,dimension (1:3)             :: Lzero,lnew,ox,oy
double precision,dimension (1:nx-1,1:ny-1)   :: fmoorx,fmoory
integer i,j,nx,ny
integer w,ww,con,lminx(1:npoint,1:3),lminy(1:npoint,1:3),tet


cvis=19.8
fmoorx=0 ; fmoory=0 
lmin=1000 ; lminx=1000  ; lminy=1000  

 

do i=1,nx-1 ; do j=1,ny-1 

 if (Ic(i,j).gt.0.9) then
  do tet=1,tetnum          !! number of tethers !!
 ! length= dsqrt( (x(i)-px(tet,4))**2 + (y(j)-py(tet,4))**2  )
 length= dsqrt( (x(i)-xbar)**2 + (y(j)-ybar)**2  )
 

    do  w=1,npoint 
      if (length.lt.lmin(w,tet) ) then 
      
       con=w
           if (con.ne.npoint) then 
             do  ww=npoint-1,con,-1
              lminx(ww+1,tet)=lminx(ww,tet) ; lminy(ww+1,tet)=lminy(ww,tet)  ;lmin(ww+1,tet)=lmin(ww,tet)
             end do 
           end if 
       lminx(w,tet)=i ; lminy(w,tet)=j  ;lmin(w,tet)=length
      
        exit 
      
      end if 
    end do 
  end do 
 end if 

end do ; end do 


  do i=1,tetnum   !! number of tethers !!
   Lnew(i)=sqrt(  (ox(i)-px(i,4))**2 + (oy(i)-py(i,4))**2  )
  end do 


  do i=1,tetnum    !! number of tethers !!
    do w=1,npoint 
     Fmoorx(lminx(w,i),lminy(w,i))=1.0/real(npoint)*Ks* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Ox(i)-Px(i,4) )/Lnew(i)
     Fmoory(lminx(w,i),lminy(w,i))=1.0/real(npoint)*Ks* ( max(Lzero(i),Lnew(i))-Lzero(i) )* ( Oy(i)-Py(i,4) )/Lnew(i)
    end do 
  end do 

!do i=1,tetnum    
!   do w=1,npoint 
!     Fmoorx(lminx(w,i),lminy(w,i))=1.0/real(npoint)*(-Ks*(xbar-xbar0)-cvis*ubar )
!     
!    end do 
!end do 


     
 
return 
end 
 
 
 subroutine impactV(fximp,fyimp,u,v,ubar,vbar,xbar,ybar,omegaz,x,y,nx,ny,Ib,tp,numimp,r2,dt)

implicit none
dOUBLE PRECISION ubar,vbar,xbar,ybar,omegaz,x(-1:nx+1),y(-1:ny+1),r2,a,b,sum,dt
double precision,dimension (0:nx+1,-1:ny+1)       ::   u,fximp
double precision,dimension (-1:nx+1,0:ny+1)       ::   v,fyimp
double precision,dimension (0:nx,0:ny)            ::   Ib,divext,divextx,divexty

integer nx,ny,i,j,tp,numimp


!divext=0 ; divextx=0 ; divexty=0 
!do i=1,nx-1 ; do j=1,ny-1 
!if ( ib(i,j)    .ne.0.0.AND.ib(i+1,j)  .ne.0.0.AND.ib(i-1,j)  .ne.0.0.AND.ib(i,j-1)  .ne.0.0.AND.ib(i,j+1).ne.0.0.AND.& 
!&    ib(i+1,j+1).ne.0.0.AND.ib(i+1,j-1).ne.0.0.AND.ib(i-1,j+1).ne.0.0.And.ib(i-1,j-1).ne.0.0.AND. &
!&    ib(i,j)    .ne.1.0.AND.ib(i+1,j)  .ne.1.0.AND.ib(i-1,j)  .ne.1.0.AND.ib(i,j-1)  .ne.1.0.AND.ib(i,j+1).ne.1.0.AND.& 
!&    ib(i+1,j+1).ne.1.0.AND.ib(i+1,j-1).ne.1.0.AND.ib(i-1,j+1).ne.1.0.And.ib(i-1,j-1).ne.1.0 ) then
! 
!
!divextx(i,j)=sqrt ( ( ( 0.25*(Ib(i,j)+Ib(i-1,j)+Ib(i,j+1)+Ib(i-1,j+1))-0.25*(Ib(i,j-1)+Ib(i-1,j-1)+Ib(i,j)+Ib(i-1,j))   ) &
!& /( 0.5*(y(j+1)-y(j-1)) ) )**2 + ( ( Ib(i,j)-Ib(i-1,j) )/( x(i)-x(i-1) ) )**2 )
!
!
!divexty(i,j)=sqrt ( ( ( 0.25*( Ib(i,j)+Ib(i+1,j)+Ib(i,j-1)+Ib(i+1,j-1) )-0.25*(Ib(i,j)+Ib(i-1,j)+Ib(i,j-1)+Ib(i-1,j-1)) )/( 0.5*(x(i+1)-x(i-1)) ) )**2 + &
!( (Ib(i,j)-Ib(i,j-1))/(y(j)-y(j-1)) )**2  ) 
! 
!
!uimp(i,j)= 50*( 1.0-0.5*(Ib(i,j)+Ib(i+1,j)) )*( ( ubar  -omegaz*(y(j)-ybar) )-u(i,j) )* &
!&          abs ( (Ib(i,j)-Ib(i-1,j))/( x(i)-x(i-1) ) )/( divextx(i,j)+0.001**2)+uimp(i,j)
!vimp(i,j)= 50*( 1.0-0.5*(Ib(i,j)+Ib(i,j+1)) )*( ( vbar  +omegaz*(x(i)-xbar) )-v(i,j) )* &
!&          abs ( (Ib(i,j)-Ib(i,j-1))/( y(j)-y(j-1) ) )/( divexty(i,j)+0.001**2)+ vimp(i,j)
!else 
!uimp(i,j)=uimp(i,j)
!vimp(i,j)=vimp(i,j)
!end if 
!
!end do ; end do 
!-(ubar-u(i,j))


if (numimp.eq.1) then 

do i=1,nx-1 ; do j=1,ny-1 

fximp(i,j)=( 1.0-ubar*dt/(x(i+1)-x(i)) )*fximp(i,j)+ ubar*dt/(x(i+1)-x(i))*fximp(i-1,j)
fyimp(i,j)=( 1.0-ubar*dt/(x(i+1)-x(i)) )*fyimp(i,j)+ ubar*dt/(x(i+1)-x(i))*fyimp(i-1,j)

end do ; end do 

end if 

a=100
b=0.03

do i=1,nx-1 ; do j=1,ny-1 
if ( ( sqrt((x(i)-xbar)**2+(y(j)-ybar)**2)-r2 )**2.lt. 0.05**2 ) then
fximp(i,j)=fximp(i,j)-(u(i,j)-ubar)*a*exp( ((sqrt((x(i)-xbar)**2+(y(j)-ybar)**2) -r2)**2 ) /-b**2)*dt
fyimp(i,j)=fyimp(i,j)-(v(i,j)-vbar)*a*exp( ((sqrt((x(i)-xbar)**2+(y(j)-ybar)**2) -r2)**2 ) /-b**2)*dt
else 
fximp(i,j)=0
fyimp(i,j)=0
end if
 
end do ; end do 



 sum=0
 do i=1,nx-1 ; do j=1,ny-1 
 
 sum=sum+sqrt(fximp(i,j)**2 + fyimp(i,j)**2) 
 
 end do ; end do 
 
 print*, "impact=",sum,numimp 
 
 
 write(145,*) numimp,sum  

if ( mod(tp,50).eq.0.AND.numimp.eq.1) then 

write(135,*) 'zone i=',nx-1,' j=',ny-1
  Do j=1,ny-1 ; Do i=1,nx-1
    write(135,*) x(i),y(j),0.5*(fximp(i,J)+fximp(i+1,j)),0.5*(fyimp(i,j)+fyimp(i,j+1)),Ib(i,j)
end do ; end do 


end if 





return 
end 
               
     
     






 
  















