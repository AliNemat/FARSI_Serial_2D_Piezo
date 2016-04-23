module M_Mesh_2D
use M_General_2D,           only: x,y,nx,ny,lx,ly,landa
use M_Platform_Constant_2D, only: len,floatx,towerb,r2
dOUBLE PRECISION        :: hx(-1:nx+1),hy(-1:ny+1)






contains 
subroutine meshgenerator()

integer          :: meshgen
dOUBLE PRECISION :: Lx1,Lx2,Lx3,Lx4,Ly1,Ly2,Ly3
dOUBLE PRECISION :: hxx1,hxx2,hxx3,hxx4,hyy1,hyy2,hyy3
dOUBLE PRECISION :: Amesh12,Amesh14,Amesh21,Amesh23,betam12,betam14,betam21,betam23
dOUBLE PRECISION :: Dmesh12,Dmesh14,Dmesh21,Dmesh23
integer          :: nx1,nx2,nx3,ny1,ny2

include 'Par_Mesh_2D.txt'
if (meshgen.eq.1) then 
    x(1)=0.5*Lx/(nx-1)
    x(-1)=-3*x(1) ; x(0)=-x(1)
    do i=1,nx+1
    x(i)=x(i-1)+Lx/(nx-1)
    end do 


    y(1)=0.5*Ly/(ny-1)
    y(-1)=-3*y(1) ; y(0)=-y(1)
    do j=1,ny+1
    y(j)=y(j-1)+Ly/(ny-1)
    end do 




else

print*,"domain size",lx,ly

nx2=nx2+nx1
nx3=nx3+nx2
 
  
ny2=ny2+ny1


    
    hxx1=1.0/dble(nx1-1) ; hxx2=1.0/dble(nx2-nx1) ;  hxx3=1.0/dble(nx3-nx2) ; hxx4=1.0/dble(nx-nx3)
    hyy1=1.0/dble(ny1-1) ; hyy2=1.0/dble(ny2-ny1) ;  hyy3=1.0/dble(ny-ny2)  
    

    Amesh12=1/(2*betam12)*log(  (  1+( exp(betam12)-1 )*Dmesh12/Lx2 )/(  1+( exp(-betam12)-1 )*Dmesh12 /Lx2   )    ) 
    Amesh14=1/(2*betam14)*log(  (  1+( exp(betam14)-1 )*Dmesh14/Lx4 )/(  1+( exp(-betam14)-1 )*Dmesh14 /Lx4   )    ) 


    Amesh21=1/(2*betam21)*log(  (  1+( exp(betam21)-1 )*Dmesh21/Ly1 )/(  1+( exp(-betam21)-1 )*Dmesh21 /Ly1   )    ) 
    Amesh23=1/(2*betam23)*log(  (  1+( exp(betam23)-1 )*Dmesh23/Ly3 )/(  1+( exp(-betam23)-1 )*Dmesh23 /Ly3   )    )

  
     do i=1,nx1-1
     x(i)=hxx1*Lx1*(dble(i)-0.5d0)
     end do 
     x(-1)=-3*x(1) ; x(0)=-x(1)
     
     do i=nx1,nx2-1  !nx-1
     kk=i-(nx1-1)
     x(i)=x(nx1-1)+Dmesh12* (  1+ (  sinh ( betam12*(dble(kk)*hxx2-Amesh12) )  )/sinh(betam12*Amesh12)  )
     end do 
     
     do i=nx2,nx3-1
     x(i)=x(nx2-1)+hxx3*Lx3*dble( i-(nx2-1) )
     end do 

     do i=nx3,nx-1
     kk=i-(nx3-1)
     x(i)=x(nx3-1)+Dmesh14* (  1+ (  sinh ( betam14*(dble(kk)*hxx4-Amesh14) )  )/sinh(betam14*Amesh14)  )
     end do 
     x(nx)=x(nx-1)+ ( x(nx-1)-x(nx-2) ) ; x(nx+1)=x(nx)+ ( x(nx)-x(nx-1) )
     
     
     

     
     
     do j=1,ny1-1 
     y(j)=Dmesh21* (  1+ (  sinh ( betam21*(dble(j-0.5)*hyy1-Amesh21) )  )/sinh(betam21*Amesh21)  )
     end do 
     y(-1)=-3*y(1) ;y(0)=-y(1)   
     
     
     do j=ny1,ny2-1
     y(j)=y(ny1-1)+hyy2*Ly2*dble(j-(ny1-1))
     end do 
     
     
     do j=ny2,ny-1
     
     kk=j-(ny2-1)
     y(j)=y(ny2-1)+Dmesh23* (  1+ (  sinh ( betam23*(dble(kk)*hyy3-Amesh23) )  )/sinh(betam23*Amesh23)  )
     end do 
     y(ny)=y(ny-1)+ ( y(ny-1)-y(ny-2) ) ; y(ny+1)=y(ny)+ ( y(ny)-y(ny-1) )
     

     
     

 end if 

hx(0)=x(1)-x(0)                           !! we assume that we use uniform grid for boundary grids!!
hx(-1)=hx(0)                        
do i=1,nx+1                          
hx(i)=2*( x(i)-x(i-1) )-hx(i-1)
end do

hy(0)=y(1)-y(0)
hy(-1)=hy(0)                   
do j=1,ny+1                       
hy(j)=2*( y(j)-y(j-1) )-hy(j-1)
end do




OPEN(205,file='meshx.plt')
OPEN(215,file='meshy.plt')
OPEN(235,file='grideval2DXY.plt')

do i=1,nx-1
write(205,*) i,x(i)
end do 
call flush (205)
do j=1,ny-1
write(215,*) j,y(j)
end do 
call flush (215)


write(235,*) 'zone i=',nx-1,' k=',ny-1
do j=1,ny-1 ;  Do i=1,nx-1

write(235,350) x(i),y(j)
end do ; end do
350 format (2(1x,e15.7))
call flush (235)


 return 
  end subroutine 
  
  


     
end module 


