module M_General_2D

     implicit none 
     save
integer, parameter :: nx=501, ny=126  !nx=143  ,ny=81   
   
double precision, parameter  :: pi=3.141592653, gy=0.0, froude=1.0, landa=2.07  
double precision,dimension (0:nx+1,-1:ny+1)       ::    u
double precision,dimension (-1:nx+1,0:ny+1)       ::    v
double precision,dimension (-1:nx+1,-1:ny+1)      ::   phi
real(8) x(-1:nx+1),y(-1:ny+1)
real(8) lx,ly

       




end module
