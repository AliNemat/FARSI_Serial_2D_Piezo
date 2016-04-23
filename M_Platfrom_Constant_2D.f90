module M_Platform_Constant_2D
 use M_General_2D, only: lx,ly,froude
 implicit none
double precision :: r2,len,TowerB,floatx
double precision :: miusolid,miucon,rosolid,rocon,percent

contains 
  subroutine Platform_Constant_Ini()
  implicit none 
    Include "Par_Platform_2D.txt"
    
 
  return 
  end subroutine 
    
end module

  
    
