module interpolation

use generalMOD, only:dp
implicit none

contains

  integer function locate(xx,x)
  implicit none
  ! Binary search algorithm
  ! Given an array xx(1:N), and given a value x, returns a value j such that x is 
  ! between xx(j) and xx(j+1). xx must be monotonic, either increasing or decreasing. 
  ! j = 0 or j = N is returned to indicate that x is out of range.
  ! Reference: Press W. H. et al. 1996, Numerical Recipes in Fortran 90, 2nd ed.,
  !            Cambridge University Press, pp. 1045
  real(dp),dimension(:),intent(in) :: xx
  real(dp),intent(in) :: x

  integer :: n,jl,jm,ju
  logical :: ascnd

  n = size(xx)  
  ascnd = (xx(n) >= xx(1))  ! True if ascending order of table, false otherwise.
  jl = 0                    ! Initialize lower
  ju = n + 1                ! and upper limits.
  
  do
     if (ju-jl <= 1) exit    ! Repeat until this condition is satisfied.
     jm = (ju+jl) / 2        ! Compute midpoint    
     if (ascnd .eqv. (x >= xx(jm))) then
        jl = jm              ! and replace either the lower limit
     else
        ju = jm              ! or the upper limit, as appropriate.
     end if
  end do
  
  if (x == xx(1)) then      ! Then set the output, being careful with the endpoints.
     locate = 1
  else if(x == xx(n)) then
     locate = n - 1
  else
     locate = jl
  end if

  end function locate
  
  ! * * * * * * * * *
    
  real(dp) function lin_interp(xx,f,x)
      implicit none
      ! This function uses linear interpolation to estimate the value
      ! of a function f at point x.
      ! f is assumed to be sampled on a regular grid, with the grid x values specified
      ! by the array xx.
      ! Reference: https://en.wikipedia.org/wiki/Linear_interpolation      
      real(dp),dimension(:),intent(in) :: xx, f
      real(dp),intent(in) :: x

      real(dp) :: denom, x1, x2
      integer :: nx,i

      nx = size(xx)
      
      i = locate(xx,x)
      if(i == 0 .or. i == nx) then
         write(*,*) 'lin_interp: selected array index out of range'
         print*,"xx", xx(1),xx(nx), "x",  x
         stop
      end if

      x1 = xx(i)
      x2 = xx(i+1)
            
      denom = (x2 - x1)
      
      lin_interp = ( f(i)*(x2-x) + f(i+1)*(x-x1) ) / denom

  end function lin_interp
  real(dp) function lin_interp_bound(xx,f,x)
      implicit none
      ! This function uses linear interpolation to estimate the value
      ! of a function f at point x.
      ! f is assumed to be sampled on a regular grid, with the grid x values specified
      ! by the array xx.
      ! Reference: https://en.wikipedia.org/wiki/Linear_interpolation      
      real(dp),dimension(:),intent(in) :: xx, f
      real(dp),intent(in) :: x

      real(dp) :: denom, x1, x2
      integer :: nx,i

      nx = size(xx)
      
      i = locate(xx,x)
      
   
      if(i == 0) then 
         i=1
         lin_interp_bound=f(i)
      elseif ( i == nx) then
         i=nx
         lin_interp_bound=f(i)
      else
         x1 = xx(i)
         x2 = xx(i+1)
               
         denom = (x2 - x1)
         lin_interp_bound = ( f(i)*(x2-x) + f(i+1)*(x-x1) ) / denom
      endif
  end function lin_interp_bound
  ! * * * * * * * * *

  real(dp) function bilin_interp(xx,yy,f,x,y)
  implicit none
  ! This function uses bilinear interpolation to estimate the value
  ! of a function f at point (x,y)
  ! f is assumed to be sampled on a regular grid, with the grid x and y values specified
  ! by the arrays xx and yy, respectively.
  ! Reference: https://en.wikipedia.org/wiki/Bilinear_interpolation         
  real(dp),dimension(:),intent(in) :: xx, yy
  real(dp),dimension(:,:),intent(in) :: f
  real(dp),intent(in) :: x, y

  real(dp) :: denom, x1, x2, y1, y2
  integer :: nx,ny,i,j

  nx = size(xx)
  ny = size(yy) 
   
  i = locate(xx,x)
  j = locate(yy,y)
  if(i == 0 .or. i == nx .or. j == 0 .or. j == ny) then
     write(*,*) 'bilin_interp: selected array index out of range'
   stop
  end if
  
  x1 = xx(i)
  x2 = xx(i+1)
  
  y1 = yy(j)
  y2 = yy(j+1) 
         
  denom = (x2 - x1) * (y2 - y1)
  
  bilin_interp = ( f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
&                f(i,j+1)*(x2-x)*(y-y1) + f(i+1,j+1)*(x-x1)*(y-y1) ) / denom

  end function bilin_interp
  

end module interpolation
