program demag_sim_serial



implicit none 


!===============================================================
!                   Variable decleration
!===============================================================

! much precision 

integer,parameter :: dp=selected_real_kind(15,300)


!################ Nature numbers ######################

! Pi
real(kind=dp),parameter :: pi = 3.141592653589793238462643_dp

! Permeability of free space
real(kind=dp),parameter :: mu_0 = 4.0_dp*pi*1.0E-7_dp

! Angstrom 

real(kind=dp),parameter :: Ang = 1.0E-10_dp

! Dipole spacing 

real(kind=dp),parameter :: l_space = 3.0_dp*Ang


! bohr magneton

real(kind=dp),parameter :: mu_b = 9.27_dp*1.0E-24_dp


!############### Numerical factors ########################


! Lets get these pre-factors out the way early 

real(kind=dp),parameter :: l3 = l_space*l_space*l_space

real(kind=dp),parameter :: pfactor_1 = mu_0*mu_b/(4.0_dp*pi)

real(kind=dp),parameter :: pfactor_2 = 2.0_dp*mu_0*mu_b/(3.0_dp*l3)

real(kind=dp),parameter :: pfactor_3 = l3/(mu_0*mu_b)

real(kind=dp),parameter :: pfactor_4 = mu_0*mu_b/(4.0_dp*pi*l3)


!############### Model stuff ########################

! lattice 

real(kind=dp),dimension(:,:,:),allocatable :: grid

! system dimensions 

integer :: dim_x,dim_y,dim_z,ierror

! vectors 

real(kind=dp),dimension(1:3) :: r_unit,m_unit,r_pos,b_field,b_av,b_0

real(kind=dp) :: no_sites

!############### Elippsoid parametisation #############

! ellipsoid origins in x,y,z

real(kind=dp) :: x_0,y_0,z_0
real(kind=dp) :: r_x,r_y,r_z
integer :: minor,major
!############### Program vars #########################

! loop indicies

integer :: i,j,k,l,m,n,count_it,p
 character(len=30) :: arg
! debug

integer :: analytic,broken
real(kind=dp) :: time_start,time_stop


!===============================================================
!                   Initialise System
!===============================================================

! Command line args ~ Note they are lower case

do
  call get_command_argument(p, arg)
  ! No args
  if (len_trim(arg) == 0) exit
  ! perform unit tests
  if (trim(arg) == '--test') then
    call unit_test()
    stop
  end if
  ! analytic case of sphere Dx=Dy=Dz=1/3
  if (trim(arg) == '--sphere') then
   analytic=1
  end if
  ! analytic case of infinite plate Dx=Dy=0;Dz=1
  if (trim(arg) == '--plate') then
   analytic=2
  end if

  p = p+1

end do



print*, 
print*, '###############################################################################'
print*, 'Serial program to calculate the components of an ellipsoid of arbitrary shape.'
print*, '###############################################################################'
print*,

!############### uniform sphere #########################
call cpu_time(time_start)

if (analytic.eq.1) then

  print*,
  print*, '---------------------------------------'
  print*, ' System Type: Uniform Magnetised Shere'
  print*, '---------------------------------------'
  print*,

  ! test dim 

  dim_x = 12
  dim_y = dim_x
  dim_z = dim_y

  ! allocate array

  allocate(grid(-int(dim_x):int(dim_x),-int(dim_y):int(dim_y),-int(dim_z):int(dim_z)),stat=ierror)
  if (ierror.ne.0) stop ' error allocating grid' 

  call init_sphere(dim_x,dim_y,dim_z,grid)

  print*, 'no atoms:',sum(grid)


   !####################  Infinite Plane #########################


else if (analytic.eq.2) then

  print*,
  print*, '---------------------------------------'
  print*, ' System Type: Uniform Magnetised plane '
  print*, '---------------------------------------'
  print*,



  dim_x = 25
  dim_y=dim_x
  dim_z = 1

  allocate(grid(-dim_x:dim_x,-dim_y:dim_y,-dim_z:dim_z),stat=ierror)
  if (ierror.ne.0) stop 'error allocating inf plane grid'

  grid(:,:,:) = 0.0_dp
  grid(:,:,0) = 1.0_dp

  print*, 'no atoms:',sum(grid)

   !####################  General elipsoid #########################

else 

  print*,
  print*, '---------------------------------------'
  print*, ' System Type: Ellipsoid'
  print*, '---------------------------------------'
  print*,

  ! now lets set up the ellipsoid. We will work with the vals form Q5
  ! rx=ry=10nm ; rz=20nm

  major = 10
  minor = 20

  dim_x = minor
  dim_y = minor
  dim_z = major


  ! allocate grid 

  allocate (grid(-dim_x:dim_x,-dim_y:dim_y,-dim_z:dim_z),stat=ierror)
  if (ierror.ne.0) stop 'problem allocating grid'

  ! populate with points satisfying ellipsoid condition

  call init_ellipsoid(major,minor,grid)

  print*, 'no atoms:',sum(grid)

end if 

!===============================================================
!                   Dipole-field Algor
!===============================================================


!######### initialise vectors ###############

! m-hat

m_unit(:) = 1.0_dp


! initialise b-field

b_field(:) = 0.0_dp
b_av(:) = 0.0_dp

!########### Begin main loop ################

! move over lattice

! i,j,k ---> x,y,z from current atom
! l,m,n ---> pos of atom of intrest on grid


do l = -dim_x,dim_x
  
  do m = -dim_y,dim_y

    do n = -dim_z,dim_z
 
      r_pos(:) = 0.0_dp
    
      if (grid(l,m,n).eq.0.0_dp) then
     
        b_field(:) = 0.0_dp  
        b_0(:) = 0.0_dp    

      else if (grid(l,m,n).eq.1.0_dp) then

        do k = -dim_x-l,dim_x-l
 
          r_pos(1) = real(k,kind=dp)
                 
          do j = -dim_y-m,dim_y-m

            r_pos(2) = real(j,kind=dp)

            do i = -dim_z-n,dim_z-n
            
              r_pos(3) = real(i,kind=dp)
                        
                      if (grid(l+k,m+j,n+i).eq.1.0_dp)then
		       		 
                        if (sum(abs(r_pos)).eq.0) then 
                         
                          b_0(:) = b_self(m_unit(:))
                       
                        else
                          ! r-hat
                         
                          r_unit(:) = r_pos(:)/sqrt(sum(r_pos**2))
                      
                       
                          b_field(:) = b_field(:) + b_dip(r_unit(:),m_unit(:),r_pos(:)) 
                         
                        end if

                      end if
         
            end do 
        
          end do 

        end do 
     
      end if  
      
    ! update average 

          
    b_field(:) = b_field(:)*pfactor_4

    b_av(:) = b_av(:) + b_field(:) + b_0(:)
  
    ! reset
    b_field(:) = 0.0_dp
    b_0(:) = 0.0_dp
      
    end do 

  end do 
  
end do 



! now average over the total number of atoms

! number of atoms in sample

no_sites = sum(grid)

if (no_sites.eq.0) stop 'something has gone v/ wrong'

b_av(:) = b_av(:)/(no_sites)

call cpu_time(time_stop)


!############## Lets have a look at what happened ########################## 


print*,'------------------------------------------------------------'
print*, 'B_field comp :','x' 
print*, 'Demagnetising factor ',1.0_dp-(pfactor_3*b_av(1))
print*,'------------------------------------------------------------'
print*, 'B_field comp :','y' 
print*, 'Demagnetising factor ',1.0_dp-(pfactor_3*b_av(2))
print*,'------------------------------------------------------------'
print*, 'B_field comp :','z' 
print*,'------------------------------------------------------------'
print*, 'Demagnetising factor ',1.0_dp-(pfactor_3*b_av(3))


deallocate(grid,stat=ierror)
if (ierror.ne.0) stop 'error deallocating grid'

print*, 
print*, 'total time:', time_stop-time_start,'seconds'




!##########################################################################
!==========================================================================
!                    Functions and Subroutines 
!==========================================================================
!##########################################################################


contains


!################ Spherical initialisation subroutine ######################

subroutine init_sphere(dx,dy,dz,A)

integer,intent(in) :: dx,dy,dz
real(kind=dp),dimension(-dx:dx,-dy:dy,-dz:dz) :: A

! locals 

real(kind=dp) :: rx,ry,rz,points,phi,dphi,theta,dtheta

!points is the number of points that a radius vector will be split into
points = real(dx,kind=dp)

A(:,:,:) = 0.0_dp

phi = 0.0_dp

dphi = pi/(2.0_dp*points)
dtheta = pi/(2.0_dp*points)


do j = 0,dx

theta=0.0_dp

  do i = 0,dy
 
    rx = cos(theta)*sin(phi)
    ry = sin(theta)*sin(phi)
    rz = cos(phi)

  

    theta = theta + dtheta

    

    grid(0:nint(points*rx),0:nint(points*ry),-nint(points*rz):nint(points*rz))   = 1.0_dp
    grid(0:nint(points*rx),-nint(points*ry):0,-nint(points*rz):nint(points*rz))  = 1.0_dp
    grid(-nint(points*rx):0,-nint(points*ry):0,-nint(points*rz):nint(points*rz)) = 1.0_dp
    grid(-nint(points*rx):0,0:nint(points*ry),-nint(points*rz):nint(points*rz))  = 1.0_dp

  end do 
  phi = phi + dphi

end do 

end subroutine 


!################ Ellipsoid generation subroutine ######################

subroutine init_ellipsoid(major,minor,A)

! alawys center ellipsoid about x=y=0. Easier to think about 

! args
integer,intent(in) :: major,minor
real(kind=dp),dimension(-minor:minor,-minor:minor,-major:major),intent(inout) :: A


! locals

real(kind=dp) :: theta,phi,dtheta,dphi,r_x,r_y,r_z,x,y,z,elipse,div
integer :: i,j,k

! initalise array

A(:,:,:) = 0.0_dp

! set r_x,r_y as the minor axis

r_x = real(minor,kind=dp)
r_y = r_x

! set r_z as major
r_z = real(major,kind=dp)

! initialise angles

theta = 0.0_dp
phi= 0.0_dp

! div controls how good the approximation is 

div = 1000.0_dp

dtheta = 2.0_dp*pi/(div) 
dphi = pi/(div)

do i = 1,int(div)
  theta = 0.0_dp
  do j = 1,int(div)

    x = r_x*cos(theta)*sin(phi)
    y = r_y*sin(theta)*sin(phi)
    z = r_z*cos(phi)

    
      elipse = (x)**2/(r_x**2) + (y)**2/(r_y**2) + (z)**2/(r_z**2)


          
      if ( (theta.le.pi).and.(phi.le.(pi/2.0_dp)) ) then 

         grid(0:int(x),0:int(y),-int(z):int(z))   = 1.0_dp
         grid(0:int(x),-int(y):0,-int(z):int(z))  = 1.0_dp
         grid(-int(x):0,-int(y):0,-int(z):int(z)) = 1.0_dp
         grid(-int(x):0,0:int(y),-int(z):int(z))  = 1.0_dp


      end if 

    theta = theta + dtheta
   
  end do 
    
  phi = phi + dphi
  
end do 



end subroutine 


!################ Dipole field function ######################

function b_dip(r_unit,m_unit,r_pos)

  ! args
  real(kind=dp),dimension(1:3) :: b_dip
  real(kind=dp),dimension(1:3),intent(in) :: m_unit,r_unit,r_pos
 
  ! local
  real(kind=dp) :: r3


  r3 = sqrt(sum(r_pos**2)) 

  ! catch 
  if (r3.eq.0) stop 'Function b_dip: Division by Zero'

  r3 = r3*r3*r3
 
  b_dip = (3.0_dp*r_unit(:)*dot_product(r_unit,m_unit)-m_unit(:))/r3
  
 
  

end function 

!---------- unit test for dipole-field function---------------

function test_b_dip() result(pass)

  implicit none

  ! Program independent Quantities for testing with 
  real(kind=dp) :: m_unit(1:3)=(/1.0_dp, 1.0_dp, 1.0_dp/)
  real(kind=dp) :: r_unit(1:3)=(/1.0_dp,0.0_dp,1.0_dp/)
  real(kind=dp) :: r_pos(1:3)=(/2.0_dp,0.0_dp,0.0_dp/)
  logical       :: pass

  ! Give the function the benifit of the doubt
  pass = .True.

  ! Check we get out a 3-vector 
  if (size(b_dip(r_unit,m_unit,r_pos)).ne.3) pass=.false.

  ! Check zero moment gives us zero field
  if (sum(abs(b_dip(r_unit,(/0.0_dp,0.0_dp,0.0_dp/),r_pos))).ne.0) pass=.false.

  ! Check non-zero moment gives us non-zero field
  if (sum(abs(b_dip(r_unit,m_unit,r_pos))).le.0) pass=.false.

  ! Check Precision
  if (precision(b_dip(r_unit,m_unit,r_pos)).ne.15) pass=.false.
   
end function 



!################ Self-interaction term ######################

function b_self(m_unit)

  real(kind=dp),dimension(1:3) :: b_self
  real(kind=dp),dimension(1:3),intent(in) :: m_unit
  

  b_self = pfactor_2*m_unit(:)


end function


!------------ Self-interaction unit test ----------------------

function test_b_self() result(pass)

  implicit none 

  real(kind=dp) :: m_unit(1:3)=(/1.0_dp,1.0_dp,1.0_dp/)

  logical :: pass

  ! dont fail me now 
  pass = .true.

  ! check the result is a 3-vector
  if (size(b_self(m_unit)).ne.3) pass=.false.
  ! check zero moment gives zero field
  if (sum(abs(b_self( (/0.0_dp,0.0_dp,0.0_dp/) ))).ne.0) pass=.false.
  ! check non-zero moment gives non-zero field 
  if (sum(abs(b_self(m_unit))).le.0) pass=.false.
  ! check precision 
  if (precision(b_self(m_unit)).ne.15) pass=.false.
  
end function

!################ Unit Test Routine ###########################

subroutine unit_test()
implicit none 


print*,
print*, 
print*, ' ------------------ Test Test Test -------------------'
print*, 
print*, 'if a test has passed ---------------------> T'
print*, 'if a test has failed ---------------------> F'
print*, 
print*, 
print*, 'Magnetic field tests for r>0:            ',test_b_dip()
print*,
print*, 'Magnetic field tests for r=0:            ',test_b_self()


end subroutine



end program 









