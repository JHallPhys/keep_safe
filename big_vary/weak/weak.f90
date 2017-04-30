program mpi_demag

use mpi

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

! loop indicies + stuff

integer :: i,j,k,l,m,n,array_size,analytic



!##############################################################
!                    MPI - Stuff
!##############################################################

! MPI variables 

integer :: num_procs,my_rank,mpifile,rank_id
integer, dimension(1:MPI_STATUS_SIZE) :: MPI_Status
integer (kind=MPI_OFFSET_KIND) :: disp
integer,dimension(1:3) :: split,final,start_it,max_it
real(kind=dp),dimension(1:3) :: b_sum

! logical for root process

logical :: root = .false.

! filename for each rank

character (len=50) :: fname

! timing 

real(kind=dp) :: time_start,time_stop
!############### Lets do this ###############################

! create the world 

call MPI_init(ierror)
if (ierror.ne.0) stop 'Error in MPI_init'

! get rank on each processor in world: rank = 0 , num procs -1

call MPI_comm_rank(MPI_comm_world,my_rank,ierror)
if (ierror.ne.0) stop 'Error in MPI_comm_rank'


! Get total number of processors

call MPI_comm_size(MPI_comm_world,num_procs,ierror)
if (ierror.ne.0) stop 'Error in MPI_comm_size'


! set root as proc 0

if (my_rank.eq.0) root = .true.




!##############################################################
!-------------------------------------------------------------
!              Generate Magnetic volume  
!-------------------------------------------------------------
!=============================================================
! The choice of system is decided by the value of analtic s.t
!-------------------------------------------------------------
! Analytic = 0 - General ellipsoid
! Analytic = 1 - Uniform magnetised sphere  
! Analytic = 2 - Infinite plate 
!=============================================================
!##############################################################

! open input file on the root process 
! input file name: " array.dat "
time_start = MPI_wtime()
if (root) then 
  
  print*, 'No procs:',num_procs
  print*,


  ! now lets set up the ellipsoid. We will work with the vals form Q5
  ! rx=ry=10nm ; rz=20nm

  major = 9
  minor = 9

  dim_x = minor
  dim_y = minor
  dim_z = major


  ! allocate grid 

  allocate (grid(-dim_x:dim_x,-dim_y:dim_y,-dim_z:dim_z),stat=ierror)
  if (ierror.ne.0) stop 'problem allocating grid'


  grid(:,:,:) = 1.0_dp
  
  fname = 'geom_init'//trim(id_str(num_procs))//'.dat'
  open(file=fname,unit=10,status='replace')

  do i = -minor,minor
    do j = -minor,minor
      do k = -major,major
      
          write(10,*) grid(i,j,k)

      end do 
    end do 
  end do 
      
  close(10)
  
 
  print*, sum(grid)
 

end if 



!##############################################################
!-------------------------------------------------------------
!                 Send copies to all ranks
!------------------------------------------------------------- 
!==============================================================
! Need to broadcast info in order of dependence 
! Array Size --> Allocatation --> Data
!==============================================================
!############################################################## 


! Send Dimensions to other procs

! x-dim
call MPI_Bcast(dim_x, 1, MPI_Integer, 0, MPI_comm_world, ierror)
if (ierror.ne.0) stop 'Error: Bcast error in dimension of x'

! y-dim

call MPI_Bcast(dim_y, 1, MPI_Integer, 0, MPI_comm_world, ierror)
if (ierror.ne.0) stop 'Error: Bcast error in dimension of y'

! z-dim

call MPI_Bcast(dim_z, 1, MPI_Integer, 0, MPI_comm_world, ierror)
if (ierror.ne.0) stop 'Error: Bcast error in in dimension of z'

!Allocate array on non root ranks
  
if (.not. root)  allocate(grid(-int(dim_x):int(dim_x),-int(dim_y):int(dim_y),-int(dim_z):int(dim_z)),stat=ierror)
if (ierror.ne.0) stop 'Error in array allocation on non root'

if (.not.root) then 
  fname = 'geom_init'//trim(id_str(num_procs))//'.dat'
  open(file=fname,unit=10,action='read')

  
  do i = -dim_x,dim_x
    do j = -dim_y,dim_y
      do k = -dim_z,dim_z
      
         read(10,*) grid(i,j,k)

      end do 
    end do 
  end do 

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

! split data across processors

split(1) = floor(real(2*dim_x+1,kind=dp)/real(num_procs,kind=dp))

split(2) = floor(real(2*dim_y+1,kind=dp)/real(num_procs,kind=dp))

split(3) = floor(real(2*dim_z+1,kind=dp)/real(num_procs,kind=dp))


! perform the algorithm for each processor
! the limits on the do loop might be different from n_proc*split for the final proc
! due to the use of floor above. So we seperate the cases

if (my_rank.lt.(num_procs-1)) then

  !???????????????????????????????????????????????????????
 ! rank_id = 20 + my_rank
 ! fname = 'points_rank'//trim(id_str(my_rank))//'.dat'
 ! open(file=fname,unit=rank_id,status='replace') ! unique file for each rank
  !???????????????????????????????????????????????????????

  print*, split(3)
  print*, sum(grid)/2.0_dp
  max_it(:) = (split(:)*(my_rank+1))
  start_it(:) = (split(:)*my_rank)
  print*, my_rank,sum(grid(-dim_x:dim_x,-dim_y:dim_y,-dim_z+start_it(3):-dim_z+max_it(3)-1))
  do l = -dim_x,dim_x
      
    do m = -dim_y,dim_y

      do n = -dim_z+start_it(3),-dim_z+max_it(3)
 
        r_pos(:) = 0.0_dp
    
        if (grid(l,m,n).eq.0.0_dp) then

          b_field(:) = 0.0_dp  
          b_0(:) = 0.0_dp
        
        else if (grid(l,m,n).eq.1.0_dp) then
        
  !        write(rank_id,*) l,m,n         

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
                   
                    !?????????????????????????????????????????????????????????????????????????????????????

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

  !close(rank_id)

else if (my_rank.eq.(num_procs-1)) then

   !???????????????????????????????????????????????????????
   !rank_id = 20 + my_rank
   !fname = 'points_rank'//trim(id_str(my_rank))//'.dat'
   !open(file=fname,unit=rank_id,status='replace') ! unique file for each rank
   !???????????????????????????????????????????????????????

    print*, my_rank,sum(grid(-dim_x:dim_x,-dim_y:dim_y,-dim_z+(split(3)*(num_procs-1)):dim_z))
    do l = -dim_x,dim_x
    
      do m = -dim_y,dim_y

        do n =-dim_z+(split(3)*(num_procs-1))+1,dim_z
 
        r_pos(:) = 0.0_dp
    
        if (grid(l,m,n).eq.0.0_dp) then

          b_field(:) = 0.0_dp  
          b_0(:) = 0.0_dp
        
        else if (grid(l,m,n).eq.1.0_dp) then
    !      write(rank_id,*) l,m,n
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
                   
                    !?????????????????????????????????????????????????????????????????????????????????????

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

!close(rank_id)
   

end if




call MPI_Barrier(MPI_comm_world, ierror)
if (ierror.ne.0) stop 'error in barrier after sphere write test' 



do i = 1,3
  call mpi_reduce(b_av(i),b_sum(i),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_comm_world,ierror)
  if (ierror.ne.0) stop 'something went wrong' 
end do 




if (root) then
  
  print*,
  print*, '++++++++++++++++++++++++++++++++++++++++++++++++++'
  print*, 'Dx:',  1.0_dp-(pfactor_3*b_sum(1)/sum(grid))
  print*, 'Dy:',  1.0_dp-(pfactor_3*b_sum(2)/sum(grid))
  print*, 'Dz:',  1.0_dp-(pfactor_3*b_sum(3)/sum(grid))
  print*, '++++++++++++++++++++++++++++++++++++++++++++++++++'
  print*,
  print*, 'Sum(D_i):', sum(1.0_dp-(pfactor_3*b_sum(:)/sum(grid)))


end if 

time_stop = MPI_wtime()

if (root) then
print*, 
print*, 'Hey Listen!'
print*,
print*, 'total time:', time_stop-time_start,'seconds'
print*,
end if

! deallocate arrays

deallocate(grid,stat=ierror)
if (ierror.ne.0) stop 'error deallocating grid'

!Shutdown MPI
call MPI_finalize(ierror)
if (ierror.ne.0) stop 'Error in mpi finalize'




























contains


!################ Spherical initialisation subroutine ######################

subroutine init_sphere(dx,dy,dz,A)
implicit none 
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

! generate approx sphere

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
implicit none
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

div = 1000.0_dp


dtheta = 2.0_dp*pi/(div) ! need to generalise this 
dphi = pi/(div)



do i = 1,nint(div)
  theta = 0.0_dp
  do j = 1,nint(div)

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
implicit none 
  ! args
  real(kind=dp),dimension(1:3) :: b_dip
  real(kind=dp),dimension(1:3),intent(in) :: m_unit,r_unit,r_pos
 
  ! local
  real(kind=dp) :: r3


  r3 = sqrt(sum(r_pos**2))

  r3 = r3*r3*r3
 
  b_dip = (3.0_dp*r_unit(:)*dot_product(r_unit,m_unit)-m_unit(:))/r3
  
 
  

end function 

!################ Self-interaction term ######################

function b_self(m_unit)
implicit none
  real(kind=dp),dimension(1:3) :: b_self
  real(kind=dp),dimension(1:3),intent(in) :: m_unit
  

  b_self = pfactor_2*m_unit(:)


end function

!############### convert ###################

function id_str(n)
  
  character(len=20) :: id_str
  integer,intent(in) :: n

  write (id_str, *) n
  id_str = adjustl(id_str)

end function 


end program 









