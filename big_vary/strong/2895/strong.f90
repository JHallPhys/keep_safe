program scale_strong_demag

implicit none

integer,parameter :: dp=selected_real_kind(15,300)

real(kind=dp) :: ideal,ideal_log,strong,strong_log,ideal_0

integer :: i,n

open(file='2895.dat',unit=20,action='read')

read(20,*) n,ideal_0

print*, ideal_0

close(20)


open(file='2895.dat',unit=20,action='read')
open(file='strong_plot.dat',unit=21,status='replace')

do i = 1,10 

  read(20,*) n,strong
  write(21,*) n,strong,ideal_0/real(i,kind=dp),log(real(n,kind=dp)),log(strong),log(ideal_0/real(i,kind=dp))
  

end do 


close(20)
close(21)
end program
