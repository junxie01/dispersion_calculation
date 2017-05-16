! get period 
subroutine getpfile(pfile)
use global_data, only: wp, nf, frequency,ierr
implicit none
character(len=80),intent(in) :: pfile
integer i
logical ext
real(wp) :: p
inquire(file=pfile,exist=ext)
if(.not.ext)then
   ierr = -1
   write(*,*)'Period file does not exist.'
   return
endif
open(8,file=pfile,status='old',form='formatted')
 i=1
100  continue
read(8,*,err=200,end=200)p
frequency(i)=1.0_wp/p
i=i+1 
goto 100
 200  continue
close(8)
nf=i-1
end subroutine getpfile
