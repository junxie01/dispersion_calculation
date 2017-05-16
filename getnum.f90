! get the earth model
! the maxmam number of the layer is 100.
! the model was in the format:
!    [layer thickness] vp vs rho
subroutine getmod(mname, mmax)
use global_data, only: wp, unit_in, mmax, thk, alpha, beta, rho,ierr
implicit none
character(len=80) :: mname,string
logical ext 
integer i,ios

!write(*,*)NL
inquire(file=mname,exist=ext)
if(.not.ext)then
   ierr = -1
   write(*,*)'Model file does not exist.'
   return
endif
open(unit_in,file=mname,status='old',form='formatted')
rewind unit_in
i=0
loop:do
      read(unit_in,'(a)',iostat=ios) string
      if (ios < 0) exit loop
      i=i+1
     enddo loop
