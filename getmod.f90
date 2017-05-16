! get the earth model
! the maxmam number of the layer is 100.
! the model was in the format:
!    [layer thickness] vp vs rho
subroutine getmod(mname)
use global_data, only: wp, unit_in, mmax, thk, alpha, beta, rho,ierr
implicit none
character(len=80), intent(in) :: mname
character(len=80) :: string
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
allocate(thk(i-1),alpha(i-1),beta(i-1),rho(i-1))
mmax=i-2
if(mmax.eq.-1)then
   ierr=-1
   write(*,*)'Error in model file.'
   return
endif
!write(*,*)'layer without half-sapce of model is',mmax
!write(*,*)'hello'
rewind unit_in
read(unit_in,'(a)') string
!write(*,*)'initial model is:'
do i=1,mmax+1
   read(unit_in,*)thk(i),alpha(i),beta(i),rho(i)
!   write(*,*)thk(i),alpha(i),beta(i),rho(i)
enddo
close(unit_in)
!write(*,*)'mmax=',mmax
!deallocate(thk,alpha,beta,rho,depth)
end subroutine getmod
