! caculate the phase velocity dispersion curve and Z/H ratio
program main
use global_data, only: wp, nf, vr, maxroot, frequency, NL, zhratio
implicit none
character(len=80) :: mname, pfile
integer i,j
interface
        subroutine getmod(mname)
                character(len=80),intent(in) :: mname
        end subroutine getmod
        subroutine getpfile(pfile)
                character(len=80),intent(in) :: pfile
        end subroutine getpfile
        subroutine flat
        end subroutine flat
        subroutine modal
        end subroutine modal
        subroutine zhr
        end subroutine zhr
end interface
allocate(frequency(NL),vr(NL,maxroot), zhratio(NL,maxroot))
mname='model'
i=iargc()
if(i.ne.2)then
   write(*,*)'usage: dsp_jx model pfile'
   return
else
   call getarg(1,mname)
   call getmod(mname)
   call getarg(2,pfile)
   call getpfile(pfile)
   call flat
   call modal
!   call flat
   call zhr
   do i=1,nf
      write(*,'(g10.5,1X,2g10.5,1X,2g10.5)')1.0/frequency(i),(vr(i,j),zhratio(i,j),j=1,maxroot)
!      write(*,*)1.0/frequency(i),(vr(i,j),zhratio(i,j),j=1,maxroot)
   enddo
endif
deallocate(frequency,vr,zhratio)
end
