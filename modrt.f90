!calculate the rdu
subroutine modrt(e11,e12,e21,e22,re11,re12,re21,re22,du,rd)
use global_data, only : wp, mmax
complex(wp), dimension(:,:,:), intent(out) :: rd
complex(wp), dimension(:,:,:), intent(in) :: e11, e12, e21, e22, re11, re12, re21, re22, du
complex(wp), dimension(:,:,:), allocatable :: td
complex(wp) :: a(4,4), b(4,4), c(4,4), d(4,2), e(2,2), f(2,2),aa(2),bb(2)
integer :: i
interface
      subroutine solve(e,bb,aa)
         complex(kind=kind(1.0d0)), dimension(2,2), intent(in) :: e
         complex(kind=kind(1.0d0)), dimension(2), intent(in) :: bb
         complex(kind=kind(1.0d0)), dimension(2), intent(out) :: aa
      end subroutine solve
end interface

allocate(td(2,2,mmax+1))
rd(:,:,mmax+1)=0.0_wp
do i=mmax,1,-1
   a(1:2,1:2) = re11(:,:,i)
   a(1:2,3:4) = re12(:,:,i)
   a(3:4,1:2) = re21(:,:,i)
   a(3:4,3:4) = re22(:,:,i)
   b(1:2,1:2) = e11(:,:,i+1)
   b(1:2,3:4) = e12(:,:,i+1)
   b(3:4,1:2) = e21(:,:,i+1)
   b(3:4,3:4) = e22(:,:,i+1)
   c=matmul(a,b)
   f=matmul(du(:,:,i+1),rd(:,:,i+1))
   !write(*,"(1I5.2,5f7.3)")i,wp,real(f(1,1)),real(f(1,2)),real(f(2,1)),real(f(2,2))
   e(1,1)=c(1,1)+c(1,3)*f(1,1)+c(1,4)*f(2,1)
   e(1,2)=c(1,2)+c(1,3)*f(1,2)+c(1,4)*f(2,2)
   e(2,1)=c(2,1)+c(2,3)*f(1,1)+c(2,4)*f(2,1)
   e(2,2)=c(2,2)+c(2,3)*f(1,2)+c(2,4)*f(2,2)
   call solve(e,du(:,1,i),td(:,1,i))
   call solve(e,du(:,2,i),td(:,2,i))
   e=matmul(f,td(:,:,i))
   d(1:2,:)=td(:,:,i)
   d(3:4,:)=e(:,:)
   rd(:,:,i)=matmul(c(3:4,:),d)
enddo
deallocate(td)
end subroutine modrt
