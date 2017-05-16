!calculate the zh ratio of rayleigh waves
subroutine zhr
use global_data, only: wp, vr, frequency, zhratio, omega, nf, maxroot, pi, mmax, beta,alpha
integer :: i,j
real(wp) :: factor, k, epsilon, k_local 
complex(wp) :: a(2,2),b(2,2),f(2,2),cd(2),cu(2),c(4)
complex(wp) :: wt(4),e(4,4),te
complex(wp), dimension(:), allocatable :: nus, nup
complex(wp), dimension(:,:,:), allocatable :: e11, e12, e21, e22, rd,du,re11,re12,re21,re22


interface
        subroutine psv(k,e11,e12,e21,e22,re11,re12,re21,re22,du,nus,nup)
		real(kind=kind(1.0d0)), intent(in) :: k
                complex(kind=kind(1.0d0)), dimension(:), intent(out) :: nus, nup
                complex(kind=kind(1.0d0)), dimension(:,:,:), intent(out) :: e11, e12, e21, e22, re11,re12,re21,re22,du
        end subroutine psv
        subroutine modrt(e11,e12,e21,e22,re11,re12,re21,re22,du,rd)
                complex(kind=kind(1.0d0)), dimension(:,:,:), intent(in) :: e11, e12, e21, e22, re11,re12,re21,re22,du
                complex(kind=kind(1.0d0)), dimension(:,:,:), intent(out) :: rd
        end subroutine modrt
        subroutine inv2x2(a,b)
                complex(kind=kind(1.0d0)), dimension(2:2),intent(in) :: a
                complex(kind=kind(1.0d0)), dimension(2),intent(in) :: b
        end subroutine inv2x2
end interface
allocate(nus(mmax+1), nup(mmax+1), e11(2,2,mmax+1), e12(2,2,mmax+1), e21(2,2,mmax+1),&
         re11(2,2,mmax+1), re12(2,2,mmax+1), re21(2,2,mmax+1),re22(2,2,mmax+1),&
         e22(2,2,mmax+1), rd(2,2,mmax+1), du(2,2,mmax))

zhratio=0.0_wp
epsilon = 0.0000001_wp

do i=1,nf
   do j=1,maxroot
!      write(*,*)vr(i,j)
      omega=2.0_wp*pi*frequency(i)
      k=omega/vr(i,j)
      k_local = k
      do while (any(abs(beta-omega/k_local) < epsilon) .or. any(abs(alpha-omega/k_local) < epsilon))
      k_local = k_local * (1.0_wp + epsilon)
      enddo
      
      call psv(k_local,e11,e12,e21,e22,re11,re12,re21,re22,du,nus,nup)
      call modrt(e11,e12,e21,e22,re11,re12,re21,re22,du,rd)

!      call inv2x2(e21(:,:,1),b)
      b=matmul(matmul(e22(:,:,1),du(:,:,1)),rd(:,:,1))
      b=b/omega
      a=e21(:,:,1)/omega
      factor=sqrt((a(1,1)+b(1,1))*(a(1,1)+b(1,1))+(a(1,2)+b(1,2))*(a(1,2)+b(1,2)))
      cd(1)=(a(1,2)+b(1,2))/factor
      cd(2)=(-a(1,1)-b(1,1))/factor
      !write(*,*)rd(1,1,1)
      !write(*,*)rd(1,2,1)
      !write(*,*)rd(2,1,1)
      !write(*,*)rd(2,2,1)
      !write(*,*)'factor=',factor
!      factor=sqrt(abs(1.0_wp-f(1,1))*abs(1.0_wp-f(1,1))+abs(f(1,2))*abs(f(1,2)))
!      te=1.0_wp/sqrt(abs(1.0_wp-f(1,1))*abs(1.0_wp-f(1,1))+abs(f(1,2))*abs(f(1,2)))
      cu=matmul(rd(:,:,1),cd)
      c(1:2)=cd
      c(3:4)=cu
!      e(1:2,1:2)=matmul(e11(:,:,1),du(:,:,1))
!      e(3:4,1:2)=matmul(e21(:,:,1),du(:,:,1))
      e(1:2,1:2)=e11(:,:,1)
      e(3:4,1:2)=e21(:,:,1)
!      e(1:2,3:4)=e12(:,:,1)
!      e(3:4,3:4)=e22(:,:,1)
      e(1:2,3:4)=matmul(e12(:,:,1),du(:,:,1))
      e(3:4,3:4)=matmul(e22(:,:,1),du(:,:,1))
      wt=matmul(e,c)
!      write(*,*)'wt='
!      write(*,*)wt
!      write(*,*)wt(2),wt(1)
      zhratio(i,j)=abs(wt(2))/abs(wt(1))
   enddo
enddo
deallocate(nus, nup, e11, e12, e21, e22,re11,re12,re21,re22,rd,du)
end subroutine zhr
