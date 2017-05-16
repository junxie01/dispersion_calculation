!calculate secular function
function secular(k)
use global_data, only : wp, mmax,rho,beta,alpha,omega
implicit none
real(wp), intent(in) :: k
real(wp) :: k_local,epsilon,factor
real(wp) :: secular
complex(wp), dimension(:,:) :: a(2,2)
complex(wp), dimension(:), allocatable :: nus,nup
complex(wp), dimension(:,:,:), allocatable :: e11,e12,e21,e22,re11,re12,re21,re22,rd,du
interface
       subroutine psv(k,e11,e12,e21,e22,re11,re12,re21,re22,du,nus,nup)
       real(kind=kind(1.0d0)), intent(in) :: k
       complex(kind=kind(1.0d0)),dimension(:), intent(out) :: nus,nup
       complex(kind=kind(1.0d0)), dimension(:,:,:), intent(out) :: e11,e12,e21,e22,re11,re12,re21,re22,du
       end subroutine psv
       subroutine modrt(e11,e12,e21,e22,re11,re12,re21,re22,du,rd)
       complex(kind=kind(1.0d0)), dimension(:,:,:), intent(in) :: e11,e12,e21,e22,re11,re12,re21,re22,du
       complex(kind=kind(1.0d0)), dimension(:,:,:), intent(out) :: rd
       end subroutine modrt
end interface
allocate(e11(2,2,mmax+1), e12(2,2,mmax+1), e21(2,2,mmax+1), e22(2,2,mmax+1), re11(2,2,mmax+1))
allocate(re12(2,2,mmax+1), re21(2,2,mmax+1), re22(2,2,mmax+1), du(2,2,mmax+1), rd(2,2,mmax+1))
allocate(nus(mmax+1),nup(mmax+1))
k_local = k

! Check to see if the trial phase velocity is equal to the shear wave velocity or
! compression wave velocity of one of the layers. epsilon is arbitrarily set equal to 0.0001
epsilon = 0.0001_wp
do while (any(abs(beta-omega/k_local) < epsilon) .or. any(abs(alpha-omega/k_local) < epsilon))
        k_local = k_local * (1.0_wp + epsilon)
end do
call psv(k_local,e11,e12,e21,e22,re11,re12,re21,re22,du,nus,nup)
call modrt(e11,e12,e21,e22,re11,re12,re21,re22,du,rd)
a=e21(:,:,1)+matmul(matmul(e22(:,:,1),du(:,:,1)),rd(:,:,1))
!factor = nus(1)*nup(1)*rho(1)*beta(1)*beta(1)*rho(1)*beta(1)*beta(1)
factor = rho(1)*beta(1)*beta(1)*rho(1)*beta(1)*beta(1)
!factor = 1.0_wp
secular=cabs((a(1,1)*a(2,2) - a(1,2)*a(2,1))/factor/omega)
deallocate(e11,e12,e21,e22,re11,re12,re21,re22,du,rd,nus,nup)
end function secular
