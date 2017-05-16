! determine parameters for psv wave
subroutine psv(k,e11,e12,e21,e22,re11,re12,re21,re22,du,nus,nup)
use global_data, only: wp, mmax, thk, alpha, beta, rho, omega
implicit none
! Internal parameters
integer :: i
real(wp), dimension(:), allocatable :: mu
complex(wp), dimension(:), allocatable :: gammas, gammap, chi, ta
! External parameters
real(wp), intent(in) :: k
complex(wp),dimension(:), intent(out) :: nus, nup
complex(wp), dimension(:,:,:), intent(out) :: e11, e12, e21, e22, re11, re12, re21, re22, du 

!write(*,*)'hello in psv stream!'
allocate(mu(mmax+1),gammas(mmax+1),gammap(mmax+1),chi(mmax+1),ta(mmax+1))
mu = rho*beta*beta

nus = sqrt(cmplx(k*k-(omega*omega)/(beta*beta),kind=wp))
where(aimag((0.0_wp,-1.0_wp)*nus) > 0.0_wp) nus = -nus
gammas = nus/k

nup = sqrt(cmplx(k*k-(omega*omega)/(alpha*alpha),kind=wp))
where(aimag((0.0_wp,-1.0_wp)*nup) > 0.0_wp) nup = -nup
gammap = nup/k
chi = 2.0_wp*k - (omega*omega)/(beta*beta)/k
ta = 1.0_wp/mu/(4.0_wp*k-2.0_wp*chi)

e11(1,1,:) = -1.0_wp
e11(1,2,:) = gammas
e12(1,1,:) = -1.0_wp
e12(1,2,:) = gammas
e11(2,1,:) = -gammap
e11(2,2,:) = 1.0_wp
e12(2,1,:) = gammap
e12(2,2,:) = -1.0_wp
e21(1,1,:) = 2.0_wp*mu*nup
e21(1,2,:) = -mu*chi
e22(1,1,:) = -2.0_wp*mu*nup
e22(1,2,:) = mu*chi
e21(2,1,:) = mu*chi
e21(2,2,:) = -2.0_wp*mu*nus
e22(2,1,:) = mu*chi
e22(2,2,:) = -2.0_wp*mu*nus

re11(1,1,:) = -2.0_wp*k*mu*ta
re11(1,2,:) = mu*chi*ta/gammap
re12(1,1,:) = ta/gammap
re12(1,2,:) = -ta 
re11(2,1,:) = -mu*chi*ta/gammas
re11(2,2,:) = 2.0_wp*k*mu*ta
re12(2,1,:) = ta 
re12(2,2,:) = -ta/gammas
re21(1,1,:) = -2*mu*k*ta
re21(1,2,:) = -mu*chi*ta/gammap
re22(1,1,:) = -ta/gammap
re22(1,2,:) = -ta
re21(2,1,:) = -mu*chi*ta/gammas
re21(2,2,:) = -2*mu*k*ta
re22(2,1,:) = -ta
re22(2,2,:) = -ta/gammas

du(1,1,:) = exp(-nup(1:mmax)*thk)
du(1,2,:) = 0.0_wp
du(2,1,:) = 0.0_wp
du(2,2,:) = exp(-nus(1:mmax)*thk)
deallocate(mu,gammas,gammap,chi,ta)
end subroutine psv
