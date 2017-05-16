subroutine modal
use global_data, only : wp, alpha, beta, mmax, frequency, vr, omega, pi, numinc, tol, maxroot,nf,rho,thk
implicit none

! Internal variables
integer :: j, m, numroot,i
real(wp) :: kmax, kmin, dk, k1, k2, k3, ktrial, kold, f1, f2, f3, ftrial, &
            vrmin, vrmax, vpmin, vsmin, vsmax, nu

interface
	real(kind=kind(1.0d0)) function secular(k)
		real(kind=kind(1.0d0)), intent(in) :: k
	end function secular
	real(kind=kind(1.0d0)) function brent(ax,bx,f,t)
		real(kind=kind(1.0d0)) :: ax, bx, t
		real(kind=kind(1.0d0)), external :: f
	end function brent
end interface

vr = 0.0_wp
vpmin = minval(alpha(1:mmax+1))
vsmin = minval(beta(1:mmax+1))
vsmax = maxval(beta(1:mmax+1))

!write(*,*)vpmin,vsmin,vsmax
! Calculate Poisson's ratio of the layer with the minimum velocity
nu = 0.5_wp*((vpmin*vpmin-2.0_wp*vsmin*vsmin)/(vpmin*vpmin-vsmin*vsmin))

! Estimate the minimum Rayleigh phase velocity using the approximate method of Achenbach (1973)
vrmin = 0.98_wp*(vsmin*(0.862_wp + 1.14_wp*nu)/(1+nu))

! Estimate (crudely) the maximum Rayleigh phase velocity
vrmax = 1.00_wp*vsmax
!do i=1,mmax+1
!write(*,*)thk(i),alpha(i),beta(i),rho(i)
!enddo

!write(*,*)nf
! Loop through the frequencies
do j = 1,nf

	numroot = 0
	omega = 2.0_wp*pi*frequency(j)
		      
	! Establish the search parameters
	kmax = omega/vrmin
	kmin = omega/vrmax
	dk = (kmax - kmin)/numinc
      
	! Establish the first and second points
	k1 = kmax
	f1 = secular(k1)
	k2 = kmax - dk
	f2 = secular(k2)
	!write(*,*)k1,f1,k2,f2

	! Establish an arbitrary high value for kold
	kold = 1.1_wp*kmax
                  
	! Loop through the remaining points
	do m = 2,numinc-1
		k3 = kmax - m*dk
		f3 = secular(k3)
			
		! Determine if a minimum is bracketed
		if ((f2 < f1) .and. (f2 < f3)) then

			! Use quadratic interpolation to refine minimum
			ktrial = brent(k3,k1,secular,1.0d-12)
			ftrial = abs(secular(ktrial))
											                           
			! Check to see if ktrial is a zero and different from the previous zero
			if (ftrial < tol .and. abs((ktrial-kold)/kold) > tol) then
				numroot = numroot + 1
				vr(j,numroot) = omega/ktrial
				!write(*,*)vr(j,numroot)
				kold = ktrial
			end if
         
		end if
                           
		if (numroot == maxroot) exit

		k1 = k2
		f1 = f2
		k2 = k3
		f2 = f3

	end do

end do  

end subroutine modal
