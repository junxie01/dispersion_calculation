module global_data
! parameter
integer, parameter :: NL = 100
! Define numerical precision
integer, parameter :: wp = kind(1.0d0)
! Pi
real(wp), parameter :: pi = 3.141592653_wp
! Define input and output unit numbers
integer, parameter :: unit_in = 5
integer, parameter :: unit_out = 7
! model parameter
integer :: mmax ! number of layer except the half space
real(wp), dimension(:), allocatable :: thk, alpha, beta, rho, depth 
! Frequency information
integer :: nf
real(wp) :: omega
real(wp), dimension(:), allocatable :: frequency 
! Phase velocity information
real(wp), dimension(:,:), allocatable :: vr
! zh ratio
real(wp), dimension(:,:), allocatable :: zhratio 
! judge
integer :: ierr
! Tolerance value with several uses throughout the program
real(wp), parameter :: tol = 1.0d-6 
! Number of modal phase velocities to find at each frequency
integer, parameter :: maxroot = 1
! Number of increments between vrmin and vrmax for finding roots of the secular function
integer, parameter :: numinc = 50

end module global_data
