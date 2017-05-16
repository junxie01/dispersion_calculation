subroutine inv2x2(a,b)
use global_data, only : wp
implicit none
! External parameters
complex(wp), dimension(2,2), intent(in) :: a
complex(wp), dimension(2,2), intent(out) :: b
! Internal parameters
complex(wp) :: a11, a12, a21, a22, c
! Assign elements of a to local variables
a11 = a(1,1); a12 = a(1,2); a21 = a(2,1); a22 = a(2,2);
! Calculate terms in the inverse
c = a11*a22 - a12*a21
b(1,1) = a22
b(1,2) = -a12
b(2,1) = -a21
b(2,2) = a11
b = b/c
end subroutine inv2x2
