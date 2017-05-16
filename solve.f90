! this program is used to solve 2X2 matrix function with complex element
subroutine solve(matr,b,a)
use global_data, only : wp
implicit none
complex(wp),dimension(2,2),intent(in) :: matr
complex(wp),dimension(2),intent(in) :: b
complex(wp),dimension(2),intent(out) :: a
complex(wp),dimension(2,2) :: m(2,2)
complex(wp) :: t
t=matr(1,1)*matr(2,2)-matr(1,2)*matr(2,1)
m(1,1)=matr(2,2)/t
m(1,2)=-matr(1,2)/t
m(2,1)=-matr(2,1)/t
m(2,2)=matr(1,1)/t
a(1)=m(1,1)*b(1)+m(1,2)*b(2)
a(2)=m(2,1)*b(1)+m(2,2)*b(2)
end subroutine solve
