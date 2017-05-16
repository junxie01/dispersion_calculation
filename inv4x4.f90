subroutine inv4x4(a,b)

! This subroutine calculates the analytical inverse of a complex 4 x 4 matrix

! Copyright 2004 by Glenn J. Rix and Carlo G. Lai

! This file is part of SWAMI.
!
! SWAMI is free software; you can redistribute it and/or
! modify it under the terms of the GNU General Public License
! as published by the Free Software Foundation; either version 2
! of the License, or (at your option) any later version.
!
! SWAMI is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

! Entry
!	a:	4 x 4 complex matrix to be inverted
!
! Exit
!	b:	analytic inverse of A

use global_data, only : wp
implicit none

! External parameters
complex(wp), dimension(4,4), intent(in) :: a
complex(wp), dimension(4,4), intent(out) :: b

! Internal parameters
complex(wp) :: a11, a12, a13, a14, a21, a22, a23, a24, a31, a32, a33, a34, a41, a42, a43, a44, c

! Assign elements of a to local variables
a11 = a(1,1); a12 = a(1,2); a13 = a(1,3); a14 = a(1,4)
a21 = a(2,1); a22 = a(2,2); a23 = a(2,3); a24 = a(2,4)
a31 = a(3,1); a32 = a(3,2); a33 = a(3,3); a34 = a(3,4)
a41 = a(4,1); a42 = a(4,2); a43 = a(4,3); a44 = a(4,4)

! Calculate terms in the inverse
c = a11*a22*a33*a44-a11*a22*a34*a43-a11*a32*a23*a44+a11*a32*a24*a43+a11*a42*a23*a34- &
    a11*a42*a24*a33-a21*a12*a33*a44+a21*a12*a34*a43+a21*a32*a13*a44-a21*a32*a14*a43- &
    a21*a42*a13*a34+a21*a42*a14*a33+a31*a12*a23*a44-a31*a12*a24*a43-a31*a22*a13*a44+ &
    a31*a22*a14*a43+a31*a42*a13*a24-a31*a42*a14*a23-a41*a12*a23*a34+a41*a12*a24*a33+ &
    a41*a22*a13*a34-a41*a22*a14*a33-a41*a32*a13*a24+a41*a32*a14*a23
    
b(1,1) = (a22*a33*a44-a22*a34*a43-a32*a23*a44+a32*a24*a43+a42*a23*a34-a42*a24*a33)
b(1,2) = -(a12*a33*a44-a12*a34*a43-a32*a13*a44+a32*a14*a43+a42*a13*a34-a42*a14*a33)
b(1,3) = (a12*a23*a44-a12*a24*a43-a22*a13*a44+a22*a14*a43+a42*a13*a24-a42*a14*a23)
b(1,4) = -(a12*a23*a34-a12*a24*a33-a22*a13*a34+a22*a14*a33+a32*a13*a24-a32*a14*a23)

b(2,1) = -(a21*a33*a44-a21*a34*a43-a31*a23*a44+a31*a24*a43+a41*a23*a34-a41*a24*a33)
b(2,2) = (a11*a33*a44-a11*a34*a43-a31*a13*a44+a31*a14*a43+a41*a13*a34-a41*a14*a33)
b(2,3) = -(a11*a23*a44-a11*a24*a43-a21*a13*a44+a21*a14*a43+a41*a13*a24-a41*a14*a23)
b(2,4) = (a11*a23*a34-a11*a24*a33-a21*a13*a34+a21*a14*a33+a31*a13*a24-a31*a14*a23)

b(3,1) = (a21*a32*a44-a21*a34*a42-a31*a22*a44+a31*a24*a42+a41*a22*a34-a41*a24*a32)
b(3,2) = -(a11*a32*a44-a11*a34*a42-a31*a12*a44+a31*a14*a42+a41*a12*a34-a41*a14*a32)
b(3,3) = (a11*a22*a44-a11*a24*a42-a21*a12*a44+a21*a14*a42+a41*a12*a24-a41*a14*a22)
b(3,4) = -(a11*a22*a34-a11*a24*a32-a21*a12*a34+a21*a14*a32+a31*a12*a24-a31*a14*a22)

b(4,1) = -(a21*a32*a43-a21*a33*a42-a31*a22*a43+a31*a23*a42+a41*a22*a33-a41*a23*a32)
b(4,2) = (a11*a32*a43-a11*a33*a42-a31*a12*a43+a31*a13*a42+a41*a12*a33-a41*a13*a32)
b(4,3) = -(a11*a22*a43-a11*a23*a42-a21*a12*a43+a21*a13*a42+a41*a12*a23-a41*a13*a22)
b(4,4) = (a11*a22*a33-a11*a23*a32-a21*a12*a33+a21*a13*a32+a31*a12*a23-a31*a13*a22)

b = b/c

end subroutine inv4x4