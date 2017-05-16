! find out solutions
function brent(ax,bx,f,t)
use global_data, only : wp
implicit none

! External variables
real(wp), intent(in) :: ax, bx, t
real(wp), external :: f
real(wp) :: brent

! Internal variables
real(wp) :: eps, a, b, c, d, e, m, p, q, r, tol, t2, u, v, w, x, fu, fv, fw, fx

! Initialize variables
a = ax
b = bx
c = (3.0_wp - sqrt(5.0_wp))/2.0_wp
eps = sqrt(epsilon(1.0_wp))
v = a + c * (b-a)
w = v
x = v
e = 0.0_wp
fv = f(x)
fw = fv
fx = fv

! Main loop
do
	m = 0.5_wp * (a+b)
	tol = eps*abs(x) + t
	t2 = 2.0_wp*tol

	! Check stopping criterion
	if (max(x-a,b-x) <= t2) exit
	
	p = 0.0_wp
	q = 0.0_wp
	r = 0.0_wp
	if (abs(e) > tol) then
		! Fit parabola
		r = (x-w) * (fx-fv)
		q = (x-v) * (fx-fw)
		p = (x-v) * q - (x-w) * r
		q = 2.0_wp * (q-r)
		if (q > 0.0_wp) then
			p = -p
		else
			q = -q
		end if
		r = e
		e = d
	end if
	
	if (abs(p) < abs(0.5_wp*q*r) .and. p < q*(a-x) .and. p < q*(b-x)) then
		! Parabolic interpolation step
		d = p/q
		u = x + d
		! f must not be evaluated too close to a or b
		if (u-a < t2 .or. b-u < t2) then
			if (x < m) then
				d = tol
			else
				d = -tol
			end if
		end if
	else
		! Golden section step
		if (x < m) then
			e = b - x
		else
			e = a - x
		end if
		d = c * e
	end if
	! f must not be evaluated too close to x
	if (abs(d) >= tol) then
		u = x + d
	else if (d > 0.0_wp) then
		u = x + tol
	else
		u = x - tol
	end if
	fu = f(u)
	
	! Update a, b, v, w, and x
	if (fu <= fx) then
		if (u < x) then
			b = x
		else
			a = x
		end if
		v = w
		fv = fw
		w = x
		fw = fx
		x = u
		fx = fu
	else
		if (u < x) then
			a = u
		else
			b = u
		end if
		if (fu <= fw .or. w == x) then
			v = w
			fv = fw
			w = u
			fw = fu
		else if (fu <= fv .or. v ==x .or. v == w) then
			v = u
			fv = fu
		end if
	end if
end do
brent = x
end function brent
