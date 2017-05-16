subroutine flat
use global_data, only: wp, thk,rho,alpha,beta,mmax
implicit none
real(wp) z0,z1,r0,r1,dr,ar,tmp,r(mmax+1),tmp1,pwr
integer i
ar=6370.0_wp
dr=0.0_wp
thk(mmax+1)=1.0_wp
pwr=2.275_wp
do i=1,mmax+1
   dr=dr+thk(i)
   r(i)=ar-dr
enddo
r0=ar
do i=1,mmax
   z0=ar*dlog(ar/r0)
   z1=ar*dlog(ar/r(i))
   thk(i)=z1-z0
   tmp=ar*(1.0/r(i+1)-1.0/r(i))/dlog(r(i)/r(i+1))
   alpha(i)=alpha(i)*tmp
   beta(i)=beta(i)*tmp
   tmp1=(r(i)**pwr-r(i+1)**pwr)/dlog(r(i)/r(i+1))/pwr/ar**pwr
   rho(i)=rho(i)*tmp1
   r0 = r(i)
enddo
thk(mmax+1)=0.0_wp
end subroutine flat
