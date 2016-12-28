SUBROUTINE tv(Np,G,m,u,v,x,y,Tkin,Vpot)
IMPLICIT NONE

INTEGER :: Np,i,j
REAL*8 :: Tkin
REAL*8 :: Vpot
REAL*8 :: G
REAL*8, DIMENSION(6) :: y
REAL*8, DIMENSION(6) :: x
REAL*8, DIMENSION(6) :: v
REAL*8, DIMENSION(6) :: u
REAL*8, DIMENSION(6) :: m


!kinetic energy
DO i=1,6
	Tkin = Tkin + 0.5*m(i)*(u(i)**2+v(i)**2)
ENDDO

!potential energy
DO i=1,6
	DO j=i+1,6
		Vpot = Vpot + (-G*m(i)*m(j)/SQRT((x(i)-x(j))**2+(y(i)-y(j))**2))
	ENDDO
ENDDO

!time average kinetic energy
!Tavg = SUM(Tkin)/n

!time average potential energy
!Vavg = SUM(Vpot)/n

! Into this module
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(in) u
!f2py intent(in) v
!f2py intent(in) dt
!f2py intent(in) Np
!f2py intent(in) m

! Out of this module
!f2py intent(out) Tkin
!f2py intent(out) Vpot

END SUBROUTINE tv