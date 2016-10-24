SUBROUTINE accel(G,Np,m,x,y,ax,ay)
IMPLICIT NONE

INTEGER :: i,j,Np
REAL*8, DIMENSION(6) :: m
REAL*8, DIMENSION(6) :: x
REAL*8, DIMENSION(6) :: y
REAL*8, DIMENSION(6) :: ax
REAL*8, DIMENSION(6) :: ay
REAL*8 :: G

ax = 0
ay = 0

!WRITE(6,*) m 
DO i=1,6
	DO j=1,6
		IF (j /= i) THEN
			ax(i) = ax(i) + (-G*m(j)/((x(i)-x(j))**2+(y(i)-y(j))**2)**(1.5))*(x(i)-x(j))
			ay(i) = ay(i) + (-G*m(j)/((x(i)-x(j))**2+(y(i)-y(j))**2)**(1.5))*(y(i)-y(j))
		END IF
	END DO
END DO

!f2py intent(out) ax
!f2py intent(out) ay
!f2py intent(in) x
!f2py intent(in) y
!f2py intent(in) m

END SUBROUTINE accel