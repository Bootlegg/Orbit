SUBROUTINE velocity(dt,Np,ax,ay,u1,v1,u2,v2)
IMPLICIT NONE

INTEGER :: Np
REAL*8 :: dt
REAL*8, DIMENSION(6) :: u1
REAL*8, DIMENSION(6) :: v1
REAL*8, DIMENSION(6) :: u2
REAL*8, DIMENSION(6) :: v2
REAL*8, DIMENSION(6) :: ax
REAL*8, DIMENSION(6) :: ay


ax = 0
ay = 0

u2(:) = u1(:) + dt*ax(:)
v2(:) = v1(:) + dt*ay(:)	

!f2py intent(in) u1
!f2py intent(in) v1
!f2py intent(in) ax
!f2py intent(in) ay
!f2py intent(in) dt
!f2py intent(in) Np
!f2py intent(out) u2
!f2py intent(out) v2

END SUBROUTINE velocity