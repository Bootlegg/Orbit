SUBROUTINE position(dt,Np,ax,ay,xnold,ynold,xnolder,ynolder,x2,y2)
IMPLICIT NONE

INTEGER :: Np
REAL*8 :: dt
REAL*8, DIMENSION(6) :: y2
REAL*8, DIMENSION(6) :: x2
REAL*8, DIMENSION(6) :: ynolder
REAL*8, DIMENSION(6) :: xnolder
REAL*8, DIMENSION(6) :: ynold
REAL*8, DIMENSION(6) :: xnold
REAL*8, DIMENSION(6) :: ax
REAL*8, DIMENSION(6) :: ay

!ax = 0
!ay = 0

x2 = 2.0*xnold - xnolder + dt*dt*ax
y2 = 2.0*ynold - ynolder + dt*dt*ay

! Into this module
!f2py intent(in) yolder
!f2py intent(in) xolder
!f2py intent(in) yold
!f2py intent(in) xold
!f2py intent(in) ax
!f2py intent(in) ay
!f2py intent(in) dt
!f2py intent(in) Np

! Out of this module
!f2py intent(out) x2
!f2py intent(out) y2

END SUBROUTINE position