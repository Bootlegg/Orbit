SUBROUTINE velocity(ax,ay,u1,v1,u2,v2,dt,np)
IMPLICIT NONE

INTEGER, intent(in) :: np
REAL*8, intent(in) :: dt
REAL*8, intent(in), DIMENSION(np) :: u1
REAL*8, intent(in), DIMENSION(np) :: v1
REAL*8, intent(in), DIMENSION(np) :: ax
REAL*8, intent(in), DIMENSION(np) :: ay
REAL*8, intent(out), DIMENSION(np) :: u2
REAL*8, intent(out), DIMENSION(np) :: v2



u2(:) = u1(:) + dt*ax(:)
v2(:) = v1(:) + dt*ay(:)

END SUBROUTINE velocity