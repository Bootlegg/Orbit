SUBROUTINE num(A,B,C)
IMPLICIT NONE

REAL*8 :: A
REAL*8, DIMENSION(6) :: B
REAL*8 :: C
B = 1
A = 300000
WRITE(6,*) A
WRITE(6,*) B
WRITE(6,*) C


!f2py intent(out) A
!f2py intent(out) B
!f2py intent(in) C

END SUBROUTINE num