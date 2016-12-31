MODULE steps
IMPLICIT NONE
CONTAINS
	SUBROUTINE accel(m,x,y,ax,ay,G,np)
		IMPLICIT NONE

		INTEGER :: i,j,np
		REAL*8, DIMENSION(np) :: m
		REAL*8, DIMENSION(np) :: x
		REAL*8, DIMENSION(np) :: y
		REAL*8, DIMENSION(np) :: ax
		REAL*8, DIMENSION(np) :: ay
		REAL*8 :: G



		DO i=1,np
			DO j=1,np
				IF (j /= i) THEN
					ax(i) = ax(i) + (-G*m(j)/((x(i)-x(j))**2+(y(i)-y(j))**2)**(1.5))*(x(i)-x(j))
					ay(i) = ay(i) + (-G*m(j)/((x(i)-x(j))**2+(y(i)-y(j))**2)**(1.5))*(y(i)-y(j))
				END IF
			ENDDO
		ENDDO

	END SUBROUTINE accel


	SUBROUTINE velocity(ax,ay,u1,v1,u2,v2,dt,np)
		IMPLICIT NONE

		INTEGER :: np
		REAL*8 :: dt
		REAL*8, DIMENSION(np) :: u1
		REAL*8, DIMENSION(np) :: v1
		REAL*8, DIMENSION(np) :: ax
		REAL*8, DIMENSION(np) :: ay
		REAL*8, DIMENSION(np) :: u2
		REAL*8, DIMENSION(np) :: v2

		u2(:) = u1(:) + dt*ax(:)
		v2(:) = v1(:) + dt*ay(:)

	END SUBROUTINE velocity

	SUBROUTINE position(ax,ay,xnold,ynold,xnolder,ynolder,x2,y2,dt,np)
		IMPLICIT NONE

		INTEGER :: np
		REAL*8 :: dt
		REAL*8, DIMENSION(np) :: y2
		REAL*8, DIMENSION(np) :: x2
		REAL*8, DIMENSION(np) :: ynolder
		REAL*8, DIMENSION(np) :: xnolder
		REAL*8, DIMENSION(np) :: ynold
		REAL*8, DIMENSION(np) :: xnold
		REAL*8, DIMENSION(np) :: ax
		REAL*8, DIMENSION(np) :: ay

		x2 = 2.0*xnold - xnolder + dt*dt*ax
		y2 = 2.0*ynold - ynolder + dt*dt*ay

	END SUBROUTINE position
END MODULE steps