SUBROUTINE mainr(ax,ay,u,v,x,y,m,Tavg,Vavg,Tkin,Vpot,dt,G,nt,np)
	USE steps
	IMPLICIT NONE


	INTEGER :: n, i, j
	REAL*8, intent(in) :: dt
	REAL*8, intent(in) :: G
	INTEGER, intent(in) :: np
	INTEGER, intent(in) :: nt

	REAL*8, intent(inout), DIMENSION(nt+1,np) :: ax
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: ay
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: u
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: v
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: x
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: y

	REAL*8, intent(in), DIMENSION(np) :: m
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: Tavg
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: Vavg
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: Tkin
	REAL*8, intent(inout), DIMENSION(nt+1,np) :: Vpot 


	! Into this module
	!f2py intent(in) G
	!f2py intent(in) dt
	!f2py intent(in) nt
	!f2py intent(in) np

	!in and out
	!f2py intent(in,out) Vpot
	!f2py intent(in,out) Tkin
	!f2py intent(in,out) Vavg
	!f2py intent(in,out) Tavg
	!f2py intent(in) m
	!f2py intent(in,out) y
	!f2py intent(in,out) x
	!f2py intent(in,out) v
	!f2py intent(in,out) u
	!f2py intent(in,out) ay
	!f2py intent(in,out) ax


	!kinetic energy
	DO i=1,np
		Tkin(1,1) = Tkin(1,1) + 0.5*m(i)*(u(1,i)**2+v(1,i)**2)
	END DO

	!potential energy
	DO i=1,np
		DO j=i+1,np
			Vpot(1,1) = Vpot(1,1) + (-G*m(i)*m(j)/SQRT((x(1,i)-x(1,j))**2+(y(1,i)-y(1,j))**2))
		END DO
	END DO

	!Basic stormer-verlet integration, first timestep
	!Accelerations(0)

	!Accel time
	DO i=1,6
		DO j=1,6
			IF (j /= i) THEN !Vent lige, det skal da ikke være ax(i) = ax(i) + Noget andet?
				ax(1,i) = ax(1,i)+(-G*m(j)/((x(1,i)-x(1,j))**2+(y(1,i)-y(1,j))**2)**(1.5))*(x(1,i)-x(1,j))
				ay(1,i) = ay(1,i)+(-G*m(j)/((x(1,i)-x(1,j))**2+(y(1,i)-y(1,j))**2)**(1.5))*(y(1,i)-y(1,j))
			END IF
		END DO
	END DO

	!positions step 1
	x(2,:) = x(1,:) + u(1,:)*dt + 0.5*ax(1,:)*dt**2
	y(2,:) = y(1,:) + v(1,:)*dt + 0.5*ay(1,:)*dt**2 


	DO n=2, nt !nt or nt+1

		!Accel time
		CALL accel(m,x(n,:),y(n,:),ax(n,:),ay(n,:),G,Np)
		
		!Velocity time
		CALL velocity(ax(n,:),ay(n,:),u(n-1,:),v(n-1,:),u(n,:),v(n,:),dt,np)

		!Position time
		CALL position(ax(n,:),ay(n,:),x(n,:),y(n,:),x(n-1,:),y(n-1,:),x(n+1,:),y(n+1,:),dt,Np)


		!Kinetic energy
		DO i=1,np
			Tkin(n,1) = Tkin(n,1) + 0.5*m(i)*(u(n,i)**2+v(n,i)**2)
		END DO
		
		!Potential energy
		DO i=1,np
			DO j=i+1,np
				Vpot(n,1) = Vpot(n,1) + (-G*m(i)*m(j)/SQRT((x(n,i)-x(n,j))**2+(y(n,i)-y(n,j))**2))
			END DO
		END DO							


		!time average kinetic energy
		Tavg(n,1) = SUM(Tkin(:n,1))/n

		!time average potential energy
		Vavg(n,1) = SUM(Vpot(:n,1))/n

	ENDDO

END SUBROUTINE mainr

