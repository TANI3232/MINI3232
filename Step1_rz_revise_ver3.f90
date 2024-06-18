!      *************************************************
!      *  Bio Heat Transfer - cyclic heating STEP1 -   *
!      *                            2-D version        *
!      *                                               *
!      *                   Programed by Tomoki TANI    *
!      *                             Feb.22.2024       *
!      *************************************************


Program cyclic_heating_Step1
    Implicit None
    integer :: i, j, n
    
! Constants
    real(8), parameter :: RMIN = 2e-3, RMAX = 50e-3, ZMAX = 100e-3, pi = 3.1415926536             ! Size of domein [mm]
    real(8), parameter :: alpha_0 = 1.374E-7, lambda_0 = 0.50         ! ROW = 1050, CP = 3465, lambda = 0.50, (Thermophysical properties)
    real(8), parameter :: h = 0.010                                   ! Heat transfer coefficient
    real(8), parameter :: Ts = 37.0, Tf = 20.0                        ! Human core temperature
!    real(8), parameter :: Heatsource = 4*pi/5                             ! Heat flux
    real(8), parameter :: Heatsource = 0.4
    real(8), parameter :: probearea = 2*pi*2E-03*0.02
    real(8), allocatable, dimension(:,:) :: T
    real(8), allocatable, dimension(:,:) :: S, Sz, deltaT, a, b, c, d, r, z
    real(8), allocatable, dimension(:) :: a_tmp,b_tmp,c_tmp,d_tmp,a_tmp2,b_tmp2,c_tmp2,d_tmp2
    real(8) :: dr, dz, dt, aa, bb, q, Area, courant_number
    real(8) :: dT_dr2, dT_dr, dT_dz2
    character(len=30) :: filename
    
    ! Get user inputs
    ! Enter the number of grid points : "
    integer :: nr = 51, nz = 101, nt = 60000, z1, z2
    ! Enter the spatial step size : "
    dr =  RMAX / (nr-1)
    dz =  ZMAX / (nz-1)
    ! Enter the time step size : "
    dt =  1.0e-3
    z1 = int(nz*0.4)
    z2 = int(nz*0.6)
!    Area = 2*pi*RMIN*20    ! 20 = Heat length 
    q=Heatsource/Area
!    q=5000
    courant_number = alpha_0*dt/dr**2
    if (courant_number >= 0.1) then
        print *,  "Warning: Courant number exceeds stability limit."
    end if
    

	allocate(T(nz,nr), S(nz,nr), Sz(nz,nr), deltaT(nz,nr), a(nz,nr), b(nz,nr), c(nz,nr), d(nz,nr), r(nz,nr), z(nz,nr))
	allocate(a_tmp(nr),b_tmp(nr),c_tmp(nr),d_tmp(nr),a_tmp2(nz),b_tmp2(nz),c_tmp2(nz),d_tmp2(nz))

! Initialization of r position
	Do j = 1, nr
	    r(:,j) = ((j-1) * dr + RMIN)
	End Do

! Set initialization condition        
    T = Ts
!    T(50,25)=100

! Time-stepping loop
    Do n= 1,nt

        ! Solve S(i,j)
        Do j = 2, nr-1
            Do i = 2, nz-1
                dT_dr2 = (T(i,j-1)-2*T(i,j)+T(i,j+1))/dr**2
                dT_dr  = (T(i,j+1)-T(i,j-1))/2/dr
                dT_dz2 = (T(i-1,j)-2*T(i,j)+T(i+1,j))/dz**2
                S(i,j) = alpha_0*dt * (dT_dr2 + (1/r(i,j))*dT_dr + dT_dz2)
            End Do
        End Do

!***********Boundary conditions of r direction***********
        ! Substitute
        aa = (alpha_0*dt)/(dr**2)

        ! left side 0-40
        a( 1:z1 , 1 )   = 0
        b( 1:z1 , 1 )   = -aa - dr*aa/r(1:z1,1)
        c( 1:z1 , 1 )   = -1 + aa + dr*aa/r(1:z1,1)
        d( 1:z1 , 1 )   = -S(1:z1,2)
        
        ! left side 41-60
        a( z1+1:z2 , 1 )  = 0
        b( z1+1:z2 , 1 )  = -aa - dr*aa/r(z1+1:z2,1)
        c( z1+1:z2 , 1 )  = -1  +aa + dr*aa/r(z1+1:z2,1)
        d( z1+1:z2 , 1 )  = dr*q*aa*dt* ((-dr/r(z1+1:z2,1)/2) -1) /lambda_0-S(z1+1:z2,2)
        
        ! left side 61-100
        a( z2+1:nz , 1 ) = 0
        b( z2+1:nz , 1 ) = -aa - dr*aa/r(z2+1:nz,1)
        c( z2+1:nz , 1 ) = -1 + aa + dr*aa/r(z2+1:nz,1)
        d( z2+1:nz , 1 ) = -S(z2+1:nz,2)

        ! right side 0-101      
        a( : , nr )     = -1 + aa + dr*aa/r(:,nr)
        b( : , nr )     = -aa + dr*aa/r(:,nr)
        c( : , nr )     = 0
        d( : , nr )     = -S(:,nr-1)

        ! i direction i=(2:N-1) 0-40
        a(:, 2:nr-1)     = -aa/2 + aa/4*dr/r(:,2:nr-1)
        b(:, 2:nr-1)     = 1 + aa
        c(:, 2:nr-1)     = -aa/2 - aa/4*dr/r(:,2:nr-1)
        d(:, 2:nr-1)     = S(:,2:nr-1)

!**************End of Boundary conditions of r direction*****************

! Solve the tridiagonal system using the Thomas Algorithm
        Do i = 1,nz
            d_tmp = d(i,:)
            a_tmp = a(i,:)
            b_tmp = b(i,:)
            c_tmp = c(i,:)
            Call solve_tridiag(a_tmp,b_tmp,c_tmp,d_tmp,Sz(i,:),Nr)
        end do

!**************Boundary conditions of z direction*****************
        ! Substitute
        bb = (alpha_0*dt)/(dz**2)

        ! Boundary conditions of z direction 
        a( 1 , : ) = 0
        b( 1 , : ) = -2*bb-dz*h*bb/lambda_0
        c( 1 , : ) = 3 + bb + 2*(1+bb)*dz*h/lambda_0
        d( 1 , : ) = (bb*dz*h*dt*(T(1,:)-Tf))/lambda_0+Sz(2,:)*(3+(2*dz*h/lambda_0))
        
        a( nz , : ) = 0
        b( nz , : ) = 1
        c( nz , : ) = 0
        d( nz , : ) = 0
        
        ! j direction j=(2:N-1) 
        a( 2:Nz-1, : ) = -bb/2
        b( 2:Nz-1, : ) = 1+bb
        c( 2:Nz-1, : ) = -bb/2
        d( 2:Nz-1, : ) = Sz(2:Nz-1,:)

!**************End of Boundary conditions of z direction*****************

! Solve the tridiagonal system using the Thomas Algorithm
        Do j = 1,nr
            d_tmp2 = d(:,j)
            a_tmp2 = a(:,j)
            b_tmp2 = b(:,j)
            c_tmp2 = c(:,j)
            Call solve_tridiag(a_tmp2,b_tmp2,c_tmp2,d_tmp2,deltaT(:,j),Nz)
        End do

! Update solution
        T(:,:) = T(:,:) + deltaT(:,:)
        T(101,:) = Ts
        write(*,*) n

!        If (MOD (n,1000).eq. 1) then    
        If (n .eq. 1) then    
            write(filename, '(A, I11.11, A)') 'Temperature_', n, '.dat'
    
            open(unit=11, file=filename, status="replace", action="write")
                Do i = 1, nz
                    write(11,*) (T(i,j), j=1,nr)
                End do
            close(11)
        End If


    End Do

    write(filename, '(A, I11.11, A)') 'Temperature_', nt, '.dat'
    open(unit=11, file=filename, status="replace", action="write")
        Do i = 1, nz
            write(11,*) (T(i,j), j=1,nr)
        End do
   close(11)

contains 

!****************************************************************    
    Subroutine solve_tridiag(a,b,c,d,x,n)
        implicit none
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: a,b,c,d
        real(8), dimension(n), intent(out) :: x
        real(8), dimension(n) :: s, e
        real(8) :: v
        integer :: i 
    
    ! Forward elimination
        s(1) = -c(1) / b(1)
        e(1) = d(1) / b(1)
        do i = 2, n
            v    = b(i) + a(i)*s(i-1)
            s(i) = -c(i)/v
            e(i) = (d(i)-a(i)*e(i-1))/v
        end do
    
    !  Back subsititution
        x(n) = e(n)
        do i = n-1, 1, -1
            x(i) = e(i) + s(i)*x(i+1)
        end do
    End subroutine solve_tridiag

End Program
