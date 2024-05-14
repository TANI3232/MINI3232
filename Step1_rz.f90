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
    real(8), parameter :: RMAX = 50, ZMAX = 100                       ! Size of domein
    real(8), parameter :: alpha_0 = 1.374E-7, lambda_0 = 0.50         ! ROW = 1050, CP = 3465, lambda = 0.50, (Thermophysical properties)
    real(8), parameter :: h = 0.010                                   ! Heat transfer coefficient
    real(8), parameter :: Ts = 37.0, Tf = 20.0                        ! Human core temperature
    real(8), parameter :: q=2.5                                       ! Heat flux
    real(8), allocatable, dimension(:,:,:) :: T
    real(8), allocatable, dimension(:,:) :: S, Sz, deltaT, a, b, c, d, r, z
    real(8), allocatable, dimension(:) :: a_tmp,b_tmp,c_tmp,d_tmp,a_tmp2,b_tmp2,c_tmp2,d_tmp2
    real(8) :: dr, dz, dt,aa,bb
    real(8) :: dT_r2, dT_r, dT_z2
    
    ! Get user inputs
    ! Enter the number og grid points : "
    integer :: nr = 51, nz = 101, nt = 10
    ! Enter the spatial step size : "
    dr =  RMAX / (nr-1)*1.0E-03
    dz =  ZMAX / (nz-1)*1.0E-03
    ! Enter the time step size : "
    dt =  1.0E-6   

	allocate(T(nz,nr,nt+1), S(nz,nr), Sz(nz,nr), deltaT(nz,nr), a(nz,nr), b(nz,nr), c(nz,nr), d(nz,nr), r(nz,nr), z(nz,nr))
	allocate(a_tmp(nr),b_tmp(nr),c_tmp(nr),d_tmp(nr),a_tmp2(nz),b_tmp2(nz),c_tmp2(nz),d_tmp2(nz) )
	

! Initialization of r position
	Do j = 1, nr
	    r(:,j) = ((j-1) * dr + 0.5*dr)*1.0E-03
	End Do

! Set initialization condition        
    T(:,:,1) = Ts

! Time-stepping loop
    Do n= 1,nt

        ! Solve S(i,j)
        Do j = 2, nr-1
            Do i = 2, nz-1
                dT_r2 = T(i,j-1,n)-2*T(i,j,n)+T(i,j+1,n)
                dT_r  = T(i,j+1,n)-T(i,j-1,n)
                dT_z2 = T(i-1,j,n)-2*T(i,j,n)+T(i+1,j,n)
                S(i,j) = (alpha_0*dt/2) * (dT_r2/dr**2 + (1/r(i,j))*(dT_r/2/dr) + dT_z2/dz**2)
            End Do
        End Do

!***********Boundary conditions of r direction***********
        ! Substitute
        aa = (alpha_0*dt)/(dr**2)

        ! left side 0-40
        a( 1:40 , 1 )   = 0
        b( 1:40 , 1 )   = -aa + dr*aa/r(1:40,1)
        c( 1:40 , 1 )   = -1 + aa - dr*aa/r(1:40,1)
        d( 1:40 , 1 )   = -S(1:40,2)
        
        ! left side 41-60
        a( 41:60 , 1 )  = 0
        b( 41:60 , 1 )  = -aa + dr*aa/r(41:60,1)
        c( 41:60 , 1 )  = -1  +aa - dr*aa/r(41:60,1)
        d( 41:60 , 1 )  = dr*q*aa* ((1/2)* (dr/r(41:60,1)) -1) /lambda_0-S(41:60,2)
        
        ! left side 61-100
        a( 61:101 , 1 ) = 0
        b( 61:101 , 1 ) = -aa + dr*aa/r(61:101,1)
        c( 61:101 , 1 ) = -1 + aa - dr*aa/r(61:101,1)
        d( 61:101 , 1 ) = -S(61:101,2)

        ! right side 0-101      
        a( : , 51 )     = -1+aa+dr*aa/r(:,51)
        b( : , 51 )     = -aa-dr*aa/r(:,51)
        c( : , 51 )     = 0
        d( : , 51 )     = -S(:,50)

        ! i direction i=(2:N-1) 0-40
        a(:, 2:50)     = -aa/2  -aa/4*dr/r(:,2:50)
        b(:, 2:50)     = 1 + aa
        c(:, 2:50)     = -aa/2 + aa/4*dr/r(:,2:50)
        d(:, 2:50)     = S(:,2:50)

!**************End of Boundary conditions of r direction*****************

! Solve the tridiagonal system using the Thomas Algorithm
        Do i = 1,nz
            d_tmp = d(i,:)
            a_tmp = a(i,:)
            b_tmp = b(i,:)
            c_tmp = c(i,:)
            Call solve_tridiag_r(a_tmp,b_tmp,c_tmp,d_tmp,Sz(i,:),Nr)
        end do

        open(unit=12, file="Sz_rz.csv", status="replace", action="write")
            do i = 1, nz
                write(12,'(*(g0,:,","))') (Sz(i,j), j = 1, nr)
            end do
        close(12)


!**************Boundary conditions of z direction*****************
        ! Substitute
        bb = (alpha_0*dt)/(dz**2)

        ! Boundary conditions of z direction 
        a( 1 , : ) = 0
        b( 1 , : ) = bb
        c( 1 , : ) = 1-bb
        d( 1 , : ) = (-2*bb*dz*h*(T(1,:,n)-Tf))/lambda_0+Sz(2,:)
        
        a( 101 , : ) = 0
        b( 101 , : ) = 1
        c( 101 , : ) = 0
        d( 101 , : ) = 0
        
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
            Call solve_tridiag_z(a_tmp2,b_tmp2,c_tmp2,d_tmp2,deltaT(:,j),Nz)
        End do

! Update solution
        T(:,:,n+1) = T(:,:,n) + deltaT(:,:)
        T(101,:,n+1) = Ts
        write(*,*) n
    
        ! Output the results
        open(unit=10, file="deltaT_rz.csv", status="replace", action="write")
            Do i = 1, nz
                write(10,'(*(g0,:,","))') (deltaT(i,j), j=1,nr)
            End do
        close(10)
    End Do

    open(unit=11, file="Temperature.csv", status="replace", action="write")
        Do i = 1, nz
!            write(11,'(f8.3(g0,:,","))') (T(i,j,n), j=1,nr)
            write(11,*) (T(i,j,n), j=1,nr)
        End do
    close(11)

contains 

!****************************************************************    
    Subroutine solve_tridiag_r(a,b,c,d,p,n)
        implicit none
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: a,b,c,d
        real(8), dimension(n), intent(out) :: p
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
        p(n) = e(n)
        do i = n-1, 1, -1
            p(i) = e(i) + s(i)*p(i+1)
        end do
    End subroutine solve_tridiag_r

!****************************************************************     
    Subroutine solve_tridiag_z(a,b,c,d,temp,n)
        Implicit none
        integer, intent(in) :: n
        real(8), dimension(n), intent(in) :: a,b,c,d
        real(8), dimension(n), intent(out) :: temp
        real(8), dimension(n) :: f, g
        real(8) :: w
        integer :: j 
    
    ! Forward elimination
        f(1) = -c(1) / b(1)
        g(1) = d(1) / b(1)
        do j = 2, n
            w    = b(j) + a(j)*f(j-1)
            f(j) = -c(j)/w
            g(j) = (d(j)-a(j)*g(j-1))/w
        end do
    
    !  Back subsititution
        temp(n) = g(n)
        do j = n-1, 1, -1
            temp(j) = g(j) + f(j)*temp(j+1)
        end do
    
    End subroutine solve_tridiag_z       

End Program
