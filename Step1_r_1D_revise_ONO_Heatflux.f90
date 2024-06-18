!      *************************************************
!      *  Bio Heat Transfer - cyclic heating STEP1 -   *
!      *                            1-D version        *
!      *                                               *
!      *                   Programed by Tomoki TANI    *
!      *                             Jun.05.2024       *
!      *************************************************

Program cyclic_heating_Step1
    Implicit None
    integer :: j, n

    ! Constants
    real(8), parameter :: RMIN = 2E-03 , RMAX = 50E-03                ! Size of domein [m]
    real(8), parameter :: alpha_0 = 1.374E-7, lambda_0 = 0.50         ! ROW = 1050, CP = 3465, lambda = 0.50, (Thermophysical properties)
    real(8), parameter :: Ts = 37.0                                   ! Human core temperature                                      
    real(8), parameter :: Heatsource = 0.4
    real(8), parameter :: probearea = 2E-03**2*pi
    real(8), allocatable, dimension(:) :: T
    real(8), allocatable, dimension(:) :: S, deltaT, a, b, c, d, r
    real(8) :: dr, dt, aa, q, courant_number
    real(8) :: dT_dr2, dT_dr
    character(len=30) :: filename

    ! Get user inputs
    ! Enter the number of grid points : "
    integer :: nr = 501, nt = 60000
    ! Enter the spatial step size : "
    dr =  RMAX / (nr-1)
    ! Enter the time step size : "
    dt =  1e-3

    ! Heat flux
    q = Heatsource / probearea

    courant_number = alpha_0*dt/dr**2
    if (courant_number >= 0.1) then
        print *,  "Warning: Courant number exceeds stability limit."
    end if

    
    allocate(T(nr), S(nr), deltaT(nr), a(nr), b(nr), c(nr), d(nr), r(nr))
    
    ! Initialization of r position
	Do j = 1, nr
	    r(j) = ((j-1) * dr + RMIN)
	End Do
    
    ! Set initialization condition
    T = Ts

    ! Time-stepping loop
    Do n = 1, nt

        ! Solve S(i,j)
        Do j = 2, nr-1
                dT_dr2 = (T(j-1)-2*T(j)+T(j+1))/dr**2
                dT_dr  = (T(j+1)-T(j-1))/2/dr
                S(j) = alpha_0*dt/2 * (dT_dr2 + (1/r(j))*dT_dr)
        End Do


!***********Boundary conditions of r direction***********

        ! Substitute
        aa = (alpha_0*dt)/(dr**2)
        
        ! left side
        a(1) = 0
        b(1) = -aa-aa*dr/r(1)
        c(1) = -1+aa+aa*dr/r(1)
        d(1) = dr/lambda_0*q*dt*(-aa-aa*dr/2/r(1))-S(2)

        ! light side
        a(nr) = -3-aa-aa*dr/r(nr)
        b(nr) = aa+aa*dr/r(nr)
        c(nr) = 0
        d(nr) = -3*S(nr-1)


        ! r direction i=(2:N-1) 
        a(2:nr-1) = -aa/2 + aa/4*dr/r(2:nr-1)
        b(2:nr-1) = 1 + aa
        c(2:nr-1) = -aa/2 - aa/4*dr/r(2:nr-1)
        d(2:nr-1) = S(2:nr-1)

!**************End of Boundary conditions of r direction*****************

        ! Solve the tridiagonal system using the Thomas Algorithm
        Call solve_tridiag(a,b,c,d,deltaT,Nr)

        ! Update solution
        T(:) = T(:) + deltaT(:)
        write(*,*) n

        If (n .eq. 1) then    
            write(filename, '(A, I11.11, A)') 'Temperature_', n, '.dat'
    
            open(unit=11, file=filename, status="replace", action="write")
                Do j = 1, nr
                    write(11,*) T(j)
                End do
            close(11)
        End If


    End Do

    write(filename, '(A, I11.11, A)') 'Temperature_', nt, '.dat'
    open(unit=11, file=filename, status="replace", action="write")
        Do j = 1, nr
            write(11,*) T(j)
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