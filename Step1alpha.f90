!      *************************************************
!      *  Bio Heat Transfer - cyclic heating STEP1 -            *
!      *                                               2-D version        *
!      *                                                                      *
!      *                           Programed by Tomoki TANI        *
!      *                                                Feb.22.2024       *
!      *************************************************


  Program cyclic_heating_Step1
        Implicit None
        integer :: i, j, n


! Constants
        real(8), parameter :: RMAX = 50, ZMAX = 100                       ! Size of domein
        real(8), parameter :: alpha_0 = 1.374E-7, lambda_0 = 0.50         !ROW = 1050, CP = 3465, lambda = 0.50, (Thermophysical properties)
        real(8), parameter :: h = 10                                      ! Heat transfer coefficient
        real(8), parameter :: Ts = 37                                     ! Human core temperature
        real(8), parameter :: q=5000                                      ! Heat flux
        real(8), allocatable, dimension(:,:) :: T, T_new, T_prev, lambda_T, S, Sz, alpha_T, deltaT, a, b, c, d, r, z, aa, bb


! Variables
        integer :: Nr = 31, Nz = 31, Nt = 101           ! Nimber of grid points
        real(8) :: dr, dz, dt, NEND
        allocate(T(Nr,Nz), T_new(Nr,Nz), T_prev(Nr,Nz), lambda_T(Nr,Nz), S(Nr,Nz), Sz(Nr,Nz), alpha_T(Nr,Nz), deltaT(Nr,Nz), a(Nr,Nz), b(Nr,Nz), c(Nr,Nz), d(Nr,Nz), r(Nr,Nz), z(Nr,Nz), aa(Nr,Nz), bb(Nr,Nz) )

! Discretization
        dr =  RMAX / (Nr-1)
        dz =  ZMAX / (Nz-1)
        dt =  1.0E-3
        

! Initialization of r position
        Do i = 1, Nr
          r(i,:) = (i-1) * dr
        End Do
        
        Do j = 1, Nz
          Do i = 1, Nr
           T(i,j) = Ts
          End Do
        End Do


 ! Time evolution
      Do n = 1, Nt
           T(i,j) = Ts
          write(*,*) n


!*********************************************************


! Substitute
          aa(i,j) = (alpha_0*dt)/(dr**2)

! Solve S(i,j)
          Do j = 1, Nz
            Do i = 1, Nr
              S(i,j) = (alpha_0(i,j)*dt/2)*((T(i-1,j)-2*T(i,j)+T(i+1,j))/dr**2+(1/r(i,j))*((T(i+1,j)-T(i-1,j))/2*dr)+(T(i,j-1)-2*T(i,j)+T(i,j+1))/dz**2)
            End Do
          End Do


! Boundary conditions of r direction

! left side 0-40
          a( 1 , 1:40 ) = 0
          b( 1 , 1:40 ) = -aa(i,j)+dr*aa(i,j)/r(i,:)
          c( 1 , 1:40 ) = -1+aa(i,j)-dr*aa(i,j)/r(i,:)
          d( 1 , 1:40 ) = -S(i,:)
       
! right side 0-40      
          a( Nr , 1:40 ) = -1+aa(i,j)+dr*aa(i,j)/r(i,:)
          b( Nr , 1:40 ) = -aa(i,j)-dr*aa(i,j)/r(i,:)
          c( Nr , 1:40 ) = 0
          d( Nr , 1:40 ) = -S(i,:)


! left side 41-60
          a( 1 , 41:60 ) = 0
          b( 1 , 41:60 ) = -aa(i,j)+dr*aa(i,j)/r(i,:)
          c( 1 , 41:60 ) = -1+aa(i,j)-dr*aa(i,j)/r(i,:)
          d( 1 , 41:60 ) = dr*q*aa(i,j)*((1/2)*(dr/r(i,:))-1)/lambda_0(i,j)-S(i,:)
       
! right side 41-60       
          a( Nr , 41:60 ) = -1+aa(i,j)+dr*aa(i,j)/r(i,:)
          b( Nr , 41:60 ) = -aa(i,j)-dr*aa(i,j)/r(i,:)
          c( Nr , 41:60 ) = 0
          d( Nr , 41:60 ) = -S(i,:)

! left side 61-100
          a( 1 , 61:100 ) = 0
          b( 1 , 61:100 ) = -aa(i,j)+dr*aa(i,j)/r(i,:)
          c( 1 , 61:100 ) = -1+aa(i,j)-dr*aa(i,j)/r(i,:)
          d( 1 , 61:100 ) = -S(i,j)

! right side 61-100       
          a( Nr , 61:100 ) = -1+aa(i,j)+dr*aa(i,j)/r(i,:)
          b( Nr , 61:100 ) = -aa(i,j)-dr*aa(i,j)/r(i,:)
          c( Nr , 61:100 ) = 0
          d( Nr , 61:100 ) = -S(i,j)
       
! i direction i=(2:N-1) 0-100
          a(2:Nr-1, j) = -aa(i,j)/4-aa(i,j)/4*dr/r(i,:)
          b(2:Nr-1, j) = 1+aa(i,j)/2
          c(2:Nr-1, j) = -aa(i,j)/4+aa(i,j)/4*dr/r(i,:)
          d(2:Nr-1, j) = S(i,j)

! Solve the tridiagonal system using the Thomas Algorithm
          Call solve_tridiag_r(a,b,c,d,Sz,Nr)

!********************************************************

! Substitute
          bb(i,j) = (alpha_0(i,j)*dt)/(dz**2)

! Boundary conditions of z direction 

          a( i , 1 ) = 0
          b( i , 1 ) = bb(i,j)+(bb(i,j)*dz*h)/lambda_0(i,j)
          c( i , 1 ) = 1-bb(i,j)
          d( i , 1 ) = (2*bb(i,j)*dz*h*Ts)/lambda_0(i,j)+Sz(i,j)
       
          a( i , Nz) = 0
          b( i , Nz) = 1
          c( i , Nz) = 0
          d( i , Nz) = Ts
       
! Solve the tridiagonal system using the Thomas Algorithm
          Call solve_tridiag_z(a,b,c,d,deltaT,Nz)

! Solve T(n+1)
          Do j = 1, Nz
            Do i = 1, Nr
             T(i,j) = deltaT(i,j) + T(i,j)
            End Do
           End Do
    
      End Do

! Output the results
!      open(unit=300,file="data N=301.csv") 
!      Do j = 1, Nz
!        Do i = 1, Nr
!         r(i,:) = (i-1) * dr
!         z(:,j) = (j-1) * dz
!         write(*,*) r(i,:), z(:,j), T(i,j)
!        End Do
!      End Do
!      Do i = 1, Nr
!         write(300,*) (T(i,j),",", j=1, Nz)  
!      End Do

! Output the results
open(unit=300, file="data_alpha.csv", status='replace', action='write') 
Do j = 1, Nz
    write(300, '(F8.3)', advance='no') T(1,j)
    Do i = 2, Nr
        write(300, '(A,F8.3)', advance='no') ",", T(i,j)
    End Do
    write(300, *)
End Do


    contains 

    
!****************************************************************    

           Subroutine solve_tridiag_r(a,b,c,d,sz,n)
            implicit none
            integer, intent(in) :: n
            real(8), dimension(n), intent(in) :: a,b,c,d
            real(8), dimension(n), intent(out) :: sz
            real(8), dimension(n) :: s, e
            real(8) :: v
            integer :: i 

! Forward elimination
            s(1) = c(1) / b(1)
            e(1) = d(1) / b(1)
             do i = 2, n
               v = b(i) - a(i)*s(i-1)
               s(i) = c(i)/v
               e(i) = (d(i)-a(i)*e(i-1))/v
             end do

!  Back subsititution
            sz(n) = e(n)
             do i = n-1, 1, -1
              sz(i) = e(i) - s(i)*sz(i+1)
             end do
           End subroutine solve_tridiag_r
    
    
!****************************************************************     


           Subroutine solve_tridiag_z(a,b,c,d,deltaT,n)
            Implicit none
            integer, intent(in) :: n
            real(8), dimension(n), intent(in) :: a,b,c,d
            real(8), dimension(n), intent(out) :: deltaT
            real(8), dimension(n) :: s, e
            real(8) :: w
            integer :: j 

! Forward elimination
             s(1) = c(1) / b(1)
             e(1) = d(1) / b(1)
              do j = 2, n
               w = b(j) - a(j)*s(j-1)
               s(j) = c(j)/w
               e(j) = (d(j)-a(j)*e(j-1))/w
              end do

!  Back subsititution
            deltaT(n) = e(n)
              do j = n-1, 1, -1
                deltaT(j) = e(j) - s(j)*deltaT(j+1)
              end do
            
           End subroutine solve_tridiag_z       
        
      


  End Program
   