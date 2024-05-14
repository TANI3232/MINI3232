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
        real(8), parameter :: alpha_0 = 1.374E-7, lambda_0 = 0.50         !ROW = 1050, CP = 3465, lambda = 0.50, (Thermophysical properties)
        real(8), parameter :: h = 0.010                                   ! Heat transfer coefficient
        real(8), parameter :: Ts = 37.0, Tf = 20.0                        ! Human core temperature
        real(8), parameter :: q=2.5                                       ! Heat flux
        real(8), allocatable, dimension(:,:,:) :: T
        real(8), allocatable, dimension(:,:) :: S, Sz, deltaT, a, b, c, d, r, z
        real(8), allocatable, dimension(:) :: dd,ddd,aaa,bbb,ccc,aaaa,bbbb,cccc
        real(8) :: dr, dz, dt,aa,bb


        ! Get user inputs
        ! Enter the number og grid points : "
        integer :: nr = 51, nz = 101, nt = 10
        ! Enter the spatial step size : "
        dr =  RMAX / (nr-1)*1.0E-03
        dz =  ZMAX / (nz-1)*1.0E-03
        ! Enter the time step size : "
        dt =  1.0E-6   

        allocate(T(nr,nz,nt+1), S(nr,nz), Sz(nr,nz), deltaT(nr,nz), a(nr,nz), b(nr,nz), c(nr,nz), d(nr,nz), r(nr,nz), z(nr,nz),  dd(nr),ddd(nz),aaa(nr),bbb(nr),ccc(nr),aaaa(nz),bbbb(nz),cccc(nz) )
        


! Initialization of r position
        Do i = 1, nr
          Do j = 1, nz
           r(i,:) = ((i-1) * dr + 0.5*dr + 1)*1.0E-03
          End Do
        End Do


! Set initialization condition        
           T(:,:,1) = Ts
!           T(25,50,1) = 2*Ts

! Time-stepping loop
    Do n= 1,nt

! Solve S(i,j)


          Do j = 2, nz-1
            Do i = 2, nr-1
               S(i,j) = (alpha_0*dt/2)*((T(i-1,j,n)-2*T(i,j,n)+T(i+1,j,n))/dr**2+(1/r(i,j))*((T(i+1,j,n)-T(i-1,j,n))/2/dr)+(T(i,j-1,n)-2*T(i,j,n)+T(i,j+1,n))/dz**2)
            End Do
          End Do


! Boundary conditions of r direction

! Substitute
        aa = (alpha_0*dt)/(dr**2)

! left side 0-40
          a( 1 , 1:40 ) = 0
          b( 1 , 1:40 ) = -aa+dr*aa/r(1,1:40)
          c( 1 , 1:40 ) = -1+aa-dr*aa/r(1,1:40)
          d( 1 , 1:40 ) = -S(2,1:40)


! left side 41-60
          a( 1 , 41:60 ) = 0
          b( 1 , 41:60 ) = -aa+dr*aa/r(1,41:60)
          c( 1 , 41:60 ) = -1+aa-dr*aa/r(1,41:60)
          d( 1 , 41:60 ) = dr*q*aa*((1/2)*(dr/r(1,:))-1)/lambda_0-S(2,41:60)


! left side 61-100
          a( 1 , 61:101 ) = 0
          b( 1 , 61:101 ) = -aa+dr*aa/r(1,61:101)
          c( 1 , 61:101 ) = -1+aa-dr*aa/r(1,61:101)
          d( 1 , 61:101 ) = -S(2,61:101)


! right side 0-101      
          a( Nr , : ) = -1+aa+dr*aa/r(Nr,:)
          b( Nr , : ) = -aa-dr*aa/r(Nr,:)
          c( Nr , : ) = 0
          d( Nr , : ) = -S(nr-1,:)


! i direction i=(2:N-1) 0-40
          a(2:Nr-1, :) = -aa/2-aa/4*dr/r(2:Nr-1,2:Nz-1)
          b(2:Nr-1, :) = 1+aa
          c(2:Nr-1, :) = -aa/2+aa/4*dr/r(2:Nr-1,2:Nz-1)
          d(2:Nr-1, :) = S(2:Nr-1,2:Nz-1)


! Solve the tridiagonal system using the Thomas Algorithm
          Do j = 1,nz
            dd = d(:,j)
            aaa = a(:,j)
            bbb = b(:,j)
            ccc = c(:,j)
          Call solve_tridiag_r(aaa,bbb,ccc,dd,Sz(:,j),Nr)
          end do


          open(unit=12, file="Sz.csv", status="replace", action="write")
          do j = 1, nz
            write(12,'(*(g0,:,","))') (Sz(i,j), i=1,nr)
          end do
            close(12)

!********************************************************

! Substitute
          bb= (alpha_0*dt)/(dz**2)


! Boundary conditions of z direction 

                a( : , 1 ) = 0
                b( : , 1 ) = bb
                c( : , 1 ) = 1-bb
                d( : , 1 ) = (-2*bb*dz*h*(T(:,1,n)-Tf))/lambda_0+Sz(:,2)

                a( : , 101) = 0
                b( : , 101) = 1
                c( : , 101) = 0
                d( : , 101) = 0


! j direction j=(2:N-1) 
                a( :, 2:Nz-1) = -bb/2
                b( :, 2:Nz-1) = 1+bb
                c( :, 2:Nz-1) = -bb/2
                d( :, 2:Nz-1) = Sz(:,2:Nz-1)


! Solve the tridiagonal system using the Thomas Algorithm
          do i = 1,nr
            ddd = d(i,:)
            aaaa = a(i,:)
            bbbb = b(i,:)
            cccc = c(i,:)
          Call solve_tridiag_z(aaaa,bbbb,cccc,ddd,deltaT(i,:),Nz)
          end do
                     


! Update solution
          T(:,:,n+1) = T(:,:,n) + deltaT
          T(:,101,n+1) = Ts
          write(*,*) n
  

          ! Output the results
          open(unit=10, file="deltaT.csv", status="replace", action="write")
          do j = 1, nz
            write(10,'(*(g0,:,","))') (deltaT(i,j), i=1,nr)
          end do
            close(10)
    End Do



    open(unit=11, file="temperature.csv", status="replace", action="write")
    do j = 1, nz
!      write(11,'(*(g0,:,","))') (T(i,j,n), i=1,nr)
      write(11,*) (T(i,j,n), i=1,nr)
    end do
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
               v = b(i) + a(i)*s(i-1)
               s(i) = -c(i)/v
               e(i) = (d(i)-a(i)*e(i-1))/v
             end do

!  Back subsititution
            p(n) = e(n)
             do i = n-1, 1, -1
              p(i) = e(i) + s(i)*p(i+1)
             end do
!             write(*, *) b
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
               w = b(j) + a(j)*f(j-1)
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
   