!This code solve the TD model equations for laminar flame in 2 dimensional cartesian coordinates
  !Additionally it includes selective frequency damping algorithm to dampen the single frequuency oscillations in order to obtain the stationaly solution
! The damping can be controlled using the paramters chi and delta
!the agorithm is described in detail in B.E. Jordi et.al. Physics of Fluid 2014 paper
module data_struct
!Defines the main data struxture used in the program.
  implicit none
  type species
     real*8  :: val, derx, der2x, der2y
!val --> Value at that point
!derx --> First x derivative
!der2x --> Second x derivative
!der2y --> Second y derivative
  end type species
      
  type TFO
     type(species) :: t,f,o
!Stores the temperature, fuel and oxidiser concentrations. Inherits the F data type to store the actual value and various derivatives at each point.
!T --> Temperature
!f --> Fuel
!o --> Oxidiser
  end type TFO
         
  type tfo1
     real*8 :: t,f,o
!Stores the temperature, fuel and oxidiser concentrations. Inherits the F data type to store the actual value and various derivatives at each point.
!T --> Temperature
!f --> Fuel
!o --> Oxidiser
  end type TFO1


  type Freal
     real  :: val, derx, der2x, der2y
            !val --> Value at that point
!derx --> First x derivative
!der2x --> Second x derivative
!der2y --> Second y derivative
  end type Freal
      
  type TFOreal
     type(Freal) :: t,f,o
!Stores the temperature, fuel and oxidiser concentrations. Inherits the F data type to store the actual value and various derivatives at each point.
!T --> Temperature
!f --> Fuel
!o --> Oxidiser
  end type TFOreal
end module data_struct
       
!module func
!contains

!end module func
module tfo1_addition
  use data_struct
  interface operator(+)
     module procedure add0,add1
  end interface
contains
  function add0(a,b)
    type(tfo1)  ::add0
    type(tfo1),intent(in)::a,b
    add0%t=a%t+b%t   !!!!addition for the components
    add0%f=a%f+b%f
    add0%o=a%o+b%o
  end function add0

  function add1(a,b) result(c)
    type(tfo1),dimension(:),intent(in)::a,b
   ! write(*,*)"size",size(a)
    type(tfo1),dimension(size(a))::c
    c%t=a%t+b%t
    c%f=a%f+b%f
    c%o=a%o+b%o
  end function add1
end module tfo1_addition

module tfo_operation
  use data_struct
  interface operator(+)
     module procedure add00,add11
  end interface
  
  interface operator(*)
     module procedure p0,p1,p3,p4 !,p3,p4
  end interface
contains
  function add00(a,b) result(c)
    type(tfo)::c
    type(tfo),intent(in)::a,b
    c%t%val=a%t%val+b%t%val   !!!!addition for the components
    c%f%val=a%f%val+b%f%val
    c%o%val=a%o%val+b%o%val
  end function add00
  
  function add11(a,b) result(c)
    type(tfo),dimension(:),intent(in)::a,b
    type(tfo),dimension(:),allocatable::c
 !    type(tfo),dimension(size(a))::c
    
    allocate(c(size(a)))
    c%t%val=a%t%val+b%t%val   !!!!addition for the components
    c%f%val=a%f%val+b%f%val
    c%o%val=a%o%val+b%o%val
    
 end function add11
 
  function p0(a,x) result(d) 
    type(tfo)  :: d
    type(tfo),intent(in)::x
    real ,intent(in)::a
    d%t%val=a*(x%t%val)
    d%f%val=a*(x%f%val)
    d%o%val=a*(x%o%val)

    end function p0

  function p1(a,x) result(d)
    type(tfo),dimension(:),intent(in)::x
    real,intent(in)::a
!    type(tfo),dimension(size(x))::d

    type(tfo),dimension(:),allocatable::d  
    allocate(d(size(x)))
    d%t%val=a*(x%t%val)
    d%f%val=a*(x%f%val)
    d%o%val=a*(x%o%val)
  end function p1  
    function p3(a,x) result(d) 
    type(tfo)  :: d
    type(tfo),intent(in)::x
    real*8 ,intent(in)::a
    d%t%val=a*(x%t%val)
    d%f%val=a*(x%f%val)
    d%o%val=a*(x%o%val)

    end function p3
        function p4(a,x) result(d) 
    type(tfo)  :: d
    type(tfo),intent(in)::x
    type(tfo1),intent(in)::a
    d%t%val=a%t*(x%t%val)
    d%f%val=a%f*(x%f%val)
    d%o%val=a%o*(x%o%val)

  end function p4
end module tfo_operation
  !addition of a tfo type  and tfo1 type
module tfo10_addition 
use data_struct
  interface operator(+)
     module procedure add01,add02
  end interface
contains

 function add01(a,b)
    type(tfo1)  ::add01
type(tfo),intent(in)::a
    type(tfo1), intent(in)::b
    add01%t=a%t%val+b%t   !!!!addition for the components
    add01%f=a%f%val+b%f
    add01%o=a%o%val+b%o
  end function add01
  
  function add02(a,b)
type(tfo),dimension(:),intent(in)::a 
    type(tfo1),dimension(:),intent(in)::b
    type(tfo1),dimension(size(a))::add02
    add02%t=a%t%val+b%t
    add02%f=a%f%val+b%f
    add02%o=a%o%val+b%o
  end function add02

end module tfo10_addition

!Start of the main program
program diff_flame
  use data_struct
  use omp_lib
  use tfo_operation  
  implicit none
  !     use func
!!!!!Problem specific variables. Not meant specifically for the program.!!!!
  !!Parameters to be specified
!Ambient Temperature
  real :: To = 288.15
  !Heat Release rate
  real :: Q
  !Cp
  real :: Cp
!Stoichiometric coeff for oxidiser
  real :: nu
  !Thermal diffusivity
  real :: Dth !That's for the sake of matalon idiot. I would prefer lambda
  !Lewis Numbers
  real :: LeF
  real :: LeO
  real :: LeF_inv, LeO_inv !It's the inverses that are always used. Better store them.
  !Fuel and Oxidiser velocity
  real :: Uo
  !Mass fractions in the inlet stream
  real :: YfO, YoO
      !Activation energy
  real :: Ea
  !Reaction rate parameter which decides the number of collisions per unit time
  real :: B
  !Gas constant
  real :: R_gas
  !The constant density assumed
  real :: rho 
  !Parameter which decides the inlet profile 
  real :: alphaparam

  !!Parameters that are to be calculated from the above parameters
  !Adiabatic Flame temperature
  real :: Ta
  !Equivalence ratio
  real :: phi
  !Zeldovich number
  real :: beta
!Heat release parameter
  real :: gamma
  !Damkohler number
  real :: Da
  
!!!!!Program specific variables!!!!
  
  !Main variable declaration
  type(TFO), dimension(:), allocatable :: arr_tfo ,arr_res,k1,arr_dt1  
  !!array to store filtered variable
  type(TFO), dimension(:), allocatable :: arr_tfo1  
  !Stores reaction rate at that point
  real*8, dimension(:), allocatable :: w 
  
  !Error, Norm
!  real*8, dimension(:), allocatable :: err
  real*8::error_rkf
      
  !Iteration variables
  integer :: i,j,l
  
  !Number of iterations
  integer :: n_iter     
!Number of steps along x axis
  integer :: nx
  !Number of steps along y axis
  integer :: ny
  integer::nn	!nx*ny
  
  
  !Fixing dx and dy 
  real :: hy, hy_inv, hy_inv2
  real :: hx, hx_inv, hx_inv2      
      !Fixing time stepping
  real :: dt 
  
  !Max reaction rate location
  integer max_locw
  
  real*8 :: w_max,w_sum
  real*8  :: e_max
  integer :: w_loc,e_loc
  
      
  !The LHS matrices for calculating the derivatives are the same every time. Declared with global scope. The LHS matrices are stored as two separate column matrices, two each for 
  
  !First derivative along x axis
  !Second derivative along x axis
  !Second derivative along y axis
  
  !a --> Diagonal matrix
!b --> Super Diagonal matrix
  !Sub Diagonal Matrix is just b's reverse      
  
  real*8, dimension(:), allocatable :: aX, bX, a2X, b2X, a2Y, b2Y
  
  ! time vars
  real t1,t2,t3,t4,t5,t6,t7,t8,t9,t10 
  integer :: nt
  
      
  !!damper parameters
  real*8:: chi,delta
  type(tfo1)::Dm11,Dm12,Dm21,Dm22
  type(tfo1)::ch ,eterm
  real*8::dqbar !!Max Norm of abs(qbar(t)-qbar(t-dt))
  real*8,parameter ::epsilon=1.0E-08
  logical:: file_exists
  character(20)::fname
  !Getting in domain data
  write(*,*)'Number of points in x direction (flow) ?'
  read(*,'(I10)')nx
  write(*,*)'Number of points in y direction (cross flow) ?'
  read(*,'(I10)')ny
  write(*,*)'Number of iterations '
  read(*,'(I10)')n_iter
  write(*,*)'Time step '
  read(*,*) dt
      
  ! set no threads
 ! nt=4
 !call omp_set_num_threads(nt)
  !get time
  write(*,*)"using ",  omp_get_max_threads()," no of threads"
  t1=omp_get_wtime()
!  call ft(t1)
  !Now some clerical work due to the stupidity of the FORTRAN language requiring all variable declarations first
  hx = 15.0/REAL(nx-1)  !Default length in x direction is 15.0
  hy = 30.0/REAL(ny-1)  !Default length in y direction is 30.0
  hx_inv = 1.0/hx
  hx_inv2 = hx_inv*hx_inv
  hy_inv = 1.0/hy
  hy_inv2 = hy_inv*hy_inv
      
  !Getting problem data
  
  !      phi = nu*YfO/YoO
  write(*,*)'Phi (Equivalence ratio) '
  read(*,*)phi
  
  !      Ta = To + Q*YfO/(Cp*(1.0+phi))
!      beta = Ea*(Ta - To)/(R_gas*Ta*Ta)
  write(*,*)'Value of Beta (Zeldovich number) '
  read(*,*)beta
  
  !      gamma = (Ta - To)/To
  write(*,*)'Value of Gamma '
  read(*,*)gamma
  
  !      Da = Dth*rho*B*YoO*exp(-Ea/(R_gas*To))/(beta*beta*beta*Uo*Uo) 
  write(*,*)'Value of Damkohler number'
  read(*,*)Da
  
  write(*,*)'Value of Fuel Lewis number '
  read(*,*)LeF
  
  write(*,*)'Value of Oxidiser Lewis number '
  read(*,*)LeO
      
  write(*,*)'Value of the parameter which decides the inlet profile'
  read(*,*)alphaparam
 
  !damper paramters
  write(*,*) 'enter value of chi damper paramter'
  read(*,*) chi

  write(*,*) 'enter vlaue of delta filter width'
  read(*,*) delta

  write(*,*) 'the value of chi is ',chi
  write(*,*) 'the value of delta is ',delta

  
  LeF_inv = 1.0/LeF
  LeO_inv = 1.0/LeO
  
!!!!Actual calculations!!!!
  
  !Array allocation for the LHS matrix vectors
  allocate(aX(nx))
  allocate(bX(nx-1))      
  allocate(a2X(nx))
  allocate(b2X(nx-1))
  allocate(a2Y(ny))
  allocate(b2Y(ny-1))      

  !Gaussian Elimination for the LHS matrix done here itself. Only once !! Great savings !!
  
  !First derivative along x axis
  aX  = 1.0
  bX(3:nx-2) = 1.0/3.0
  bX(1) = 3.0
  bX(2) = 0.25
  bX(nx-1) = 0.25
  do i = 2,nx
     aX(i) = aX(i) - bX(i-1)*bX(nx+1-i)/aX(i-1)         
  end do
      
  !Second derivative along x axis
  a2X  = 1.0
  b2X(3:nx-2) = 2.0/11.0
  b2X(1) = 10.0
  b2X(2) = 0.1
  b2X(nx-1) = 0.1
  do i = 2,nx
     a2X(i) = a2X(i) - b2X(i-1)*b2X(nx+1-i)/a2X(i-1)         
  end do
  
      !Second derivative along y axis
  a2Y  = 1.0
  b2Y(3:ny-2) = 2.0/11.0
  b2Y(1) = 10.0
  b2Y(2) = 0.1
  b2Y(ny-1) = 0.1
  do i = 2,ny
     a2Y(i) = a2Y(i) - b2Y(i-1)*b2Y(ny+1-i)/a2Y(i-1)         
  end do
  
  !calculate elements of damping matrix Dm
  ch%t=chi
  ch%f=chi
  ch%o=chi
  
  eterm%t=exp((-ch%t-1.0/delta)*dt)
  eterm%f=exp((-ch%f-1.0/delta)*dt)
  eterm%o=exp((-ch%o-1.0/delta)*dt)
 
  Dm11%t=(1+ch%t*delta*eterm%t)/(1+ch%t*delta)
  Dm12%t=(ch%t*delta*(1-eterm%t))/(1+ch%t*delta)
  Dm21%t=(1-eterm%t)/(1+ch%t*delta)
  Dm22%t=(ch%t*delta+eterm%t)/(1+ch%t*delta)
  
  Dm11%f=(1+ch%f*delta*eterm%f)/(1+ch%f*delta)
  Dm12%f=(ch%f*delta*(1-eterm%f))/(1+ch%f*delta)
  Dm21%f=(1-eterm%f)/(1+ch%f*delta)
  Dm22%f=(ch%f*delta+eterm%f)/(1+ch%f*delta)
  
    
  Dm11%o=(1+ch%o*delta*eterm%o)/(1+ch%o*delta)
  Dm12%o=(ch%o*delta*(1-eterm%o))/(1+ch%o*delta)
  Dm21%o=(1-eterm%o)/(1+ch%o*delta)
  Dm22%o=(ch%o*delta+eterm%o)/(1+ch%o*delta)

  write(*,*)'Damper matrix is (',Dm11,',',Dm12,',',Dm21,',',Dm22,')'
  !Array allocations the values in the domain
  allocate(arr_tfo(nx*ny))
  if(allocated(arr_tfo)) write(*,*)'First array allocated'
  allocate(arr_tfo1(nx*ny))
  if(allocated(arr_tfo1)) write(*,*)'filter array allocated'
  allocate(w(nx*ny)) 
  if(allocated(w)) write(*,*)'Reaction rate array allocated'
!  allocate(err(nx*ny))
!  if(allocated(err)) write(*,*)'Error norm array allocated'
      
      w(1:nx*ny) = 0.0
      
      arr_tfo(1:nx*ny)%f%val = 0.0
      arr_tfo(1:nx*ny)%o%val = 0.0
      arr_tfo(1:nx*ny)%t%val = 0.0
      

      !Restart code
      call restart(arr_tfo,arr_tfo1) 
    
!      call boundary(arr_tfo)  !!!not necessary for restart code
	nn=nx*ny
      allocate(k1(1:nn),arr_res(1:nn),arr_dt1(1:nn)) !, working arrays for time stepping
!!initiliaze tfo1
       fname="output_filter.dat"
       INQUIRE(FILE=fname, EXIST=file_exists)
       if(.not. file_exists) then
          write(6,*)" file ", fname, " does not exists. hence initializing filtered variable"
          !!initialize the filtered variable arr_tfo1
          arr_tfo1(1:nn)%t%val= arr_tfo(1:nn)%t%val
          arr_tfo1(1:nn)%f%val= arr_tfo(1:nn)%f%val
          arr_tfo1(1:nn)%o%val= arr_tfo(1:nn)%o%val
      end if
      
      open(10,file="err.dat", access="append", action="write")
      open(20,file="rrate.dat", access="append", action="write")
      open(50,file="norm.dat", access="append", action="write")


!      call ft(t2)
        t2=omp_get_wtime()

        !!write initial norms
        write(50,'(6F15.8)') sum(arr_tfo(1:nn)%t%val),sum(arr_tfo(1:nn)%f%val),sum(arr_tfo(1:nn)%o%val),&
             sum(arr_tfo1(1:nn)%t%val),sum(arr_tfo1(1:nn)%f%val),sum(arr_tfo1(1:nn)%o%val)

!!!!!!!!!!!!!iterations start      
      write(*,*)"iterations start will  be running ",n_iter,"iterations" 
     ! t8=0.0 
         
    !  dqbar=10.0
      do l = 1,n_iter
         
         !!call damped version of time stepping for q
         call time_step_damp(arr_tfo, dt,k1,arr_res,arr_dt1,arr_tfo1)
         
        !!call time stepping for filtered variable  qbar here arr_tfo1
         call time_step_filter(arr_tfo1, dt,k1,arr_res,arr_dt1,arr_tfo,dqbar)    
!$omp barrier
         !$omp master
        ! call ft(t10)
    !     call correct(arr_tfo)
    !     call correct(arr_tfo1)
         
     !calculates reaction rate
         call rrate_max(arr_tfo1,w_max,w_loc,w_sum)

         !!calculate new residue based on filtered varaibles to verify if f(q)=0
         call residual(arr_tfo1,e_max,e_loc)
         !writess err rrate to file
         write(10,'(ES25.16,1X,F15.8,1X,F15.8,1X,ES25.16)') e_max,mod(e_loc-1,nx)*hx, int(e_loc/nx)*hy - (ny-1)*hy/2.0,dqbar !, error_rkf
         
         write(20,'(F20.10,1X,F10.5,1X,F20.10)') w_max, mod(w_loc-1,nx)*hx, w_sum
         write(50,'(6F15.8)') sum(arr_tfo(1:nn)%t%val),sum(arr_tfo(1:nn)%f%val),sum(arr_tfo(1:nn)%o%val),&
              sum(arr_tfo1(1:nn)%t%val),sum(arr_tfo1(1:nn)%f%val),sum(arr_tfo1(1:nn)%o%val)

      !   write(*,*) e_max, w_max, w_sum,dqbar

         !time
       !  call ft(t9)
       !  t9=omp_get_wtime()

      !   t8=t8+t9-t10
         !$omp end master
         !check for stopping criteria
         if(l>100000 .and. dqbar <epsilon ) then
            write(*,*) "convergence criteria has met, hence stopping the simulation"
            write(*,*) "max precision error is ", dqbar
            write(*,*) "max accuracy error is ", e_max
            exit
         end if
      end do
     if(dqbar >epsilon ) write(*,*) "Max iterations reached but convergence goal not reached hence stopping the simulations"
      write(*,*) "max precision error is ", dqbar
      write(*,*) "max accuracy error is ", e_max

      ! ending time calculations
      !t3=omp_get_wtime()
      !call ft(t3)
      t3=omp_get_wtime()

      write(*,*) "total time for ",n_iter," iterations",t3-t2
      !write(*,*) 'total time for err',t8
      
      
      close(10)
      close(20)
      close(50)
   
      !Dump to file
      call  rrate_calc(arr_tfo,w)
      call dump(arr_tfo,w)
      call  rrate_calc(arr_tfo1,w)
      call dump_filtered(arr_tfo1,w)
      
      deallocate(arr_tfo,w,arr_tfo1,arr_dt1,arr_res)
      deallocate(aX, bX, a2X, b2X, a2Y, b2Y)
      t4=omp_get_wtime()
     ! call ft(t4)
      write(*,*) 'total time for prg=',t4-t2
      
      !Subroutine definitions
    contains
      
      subroutine read_input(nl_iter,dt)
        integer :: nl_iter
        integer :: nn1
        real::dt
        open(1,file='le1p7inputs.txt',status="old",action="read",access="sequential")
        
        read(1,*)nn1
        read(1,*)nn1
        read(1,*)nl_iter
        read(1,*)dt
      end subroutine read_input
      
     subroutine restart(arr_in,arr_filt)
       !Code restarting
       type(TFO), dimension(:), intent(out) :: arr_in,arr_filt
       logical :: file_exists
       character(20)::fname
       real :: x1,y1 ! Dummy real variables.. to read those x and y values from the restart file
       real :: r1 !Dummy real variable to read reaction rate
       real dai   !dummy check damkoler no
       character :: a !Dummy variable to read blank line
       integer ::p,q
       open(25,file="output_inc.dat", status="old", action="read", access="sequential")
       
       read(25,*) x1,y1,dai
       write(*,*)"output file supplied for Le=", x1,"and alphaparam=",y1,"da=",dai
       do p = 0,ny-1
          do q = 1,nx
             read(25,'(F15.8,1X,F15.8,1X,F25.16,1X,F25.16,1X,F25.16,1X,F25.16)') &
                  x1, y1, arr_in(p*nx+q)%t%val, arr_in(p*nx+q)%f%val, arr_in(p*nx+q)%o%val, r1
          end do
          read(25,'(a)')a
       end do
       close(25)
       !Fixing the boundary conditions
       call read_inlet(arr_in)
       !!check if the filtered data is avialble in a file
       
       fname="output_filter.dat"
       INQUIRE(FILE=fname, EXIST=file_exists)
       if(file_exists) then
          write(6,*) "file ",fname," exists  reading filtered variable data..."
          
          !!read filtered array
          open(205,file=fname, status="old", action="read", access="sequential")
          read(205,*) x1,y1,dai
          write(*,*)"output file supplied for Le=", x1,"and alphaparam=",y1,"da=",dai
          do p = 0,ny-1
             do q = 1,nx
                read(205,'(F15.8,1X,F15.8,1X,F25.16,1X,F25.16,1X,F25.16,1X,F25.16)') &
                     x1, y1, arr_filt(p*nx+q)%t%val, arr_filt(p*nx+q)%f%val, arr_filt(p*nx+q)%o%val, r1
             end do
             read(205,'(a)')a
          end do
          close(205)
       end if
       
       
       write(*,*)'code restarted.. hope for speedy convergence !!'
       
     end subroutine restart
     
    subroutine dump(arr_out,rrate)
          !Writes the given array object to file "output_inc.dat"
          type(TFO), dimension(:), intent(in) :: arr_out
          real*8, dimension(:), intent(in) :: rrate
          integer :: p,q !Iteration variables     

          open(17,file="output_inc.dat", status="replace", action="write", access="sequential")
          write(17,'(F15.8,1X,F15.8,1X,F15.8)')LeF,alphaparam,da
          
          do p = 0,ny-1
             do q = 1,nx
                write(17,'(F15.8,1X,F15.8,1X,F25.16,1X,F25.16,1X,F25.16,1X,F25.16)')(q-1)*hx, p*hy - (ny-1)*hy/2, &
                     arr_out(p*nx+q)%t%val, arr_out(p*nx+q)%f%val, arr_out(p*nx+q)%o%val, rrate(p*nx+q)
            ! write(*,*)(q-1)*hx, p*hy - (ny-1)*hy/2, &
             !        arr_out(p*nx+q)%t%val, arr_out(p*nx+q)%f%val, arr_out(p*nx+q)%o%val, rrate(p*nx+q)
             end do
             write(17,*)''
          end do
          close(17)
        end subroutine dump
    subroutine dump_filtered(arr_out,rrate)
          !Writes the given array object to file "output_inc.dat"                                                                                                                                                  
          type(TFO), dimension(:), intent(in) :: arr_out
          real*8, dimension(:), intent(in) :: rrate
          integer :: p,q !Iteration variables                                                                                                                                                                      

          open(107,file="output_filter.dat", status="replace", action="write", access="sequential")
          write(107,'(F15.8,1X,F15.8,1X,F15.8)')LeF,alphaparam,da

          do p = 0,ny-1
             do q = 1,nx
                write(107,'(F15.8,1X,F15.8,1X,F25.16,1X,F25.16,1X,F25.16,1X,F25.16)')(q-1)*hx, p*hy - (ny-1)*hy/2, &
                     arr_out(p*nx+q)%t%val, arr_out(p*nx+q)%f%val, arr_out(p*nx+q)%o%val, rrate(p*nx+q)
            ! write(*,*)(q-1)*hx, p*hy - (ny-1)*hy/2, &                                                                                                                                                            
             !        arr_out(p*nx+q)%t%val, arr_out(p*nx+q)%f%val, arr_out(p*nx+q)%o%val, rrate(p*nx+q)                                                                                                          
             end do
	     write(107,*)''
          end do
          close(107)
        end subroutine dump_filtered

subroutine time_step(arr_dt, delta_t,k1,arr_res,arr_dt1)
  
  !Time stepping using RKF method
  !          use func
  use tfo_operation
  type(TFO), dimension(:), intent(inout) :: arr_dt !The current time array

          real, intent(in) :: delta_t
          type(TFO), dimension(:),intent(inout) :: k1,arr_res,arr_dt1 !,arr_res !,arr_err !! !The intermediate working arrays functrion ealuation for rkf methods
!         real*8::err_res          !!local variables

          !Iteration variables
integer :: i,j
!      write(*,*)nx,ny,nz
!assigning arr_dt0
!!!!function evaluaton
          call deriv_boundary(arr_dt)
          call fd(Arr_dt,arr_res) 

!$omp parallel default(shared)
!$omp do private(i) schedule(static)           
          do i=1,nn
             k1(i)=delta_t*arr_res(i)
             arr_res(i)=k1(i)
             arr_dt1(i)=arr_dt(i)+0.5*k1(i)
          end do
!$omp end do
!$omp end parallel

          call deriv_boundary(arr_dt1)
          call fd(Arr_dt1,k1)

!$omp parallel default(shared)
!$omp do private(i) schedule(static)
          do i=1,nn
             k1(i)=delta_t*k1(i)
             arr_res(i)=arr_res(i)+2.0*k1(i)
             arr_dt1(i)= arr_dt(i)+0.5*k1(i)
          end do
!$omp end do
!$omp end parallel                  

          call deriv_boundary(arr_dt1)
          call fd(Arr_dt1,k1)

!$omp parallel default(shared)
!$omp do private(i) schedule(static)
          do i=1,nn
             k1(i)=delta_t*k1(i)    !k3
             arr_res(i)=arr_res(i)+2.0*k1(i) 
             arr_dt1(i)= arr_dt(i)+k1(i)
          end do
!$omp end do
!$omp end parallel                                                                                                           

          call deriv_boundary(arr_dt1)
          call fd(Arr_dt1,k1)
          !last iteration 


!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(static)
          do i=1,nn
             k1(i)=delta_t*k1(i)
             arr_res(i)=arr_res(i)+k1(i)          
             arr_res(i)=(1.0/6.0)*arr_res(i)  !
             arr_dt(i) =arr_dt(i)+arr_res(i)
          end do
!$omp end parallel do

!$omp barrier

         
         
        end subroutine time_step
        
        subroutine time_step_damp(arr_dt, delta_t,k1,arr_res,arr_dt1,arr_dtbar)
  
  !Time stepping using RKF method
  !          use func
  use tfo_operation
  type(TFO), dimension(:), intent(inout) :: arr_dt !The current time array
  type(TFO), dimension(:), intent(in):: arr_dtbar  !steady state array

          real, intent(in) :: delta_t
          type(TFO), dimension(:),intent(inout) :: k1,arr_res,arr_dt1 !,arr_res !,arr_err !! !The intermediate working arrays functrion ealuation for rkf methods
!         real*8::err_res          !!local variables

          !Iteration variables
integer :: i,p,q
!      write(*,*)nx,ny,nz
!assigning arr_dt0
!!!!function evaluaton
          call deriv_boundary(arr_dt)
          call fd(Arr_dt,arr_res) 
          call fdamp(Arr_dt,arr_dtbar,arr_res) 
!$omp parallel default(shared)
!$omp do private(i,p,q) schedule(static)           
         do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*arr_res(i)
             arr_res(i)=k1(i)
             arr_dt1(i)=arr_dt(i)+0.5*k1(i)
          end do
         end do
!$omp end do
!$omp end parallel

          call deriv_boundary(arr_dt1)
          call fd(Arr_dt1,k1)
          call fdamp(Arr_dt1,arr_dtbar,k1) 
          call correct(arr_dt1)
!$omp parallel default(shared)
!$omp do private(i,p,q) schedule(static)
          do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*k1(i)
             arr_res(i)=arr_res(i)+2.0*k1(i)
             arr_dt1(i)= arr_dt(i)+0.5*k1(i)
          end do
         end do
          !$omp end do
          !$omp end parallel
          call correct(arr_dt1)


          call deriv_boundary(arr_dt1)
          call fd(Arr_dt1,k1)
                    call fdamp(Arr_dt1,arr_dtbar,k1) 
         !$omp parallel default(shared)
!$omp do private(i,p,q) schedule(static)
         do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*k1(i)    !k3
             arr_res(i)=arr_res(i)+2.0*k1(i) 
             arr_dt1(i)= arr_dt(i)+k1(i)
          end do
         end do
!$omp end do
!$omp end parallel                                                                                                           
          call correct(arr_dt1)

          call deriv_boundary(arr_dt1)
          call fd(Arr_dt1,k1)
          call fdamp(Arr_dt1,arr_dtbar,k1) 
          !last iteration 


!$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i,p,q) SCHEDULE(static)
     do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*k1(i)
             arr_res(i)=arr_res(i)+k1(i)          
             arr_res(i)=(1.0/6.0)*arr_res(i)  !
             arr_dt(i) =arr_dt(i)+arr_res(i)

          end do
         end do
!$omp end parallel do
          call correct(arr_dt)

!$omp barrier

         
         
        end subroutine time_step_damp
        
      subroutine time_step_filter(arr_dt, delta_t,k1,arr_res,arr_dt1,arr_unfiltered,err_prec)
  
  !Time stepping using RKF method
  !          use func
  use tfo_operation
  type(TFO), dimension(:), intent(inout) :: arr_dt !The current time array to be marched in this case qbar
  type(TFO), dimension(:), intent(in)::arr_unfiltered  !!unfiltered array q
          real, intent(in) :: delta_t
          real*8, intent(out)::err_prec
          type(TFO), dimension(:),intent(inout) :: k1,arr_res,arr_dt1 !,arr_res !,arr_err !! !The intermediate working arrays functrion ealuation for rkf methods
         real*8::err_res          !!local variables

          !Iteration variables
integer :: i,p,q
!      write(*,*)nx,ny,nz
!assigning arr_dt0
!!!!function evaluaton
          !call deriv_boundary(arr_dt)
    !      call fd(Arr_dt,arr_res) 
          call fsdamp(arr_dt,arr_unfiltered,arr_res)
!$omp parallel default(shared)
!$omp do private(i,p,q) schedule(static)           
      do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*arr_res(i)
             arr_res(i)=k1(i)
             arr_dt1(i)=arr_dt(i)+0.5*k1(i)
          end do
         end do
!$omp end do
          call correct(arr_dt1)
!$omp end parallel

         ! call deriv_boundary(arr_dt1)
        !  call fd(Arr_dt1,k1)
          call fsdamp(arr_dt1,arr_unfiltered,k1)
!$omp parallel default(shared)
!$omp do private(i,p,q) schedule(static)
        do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*k1(i)
             arr_res(i)=arr_res(i)+2.0*k1(i)
             arr_dt1(i)= arr_dt(i)+0.5*k1(i)
          end do
         end do
!$omp end do
!$omp end parallel 
          call correct(arr_dt1)                 

        !  call deriv_boundary(arr_dt1)
        !  call fd(Arr_dt1,k1)
                    call fsdamp(arr_dt1,arr_unfiltered,k1)

!$omp parallel default(shared)
!$omp do private(i,p,q) schedule(static)
          do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*k1(i)    !k3
             arr_res(i)=arr_res(i)+2.0*k1(i) 
             arr_dt1(i)= arr_dt(i)+k1(i)
          end do
         end do
!$omp end do
!$omp end parallel 
          call correct(arr_dt1)                                                                                                          

       !   call deriv_boundary(arr_dt1)
        !  call fd(Arr_dt1,k1)
          call fsdamp(arr_dt1,arr_unfiltered,k1)
          !last iteration 

err_prec=0.0  !! minimum value of err_res is zero due to abs

          !$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i,err_res,p,q) SCHEDULE(static) REDUCTION(MAX:err_prec)
        do p = 0,ny-1
             do q = 1,nx
             i=p*nx+q
             k1(i)=delta_t*k1(i)
             arr_res(i)=arr_res(i)+k1(i)          
             arr_res(i)=(1.0/6.0)*arr_res(i)  
             arr_dt(i) =arr_dt(i)+arr_res(i)
             err_res=abs(arr_res(i)%t%val)+abs(arr_res(i)%f%val)+abs(arr_res(i)%o%val)
             if(err_res>err_prec) then 
             err_prec=err_res
              end if
          end do
                  end do
!$omp end parallel do
          call correct(arr_dt)

!$omp barrier

        end subroutine time_step_filter  
        
        subroutine correct(arr_dt)
          use data_struct
          type(tfo),dimension(:),intent(inout)::arr_dt
          !$OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(static) 
          do i=1,nn
             if(arr_dt(i)%t%val < 0.0) then
                arr_dt(i)%t%val=0.0
             end if
             if(arr_dt(i)%f%val <0.0) then
                arr_dt(i)%f%val=0.0
             end if
             if(arr_dt(i)%o%val <0.0) then
                arr_dt(i)%o%val=0.0
             end if
                          if(arr_dt(i)%t%val >1) then
                arr_dt(i)%t%val=1.0
             end if
             if(arr_dt(i)%f%val >1) then
                arr_dt(i)%f%val=1.0
             end if
             if(arr_dt(i)%o%val >1) then
                arr_dt(i)%o%val=1.0
             end if
          end do
          !$omp end parallel do                                                                                                                                                                                   
        end   subroutine correct
!!function to dampen the oscillations
        subroutine damper(arr_tfo,arr_tfo1)
  use data_struct
  type(TFO), dimension(:), intent(inout) :: arr_tfo,arr_tfo1 !The current time array and the filtered array                                                                                                                                            
  type(TFO) :: q1,q1bar !the intermediate working arrays 
  !store first the newqly calucate varibales q and qbar in q1 and q1bar and then assign to array
!    write(*,*)'Damper matrix is (',Dm11,',',Dm12,',',Dm21,',',Dm22,')'

            do i=1,nn
            !calculate new values by multiplying the damper matrix to (q,qbar)
            q1=Dm11*arr_tfo(i)+Dm12*arr_tfo1(i)
            q1bar=Dm12*arr_tfo(i)+Dm22*arr_tfo1(i)
!!update the variables
            arr_tfo(i)=q1
            arr_tfo1(i)=q1bar            
            end	do

end subroutine damper
!!serial subroutine to find max residue absolute value

        subroutine residual(arr_dt,err_max,err_loc)
                  type(TFO), dimension(:), intent(inout) :: arr_dt !The current time array
                  real*8,intent(out)::err_max  !!max residue
                            integer,intent(out)::err_loc  !!location of max residue
             real*8::err_res  !!local variable
             !!calculate derivatives
             call deriv_boundary(arr_dt)
             !!calculate residual
                  call fd(arr_dt,arr_res)
                  err_max=0.0
          	err_loc=0

! $OMP  PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(static) REDUCTION(MAX:err_res)
          do i=1,nn  
             err_res=abs(arr_res(i)%t%val)+abs(arr_res(i)%f%val)+abs(arr_res(i)%o%val)
             if (err_res.gt.err_max) then
                err_max=err_res
                err_loc=i
             end if
          end do
! $omp end parallel do

 ! $omp barrier
        end subroutine residual
        
!parallel subroutine for calculatinfg max err and max reaction rate

        subroutine max_err1(arr_dt,arr_dt1,w,w_max,err_max,w_loc,e_loc)
          implicit none
          type(TFO), dimension(:), intent(in) :: arr_dt !The current time array
          type(TFO), dimension(:), intent(in) :: arr_dt1 !The updated array
          real*8, dimension(:), intent(in) :: w !reaction rate array
          real*8,intent(out) :: w_max
          real*8,intent(out) :: err_max
          integer, intent(out) :: w_loc,e_loc        
          !!local variables
          !Iteration variables
          integer :: i,j,nn1
          real*8, dimension(:), allocatable :: w1
          real*8, dimension(:), allocatable :: e1, e2,errd
          integer, dimension(:), allocatable :: el, wl
          real*8 :: err
          allocate(w1(nt))
          allocate(e1(nt))
          allocate(wl(nt))
          allocate(el(nt))
          allocate(errd(nx*ny))
          !write(*,*) 'array allocated succ'
          w1(1:nt)=-1000.0
          e1(1:nt)=-1000.0  
          el(1:nt)=0
          wl(1:nt)=0
          nn1=int(nx*ny/nt)
          !$omp parallel default(shared)
          !$omp do private(j,i) schedule(static,1)
          do j = 1,nt
             do i = nn1*(j-1)+1,nn1*j
                If (w1(j) .lt. w(i)) then
                   wl(j)=i 
                   w1(j)=w(i) 
                end if
     
                err=abs(arr_dt1(i)%t%val - arr_dt(i)%t%val)+           &
                     abs(arr_dt1(i)%o%val-arr_dt(i)%o%val) + abs(arr_dt1(i)%f%val - arr_dt(i)%f%val)
                
                If (e1(j) .lt. err) then 
   el(j)=i 
    e1(j)=err 
!                        write(*,*) i,j,err,e1(j),el(j)
                        !                        write(*,*) i,j,w(i),w1(j),wl(j)
                     end if
                  end do
               end do
               !$omp end do 
               !$omp end parallel
               !updating the wmax        
               w_max=-100
               err_max=-100
               e_loc=0
               w_loc=0
               
               !calculating max among the threads
               do i=1,nt
                  If (w_max .lt. w1(i)) then
                     w_max=w1(i)
                     w_loc=wl(i)
                  end if
                  
                  If (err_max .lt. e1(i)) then
                     err_max=e1(i)
                     e_loc=el(i)
                  end if
               end do
               
               ! calculating max for remaining terms       
               if (mod(nx*ny,nt) .ne. 0) then
                  do i=nn*nt+1,nx*ny
                     If (w_max .lt. w(i)) then
                        w_max=w(i)
                        w_loc=i
                     end if
                     
                     err=abs(arr_dt1(i)%t%val - arr_dt(i)%t%val)+ &
             abs(arr_dt1(i)%o%val-arr_dt(i)%o%val) + abs(arr_dt1(i)%f%val - arr_dt(i)%f%val)
                     
                     If (err_max .lt. err) then
                        err_max=err
                        e_loc=i
                     end if
                  end do
               end if
               
  errd=abs(arr_dt1%t%val - arr_dt%t%val)+ abs(arr_dt1%o%val-arr_dt%o%val) + abs(arr_dt1%f%val - arr_dt%f%val)
  write(*,*) 'actual max error and loc =',maxval(errd),mod(maxloc(errd)-1,nx)*hx, int(maxloc(errd)/nx)*hy - (ny-1)*hy/2.0
  write(*,'(F25.16,1X,F15.8,1X,F15.8)')err_max, int(e_loc/nx)*hy - (ny-1)*hy/2.0, mod(e_loc-1,nx)*hx

  write(*,*) 'actual max w and loc', maxval(w),mod(maxloc(w)-1,nx)*hx, int(maxloc(w)/nx)*hy - (ny-1)*hy/2.0
  write(*,'(F25.16,1X,F15.8,1X,F15.8)')w_max, mod(w_loc-1,nx)*hx

 deallocate(e1,w1,el,wl,errd)

             end subroutine max_err1

!subroutine for reading the inlet concentration profile   
             subroutine read_inlet(arr_in)
               !Code restarting
               type(TFO), dimension(:), intent(inout) :: arr_in
               integer ::p               
               open(205,file="inlet.dat", status="old", action="read", access="sequential")
               !         read(25,*) x1,y1,dai
               !        write(*,*)"output file supplied for Le=", x1,"and alphaparam=",y1,"da=",dai
         
               do p = 1, (ny-1)*nx + 1, nx
                  read(205,*)arr_in(p)%f%val  !,  arr_in(p)%o%val
                  arr_in(p)%o%val = 1-arr_in(p)%f%val
                  arr_in(p)%t%val=0.0 
                 ! write(*,*)arr_in(p)%f%val,arr_in(p)%o%val
               end do
               
               close(205)
          !convert input values to reqd precision
  !        arr_in%t%val=real(arr_in%t%val,8)
   !       arr_in%f%val=real(arr_in%f%val,8)
 !         arr_in%f%val=real(arr_in%f%val,8)
 
 !write(*,*)'code restarted.. hope for speedy convergence !!'
          
             end subroutine read_inlet

   
             subroutine boundary(bc_ar)
               type(TFO), dimension(:),intent(inout) :: bc_ar
               ! real, intent(in) :: xparam
               
               integer i !Iteration variable
               ! real :: yloc,yloc0  !Used for defining the y location in fixing the in flux BC
               
!!!Derivative zero bc's SIXTH ORDER ACCURATE FOR NOW
               !!Left end bc's  in y direction  at y= -h
               !$omp parallel default(shared)
               !$omp do private(i) schedule(static,25)  
               
               do i = 2, nx-1 
                  bc_ar(i)%t%val = (360.0*bc_ar(i+nx)%t%val - 450.0*bc_ar(i+2*nx)%t%val + 400.0*bc_ar(i+3*nx)%t%val&
                       - 225.0*bc_ar(i+4*nx)%t%val + 72.0*bc_ar(i+5*nx)%t%val - 10.0*bc_ar(i+6*nx)%t%val  )/147.0
                  !            bc_ar(i)%t%val = bc_ar(i+nx)%t%val
                  bc_ar(i)%o%val = (360.0*bc_ar(i+nx)%o%val - 450.0*bc_ar(i+2*nx)%o%val + 400.0*bc_ar(i+3*nx)%o%val&
                       - 225.0*bc_ar(i+4*nx)%o%val + 72.0*bc_ar(i+5*nx)%o%val - 10.0*bc_ar(i+6*nx)%o%val  )/147.0
                  !            bc_ar(i)%o%val = bc_ar(i+nx)%o%val
                  bc_ar(i)%f%val = (360.0*bc_ar(i+nx)%f%val - 450.0*bc_ar(i+2*nx)%f%val + 400.0*bc_ar(i+3*nx)%f%val&
                       - 225.0*bc_ar(i+4*nx)%f%val + 72.0*bc_ar(i+5*nx)%f%val - 10.0*bc_ar(i+6*nx)%f%val  )/147.0
                  !            bc_ar(i)%f%val = bc_ar(i+nx)%f%val

               end do
               !$omp end do
               !!Right end bc's   y bc 2  at y =h
               !$omp do private(i) schedule(static,25)
               do i = (ny-1)*nx + 2, ny*nx - 1 
                  bc_ar(i)%t%val = (360.0*bc_ar(i-nx)%t%val - 450.0*bc_ar(i-2*nx)%t%val + 400.0*bc_ar(i-3*nx)%t%val&
                       - 225.0*bc_ar(i-4*nx)%t%val + 72.0*bc_ar(i-5*nx)%t%val - 10.0*bc_ar(i-6*nx)%t%val  )/147.0
                  !            bc_ar(i)%t%schedule(static,25)val = bc_ar(i-nx)%t%val
                  bc_ar(i)%o%val = (360.0*bc_ar(i-nx)%o%val - 450.0*bc_ar(i-2*nx)%o%val + 400.0*bc_ar(i-3*nx)%o%val&
                       - 225.0*bc_ar(i-4*nx)%o%val + 72.0*bc_ar(i-5*nx)%o%val - 10.0*bc_ar(i-6*nx)%o%val  )/147.0
                  !            bc_ar(i)%o%val = bc_ar(i-nx)%o%val
                  bc_ar(i)%f%val = (360.0*bc_ar(i-nx)%f%val - 450.0*bc_ar(i-2*nx)%f%val + 400.0*bc_ar(i-3*nx)%f%val&
                       - 225.0*bc_ar(i-4*nx)%f%val + 72.0*bc_ar(i-5*nx)%f%val - 10.0*bc_ar(i-6*nx)%f%val  )/147.0
!            bc_ar(i)%f%val = bc_ar(i-nx)%f%val
               end do
               !$omp end do
               !!Top end (Far field bc's)   at x=L
               !$omp do private(i) schedule(static,25)
               do i = 2*nx, nx*(ny-1), nx             
                  bc_ar(i)%t%val = (360.0*bc_ar(i-1)%t%val - 450.0*bc_ar(i-2)%t%val + 400.0*bc_ar(i-3)%t%val&
                       - 225.0*bc_ar(i-4)%t%val + 72.0*bc_ar(i-5)%t%val - 10.0*bc_ar(i-6)%t%val  )/147.0
                  !            bc_ar(i)%t%val = bc_ar(i-1)%t%val
                  bc_ar(i)%o%val = (360.0*bc_ar(i-1)%o%val - 450.0*bc_ar(i-2)%o%val + 400.0*bc_ar(i-3)%o%val&
                       - 225.0*bc_ar(i-4)%o%val + 72.0*bc_ar(i-5)%o%val - 10.0*bc_ar(i-6)%o%val  )/147.0
!            bc_ar(i)%o%val = bc_ar(i-1)%o%val
                  bc_ar(i)%f%val = (360.0*bc_ar(i-1)%f%val - 450.0*bc_ar(i-2)%f%val + 400.0*bc_ar(i-3)%f%val&
                       - 225.0*bc_ar(i-4)%f%val + 72.0*bc_ar(i-5)%f%val - 10.0*bc_ar(i-6)%f%val  )/147.0
                  !            bc_ar(i)%f%val = bc_ar(i-1)%f%val
                  
               end do
               !$omp end do 
               !$omp end parallel
               
               
               close(3)
             end subroutine boundary
             
             subroutine deriv_x(arr_derx, h_inv, n)
               !This function fixes the value of the first 'x' derivative of the given array. The derivative is calculated using a Tridiagonal scheme. It is 6 order accurate at the centre and 4th accurate at the boundaries

!    Inputs:
!    arr_derx --> The array for which the first 'x' derivative is to be found. The input array is of type TFO. Since the arrays are passed by reference, the derivatives are stored in the array itself in an appropriate manner.
!    n --> The length of the input array
!    h_inv --> (1/h) where h is the space discretization
          
          !!Declaration of input variables
          type(TFO), dimension(:)  :: arr_derx
          real, intent(in) :: h_inv
          integer, intent(in) :: n

          !!Variables declared within the function
          !Iteration variables
          integer :: k
          !The coefficients of the derivatives
          real alpha, a, b
          !The RHS matrices used to solve the tridiagonal equations. Made allocatable and later made of dimension 'n' --> the input array size
          real*8, dimension(:), allocatable :: BT, BO, BF

          allocate(BT(n))
          allocate(BO(n))
          allocate(BF(n))

!!! choice of alpha  
          alpha = 1.0/3.0
    
!!Deciding other parameters from alpha
!!beta = 0.. tridiagonal scheme
          a = (alpha + 9.0)/6.0
          b = (32.0*alpha - 9.0)/15.0
!!c = (-3*alpha +1)/10
          
!Calculating the RHS matrices
          !Inner nodes
          BT(3:n-2) = ( b*(arr_derx(5:)%t%val - arr_derx(:n-4)%t%val)*0.5 + &
               a*(arr_derx(4:n-1)%t%val - arr_derx(2:n-3)%t%val) )*0.5*h_inv
          BO(3:n-2) = ( b*(arr_derx(5:)%o%val - arr_derx(:n-4)%o%val)*0.5 + &
               a*(arr_derx(4:n-1)%o%val - arr_derx(2:n-3)%o%val) )*0.5*h_inv
          BF(3:n-2) = ( b*(arr_derx(5:)%f%val - arr_derx(:n-4)%f%val)*0.5 + &
               a*(arr_derx(4:n-1)%f%val - arr_derx(2:n-3)%f%val) )*0.5*h_inv

          !Penultimate nodes 4th order accurate
          BT(2) = 3*(arr_derx(3)%t%val - arr_derx(1)%t%val)*0.25*h_inv
          BT(n-1) = 3*(arr_derx(n)%t%val - arr_derx(n-2)%t%val)*0.25*h_inv
          BO(2) = 3*(arr_derx(3)%o%val - arr_derx(1)%o%val)*0.25*h_inv 
          BO(n-1) = 3*(arr_derx(n)%o%val - arr_derx(n-2)%o%val)*0.25*h_inv
          BF(2) = 3*(arr_derx(3)%f%val - arr_derx(1)%f%val)*0.25*h_inv
          BF(n-1) = 3*(arr_derx(n)%f%val - arr_derx(n-2)%f%val)*0.25*h_inv

          !Boundary nodes 4th order accurate
          BT(1) = (-17.0*arr_derx(1)%t%val + 9.0*(arr_derx(2)%t%val + arr_derx(3)%t%val)  - arr_derx(4)%t%val)*h_inv/6.0
          BO(1) = (-17.0*arr_derx(1)%o%val + 9.0*(arr_derx(2)%o%val + arr_derx(3)%o%val)  - arr_derx(4)%o%val)*h_inv/6.0
          BF(1) = (-17.0*arr_derx(1)%f%val + 9.0*(arr_derx(2)%f%val + arr_derx(3)%f%val)  - arr_derx(4)%f%val)*h_inv/6.0
          BT(n) = (17.0*arr_derx(n)%t%val - 9.0*(arr_derx(n-1)%t%val + arr_derx(n-2)%t%val)  + arr_derx(n-3)%t%val)*h_inv/6.0
          BO(n) = (17.0*arr_derx(n)%o%val - 9.0*(arr_derx(n-1)%o%val + arr_derx(n-2)%o%val)  + arr_derx(n-3)%o%val)*h_inv/6.0
          BF(n) = (17.0*arr_derx(n)%f%val - 9.0*(arr_derx(n-1)%f%val + arr_derx(n-2)%f%val)  + arr_derx(n-3)%f%val)*h_inv/6.0


          !Gaussian Elimination for the RHS matrices alone now.

          do k=2,n
             BT(k) = BT(k) - BT(k-1)*bX(n+1-k)/aX(k-1)
             BO(k) = BO(k) - BO(k-1)*bX(n+1-k)/aX(k-1)
             BF(k) = BF(k) - BF(k-1)*bX(n+1-k)/aX(k-1)             
          end do
           
          !Solving backwards.. Final step

          arr_derx(nx)%t%derx = BT(nx)/aX(nx)
          arr_derx(nx)%o%derx = BO(nx)/aX(nx)
          arr_derx(nx)%f%derx = BF(nx)/aX(nx)
          
          do k=n-1,1,-1
             arr_derx(k)%t%derx = (BT(k) - bX(k)*arr_derx(k+1)%t%derx)/aX(k)
             arr_derx(k)%o%derx = (BO(k) - bX(k)*arr_derx(k+1)%o%derx)/aX(k)
             arr_derx(k)%f%derx = (BF(k) - bX(k)*arr_derx(k+1)%f%derx)/aX(k)
          end do

          deallocate(BT, BO, BF)


        end subroutine deriv_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine deriv2_x(arr_der2x, h_inv2, n)
!This function fixes the value of the second 'x' derivative of the given array. The derivative is calculated using a Tridiagonal scheme. It is 6 order accurate at the centre and 4th accurate at the boundaries

!    Inputs:
!    arr_der2x --> The array for which the second 'x' derivative is to be found. The input array is of type TFO. Since the arrays are passed by reference, the derivatives are stored in the array itself in an appropriate manner.
!    n --> The length of the input array
!    h_inv2 --> (1/h)^2 where h is the space discretization
          
          !!Declaration of input variables
          type(TFO), dimension(:)  :: arr_der2x
          real, intent(in) :: h_inv2
          integer, intent(in) :: n

          !!Variables declared within the function
          !Iteration variables
          integer :: k
          !The coefficients of the derivatives
          real alpha, a, b
          !The RHS matrices used to solve the tridiagonal equations. Made allocatable and later made of dimension 'n' --> the input array size
          real*8, dimension(:), allocatable :: BT, BO, BF

          allocate(BT(n))
          allocate(BO(n))
          allocate(BF(n))

!!! choice of alpha
          alpha = 2.0/11.0
    
!!Deciding other parameters from alpha
!!beta = 0.. tridiagonal scheme
          a = 4.0*(1.0 - alpha)/3.0
          b = (10.0*alpha - 1.0)/3.0
!c = 0
          
!Calculating the RHS matrices
          !Inner nodes
          BT(3:n-2) = ( b*(arr_der2x(5:)%t%val - 2.0*arr_der2x(3:n-2)%t%val + arr_der2x(:n-4)%t%val)*0.25 + &
               a*(arr_der2x(4:n-1)%t%val  - 2.0*arr_der2x(3:n-2)%t%val + arr_der2x(2:n-3)%t%val) )*h_inv2
          BO(3:n-2) = ( b*(arr_der2x(5:)%o%val - 2.0*arr_der2x(3:n-2)%o%val + arr_der2x(:n-4)%o%val)*0.25 + &
               a*(arr_der2x(4:n-1)%o%val  - 2.0*arr_der2x(3:n-2)%o%val + arr_der2x(2:n-3)%o%val) )*h_inv2
          BF(3:n-2) = ( b*(arr_der2x(5:)%f%val - 2.0*arr_der2x(3:n-2)%f%val + arr_der2x(:n-4)%f%val)*0.25 + &
               a*(arr_der2x(4:n-1)%f%val  - 2.0*arr_der2x(3:n-2)%f%val + arr_der2x(2:n-3)%f%val) )*h_inv2

          !Penultimate nodes 4th order accurate
          BT(2) = 1.2*(arr_der2x(3)%t%val - 2.0*arr_der2x(2)%t%val + arr_der2x(1)%t%val)*h_inv2
          BT(n-1) =  1.2*(arr_der2x(n)%t%val - 2.0*arr_der2x(n-1)%t%val + arr_der2x(n-2)%t%val)*h_inv2
          BO(2) =  1.2*(arr_der2x(3)%o%val - 2.0*arr_der2x(2)%o%val + arr_der2x(1)%o%val)*h_inv2
          BO(n-1) =  1.2*(arr_der2x(n)%o%val - 2.0*arr_der2x(n-1)%o%val + arr_der2x(n-2)%o%val)*h_inv2
          BF(2) =  1.2*(arr_der2x(3)%f%val - 2.0*arr_der2x(2)%f%val + arr_der2x(1)%f%val)*h_inv2
          BF(n-1) =  1.2*(arr_der2x(n)%f%val - 2.0*arr_der2x(n-1)%f%val + arr_der2x(n-2)%f%val)*h_inv2

          !Boundary nodes 4th order accurate
          BT(1) = (145.0*arr_der2x(1)%t%val - 304.0*arr_der2x(2)%t%val + 174.0*arr_der2x(3)%t%val  - &
               16.0*arr_der2x(4)%t%val + arr_der2x(5)%t%val)*h_inv2/12.0
          BO(1) = (145.0*arr_der2x(1)%o%val - 304.0*arr_der2x(2)%o%val + 174.0*arr_der2x(3)%o%val  - &
               16.0*arr_der2x(4)%o%val + arr_der2x(5)%o%val)*h_inv2/12.0
          BF(1) = (145.0*arr_der2x(1)%f%val - 304.0*arr_der2x(2)%f%val + 174.0*arr_der2x(3)%f%val  - &
               16.0*arr_der2x(4)%f%val + arr_der2x(5)%f%val)*h_inv2/12.0

          BT(n) = (145.0*arr_der2x(n)%t%val - 304.0*arr_der2x(n-1)%t%val + 174.0*arr_der2x(n-2)%t%val  - &
               16.0*arr_der2x(n-3)%t%val + arr_der2x(n-4)%t%val )*h_inv2/12.0
          BO(n) = (145.0*arr_der2x(n)%o%val - 304.0*arr_der2x(n-1)%o%val + 174.0*arr_der2x(n-2)%o%val  - &
               16.0*arr_der2x(n-3)%o%val + arr_der2x(n-4)%o%val )*h_inv2/12.0
          BF(n) = (145.0*arr_der2x(n)%f%val - 304.0*arr_der2x(n-1)%f%val + 174.0*arr_der2x(n-2)%f%val  - &
               16.0*arr_der2x(n-3)%f%val + arr_der2x(n-4)%f%val )*h_inv2/12.0


          !Gaussian Elimination for the RHS matrices alone now.
          
          do k=2,n
             BT(k) = BT(k) - BT(k-1)*b2X(n+1-k)/a2X(k-1)
             BO(k) = BO(k) - BO(k-1)*b2X(n+1-k)/a2X(k-1)
             BF(k) = BF(k) - BF(k-1)*b2X(n+1-k)/a2X(k-1)             
          end do
          
          !Solving backwards.. Final step

          arr_der2x(n)%t%der2x = BT(n)/a2X(n)
          arr_der2x(n)%o%der2x = BO(n)/a2X(n)
          arr_der2x(n)%f%der2x = BF(n)/a2X(n)
          
          do k=n-1,1,-1
             arr_der2x(k)%t%der2x = (BT(k) - b2X(k)*arr_der2x(k+1)%t%der2x)/a2X(k)
             arr_der2x(k)%o%der2x = (BO(k) - b2X(k)*arr_der2x(k+1)%o%der2x)/a2X(k)
             arr_der2x(k)%f%der2x = (BF(k) - b2X(k)*arr_der2x(k+1)%f%der2x)/a2X(k)
          end do


          deallocate(BT, BO, BF)


        end subroutine deriv2_x


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        subroutine deriv2_y(arr_der2y, h_inv2, n)
!This function fixes the value of the second 'y' derivative of the given array. The derivative is calculated using a Tridiagonal scheme. It is 6 order accurate at the centre and 4th accurate at the boundaries

!    Inputs:
!    arr_der2y --> The array for which the second 'y' derivative is to be found. The input array is of type TFO. Since the arrays are passed by reference, the derivatives are stored in the array itself in an appropriate manner.
!    n --> The length of the input array
!    h_inv2 --> (1/h)^2 where h is the space discretization
          
          !!Declaration of input variables
          type(TFO), dimension(:)  :: arr_der2y
          real, intent(in) :: h_inv2
          integer, intent(in) :: n

          !!Variables declared within the function
          !Iteration variables
          integer :: k
          !The coefficients of the derivatives
          real alpha, a, b
          !The RHS matrices used to solve the tridiagonal equations. Made allocatable and later made of dimension 'n' --> the input array size
          real*8, dimension(:), allocatable :: BT, BO, BF

          allocate(BT(n))
          allocate(BO(n))
          allocate(BF(n))

!!! choice of alpha
          alpha = 2.0/11.0
    
!!Deciding other parameters from alpha
!!beta = 0.. tridiagonal scheme
          a = 4.0*(1.0 - alpha)/3.0
          b = (10.0*alpha - 1.0)/3.0
!c = 0
          
!Calculating the RHS matrices
          !Inner nodes
          BT(3:n-2) = ( b*(arr_der2y(5:)%t%val - 2.0*arr_der2y(3:n-2)%t%val + arr_der2y(:n-4)%t%val)*0.25 + &
               a*(arr_der2y(4:n-1)%t%val  - 2.0*arr_der2y(3:n-2)%t%val + arr_der2y(2:n-3)%t%val) )*h_inv2
          BO(3:n-2) = ( b*(arr_der2y(5:)%o%val - 2.0*arr_der2y(3:n-2)%o%val + arr_der2y(:n-4)%o%val)*0.25 + &
               a*(arr_der2y(4:n-1)%o%val  - 2.0*arr_der2y(3:n-2)%o%val + arr_der2y(2:n-3)%o%val) )*h_inv2
          BF(3:n-2) = ( b*(arr_der2y(5:)%f%val - 2.0*arr_der2y(3:n-2)%f%val + arr_der2y(:n-4)%f%val)*0.25 + &
               a*(arr_der2y(4:n-1)%f%val  - 2.0*arr_der2y(3:n-2)%f%val + arr_der2y(2:n-3)%f%val) )*h_inv2

          !Penultimate nodes 4th order accurate
          BT(2) = 1.2*(arr_der2y(3)%t%val - 2.0*arr_der2y(2)%t%val + arr_der2y(1)%t%val)*h_inv2
          BT(n-1) =  1.2*(arr_der2y(n)%t%val - 2.0*arr_der2y(n-1)%t%val + arr_der2y(n-2)%t%val)*h_inv2
          BO(2) =  1.2*(arr_der2y(3)%o%val - 2.0*arr_der2y(2)%o%val + arr_der2y(1)%o%val)*h_inv2
          BO(n-1) =  1.2*(arr_der2y(n)%o%val - 2.0*arr_der2y(n-1)%o%val + arr_der2y(n-2)%o%val)*h_inv2
          BF(2) =  1.2*(arr_der2y(3)%f%val - 2.0*arr_der2y(2)%f%val + arr_der2y(1)%f%val)*h_inv2
          BF(n-1) =  1.2*(arr_der2y(n)%f%val - 2.0*arr_der2y(n-1)%f%val + arr_der2y(n-2)%f%val)*h_inv2

          !Boundary nodes 4th order accurate
          BT(1) = (145.0*arr_der2y(1)%t%val - 304.0*arr_der2y(2)%t%val + 174.0*arr_der2y(3)%t%val  - &
               16.0*arr_der2y(4)%t%val + arr_der2y(5)%t%val)*h_inv2/12.0
          BO(1) = (145.0*arr_der2y(1)%o%val - 304.0*arr_der2y(2)%o%val + 174.0*arr_der2y(3)%o%val  - &
               16.0*arr_der2y(4)%o%val + arr_der2y(5)%o%val)*h_inv2/12.0
          BF(1) = (145.0*arr_der2y(1)%f%val - 304.0*arr_der2y(2)%f%val + 174.0*arr_der2y(3)%f%val  - &
               16.0*arr_der2y(4)%f%val + arr_der2y(5)%f%val)*h_inv2/12.0

          BT(n) = (145.0*arr_der2y(n)%t%val - 304.0*arr_der2y(n-1)%t%val + 174.0*arr_der2y(n-2)%t%val  - &
               16.0*arr_der2y(n-3)%t%val + arr_der2y(n-4)%t%val )*h_inv2/12.0
          BO(n) = (145.0*arr_der2y(n)%o%val - 304.0*arr_der2y(n-1)%o%val + 174.0*arr_der2y(n-2)%o%val  - &
               16.0*arr_der2y(n-3)%o%val + arr_der2y(n-4)%o%val )*h_inv2/12.0
          BF(n) = (145.0*arr_der2y(n)%f%val - 304.0*arr_der2y(n-1)%f%val + 174.0*arr_der2y(n-2)%f%val  - &
               16.0*arr_der2y(n-3)%f%val + arr_der2y(n-4)%f%val )*h_inv2/12.0


          !Gaussian Elimination for the RHS matrices alone now.
          
          do k=2,n
             BT(k) = BT(k) - BT(k-1)*b2Y(n+1-k)/a2Y(k-1)
             BO(k) = BO(k) - BO(k-1)*b2Y(n+1-k)/a2Y(k-1)
             BF(k) = BF(k) - BF(k-1)*b2Y(n+1-k)/a2Y(k-1)             
          end do
          
          !Solving backwards.. Final step

          arr_der2y(n)%t%der2y = BT(n)/a2Y(n)
          arr_der2y(n)%o%der2y = BO(n)/a2Y(n)
          arr_der2y(n)%f%der2y = BF(n)/a2Y(n)

          do k=n-1,1,-1
             arr_der2y(k)%t%der2y = (BT(k) - b2Y(k)*arr_der2y(k+1)%t%der2y)/a2Y(k)
             arr_der2y(k)%o%der2y = (BO(k) - b2Y(k)*arr_der2y(k+1)%o%der2y)/a2Y(k)
             arr_der2y(k)%f%der2y = (BF(k) - b2Y(k)*arr_der2y(k+1)%f%der2y)/a2Y(k)
          end do

          deallocate(BT, BO, BF)


        end subroutine deriv2_y

!subroutine to calculate the absolute time in seconds
        subroutine ft(abstime)
          character*8 date
          character*10 time
          character*5 zone
          INTEGER*4 VALUES(8)
          real abstime
          call date_and_time( date, time, zone, values )
          abstime=values(5)*3600+values(6)*60+values(7)+values(8)*0.0001
        end subroutine ft



!function for caluating the time derivative 

 subroutine fd(arr_dt,fo)
 ! calculates time derivative
  use data_struct  

  type(TFO), dimension(:), intent(in) :: arr_dt !The current time array
  type(TFO), dimension(:),intent(out)::fo   !function ar arr_st
 
   !!local variables
  REAL*8 :: wr
  !Iteration variables
  integer :: i,j
 ! allocate(fo(size(arr_dt)))

  do j = 1,ny-2
     do i = 2, nx-1
        
        !Calculating reaction rate
        wr = (beta*arr_dt(j*nx+i)%f%val)*(beta*arr_dt(j*nx+i)%o%val)*(Da*beta*&
             exp(beta*(arr_dt(j*nx+i)%t%val - 1.0)*(1+gamma)/(1.0 + gamma*arr_dt(j*nx+i)%t%val)))
        
        !Updating t values
        fo(j*nx+i)%t%val =(arr_dt(j*nx+i)%t%der2x + arr_dt(j*nx+i)%t%der2y - arr_dt(j*nx+i)%t%derx + (1.0 + phi)*wr )
        
        !Updating o values
        fo(j*nx+i)%o%val =(LeO_inv*(arr_dt(j*nx+i)%o%der2x + arr_dt(j*nx+i)%o%der2y) -  arr_dt(j*nx+i)%o%derx - phi*wr )
        
        !Updating f values
        fo(j*nx+i)%f%val =(LeF_inv*(arr_dt(j*nx+i)%f%der2x + arr_dt(j*nx+i)%f%der2y) -  arr_dt(j*nx+i)%f%derx - wr )
!        write(*,*)i," ",j
     end do
  end do
!omp parallel end do  

!deallocate(fo)
end subroutine fd
 subroutine fdamp(arr_dt,arr_dt1,fo)
 
 ! adds damping to residual arrray
  use data_struct  

  type(TFO), dimension(:), intent(in) :: arr_dt,arr_dt1 !The current time array and filtered array
  type(TFO), dimension(:),intent(out)::fo   !residual aray 
   !!local variables
  REAL*8 :: wr
  type(tfo1):: qdiff
  !Iteration variables
  integer :: i,j
 ! allocate(fo(size(arr_dt)))
      ! write(*,*)"value of chi si ",chi

  do j = 1,ny-2
     do i = 2, nx-1
      !Calculating fictitious damping reaction rate   based on difference                                                                                                                                                                              
       
        qdiff%t=(arr_dt(j*nx+i)%t%val-arr_dt1(j*nx+i)%t%val)
        qdiff%f=(arr_dt(j*nx+i)%f%val-arr_dt1(j*nx+i)%f%val)
        qdiff%o=(arr_dt(j*nx+i)%o%val-arr_dt1(j*nx+i)%o%val)
       ! wr = (Da*beta*&
       !      exp(beta*( abs(qdiff%t) - 1.0)*(1+gamma)/(1.0 + gamma*abs(qdiff%t))))
        !Updating t values
        fo(j*nx+i)%t%val =  fo(j*nx+i)%t%val-chi*qdiff%t
        
        !Updating o values
        fo(j*nx+i)%o%val =fo(j*nx+i)%o%val-chi*qdiff%o
        
        !Updating f values
        fo(j*nx+i)%f%val =fo(j*nx+i)%f%val-chi*qdiff%f
        !        write(*,*)i," ",j
     end do
  end do
  !omp parallel end do  
  
  !deallocate(fo)
end subroutine fdamp

 subroutine fsdamp(arr_dt,arr_unfiltered,fo)
 
!calculates the filtered variable used for time marching of filtered variable
  use data_struct  

  type(TFO), dimension(:), intent(in) :: arr_unfiltered,arr_dt !The current time array and filtered array
  type(TFO), dimension(:),intent(out)::fo   !residual aray 
   !!local variables
  REAL*8 :: wr
  !Iteration variables
  integer :: i,j
 ! allocate(fo(size(arr_dt)))
wr=1.0/delta
  do j = 1,ny-2
     do i = 2, nx-1
        !Updating t values
        fo(j*nx+i)%t%val =  wr*(arr_unfiltered(j*nx+i)%t%val-arr_dt(j*nx+i)%t%val)
        
        !Updating o values
        fo(j*nx+i)%o%val =wr*(arr_unfiltered(j*nx+i)%o%val-arr_dt(j*nx+i)%o%val)
        
        !Updating f values
        fo(j*nx+i)%f%val =wr*(arr_unfiltered(j*nx+i)%f%val-arr_dt(j*nx+i)%f%val)
!        write(*,*)i," ",j
     end do
  end do
!omp parallel end do  

!deallocate(fo)
end subroutine fsdamp


function fd1(arr_dt) result(fo)
 ! calculates time derivative
  use data_struct
  type(TFO), dimension(:), intent(in) :: arr_dt !The current time array
  type(TFO), dimension(:),allocatable::fo   !function ar arr_st

   !!local variables
  REAL*8      :: wr
  !Iteration variables
  integer :: i,j
  allocate(fo(size(arr_dt)))

  do j = 0,ny-1
     do i = 1, nx

        !Calculating reaction rate
        wr = (beta*arr_dt(j*nx+i)%f%val)*(beta*arr_dt(j*nx+i)%o%val)*(Da*beta*&
             exp(beta*(arr_dt(j*nx+i)%t%val - 1.0)*(1+gamma)/(1.0 + gamma*arr_dt(j*nx+i)%t%val)))

        !Updating t values
        fo(j*nx+i)%t%val =(arr_dt(j*nx+i)%t%der2x + arr_dt(j*nx+i)%t%der2y - arr_dt(j*nx+i)%t%derx + (1.0 + phi)*wr )

        !Updating o values
        fo(j*nx+i)%o%val =(LeO_inv*(arr_dt(j*nx+i)%o%der2x + arr_dt(j*nx+i)%o%der2y) -  arr_dt(j*nx+i)%o%derx - phi*wr )

        !Updating f values
        fo(j*nx+i)%f%val =(LeF_inv*(arr_dt(j*nx+i)%f%der2x + arr_dt(j*nx+i)%f%der2y) -  arr_dt(j*nx+i)%f%derx - wr )
!        write(*,*)i," ",j
     end do
  end do
!omp parallel end do

!deallocate(fo)
end function fd1




subroutine deriv_boundary(arr_dt)
 type(TFO), dimension(:), intent(inout) :: arr_dt !The current time array
integer:: i,k
!y derivatives
!$omp parallel shared(arr_dt) default(shared)
!$omp do private(i) schedule(static,25)
         do i = 2,nx-1 
            call deriv2_y(arr_dt(i:(ny-1)*nx+i:nx), hy_inv2, ny)
         end do
!$omp end do 

         !x derivatives

!$omp do private(i) schedule(static,25)
         do i = 1,ny-2            
            call deriv_x(arr_dt(i*nx+1:nx*(i+1)), hx_inv, nx)
            call deriv2_x(arr_dt(i*nx+1:nx*(i+1)), hx_inv2, nx)
         end do
!$omp end do
!$omp end parallel 
        call boundary(arr_dt(1:nx*ny))
!end do


end subroutine deriv_boundary

!calcukates reaction rate from array
subroutine rrate_calc(arr_dt,w)
 ! calculates time derivative
  use data_struct  
  type(TFO), dimension(:), intent(in) :: arr_dt !The current time array
  real*8, dimension(:),intent(out)::w   !function ar arr_st

  !Iteration variables
  integer :: i,j

  do j = 1,ny*nx
     !Calculating reaction rate
     w(j) = (beta*arr_dt(j)%f%val)*(beta*arr_dt(j)%o%val)*(Da*beta*&
          exp(beta*(arr_dt(j)%t%val - 1.0)*(1+gamma)/(1.0 + gamma*arr_dt(j)%t%val)))
  end do
     
end subroutine rrate_calc

subroutine rrate_max(arr_dt,wmax,wloc,wsum)
 ! calculates time derivative
  use data_struct
  type(TFO), dimension(:), intent(in) :: arr_dt !The current time array
  real*8, intent(out)::wmax,wsum   !function ar arr_st
  integer,intent(out)::wloc
  !Iteration variables
  integer :: i,j
  real*8::wr
  wmax=0.0
  wloc=0
  wsum=0.0
  do j = 1,ny*nx
     !Calculating reaction rate
     wr = (beta*arr_dt(j)%f%val)*(beta*arr_dt(j)%o%val)*(Da*beta*&
          exp(beta*(arr_dt(j)%t%val - 1.0)*(1+gamma)/(1.0 + gamma*arr_dt(j)%t%val)))
     wsum=wsum+wr
     if (wr .gt. wmax)then
        wmax=wr
        wloc=j  
     end if
  end do
  
end subroutine rrate_max




end program diff_flame

      !      subroutine OMP_set_num_threads(number_of_threads)
 !       integer(kind = OMP_integer_kind), intent(in) :: number_of_threads
  !    end subroutine OMP_set_num_threads

