module setup

  !
  ! Possible values for n < 5000
  !
  ! 0,1,2,3,4,5,6,7,8,9,10,12,14,15,16,18,20,21,24,25,27,28,30,32,35
  ! ,36,40,42,45,48,49,50,54,56,60,63,64,70,72,75,80,81,84,90,96,98
  ! ,100,105,108,112,120,125,126,128,135,140,144,147,150,160,162,168
  ! ,175,180,189,192,196,200,210,216,224,225,240,243,245,250,252,256
  ! ,270,280,288,294,300,315,320,324,336,343,350,360,375,378,384,392
  ! ,400,405,420,432,441,448,450,480,486,490,500,504,512,525,540,560
  ! ,567,576,588,600,625,630,640,648,672,675,686,700,720,729,735,750
  ! ,756,768,784,800,810,840,864,875,882,896,900,945,960,972,980,1000
  ! ,1008,1024,1029,1050,1080,1120,1125,1134,1152,1176,1200,1215,1225
  ! ,1250,1260,1280,1296,1323,1344,1350,1372,1400,1440,1458,1470,1500
  ! ,1512,1536,1568,1575,1600,1620,1680,1701,1715,1728,1750,1764,1792
  ! ,1800,1875,1890,1920,1944,1960,2000,2016,2025,2048,2058,2100,2160
  ! ,2187,2205,2240,2250,2268,2304,2352,2400,2401,2430,2450,2500,2520
  ! ,2560,2592,2625,2646,2688,2700,2744,2800,2835,2880,2916,2940,3000
  ! ,3024,3072,3087,3125,3136,3150,3200,3240,3360,3375,3402,3430,3456
  ! ,3500,3528,3584,3600,3645,3675,3750,3780,3840,3888,3920,3969,4000
  ! ,4032,4050,4096,4116,4200,4320,4374,4375,4410,4480,4500,4536,4608
  ! ,4704,4725,4800,4802,4860,4900,5000
  !

  use precision
  use mpi
  use globals
  use conversion
  use misc
  use fft

  implicit none


  
contains



  subroutine s_load_data()

    character(len=20) :: object
    character(len=256) :: func
    character(len=maxstringlength) :: inputfilename
    real(fdp), dimension(dim) :: t_point1, t_point2
    real(fdp) :: t_real1, t_real2, t_real3
    integer(idp) :: t_int1
    real(fdp), dimension(2) :: potrange
    real(fdp), dimension(4) :: potbound

    
#ifdef DEBUG
    open(unit=potlog, file='log.input', action='write')
#endif

    call getarg(1, inputfilename)

    open(unit=inputfile, file=trim(inputfilename), action='read',&
         & form='formatted')

    potrange = (/ 0.0d0, 0.0d0 /)
    potbound = (/ 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

    do

       !
       ! Find out what the next input object it
       !

       call s_read(object=object)

#ifdef DEBUG
       if (trim(object) /= '#') then
          
          write(potlog, *) 'Loading object ', trim(object)

       end if
#endif

       select case (object)

          !
          ! Load in the string parameters
          !

       case ('datadir')
          call s_read(string=parameters%datadir)

          defined_parameters%datadir = .true.

       case ('chkpntfile')
          call s_read(string=parameters%chkpntfile)

          defined_parameters%chkpntfile = .true.

       case ('basename')
          call s_read(string=parameters%basename)

          defined_parameters%basename = .true.

       case('initialwavefn')
          call s_read(string=parameters%initialwavefn)

          defined_parameters%initialwavefn = .true.

       case('optimisation')
          call s_read(string=parameters%optimisation)

          defined_parameters%optimisation = .true.
          defined_parameters%potential = .true.

       case ('potsmoothtype')
          call s_read(string=parameters%potsmoothtype)

          defined_parameters%potsmoothtype = .true.

       case('gauge')
          call s_read(string=parameters%gauge)

          defined_parameters%gauge = .true.

          !
          ! Load in the single real parameters
          !

       case ('mass')
          call s_read(realnum=parameters%mass)

          defined_parameters%mass = .true.

       case ('fermienergy')
          call s_read(realnum=parameters%fermienergy)

          defined_parameters%fermienergy = .true.

       case ('temperature')
          call s_read(realnum=parameters%temperature)

          defined_parameters%temperature = .true.

       case ('initialenergy')
          call s_read(realnum=parameters%initialenergy)

          defined_parameters%initialenergy = .true.

       case ('t')
          call s_read(realnum=parameters%t)

          defined_parameters%t = .true.

       case ('tstep')
          call s_read(realnum=parameters%tstep)

          defined_parameters%tstep = .true.

       case ('distance')
          call s_read(realnum=parameters%distance)

          defined_parameters%distance = .true.

       case ('speccutoff')
          call s_read(realnum=parameters%speccutoff)

          defined_parameters%speccutoff = .true.

       case ('externale')
          call s_read(realnum=parameters%externale)
          
          defined_parameters%externale = .true.
          defined_parameters%potential = .true.

       case ('externalb')
          call s_read(realnum=parameters%externalb)

          defined_parameters%externalb = .true.
          defined_parameters%potential = .true.

       case ('omegal')
          call s_read(realnum=parameters%omegal)

          defined_parameters%omegal = .true.

       case('eigenfactor')
          call s_read(realnum=parameters%eigenfactor)

          defined_parameters%eigenfactor = .true.

          !
          ! Load in the single integer parameters
          !

       case ('timesteps')
          call s_read(intnum=parameters%timesteps)

          defined_parameters%timesteps = .true.

       case ('chkpntfreq')
          call s_read(intnum=parameters%chkpntfreq)

          defined_parameters%chkpntfreq = .true.

       case ('plotnumber')
          call s_read(intnum=parameters%plotnumber)

          defined_parameters%plotnumber = .true.

          !
          ! Load in the single logical parameters
          !

       case ('plotcomplex')
          call s_read(logic=parameters%plotcomplex)

          defined_parameters%plotcomplex = .true.

       case ('plotwfn')
          call s_read(logic=parameters%plotwfn)

          defined_parameters%plotwfn = .true.

       case ('plotwfnp')
          call s_read(logic=parameters%plotwfnp)

          defined_parameters%plotwfnp = .true.

       case ('plotpot')
          call s_read(logic=parameters%plotpot)

          defined_parameters%plotpot = .true.

       case ('tranrefl')
          call s_read(logic=parameters%tranrefl)

          defined_parameters%tranrefl = .true.

       case ('espectrum')
          call s_read(logic=parameters%espectrum)

          defined_parameters%espectrum = .true.

       case ('pspectrum')
          call s_read(logic=parameters%pspectrum)

          defined_parameters%pspectrum = .true.

       case ('tspectrum')
          call s_read(logic=parameters%tspectrum)

          defined_parameters%tspectrum = .true.

       case ('conductance')
          call s_read(logic=parameters%conductance)

          defined_parameters%conductance = .true.

       case('primefactorn')
          call s_read(logic=parameters%primefactorn)

          defined_parameters%primefactorn = .true.

          !
          ! Load in multi dimensional real parameters
          !
          
       case ('p')
          call s_read(realnums=parameters%p)
          
          defined_parameters%p = .true.
          
       case ('pspread')
          call s_read(realnums=parameters%pspread)
          
          defined_parameters%pspread = .true.

       case ('pmax')
          call s_read(realnums=parameters%pmax)

          defined_parameters%pmax = .true.

       case ('accuracyfactor')
          call s_read(realnums=parameters%accuracyfactor)

          defined_parameters%accuracyfactor = .true.

       case ('position')
          call s_read(realnums=parameters%position)

          defined_parameters%position = .true.

       case ('wavefnwidth')
          call s_read(realnums=parameters%wavefnwidth)

          defined_parameters%wavefnwidth = .true.

       case ('gridspacing')
          call s_read(realnums=parameters%gridspacing)
        
          defined_parameters%gridspacing = .true.

       case ('maxgridspacing')
          call s_read(realnums=parameters%maxgridspacing)
        
          defined_parameters%maxgridspacing = .true.

       case ('spacesize')
          call s_read(realnums=parameters%spacesize)
          
          defined_parameters%spacesize = .true.
          
       case ('pspacing')
          call s_read(realnums=parameters%pspacing)

          defined_parameters%pspacing = .true.

       case ('potrange')
          call s_read(realnums=parameters%potrange)
          
          defined_parameters%potrange = .true.
          
       case ('potbound')
          call s_read(realnums=parameters%potbound)
          
          defined_parameters%potbound = .true.
          
          !
          ! Load in multi dimensional integer parameters
          !

       case ('smooth')
          call s_read(intnums=parameters%smooth)

          defined_parameters%smooth = .true.

       case ('n')
          call s_read(intnums=parameters%n)
          
          defined_parameters%n = .true.
          
       case ('plotpointnum')
          call s_read(intnums=parameters%plotpointnum)
          
          defined_parameters%plotpointnum = .true.

       case ('extraboundary')
          call s_read(intnums=parameters%extraboundary)
          
          defined_parameters%extraboundary = .true.
          
          !
          ! Load in multi dimensional logical parameters
          !

       case('comppotpos')
          call s_read(logics=parameters%comppotpos)

          defined_parameters%comppotpos = .true.

          !
          ! Load in potential components
          !

       case('potfn')

          !
          ! Read a functional potential
          !

          call s_read(intnum=t_int1)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)
          call s_read(string=func)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, a, a)', advance='no') object, t_int1&
               & , ' ', func

          potrange(1) = min(t_real1, potrange(1))
          potrange(2) = max(t_real2, potrange(2))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('circle')

          !
          ! Read a circle
          !

          call s_read(intnum=t_int1)
          call s_read(realnums=t_point1)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)
          call s_read(realnum=t_real3)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, 2e, e, e, e)', advance='no') object,&
               & t_int1, t_point1, t_real1, t_real2, t_real3

          potrange(1) = min(t_real3, potrange(1))
          potrange(2) = max(t_real3, potrange(2))

          potbound(1) = min(t_point1(1)-t_real1-t_real2, potbound(1))
          potbound(2) = max(t_point1(1)+t_real1+t_real2, potbound(2))
          potbound(3) = min(t_point1(2)-t_real1-t_real2, potbound(3))
          potbound(4) = max(t_point1(2)+t_real1+t_real2, potbound(4))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('fillcircle')
          call s_read(intnum=t_int1)
          call s_read(realnums=t_point1)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, 2e, e, e)', advance='no') object,&
               & t_int1, t_point1, t_real1, t_real2

          potrange(1) = min(t_real2, potrange(1))
          potrange(2) = max(t_real2, potrange(2))

          potbound(1) = min(t_point1(1)-t_real1, potbound(1))
          potbound(2) = max(t_point1(1)+t_real1, potbound(2))
          potbound(3) = min(t_point1(2)-t_real1, potbound(3))
          potbound(4) = max(t_point1(2)+t_real1, potbound(4))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('ellipse')
          call s_read(intnum=t_int1)
          call s_read(realnums=t_point1)
          call s_read(realnums=t_point2)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, 2e, 2e, e, e)', advance='no') object&
               & , t_int1, t_point1, t_point2, t_real1, t_real2

          potrange(1) = min(t_real2, potrange(1))
          potrange(2) = max(t_real2, potrange(2))

          potbound(1) = min(t_point1(1)-t_point2(1)-t_real1,&
               & potbound(1))
          potbound(2) = max(t_point1(1)+t_point2(1)+t_real1,&
               & potbound(2))
          potbound(3) = min(t_point1(2)-t_point2(2)-t_real1,&
               & potbound(3))
          potbound(4) = max(t_point1(2)+t_point2(2)+t_real1,&
               & potbound(4))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('fillellipse')
          call s_read(intnum=t_int1)
          call s_read(realnums=t_point1)
          call s_read(realnums=t_point2)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, 2e, 2e, e, e)', advance='no') object&
               & , t_int1, t_point1, t_point2, t_real1, t_real2

          potrange(1) = min(t_real2, potrange(1))
          potrange(2) = max(t_real2, potrange(2))

          potbound(1) = min(t_point1(1)-t_point2(1)-t_real1,&
               & potbound(1))
          potbound(2) = max(t_point1(1)+t_point2(1)+t_real1,&
               & potbound(2))
          potbound(3) = min(t_point1(2)-t_point2(2)-t_real1,&
               & potbound(3))
          potbound(4) = max(t_point1(2)+t_point2(2)+t_real1,&
               & potbound(4))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('line')
          call s_read(intnum=t_int1)
          call s_read(realnums=t_point1)
          call s_read(realnums=t_point2)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, 2e, 2e, e, e)', advance='no') object&
               & , t_int1, t_point1, t_point2, t_real1, t_real2

          potrange(1) = min(t_real2, potrange(1))
          potrange(2) = max(t_real2, potrange(2))

          potbound(1) = min(t_point1(1)-t_real1, t_point2(1)-t_real1,&
               & potbound(1))
          potbound(2) = max(t_point1(1)+t_real1, t_point2(1)+t_real1,&
               & potbound(2))
          potbound(3) = min(t_point1(2)-t_real1, t_point2(2)-t_real1,&
               & potbound(3))
          potbound(4) = max(t_point1(2)+t_real1, t_point2(2)+t_real1,&
               & potbound(4))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('ramp')
          call s_read(intnum=t_int1)
          call s_read(realnums=t_point1)
          call s_read(realnums=t_point2)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)
          call s_read(realnum=t_real3)

          defined_parameters%potential = .true.

          write(potdata, '(a, i, 2e, 2e, e, e, e)', advance='no')&
               & object, t_int1, t_point1, t_point2, t_real1, t_real2&
               & , t_real3

          potrange(1) = min(t_real2, t_real3, potrange(1))
          potrange(2) = max(t_real2, t_real3, potrange(2))

          potbound(1) = min(t_point1(1)-t_real1, t_point2(1)-t_real1,&
               & potbound(1))
          potbound(2) = max(t_point1(1)+t_real1, t_point2(1)+t_real1,&
               & potbound(2))
          potbound(3) = min(t_point1(2)-t_real1, t_point2(2)-t_real1,&
               & potbound(3))
          potbound(4) = max(t_point1(2)+t_real1, t_point2(2)+t_real1,&
               & potbound(4))

          parameters%potcomps(t_int1) = .true.

          write(potdata, '(a)') ''

       case ('polygon')
          call s_read(intnum=t_int1)
          call s_read(realnum=t_real1)
          call s_read(realnum=t_real2)
          call s_read(realnums=t_point1)

          defined_parameters%potential = .true.

          potrange(1) = min(t_real2, potrange(1))
          potrange(2) = max(t_real2, potrange(2))

          write(potdata, '(a, i, e, e, 2e)', advance='no') object,&
               & t_int1, t_real1, t_real2, t_point1

          call s_read(realnums=t_point2)
          write(potdata, '(2e)', advance='no') t_point2
          
          potbound(1) = min(t_point2(1)-t_real1,  potbound(1))
          potbound(2) = max(t_point2(1)+t_real1,  potbound(2))
          potbound(3) = min(t_point2(2)-t_real1,  potbound(3))
          potbound(4) = max(t_point2(2)+t_real1,  potbound(4))
          
          do while(sum(abs(t_point1-t_point2))/minval(parameters &
               &%spacesize) > cutoff)
             
             call s_read(realnums=t_point2)
             write(potdata, '(2e)', advance='no') t_point2

             potbound(1) = min(t_point2(1)-t_real1,  potbound(1))
             potbound(2) = max(t_point2(1)+t_real1,  potbound(2))
             potbound(3) = min(t_point2(2)-t_real1,  potbound(3))
             potbound(4) = max(t_point2(2)+t_real1,  potbound(4))

          end do

          write(potdata, '(a)') ''

          parameters%potcomps(t_int1) = .true.

       case ('fillpolygon')
          call s_read(intnum=t_int1)
          call s_read(realnum=t_real1)
          call s_read(realnums=t_point1)

          defined_parameters%potential = .true.

          potrange(1) = min(t_real1, potrange(1))
          potrange(2) = max(t_real1, potrange(2))

          write(potdata, '(a, i, e, 2e)', advance='no') object,&
               & t_int1, t_real1, t_point1
          
          call s_read(realnums=t_point2)
          write(potdata, '(2e)', advance='no') t_point2

          potbound(1) = min(t_point2(1),  potbound(1))
          potbound(2) = max(t_point2(1),  potbound(2))
          potbound(3) = min(t_point2(2),  potbound(3))
          potbound(4) = max(t_point2(2),  potbound(4))

          do while(sum(abs(t_point1-t_point2))/minval(parameters &
               &%spacesize) > cutoff)

             call s_read(realnums=t_point2)
             write(potdata, '(2e)', advance='no') t_point2

             potbound(1) = min(t_point2(1),  potbound(1))
             potbound(2) = max(t_point2(1),  potbound(2))
             potbound(3) = min(t_point2(2),  potbound(3))
             potbound(4) = max(t_point2(2),  potbound(4))

          end do

          write(potdata, '(a)') ''

          parameters%potcomps(t_int1) = .true.

       case ('end')
          exit

       case default

       end select

#ifdef DEBUG
       if (trim(object) /= '#') then
          
          write(potlog, *) 'Finished object ', trim(object)

       end if
#endif

    end do

#ifdef DEBUG
    close(potlog)
#endif

    if (.not. defined_parameters%potrange) then

       parameters%potrange = potrange

    end if

    if (.not. defined_parameters%potbound) then

       parameters%potbound = potbound

    end if

    close(inputfile)

  end subroutine s_load_data



  subroutine s_derive_parameters()

    integer(isp), dimension(2*dim) :: potposition
    real(fdp), dimension(dim) :: tstepmin


    where (parameters%comppotpos)
       
       potposition = 1
       
    end where
    
    if (.not. defined_parameters%extraboundary) then
       
       select case(parameters%optimisation)
          
       case ('exponential', 'gaussian', 'linear')
          parameters%extraboundary = parameters%extraboundary + (/&
               & potposition(1)+potposition(2), potposition(3) &
               &+potposition(4) /)
          
       case default
          
       end select
       
    end if
    
    if (.not. defined_parameters%p) then
       
       parameters%p = (/ sqrt(2.0d0*parameters%mass*parameters &
            &%initialenergy), 0.0d0 /)
       
    end if
    
    if (.not. defined_parameters%initialenergy) then
       
       parameters%initialenergy = sum(parameters%p**2)/2.0d0 &
            &/parameters%mass
       
    end if

    parameters%pspread = parameters%pspread*sqrt(sum(parameters%p**2))
       
    if (.not. defined_parameters%omegal) then

       parameters%omegal = abs(chargee*parameters%externalb/2.0d0 &
            &/parameters%mass)

    end if

    if (.not. defined_parameters%wavefnwidth) then

       if (defined_parameters%externalb .and. (.not.&
            & defined_parameters%pspread)) then
          
          parameters%wavefnwidth = sqrt(hbar/parameters%mass &
               &/parameters%omegal)

       else
          
          parameters%wavefnwidth = hbar/parameters%pspread

       end if

    end if

    if (.not. defined_parameters%pmax) then

       !
       ! compute pmax from the width of the fermi or gauss packet and
       ! the depth of the potential well
       !

       parameters%pmax = abs(parameters%p)
       
       select case(parameters%initialwavefn)

       case ('gaussian')
          parameters%pmax = parameters%pmax + sqrt(-2.0d0*log(cutoff)&
               & /parameters%wavefnwidth**2*hbar**2)
          
       case ('fermi')

          !
          ! The fermi-dirac function is roughly bound by a gaussian
          ! with width
          ! w=sqrt(hbar^2/bk/m/t)
          !

          parameters%pmax = parameters%pmax + sqrt(-2.0d0*log(cutoff)&
               & *parameters%mass*parameters%temperature*boltzk)

       case default

       end select
       
       parameters%pmax = parameters%pmax + sqrt(2.0d0*parameters%mass&
            & *(abs(parameters%potrange(2)-parameters%potrange(1))&
            & +abs(parameters %externale)))
       
       parameters%pmax = parameters%pmax * parameters%accuracyfactor

       parameters%pmax = maxval(parameters%pmax)
       
    end if

    if (.not. defined_parameters%gridspacing) then
       
       parameters%gridspacing = f_gridspacing(parameters%pmax)

    end if

    if (defined_parameters%maxgridspacing) then

       parameters%gridspacing = min(parameters%gridspacing,&
            & parameters%maxgridspacing)

    end if

    if (.not. defined_parameters%n) then

       parameters%n = ceiling(dble(int(parameters%spacesize &
            &/parameters%gridspacing)+1+parameters%extraboundary &
            &*parameters%smooth)/2.0d0)*2

       if (parameters%primefactorn) then
          
          parameters%n(1) = f_factorprime(parameters%n(1))
          parameters%n(2) = f_factorprime(parameters%n(2))

       end if
       
    end if

    if (.not. defined_parameters%gridspacing) then

       parameters%gridspacing = parameters%spacesize/dble(parameters &
            &%n-1-parameters%extraboundary*parameters%smooth)

    end if

    if (.not. defined_parameters%pmax) then

       parameters%pmax = f_pmax(parameters%gridspacing)

    end if

    if (.not. defined_parameters%pspacing) then
       
       !
       ! Must align the grid to the FFT computed momentum
       ! thus the n-1 term is n only
       !

       parameters%pspacing = 2.0d0*parameters%pmax / dble(parameters &
            &%n)

    end if

    if (.not. defined_parameters%t) then

       parameters%t = abs(parameters%distance)*parameters%mass &
            &/parameters%p(1)

    end if

    if (.not. defined_parameters%tstep) then

       select case(parameters%optimisation)

       case ('exponential', 'gaussian', 'linear')

          tstepmin = parameters%t

          ! tstepmin = sqrt(2.0d0*parameters%mass/(parameters%pmax
          ! **2/2.0d0/parameters%mass)**(5.0d0/4.0d0)*dble(parameters
          ! %smooth)*parameters%gridspacing *(-parameters%wavefnwidth
          ! **2*parameters%pmax/hbar +sqrt(parameters%pmax**2
          ! *parameters%wavefnwidth**4/hbar**2 +2.0d0*log(2.0d0*pi
          ! *product(parameters%wavefnwidth))
          ! *parameters%wavefnwidth**2)))

          parameters%tstep = maxval(tstepmin)
          
          if (parameters%comppotpos(1) .or. parameters%comppotpos(2))&
               & then
             
             parameters%tstep = min(tstepmin(1), parameters%tstep)
             
          end if
          
          if (parameters%comppotpos(3) .or. parameters%comppotpos(4))&
               & then
             
             parameters%tstep = min(tstepmin(2), parameters%tstep)
             
          end if

       case default
          parameters%tstep = parameters%t

       end select

    end if

    if (.not. defined_parameters%timesteps) then
       
       parameters%timesteps = ceiling(parameters%t/parameters%tstep)
       
    end if

    if (.not. defined_parameters%tstep) then

       parameters%tstep = parameters%t/dble(parameters%timesteps)
       
    end if

    if (.not. defined_parameters%plotnumber) then
       
       parameters%plotnumber = parameters%timesteps
       
    end if

    if (.not. defined_parameters%plotpointnum) then

       parameters%plotpointnum = parameters%n

    end if

    if (defined_parameters%externalb .and. .not. defined_parameters&
         &%eigenfactor) then

       parameters%eigenfactor = parameters%eigenfactor*2.0d0

    end if
    
    if (.not. defined_parameters%pspread) then
       
      parameters%pspread = hbar/parameters%wavefnwidth
       
    end if

  end subroutine s_derive_parameters


  
  subroutine s_display_parameters()

    call s_display('mass= ', realnum=parameters%mass)
    call s_display('initialenergy= ', realnum=parameters &
         &%initialenergy)
    call s_display('p= ', realnums=parameters%p)
    call s_display('pspread= ', realnums=parameters%pspread)
    call s_display('pspread%= ', realnums=parameters%pspread&
         &/sqrt(sum(parameters%p**2))*100.0d0)
    call s_display('pmax= ', realnums=parameters%pmax)
    call s_display('gridspacing= ', realnums=parameters%gridspacing)
    call s_display('start position= ', realnums=parameters%position)
    call s_display('distance= ', realnum=parameters%distance)
    call s_display('pspacing= ', realnums=parameters%pspacing)
    call s_display('wavefnwidth= ', realnums=parameters%wavefnwidth)
    call s_display('temperature= ', realnum=parameters%temperature)
    call s_display('fermi energy= ', realnum=parameters%fermienergy)
    call s_display('external e= ',realnum=parameters%externale)
    call s_display('external b= ',realnum=parameters%externalb)
    call s_display('lomar frequency= ', realnum=parameters%omegal)
    call s_display('potrange= ', realnums=parameters%potrange)
    call s_display('t= ', realnum=parameters%t)
    call s_display('n= ', intnums=parameters%n)
    call s_display('potbound xmin, xmax= ', realnums=parameters &
         &%potbound(1:2))
    call s_display('potbound ymin, ymax= ', realnums=parameters &
         &%potbound(3:4))
    call s_display('space size= ', realnums=parameters%spacesize)
    call s_display('time steps= ', intnum=parameters%timesteps)
    call s_display('tstep= ', realnum=parameters%tstep)

  end subroutine s_display_parameters



  subroutine s_setup()
    
    call s_load_data()
    
    call s_derive_parameters()
    
    call s_display_parameters()
    
  end subroutine s_setup



  subroutine s_initialwavefn(wavefn, nl, nu)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(out) ::&
         & wavefn

    real(fdp), dimension(dim) :: t_wavefnwidth
    real(fdp) :: x, y
    integer(idp) :: i, j

    
    select case(parameters%initialwavefn)
       
    case ('gaussian')
 
       t_wavefnwidth = sqrt(-2.0d0*parameters%wavefnwidth**2&
            &*log(2.0d0*pi*parameters%wavefnwidth(1)*parameters&
            &%wavefnwidth(2)*cutoff))
       
       call s_gaussian(wavefn, nl, nu, t_wavefnwidth)
       
    case ('fermi')
       call s_fermidistribution(wavefn, nl, nu, parameters%position,&
            & parameters%p, parameters%gridspacing)

    case default

    end select

    !
    !  The phase needs to be adjusted for some of the gauges
    !  to give the same results as the (-y, x, 0) gauge
    !

    if (defined_parameters%externalb) then
       
       select case (parameters%gauge)
          
       case ('(-y,0,0)')
          do j = nl(2), nu(2)
             y = f_positiony(j) - parameters%position(2)
             do i = nl(1), nu(1)
                x = f_positionx(i) - parameters%position(1)
                wavefn(i, j) = wavefn(i, j) * exp(-ii*chargee*parameters%externalb*x*y/2.0d0/hbar)
             end do
          end do
          
       case ('(x-y,x-y,0)')
          do j = nl(2), nu(2)
             y = f_positiony(j) - parameters%position(2)
             do i = nl(1), nu(1)
                x = f_positionx(i) - parameters%position(1)
                wavefn(i, j) = wavefn(i, j) * exp(ii*chargee*parameters%externalb*(x**2-y**2)/4.0d0/hbar)
             end do
          end do
          
       end select
       
    end if
    
  end subroutine s_initialwavefn
  
  
  
end module setup
