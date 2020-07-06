module globals

  !
  ! This code uses atomic units in all calculations.  These are not
  ! always the way you want to quote them, so the conversions are:
  !
  ! Time:        1 a.u. = 2.41888E-17 s
  ! Energy:      1 a.u. = 27.2114 eV
  ! Distance:    1 a.u. = 5.29177E-11 m
  ! Momentum:    1 a.u. = 1.99285E-24 kg m/s
  ! Magnetism:   1 a.u. = 2.35050E+5 T
  !



  use precision

  public

  !
  ! Define some global parameters
  !

  integer(idp), parameter :: dim = 2, potdata = 99, potlog = 98,&
       & chebwfn = 97, inputfile=96, chkptfile=95, maxpotcomps = 50,&
       & maxstringlength = 256, fftfor = -1, fftbac = 1,&
       & maxfilenamelen = 256, numcharlen = 6

  real(fdp), parameter :: pi =&
       & 3.141592653589793238462643383279502884197, cutoff = 1.0d-14,&
       & boltzk = 1.380658d-23, h = 6.6260755d-34, hbar = h/2.0d0/pi,&
       & masse = 9.1093898d-31, nm = 1.0d-9, chargee = 1.60217733d-19&
       & , eV = chargee, aut = 2.41888d-17, aue = 27.2114*eV, aus =&
       & 5.29177d-11, aup = 1.99285d-24, aub = 2.35050d5, c =&
       & 2.99792458d8, mu0 = 4.0d0*pi*1.0d-7

  complex(cdp), parameter :: ii = (0.0d0, 1.0d0)

  real(fdp), dimension(15), parameter :: constantvalues = (/ pi, c,&
       & masse , ev, nm, chargee, h, hbar, boltzk, mu0, aut, aue, aus&
       &, aup, aub /)

  character(LEN = *), parameter :: constants =&
       & 'pi c masse ev nm chargee h hbar boltzk mu0 aut aue aus aup&
       & aub'


  !
  ! Define the type and default values for all the user-defineable
  ! global variables
  !

  type params
     
     character(LEN = 50) :: datadir = 'data/', chkpntfile = 'chkpt.',&
          & basename = 'wavefn.', optimisation = 'none' ,&
          & potsmoothtype = 'none', initialwavefn = 'gaussian', gauge&
          & = '(-y,x,0)'
     
     real(fdp) :: mass = 6.67d-2*masse, fermienergy = 0.0d0,&
          & temperature = 0.01d0, initialenergy = 0.0d0, t = 0.0d0,&
          & distance = 1.0d4*aus, tstep = 0.0d0, speccutoff = 0.01d0,&
          & externale = 0.0d0, externalb = 0.0d0, omegal = 0.0d0,&
          & eigenfactor = 1.25d0

     real(fdp), dimension(dim) :: p = 0.0d0, pspread = 0.05d0, pmax =&
          & 0.0d0, accuracyfactor = 1.0d0, position = 0.0d0,&
          & wavefnwidth = 0.0d0, gridspacing = 0.0d0, pspacing =&
          & 0.0d0, spacesize = 2.0d4*aus, maxgridspacing = 0.0d0

     real(fdp), dimension(2) :: potrange = 0.0d0

     integer(idp), dimension(dim) :: n = 0, extraboundary = 2, smooth&
          & = 32, plotpointnum = 0

     integer(idp) :: timesteps = 1, plotnumber = 0, chkpntfreq = 1

     logical, dimension(maxpotcomps) :: potcomps = .false.

     logical, dimension(2*dim) :: comppotpos = .true.

     real(fdp), dimension(4) :: potbound = 0.0d0

     logical :: plotwfn = .true., plotwfnp = .true., plotpot = .true.&
          & , tranrefl = .true., plotcomplex = .true., espectrum =&
          & .true., pspectrum = .true., tspectrum = .true.,&
          & conductance = .true. , primefactorn = .true.

  end type params


  type def_params

     !
     ! Defines parameters which are avalaible for all routines
     !

     logical :: datadir = .false., chkpntfile = .false., basename = .false., optimisation = .false.,&
          & potsmoothtype = .false., initialwavefn = .false., mass =&
          & .false., fermienergy = .false., temperature = .false.,&
          & initialenergy = .false., t = .false., distance = .false.,&
          & tstep = .false., speccutoff = .false., p = .false.,&
          & pspread = .false., pmax = .false., accuracyfactor =&
          & .false., position = .false., wavefnwidth = .false.,&
          & gridspacing = .false., pspacing = .false., spacesize =&
          & .false., potrange = .false., n = .false., extraboundary =&
          & .false., plotpointnum = .false., timesteps = .false., chkpntfreq = .false.,&
          & smooth = .false., plotnumber = .false., potcomps =&
          & .false., comppotpos = .false., potbound = .false.,&
          & plotwfn = .false., plotwfnp = .false., plotpot = .false.,&
          & tranrefl = .false., plotcomplex = .false., espectrum =&
          & .false., pspectrum = .false., tspectrum = .false.,&
          & conductance = .false., externale = .false., externalb =&
          & .false., potential = .false., primefactorn = .false.,&
          & omegal = .false., gauge = .false., eigenfactor = .false.,&
          & maxgridspacing = .false.

  end type def_params


  type(params) :: parameters
  type(def_params) :: defined_parameters

end module globals
