program cpssetup

  use precision
  use mpi
  use globals
  use conversion
  use misc
  use setup
  use hamiltonian
  use propagate

  implicit none

  
  complex(cdp), dimension(:, :), allocatable :: wavefn, wavefn_t
  complex(cdp), dimension(:), allocatable :: espectrum
  real(fdp) :: e_start
  integer(idp), dimension(dim) :: nl, nu, trnl, trnu


#ifdef MPI  
  call s_mpiinit()
#endif

  open(unit=potdata, status='scratch', form="formatted")

  !
  !  Read the input file and setup initialise global variables
  !

  call s_setup()

#ifdef MPI
  call s_mpiarraysize(nl, nu, trnl, trnu)
#else
  nl = (/ 1, 1 /)
  nu = (/ parameters%n(1), parameters%n(2) /)
  trnl = (/ 1, 1 /)
  trnu = (/ parameters%n(2), parameters%n(1) /)
#endif

  allocate(wavefn(nl(1):nu(1), nl(2):nu(2)), wavefn_t(nl(1):nu(1), nl(2):nu(2)), espectrum(parameters%timesteps))

  !
  !  Initialise the wavefunction
  !

  call s_initialwavefn(wavefn, nl, nu)

  call s_display('================================================================================')
  call s_display('time step= ', char=f_numchar(0)//' / '//f_numchar(parameters%timesteps))

  !
  !  Propagate the wavefunction for t=0s.  This is done to initialise
  !  the potential and various components required for the full
  !  propagation
  !

  call s_propagatewavefn(wavefn, nl, nu, 0.0d0, 0.0d0)

  !
  !  Save the initial wave function to disk in a converted format
  !  for loading into Mathematica
  !

  call s_display('norm= ', realnum=f_norm(real(wavefn*conjg(wavefn)), nl, nu))
  call s_display('system energy= ', realnum=f_wavefnenergy(wavefn, nl, nu))

  if (parameters%plotwfn) then

     call s_write2drarray(real(wavefn*conjg(wavefn)), nl, nu, trim(f_numchar(0)))
#if 1
     call s_write2drarray(real(wavefn), nl, nu, 'real.'//trim(f_numchar(0)))
     call s_write2drarray(aimag(wavefn), nl, nu, 'comp.'//trim(f_numchar(0)))
#endif

  end if

  wavefn_t = wavefn

  call s_c2dfft(wavefn_t, nl, nu, fftfor)

#ifdef MPI
  call s_mpitranspose(wavefn_t, trnl, trnu)
#endif

  call s_cshift(wavefn_t, nl, nu, parameters%n(1)/2, parameters%n(2)/2)

  if (parameters%plotwfnp) then

     call s_write2drarray(real(wavefn_t*conjg(wavefn_t)), nl, nu,&
          & 'p.'//trim(f_numchar(0)), incboundary=.true.)

  end if

  e_start = f_wavefnenergy(wavefn, nl, nu)

  espectrum = 0.0d0

  call s_checkpoint(0, wavefn, wavefn, nl, nu, e_start, espectrum)

  deallocate(wavefn, wavefn_t, espectrum)

  close(potdata)

#ifdef MPI
  call s_mpifinalise()
#endif

end program cpssetup
