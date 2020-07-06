program cpsfinish

  use precision
  use mpi
  use globals
  use conversion
  use misc
  use setup
  use hamiltonian
  use propagate

  implicit none

  
  complex(cdp), dimension(:, :), allocatable :: wavefn, wavefn0
  complex(cdp), dimension(:), allocatable :: espectrum
  real(fdp) :: e_before, e_after, e_start, e_finish, plotfrequency, trans, refl, trap, norm
  integer(idp) :: step
  integer(idp), dimension(4) :: potbound
  integer(idp), dimension(dim) :: nl, nu, trnl, trnu

#ifdef MPI  
  call s_mpiinit()
#endif

  open(unit=potdata, status='scratch', form="formatted")

  call s_setup()

#ifdef MPI
  call s_mpiarraysize(nl, nu, trnl, trnu)
#else
  nl = (/ 1, 1 /)
  nu = (/ parameters%n(1), parameters%n(2) /)
  trnl = (/ 1, 1 /)
  trnu = (/ parameters%n(2), parameters%n(1) /)
#endif

  call s_checkrestart(nl, nu)

  allocate(wavefn(nl(1):nu(1), nl(2):nu(2)), wavefn0(nl(1):nu(1), nl(2):nu(2)), espectrum(parameters%timesteps))

  call s_restart(step, wavefn, wavefn0, nl, nu, e_start, espectrum)

  e_finish = f_wavefnenergy(wavefn, nl, nu)

  call s_display('===================================================&
       &=============================')
  call s_display('final norm= ', realnum=f_norm(real(wavefn&
       &*conjg(wavefn)), nl, nu))
  call s_display('final energy= ', realnum=e_finish)
  call s_display('tot. change in energy= ', realnum=e_start -&
       & e_finish)
  call s_display('===================================================&
       &=============================')

  !
  !  Do various post processing operations
  !

#ifndef MPI
  if (parameters%pspectrum) then

     call s_pspectrum(wavefn, wavefn0, nl, nu)

  end if

  if (parameters%espectrum) then

     call s_espectrum(espectrum)

  end if

  if (parameters%tspectrum) then

     call s_tspectrum(wavefn, wavefn0, nl, nu)

  end if

!!$  if (parameters%conductance) then
!!$
!!$     call s_conductance(wavefn, nl, nu)
!!$
!!$  end if
#endif

  deallocate(wavefn, wavefn0, espectrum)

  close(potdata)

#ifdef MPI
  call s_mpifinalise()
#endif

end program cpsfinish
