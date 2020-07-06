program cpspropagate

  use precision
  use mpi
  use globals
  use conversion
  use misc
  use setup
  use hamiltonian
  use propagate

  implicit none

  
  complex(cdp), dimension(:, :), allocatable :: wavefn, wavefn_t, wavefn0
  complex(cdp), dimension(:), allocatable :: espectrum
  real(fdp) :: e_before, e_after, e_start, e_finish, plotfrequency, trans, refl, trap, norm
  integer(idp) :: step, tstep
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

  allocate(wavefn(nl(1):nu(1), nl(2):nu(2)), &
       & wavefn0(nl(1):nu(1), nl(2):nu(2)), &
       & wavefn_t(nl(1):nu(1), nl(2):nu(2)), &
       & espectrum(parameters%timesteps))

  call s_restart(tstep, wavefn, wavefn0, nl, nu, e_start, espectrum)

  e_after = f_wavefnenergy(wavefn, nl, nu)

  plotfrequency = dble(parameters%timesteps)/dble(parameters&
       &%plotnumber)

  do step = tstep+1, parameters%timesteps

     e_before = e_after

     call s_display('================================================&
          &================================')
     call s_display('time step= ', char=f_numchar(step)//' / '&
          &//f_numchar(parameters%timesteps))

     !
     !  Propagate the wavefunction the required amount
     !

     call s_propagatewavefn(wavefn, nl, nu, dble(step-1)*parameters&
          &%tstep, dble(step)*parameters%tstep)

     e_after = f_wavefnenergy(wavefn, nl, nu)

     !
     !  Save the wave function to disk in a converted format for
     !  loading into Mathematica
     !
#ifndef MPI
     if (parameters%espectrum) then

        espectrum(step) = sum(wavefn*conjg(wavefn0))&
             &*product(parameters%gridspacing)

     end if
#endif

     if (parameters%plotwfn .and. ceiling(aint(dble(step)&
          &/plotfrequency)*plotfrequency) == step) then

        call s_write2drarray(real(wavefn*conjg(wavefn)), nl, nu,&
             & trim(f_numchar(step)))

     end if

     wavefn_t = wavefn

     call s_c2dfft(wavefn_t, nl, nu, fftfor)

#ifdef MPI
     call s_mpitranspose(wavefn_t, trnl, trnu)
#endif

     call s_cshift(wavefn_t, nl, nu, parameters%n(1)/2, parameters&
          &%n(2)/2)

     if (parameters%plotwfnp .and. ceiling(aint(dble(step)&
          &/plotfrequency)*plotfrequency) == step) then

        call s_write2drarray(real(wavefn_t*conjg(wavefn_t)), nl, nu,&
             & 'p.'//trim(f_numchar(step)), incboundary=.true.)

     end if

     call s_display('norm= ', realnum=f_norm(real(wavefn&
          &*conjg(wavefn)), nl, nu))
     call s_display('system energy= ', realnum=e_after)
     call s_display('change in system energy= ', realnum=e_start&
          &-e_after)

     if (parameters%tranrefl) then
        
        !
        !  Compute and display the amount transmitted and reflected
        !

        potbound(1) = f_gridpointx(parameters%potbound(1))
        potbound(2) = f_gridpointx(parameters%potbound(2))
        potbound(3) = f_gridpointy(parameters%potbound(3))
        potbound(4) = f_gridpointy(parameters%potbound(4))

        !
        ! Compute a normal transmission
        !
 
        refl = f_reflection(real(wavefn*conjg(wavefn)), nl, nu,&
             & potbound)
        trap = f_trapped(real(wavefn*conjg(wavefn)), nl, nu, potbound)
        trans = f_transmission(real(wavefn*conjg(wavefn)), nl, nu,&
             & potbound)
        
        norm = refl + trap + trans
        
        call s_display('====  Real Space  ====         Absolute      &
             &    Relative')
        call s_display('rs. norm= ', realnums=(/ norm, norm/norm /))
        call s_display('rs. reflected= ', realnums=(/ refl, refl/norm&
             & /))
        call s_display('rs. trapped= ', realnums=(/ trap, trap/norm/))
        call s_display('rs. transmitted= ', realnums=(/ trans, trans&
             &/norm /))

        !
        ! Compute transmission in momentum space
        !
        
        refl = f_preflection(real(wavefn_t*conjg(wavefn_t)), nl, nu)
        trans = f_ptransmission(real(wavefn_t*conjg(wavefn_t)), nl,&
             & nu)
        
        norm = refl + trans
        
        call s_display('====  Momentum Space  ====     Absolute      &
             &    Relative')
        call s_display('ms. norm= ', realnums=(/ norm, norm/norm /))
        call s_display('ms. backward= ', realnums=(/ refl, refl/norm&
             &/))
        call s_display('ms. forward= ', realnums=(/ trans, trans/norm&
             & /))
        
     end if

     if ((step/parameters%chkpntfreq)*parameters%chkpntfreq == step) then
        
        call s_checkpoint(step, wavefn, wavefn, nl, nu, e_start, espectrum)

     end if

  end do

  deallocate(wavefn, wavefn0, wavefn_t, espectrum)

  close(potdata)

#ifdef MPI
  call s_mpifinalise()
#endif

end program cpspropagate
