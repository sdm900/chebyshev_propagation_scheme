module hamiltonian

  use precision
  use mpi
  use globals
  use fft
  use misc
  use potential

  implicit none

  

contains



  function f_hamiltonian(wavefn, nl, nu) result(hamiltonian)

    !
    ! This function computes the hamiltonian of a wave function
    !
    
    integer(idp), dimension(2), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & wavefn
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: hamiltonian

    complex(cdp), dimension(:, :), allocatable, save :: del2, pot,&
         & delwavefnx, delwavefny, delwavefn, delx, dely, delxy,&
         & vecpotfn2
    real(fdp), dimension(:, :), allocatable, save :: vecpotfnx,&
         & vecpotfny, vecpotfn
    integer(idp), dimension(dim), save :: nlold = 0, nuold = 0
    logical, save :: defvecpot = .false., defpot = .false., defdel =&
         & .false.


    !
    ! If the space size changes, deallocate all the stored function
    !

    if (nlold(1) /= nl(1) .or. nlold(2) /= nl(2) .or. nuold(1) /=&
         & nu(1) .or. nuold(2) /= nu(2)) then
       
       if (defdel) then
          
          deallocate(del2)

          defdel = .false.
          
       end if
       
       if (defvecpot) then

          select case (parameters%gauge)

          case('(-y,0,0)')
             deallocate(delx, vecpotfny, vecpotfn2, delwavefnx)

          case('(-y,x,0)')
             
             deallocate(delx, dely, vecpotfnx,  vecpotfny, delwavefnx&
                  &, delwavefny, vecpotfn2)
             
          case('(x-y,x-y,0)', 'fft')
             deallocate(delxy, vecpotfn,  delwavefn, vecpotfn2)
             
          end select
          
          defvecpot = .false.
          
       end if
       
       if (defpot) then
          
          deallocate(pot)
          
          defpot = .false.
          
       end if
       
       nlold = nl
       nuold = nu

    end if

    !
    ! If the required function are not allocated and computed,
    ! allocate them and compute them
    !

    if (.not. defdel) then

       allocate(del2(nl(1):nu(1), nl(2):nu(2)))
       
       call s_momentumfac(del2, nl, nu, (/ 1.0d0, 1.0d0 /)&
            &/dble(product(parameters%n)), (/ 2, 2 /))

       !
       ! Include terms for speed
       !

       del2 = - del2 * hbar**2 / 2.0d0 / parameters%mass

       defdel = .true.
       
    end if

    if (defined_parameters%potential .and. (.not. defpot) .and.&
         & (.not. defvecpot)) then
       
       !
       ! Define the potential
       !

       allocate(pot(nl(1):nu(1), nl(2):nu(2)))

       call s_setpotential(pot, nl, nu)

       !
       ! Save the potential to disk for plotting
       !
       
       if (parameters%plotpot) then
          
          call s_write2drarray(real(pot), nl, nu, filename='rpot',&
               & incboundary=.true.)
          call s_write2drarray(aimag(pot), nl, nu, filename='ipot',&
               & incboundary=.true.)
          
       end if

       defpot = .true.
       
    end if

    if (defined_parameters%externalb .and. (.not. defvecpot) ) then

       select case (parameters%gauge)
          
       case('(-y,0,0)')
          
          allocate(delx(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfny(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfn2(nl(1):nu(1), nl(2):nu(2)),&
               & delwavefnx(nl(1):nu(1), nl(2):nu(2)))
          
          call s_bpoty(delx, vecpotfny, vecpotfn2, nl, nu)
          
       case('(-y,x,0)')
          
          allocate(delx(nl(1):nu(1), nl(2):nu(2)), dely(nl(1):nu(1),&
               & nl(2):nu(2)), vecpotfnx(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfny(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfn2(nl(1):nu(1), nl(2):nu(2)),&
               & delwavefnx(nl(1):nu(1), nl(2):nu(2)),&
               & delwavefny(nl(1):nu(1), nl(2):nu(2)))

          call s_bpotyx(delx, dely, vecpotfnx, vecpotfny, vecpotfn2,&
               & nl, nu)
          
       case('(x-y,x-y,0)')
          allocate(delxy(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfn(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfn2(nl(1):nu(1), nl(2):nu(2)),&
               & delwavefn(nl(1):nu(1), nl(2):nu(2)))

          call s_bpotxyxy(delxy, vecpotfn, vecpotfn2, nl, nu)
                    
       case('fft')
          allocate(delxy(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfn(nl(1):nu(1), nl(2):nu(2)),&
               & vecpotfn2(nl(1):nu(1), nl(2):nu(2)),&
               & delwavefn(nl(1):nu(1), nl(2):nu(2)))

          call s_bpotfft(delxy, vecpotfn, vecpotfn2, nl, nu)
                    
       end select
    
       defvecpot = .true.

       if (defpot) then

          vecpotfn2 = vecpotfn2 + pot
          
          deallocate(pot)
          
          defpot = .false.

       end if
       
    end if
    
    !
    ! Compute the derivatives of the wavefn
    !

    hamiltonian = wavefn
    
    call s_c2dfft(hamiltonian, nl, nu, fftfor)

    if (defined_parameters%externalb) then
       
       select case (parameters%gauge)
          
       case('(-y,0,0)')
          delwavefnx = hamiltonian*delx
          
          call s_c2dfft(delwavefnx, nl, nu, fftbac)
          
       case('(-y,x,0)')
          delwavefnx = hamiltonian*delx
          delwavefny = hamiltonian*dely
          
          call s_c2dfft(delwavefnx, nl, nu, fftbac)
          call s_c2dfft(delwavefny, nl, nu, fftbac)
          
       case('(x-y,x-y,0)', 'fft')
          delwavefn = hamiltonian*delxy

          call s_c2dfft(delwavefn, nl, nu, fftbac)
          
       end select

    end if

    hamiltonian = hamiltonian*del2
    
    call s_c2dfft(hamiltonian, nl, nu, fftbac)

    !
    ! Build up the hamiltonian from its various bits
    !

    if (defined_parameters%externalb .and. defvecpot) then
       
       select case (parameters%gauge)
          
       case('(-y,0,0)')
          hamiltonian = hamiltonian + vecpotfny*delwavefnx +&
               & vecpotfn2*wavefn

       case('(-y,x,0)')
          hamiltonian = hamiltonian + vecpotfnx*delwavefny +&
               & vecpotfny*delwavefnx + vecpotfn2*wavefn

       case('(x-y,x-y,0)', 'fft')
          hamiltonian = hamiltonian + vecpotfn*delwavefn + vecpotfn2&
               &*wavefn

       end select
       
    end if
    
    if (defined_parameters%potential .and. defpot) then
       
       hamiltonian = hamiltonian + pot*wavefn

    end if

    !
    ! If we have a magnetic field, ensure the boundary of the wave
    ! function
    ! is zero.  This is a nasty hack, but I think required. 
    ! Hopefully it
    ! will make this work well for magnetic fields.
    !

    if(defined_parameters%externalb) then

       call s_smoothedge(hamiltonian, nl, nu, parameters%smooth)

    end if

  end function f_hamiltonian



  function f_wavefnenergy(wavefn, nl, nu) result(energy)

    !
    ! this module calculates the energy of a function
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & wavefn
    real(fdp) :: energy

    !
    !  Compute the hamiltonian
    !

    energy = f_norm(real(conjg(wavefn)*f_hamiltonian(wavefn, nl, nu))&
         &, nl, nu)

  end function f_wavefnenergy



end module hamiltonian
