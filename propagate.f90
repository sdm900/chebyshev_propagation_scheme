module propagate

  !
  ! This module contains the routines which actually propagate the
  ! electron through the device.
  !

  use precision
  use globals
  use conversion
  use misc
  use hamiltonian

  implicit none

  interface

     !
     ! Interfaces to the external bessel function routines.  This
     ! improves speed and reliability
     !

     function dbesj0(alpha) result(bes)
       real(8), intent(in) :: alpha
       real(8) :: bes
     end function dbesj0

     function dbesj1(alpha) result(bes)
       real(8), intent(in) :: alpha
       real(8) :: bes
     end function dbesj1

     function dbesjn(i, alpha) result(bes)
       integer(4), intent(in) :: i
       real(8), intent(in) :: alpha
       real(8) :: bes
     end function dbesjn

  end interface

  !
  ! These variables are required from run to run and between modules
  !

  real(fdp), save :: deltaeigen, eigenmin, eigenmax, alpha



contains



  function f_modhamiltonian(wavefn, nl, nu) result(hamwavefn)

    !
    ! Calculate the modified hamiltonian of a wave function as
    ! required by the Chebyshev scheme
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & wavefn
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: hamwavefn


    hamwavefn = 2.0d0/deltaeigen*f_hamiltonian(wavefn, nl, nu) -&
         & (eigenmax + eigenmin)/deltaeigen*wavefn

  end function f_modhamiltonian



!!$  subroutine s_splitpropagation(wavefn, t)
!!$    
!!$    complex(cdp), dimension(:, :), intent(inout) :: wavefn
!!$    real(fdp), intent(in) :: t
!!$    
!!$    complex(cdp), dimension(:, :), allocatable :: splitfn
!!$    complex(cdp), dimension(:, :), allocatable :: wavefn_i
!!$    real(fdp), dimension(dim) :: splitsmooth
!!$    real(fdp), dimension(4) :: split, edge
!!$    integer(idp), dimension(4) :: gridsplit, gridedge
!!$    integer(idp), dimension(dim) :: n
!!$    logical :: cont
!!$    
!!$    
!!$    n = (/ size(wavefn, 1), size(wavefn, 2) /)
!!$    
!!$    !
!!$    ! Setup all the variables needed to use the split operator
!!$    ! propagation
!!$    !
!!$    
!!$    splitsmooth = dble(parameters%smooth) * parameters%gridspacing
!!$    
!!$    split(1) = parameters%potbound(1) - parameters%smooth(1) !!$  
  !!       &*parameters%gridspacing(1) - splitsmooth(1)
!!$    split(2) = parameters%potbound(2) + parameters%smooth(1) !!$  
  !!       &*parameters%gridspacing(1) + splitsmooth(1)
!!$    split(3) = parameters%potbound(3) - parameters%smooth(2) !!$  
  !!       &*parameters%gridspacing(2) - splitsmooth(2)
!!$    split(4) = parameters%potbound(4) + parameters%smooth(2) !!$  
  !!       &*parameters%gridspacing(2) + splitsmooth(2)
!!$    
!!$    gridsplit(1) = max(1, f_gridpointx(split(1))-1)
!!$    gridsplit(2) = min(n(1), f_gridpointx(split(2))+1)
!!$    gridsplit(3) = max(1, f_gridpointy(split(3))-1)
!!$    gridsplit(4) = min(n(2), f_gridpointy(split(4))+1)
!!$    
!!$    edge(1) = split(1) - splitsmooth(1) - abs(t)*parameters
  !!%pmax(1) !!$         &/parameters%mass
!!$    edge(2) = split(2) + splitsmooth(1) + abs(t)*parameters
  !!%pmax(1) !!$         &/parameters%mass
!!$    edge(3) = split(3) - splitsmooth(2) - abs(t)*parameters
  !!%pmax(2) !!$         &/parameters%mass
!!$    edge(4) = split(4) + splitsmooth(2) + abs(t)*parameters
  !!%pmax(2) !!$         &/parameters%mass
!!$    
!!$    gridedge(1) = max(1, f_gridpointx(edge(1))-1)
!!$    gridedge(2) = min(n(1), f_gridpointx(edge(2))+1)
!!$    gridedge(3) = max(1, f_gridpointy(edge(3))-1)
!!$    gridedge(4) = min(n(2), f_gridpointy(edge(4))+1)
!!$    
!!$    cont = .false.
!!$    
!!$    !
!!$    ! Ensure that the new array sizes are actually powers of 2, 3,
  !! 5
!!$    ! or 7
!!$    !
!!$    
!!$    do while ((.not. cont) .and. parameters%primefactorn)
!!$       
!!$       cont = .true.
!!$       
!!$       if (f_factorprime(gridedge(2)-gridedge(1)+1) /=
  !! (gridedge(2) !!$            &-gridedge(1)+1)) then
!!$          
!!$          gridedge(1) = max(gridedge(1)-1, 1)
!!$          gridedge(2) = min(gridedge(2)+1, n(1))
!!$          
!!$          cont = .false.
!!$          
!!$       end if
!!$       
!!$       if (f_factorprime(gridedge(4)-gridedge(3)+1) /=
  !! (gridedge(4) !!$            &-gridedge(3)+1)) then
!!$          
!!$          gridedge(3) = max(gridedge(3)-1, 1)
!!$          gridedge(4) = min(gridedge(4)+1, n(2))
!!$          
!!$          cont = .false.
!!$          
!!$       end if
!!$       
!!$    end do
!!$    
!!$    allocate(splitfn(n(1), n(2)), wavefn_i(gridedge(1):gridedge(2)
  !!, !!$         & gridedge(3):gridedge(4)))
!!$    
!!$    splitfn = 0.0d0
!!$    
!!$    splitfn(gridsplit(1):gridsplit(2), gridsplit(3):gridsplit(4)) 
  !!= !!$         & 1.0d0
!!$    
!!$    call s_blur(splitfn, splitsmooth)
!!$    
!!$    splitfn = real(splitfn)-f_arraymin(real(splitfn))
!!$    splitfn = real(splitfn)/f_arraymax(real(splitfn))
!!$    
!!$    wavefn_i = real(splitfn(gridedge(1):gridedge(2), !!$         &
  !! gridedge(3):gridedge(4))) * wavefn(gridedge(1):gridedge(2) !!$  
  !!       &, gridedge(3):gridedge(4))
!!$    wavefn = (1.0d0 - real(splitfn)) * wavefn
!!$    
!!$    call s_chebpropagation(wavefn_i, t)
!!$    call s_freepropagation(wavefn, t)
!!$    
!!$    wavefn(gridedge(1):gridedge(2), gridedge(3):gridedge(4)) = !!$
  !!         & wavefn(gridedge(1):gridedge(2),
  !! gridedge(3):gridedge(4)) + !!$         & wavefn_i
!!$    
!!$    deallocate(splitfn, wavefn_i)
!!$    
!!$  end subroutine s_splitpropagation



  subroutine s_chebpropagation(wavefn, nl, nu, t)

    !
    ! Propagate the electron via the Chebyshev scheme
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: wavefn
    real(fdp), intent(in) :: t

    complex(cdp), dimension(:, :), allocatable :: chebjm1_t, chebjm0_t
    integer(idp) :: j
    real(fdp) :: tmp, norm, normold, prevnorm
    logical :: lastloop, finished, firstrun, chebwfnopened
    logical, save :: convergent = .true.
    real(fdp), save :: eigenfactor = -1.0d0, highnorm = -1.0d0


    allocate(chebjm1_t(nl(1):nu(1), nl(2):nu(2)),&
         & chebjm0_t(nl(1):nu(1), nl(2):nu(2)))

    inquire(unit=chebwfn, opened=chebwfnopened)

    if (.not. chebwfnopened) then
       open(unit=chebwfn, status='scratch', form='unformatted')
    end if

    rewind(chebwfn)
    write(chebwfn) wavefn

    finished = .false.
    prevnorm = huge(prevnorm)
    firstrun = .true.

    if (highnorm < 0.0d0) then

       highnorm = 1.0d0+sqrt(cutoff)

    end if

    if (eigenfactor < 0.0d0) then

       eigenfactor = parameters%eigenfactor

    end if

    do while (.not. finished)

       rewind(chebwfn)
       read(chebwfn) wavefn

       if (.not. firstrun) then

          eigenfactor = eigenfactor*max(sqrt(2.0), dble(j)/alpha)

       end if

       eigenmin = min(parameters%potrange(1), 0.0d0)*eigenfactor
       eigenmax = sum(parameters%pmax**2)/2.0d0/parameters%mass&
            &*eigenfactor
       deltaeigen = eigenmax-eigenmin

       alpha = eigenmax*t/2.0d0/hbar

       call s_display('eigenmax= ', realnum=eigenmax)
       call s_display('alpha= ', realnum=alpha)
       call s_display('size chebyshev region= ', intnums=nu-nl+1)

       chebjm0_t = wavefn

       chebjm1_t = - f_modhamiltonian(wavefn, nl, nu)

       tmp = dbesj1(alpha)

       wavefn = ii*2.0d0*tmp*chebjm1_t + dbesj0(alpha)*chebjm0_t

       j = 1

       !
       ! The sum of polynomials ends when the bessel functions fall
       ! to zero
       !

       normold = 0.0d0

       norm = f_norm(real(wavefn*conjg(wavefn)), nl, nu)

       if (abs(alpha) < cutoff) then

          lastloop = .true.

       else

          lastloop = .false.

       end if

       do while (.not. lastloop)

          j = j+1

          !
          ! To create faster code and to save memory (one less
          ! temporary variable) the summation is done two terms at a
          ! time.  This allows the reuse of a temporary variable
          !

          chebjm0_t = -2.0d0*f_modhamiltonian(chebjm1_t, nl, nu) -&
               & chebjm0_t

          tmp = dbesjn(j, alpha)

          wavefn = wavefn + (ii**j)*2.0d0*tmp*chebjm0_t

          j = j+1

          chebjm1_t = -2.0d0*f_modhamiltonian(chebjm0_t, nl, nu) -&
               & chebjm1_t

          tmp = dbesjn(j, alpha)

          wavefn = wavefn + (ii**j)*2.0d0*tmp*chebjm1_t

          !
          !  Do some checks to ensure the sum is working
          !

          normold = norm

          norm = f_norm(real(wavefn*conjg(wavefn)), nl, nu)

          if (norm > huge(norm)*(cutoff**2)/dble(product(parameters&
               &%n))) then

             lastloop = .true.

          end if

          if (abs(normold-norm) < cutoff) then

             lastloop = .true.

          end if

       end do

       wavefn = exp(-ii*(eigenmax+eigenmin)*t/2.0d0/hbar)*wavefn

       call s_display('terms in sum= ', intnum=j)

       !
       ! Check whether the sum has converged
       !

       finished = .true.

       norm = f_norm(real(wavefn*conjg(wavefn)), nl, nu)

       if ((norm > highnorm+min(1.0d0, abs(highnorm-1.0d0))&
            &+sqrt(cutoff)) .and. convergent) then

          call s_display('****  It looks like the propagation is non&
               &-convergent')

          call s_display('norm= ', realnum=norm)

          convergent = .false.

          if (abs(1.0d0-prevnorm) > sqrt(2.0d0)*abs(1.0d0-norm)) then

             call s_display('****  Will try again with an increased&
                  & maximum eigen energy')

             prevnorm = norm

             finished = .false.

             convergent = .true.

          end if

       end if

       firstrun = .false.

    end do

!    close(chebwfn)

    deallocate(chebjm1_t, chebjm0_t)

    highnorm = max(norm, highnorm)

  end subroutine s_chebpropagation



  subroutine s_freepropagation(wavefn, nl, nu, t)

    !
    ! Propagate the wave function in free space
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: wavefn
    real(fdp), intent(in) :: t

    complex(cdp), dimension(:, :), allocatable :: psquared
    integer(idp) :: i


    allocate(psquared(nl(1):nu(1), nl(2):nu(2)))

    call s_momentumfac(psquared, nl, nu, (/ 1.0d0, 1.0d0 /), (/ 2, 2 &
         &/))

    !
    ! Compute the laplacian.  Include the normalisation of the
    ! fourier transform in the momentum operator for speed
    !

    psquared = exp(-ii*hbar*(-psquared)/2.0d0/parameters%mass*t) &
         &/dble(parameters%n(1)*parameters%n(2))

    call s_c2dfft(wavefn, nl, nu, fftfor)

    wavefn = wavefn*psquared

    call s_c2dfft(wavefn, nl, nu, fftbac)

    deallocate(psquared)

  end subroutine s_freepropagation



  subroutine s_theobpropagation(wavefn, nl, nu, t)

    !
    ! Propagate the wave function in free space with a uniform
    ! constant magnetic field
    !
    ! This assumes a gaussian wave packet of a certain shape and size
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: wavefn
    real(fdp), intent(in) :: t

    real(fdp) :: invw, x, y, xx, yy, C
    complex(cdp) :: lx, ly
    integer(idp) :: i, j


    invw = 1.0d0/parameters%wavefnwidth(1)
    C = invw/sqrt(pi)*exp(-(parameters%p(1)/2.0d0/hbar/invw)**2 &
               &-(parameters%p(2)/2.0d0/hbar/invw)**2)

    lx = ii*parameters%p(1)/2.0d0/hbar/invw*exp(-ii*parameters%omegal&
         &*t)

    ly = ii*parameters%p(2)/2.0d0/hbar/invw*exp(-ii*parameters%omegal&
         &*t)

    do j = nl(2), nu(2)

       y = f_positiony(j) - parameters%position(2)

       do i = nl(1), nu(1)

          x = f_positionx(i) - parameters%position(1)

          xx = x*cos(parameters%omegal*t) - y*sin(parameters%omegal*t)
          yy = y*cos(parameters%omegal*t) + x*sin(parameters%omegal*t)

          wavefn(i, j) = C*exp(-lx**2-ly**2+2.0d0*invw*lx*xx+2.0d0*invw*ly*yy&
               &-invw**2*xx**2/2.0d0 - invw**2*yy**2/2.0d0-ii*parameters&
               &%omegal*t)
          
       end do

    end do

  end subroutine s_theobpropagation



  subroutine s_theopropagation(wavefn, nl, nu, t)

    !
    ! Propagate the wave function in free space by they theoretical
    ! equation
    !
    ! This assumes a gaussian wave packet
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: wavefn
    real(fdp), intent(in) :: t

    real(fdp) :: theta1, theta2, phi1, phi2
    complex(cdp) :: eiphi1, eiphi2
    real(fdp) :: x, y
    integer(idp) :: i, j


    theta1 = atan(hbar*t/parameters%mass/parameters%wavefnwidth(1) &
         &**2)/2.0d0
    theta2 = atan(hbar*t/parameters%mass/parameters%wavefnwidth(2) &
         &**2)/2.0d0

    phi1 = -theta1-parameters%p(1)**2*t/2.0d0/parameters%mass/hbar
    phi2 = -theta2-parameters%p(2)**2*t/2.0d0/parameters%mass/hbar

    eiphi1 = exp(ii*phi1)/(4.0d0*parameters%wavefnwidth(1)**4+4.0d0 &
         &*hbar**2*t**2/parameters%mass**2)**(1.0d0/4.0d0)
    eiphi2 = exp(ii*phi2)/(4.0d0*parameters%wavefnwidth(2)**4+4.0d0 &
         &*hbar**2*t**2/parameters%mass**2)**(1.0d0/4.0d0)

    do j=nl(2), nu(2)

       y = f_positiony(j)

       do i=nl(1), nu(1)

          x = f_positionx(i)

          wavefn(i, j) = (16.0d0*parameters%wavefnwidth(1)**2 &
               &*parameters%wavefnwidth(2)**2/pi**2)**(1.0d0/4.0d0) &
               &*eiphi1*eiphi2* exp(ii*parameters%p(1)*x/hbar+ii &
               &*parameters%p(2)*y/hbar- (x-parameters%position(1) &
               &-parameters%p(1)/parameters%mass*t)**2/(2.0d0 &
               &*parameters%wavefnwidth(1)**2+2.0d0*ii*hbar*t &
               &/parameters%mass)- (y-parameters%position(2) &
               &-parameters%p(2)/parameters%mass*t)**2/(2.0d0 &
               &*parameters%wavefnwidth(2)**2+2.0d0*ii*hbar*t &
               &/parameters%mass))

       end do

    end do


  end subroutine s_theopropagation



  subroutine s_propagatewavefn(wavefn, nl, nu, t1, t2)

    !
    ! Setup the various parts of the wave function for propagation
    ! under free space and the Chebyshev scheme
    !
    ! There are three possibilities:
    !
    !     entirely chebyshev
    !     split operator
    !     entirely free space
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: wavefn
    real(fdp), intent(in) :: t1, t2


    if (defined_parameters%potential .or. defined_parameters &
         &%externalb .or. defined_parameters%optimisation) then

       select case (parameters%optimisation)

!!$       case ('splitwavefn')
!!$
!!$           call s_splitpropagation(wavefn, t)
!!$
!!$          call s_display('****  WARNING: The split wave function
          !! !!$               & propagation is not working at the
          !! moment.')          


       case ('exponential', 'gaussian', 'linear', 'none')
          call s_chebpropagation(wavefn, nl, nu, t2-t1)

       case ('freespace')
          call s_freepropagation(wavefn, nl, nu, t2-t1)

       case ('theoreticalfreespace')
          call s_theopropagation(wavefn, nl, nu, t2)

       case ('theoreticalbfield')
          call s_theobpropagation(wavefn, nl, nu, t2)

       case default
          call s_display('****  WARNING: no acceptable optimisation&
               & defined, wavefn not propagated')          

       end select

    else

       !
       ! No potential, so just propagate in free space
       !

       call s_freepropagation(wavefn, nl, nu, t2-t1)

    end if

  end subroutine s_propagatewavefn



end module propagate
