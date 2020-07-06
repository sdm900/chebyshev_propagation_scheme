module potential

  use precision
  use globals
  use mpi
  use interpreter
  use misc

  implicit none



contains



  subroutine s_setpotential(pot, nl, nu)

    !
    ! Define the electronic potential for use in propagation.  It
    ! reads in a scratch file which contains all the layout
    ! information for the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(out) ::&
         & pot

    logical, dimension(:, :), allocatable :: mask
    character(len=20) :: object
    integer(idp) :: t_int1, i, j, nxmin, nxmax, nymin, nymax, nmin,&
         & smooth
    logical :: donecomponent
    complex(cdp) :: cmplxboundary, boundaryheight


    allocate(mask(nl(1):nu(1), nl(2):nu(2)))

    mask = .false.
    pot = 0.0d0

    do i=1, maxpotcomps

       if (parameters%potcomps(i)) then

          rewind(potdata)

          donecomponent = .false.

          do while(.not. donecomponent)

             read(potdata, '(a)', advance='no') object

             select case (trim(object))

             case('potfn')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_potfn(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('circle')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_circle(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('fillcircle')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_fillcircle(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('ellipse')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_ellipse(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('fillellipse')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_fillellipse(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('line')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_line(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('ramp')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_ramp(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('polygon')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_polygon(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case ('fillpolygon')
                read(potdata, '(i)', advance='no') t_int1

                if (t_int1 == i) then

                   call s_fillpolygon(pot, mask, nl, nu)

                   donecomponent = .true.

                end if

             case default

             end select

             read(potdata, '(a)')

          end do

       end if

    end do

    if (defined_parameters%potsmoothtype .and. f_arraymax(abs(pot),&
         & nl, nu) > 0.0d0) then

       select case (parameters%potsmoothtype)

!!$       case ('poisson')
!!$
!!$          !
!!$          ! Apply Poisson's equation
!!$          !
!!$
!!$          call s_display('applying poisson equation')
!!$
!!$          pot = real(pot)
!!$
!!$          call s_relaxpoisson(pot, mask, nl, nu, parameters
          !!%gridspacing)
!!$
!!$          pot = real(pot)

       case ('gaussblur')

          !
          ! Blur the potential
          !

          call s_display('applying gauss blur, smooth= ', intnums&
               &=parameters%smooth)

          pot = real(pot)

          call s_blur(pot, nl, nu)

          pot = real(pot)

       case('smooth')

          !
          ! Apply an averaging type smooth
          !

          call s_display('applying smooth, smooth= ', intnum&
               &=maxval(parameters%smooth))

          pot = real(pot)

          call s_smooth(pot, nl, nu, maxval(parameters%smooth))

          pot = real(pot)

       case default

       end select

    end if

    if (defined_parameters%optimisation) then

       cmplxboundary = cmplx(0.0d0, -(sum(parameters%pmax**2)/2.0d0 &
            &/parameters%mass)**(5.0d0/4.0d0))

       do j = nl(2), nu(2)

          do i = nl(1), nu(1)

             nxmin = abs(i-1)
             nxmax = abs(parameters%n(1)-i)
             nymin = abs(j-1)
             nymax = abs(parameters%n(2)-j)
             nmin = min(nxmin, nxmax, nymin, nymax)

             if (nmin == nxmin .or. nmin == nxmax) then

                smooth = parameters%smooth(1)

             else

                smooth = parameters%smooth(2)

             end if

             if (nmin < smooth) then

                if ((nmin == nxmin .and. parameters%comppotpos(1))&
                     & .or. (nmin == nxmax .and. parameters&
                     &%comppotpos(2)) .or. (nmin == nymin .and.&
                     & parameters%comppotpos(3)) .or. (nmin == nymax&
                     & .and. parameters%comppotpos(4))) then

                   select case (parameters%optimisation)
                      
                   case ('exponential')
                      
                      !      
                      ! The optimal potential
                      !
                      
                      boundaryheight = cmplxboundary*exp(&
                           &-dble(smooth)/dble(smooth - nmin))&
                           &*exp(1.0d0)
                      
                   case ('gaussian')
                      
                      !
                      ! Gaussian potential
                      !
                      
                      boundaryheight = -cmplxboundary*(exp(-(2.0d0&
                           &*dble(smooth - nmin)/dble(smooth))**2)&
                           &-1.0d0)
                      
                   case ('linear')
                      
                      !
                      ! The linear ramp
                      !
                      
                      boundaryheight = cmplxboundary*dble(smooth -&
                           & nmin)/dble(smooth)
                      
                   case default
                      exit
                      
                   end select
                   
                   pot(i, j) = boundaryheight
                   
                end if
                
             end if

          end do

       end do
       
    end if
    
    if (defined_parameters%externale) then
       
       do i = nl(1), nu(1)
          
          if (i <= f_gridpointx(parameters%potbound(1))) then
             
             pot(i, :) = pot(i, :) + parameters%externale
             
          elseif (i <= f_gridpointx(parameters%potbound(2))) then
             
             pot(i, :) = pot(i, :) + (1.0d0 - dble(i -&
                  & f_gridpointx(parameters%potbound(1)))&
                  &/dble(f_gridpointx(parameters%potbound(2)) -&
                  & f_gridpointx(parameters%potbound(1))))*parameters&
                  &%externale

          end if

       end do

    end if

    deallocate(mask)

  end subroutine s_setpotential



  subroutine s_fillmask(mask, nl, nu)

    !
    ! Fill in an object, given a mask
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    integer(idp), dimension(:, :), allocatable :: potxl, potxr, potyu&
         &, potyd
    integer(idp) :: i, j
#ifdef MPI
    integer(idp), dimension(:), allocatable :: buffer
    logical :: yusent, ydsent
    integer, dimension(4) :: mpirequests
    integer, dimension(MPI_STATUS_SIZE, 4) :: mpistatuss
    integer :: mpiindex
#endif


    allocate(potxl(nl(1):nu(1), nl(2):nu(2)), potxr(nl(1):nu(1),&
         & nl(2):nu(2)), potyu(nl(1):nu(1), nl(2):nu(2)),&
         & potyd(nl(1):nu(1), nl(2):nu(2)))

    potxl = 0
    potxr = 0
    potyu = 0
    potyd = 0

    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          if (i > nl(1)) then
             
             if (mask(i, j) /= mask(i-1, j)) then
                
                potxl(i, j) = modulo(potxl(i-1, j)+1, 4)
                
             else
                
                potxl(i, j) = potxl(i-1, j)
                
             end if
             
          elseif (mask(i, j)) then
             
             potxl(i, j) = 1
             
          end if

          if (nu(1)-i+nl(1) < nu(1)) then
             
             if (mask(nu(1)-i+nl(1), j) /= mask(nu(1)-i+nl(1)+1, j))&
                  & then
                
                potxr(nu(1)-i+nl(1), j) = modulo(potxr(nu(1)-i+nl(1)&
                     &+1, j)+1, 4)
                
             else
                
                potxr(nu(1)-i+nl(1), j) = potxr(nu(1)-i+nl(1)+1, j)
                
             end if

          elseif (mask(nu(1)-i+nl(1), j)) then
             
             potxr(nu(1)-i+nl(1), j) = 1
             
          end if

          if (j > nl(2)) then
             
             if(mask(i, j) /= mask(i, j-1)) then
                
                potyu(i, j) = modulo(potyu(i, j-1)+1, 4)
                
             else
                
                potyu(i, j) = potyu(i, j-1)
                
             end if

          elseif (mask(i, j)) then
                
             potyu(i, j) = 1

          end if
          
          if (nu(2)-j+nl(2) < nu(2)) then
             
             if (mask(i, nu(2)-j+nl(2)) /= mask(i, nu(2)-j+nl(2)+1))&
                  & then
                
                potyd(i, nu(2)-j+nl(2)) = modulo(potyd(i, nu(2)-j&
                     &+nl(2)+1)+1, 4)
                
             else
                
                potyd(i, nu(2)-j+nl(2)) = potyd(i, nu(2)-j+nl(2)+1)
                
             end if

          elseif (mask(i, nu(2)-j+nl(2))) then
             
             potyd(i, nu(2)-j+nl(2)) = 1
             
          end if

        end do
       
    end do

#ifdef MPI
    mpirequests = MPI_REQUEST_NULL
    mpistatuss = 0
    yusent = .false.
    ydsent = .false.

    allocate(buffer(nl(1):nu(1)))

    if (mpirank == mpistart) then

       call s_mpisendint(potyu(nl(1), nu(2)), size(potyu(:, nu(2))),&
            & mpirank+1, mpirequests(1))

       yusent = .true.

    elseif (mpirank == mpiend) then
       
       call s_mpisendint(potyd(nl(1), nl(2)), size(potyd(:, nl(2))),&
            & mpirank-1, mpirequests(1))
       
       ydsent = .true.

    end if

    do while (.not. (yusent .and. ydsent)) 

       call s_mpirecvint(buffer(1), size(buffer), MPI_ANY_SOURCE,&
            & mpirequests(2))

       call s_mpiwaitall(size(mpirequests), mpirequests(1),&
            & mpistatuss(1, 1))

       if (mpistatuss(MPI_SOURCE, 2) == mpirank-1) then

          do j = nl(2), nu(2)

             do i = nl(1), nu(1)

                potyu(i, j) = modulo(potyu(i, j)+buffer(i), 4)

             end do

          end do

          if (mpirank < mpiend) then
             
             call s_mpisendint(potyu(nl(1), nu(2)), size(potyu(:,&
                  & nu(2))), mpirank+1, mpirequests(3))
             
          end if

          yusent = .true.
          
       elseif (mpistatuss(MPI_SOURCE, 2) == mpirank+1) then

          do j = nl(2), nu(2)
             
             do i = nl(1), nu(1)
                
                potyd(i, j) = modulo(potyd(i, j)+buffer(i), 4)
                
             end do

          end do

          if (mpirank > mpistart) then
             
             call s_mpisendint(potyd(nl(1), nl(2)), size(potyd(:,&
                  & nl(2))), mpirank-1, mpirequests(4))

          end if
          
          ydsent = .true.
          
       end if

    end do

    call s_mpiwaitall(size(mpirequests), mpirequests(1), mpistatuss(1&
         &, 1))

    deallocate(buffer)
#endif

    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          if (mask(i, j) .or. (potxl(i, j) == potxr(i, j) .and.&
               & potxl(i, j) > 0 .and. (potyu(i, j) > 0 .or. potyd(i,&
               & j) > 0)) .or. (potyu(i, j) == potyd(i, j) .and.&
               & potyu(i, j) > 0 .and. (potxl(i, j) > 0 .or. potxr(i,&
               & j) > 0))) then

             mask(i, j) = .true.

          else

             mask(i, j) = .false.

          end if

       end do

    end do

    deallocate(potxl, potxr, potyu, potyd)

 end subroutine s_fillmask



  subroutine s_circlemask(mask, nl, nu, centre, radius, width)

    !
    ! Define a circle mask
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask
    real(fdp), dimension(dim), intent(in) :: centre
    real(fdp), intent(in) :: radius, width

    real(fdp), dimension(dim) :: position
    real(fdp) :: d
    integer(idp) :: i, j


    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          position = f_position( (/ i, j /) )

          d = abs(sqrt(sum((centre-position)**2))-radius)

          if (d <= width) then

             mask(i, j) = .true.

          end if

       end do

    end do

  end subroutine s_circlemask



  subroutine s_ellipsemask(mask, nl, nu, centre, xylen, width)

    !
    ! Define an ellipse mask
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask
    real(fdp), dimension(dim), intent(in) :: centre, xylen
    real(fdp), intent(in) :: width

    real(fdp), dimension(dim) :: position
    real(fdp) :: d
    integer(idp) :: i, j


    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          position = f_position( (/ i, j /) )

          d =  abs(sqrt(sum((centre-position)**2))*(1.0d0-xylen(1)&
               &*xylen(2)/sqrt(sum(xylen**2*(centre-position)**2))))

          if (d <= width) then

             mask(i, j) = .true.

          end if

       end do

    end do

  end subroutine s_ellipsemask



  subroutine s_linemask(mask, nl, nu, point1, point2, width)

    !
    ! Define a line mask
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask
    real(fdp), dimension(dim), intent(in) :: point1, point2
    real(fdp), intent(in) :: width

    real(fdp), dimension(dim) :: position
    real(fdp) :: a, b, c, d
    integer(idp) :: i, j


    c = sqrt(sum((point1-point2)**2))

    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          position = f_position( (/ i, j /) )

          a = sqrt(sum((point1-position)**2))

          b = sqrt(sum((point2-position)**2))

          d = sqrt(abs(-a**4+2*(b**2+c**2)*a**2-(b**2-c**2)**2))&
               &/2.0d0/c

          if (abs(a-d)/nm < sqrt(cutoff)) then

             d = a

          end if

          if (abs(b-d)/nm < sqrt(cutoff)) then

             d = b

          end if

          if ( abs(sqrt(a**2-d**2)+sqrt(b**2-d**2)-c)/nm <=&
               & sqrt(cutoff) .and. d <= width .or. a <= width .or. b&
               & <= width) then

             mask(i, j) = .true.

          end if

       end do

    end do

  end subroutine s_linemask



  subroutine s_circle(pot, mask, nl, nu)

    !
    ! Define a circle component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp), dimension(dim) :: centre, position
    real(fdp) :: radius, width, height


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(2f, f, f, f)', advance='no') centre, radius,&
         & width, height

    call s_circlemask(t_mask, nl, nu, centre, radius, width)

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_circle



  subroutine s_fillcircle(pot, mask, nl, nu)

    !
    ! Define a filled circle component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp), dimension(dim) :: centre
    real(fdp) :: radius, height


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(f, f, 2f)', advance='no') centre, radius, height

    call s_circlemask(t_mask, nl, nu, centre, radius,&
         & sqrt(sum(parameters%gridspacing**2)))

    call s_fillmask(t_mask, nl, nu)

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_fillcircle



  subroutine s_ellipse(pot, mask, nl, nu)

    !
    ! Define an ellipse component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp), dimension(dim) :: centre, xylen, position
    real(fdp) :: width, height


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(2f, 2f, f, f)', advance='no') centre, xylen,&
         & width, height

    call s_ellipsemask(t_mask, nl, nu, centre, xylen, width)

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_ellipse



  subroutine s_fillellipse(pot, mask, nl, nu)

    !
    ! Define a filled ellipse component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp), dimension(dim) :: centre, xylen
    real(fdp) :: height


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(2f, 2f, f)', advance='no') centre, xylen, height

    call s_ellipsemask(t_mask, nl, nu, centre, xylen,&
         & sqrt(sum(parameters%gridspacing**2)))

    call s_fillmask(t_mask, nl, nu)

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_fillellipse



  subroutine s_line(pot, mask, nl, nu)

    !
    ! Define a line component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp), dimension(dim) :: point1, point2, position
    real(fdp) :: width, height


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(2f, 2f, f, f)', advance='no') point1, point2,&
         & width, height

    call s_linemask(t_mask, nl, nu, point1, point2, width)

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_line



  subroutine s_ramp(pot, mask, nl, nu)

    !
    ! Define a straight line ramp component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    real(fdp), dimension(dim) :: point1, point2, position
    real(fdp) :: width, a, b, c, d, height_to, height_from
    integer(idp) :: i, j


    read(potdata, '(2f, 2f, f, f, f)', advance='no') point1, point2,&
         & width, height_from, height_to

    c = sqrt(sum((point1-point2)**2))

    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          position = f_position( (/ i, j /) )

          a = sqrt(sum((point1-position)**2))

          b = sqrt(sum((point2-position)**2))

          d = sqrt(abs(-a**4+2*(b**2+c**2)*a**2-(b**2-c**2)**2))&
               &/2.0d0/c

          if (abs(a-d)/nm < sqrt(cutoff)) then

             d = a

          end if

          if (abs(b-d)/nm < sqrt(cutoff)) then

             d = b

          end if

          if ( abs(sqrt(a**2-d**2)+sqrt(b**2-d**2)-c)/nm <=&
               & sqrt(cutoff) .and. (d <= width .or. a <= width .or.&
               & b <= width)) then

             pot(i, j) = (height_to-height_from)*sqrt(a**2-d**2)/c&
                  &+height_from

             mask(i, j) = .true.

          end if

       end do

    end do

  end subroutine s_ramp



  subroutine s_polygon(pot, mask, nl, nu)
    
    !
    ! Define an arbitrary shaped polygone component in the potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp) :: height, width
    real(fdp), dimension(dim) :: point1, point2, point3


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(f, f, 2f)', advance='no') width, height, point1

    point2 = point1

    point3 = point1 + 1.0d0

    do while(sum(abs(point1-point3))/nm > sqrt(cutoff))

       read(potdata, '(2f)', advance='no') point3

       call s_linemask(t_mask, nl, nu, point2, point3, width)

       point2 = point3

    end do

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_polygon



  subroutine s_fillpolygon(pot, mask, nl, nu)

    !
    ! Define a filled arbitrary shaped polygon component in the
    ! potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    logical, dimension(:, :), allocatable :: t_mask
    real(fdp) :: height
    real(fdp), dimension(dim) :: point1, point2, point3


    allocate(t_mask(nl(1):nu(1), nl(2):nu(2)))

    t_mask = .false.

    read(potdata, '(f, 2f)', advance='no') height, point1

    point2 = point1

    point3 = point1 + 1.0d0

    do while(sum(abs(point1-point3))/nm > sqrt(cutoff))

       read(potdata, '(2f)', advance='no') point3

       call s_linemask(t_mask, nl, nu, point2, point3,&
            & sqrt(sum(parameters%gridspacing**2)))

       point2 = point3

    end do

    call s_fillmask(t_mask, nl, nu)

    where (t_mask)

       pot = height

       mask = .true.

    end where

    deallocate(t_mask)

  end subroutine s_fillpolygon



  subroutine s_potfn(pot, mask, nl, nu)

    !
    ! Define an arbitrary mathematical function component in the
    ! potential
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: pot
    logical, dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout) ::&
         & mask

    integer(idp) :: i, j
    character(len = maxstringlength) :: func, func_t
    real(fdp), dimension(2+size(constantvalues)) :: t_vals


    read(potdata, '(a)', advance='no') func

    t_vals(3:) = constantvalues

    do j=nl(2), nu(2)

       do i=nl(1), nu(1)

          t_vals(1:2) = f_position( (/ i, j/) )

          pot(i, j) = f_evaluatefn(func, 'X Y'//' '//constants,&
               & t_vals )

       end do

    end do

    where (real(pot*conjg(pot))/nm > sqrt(cutoff))

       mask = .true.

    end where

  end subroutine s_potfn



  subroutine s_bpoty(delx, vecpotfny, vecpotfn2, nl, nu)

    !
    ! Define the vector potential function such that 
    !
    ! A = b*(-y, 0, 0)
    ! B = curl A = (0, 0, b)
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: delx,&
         & vecpotfn2
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: vecpotfny

    complex(cdp), dimension(:, :), allocatable :: tmp
    integer(idp) :: i, j
    real(fdp) :: x, y


    allocate(tmp(nl(1):nu(1), nl(2):nu(2)))

    vecpotfny = 0.0d0

    do j = nl(2), nu(2)

       y = f_positiony(j) - parameters%position(2)
       
       do i = nl(1), nu(1)

          x = f_positionx(i) - parameters%position(1)
          
          vecpotfny(i, j) = - y*parameters%externalb
          
       end do
       
    end do

    !
    ! Blur the vector potential function
    !

    tmp = vecpotfny
    
    call s_smoothedge(tmp, nl, nu, parameters%smooth)

    vecpotfny = real(tmp)

    !
    ! Define the vector potential function squared
    !

    vecpotfn2 = vecpotfny**2

    !
    ! Compute the actual magnetic field from the above vector
    ! potential.  This is done just for a check.
    !

    if (parameters%plotpot) then

       !
       ! Define the first derivatives d/dy and temporarily store it
       ! in delx
       !
       
       call s_momentumfac(delx, nl, nu, (/ 0.0d0, 1.0d0 /)&
            &/dble(product(parameters%n)), (/ 0, 1 /))
       
       tmp = vecpotfny

       call s_c2dfft(tmp, nl, nu, fftfor)

       tmp = tmp*delx

       call s_c2dfft(tmp, nl, nu, fftbac)

       call s_write2drarray(real(-tmp), nl, nu, filename='bfield',&
            & incboundary=.true.)
       call s_write2drarray(vecpotfny, nl, nu, filename='vecpoty',&
            & incboundary=.true.)

    end if

    !
    ! Define the first derivatives d/dx 
    !

    call s_momentumfac(delx, nl, nu, (/ 1.0d0, 0.0d0 /)&
         &/dble(product(parameters%n)), (/ 1, 0 /))

    !
    ! Compute the cross terms
    !

    tmp = vecpotfny

    call s_c2dfft(tmp, nl, nu, fftfor)
    
    tmp = tmp*delx
    
    call s_c2dfft(tmp, nl, nu, fftbac)

    !
    ! Include terms for speed
    !

    delx = delx*ii*hbar
    vecpotfny = vecpotfny/parameters%mass*chargee
    vecpotfn2 = vecpotfn2/2.0d0/parameters%mass*chargee**2+tmp/2.0d0/parameters%mass*ii*chargee*hbar

    deallocate(tmp)

  end subroutine s_bpoty



  subroutine s_bpotyx(delx, dely, vecpotfnx, vecpotfny, vecpotfn2, nl&
       &, nu)

    !
    ! Define the vector potential function such that 
    !
    ! A = b/2*(-y, x, 0)
    ! B = curl A = (0, 0, b)
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: delx, dely,&
         & vecpotfn2
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: vecpotfnx, &
         & vecpotfny

    complex(cdp), dimension(:, :), allocatable :: tmpx, tmpy
    integer(idp) :: i, j
    real(fdp) :: norm, x, y


    allocate(tmpx(nl(1):nu(1), nl(2):nu(2)), tmpy(nl(1):nu(1),&
         & nl(2):nu(2)))

    vecpotfnx = 0.0d0
    vecpotfny = 0.0d0

    do j = nl(2), nu(2)

       y = f_positiony(j) - parameters%position(2)
       
       do i = nl(1), nu(1)

          x = f_positionx(i) - parameters%position(1)
          
          vecpotfny(i, j) = -y*parameters%externalb/2.0d0
          vecpotfnx(i, j) = x*parameters%externalb/2.0d0
          
       end do
       
    end do

    !
    ! Blur the vector potential function
    !

    tmpy = vecpotfny
    tmpx = vecpotfnx

    call s_smoothedge(tmpy, nl, nu, parameters%smooth)
    call s_smoothedge(tmpx, nl, nu, parameters%smooth)

    vecpotfny = real(tmpy)
    vecpotfnx = real(tmpx)
    
    !
    ! Define the vector potential function squared
    !

    vecpotfn2 = vecpotfnx**2 + vecpotfny**2
    
    !
    ! Define the first derivatives d/dx and d/dy
    !
    
    call s_momentumfac(delx, nl, nu, (/ 1.0d0, 0.0d0 /)&
         &/dble(product(parameters%n)), (/ 1, 0 /))
    call s_momentumfac(dely, nl, nu, (/ 0.0d0, 1.0d0 /)&
         &/dble(product(parameters%n)), (/ 0, 1 /))
       
    !
    ! Compute the actual magnetic field from the above vector
    ! potential.  This is done just for a check.
    !
    
    if (parameters%plotpot) then

       tmpy = vecpotfny
       tmpx = vecpotfnx

       call s_c2dfft(tmpy, nl, nu, fftfor)
       call s_c2dfft(tmpx, nl, nu, fftfor)

       tmpy = tmpy*dely
       tmpx = tmpx*delx

       call s_c2dfft(tmpy, nl, nu, fftbac)
       call s_c2dfft(tmpx, nl, nu, fftbac)

       call s_write2drarray(real(tmpx-tmpy), nl, nu, filename&
            &='bfield', incboundary=.true.)
       call s_write2drarray(vecpotfny, nl, nu, filename='vecpoty',&
            & incboundary=.true.)
       call s_write2drarray(vecpotfnx, nl, nu, filename='vecpotx',&
            & incboundary=.true.)

    end if

    !
    ! Compute the cross terms
    !

    tmpy = vecpotfny
    tmpx = vecpotfnx

    call s_c2dfft(tmpy, nl, nu, fftfor)
    call s_c2dfft(tmpx, nl, nu, fftfor)
    
    tmpy = tmpy*delx
    tmpx = tmpx*dely
    
    call s_c2dfft(tmpy, nl, nu, fftbac)
    call s_c2dfft(tmpx, nl, nu, fftbac)

    !
    ! Include terms for speed
    !

    delx = delx*ii*hbar
    dely = dely*ii*hbar
    vecpotfnx = vecpotfnx/parameters%mass*chargee
    vecpotfny = vecpotfny/parameters%mass*chargee
    vecpotfn2 = vecpotfn2/2.0d0/parameters%mass*chargee**2 + (tmpx&
         &+tmpy)/2.0d0/parameters%mass*ii*chargee*hbar

    deallocate(tmpx, tmpy)

  end subroutine s_bpotyx



  subroutine s_bpotxyxy(delxy, vecpotfn, vecpotfn2, nl, nu)

    !
    ! Define the vector potential function such that 
    !
    ! A = b/2*(x-y, x-y, 0)
    ! B = curl A = (0, 0, b)
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: delxy,&
         & vecpotfn2
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: vecpotfn

    complex(cdp), dimension(:, :), allocatable :: tmp
    integer(idp) :: i, j
    real(fdp) :: x, y


    allocate(tmp(nl(1):nu(1), nl(2):nu(2)))

    vecpotfn = 0.0d0
    
    do j = nl(2), nu(2)

       y = f_positiony(j) - parameters%position(2)
       
       do i = nl(1), nu(1)

          x = f_positionx(i) - parameters%position(1)
          
          vecpotfn(i, j) = (x - y)*parameters%externalb/2.0d0
          
       end do
       
    end do
    
    !
    ! Blur the vector potential function
    !

    tmp = vecpotfn

    call s_smoothedge(tmp, nl, nu, parameters%smooth)

    vecpotfn = real(tmp)
    
    !
    ! Define the vector potential function squared
    !

    vecpotfn2 = vecpotfn**2

    !
    ! Compute the actual magnetic field from the above vector
    ! potential.  This is done just for a check.
    !

    if (parameters%plotpot) then

       !
       ! Define the first derivatives d/dx-d/dy
       !
       
       call s_momentumfac(delxy, nl, nu, (/ 1.0d0, -1.0d0 /)&
            &/dble(product(parameters%n)), (/ 1, 1 /))
       
       !
       ! Include normlaisation for speed
       !
       
       tmp = vecpotfn

       call s_c2dfft(tmp, nl, nu, fftfor)

       tmp = tmp*delxy

       call s_c2dfft(tmp, nl, nu, fftbac)

       call s_write2drarray(real(tmp), nl, nu, filename='bfield',&
            & incboundary=.true.)
       call s_write2drarray(vecpotfn, nl, nu, filename='vecpot',&
            & incboundary=.true.)

    end if

    !
    ! Define the first derivatives d/dx+d/dy
    !

    call s_momentumfac(delxy, nl, nu, (/ 1.0d0, 1.0d0 /)&
         &/dble(product(parameters%n)), (/ 1, 1 /))

    !
    ! Compute derivative of vector potential function
    !

    tmp = vecpotfn
    
    call s_c2dfft(tmp, nl, nu, fftfor)
    
    tmp = tmp*delxy
    
    call s_c2dfft(tmp, nl, nu, fftbac)
    
    !
    ! Include terms for speed
    !

    delxy = delxy*ii*hbar
    vecpotfn = vecpotfn/parameters%mass*chargee
    vecpotfn2 = vecpotfn2/parameters%mass*chargee**2 + tmp*ii*hbar&
         &*chargee/2.0d0/parameters%mass

    deallocate(tmp)

  end subroutine s_bpotxyxy



  subroutine s_bpotfft(delf, vecpotfn, vecpotfn2, nl, nu)

    !
    ! Define the vector potential function such that 
    !
    ! A = b/2*(f, f, 0)
    ! B = curl A = (0, 0, b)
    !
    ! by direct solution using fft's.
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: delf,&
         & vecpotfn2
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: vecpotfn

    complex(cdp), dimension(:, :), allocatable :: tmp


    allocate(tmp(nl(1):nu(1), nl(2):nu(2)))

    !
    ! Compute the actual magnetic field with the smoothed boundary
    !

    tmp = parameters%externalb

    call s_smoothedge(tmp, nl, nu, parameters%smooth)

    !
    ! Compute the momentum factor representation for grad operator
    !

    call s_momentumfac(delf, nl, nu, (/ 1.0d0, -1.0d0 /)&
         &/dble(product(parameters%n)), (/ 1, 1 /))

    !
    ! Solve for vector potential
    !

    call s_c2dfft(tmp, nl, nu, fftfor)

    where(abs(delf) < cutoff)

       tmp = cutoff

    elsewhere
       
       tmp = tmp/delf

    end where
    
    call s_c2dfft(tmp, nl, nu, fftbac)

    vecpotfn = real(tmp)
    vecpotfn2 = vecpotfn**2

    !
    ! Plot potential if required
    !

    if (parameters%plotpot) then

       !
       ! Define the first derivatives d/dx-d/dy
       !
       
       call s_momentumfac(delf, nl, nu, (/ 1.0d0, -1.0d0 /)&
            &/dble(product(parameters%n)), (/ 1, 1 /))
       
       !
       ! Include normlaisation for speed
       !
       
       tmp = vecpotfn

       call s_c2dfft(tmp, nl, nu, fftfor)

       tmp = tmp*delf

       call s_c2dfft(tmp, nl, nu, fftbac)

       call s_write2drarray(real(tmp), nl, nu, filename='bfield',&
            & incboundary=.true.)
       call s_write2drarray(vecpotfn, nl, nu, filename='vecpot',&
            & incboundary=.true.)

    end if

    !
    ! Define the first derivatives d/dx+d/dy
    !

    call s_momentumfac(delf, nl, nu, (/ 1.0d0, 1.0d0 /)&
         &/dble(product(parameters%n)), (/ 1, 1 /))

    !
    ! Compute derivative of vector potential function
    !

    tmp = vecpotfn
    
    call s_c2dfft(tmp, nl, nu, fftfor)
    
    tmp = tmp*delf
    
    call s_c2dfft(tmp, nl, nu, fftbac)
    
    !
    ! Include terms for speed
    !

    delf = delf*ii*hbar
    vecpotfn = vecpotfn/parameters%mass*chargee
    vecpotfn2 = vecpotfn2/parameters%mass*chargee**2 + tmp*ii*hbar&
         &*chargee/2.0d0/parameters%mass

    deallocate(tmp)

  end subroutine s_bpotfft



end module potential
