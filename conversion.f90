module conversion

  !
  ! This modules handles conversions of all sorts
  !

  use precision
  use globals

  implicit none



contains



  function f_numchar(num) result(numchar)

    !
    ! Integer to a 4 digit ASCII string
    !

    integer(idp), intent(in) :: num
    character(len=numcharlen) :: numchar


    write(numchar, '(i6.6)') num

  end function f_numchar



  function f_pmax(gridspacing) result(pmax)

    !
    ! Maximum momentum from gridspacing
    !

    real(fdp), dimension(dim), intent(in) :: gridspacing
    real(fdp), dimension(dim) :: pmax


    pmax = pi/gridspacing*hbar

  end function f_pmax



  function f_gridspacing(pmax) result(gridspacing)

    !
    ! Gridspacing from maximum momentum
    !

    real(fdp), dimension(dim), intent(in) :: pmax
    real(fdp), dimension(dim) :: gridspacing


    gridspacing = pi/pmax*hbar

  end function f_gridspacing



  function f_gridpointx(x) result(gridpoint)

    !
    ! Y grid point corresponding to x coordinate
    !

    real(fdp), intent(in) :: x
    integer(idp) :: gridpoint


    gridpoint = nint(x/parameters%gridspacing(1) + dble(parameters&
         &%n(1)+1)/2.0d0)

  end function f_gridpointx



  function f_gridpointy(y) result(gridpoint)

    !
    ! Y gridpoint corresponding to y coordinate
    !

    real(fdp), intent(in) :: y
    integer(idp) :: gridpoint


    gridpoint = nint(y/parameters%gridspacing(2) + dble(parameters&
         &%n(2)+1)/2.0d0)

  end function f_gridpointy



  function f_gridpoint(position) result(gridpoint)

    !
    ! Gridpoint corresponding to coordinate

    real(fdp), dimension(dim), intent(in) :: position
    integer(idp), dimension(dim) :: gridpoint


    gridpoint = (/ f_gridpointx(position(1)),&
         & f_gridpointy(position(2)) /)

  end function f_gridpoint



  function f_positionx(gridpointx) result(x)

    !
    ! X position corresponding to x gripoint
    !

    integer(idp), intent(in) :: gridpointx
    real(fdp) :: x


    x = (dble(gridpointx)-dble(parameters%n(1)+1)/2.0d0)*parameters&
         &%gridspacing(1)

  end function f_positionx



  function f_positiony(gridpointy) result(y)

    !
    ! Y position corresponding to y gripoint
    !

    integer(idp), intent(in) :: gridpointy
    real(fdp) :: y


    y = (dble(gridpointy)-dble(parameters%n(2)+1)/2.0d0)*parameters&
         &%gridspacing(2)

  end function f_positiony



  function f_position(gridpoint) result(position)

    !
    ! Position corresponding to a gripoint
    !

    integer(idp), dimension(dim), intent(in) :: gridpoint
    real(fdp), dimension(dim) :: position


    position = (/ f_positionx(gridpoint(1)),&
         & f_positiony(gridpoint(2)) /)

  end function f_position



  function f_momentumx(x) result(p)

    !
    ! X momentum corresponding to x gridpoint
    !

    integer(idp), intent(in) :: x
    real(fdp) :: p

    !
    ! Must align the grid to what FFT computes
    !

    p = dble(x-1-parameters%n(1)/2)*parameters%pspacing(1)

  end function f_momentumx



  function f_momentumy(y) result(p)

    !
    ! Y momentum corresponding to y gridpoint
    !

    integer(idp), intent(in) :: y
    real(fdp) :: p


    !
    ! Must align the grid to what FFT computes
    !

    p = dble(y-1-parameters%n(2)/2)*parameters%pspacing(2)

  end function f_momentumy



  function f_momentum(point) result(p)

    !
    ! Momentum corresponding to a gridpoint
    !

    integer(idp), dimension(dim), intent(in) :: point
    real(fdp), dimension(dim) :: p


    p = (/ f_momentumx(point(1)), f_momentumy(point(2)) /)

  end function f_momentum



  function f_pgridpointx(x) result(gridpoint)

    !
    ! X gridpoint corresponding to x momentum
    !

    real(fdp), intent(in) :: x
    integer(idp) :: gridpoint


    gridpoint = nint(x/parameters%pspacing(1))+parameters%n(1)/2+1

  end function f_pgridpointx



  function f_pgridpointy(y) result(gridpoint)

    !
    ! Y gridpoint corresponding to y momentum
    !

    real(fdp), intent(in) :: y
    integer(idp) :: gridpoint


    gridpoint = nint(y/parameters%pspacing(2))+parameters%n(2)/2+1

  end function f_pgridpointy



  function f_pgridpoint(position) result(gridpoint)

    !
    ! Gridpoint corresponding to momentum
    !

    real(fdp), dimension(dim), intent(in) :: position
    integer(idp), dimension(dim) :: gridpoint


    gridpoint = (/ f_gridpointx(position(1)),&
         & f_gridpointy(position(2)) /)

  end function f_pgridpoint



end module conversion
