module bessel

  use precision
  use globals
  use nswclib

  implicit none


  interface f_besselj
     
     module procedure f_dbesselj
     module procedure f_zbesselj

  end interface
     

  interface f_besseli
     
     module procedure f_dbesseli
     module procedure f_zbesseli

  end interface
     

contains



  function f_dbesselj(n, z) result(bessel)

    integer(idp), intent(in) :: n
    real(fdp), intent(in) :: z

    real(fdp) :: bessel
    complex(cdp) :: bessel_t, z_t


    z_t = z

    call bsslj(z_t, n, bessel_t)

    bessel = real(bessel_t)

  end function f_dbesselj



   function f_zbesselj(n, z) result(bessel)

    integer(idp), intent(in) :: n
    complex(cdp), intent(in) :: z

    complex(cdp) :: bessel


    call bsslj(z, n, bessel)

  end function f_zbesselj



  function f_dbesseli(n, z) result(bessel)

    integer(idp), intent(in) :: n
    real(fdp), intent(in) :: z

    real(fdp) :: bessel
    complex(cdp) :: bessel_t, z_t


    z_t = z

    call bssli(0, z_t, n, bessel_t)

    bessel = real(bessel_t)

  end function f_dbesseli



  function f_zbesseli(n, z) result(bessel)

    integer(idp), intent(in) :: n
    complex(cdp), intent(in) :: z

    complex(cdp) :: bessel, z_t


    call bssli(0, z, n, bessel)

  end function f_zbesseli



end module bessel
