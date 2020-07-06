program test

  use precision
  use globals
  use fft
  use misc


  integer(idp) :: i, j, k
  integer(idp), dimension(2) :: n
  real(fdp), dimension(2) :: t
  complex(cdp), dimension(:,:), allocatable :: a, dely, delx, xdelydel, array_t, x, y, t1, t2


  read(*,*) n

  parameters%n = n
  parameters%pspacing = 1.0d0/dble(n)
  parameters%gridspacing = 1.0

  allocate(xdelydel(n(1), n(2)), dely(n(1), n(2)), delx(n(1), n(2)), array_t(n(1), n(2)), a(n(1), n(2)), x(n(1), n(2)), y(n(1), n(2)), t1(n(1), n(2)), t2(n(1), n(2)))
  
  do j=1, n(2)
     do i=1, n(1)

        a(i, j) = exp(-dble(i-n(1)/3)**2/50.0d0-dble(j-n(2)/3)**2/50)

     end do
  end do

  call s_momentumfac(xdelydel, (/ 1.0d0, 1.0d0 /), (/ 1, 1 /()
  
  call s_momentumfac(array_t, (/ 1.0d0, 1.0d0 /), (/ 1, 1 /))
  
  call s_momentumfac(dely, (/ 1.0d0, 1.0d0 /), (/ 1, 1 /))
  
  call s_momentumfac(delx, (/ 1.0d0, 1.0d0 /), (/ 1, 1 /))
  
  write(*,*) "a ",sum(abs(a))
  write(*,*) "array_t ",sum(abs(array_t))
  write(*,*) "xdelydel ",sum(abs(xdelydel))
  write(*,*) "delx ",sum(abs(delx))
  write(*,*) "dely ",sum(abs(dely))

  do i=1, n(1)
     
     xdelydel(i, :) = f_positionx(i)*xdelydel(i, :)
     x(i, :) = f_positionx(i)
     
  end do
  
  do i=1, n(2)
     
     xdelydel(:, i) = xdelydel(:, i)-f_positiony(i)*array_t(:, i)
     y(:, i) = f_positiony(i)
     
  end do
  
  xdelydel = xdelydel / dble(n(1)*n(2))
  delx = delx /  dble(n(1)*n(2))
  dely = dely /  dble(n(1)*n(2))
   
  call s_c2dfft(a, n, fftfor)

  array_t = a*xdelydel
  call s_c2dfft(array_t, n, fftbac)

  t1 = dely*a
  call s_c2dfft(t1, n, fftbac)

  t2 = delx*a
  call s_c2dfft(t2, n, fftbac)

  t2 = x*dely*a-y*delx*a
  call s_c2dfft(t2, n, fftbac)

  t1 = x*t1-y*t2

  write(*,*) "t1-t2 ",sum(abs(t2-t1))

  call s_c2dfft(a, n, fftbac)
  a = a / dble(n(1)*n(2))

  write(*,*) "a ",sum(abs(a))
  write(*,*) "array_t ",sum(abs(array_t))
  write(*,*) "t1 ",sum(abs(t1))
  write(*,*) "array_t-t1",sum(abs(array_t-t1))
     
end program test

