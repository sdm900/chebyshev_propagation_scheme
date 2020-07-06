module fft

  use precision
  use mpi
  use globals

  implicit none

  interface

#ifdef FFTW
     subroutine fftwnd_f77_create_plan(plan, dim, n, dir, flags)
       integer(PNTSIZE), intent(inout) :: plan
       integer(4), intent(in) :: dim, dir, flags
       integer(4), dimension(*), intent(in) :: n
     end subroutine fftwnd_f77_create_plan

     subroutine fftwnd_f77_one(plan, in, out)
       integer(PNTSIZE), intent(in) :: plan
       complex(8), dimension(*), intent(in) :: in
       complex(8), dimension(*), intent(out) :: out
     end subroutine fftwnd_f77_one

     subroutine fftwnd_f77_destroy_plan(plan)
       integer(PNTSIZE) :: plan
     end subroutine fftwnd_f77_destroy_plan
#endif

#ifdef MPI
     subroutine fftwnd_f77_mpi_create_plan(plan, mpicomm, dim, n, dir&
          &, flags)
       integer(PNTSIZE), intent(inout) :: plan
       integer(4), intent(in) :: mpicomm, dim, dir, flags
       integer(4), dimension(*), intent(in) :: n
     end subroutine fftwnd_f77_mpi_create_plan
       
     subroutine fftwnd_f77_mpi(plan, nfields, data, work, usework,&
          & ioorder)
       integer(PNTSIZE), intent(in) :: plan
       integer(4), intent(in) :: nfields, usework, ioorder
       complex(8), dimension(*), intent(inout) :: data, work
     end subroutine fftwnd_f77_mpi

     subroutine fftwnd_f77_mpi_destroy_plan(plan)
       integer(PNTSIZE), intent(in) :: plan
     end subroutine fftwnd_f77_mpi_destroy_plan
#endif

#ifdef SSL2
     subroutine dvmcft(datareal, dataimag, n, dim, isn, w, lenw, icon)
       real(8), dimension(*), intent(inout) :: datareal, dataimag, w
       integer(4), dimension(*), intent(in) :: n, isn
       integer(4), intent(in) :: dim, lenw, icon
     end subroutine dvmcft
#endif

#ifdef CXML
     function zfft_2d(infmt, outfmt, dir, in, out, nx, ny, lda, nxs,&
          & nys) result(fft)
       integer(4) :: fft
       integer(4), intent(in) :: nx, ny, nxs, nys, lda
       complex(8), dimension(*), intent(in) :: in
       complex(8), dimension(*), intent(out) :: out
       character(*), intent(in) :: infmt, outfmt, dir
     end function zfft_2d
#endif

  end interface



contains



#ifdef SSL2
  subroutine s_c2dfft(data, n, isign)

    ! This routine implements an optimised fourier transform on the
    ! VPP Super Computer
    !
    ! The result is an un-normalised fourier transform.

    integer(idp), dimension(dim), intent(in) :: n
    complex(cdp), dimension(n(1), n(2)), intent(inout) :: data
    integer(idp), intent(in) :: isign

    integer(idp) :: lenw, icon, i, j
    integer(idp), dimension(dim) :: t_n, isn
    real(fdp), dimension(n(1), n(2)) :: datareal1, dataimag1
    real(fdp), dimension(n(2), n(1)) :: datareal2, dataimag2
    real(fdp), dimension(:), allocatable :: w
    real(fdp), dimension(1) :: wdummy


    t_n = (/ n(2), n(1) /)
    isn = (/ 0 , -isign /)

    datareal1 = real(data)
    dataimag1 = aimag(data)

    lenw = 0

    call dvmcft(datareal1,dataimag1,n,dim,isn,wdummy,lenw,icon)

    allocate(w(lenw))

    call dvmcft(datareal1,dataimag1,n,dim,isn,w,lenw,icon)

    deallocate(w)

    if (icon /= 0) then

       call s_ffterror(icon)

    end if

    do i = 1, t_n(1)
       do j = 1, t_n(2)

          datareal2(i, j) = datareal1(j, i)
          dataimag2(i, j) = dataimag1(j, i)

       end do
    end do

    lenw = 0

    call dvmcft(datareal2,dataimag2,t_n,dim,isn,wdummy,lenw,icon)

    allocate(w(lenw))

    call dvmcft(datareal2,dataimag2,t_n,dim,isn,w,lenw,icon)

    deallocate(w)

    call s_vppffterror(icon)

    do j = 1, n(2)
       do i = 1, n(1)

          data(i, j) = cmplx(datareal2(j, i), dataimag2(j, i))

       end do
    end do

  end subroutine s_c2dfft



  subroutine s_c1dfft(data, n, isign)

    ! This routine implements an optimised fourier transform on the
    ! VPP Super Computer
    !
    ! The result is an un-normalised fourier transform.

    integer(idp), intent(in) :: n, isign
    complex(cdp), dimension(n), intent(inout) :: data

    integer(idp) :: lenw, icon
    real(fdp), dimension(n) :: datareal, dataimag
    real(fdp), dimension(:), allocatable :: w
    real(fdp), dimension(1) :: wdummy


    datareal1 = real(data)
    dataimag1 = aimag(data)

    lenw = 0

    call dvmcft(datareal,dataimag,n,dim,-isign,wdummy,lenw,icon)

    allocate(w(lenw))

    call dvmcft(datareal,dataimag,n,dim,-isign,w,lenw,icon)

    deallocate(w)

    if (icon /= 0) then

       call s_vppffterror(icon)

    end if

  end subroutine s_c1dfft



  subroutine s_vppffterror(num)

    integer(idp), intent(in) :: num


    if (icon /= 0) then

       select case (icon)
       case(30000)
          write(*,*) 'Fourier transform ended abnormally (M=<0)'
       case(30001)
          write(*,*) 'Fourier transform ended abnormally (The size of&
               & the work array is too small)'
       case(30002)
          write(*,*) 'Fourier transform ended abnormally (For a i,&
               & ISN(i)>1 or ISN(i)<-1)'
       case(30003)
          write(*,*) 'Fourier transform ended abnormally (For a i&
               &,N(i)<1)'
       case(30004)
          write(*,*) 'Fourier transform ended abnormally (For any i,&
               & ISN(i)=0)'
       end select

    end if

  end subroutine s_vppffterror
#endif



#ifdef FFTW
  subroutine s_c2dfft(data, nl, nu, isign)

    !
    ! This routine implements an optimised complex fft for a DEC
    ! Alpha workstation
    !
    ! The result is an unnormalised fourier transform.
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: data
    integer(idp), intent(in) :: isign

    integer(pnt), target, save :: planfftfor = 0, planfftbac = 0
    integer(pnt), pointer :: plan
    integer(idp), dimension(2), save :: n = 0, nlold = 0, nuold = 0
    integer(idp), parameter :: FFTW_ESTIMATE=0, FFTW_MEASURE=1,&
         & FFTW_IN_PLACE=8, FFTW_USE_WISDOM=16


    n = nu - nl + 1

    if (isign == fftfor) then

       plan => planfftfor

    elseif (isign == fftbac) then

       plan => planfftbac

    end if

    if ((nlold(1) /= nl(1) .or. nuold(1) /= nu(1) .or. nlold(2) /=&
         & nl(2) .or. nuold(2) /= nu(2)) .and. plan /= 0) then

       call fftwnd_f77_destroy_plan(plan)

       plan = 0

    end if

    if (plan == 0) then

       call fftwnd_f77_create_plan(plan, 2, n, isign, FFTW_MEASURE&
            &+FFTW_IN_PLACE+FFTW_USE_WISDOM)

       nlold = nl
       nuold = nu

    end if

    call fftwnd_f77_one(plan, data, data)

  end subroutine s_c2dfft



  subroutine s_c1dfft(data, n, isign)

    ! This routine implements an optimised complex fft for a DEC
    ! Alpha workstation
    !
    ! The result is an unnormalised fourier transform.

    integer(idp), intent(in) :: n
    complex(cdp), dimension(n), intent(inout) :: data
    integer(idp), intent(in) :: isign

    integer(pnt), target, save :: planfftfor = 0, planfftbac = 0
    integer(pnt), pointer :: plan
    integer(idp), save :: nold
    integer(idp), parameter :: FFTW_ESTIMATE=0, FFTW_MEASURE=1,&
         & FFTW_IN_PLACE=8, FFTW_USE_WISDOM=16


    if (isign == fftfor) then

       plan => planfftfor

    elseif (isign == fftbac) then

       plan => planfftbac

    end if

    if (nold /= n .and. plan /= 0) then

       call fftwnd_f77_destroy_plan(plan)

       plan = 0

    end if

    if (plan == 0) then

       call fftwnd_f77_create_plan(plan, 1, (/ n /), isign,&
            & FFTW_MEASURE+FFTW_IN_PLACE+FFTW_USE_WISDOM)

       nold = n

    end if

    call fftwnd_f77_one(plan, data, data)

  end subroutine s_c1dfft
#endif



#ifdef CXML
  subroutine s_c2dfft(data, n, isign)

    ! This routine implements a complex fourier transform in 2d
    ! using Digitals Extended Math Library (dxml).
    !
    ! the results is an unnormalised fourier transform

    integer(idp), dimension(2), intent(in) :: n
    complex(cdp), dimension(n(1), n(2)), intent(inout) :: data
    integer(idp), intent(in) :: isign

    integer(idp) :: status


    status = -1

    if (isign == 1) then

       status = zfft_2d('C', 'C', 'B', data, data, n(1), n(2), n(1) ,&
            & 1, 1)

    else if (isign == -1) then

       status = zfft_2d('C', 'C', 'F', data, data, n(1), n(2), n(1) ,&
            & 1, 1)

    end if

    if (status /= 0) then

       write(*,*) 'Error with the zfft_2d = ', status

    end if

  end subroutine s_c2dfft



  subroutine s_c1dfft(data, n, isign)

    ! This routine implements a complex fourier transform in 2d
    ! using Digitals Extended Math Library (dxml).
    !
    ! the results is an unnormalised fourier transform

    integer(idp), intent(in) :: n
    complex(cdp), dimension(n), intent(inout) :: data
    integer(idp), intent(in) :: isign

    integer(idp) :: status


    status = -1

    if (isign == 1) then

       status = zfft('C', 'C', 'B', data, data, n, 1)

    else if (isign == -1) then

       status = zfft('C', 'C', 'F', data, data, n, 1)

    end if

    if (status /= 0) then

       write(*,*) 'Error with the zfft = ', status

    end if

  end subroutine s_c1dfft
#endif



#ifdef MPI
  subroutine s_c2dfft(data, nl, nu, isign)

    !
    ! This routine implements an optimised complex fft for a DEC
    ! Alpha workstation
    !
    ! The result is an unnormalised fourier transform.
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: data
    integer(idp), intent(in) :: isign

    integer(pnt), target, save :: planfor = 0, planbac = 0
    integer(pnt), pointer :: plan
    integer(idp), dimension(2), save :: nlold = 0, nuold = 0, n = 0
    integer(idp), parameter :: FFTW_ESTIMATE=0, FFTW_MEASURE=1,&
         & FFTW_USE_WISDOM=16, FFTW_TRANSPOSED_ORDER=1,&
         & FFTW_NORMAL_ORDER=0
    complex(cdp), dimension(:), allocatable, save :: work


    if (isign == fftfor) then

       plan => planfor
       n = parameters%n

    elseif (isign == fftbac) then

       plan => planbac
       n = (/ parameters%n(2), parameters%n(1) /)

    end if

    if ((nlold(1) /= nl(1) .or. nuold(1) /= nu(1) .or. nlold(2) /=&
         & nl(2) .or. nuold(2) /= nu(2)) .and. plan /= 0) then

       call fftwnd_f77_mpi_destroy_plan(plan)

       plan = 0

       deallocate(work)

    end if

    if (plan == 0) then

       call fftwnd_f77_mpi_create_plan(plan, mpicomm, 2, n, isign,&
            & FFTW_MEASURE+FFTW_USE_WISDOM)

       nlold = nl
       nuold = nu

    end if

    if (.not. allocated(work)) then

       allocate(work(size(data)))

    end if

    call fftwnd_f77_mpi(plan, 1, data, work, 1, FFTW_TRANSPOSED_ORDER)

  end subroutine s_c2dfft
#endif


end module fft
