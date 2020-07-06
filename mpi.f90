module mpi

  use precision
  use globals

  implicit none

#ifdef MPI
  include 'mpif.h'

  integer, parameter :: inttag=1000, realtag=2000, comptag=3000
  integer, save :: mpicomm=-1, mpisize = 1, mpirank=0, mpistart = 0, mpiend = 0, mpierror=0


contains


  subroutine s_mpierror (error)

    character(*) :: error

    
    if (mpierror /= 0) then
       
       write(*,*) '**********'
       write(*,*) 'Error: ', error
       write(*,*) 'MPI Error: ', mpierror
       write(*,*) 'Comm: ', mpicomm
       write(*,*) 'Size: ', mpisize
       write(*,*) 'Rank: ', mpirank
       write(*,*) '**********'

    end if

  end subroutine s_mpierror

  

  subroutine s_mpiinit()

    call mpi_init(mpierror)
    
    call s_mpierror('mpi_init')
    
    mpicomm = MPI_COMM_WORLD
    
    call mpi_comm_size(mpicomm, mpisize, mpierror)
    
    call s_mpierror('mpi_comm_size')
    
    call mpi_comm_rank(mpicomm, mpirank, mpierror)
    
    call s_mpierror('mpi_comm_rank')

    mpistart = 0

    mpiend = mpisize - 1
    
  end subroutine s_mpiinit

  

  subroutine s_mpifinalise()

    call mpi_finalize(mpierror)

    call s_mpierror('mpi_finalize')

  end subroutine s_mpifinalise



  subroutine s_mpiwait(mpirequest, mpistatus)

    integer :: mpirequest
    integer, dimension(MPI_STATUS_SIZE) :: mpistatus

    call mpi_wait(mpirequest, mpistatus, mpierror)

    call s_mpierror('Failed wait')

  end subroutine s_mpiwait



  subroutine s_mpiwaitany(sizempirequests, mpirequests, mpiindex, mpistatus)
    
    integer :: sizempirequests
    integer, dimension(sizempirequests) :: mpirequests
    integer, dimension(MPI_STATUS_SIZE) :: mpistatus
    integer :: mpiindex

    call mpi_waitany(sizempirequests, mpirequests, mpiindex, mpistatus, mpierror)

    call s_mpierror('Failed waitany')

  end subroutine s_mpiwaitany



  subroutine s_mpiwaitall(sizempirequests, mpirequests, mpistatuss)

    integer :: sizempirequests
    integer, dimension(sizempirequests) :: mpirequests
    integer, dimension(MPI_STATUS_SIZE, sizempirequests) :: mpistatuss

    call mpi_waitall(sizempirequests, mpirequests, mpistatuss, mpierror)

    call s_mpierror('Failed waitall')

  end subroutine s_mpiwaitall



  subroutine s_mpisendint(ints, sizeints, sendto, mpirequest)

    integer :: sizeints
    integer, dimension(sizeints) :: ints
    integer :: sendto
    integer :: mpirequest

    
    call mpi_issend(ints, sizeints, MPI_INTEGER, sendto, inttag, mpicomm, mpirequest, mpierror)
    
    call s_mpierror('Failed send')

  end subroutine s_mpisendint



  subroutine s_mpirecvint(ints, sizeints, recvfrom, mpirequest)

    integer :: sizeints
    integer, dimension(sizeints) :: ints
    integer :: recvfrom
    integer :: mpirequest


    call mpi_irecv(ints, sizeints, MPI_INTEGER, recvfrom, inttag, mpicomm, mpirequest, mpierror)
    
    call s_mpierror('Failed recv')

  end subroutine s_mpirecvint



  subroutine s_mpisendreal(reals, sizereals, sendto, mpirequest)

    integer :: sizereals
    real(8), dimension(sizereals) :: reals
    integer :: sendto
    integer :: mpirequest

    
    call mpi_issend(reals, sizereals, MPI_DOUBLE_PRECISION, sendto, realtag, mpicomm, mpirequest, mpierror)
    
    call s_mpierror('Failed send')

  end subroutine s_mpisendreal



  subroutine s_mpirecvreal(reals, sizereals, recvfrom, mpirequest)

    integer :: sizereals
    real(8), dimension(sizereals) :: reals
    integer :: recvfrom
    integer :: mpirequest


    call mpi_irecv(reals, sizereals, MPI_DOUBLE_PRECISION, recvfrom, realtag, mpirequest, mpierror)
    
    call s_mpierror('Failed recv')

  end subroutine s_mpirecvreal



  subroutine s_mpisendcomp(complexs, sizecomplexs, sendto, mpirequest)

    integer :: sizecomplexs
    complex(8), dimension(sizecomplexs) :: complexs
    integer :: sendto
    integer :: mpirequest

    
    call mpi_issend(complexs, sizecomplexs, MPI_DOUBLE_COMPLEX, sendto, comptag, mpicomm, mpirequest, mpierror)
    
    call s_mpierror('Failed send')

  end subroutine s_mpisendcomp



  subroutine s_mpirecvcomp(complexs, sizecomplexs, recvfrom, mpirequest)

    integer :: sizecomplexs
    complex(8), dimension(sizecomplexs) :: complexs
    integer :: recvfrom
    integer :: mpirequest


    call mpi_irecv(complexs, sizecomplexs, MPI_DOUBLE_COMPLEX, recvfrom, comptag, mpicomm, mpirequest, mpierror)
    
    call s_mpierror('Failed recv')

  end subroutine s_mpirecvcomp



  subroutine s_mpiarraysize(nl, nu, trnl, trnu)

    integer(4), dimension(dim) :: nl, nu, trnl, trnu


    if ((parameters%n(1)/mpisize)*mpisize == parameters%n(1) .and. (parameters%n(2)/mpisize)*mpisize == parameters%n(2)) then
       
       nu(1) = parameters%n(1)
       nl(1) = 1
       trnu(1) = parameters%n(2)
       trnl(1) = 1

       nu(2) = parameters%n(2) / mpisize * mpirank + parameters%n(2) / mpisize
       nl(2) = parameters%n(2) / mpisize * mpirank + 1
       trnu(2) = parameters%n(1) / mpisize * mpirank + parameters%n(1) / mpisize
       trnl(2) = parameters%n(1) / mpisize * mpirank + 1

    else

       call s_mpierror('Can not partition array')

    end if

  end subroutine s_mpiarraysize



  subroutine s_mpiallreducesumint(number)

    integer :: number

    integer :: tmp


    call mpi_allreduce(number, tmp, 1, MPI_INTEGER, MPI_SUM, mpicomm, mpierror)

    call s_mpierror('Failed reducesum')

    number = tmp

  end subroutine s_mpiallreducesumint



  subroutine s_mpiallreducesumreal(number)

    real(8) :: number

    real(8) :: tmp


    call mpi_allreduce(number, tmp, 1, MPI_DOUBLE_PRECISION, MPI_SUM, mpicomm, mpierror)

    call s_mpierror('Failed reducesum')

    number = tmp

  end subroutine s_mpiallreducesumreal



  subroutine s_mpiallreducemaxreal(number)

    real(8) :: number

    real(8) :: tmp


    call mpi_allreduce(number, tmp, 1, MPI_DOUBLE_PRECISION, MPI_MAX, mpicomm, mpierror)

    call s_mpierror('Failed reducesum')

    number = tmp

  end subroutine s_mpiallreducemaxreal



  subroutine s_mpiallreduceminreal(number)

    real(8) :: number

    real(8) :: tmp


    call mpi_allreduce(number, tmp, 1, MPI_DOUBLE_PRECISION, MPI_MIN, mpicomm, mpierror)

    call s_mpierror('Failed reducesum')

    number = tmp

  end subroutine s_mpiallreduceminreal



  subroutine s_mpiallgatherint(intv, sizeintv, inta, sizeinta)
    
    integer :: sizeintv, sizeinta
    integer, dimension(sizeintv) :: intv
    integer, dimension(sizeinta) :: inta


    call mpi_allgather(intv, sizeintv, MPI_INTEGER, inta, sizeintv, MPI_INTEGER, mpistart, mpicomm, mpierror)
    
    call s_mpierror('Failed gatherreal')
    
  end subroutine s_mpiallgatherint



  subroutine s_mpigatherreal(array, sizearray, totarray, sizetotarray)
    
    integer :: sizearray, sizetotarray
    real(8), dimension(sizearray) :: array
    real(8), dimension(sizetotarray) :: totarray


    call mpi_gather(array, sizearray, MPI_DOUBLE_PRECISION, totarray, sizearray, MPI_DOUBLE_PRECISION, mpistart, mpicomm, mpierror)
    
    call s_mpierror('Failed gatherreal')
    
  end subroutine s_mpigatherreal



  subroutine s_mpigathercomp(array, sizearray, totarray, sizetotarray)
    
    integer :: sizearray, sizetotarray
    complex(8), dimension(sizearray) :: array
    complex(8), dimension(sizetotarray) :: totarray


    call mpi_gather(array, sizearray, MPI_DOUBLE_COMPLEX, totarray, sizearray, MPI_DOUBLE_COMPLEX, mpistart, mpicomm, mpierror)
    
    call s_mpierror('Failed gathercomp')
    
  end subroutine s_mpigathercomp



  subroutine s_mpiscatterreal(totarray, sizetotarray, array, sizearray)
    
    integer :: sizearray, sizetotarray
    real(8), dimension(sizearray) :: array
    real(8), dimension(sizetotarray) :: totarray


    call mpi_scatter(totarray, sizearray, MPI_DOUBLE_PRECISION, array, sizearray, MPI_DOUBLE_PRECISION, mpistart, mpicomm, mpierror)
    
    call s_mpierror('Failed scatterreal')
    
  end subroutine s_mpiscatterreal



  subroutine s_mpiscattercomp(totarray, sizetotarray, array, sizearray)
    
    integer :: sizearray, sizetotarray
    complex(8), dimension(sizearray) :: array
    complex(8), dimension(sizetotarray) :: totarray


    call mpi_scatter(totarray, sizearray, MPI_DOUBLE_COMPLEX, array, sizearray, MPI_DOUBLE_COMPLEX, mpistart, mpicomm, mpierror)
    
    call s_mpierror('Failed scattercomplex')
    
  end subroutine s_mpiscattercomp



  subroutine s_mpifileopenwrite(filename, handle)

    character(len=*) :: filename
    integer :: handle
    
    integer :: mpiinfo

    
    call mpi_file_open(mpicomm, filename, ior(MPI_MODE_WRONLY, MPI_MODE_CREATE), mpiinfo, handle, mpierror)

    call s_mpierror('Failed fileopen')

  end subroutine s_mpifileopenwrite



  subroutine s_mpifileclose(handle)

    integer :: handle

    
    call mpi_file_close(handle, mpierror)

    call s_mpierror('Failed fileclose')

  end subroutine s_mpifileclose



  subroutine s_mpifilewriteint(handle, intv, sizeintv)
    
    integer :: sizeintv
    integer :: handle
    integer, dimension(sizeintv) :: intv

    integer, dimension(MPI_STATUS_SIZE) :: mpistatus


    call mpi_file_write_shared(handle, intv, sizeintv, MPI_INTEGER, mpistatus, mpierror)

    call s_mpierror('Filed writeshared')

  end subroutine s_mpifilewriteint



  subroutine s_mpifilewriteorderedreal(handle, array, sizearray)

    integer :: sizearray
    integer :: handle
    real(8), dimension(sizearray) :: array

    integer, dimension(MPI_STATUS_SIZE) :: mpistatus
 

    call mpi_file_write_ordered(handle, array, sizearray, MPI_DOUBLE_PRECISION, mpistatus, mpierror)

    call s_mpierror('Failed writeordered')

  end subroutine s_mpifilewriteorderedreal



  subroutine s_mpifilewriteorderedcomp(handle, array, sizearray)

    integer :: sizearray
    integer :: handle
    complex(8), dimension(sizearray) :: array

    integer, dimension(MPI_STATUS_SIZE) :: mpistatus
 

    call mpi_file_write_ordered(handle, array, sizearray, MPI_DOUBLE_COMPLEX, mpistatus, mpierror)

    call s_mpierror('Failed writeordered')

  end subroutine s_mpifilewriteorderedcomp



  subroutine s_mpialltoallcomp(arraysend, arrayrecv, sizearray, chunksize)

    integer(4) :: sizearray
    complex(8), dimension(sizearray) :: arraysend, arrayrecv
    integer :: chunksize


    call mpi_alltoall(arraysend, chunksize, MPI_DOUBLE_COMPLEX, arrayrecv, chunksize, MPI_DOUBLE_COMPLEX, mpicomm, mpierror)

    call s_mpierror('Failed alltoall')

  end subroutine s_mpialltoallcomp



  subroutine s_mpibarrier()

    call mpi_barrier(mpicomm, mpierror)

    call s_mpierror('Failed barrier')

  end subroutine s_mpibarrier


  subroutine s_mpiabort()

    call mpi_abort(mpicomm, mpierror)

    call s_mpierror('Failed abort')

  end subroutine s_mpiabort




  subroutine s_mpitranspose(data, nl, nu)

    !
    ! This routine implements an optimised MPI tranpose, while
    ! leaving the fortran array the same size and dimensions
    !

    integer(4), dimension(dim), intent(in) :: nl, nu
    complex(8), dimension((nu(1)-nl(1)+1)*(nu(2)-nl(2)+1)), intent(inout) :: data

    complex(8), dimension(:), allocatable :: work
    integer(4) :: i, j, k, l, chunksize
    integer(4), dimension(dim) :: n

    
    n = (/ nu(1)-nl(1)+1, nu(2)-nl(2)+1 /)

    allocate(work(size(data)))

    chunksize = size(data)/mpisize

    l = 0

    do i = 1, n(1)
       do j = 1, n(2)

          l = l + 1

          work(l) = data(i+n(1)*(j-1))

       end do
    end do

    call s_mpialltoallcomp(work(1), data(1), size(data), chunksize)

    l = 0

    do i = 1, n(1)/mpisize
       do k = 1, mpisize
          do j = 1, n(2)
             
             l = l + 1
             
             work(l) = data(j+chunksize*(k-1)+n(2)*(i-1))
             
          end do
       end do
    end do

    data = work

    deallocate(work)

  end subroutine s_mpitranspose
#endif


end module mpi
