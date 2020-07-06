module misc

  use precision
  use globals
  use interpreter
  use fft
  use mpi
  use conversion

  implicit none



contains



  subroutine s_checkpoint(tstep, wavefn, wavefn0, nl, nu, energy0, espectrum)

    integer(idp), intent(in) :: tstep
    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) :: wavefn, wavefn0
    complex(cdp), dimension(:), intent(in) :: espectrum
    real(fdp), intent(in) :: energy0

    character(len=numcharlen) :: numchar
    character(len=maxstringlength) :: potline
    integer :: readerr


#ifdef MPI
    numchar = f_numchar(mpirank)
#else
    numchar = f_numchar(0)
#endif

    open(unit=chkptfile, file=trim(parameters%chkpntfile)//numchar, form="unformatted", action="write")
#ifdef MPI
    write(chkptfile) mpisize
#else
    write(chkptfile) 1
#endif
    write(chkptfile) nl
    write(chkptfile) nu
    write(chkptfile) tstep
    write(chkptfile) energy0
    write(chkptfile) wavefn
    write(chkptfile) wavefn0
    write(chkptfile) espectrum
    close(chkptfile)

  end subroutine s_checkpoint



  subroutine s_checkrestart(nl, nu)

    integer(idp), dimension(dim), intent(inout) :: nl, nu

    character(len=numcharlen) :: numchar
    integer :: tmp
    integer(idp), dimension(dim) :: tmp1, tmp2


#ifdef MPI
    numchar = f_numchar(mpirank)
#else
    numchar = f_numchar(0)
#endif

    open(unit=chkptfile, file=trim(parameters%chkpntfile)//numchar, form="unformatted", action="read")

    read(chkptfile) tmp

#ifdef MPI
    if (tmp /= mpisize) then
#else
    if (tmp /= 1) then
#endif
       
       call s_display(text="**** Error, trying to restart with incompatible number of processors")

#ifdef MPI
       call s_mpiabort()
#endif

       stop

    end if
    

    read(chkptfile) tmp1
    read(chkptfile) tmp2

    if (tmp1(1) /= nl(1) .or. tmp1(2) /= nl(2) .or. tmp2(1) /= nu(1) .or. tmp2(2) /= nu(2)) then
       
       call s_display(text="**** Error, restart size incompatible with input file")

#ifdef MPI
       call s_mpiabort()
#endif

       stop

    end if

    close(chkptfile)

  end subroutine s_checkrestart



  subroutine s_restart(tstep, wavefn, wavefn0, nl, nu, energy0, espectrum)

    integer(idp), intent(out) :: tstep
    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(out) :: wavefn, wavefn0
    complex(cdp), dimension(:), intent(out) :: espectrum
    real(fdp), intent(out) :: energy0

    character(len=numcharlen) :: numchar
    character(len=maxstringlength) :: potline
    integer :: readerr, tmp
    integer(idp), dimension(dim) :: tmp1

#ifdef MPI
    numchar = f_numchar(mpirank)
#else
    numchar = f_numchar(0)
#endif

    open(unit=chkptfile, file=trim(parameters%chkpntfile)//numchar, form="unformatted", action="read")
    read(chkptfile) tmp
    read(chkptfile) tmp1
    read(chkptfile) tmp1
    read(chkptfile) tstep
    read(chkptfile) energy0
    read(chkptfile) wavefn
    read(chkptfile) wavefn0
    read(chkptfile) espectrum
    close(chkptfile)

  end subroutine s_restart



  subroutine s_write2drarray(array, nl, nu, filename, basename, skip,&
       & incboundary)
   
    integer(idp), dimension(dim), intent(in) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in),&
         & target :: array
    character(len = *), intent(in) :: filename
    character(len = *), intent(in), optional :: basename
    integer(idp), dimension(dim), intent(in), optional :: skip
    logical, intent(in), optional :: incboundary

    integer(idp), dimension(dim) :: skipt
    logical :: incboundaryt
    integer(idp), dimension(dim) :: n1, nut, nlt
    real(fdp), dimension(:, :), allocatable :: writearray
    character(len = maxfilenamelen) :: fullfilename
#ifdef MPI
    integer :: handle
    integer(idp), dimension(dim) :: effectiven
#endif


    if (present(skip)) then

       skipt = skip

    else

       skipt = max(floor(dble(parameters%n)/dble(parameters&
            &%plotpointnum)), 1)

    end if

    if (present(incboundary)) then

       incboundaryt = incboundary

    else

       incboundaryt = parameters%plotcomplex

    end if

    if (incboundaryt) then

       n1 = 0

    else

       n1 = parameters%extraboundary*parameters%smooth / 2

    end if

    if (present(basename)) then
       
       fullfilename = trim(parameters%datadir)//trim(basename)//trim(filename)

    else

       fullfilename = trim(parameters%datadir)//trim(parameters%basename)&
            &//trim(filename)
       
    end if

#ifdef MPI
    if (mpirank == mpistart) then

       nlt = (/ nl(1)+n1(1), nl(2)+n1(2) /)
       nut = (/ nu(1)-n1(1), nu(2) /)

    elseif (mpirank == mpiend) then

       nlt = (/ nl(1)+n1(1), nl(2)+mod(nl(2)+n1(2)+1, skipt(2)) /)
       nut = (/ nu(1)-n1(1), nu(2)-n1(2) /)

    else

       nlt = (/ nl(1)+n1(1), nl(2)+mod(nl(2)+n1(2)+1, skipt(2)) /)
       nut = (/ nu(1)-n1(1), nu(2) /)

    end if

    allocate(writearray((nut(1)-nlt(1))/skipt(1)+1, (nut(2)-nlt(2))&
         &/skipt(2)+1))

    writearray = array(nlt(1):nut(1):skipt(1), nlt(2):nut(2):skipt(2))

    effectiven = (/ size(writearray, 1), size(writearray, 2) /)
    
    call s_mpiallreducesumint(effectiven(2))
    
    call s_mpifileopenwrite(fullfilename, handle)
    
    if (mpirank == mpistart) then
       
       call s_mpifilewriteint(handle, (/ size(effectiven)*idp /), 1)
       call s_mpifilewriteint(handle, effectiven(1), size(effectiven))
       call s_mpifilewriteint(handle, (/ size(effectiven)*idp /), 1)
       call s_mpifilewriteint(handle, (/ product(effectiven)*fdp /),&
            & 1)
       
    end if
    
    call s_mpifilewriteorderedreal(handle, writearray(1, 1),&
         & size(writearray))
    
    if (mpirank == mpistart) then
       
       call s_mpifilewriteint(handle, (/ product(effectiven)*fdp /),&
            & 1)

    end if
       
    call s_mpifileclose(handle)

    deallocate(writearray)
#else
    nlt = (/ nl(1)+n1(1), nl(2)+n1(2) /)
    nut = (/ nu(1)-n1(1), nu(2)-n1(2) /)

    allocate(writearray((nut(1)-nlt(1))/skipt(1)+1, (nut(2)-nlt(1))&
         &/skipt(2)+1))

    writearray = array(nlt(1):nut(1):skipt(1), nlt(2):nut(2):skipt(2))

    open(unit = 1, file = fullfilename, action = 'write', form&
         &='unformatted')

    write(1) (/ size(writearray, 1), size(writearray, 2) /)
    
    write(1) writearray
    
    close(1)

    deallocate(writearray)
#endif
       
  end subroutine s_write2drarray



  subroutine s_write2dcarray(array, nl, nu, filename, basename, skip,&
       & incboundary)
   
    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in),&
         & target :: array
    character(len = *), intent(in) :: filename
    character(len = *), intent(in), optional :: basename
    integer(idp), dimension(dim), intent(in), optional :: skip
    logical, intent(in), optional :: incboundary

    integer(idp), dimension(dim) :: skipt
    logical :: incboundaryt
    integer(idp), dimension(dim) :: n1, nut, nlt
    complex(cdp), dimension(:, :), allocatable :: writearray
    character(len = maxfilenamelen) :: fullfilename
#ifdef MPI
    integer :: handle
    integer(idp), dimension(dim) :: effectiven
#endif


    if (present(skip)) then

       skipt = skip

    else

       skipt = max(floor(dble(parameters%n)/dble(parameters&
            &%plotpointnum)), 1)

    end if

    if (present(incboundary)) then

       incboundaryt = incboundary

    else

       incboundaryt = parameters%plotcomplex

    end if

    if (incboundaryt) then

       n1 = 0

    else

       n1 = parameters%extraboundary*parameters%smooth / 2

    end if

    if (present(basename)) then
       
       fullfilename = trim(parameters%datadir)//trim(basename)//trim(filename)

    else

       fullfilename = trim(parameters%datadir)//trim(parameters%basename)&
            &//trim(filename)
       
    end if

#ifdef MPI
    if (mpirank == mpistart) then

       nlt = (/ nl(1)+n1(1), nl(2)+n1(2) /)
       nut = (/ nu(1)-n1(1), nu(2) /)

    elseif (mpirank == mpiend) then

       nlt = (/ nl(1)+n1(1), nl(2)+mod(nl(2)+n1(2)+1, skipt(2)) /)
       nut = (/ nu(1)-n1(1), nu(2)-n1(2) /)

    else

       nlt = (/ nl(1)+n1(1), nl(2)+mod(nl(2)+n1(2)+1, skipt(2)) /)
       nut = (/ nu(1)-n1(1), nu(2) /)

    end if

    allocate(writearray((nut(1)-nlt(1))/skipt(1)+1, (nut(2)-nlt(1))&
         &/skipt(2)+1))

    writearray = array(nlt(1):nut(1):skipt(1), nlt(2):nut(2):skipt(2))

    call s_mpifileopenwrite(fullfilename, handle)
    
    if (mpirank == mpistart) then
       
       effectiven = (/ size(writearray, 1), size(writearray, 2) /)
       
       call s_mpiallreducesumint(effectiven(2))

       call s_mpifilewriteint(handle, (/ size(effectiven)*idp /), 1)
       call s_mpifilewriteint(handle, effectiven(1), size(effectiven))
       call s_mpifilewriteint(handle, (/ size(effectiven)*idp /), 1)
       call s_mpifilewriteint(handle, (/ product(effectiven)*cdp*2 /)&
            &, 1)

    end if
       
    call s_mpifilewriteorderedcomp(handle, writearray(1, 1),&
         & size(writearray))

    if (mpirank == mpistart) then

       call s_mpifilewriteint(handle, (/ product(effectiven)*cdp*2 /)&
            &, 1)

    end if
       
    call s_mpifileclose(handle)

    deallocate(writearray)
#else
    nlt = (/ nl(1)+n1(1), nl(2)+n1(2) /)
    nut = (/ nu(1)-n1(1), nu(2)-n1(2) /)

    allocate(writearray((nut(1)-nlt(1))/skipt(1)+1, (nut(2)-nlt(1))&
         &/skipt(2)+1))

    writearray = array(nlt(1):nut(1):skipt(1), nlt(2):nut(2):skipt(2))

    open(unit = 1, file = fullfilename, action = 'write', form&
         &='unformatted')

    write(1) (/ size(writearray, 1), size(writearray, 2) /)
    
    write(1) writearray
    
    close(1)

    deallocate(writearray)
#endif
       
  end subroutine s_write2dcarray



  subroutine s_write1darray(array, filename, basename, skip1,&
       & incboundary)

    real(fdp), dimension(:), intent(in), target :: array
    character(len = *), intent(in) :: filename
    character(len = *), intent(in), optional :: basename
    integer(idp), intent(in), optional :: skip1
    logical, intent(in), optional :: incboundary
    real(fdp), dimension(:), pointer :: writearray

    integer(idp) :: skip1t
    logical :: incboundaryt
    integer(idp) :: n, n1


#ifdef MPI
    if (mpirank == mpistart) then
#endif

       if (present(skip1)) then
          
          skip1t = skip1
          
       else
          
          skip1t = 1
          
       end if
       
       if (present(incboundary)) then
          
          incboundaryt = incboundary
          
       else
          
          incboundaryt = .true.
          
       end if
       
       n = size(array)
       
       if (incboundaryt) then
          
          n1 = n
          
       else
          
          n1 = size(array(parameters%extraboundary(1)/2*parameters &
               &%smooth(1)+1:n - parameters%extraboundary(1)/2 &
               &*parameters%smooth(1)))
          
       end if
       
       if (present(basename)) then
          
          open(unit = 1, file = trim(parameters%datadir)//trim(basename) &
               &//trim(filename), action = 'write', form&
               &='unformatted')
          
       else
          
          open(unit = 1, file = trim(parameters%datadir)//trim(parameters &
               &%basename)//trim(filename), action = 'write', form &
               &='unformatted')
          
       end if
       
       writearray => array((n - n1)/2+1:(n + n1)/2:skip1t)
       
       write(1) size(writearray)
       
       write(1) writearray
       
       close(1)

#ifdef MPI
    endif
#endif

  end subroutine s_write1darray



  subroutine s_display(text, realnum, intnum, realnums, intnums, char)
    
    character(len=*), intent(in) :: text
    real(fdp), intent(in), optional :: realnum
    integer(idp), intent(in), optional :: intnum
    real(fdp), dimension(:), intent(in), optional :: realnums
    integer(idp), dimension(:), intent(in), optional :: intnums
    character(len=*), intent(in), optional :: char

    integer(idp) :: n, i
    character(len=maxstringlength) :: t_string
    

#ifdef MPI
    if (mpirank == mpistart) then
#endif
       
       if (present(realnum)) then
          
          write(*, '(a30, es18.10e3)') text, realnum
          
       else if (present(intnum)) then
          
          write(*, '(a30, i16)') text, intnum
          
       else if (present(char)) then
          
          write(*, '(a30, a)') text, char
          
       else if (present(realnums)) then
          
          n = size(realnums)
          
          t_string = ''
          
          do i=1, n
             
             t_string = trim(t_string) // ', es18.10e3'
             
          end do
          
          t_string = '(a30'//trim(t_string)//')'
          
          write(*, t_string) text, realnums
          
       else if (present(intnums)) then
          
          n = size(intnums)
          
          t_string = ''
          
          do i=1, n
             
             t_string = trim(t_string) // ', i16'
             
          end do
          
          t_string = '(a30'//trim(t_string)//')'
          
          write(*, t_string) text, intnums
          
       else
          
          write(*, '(a)') text
          
       end if

#ifdef MPI
    endif
#endif
    
  end subroutine s_display



  subroutine s_readline(string)
    
    character(len=*), intent(out) :: string

    
    string = ''

    read(inputfile, '(a)') string

    string = trim(adjustl(string))

  end subroutine s_readline



  subroutine s_getatom(string, atom)

    character(len=*), intent(inout) :: string
    character(len=*), intent(out) :: atom

    integer(idp) :: n

    string = trim(adjustl(string))
    atom = ''
    
    n = scan(string, ' ')

    if (n == 0) then
       
       read(string, '(a)') atom
       
    else
       
       read(string(:n), '(a)') atom
       
    end if

    atom = trim(adjustl(atom))

    string = trim(adjustl(string(n:)))

  end subroutine s_getatom
  



  subroutine s_read(object, string, realnum, intnum, logic, realnums,&
       & intnums, logics)
    
    character(len=*), intent(out), optional :: object, string
    real(fdp), intent(out), optional :: realnum
    integer(idp), intent(out), optional :: intnum
    logical, intent(out), optional :: logic
    real(fdp), dimension(:), intent(out), optional :: realnums
    integer(idp), dimension(:), intent(out), optional :: intnums
    logical, dimension(:), intent(out), optional :: logics

    character(len=maxstringlength) :: inputline, atom
    real(fdp) :: t_real
    integer(idp) :: n, i

    
    if (present(object)) then
       
       call s_readline(object)
       
    end if
    
    if (present(string)) then
       
       call s_readline(string)
       
    end if
    
    if (present(realnum)) then
       
       call s_readline(inputline)
       
       realnum = f_evaluatefn(inputline, constants, constantvalues)
       
    end if
    
    if (present(intnum)) then
       
       call s_readline(inputline)
       
       t_real = f_evaluatefn(inputline, constants, constantvalues)
       
       intnum = aint(t_real)

    end if

    if (present(logic)) then

       call s_readline(inputline)

       read(inputline, *) logic

    end if

    if (present(realnums)) then

       n = size(realnums)

       call s_readline(inputline)

       do i = 1, n

          call s_getatom(inputline, atom)

          realnums(i) = f_evaluatefn(atom, constants, constantvalues)

       end do

    end if

    if (present(intnums)) then

       n = size(intnums)

       call s_readline(inputline)
       
       do i = 1,n
          
          call s_getatom(inputline, atom)

          t_real = f_evaluatefn(atom, constants, constantvalues)
          
          intnums(i) = aint(t_real)

       end do

    end if

    if (present(logics)) then
       
       n = size(logics)

       call s_readline(inputline)
       
       do i = 1, n

          call s_getatom(inputline, atom)
          
          read(atom, *) logics(i)
          
       end do
       
    end if

  end subroutine s_read
  
  
  
  function f_norm(array, nl, nu, gridspacing) result(area)
    
    !
    ! Compute the norm
    !

    integer(idp), dimension(dim) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & array
    real(fdp), dimension(dim), optional, intent(in) :: gridspacing
    real(fdp), dimension(dim) :: t_gridspacing
    real(fdp) :: area


    if (present(gridspacing)) then

       t_gridspacing = gridspacing

    else

       t_gridspacing = parameters%gridspacing

    end if

    area = sum(array)*product(t_gridspacing)

#ifdef MPI
    call s_mpiallreducesumreal(area)
#endif

  end function f_norm



  function f_arraymax(array, nl, nu) result(arraymax)
    
    !
    ! Compute the maximum value
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & array
    real(fdp) :: arraymax


    arraymax = maxval(array)

#ifdef MPI
    call s_mpiallreducemaxreal(arraymax)
#endif

  end function f_arraymax



  function f_arraymin(array, nl, nu) result(arraymin)
    
    !
    ! Compute the maximum value
    !

    integer(idp), dimension(dim), intent(in) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & array
    real(fdp) :: arraymin


    arraymin = minval(array)

#ifdef MPI
    call s_mpiallreduceminreal(arraymin)
#endif

  end function f_arraymin



  subroutine s_cshift(array, nl, nu, n1, n2)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: array
    integer(idp), intent(in) :: n1, n2
#ifdef MPI
    integer(idp) :: rotation, offset
    integer, dimension(2) :: sendto, recvfrom
    complex(cdp), dimension(:, :), allocatable :: buffer
    integer, dimension(4) :: mpirequests
    integer, dimension(MPI_STATUS_SIZE, 4) :: mpistatuss
#endif


#ifdef MPI
    mpirequests = MPI_REQUEST_NULL
    mpistatuss=0

    array = cshift(array, n1, 1)

    rotation = n2*mpisize/parameters%n(2)

    sendto(1) = modulo(mpirank+rotation, mpisize)
    sendto(2) = modulo(sendto(1)+1, mpisize)
    recvfrom(2) = modulo(mpirank-rotation, mpisize)
    recvfrom(1) = modulo(recvfrom(2)-1, mpisize)

    offset = modulo(recvfrom(2)*parameters%n(2)/mpisize + n2,&
         & parameters%n(2)) - nl(2) + 1

    allocate(buffer(nl(1):nu(1), nl(2):nu(2)))

    if (offset < 1) then

       call s_mpirecvcomp(buffer(nl(1), nl(2)), size(buffer(:,&
            & nl(2):nu(2))), recvfrom(2), mpirequests(1))
       call s_mpisendcomp(array(nl(1), nl(2)), size(array(:,&
            & nl(2):nu(2))), sendto(1), mpirequests(2))
       
    elseif (offset > nu(2)-nl(2)) then
       
       call s_mpirecvcomp(buffer(nl(1), nl(2)), size(buffer(:,&
            & nl(2):nu(2))), recvfrom(1), mpirequests(1))
       call s_mpisendcomp(array(nl(1), nl(2)), size(array(:,&
            & nl(2):nu(2))), sendto(2), mpirequests(2))
       
    else
       
       call s_mpirecvcomp(buffer(nl(1), nl(2)), size(buffer(:,&
            & nl(2):nl(2)+offset-1)), recvfrom(1), mpirequests(1))
       call s_mpirecvcomp(buffer(nl(1), nl(2)+offset), size(buffer(:,&
            & nl(2)+offset:nu(2))), recvfrom(2), mpirequests(2))
       call s_mpisendcomp(array(nl(1), nl(2)), size(array(:,&
            & nl(2):nu(2)-offset)), sendto(1), mpirequests(3))
       call s_mpisendcomp(array(nl(1), nu(2)-offset+1), size(array(:,&
            & nu(2)-offset+1:nu(2))), sendto(2), mpirequests(4))
       
    end if
    
    call s_mpiwaitall(size(mpirequests), mpirequests(1), mpistatuss(1&
         &, 1))
    
    array = buffer

    deallocate(buffer)
#else
    array = cshift(cshift(array, n1, 1), n2, 2)
#endif

  end subroutine s_cshift



  function f_transmission(array, nl, nu, potbound)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in),&
         & target :: array
    integer(idp), dimension(4), intent(in) :: potbound

    real(fdp) :: f_transmission


    f_transmission = 0.0d0

#ifdef MPI
    if (max(potbound(2), nl(1)) <= nu(1)) then
       
       f_transmission = sum(array(max(potbound(2), nl(1)):, :))&
            &*product(parameters%gridspacing)

    end if

    call s_mpiallreducesumreal(f_transmission)
#else
    f_transmission = sum(array(max(potbound(2), nl(1)):, :))&
         &*product(parameters%gridspacing)
#endif
    
  end function f_transmission



  function f_reflection(array, nl, nu, potbound)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in),&
         & target :: array
    integer(idp), dimension(4), intent(in) :: potbound
    
    real(fdp) :: f_reflection


    f_reflection = 0.0d0

#ifdef MPI
    if (min(potbound(1), nu(1)) >= nl(1)) then
       
       f_reflection = sum(array(:min(potbound(1), nu(1)), :))&
            &*product(parameters%gridspacing)
       
    end if

    call s_mpiallreducesumreal(f_reflection)
#else
    f_reflection = sum(array(:min(potbound(1), nu(1)), :))&
         &*product(parameters%gridspacing)
#endif

  end function f_reflection



  function f_trapped(array, nl, nu, potbound)
        
    integer(idp), dimension(dim), intent(in) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in),&
         & target :: array
    integer(idp), dimension(4), intent(in) :: potbound

    real(fdp) :: f_trapped


    f_trapped = 0.0d0

#ifdef MPI
    if (max(potbound(1), nl(1)) <= nu(1) .and. min(potbound(2),&
         & nu(1)) >= nl(1)) then
       
       f_trapped = sum(array(max(potbound(1), nl(1)):min(potbound(2),&
            & nu(1)), :))*product(parameters%gridspacing)
       
    end if

    call s_mpiallreducesumreal(f_trapped)
#else
    f_trapped = sum(array(max(potbound(1), nl(1)):min(potbound(2),&
         & nu(1)), :))*product(parameters%gridspacing)
#endif

  end function f_trapped



  function f_ptransmission(array, nl, nu)

    integer(idp), dimension(dim) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)) :: array

    real(fdp) :: f_ptransmission


    f_ptransmission = 0.0d0

#ifdef MPI
    if (max(parameters%n(1)/2, nl(1)) <= nu(1)) then
       
       f_ptransmission = sum(array(max(parameters%n(1)/2, nl(1)):,&
            & :))*product(parameters%pspacing)

    end if

    call s_mpiallreducesumreal(f_ptransmission)
#else
    f_ptransmission = sum(array(max(parameters%n(1)/2, nl(1)):, :))&
         &*product(parameters%pspacing)
#endif

  end function f_ptransmission



  function f_preflection(array, nl, nu)

    integer(idp), dimension(dim) :: nl, nu
    real(fdp), dimension(nl(1):nu(1), nl(2):nu(2)), target :: array

    real(fdp) :: f_preflection


    f_preflection = 0.0d0

#ifdef MPI
    if (min(parameters%n(1)/2, nu(1)) >= nl(1)) then
       
       f_preflection = sum(array(:min(parameters%n(1)/2, nu(1)), :))&
            &*product(parameters%pspacing)
       
    end if

    call s_mpiallreducesumreal(f_preflection)
#else
    f_preflection = sum(array(:min(parameters%n(1)/2, nu(1)), :))&
         &*product(parameters%pspacing)
#endif

  end function f_preflection



  subroutine s_momentumfac(factor, nl, nu, derdims, order)

    integer(idp), dimension(dim), intent(in) :: nu, nl
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(out) ::&
         & factor
    real(fdp), dimension(dim), intent(in) :: derdims
    integer(idp), dimension(dim), intent(in) :: order

    integer(idp) :: i, j


    do j = nl(2), nu(2)
       do i = nl(1), nu(1)
          
          factor(i, j) = derdims(1)*(ii*f_momentumx(i)/hbar)&
               &**order(1) + derdims(2)*(ii*f_momentumy(j)/hbar)&
               &**order(2)
       
       end do
    end do

    !
    ! Rotate the factor so that it aligns with the fourier transformed
    ! wave function (zero momentum state first)
    !

    call s_cshift(factor, nl, nu, parameters%n(1)/2, parameters%n(2)&
         &/2)

#ifdef MPI
    call s_mpitranspose(factor, nl, nu)
#endif

  end subroutine s_momentumfac



  function f_factorprime(n) result(n_t)

    !
    ! The fourier transform routine FFTW has been hand optimised for
    ! FFT's with factors given in optimalsize.  This routine factors
    ! an integer into these values for optimal performance.
    !

    integer(idp), parameter :: numoptimal = 5
    integer(idp), dimension(numoptimal), parameter :: optimalsize = (&
         &/ 2, 3, 5, 7, 11 /)

    integer(idp), intent(in) :: n
    integer(idp) :: n_t

    integer(idp) :: tmp, i
    logical :: cont


    n_t = n

    n_t = ceiling(dble(n_t) / 2.0d0 )
    
#ifdef MPI
    n_t = ceiling(dble(n_t) / dble(mpisize))
#endif
    
    tmp = n_t
    
    do while (tmp /= 1)
       
       tmp = n_t
       
       cont = .false.
       
       do while ((.not. cont ) .and. tmp /= 1)
          
          cont = .true.
          
          do i=1, numoptimal
             
             if (mod(tmp, optimalsize(i)) == 0) then
                
                tmp = tmp/optimalsize(i)
                
                cont = .false.
                
             end if
             
          end do

       end do

       if (tmp /= 1) then

          n_t = n_t + 1

       end if

    end do

    n_t = n_t * 2

#ifdef MPI
    n_t = n_t * mpisize
#endif

  end function f_factorprime



  subroutine s_blur(array, nl, nu, smooth)

    integer(idp), dimension(dim), intent(in) :: nu, nl
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: array
    real(fdp), dimension(dim), optional, intent(in) :: smooth

    real(fdp), dimension(dim) :: smoothness
    complex(cdp), dimension(:, :), allocatable :: blur


    allocate(blur(nl(1):nu(1), nl(2):nu(2)))

    if (present(smooth)) then
       
       smoothness = smooth
       
    else

       smoothness = dble(parameters%smooth)*parameters%gridspacing

    end if

    call s_gaussian(blur, nl, nu, smoothness, (/ 0.0d0, 0.0d0 /), (/&
         & 0.0d0, 0.0d0 /))

    blur = real(blur) / f_norm(real(blur), nl, nu, (/ 1.0d0, 1.0d0 /))

    !
    ! Perform the gaussian blur
    !

    if (f_arraymax(abs(array), nl, nu) > 0.0d0) then

       call s_c2dfft(array, nl, nu, fftfor)

       call s_c2dfft(blur, nl, nu, fftfor)

       array = array*blur
       
       call s_c2dfft(array, nl, nu, fftbac)

       array = array/dble(product(parameters%n))

       call s_cshift(array, nl, nu, parameters%n(1)/2, parameters&
            &%n(2)/2)

    else

       array = 0.0d0

    end if

    deallocate(blur)

  end subroutine s_blur



  subroutine s_smooth(array, nl, nu, smooth)

    !
    ! Smooth the array using an averaging technique
    !

    integer(idp), dimension(dim), intent(in) :: nu, nl
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: array
    integer(idp), intent(in) :: smooth

    complex(cdp), dimension(:, :), allocatable :: t_array
    integer(idp) :: i, j, ii, jj, k, nn
    real(fdp) :: s
#ifdef MPI
    integer, dimension(4) :: mpirequests
    integer, dimension(MPI_STATUS_SIZE, 4) :: mpistatuss
#endif


    nn = 2

    allocate(t_array(nl(1)-nn:nu(1)+nn, nl(2)-nn:nu(2)+nn))
    
    t_array = 0.0d0

    do k = 1, smooth
       
       t_array(nl(1):nu(1), nl(2):nu(2)) = array

#ifdef MPI
       mpirequests = MPI_REQUEST_NULL
       mpistatuss = 0

       if (mpirank == mpistart) then

          call s_mpirecvcomp(t_array(nl(1)-nn, nu(2)+1),&
               & size(t_array(:, nu(2)+1:nu(2)+nn)), mpirank+1,&
               & mpirequests(1))

          call s_mpisendcomp(t_array(nl(1)-nn, nu(2)-nn+1),&
               & size(t_array(:, nu(2)-nn+1:nu(2))), mpirank+1,&
               & mpirequests(2))
          
       elseif (mpirank == mpiend) then

          call s_mpirecvcomp(t_array(nl(1)-nn, nl(2)-nn),&
               & size(t_array(:, nl(2)-nn:nl(2)-1)), mpirank-1,&
               & mpirequests(1))

          call s_mpisendcomp(t_array(nl(1)-nn, nl(2)), size(t_array(:&
               &, nl(2):nl(2)+nn-1)), mpirank-1, mpirequests(2))
          
       else

          call s_mpirecvcomp(t_array(nl(1)-nn, nu(2)+1),&
               & size(t_array(:, nu(2)+1:nu(2)+nn)), mpirank+1,&
               & mpirequests(1))
          call s_mpirecvcomp(t_array(nl(1)-nn, nl(2)-nn),&
               & size(t_array(:, nl(2)-nn:nl(2)-1)), mpirank-1,&
               & mpirequests(2))

          call s_mpisendcomp(t_array(nl(1)-nn, nu(2)-nn+1),&
               & size(t_array(:, nu(2)-nn+1:nu(2))), mpirank+1,&
               & mpirequests(3))
          call s_mpisendcomp(t_array(nl(1)-nn, nl(2)), size(t_array(:&
               &, nl(2):nl(2)+nn-1)), mpirank-1, mpirequests(4))
          
       end if

       call s_mpiwaitall(size(mpirequests), mpirequests(1),&
            & mpistatuss(1, 1))
#endif

       do j = max(nl(2), nn+1), min(nu(2), parameters%n(2)-nn)
          
          do i = max(nl(1), nn+1), min(nu(1), parameters%n(1)-nn)
             
             s = 0.0d0
             
             do jj = -nn, nn
                
                do ii = -nn, nn
                   
                   s = s + t_array(i+ii, j+jj)/max(1.0d0, 8.0d0&
                        &*abs(dble(ii)), 8.0d0*abs(dble(jj)))
                   
                end do
                
             end do

             array(i, j) = s/dble(nn+1)
             
          end do
          
       end do

    end do

    deallocate(t_array)

  end subroutine s_smooth



  subroutine s_smoothedge(array, nl, nu, smooth)

    !
    ! Produce a nice smooth array using both the smooth and
    ! the gauss blur
    !

    integer(idp), dimension(dim), intent(in) :: nu, nl
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: array
    integer(idp), dimension(dim), intent(in) :: smooth

    complex(cdp), dimension(:, :), allocatable, save :: smoothedge
    integer(idp) :: i, j
    integer(idp), dimension(dim), save :: nuold=0, nlold=0


    if ((nlold(1) /= nl(1) .or. nuold(1) /= nu(1) .or. nlold(2) /=&
         & nl(2) .or. nuold(2) /= nu(2)) .and. allocated(smoothedge))&
         & then

       deallocate(smoothedge)

    end if

    if (.not. allocated(smoothedge)) then
       
       allocate(smoothedge(nl(1):nu(1), nl(2):nu(2)))
       
       smoothedge = 0.0d0

       do j = max(nl(2), smooth(2)+1), min(nu(2), nu(2)-smooth(2))
          
          do i = max(nl(1), smooth(1)+1), min(nu(1), nu(1)-smooth(1))
             
             smoothedge(i, j) = 1.0d0
             
          end do
          
       end do
       
       call s_smooth(smoothedge, nl, nu, 4)
       
       smoothedge = real(smoothedge)
       
       call s_blur(smoothedge, nl, nu, dble(smooth*2)*parameters&
            &%gridspacing)
       
       smoothedge = real(smoothedge)
       
       smoothedge = (smoothedge - minval(real(smoothedge)))&
            &/(maxval(real(smoothedge))-minval(real(smoothedge)))
       
       !    call s_write2drarray(real(smoothedge), nl, nu, filename
       !='smoothedge', incboundary=.true.)

    end if
       
    array = array * smoothedge

  end subroutine s_smoothedge



  subroutine s_gaussian(array, nl, nu, width, position, p,&
       & gridspacing)

    integer(idp), dimension(dim), intent(in) :: nu, nl
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(out) ::&
         & array
    real(fdp), dimension(dim), optional, intent(in) :: width,&
         & position, p, gridspacing
    
    real(fdp), dimension(dim) :: t_width, t_position, t_p,&
         & t_gridspacing
    real(fdp), dimension(dim) :: halfwidth, amplitude, halfwidth_t
    integer(idp) :: i, j
    real(fdp) :: x, y


    if (present(width)) then

       t_width = width

    else

       t_width = parameters%wavefnwidth

    end if

    if (present(position)) then

       t_position = position

    else

       t_position = parameters%position

    end if

    if (present(p)) then

       t_p = p

    else

       t_p = parameters%p

    end if

    if (present(gridspacing)) then

       t_gridspacing = gridspacing

    else

       t_gridspacing = parameters%gridspacing

    end if

    halfwidth_t = 0.0d0
    halfwidth = t_width
    amplitude = 1.0d0/2.0d0/pi/halfwidth(1)/halfwidth(2)

    do while (abs(halfwidth_t(1)-halfwidth(1))/minval(parameters&
         &%spacesize) > cutoff .or. abs(halfwidth_t(2)-halfwidth(2))&
         &/minval(parameters%spacesize) > cutoff)

       halfwidth_t = halfwidth

       halfwidth = sqrt(abs(-t_width**2/2.0d0/log(cutoff/amplitude)))

       amplitude = 1.0d0/2.0d0/pi/halfwidth(1)/halfwidth(2)

    end do

    array = 0.0d0

    do j = nl(2), nu(2)

       y = f_positiony(j) - t_position(2)

       do i = nl(1), nu(1)

          x = f_positionx(i) - t_position(1)

          array(i, j) = exp(-x**2/2.0d0/halfwidth(1)**2-y**2/2.0d0&
               &/halfwidth(2)**2+ii*t_p(1)*x/hbar+ii*t_p(2)*y/hbar)

       end do

    end do

    array = array/sqrt(f_norm(real(array*conjg(array)), nl, nu,&
         & t_gridspacing))

  end subroutine s_gaussian



!!$  subroutine s_relaxpoisson(pot, mask, gridspacing)
!!$
!!$    complex(cdp), dimension(:, :), intent(inout) :: pot
!!$    logical, dimension(:, :), intent(in) :: mask
!!$    real(fdp), dimension(dim), intent(in) :: gridspacing
!!$
!!$    real(fdp) :: norm, oldnorm, relaxation
!!$    integer(idp) :: i, j
!!$
!!$
!!$    oldnorm = 0.0d0
!!$
!!$    norm = f_norm(real(pot*conjg(pot)), gridspacing)
!!$
!!$    !
!!$    ! Solves poisson's equation by successive over relaxation
!!$    !
!!$
!!$    relaxation = 2.0d0/(1.0d0+pi/sqrt(sum(gridspacing**2)))
!!$
!!$    do while (abs(norm - oldnorm) > sqrt(cutoff)*parameters !!$   
  !!      &%potrange(2))
!!$
!!$       do j = f_gridpointy(parameters%potbound(3)), !!$           
  !! & f_gridpointy(parameters%potbound(4))
!!$
!!$          do i = f_gridpointx(parameters%potbound(1)), !!$        
  !!       & f_gridpointx(parameters%potbound(2))
!!$
!!$             if (.not. mask(i, j)) then
!!$
!!$                pot(i, j) = pot(i, j) + relaxation*(gridspacing(1)
  !!**2 !!$                     &/2.0d0/sum(gridspacing**2)*(pot(i+1,
  !! j)+pot(i-1 !!$                     &, j)) + gridspacing(2)**2
  !!/2.0d0/sum(gridspacing !!$                     &**2)*(pot(i, j+1)
  !!+pot(i, j-1)) - pot(i, j))
!!$
!!$             end if
!!$
!!$          end do
!!$
!!$       end do
!!$
!!$       do j = f_gridpointy(parameters%potbound(3)), !!$           
  !! & f_gridpointy(parameters%potbound(4))
!!$
!!$          do i = f_gridpointx(parameters%potbound(1)), !!$        
  !!       & f_gridpointx(parameters%potbound(2))
!!$
!!$             if (.not. mask(i, j)) then
!!$
!!$                pot(i, j) = pot(i, j) + relaxation*(gridspacing(1)
  !!**2 !!$                     &/2.0d0/sum(gridspacing**2)*(pot(i+1,
  !! j)+pot(i-1 !!$                     &, j)) + gridspacing(2)**2
  !!/2.0d0/sum(gridspacing !!$                     &**2)*(pot(i, j+1)
  !!+pot(i, j-1)) - pot(i, j))
!!$
!!$             end if
!!$
!!$          end do
!!$
!!$       end do
!!$
!!$       oldnorm = norm
!!$
!!$       norm = f_norm(real(pot*conjg(pot)), gridspacing)
!!$
!!$    end do
!!$
!!$  end subroutine s_relaxpoisson



  function f_fermifunction(energy) result(amp)

    real(fdp), intent(in) :: energy
    real(fdp) :: amp, e


    e = (energy - parameters%fermienergy)/boltzk/parameters &
         &%temperature

    amp = 1.0d0/(exp(min(e, 200.0d0)) + 1.0d0)
    
  end function f_fermifunction



  subroutine s_fermidistribution(array, nl, nu, position, p,&
       & gridspacing)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(inout)&
         & :: array
    real(fdp), dimension(dim), optional, intent(in) :: position, p,&
         & gridspacing

    real(fdp), dimension(dim) :: t_position, t_p, t_gridspacing
    real(fdp) :: energy
    integer(idp), dimension(dim) :: gridposition
    integer(idp) :: i, j


    if (present(position)) then

       t_position = position

    else

       t_position = parameters%position

    end if

    if (present(p)) then

       t_p = p

    else

       t_p = parameters%p

    end if

    if (present(gridspacing)) then

       t_gridspacing = gridspacing

    else

       t_gridspacing = parameters%gridspacing

    end if

    do j = nl(2), nu(2)

       do i = nl(1), nu(1)

          energy = sum((f_momentum((/ i, j /))-t_p)**2)/2.0d0&
               &/parameters%mass

          array(i, j) = f_fermifunction(energy)

       end do

    end do

    call s_cshift(array, nl, nu, -parameters%n(1)/2, -parameters%n(2)&
         &/2)

#ifdef MPI
    call s_mpitranspose(array, nl, nu)
#endif

    call s_c2dfft(array, nl, nu, fftbac)

    call s_cshift(array, nl, nu, parameters%n(1)/2, parameters%n(2)/2)

    gridposition = f_gridpoint(t_position)

    call s_cshift(array, nl, nu, parameters%n(1)/2 - gridposition(1),&
         & parameters%n(2)/2-gridposition(2))

    array = array/sqrt(f_norm(real(array*conjg(array)), nl, nu,&
         & t_gridspacing))

  end subroutine s_fermidistribution


#ifndef MPI
  subroutine s_tspectrum(wavefn, wavefn0, nl, nu)

    integer(idp), dimension(dim), intent(in) :: nu, nl
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & wavefn, wavefn0

    complex(cdp), dimension(:, :), allocatable :: wavefn1_t, wavefn2_t
    real(fdp) :: norm1_t
    integer(idp) :: i


    allocate(wavefn1_t(nl(1):nu(1), nl(2):nu(2)),&
         & wavefn2_t(nl(1):nu(1), nl(2):nu(2)))

    wavefn1_t = wavefn0
    wavefn2_t = wavefn
    
    call s_c2dfft(wavefn1_t, nl, nu, fftfor)
    call s_c2dfft(wavefn2_t, nl, nu, fftfor)
    
    wavefn1_t = real(wavefn1_t*conjg(wavefn1_t))
    wavefn2_t = real(wavefn2_t*conjg(wavefn2_t))

    norm1_t = maxval(abs(wavefn1_t))
    
    where(real(wavefn1_t) > parameters%speccutoff*norm1_t)
       
       wavefn2_t = wavefn2_t / wavefn1_t
       
    elsewhere
       
       wavefn2_t = 0.0d0
       
    end where
    
    call s_cshift(wavefn2_t, nl, nu, parameters%n(1)/2, parameters&
         &%n(2)/2)

    call s_write2drarray(real(wavefn2_t), nl, nu, filename='t.spec'&
         &,skip=(/ 1, 1 /), incboundary=.true.)

    deallocate(wavefn1_t, wavefn2_t)
    
  end subroutine s_tspectrum



  subroutine s_pspectrum(wavefn, wavefn0, nl, nu)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & wavefn, wavefn0

    complex(cdp), dimension(:, :), allocatable :: wavefn1_t, wavefn2_t
    real(fdp), dimension(:), allocatable :: pspectrum1, pspectrum2
    integer(idp) :: i


    allocate(wavefn1_t(nl(1):nu(1), nl(2):nu(2)),&
         & wavefn2_t(nl(1):nu(1), nl(2):nu(2)),&
         & pspectrum1(nl(1):nu(1)), pspectrum2(nl(1):nu(1)))

    wavefn1_t = wavefn0
    wavefn2_t = wavefn
    
    call s_c2dfft(wavefn1_t, nl, nu, fftfor)
    call s_c2dfft(wavefn2_t, nl, nu, fftfor)
    
    wavefn1_t = abs(real(wavefn1_t*conjg(wavefn1_t)))
    wavefn2_t = abs(real(wavefn2_t*conjg(wavefn2_t)))

    do i = nl(1), nu(1)
       
       pspectrum2(i) = sum(real(wavefn2_t(i, :)))
       pspectrum1(i) = sum(real(wavefn1_t(i, :)))
       
    end do
    
    where(real(pspectrum1) > parameters%speccutoff*maxval(pspectrum1))
       
       pspectrum2 = pspectrum2/pspectrum1
       
    elsewhere
       
       pspectrum2 = 0.0d0
       
    end where
    
    pspectrum2 = cshift(pspectrum2, parameters%n(1)/2, 1)
    
    call s_write1darray(pspectrum2, 'p.spec')

    deallocate(wavefn1_t, wavefn2_t, pspectrum1, pspectrum2)
    
  end subroutine s_pspectrum



  subroutine s_espectrum(espectrum)

    complex(cdp), dimension(:), intent(in) :: espectrum

    complex(cdp), dimension(:), allocatable :: espectrum_t
    

    allocate(espectrum_t(size(espectrum)))

    espectrum_t = espectrum
    
    call s_write1darray(sqrt(real(espectrum_t*conjg(espectrum_t)))&
         &,'e')
    call s_write1darray(real(espectrum_t),'e.r')
    call s_write1darray(aimag(espectrum_t),'e.c')

    call s_c1dfft(espectrum_t, parameters%timesteps, fftbac)
    
    call s_write1darray(sqrt(real(espectrum_t*conjg(espectrum_t)))&
         &,'e.spec')
    call s_write1darray(real(espectrum_t),'e.r.spec')
    call s_write1darray(aimag(espectrum_t),'e.c.spec')

    deallocate(espectrum_t)

  end subroutine s_espectrum



  subroutine s_conductance(wavefn, nl, nu)

    integer(idp), dimension(dim), intent(in) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), intent(in) ::&
         & wavefn

    complex(cdp), dimension(:, :), allocatable :: wavefn_t, cond
    integer(idp) :: i, j
    complex(cdp) :: num1, num2


    allocate(cond(parameters%n(1)/2, parameters%n(2)),&
         & wavefn_t(nl(1):nu(1), nl(2):nu(2)))

    wavefn_t = wavefn

    call s_c2dfft(wavefn_t, nl, nu, fftfor)

    wavefn_t = cshift(cshift(wavefn_t, parameters%n(1)/2, 1),&
         & parameters%n(2)/2, 2)

    cond = 0.0d0

    do j=1, parameters%n(2)
       
       do i=1, parameters%n(1)/2
        
          num1 = wavefn_t(i, j)
          num2 = wavefn_t(parameters%n(1)/2+1-i, j)
          
          if (real(num1*conjg(num1)) > sqrt(cutoff)) then
             
             cond(parameters%n(1)/2+1-i, j) = num2/num1
             
          end if
          
       end do

    end do
    
    call s_write2drarray(sqrt(real(cond*conjg(cond))), nl, nu,&
         & filename ='cond', incboundary=.true.)

    deallocate(cond, wavefn_t)

  end subroutine s_conductance
#endif


  subroutine s_arraymatch(array1, array2, nl, nu)

    !
    ! Adjust array1 to as close as possible to array2
    ! Using Newton's Method
    !

    integer(idp), dimension(dim) :: nl, nu
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), target,&
         & intent(inout) :: array1
    complex(cdp), dimension(nl(1):nu(1), nl(2):nu(2)), target,&
         & intent(in) :: array2

    real(fdp) :: c, d, dd, cmin, dmin
    integer(idp) :: count, i, j


    c = 1.0d0

    d = 0.0d0

    do j = max(nl(2), 1+parameters%smooth(2)), min(nu(2), parameters&
         &%n(2)-parameters%smooth(2))

       do i = max(nl(1), 1+parameters%smooth(1)), min(nu(1),&
            & parameters%n(1)-parameters%smooth(1))

          if (abs(real(array2(i, j))) > cutoff**2) then
             
             d = d + sqrt((c*real(array1(i, j))/real(array2(i, j)) -&
                  & 1.0d0)**2)

          end if

       end do

    end do

#ifdef MPI
    call s_mpiallreducesumreal(d)
#endif

    count = 1

    dmin = d
    
    cmin = c

    do while (abs(d) > cutoff .and. count < 20)
       
       count = count + 1
       
       dd = 0.0d0
       
       do j = max(nl(2), 1+parameters%smooth(2)), min(nu(2),&
            & parameters%n(2)-parameters%smooth(2))
          
          do i = max(nl(1), 1+parameters%smooth(1)), min(nu(1),&
               & parameters%n(1)-parameters%smooth(1))
             
             if (abs(real(array2(i, j))) > cutoff**2) then
                
                dd = dd + real(array1(i, j))/real(array2(i,j))*(c &
                     &*real(array1(i, j))/real(array2(i, j))-1.0d0) &
                     &/sqrt((c*real(array1(i, j))/real(array2(i, j)) &
                     &-1.0d0)**2)
                
             end if
             
          end do
          
       end do

#ifdef MPI
       call s_mpiallreducesumreal(dd)
#endif
       
       c = c - d/dd
       
       d = 0.0d0
       
       do j = max(nl(2), 1+parameters%smooth(2)), min(nu(2),&
            & parameters%n(2)-parameters%smooth(2))
          
          do i = max(nl(1), 1+parameters%smooth(1)), min(nu(1),&
               & parameters%n(1)-parameters%smooth(1))
             
             if (abs(real(array2(i, j))) > cutoff**2) then
                
                d = d + sqrt((c*real(array1(i, j))/real(array2(i, j))&
                     & - 1.0d0)**2)
                
             end if
             
          end do
          
       end do

#ifdef MPI
       call s_mpiallreducesumreal(d)
#endif
       
       if (abs(d) < abs(dmin)) then
          
          dmin = d
          cmin = c
          
       end if
       
    end do
    
    array1 = cmin*array1
    
  end subroutine s_arraymatch
  


end module misc
