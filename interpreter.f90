module interpreter

  !
  ! This module implements a function interpreter so that function
  ! form potentials can be used
  !

  use precision

  implicit none

  integer(idp), private, parameter :: commandlength = 10
  integer(idp), private, parameter :: variablelength = 10
  integer(idp), private, parameter :: errorlength = 255
  character(len = variablelength), private, dimension(:), allocatable :: variables
  real(fdp), private, dimension(:), allocatable :: varvalues
  character(len=errorlength), private :: error='OK'


contains



  subroutine s_seterror(t_error)

    character(len=*), intent(in) :: t_error


    if (len_trim(t_error) > 0) then

       error = trim(t_error)

    end if

  end subroutine s_seterror



  subroutine s_assignvariables(vars, varvals)

    character(len = *), intent(in) :: vars
    real(fdp), dimension(:), intent(in) :: varvals
    
    integer(idp) :: n


    n = size(varvals)

    allocate(variables(n), varvalues(n))

    read(vars, *) variables

    varvalues = varvals

  end subroutine s_assignvariables



  subroutine s_erasevariables()

    if (allocated(variables)) then

       deallocate(variables)

    end if

    if (allocated(varvalues)) then

       deallocate(varvalues)

    end if

  end subroutine s_erasevariables



  function f_readnumber(func) result(number)

    character(len = *), intent(inout) :: func
    real(fdp) :: number

    integer(idp) :: pos
    logical :: dflag, exponent, finish


    exponent = .false.
    dflag = .false.
    finish = .false.
    pos = 0
    number = 0.0d0

    do while(pos < len_trim(func) .and. .not. finish)

       pos = pos + 1

       select case(func(pos:pos))

       case ('0':'9')
          if (dflag) then

             exponent = .true.
             dflag = .false.

          end if

       case ('d','e')
          if (dflag .or. exponent) then

             finish = .true.

          else

             dflag = .true.

          end if

       case('.')
          if (dflag .or. exponent) then

             finish = .true.

          end if

       case('+', '-')
          if (dflag) then
             
             dflag = .false.
             exponent = .true.

          else

             finish = .true.

          end if

       case default
          finish = .true.

       end select

    end do

    if (pos == len_trim(func) .and. .not. finish) then

       read(func, *) number

       func = ''

    else

       read(func(:pos-1), *) number

       func = func(pos:)

    end if

  end function f_readnumber



  subroutine s_readtext(func, command, number)

    character(len = *), intent(inout) :: func
    character(len = commandlength), intent(out) :: command
    real(fdp), intent(out) :: number

    integer(idp) :: i, funccut
    integer(idp), parameter :: numcommands = 16
    character(len = commandlength), dimension(numcommands) :: commands


    commands(1) = 'cos'
    commands(2) = 'sin'
    commands(3) = 'tan'
    commands(4) = 'exp'
    commands(5) = 'log'
    commands(6) = 'log10'
    commands(7) = 'sqrt'
    commands(8) = 'acos'
    commands(9) = 'asin'
    commands(10) = 'atan'
    commands(11) = 'cosh'
    commands(12) = 'sinh'
    commands(13) = 'tanh'
    commands(14) = 'anint'
    commands(15) = 'aint'
    commands(16) = 'abs'

    command = ''
    number = 0.0d0
    funccut = 0

    do i = 1, numcommands

       if (index(trim(func), trim(commands(i))//'(') == 1 .and.&
            & len_trim(commands(i)) > funccut) then

          command = commands(i)
          
          number = 0.0d0

          funccut = len_trim(command)

       end if
       
    end do

    if (command == '') then

       do i = 1, size(varvalues)
          
          if (index(trim(func), trim(variables(i))) == 1 .and.&
               & len_trim(variables(i)) > funccut) then
             
             number = varvalues(i)

             command = 'variable'

             funccut = len_trim(variables(i))

          end if
          
       end do
       
    end if

    if (command == '') then

       call s_seterror('Unknown command or variable: '//trim(func))

    end if
    
    if (len_trim(func) > funccut) then

       func = func(funccut+1:)

    else

       func = ''

    end if

    if (funccut == 0) then

       func = func(2:)

    end if

  end subroutine s_readtext



  function f_readoperator(func) result(operator)

    character(len = *), intent(inout) :: func
    character(len = commandlength) :: operator

    integer(idp) :: i
    integer(idp), parameter :: numoperators = 6
    character(len = 2), dimension(numoperators) :: operators

    
    operators(1) = '+'
    operators(2) = '-'
    operators(3) = '/'
    operators(4) = '*'
    operators(5) = '**'
    operators(6) = '^'

    operator = ''

    do i = 1, numoperators

       if (index(trim(func), trim(operators(i))) == 1 .and.&
            & len_trim(operators(i)) > len_trim(operator)) then

          operator = operators(i)

       end if
       
    end do

    if (operator == '') then

       call s_seterror('Unknown operator: '//operator)

    end if
    
    if (len_trim(func) > len_trim(operator)) then

       func = func(len_trim(operator)+1:)

    else

       func = ''

    end if

  end function f_readoperator



  function f_readbracket(func) result(subfunc)

    character(len = *), intent(inout) :: func
    character(len = len(func)) :: subfunc

    integer(idp) :: pos, count

    
    pos = 0
    count = 0

    do while((pos < len_trim(func)) .and. (count >= 0))

       pos = pos + 1

       select case(func(pos:pos))

       case('(')
          count = count + 1

       case(')')
          count = count - 1

       case default

       end select

       if (count == 0) then

          exit

       end if

    end do

    subfunc = trim(func(2:pos-1))

    if (count /= 0) then

       call s_seterror('Unbalanced brackets: '//trim(func))

    end if

    if (len_trim(subfunc) == 0) then

       call s_seterror('Empty brackets encountered: '//trim(func))

    end if

    if (len_trim(func) < pos+1) then

       func = ''

    else
       
       func = func(pos+1:)

    end if

  end function f_readbracket



  recursive function f_compute(func, value, command) result(number)

    character(len = *), intent(inout) :: func
    real(fdp), intent(in) :: value
    character(len = commandlength), intent(in) :: command
    real(fdp) :: number

    character(len = commandlength) :: nextcommand
    character(len = len(func)) :: func_t
    real(fdp) :: nextnumber, real_t


    select case(command)

    case('number')
       number = value

    case('variable')
       number = value

    case('bracket')
       func_t = f_readbracket(func)
       number = f_evaluate(func_t)

    case('nop')
       number = value

    case('end')
       number = value

    case('+')
       number = value + f_evaluateblock(func, command)

    case('-')
       number = value - f_evaluateblock(func, command)

    case('*')
       number = value * f_evaluateblock(func, command)

    case('/')
       number = value / f_evaluateblock(func, command)

    case('**', '^')
       number = value ** f_evaluateblock(func, command)

    case('cos')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = cos(f_compute(func, nextnumber, nextcommand))

    case('sin')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = sin(f_compute(func, nextnumber, nextcommand))

    case('tan')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = tan(f_compute(func, nextnumber, nextcommand))

    case('exp')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = exp(f_compute(func, nextnumber, nextcommand))

    case('log')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = log(f_compute(func, nextnumber, nextcommand))

    case('log10')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = log10(f_compute(func, nextnumber, nextcommand))

    case('sqrt')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = sqrt(f_compute(func, nextnumber, nextcommand))

    case('acos')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = acos(f_compute(func, nextnumber, nextcommand))

    case('asin')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = asin(f_compute(func, nextnumber, nextcommand))

    case('atan')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = atan(f_compute(func, nextnumber, nextcommand))

    case('cosh')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = cosh(f_compute(func, nextnumber, nextcommand))

    case('sinh')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = sinh(f_compute(func, nextnumber, nextcommand))

    case('tanh')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = tanh(f_compute(func, nextnumber, nextcommand))

    case('anint')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = anint(f_compute(func, nextnumber, nextcommand))

    case('aint')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = aint(f_compute(func, nextnumber, nextcommand))

    case('abs')
       call s_getnextatom(func, nextnumber, nextcommand)

       number = abs(f_compute(func, nextnumber, nextcommand))

    case default
       call s_seterror('Unknown input: '//command)

    end select

  end function f_compute



  subroutine s_getnextatom(func, number, command)

    character(len = *), intent(inout) :: func
    real(fdp), intent(inout) :: number
    character(len = commandlength), intent(inout) :: command


    command = ''

    do while (command == '' .or. command == 'nop')

       command = ''
       
       select case(func(1:1))
          
       case('0':'9', '.')
          number = f_readnumber(func)
          command = 'number'
          
       case('+', '-', '/', '*', '^')
          command = f_readoperator(func)
          
       case('a':'z','A':'Z')
          call s_readtext(func, command, number)
          
       case('(')
          command = 'bracket'

       case(' ')
          if (len_trim(func) > 1) then
             
             command = 'nop'
             
             func = func(2:)

          else

             func = ''

             command = 'end'

          end if

       case default
          if (len_trim(func) > 1) then
             
             call s_seterror('Unknown input: '//trim(func))
             
             func = func(2:)
             
             command = 'nop'

          else
             
             func = ''
             
             command = 'end'
             
          end if
          
       end select

    end do

  end subroutine s_getnextatom



  function f_order(command) result(order)

    character(len = *), intent(in) :: command
    real(fdp) :: order


    order = 0.0d0
    
    select case(trim(command))

    case('+', '-')
       order = 1.0d0

    case('*', '/')
       order = 2.0d0

    case('**', '^')
       order = 3.0d0

    case default
       order = 0.0d0

    end select
       
  end function f_order
  
  
    
  subroutine s_peeknextatom(func, number, command)

    character(len = *), intent(in) :: func
    real(fdp), intent(inout) :: number
    character(len = commandlength), intent(inout) :: command

    character(len = len(func)) :: func_t


    func_t = func
    
    call s_getnextatom(func_t, number, command)
    
  end subroutine s_peeknextatom


  recursive function f_evaluateblock(func, blockcommand)&
       & result(number)

    character(len = *), intent(inout) :: func
    character(len = *), intent(in) :: blockcommand
    real(fdp) :: number

    character(len = commandlength) :: command, nextcommand
    real(fdp) :: nextnumber


    call s_peeknextatom(func, nextnumber, nextcommand)

    if (trim(nextcommand) == 'end') then

       call s_seterror('Premature end of function')

    end if

    do while((f_order(nextcommand) == 0 .or. (f_order(nextcommand) >&
         & f_order(blockcommand))) .and. trim(nextcommand) /= 'end')
       
       call s_getnextatom(func, number, command)
          
       number = f_compute(func, number, command)
       
       call s_peeknextatom(func, nextnumber, nextcommand)

    end do

  end function f_evaluateblock



  recursive function f_evaluate(func) result(number)

    character(len = *), intent(inout) :: func
    real(fdp) :: number

    character(len = commandlength) :: command
    character(len = len(func)) :: func_t


    number = 0
    command = ''

    do while(command /= 'end')
       
       call s_getnextatom(func, number, command)
       
       number = f_compute(func, number, command)

    end do

  end function f_evaluate



  function f_evaluatefn(func, vars, varvals, t_error) result(number)

    character(len = *), intent(in) :: func, vars
    real(fdp), dimension(:), intent(in) :: varvals
    character(len=*), intent(out), optional :: t_error
    real(fdp) :: number

    character(len = len(func)) :: func_t
    integer(idp) :: i

    
    call s_assignvariables(vars, varvals)

    func_t = trim(func)

    number = f_evaluate(func_t)

    if (present(t_error)) then

       t_error = error

    end if

    call s_erasevariables()

  end function f_evaluatefn



end module interpreter
