!*******************************************************************************
!> license: BSD
!
!  Modern Fortran function parser.
!
!### Original Description
!  This function parser module is intended for applications where a set of mathematical
!  fortran-style expressions is specified at runtime and is then evaluated for a large
!  number of variable values. This is done by compiling the set of function strings
!  into byte code, which is interpreted efficiently for the various variable values.
!
!### History
!  * Original code (Fortran 90 function parser v1.1) by:
!    Roland Schmehl <roland.schmehl@alumni.uni-karlsruhe.de>.
!    The source code is available from http://fparser.sourceforge.net
!  * The function parser concept is based on a C++ class library written by Juha
!    Nieminen <warp@iki.fi> available from http://warp.povusers.org/FunctionParser/
!  * Expanded and updated to modern Fortran by Jacob Williams, 2017.
!
!### License
!  * Copyright (c) 2000-2008, Roland Schmehl. All rights reserved.
!  * Copyright (c) 2017, Jacob Williams. All rights reserved.
!  * This software is distributable under the BSD license. See the terms of the
!    BSD license in the documentation provided with this software.
!
!@note The default real kind (`wp`) can be
!      changed using optional preprocessor flags.
!      This library was built with real kind:
#ifdef REAL32
!      `real(kind=real32)` [4 bytes]
#elif REAL64
!      `real(kind=real64)` [8 bytes]
#elif REAL128
!      `real(kind=real128)` [16 bytes]
#else
!      `real(kind=real64)` [8 bytes]
#endif

    module function_parser

    use error_module,    only: list_of_errors
    use iso_fortran_env

    implicit none

    private

#ifdef REAL32
    integer,parameter,public :: fparser_rk = real32   !! real kind used by this module [4 bytes]
#elif REAL64
    integer,parameter,public :: fparser_rk = real64   !! real kind used by this module [8 bytes]
#elif REAL128
    integer,parameter,public :: fparser_rk = real128  !! real kind used by this module [16 bytes]
#else
    integer,parameter,public :: fparser_rk = real64   !! real kind used by this module [8 bytes]
#endif

    integer,parameter :: wp = fparser_rk  !! local copy of `fparser_rk` with a shorter name

    !parameters:
    real(wp), parameter :: zero = 0.0_wp
    real(wp), parameter :: one  = 1.0_wp
    real(wp), parameter :: pi = acos(-one)

    ! Note: these should be continuous, unique integers:
    ! [they must have the values that correspond to the array indices below]
    integer, parameter ::   cImmed   = 1
    integer, parameter ::   cNeg     = 2
    integer, parameter ::   cAdd     = 3, &    ! Operators
                            cSub     = 4, &
                            cMul     = 5, &
                            cDiv     = 6, &
                            cPow     = 7
    integer, parameter ::   cAbs     = 8,  &   ! Functions
                            cExp     = 9,  &
                            cLog10   = 10, &
                            cLog     = 11, &
                            cSqrt    = 12, &
                            cSinh    = 13, &
                            cCosh    = 14, &
                            cTanh    = 15, &
                            cSin     = 16, &
                            cCos     = 17, &
                            cTan     = 18, &
                            cAsin    = 19, &
                            cAcos    = 20, &
                            cAtan2   = 21, &    ! atan2 must precede atan to prevent aliasing.
                            cAtan    = 22, &
                            cPi      = 23, &    ! Pi (function with zero arguments)
                            cCeil    = 24, &
                            cFloor   = 25, &
                            cGamma   = 26, &
                            cHypot   = 27, &
                            cMax     = 28, &
                            cMin     = 29, &
                            cModulo  = 30, &
                            cMod     = 31, &
                            cSign    = 32, &
                            cIf      = 33     ! if (three arguments)
    integer, parameter ::   VarBegin = 34

    character(len=1), dimension(cAdd:cPow), parameter ::  operators = [ '+', &  ! plus
                                                                        '-', &  ! minus
                                                                        '*', &  ! multiply
                                                                        '/', &  ! divide
                                                                        '^'  ]  ! power

    character(len=7), dimension(cAbs:cIf), parameter :: functions = [    'abs    ', &
                                                                         'exp    ', &
                                                                         'log10  ', &
                                                                         'log    ', &
                                                                         'sqrt   ', &
                                                                         'sinh   ', &
                                                                         'cosh   ', &
                                                                         'tanh   ', &
                                                                         'sin    ', &
                                                                         'cos    ', &
                                                                         'tan    ', &
                                                                         'asin   ', &
                                                                         'acos   ', &
                                                                         'atan2  ', &
                                                                         'atan   ', &
                                                                         'pi     ', &
                                                                         'ceiling', &
                                                                         'floor  ', &
                                                                         'gamma  ', &
                                                                         'hypot  ', &
                                                                         'max    ', &
                                                                         'min    ', &
                                                                         'modulo ', &
                                                                         'mod    ', &
                                                                         'sign   ', &
                                                                         'if     ' ]

    ! Specify the number of required arguments each `functions` element must have.
    integer, dimension(cAbs:cIf), parameter :: required_args = [ 1, & ! abs
                                                                 1, & ! exp
                                                                 1, & ! log10
                                                                 1, & ! log
                                                                 1, & ! sqrt
                                                                 1, & ! sinh
                                                                 1, & ! cosh
                                                                 1, & ! tanh
                                                                 1, & ! sin
                                                                 1, & ! cos
                                                                 1, & ! tan
                                                                 1, & ! asin
                                                                 1, & ! acos
                                                                 2, & ! atan2
                                                                 1, & ! atan
                                                                 0, & ! pi
                                                                 1, & ! Ceiling
                                                                 1, & ! Floor
                                                                 1, & ! Gamma
                                                                 2, & ! Hypot
                                                                 2, & ! Max
                                                                 2, & ! Min
                                                                 2, & ! Modulo
                                                                 2, & ! Mod
                                                                 2, & ! Sign
                                                                 3  ] ! if

    ! Specify the number of optional arguments each `functions` element might have.
    integer, dimension(cAbs:cIf), parameter :: optional_args = [ 0, & ! abs
                                                                 0, & ! exp
                                                                 0, & ! log10
                                                                 0, & ! log
                                                                 0, & ! sqrt
                                                                 0, & ! sinh
                                                                 0, & ! cosh
                                                                 0, & ! tanh
                                                                 0, & ! sin
                                                                 0, & ! cos
                                                                 0, & ! tan
                                                                 0, & ! asin
                                                                 0, & ! acos
                                                                 0, & ! atan2
                                                                 1, & ! atan
                                                                 0, & ! pi
                                                                 0, & ! Ceiling
                                                                 0, & ! Floor
                                                                 0, & ! Gamma
                                                                 0, & ! Hypot
                                                                 0, & ! Max
                                                                 0, & ! Min
                                                                 0, & ! Modulo
                                                                 0, & ! Mod
                                                                 0, & ! Sign
                                                                 0  ] ! if

    ! The maximum number of arguments any `functions` element might have.
    integer, parameter :: max_func_args = maxval(required_args + optional_args)

    !the list of error messages:
    integer,parameter :: error_div_by_zero       = 1
    integer,parameter :: error_sqrt_arg_neg      = 2
    integer,parameter :: error_log_arg_neg       = 3
    integer,parameter :: error_asin_arg_illegal  = 4
    integer,parameter :: error_acos_arg_illegal  = 5
    integer,parameter :: error_invalid_operation = 6
    character(len=25),dimension(6),parameter :: error_messages = &
        [ 'Division by zero         ', &    ! 1
          'Argument of SQRT negative', &    ! 2
          'Argument of LOG negative ', &    ! 3
          'Argument of ASIN illegal ', &    ! 4
          'Argument of ACOS illegal ', &    ! 5
          'Invalid operation        '  ]

    type stack_func_container
        !! to create an array of the function pointers in the fparser
        procedure(stack_func),pointer,nopass :: f => null()
    end type stack_func_container

  !****************************************************************
  !>
  !  The function parser class.

    type,public :: fparser

        private

        type(stack_func_container),dimension(:),allocatable :: bytecode_ops  !! array of function pointers
        integer,dimension(:),allocatable :: bytecode  !! array of integers
        integer :: bytecodesize = 0

        real(wp),dimension(:),allocatable :: immed
        integer :: immedsize = 0

        real(wp),dimension(:),allocatable :: stack
        integer :: stacksize = 0
        integer :: stackptr = 0

        type(list_of_errors) :: error_msg       !! list of error messages

    contains

        private

        procedure,public :: parse    => parse_function
        procedure,public :: evaluate => evaluate_function
        procedure,public :: destroy  => destroy_parser
        procedure,public :: error
        procedure,public :: print_errors
        procedure,public :: clear_errors

        procedure :: compile_substr
        procedure :: compile
        procedure :: check_syntax
        procedure :: add_error

    end type fparser
    !****************************************************************

    !****************************************************************
    !>
    !  A wrapper to [[fparser]] to evaluate an array of functions.
    !
    !@note Each parser has the same variables.

      type,public :: fparser_array
          private
          type(fparser),dimension(:),allocatable :: f  !! An array of parsers.
                                                       !! Each one has the
                                                       !! same variables.
      contains
          private
          procedure,public :: parse        => parse_function_array
          procedure,public :: evaluate     => evaluate_function_array
          procedure,public :: destroy      => destroy_parser_array
          procedure,public :: error        => error_array
          procedure,public :: print_errors => print_errors_array
          procedure,public :: clear_errors => clear_errors_array
      end type fparser_array
      !****************************************************************

    !interface to functions:
    abstract interface
        subroutine stack_func(me,ip,dp,sp,val,ierr)
            !! a function that operates on the stack
            import :: wp,fparser
            implicit none
            class(fparser),intent(inout)     :: me
            integer,intent(in)               :: ip    !! instruction pointer
            integer,intent(inout)            :: dp    !! data pointer
            integer,intent(inout)            :: sp    !! stack pointer
            real(wp),dimension(:),intent(in) :: val   !! variable values
            integer,intent(out)              :: ierr  !! error flag
        end subroutine stack_func
    end interface

    contains
!*******************************************************************************

!*******************************************************************************
!>
!  Returns true if there are any errors in the class.

    pure elemental function error(me)

    implicit none

    class(fparser),intent(in) :: me
    logical :: error !! true if there are any errors in the class

    error = me%error_msg%has_errors()

    end function error
!*******************************************************************************

!*******************************************************************************
!>
!  Prints the error messages (if any) in the class.

    subroutine print_errors(me,iunit)

    implicit none

    class(fparser),intent(inout) :: me
    integer,intent(in) :: iunit  !! unit number for printing
                                 !! (assumed to be open)

    call me%error_msg%print(iunit)

    end subroutine print_errors
!*******************************************************************************

!*******************************************************************************
!>
!  Clears any error messages in the class.

    pure elemental subroutine clear_errors(me)

    implicit none

    class(fparser),intent(inout) :: me

    call me%error_msg%destroy()

    end subroutine clear_errors
!*******************************************************************************

!*******************************************************************************
!>
!  Returns true if there are any errors in the class.

    pure elemental function error_array(me)

    implicit none

    class(fparser_array),intent(in) :: me
    logical :: error_array !! true if there are any errors in the class

    if (allocated(me%f)) then
        error_array = any(me%f%error())
    else
        error_array = .false.
    end if

    end function error_array
!*******************************************************************************

!*******************************************************************************
!>
!  Prints the error messages (if any) in the class.

    subroutine print_errors_array(me,iunit)

    implicit none

    class(fparser_array),intent(inout) :: me
    integer,intent(in) :: iunit  !! unit number for printing
                                 !! (assumed to be open)

    integer :: i  !! counter

    if (allocated(me%f)) then
        do i=1, size(me%f)
            call me%f(i)%print_errors(iunit)
        end do
    end if

    end subroutine print_errors_array
!*******************************************************************************

!*******************************************************************************
!>
!  Clears any error messages in the class.

    pure elemental subroutine clear_errors_array(me)

    implicit none

    class(fparser_array),intent(inout) :: me

    integer :: i  !! counter

    if (allocated(me%f)) then
        do i=1, size(me%f)
            call me%f(i)%clear_errors()
        end do
    end if

    end subroutine clear_errors_array
!*******************************************************************************

!*******************************************************************************
!>
!  [[fparser]] destructor.
!
!  This can be called manually, and it is also called in [[parse_function]].

    pure elemental subroutine destroy_parser(me)

    implicit none

    class(fparser),intent(inout) :: me

    if (allocated(me%bytecode))     deallocate(me%bytecode)
    if (allocated(me%immed))        deallocate(me%immed)
    if (allocated(me%stack))        deallocate(me%stack)
    if (allocated(me%bytecode_ops)) deallocate(me%bytecode_ops)

    call me%error_msg%destroy()

    end subroutine destroy_parser
!*******************************************************************************

!*******************************************************************************
!>
!  [[fparser_array]] destructor.

    pure elemental subroutine destroy_parser_array(me)

    implicit none

    class(fparser_array),intent(inout) :: me

    integer :: i !! counter

    if (allocated(me%f)) then
        do i=1,size(me%f)
            call me%f(i)%destroy()
        end do
        deallocate(me%f)
    end if

    end subroutine destroy_parser_array
!*******************************************************************************

!*******************************************************************************
!>
!  Parse the function string `funcstr` and compile it into bytecode

    subroutine parse_function (me, funcstr, var, case_sensitive)

    implicit none

    class(fparser),intent(inout)               :: me
    character(len=*),intent(in)                :: funcstr         !! function string
    character(len=*), dimension(:), intent(in) :: var             !! array with variable names
    logical,intent(in),optional                :: case_sensitive  !! are the variables case sensitive?
                                                                  !! [default is false]

    character (len=:),allocatable                :: func     !! function string, local use
    character(len=len(var)),dimension(size(var)) :: tmp_var  !! variable list, local use
    integer,dimension(:),allocatable             :: ipos

    logical :: is_case_sensitive
    integer :: k

    if (present(case_sensitive)) then
        is_case_sensitive = case_sensitive
    else
        is_case_sensitive = .false.
    end if

    !first, initialize:
    call me%destroy()

    !if is case insensitive, then convert both to lowercase:
    func = trim(adjustl(funcstr))   ! local copy of function string
    if (is_case_sensitive) then
        tmp_var = var
    else
        call to_lowercase (func, func)    ! local copy of function string
        call to_lowercase (var, tmp_var)  !
    end if

    ! Allocate and initialize map from preprocessed function string to user input
    ! function string (used to print useful error messages).
    allocate (ipos(len_trim(funcstr)))
    ipos = [ (k,k=1,size(ipos)) ]

    !preprocess and check syntax:
    call replace_string ('**','^ ',func)      ! exponent into 1-char. format
    call remove_spaces (func,ipos)            ! condense function string
    call me%check_syntax(func,funcstr,tmp_var,ipos)

    ! Do not compile if `check_syntax` failed.
    if (.not. me%error()) then
        call me%compile(func,tmp_var)         ! compile into bytecode
    end if

    deallocate (ipos)

    end subroutine parse_function
!*******************************************************************************

!*******************************************************************************
!>
!  Evaluate bytecode of function for the values passed in array `val`.

    subroutine evaluate_function (me, val, res)

    implicit none

    class(fparser),intent(inout)        :: me
    real(wp), dimension(:), intent(in)  :: val  !! variable values
    real(wp),intent(out)                :: res  !! result

    integer :: ip   !! instruction pointer
    integer :: dp   !! data pointer
    integer :: sp   !! stack pointer
    integer :: ierr !! error flag

    !initialize:
    dp = 1
    sp = 0

    !do all the operations:
    do ip=1,me%bytecodesize
        call me%bytecode_ops(ip)%f(me,ip,dp,sp,val,ierr)
        if (ierr/=0) then
            call me%error_msg%add(trim(get_error_message_string(ierr)))
            res = zero
            return
        end if
    end do

    !the result:
    res = me%stack(1)

    end subroutine evaluate_function
!*******************************************************************************

!*******************************************************************************
!>
!  Alternate version of [[parse_function]] for the [[fparser_array]] class.

    subroutine parse_function_array (me, funcstr, var, case_sensitive)

    implicit none

    class(fparser_array),intent(inout)         :: me
    character(len=*),dimension(:),intent(in)   :: funcstr         !! function string array
    character(len=*),dimension(:),intent(in)   :: var             !! array with variable names
    logical,intent(in),optional                :: case_sensitive  !! are the variables case sensitive?
                                                                  !! [default is false]

    integer :: i       !! counter
    integer :: n_funcs !! number of functions in the class

    !first, initialize:
    call me%destroy()

    n_funcs = size(funcstr)
    allocate(me%f(n_funcs))

    do i=1,n_funcs
        call me%f(i)%parse(funcstr(i),var,case_sensitive)
        if (me%f(i)%error_msg%has_errors()) exit  ! stop if there are any errors
    end do

    end subroutine parse_function_array
!*******************************************************************************

!*******************************************************************************
!>
!  Alternate version of [[evaluate_function]] for the [[fparser_array]] class.

    subroutine evaluate_function_array (me, val, res)

    implicit none

    class(fparser_array),intent(inout)  :: me
    real(wp), dimension(:), intent(in)  :: val  !! variable values
    real(wp),dimension(:),intent(out)   :: res  !! result. Should be `size(me%f)`

    integer :: i       !! counter
    integer :: n_funcs !! number of functions in the class

    if (allocated(me%f)) then
        n_funcs = size(me%f)
        if (n_funcs == size(res)) then
            do i=1,n_funcs
                call me%f(i)%evaluate(val,res(i))
                if (me%f(i)%error_msg%has_errors()) exit  ! stop if there are any errors
            end do
        else
            call me%f(i)%error_msg%add('Error: the res vector is not the correct size.')
            res = zero
        end if
    else
        call me%f(i)%error_msg%add('Error: the fparser_array has not been initialized.')
        res = zero
    end if

    end subroutine evaluate_function_array
!*******************************************************************************

!******************************************************************
    subroutine cimmed_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    sp=sp+1
    me%stack(sp) = me%immed(dp)
    dp=dp+1
    ierr = 0

    end subroutine cimmed_func
!******************************************************************

!******************************************************************
!>
!  Negative function

    subroutine cneg_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = -me%stack(sp)
    ierr = 0

    end subroutine cneg_func
!******************************************************************

!******************************************************************
!>
!  Add function

    subroutine cadd_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = me%stack(sp-1) + me%stack(sp)
    sp=sp-1
    ierr = 0

    end subroutine cadd_func
!******************************************************************

!******************************************************************
!>
!  Subtract function

    subroutine csub_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = me%stack(sp-1) - me%stack(sp)
    sp=sp-1
    ierr = 0

    end subroutine csub_func
!******************************************************************

!******************************************************************
!>
!  Multiply function

    subroutine cmul_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = me%stack(sp-1) * me%stack(sp)
    sp=sp-1
    ierr = 0

    end subroutine cmul_func
!******************************************************************

!******************************************************************
!>
!  Division function

    subroutine cdiv_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if (me%stack(sp)==zero) then

        ierr = error_div_by_zero    !divide by zero error

    else

        me%stack(sp-1) = me%stack(sp-1)/me%stack(sp)
        sp=sp-1
        ierr = 0

    end if

    end subroutine cdiv_func
!******************************************************************

!******************************************************************
!>
!  Power

    subroutine cpow_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = me%stack(sp-1)**me%stack(sp)
    sp=sp-1
    ierr = 0

    end subroutine cpow_func
!******************************************************************

!******************************************************************
!>
!  Absolute value function

    subroutine cabs_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = abs(me%stack(sp))
    ierr = 0

    end subroutine cabs_func
!******************************************************************

!******************************************************************
!>
!  Exponential function

    subroutine cexp_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = exp(me%stack(sp))
    ierr = 0

    end subroutine cexp_func
!******************************************************************

!******************************************************************
!>
!  log10 function

    subroutine clog10_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if (me%stack(sp)<=zero) then
        ierr = error_log_arg_neg
    else
        me%stack(sp) = log10(me%stack(sp))
        ierr = 0
    end if

    end subroutine clog10_func
!******************************************************************

!******************************************************************
!>
!  log function

    subroutine clog_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if (me%stack(sp)<=zero) then
        ierr = error_log_arg_neg
    else
        me%stack(sp)=log(me%stack(sp))
        ierr = 0
    end if

    end subroutine clog_func
!******************************************************************

!******************************************************************
!>
!  square root function

    subroutine csqrt_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if (me%stack(sp)<zero) then
        ierr = error_sqrt_arg_neg
    else
        me%stack(sp) = sqrt(me%stack(sp))
        ierr = 0
    end if

    end subroutine csqrt_func
!******************************************************************

!******************************************************************
!>
!  sinh function

    subroutine csinh_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = sinh(me%stack(sp))
    ierr = 0

    end subroutine csinh_func
!******************************************************************

!******************************************************************
!>
!  cosh function

    subroutine ccosh_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = cosh(me%stack(sp))
    ierr = 0

    end subroutine ccosh_func
!******************************************************************

!******************************************************************
!>
!  tanh function

    subroutine ctanh_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = tanh(me%stack(sp))
    ierr = 0

    end subroutine ctanh_func
!******************************************************************

!******************************************************************
!>
!  sin function

    subroutine csin_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = sin(me%stack(sp))
    ierr = 0

    end subroutine csin_func
!******************************************************************

!******************************************************************
!>
!  cos function

    subroutine ccos_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = cos(me%stack(sp))
    ierr = 0

    end subroutine ccos_func
!******************************************************************

!******************************************************************
!>
!  tan function

    subroutine ctan_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = tan(me%stack(sp))
    ierr = 0

    end subroutine ctan_func
!******************************************************************

!******************************************************************
!>
!  asin function

    subroutine casin_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if ((me%stack(sp)<-one).or.(me%stack(sp)>one)) then
        ierr = error_asin_arg_illegal
    else
        me%stack(sp)=asin(me%stack(sp))
        ierr = 0
    end if

    end subroutine casin_func
!******************************************************************

!******************************************************************
!>
!  acos function

    subroutine cacos_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if ((me%stack(sp)<-one).or.(me%stack(sp)>one)) then
        ierr = error_acos_arg_illegal
    else
        me%stack(sp) = acos(me%stack(sp))
        ierr = 0
    end if

    end subroutine cacos_func
!******************************************************************

!******************************************************************
!>
!  atan function

    subroutine catan_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp)=atan(me%stack(sp))
    ierr = 0

    end subroutine catan_func
!******************************************************************

!******************************************************************
!>
!  atan2 function

    subroutine catan2_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = atan2(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine catan2_func
!******************************************************************

!******************************************************************
!>
!  Pi. A function with zero arguments.

    subroutine cPi_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    sp = sp + 1
    me%stack(sp) = pi
    ierr = 0

    end subroutine cPi_func
!******************************************************************

!******************************************************************
!>
!  ceiling function

    subroutine cceil_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = ceiling(me%stack(sp))
    ierr = 0

    end subroutine cceil_func
!******************************************************************

!******************************************************************
!>
!  floor function

    subroutine cfloor_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = floor(me%stack(sp))
    ierr = 0

    end subroutine cfloor_func
!******************************************************************

!******************************************************************
!>
!  gamma function

    subroutine cgamma_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp) = gamma(me%stack(sp))
    ierr = 0

    end subroutine cgamma_func
!******************************************************************

!******************************************************************
!>
!  hypot function

    subroutine chypot_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = hypot(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine chypot_func
!******************************************************************

!******************************************************************
!>
!  max function

    subroutine cmax_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = max(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine cmax_func
!******************************************************************

!******************************************************************
!>
!  min function

    subroutine cmin_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = min(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine cmin_func
!******************************************************************

!******************************************************************
!>
!  mod function

    subroutine cmod_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = mod(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine cmod_func
!******************************************************************

!******************************************************************
!>
!  modulo function

    subroutine cmodulo_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = modulo(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine cmodulo_func
!******************************************************************

!******************************************************************
!>
!  sign function

    subroutine csign_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    me%stack(sp-1) = sign(me%stack(sp-1), me%stack(sp))
    sp = sp - 1
    ierr = 0

    end subroutine csign_func
!******************************************************************

!******************************************************************
!>
!  If function with three arguments.
!
!  `If(expression, value is true, value if false)`
!
!  Where: 0 is false and /=0 is true.

    subroutine cif_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    if (me%stack(sp-2) /= zero) then ! true
        me%stack(sp-2) = me%stack(sp-1)
    else ! false
        me%stack(sp-2) = me%stack(sp)
    end if

    sp = sp - 2
    ierr = 0

    end subroutine cif_func
!******************************************************************

!******************************************************************
    subroutine cdefault_func(me,ip,dp,sp,val,ierr)

    implicit none

    class(fparser),intent(inout)     :: me    !! parser class
    integer,intent(in)               :: ip    !! instruction pointer
    integer,intent(inout)            :: dp    !! data pointer
    integer,intent(inout)            :: sp    !! stack pointer
    real(wp),dimension(:),intent(in) :: val   !! variable values
    integer,intent(out)              :: ierr  !! error flag

    sp=sp+1
    me%stack(sp) = val(me%bytecode(ip)-varbegin+1)
    ierr = 0

    end subroutine cdefault_func
!******************************************************************

!*******************************************************************************
!>
!  Find the end positions of each function argument.
!
!                12345678901234567890123456789012
!  For example, "funcWithThreeArgs(12.4,x,abs(x))" will result in an
!  'arg_pos' of [22, 24, 31, 0, 0].
!
! Errors returned:
!   1 : Function's closing parenthesis is missing.
!   2 : Function has more than `max_func_args` arguments.
!   3 : Function has an empty argument.

    subroutine find_arg_positions(paren_start, func, num_args, arg_pos, ierr, err_pos)

    implicit none

    integer,                           intent(in)    :: paren_start   !! The position of the function's opening parenthesis.
    character(len=*),                  intent(in)    :: func          !! Pre-processed function string (no spaces or '**').
    integer,                           intent(out)   :: num_args      !! The number of arguments found.
    integer, dimension(max_func_args), intent(out)   :: arg_pos       !! The end positions of each argument.
    integer, optional,                 intent(out)   :: ierr          !! Error code; non-zero if an error occurred.
    integer, optional,                 intent(out)   :: err_pos       !! The position in `func` of an error; 0 unless `ierr` /= 0.

    integer, parameter :: missing_closing_paren = 1, &
                          too_many_args         = 2, &
                          empty_arg             = 3

    character(len=1) :: cur_char                !! The current character in `func` being processed.

    integer          :: cur_pos,         &      !! The current position in `func` being processed.
                        func_len,        &      !! The length of `func`.
                        open_parens,     &      !! The number of open parentheses.
                        arg_len,         &      !! The length of an argument.
                        iarg                    !! Argument index.

    ! Initialize outputs.
    num_args = 1
    arg_pos = 0
    if (present(ierr)) ierr = 0
    if (present(err_pos)) err_pos = 0

    func_len = len_trim(func)
    open_parens = 1

    cur_pos = paren_start + 1
    func_len = len_trim(func)

    ! Step through the function string until we find the function's closing parenthesis.
    ! Every time we find a comma character at the same parentheses level as the function's
    ! opening parenthesis, increment the number of arguments and record the previous
    ! argument's last character.
    do while (open_parens > 0)
        if (cur_pos > func_len) then
            ! The function did not have a closing parenthesis.
            if (present(ierr)) ierr = missing_closing_paren
            if (present(err_pos)) err_pos = func_len
            return
        end if

        cur_char = func(cur_pos:cur_pos)

        if (cur_char == '(') then
            open_parens = open_parens + 1

        elseif (cur_char == ')') then
            ! For a "functionWithNoArgs()" set num_args = 0 and return.
            if (cur_pos == paren_start + 1) then
                num_args = 0
                return
            end if

            open_parens = open_parens - 1

            ! We have arrived at the function's closing parenthesis.
            if (open_parens == 0) arg_pos(num_args) = cur_pos - 1

        elseif (cur_char == ',') then
            if (open_parens == 1) then
                ! Check if we have exceeded the maximum allowed number of function
                ! arguments.  This is not checking that 'num_args' is necessarily
                ! valid for a specific function.
                if (num_args + 1 > max_func_args) then
                    if (present(ierr)) ierr = too_many_args
                    if (present(err_pos)) err_pos = cur_pos
                    return
                endif

                arg_pos(num_args) = cur_pos - 1
                num_args = num_args + 1
            end if
        end if

        cur_pos = cur_pos + 1
    end do

    ! Check for empty arguments.
    do iarg = 1, num_args
        if (iarg == 1) then
            arg_len = arg_pos(iarg) - paren_start
        else
            arg_len = arg_pos(iarg) - arg_pos(iarg - 1) - 1
        endif

        if (arg_len == 0) then
            if (present(ierr)) ierr = empty_arg
            if (present(err_pos)) err_pos = arg_pos(iarg)
            return
        end if
    end do

    end subroutine find_arg_positions
!*******************************************************************************

!*******************************************************************************
!>
!  Check syntax of function string.

    recursive subroutine check_syntax (me,func,funcstr,var,ipos)

    implicit none

    class(fparser),intent(inout)              :: me
    character(len=*),intent(in)               :: func       !! function string without spaces
    character(len=*),intent(in)               :: funcstr    !! original function string
    character(len=*), dimension(:),intent(in) :: var        !! array with variable names
    integer,dimension(:),intent(in)           :: ipos

    integer          :: n
    character(len=1) :: c
    real(wp)         :: r
    logical          :: err
    integer          :: parcnt        !! parenthesis counter
    integer          :: j,ib,in,lFunc
    integer          :: arg_pos(max_func_args)
    integer          :: num_args, iarg
    integer          :: ierr, err_pos
    integer          :: func_arg_ipos(size(ipos))
    integer          :: arg_start, arg_end
    logical          :: end_of_function

    lFunc = len_trim(func)

    ! Initial loop over function to verify matched parentheses.  This slightly
    ! increases runtime but allows for improved error reporting.
    parcnt = 0
    do j = 1, lFunc
         c = func(j:j)
         if (c == '(') then
            parcnt = parcnt + 1
         elseif (c == ')') then
            parcnt = parcnt - 1
         end if

        if (parcnt < 0) then
            call me%add_error(j, ipos, funcstr, 'Missing opening parenthesis')
            return
        end if
    end do

    if (parcnt /= 0) then
        call me%add_error(lFunc, ipos, funcstr, 'Missing closing parenthesis')
        return
    endif

    j = 1
    step: do
        if (j > lFunc) then
            call me%add_error(j, ipos, funcstr)
            return
        end if

        c = func(j:j)
        ! Check for valid operand (must appear)
        if (c == '-' .or. c == '+') then                      ! Check for leading - or +
            j = j+1
            if (j > lFunc) then
                call me%add_error(j, ipos, funcstr, 'Missing operand')
                return
            end if
            c = func(j:j)
            if (any(c == operators)) then
                call me%add_error(j, ipos, funcstr, 'Multiple operators')
                return
            end if
        end if

        end_of_function = .false.
        n = mathfunction_index (func(j:), var)
        if (n > 0) then                                       ! Check for math function
            j = j+len_trim(functions(n))
            if (j > lFunc) then
                call me%add_error(j, ipos, funcstr, 'Missing function argument')
                return
            end if
            c = func(j:j)
            if (c /= '(') then
                call me%add_error(j, ipos, funcstr, 'Missing opening parenthesis')
                return
            end if

            ! Find the number of function arguments and argument substring positions
            ! in `func`.
            call find_arg_positions(j, func, num_args, arg_pos, ierr, err_pos)
            if (ierr /= 0) then
                select case (ierr)
                    case     (1); call me%add_error(err_pos, ipos, funcstr, 'Missing function closing parenthesis')
                    case     (2); call me%add_error(err_pos, ipos, funcstr, 'Function has too many arguments')
                    case     (3); call me%add_error(err_pos, ipos, funcstr, 'Function has an empty argument')
                    case default; call me%add_error(err_pos, ipos, funcstr, 'Unknown find argument position error')
                end select
                return
            end if

            ! Verify that the number of function arguments present is consistent
            ! with the specified function.
            if (num_args < required_args(n)) then
                call me%add_error(j, ipos, funcstr, 'Missing required function argument')
                return
            elseif (num_args > required_args(n) + optional_args(n)) then
                call me%add_error(j, ipos, funcstr, 'Too many function arguments')
                return
            end if

            ! Recursively check each argument substring.
            if (num_args == 0) then
                j = j + 2
            else
                do iarg = 1, num_args
                    if (iarg == 1) then
                        arg_start = j + 1
                    else
                        arg_start = arg_pos(iarg-1) + 2
                    endif

                    arg_end = arg_pos(iarg)
                    func_arg_ipos(1:(arg_end-arg_start+1)) = ipos(arg_start:arg_end)

                    call me%check_syntax(func(arg_start:arg_end), funcstr, var, func_arg_ipos)
                    if (me%error()) return
                end do

                j = arg_pos(num_args) + 2
            endif

            if (j > lFunc) exit
            c = func(j:j)

            ! We've moved past the closing parenthesis of "someFunction(...)", so
            ! next we either need an operator or a closing parenthesis.
            end_of_function = .true.
        end if

        if (.not. end_of_function) then
            if (c == '(') then                                ! Check for opening parenthesis
                parcnt = parcnt+1
                j = j+1
                cycle step
            end if

            if (scan(c,'0123456789.') > 0) then               ! Check for number
                r = string_to_real (func(j:),ib,in,err)
                if (err) then
                    call me%add_error(j, ipos, funcstr, 'Invalid number format:  '//func(j+ib-1:j+in-2))
                    return
                end if
                j = j+in-1
                if (j > lFunc) exit
                c = func(j:j)
            else                                              ! Check for variable
                n = variable_index (func(j:),Var,ib,in)
                if (n == 0) then
                    call me%add_error(j, ipos, funcstr, 'Invalid element: '//func(j+ib-1:j+in-2))
                    return
                end if
                j = j+in-1
                if (j > lFunc) exit
                c = func(j:j)
            end if
        end if

        do while (c == ')')                                   ! Check for closing parenthesis
            parcnt = parcnt-1
            if (parcnt < 0) then
                call me%add_error(j, ipos, funcstr, 'Mismatched parenthesis')
                return
            end if
            if (func(j-1:j-1) == '(') then
                call me%add_error(j-1, ipos, funcstr, 'Empty parentheses')
                return
            end if
            j = j+1
            if (j > lFunc) exit
            c = func(j:j)
        end do

        ! Now, we have a legal operand: A legal operator or end of string must follow
        if (j > lFunc) exit
        if (any(c == operators)) then                         ! Check for multiple operators
            if (j+1 > lFunc) then
                call me%add_error(j, ipos, funcstr)
                return
            end if
            if (any(func(j+1:j+1) == operators)) then
                call me%add_error(j+1, ipos, funcstr, 'Multiple operators')
                return
            end if
        else                                                  ! Check for next operand
            call me%add_error(j, ipos, funcstr, 'Missing operator')
            return
        end if

        ! Now, we have an operand and an operator: the next loop will check for another
        ! operand (must appear)
        j = j+1
    end do step

    end subroutine check_syntax
!*******************************************************************************

!*******************************************************************************
!>
!  return error message string, given the error code.

    function get_error_message_string (ierr) result (msg)

    implicit none

    integer,intent(in)                   :: ierr   !! error message number
    character (len=len(error_messages))  :: msg    !! the error message string

    if (ierr == 0 .or. ierr > size(error_messages)) then
        msg = ''
    else
        msg = error_messages(ierr)
    endif

    end function get_error_message_string
!*******************************************************************************

!*******************************************************************************
!>
!  add error message to the class.

    subroutine add_error (me, j, ipos, funcstr, msg)

    implicit none

    class(fparser),intent(inout)         :: me
    integer,intent(in)                   :: j
    integer,dimension(:),intent(in)      :: ipos
    character(len=*),intent(in)          :: funcstr     !! original function string
    character(len=*),optional,intent(in) :: msg

    character(len=:),allocatable :: tmp !! to indicate where on
                                        !! the line the error occurs

    if (present(msg)) then
        call me%error_msg%add('*** Error in syntax of function string: '//Msg)
    else
        call me%error_msg%add('*** Error in syntax of function string:')
    endif

    call me%error_msg%add(' '//trim(funcstr))

    tmp = repeat(' ',ipos(min(j,size(ipos))))//'?'     ! Advance to the jth position
    call me%error_msg%add(tmp)
    deallocate(tmp)

    end subroutine add_error
!*******************************************************************************

!*******************************************************************************
!>
!  return operator index

    function operator_index (c) result (n)

    implicit none

    character(len=1), intent(in) :: c
    integer :: n

    integer :: j  !! counter

    n = 0
    do j=cadd,cpow
        if (c == operators(j)) then
            n = j
            exit
        end if
    end do

    end function operator_index
!*******************************************************************************

!*******************************************************************************
!>
!  Return index of math function beginning at 1st position of string `str`

    function mathfunction_index (str, var) result (n)

    implicit none

    character(len=*), intent(in) :: str
    character(len=*), dimension(:),intent(in) :: var        !! array with variable names
    integer :: n

    integer :: j
    integer :: k
    character (len=len(functions)) :: fun

    n = 0
    do j=cAbs,cIf                          ! check all math functions
       k = min(len_trim(functions(j)), len(str))
       call to_lowercase (str(1:k), fun)
       if (fun == functions(j)) then              ! compare lower case letters
          n = j                               ! found a matching function
          exit
       end if
    end do

    if (n>0) then
        if (any(functions(n) == var)) then
            ! in this case, there is a variable with the same
            ! name as this function. So, check to make sure this
            ! is really the function.
            if (k+1<=len(str)) then
                if (str(k+1:k+1) /= '(') then  ! this assumes that spaces have been removed
                    n = 0  ! assume it is the variable
                end if
            end if
        end if
    end if

    end function mathfunction_index
!*******************************************************************************

!*******************************************************************************
!>
!  return index of variable at begin of string str (returns 0 if no variable found)

    function variable_index (str, var, ibegin, inext) result (n)

    implicit none

    character(len=*),intent(in)              :: str     !! string
    character(len=*),dimension(:),intent(in) :: var     !! array with variable names
    integer                                  :: n       !! index of variable
    integer, optional,intent(out)            :: ibegin  !! start position of variable name
    integer, optional,intent(out)            :: inext   !! position of character after name

    integer :: j,ib,in,lstr

    n = 0
    ib = 0
    in = 0
    lstr = len_trim(str)
    if (lstr > 0) then
        do ib=1,lstr                                   ! search for first character in str
            if (str(ib:ib) /= ' ') exit                ! when lstr>0 at least 1 char in str
        end do
        do in=ib,lstr                                  ! search for name terminators
            if (scan(str(in:in),'+-*/^) ') > 0) exit   ! NOTE: all the operators must be here [cAdd,cSub,cMul,cDiv,cPow]
        end do
        do j=1,size(var)
            if (str(ib:in-1) == var(j)) then
                n = j                                  ! variable name found
                exit
            end if
        end do
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in

    end function variable_index
!*******************************************************************************

!*******************************************************************************
!>
!  Remove Spaces from string, remember positions of characters in old string

    subroutine remove_spaces (str, ipos)

    implicit none

    character(len=*),intent(inout)     :: str
    integer,dimension(:),intent(inout) :: ipos

    integer :: k,lstr

    lstr = len_trim(str)
    k = 1
    do while (str(k:lstr) /= ' ')
        if (str(k:k) == ' ') then
            str(k:lstr)  = str(k+1:lstr)//' '    ! move 1 character to left
            ipos(k:lstr) = [ ipos(k+1:lstr), 0 ] ! move 1 element to left
            k = k-1
        end if
        k = k+1
    end do

  end subroutine remove_spaces
!*******************************************************************************

!*******************************************************************************
!>
!  replace ALL appearances of character set `ca` in
!  string `str` by character set `cb`

    subroutine replace_string(ca,cb,str)

    implicit none

    character(len=*),intent(in)       :: ca
    character(len=len(ca)),intent(in) :: cb    !! `len(ca)` must be `len(cb)`
    character(len=*),intent(inout)    :: str

    integer :: j,lca

    lca = len(ca)
    do j=1,len_trim(str)-lca+1
        if (str(j:j+lca-1) == ca) str(j:j+lca-1) = cb
    end do

    end subroutine replace_string
!*******************************************************************************

!*******************************************************************************
!>
!  Compile function string `f` into bytecode
!
!@note This is not very efficient since it is parsing it twice
!      just to get the size of all the arrays.

    subroutine compile (me, f, var)

    implicit none

    class(fparser),intent(inout)             :: me
    character(len=*),intent(in)              :: f      !! function string
    character(len=*),dimension(:),intent(in) :: var    !! array with variable names

    integer :: istat  !! allocation status flag

    me%bytecodesize = 0
    me%immedsize    = 0
    me%stacksize    = 0
    me%stackptr     = 0

    ! compile string to determine size:
    call me%compile_substr (f,1,len_trim(f),var)

    allocate ( me%bytecode(me%bytecodesize),      &
               me%bytecode_ops(me%bytecodesize),  &
               me%immed(me%immedsize),            &
               me%stack(me%stacksize),            &
               stat = istat                       )

    if (istat /= 0) then
        call me%error_msg%add('*** Parser error: Memory allocation for byte code failed')
    else
        me%bytecodesize = 0
        me%immedsize    = 0
        me%stacksize    = 0
        me%stackptr     = 0
        call me%compile_substr (f,1,len_trim(f),var) ! compile string into bytecode
    end if

    end subroutine compile
!*******************************************************************************

!*******************************************************************************
!>
!  Add compiled byte to bytecode

    subroutine add_compiled_byte (me, b, num_args)

    implicit none

    class(fparser),intent(inout)  :: me
    integer,intent(in)            :: b          !! value of byte to be added
    integer, optional, intent(in) :: num_args

    integer :: args

    if (present(num_args)) then
        args = num_args
    else
        ! The required_args parameter array is not indexed from 1.
        if ( (b >= lbound(required_args, 1)) .and. (b <= ubound(required_args, 1)) ) then
            args = required_args(b)
        endif
    endif

    me%bytecodesize = me%bytecodesize + 1

    if (allocated(me%bytecode)) then

        !integer:
        me%bytecode(me%bytecodesize) = b

        !set the function pointer:
        ! [this replaces the original code which used
        !  a case statement during the evaluation]
        select case (b)
        case (cImmed);          me%bytecode_ops(me%bytecodesize)%f => cimmed_func
        case   (cNeg);          me%bytecode_ops(me%bytecodesize)%f => cneg_func
        case   (cAdd);          me%bytecode_ops(me%bytecodesize)%f => cadd_func
        case   (cSub);          me%bytecode_ops(me%bytecodesize)%f => csub_func
        case   (cMul);          me%bytecode_ops(me%bytecodesize)%f => cmul_func
        case   (cDiv);          me%bytecode_ops(me%bytecodesize)%f => cdiv_func
        case   (cPow);          me%bytecode_ops(me%bytecodesize)%f => cpow_func
        case   (cabs);          me%bytecode_ops(me%bytecodesize)%f => cabs_func
        case   (cExp);          me%bytecode_ops(me%bytecodesize)%f => cexp_func
        case (cLog10);          me%bytecode_ops(me%bytecodesize)%f => clog10_func
        case   (cLog);          me%bytecode_ops(me%bytecodesize)%f => clog_func
        case  (cSqrt);          me%bytecode_ops(me%bytecodesize)%f => csqrt_func
        case  (cSinh);          me%bytecode_ops(me%bytecodesize)%f => csinh_func
        case  (cCosh);          me%bytecode_ops(me%bytecodesize)%f => ccosh_func
        case  (cTanh);          me%bytecode_ops(me%bytecodesize)%f => ctanh_func
        case   (cSin);          me%bytecode_ops(me%bytecodesize)%f => csin_func
        case   (cCos);          me%bytecode_ops(me%bytecodesize)%f => ccos_func
        case   (cTan);          me%bytecode_ops(me%bytecodesize)%f => ctan_func
        case  (cAsin);          me%bytecode_ops(me%bytecodesize)%f => casin_func
        case  (cAcos);          me%bytecode_ops(me%bytecodesize)%f => cacos_func
        case (cAtan2);          me%bytecode_ops(me%bytecodesize)%f => catan2_func
        case  (cAtan)
            select case (args)
            case (1);           me%bytecode_ops(me%bytecodesize)%f => catan_func
            case (2);           me%bytecode_ops(me%bytecodesize)%f => catan2_func
            end select
        case (cPi);             me%bytecode_ops(me%bytecodesize)%f => cPi_func

        case(cCeil);            me%bytecode_ops(me%bytecodesize)%f => cceil_func
        case(cFloor);           me%bytecode_ops(me%bytecodesize)%f => cfloor_func
        case(cGamma);           me%bytecode_ops(me%bytecodesize)%f => cgamma_func
        case(cHypot);           me%bytecode_ops(me%bytecodesize)%f => chypot_func
        case(cMax);             me%bytecode_ops(me%bytecodesize)%f => cmax_func
        case(cMin);             me%bytecode_ops(me%bytecodesize)%f => cmin_func
        case(cMod);             me%bytecode_ops(me%bytecodesize)%f => cmod_func
        case(cModulo);          me%bytecode_ops(me%bytecodesize)%f => cmodulo_func
        case(cSign);            me%bytecode_ops(me%bytecodesize)%f => csign_func

        case (cIf);             me%bytecode_ops(me%bytecodesize)%f => cif_func
        case  default;          me%bytecode_ops(me%bytecodesize)%f => cdefault_func
        end select

    end if

    end subroutine add_compiled_byte
!*******************************************************************************

!*******************************************************************************
!>
!  return math item index, if item is real number, enter it into Comp-structure

    function mathitem_index (me, f, var) result (n)

    implicit none

    class(fparser),intent(inout)             :: me
    character(len=*),intent(in)              :: f     !! function substring
    character(len=*),dimension(:),intent(in) :: var   !! array with variable names

    integer :: n  !! byte value of math item

    n = 0
    if (len(f)==0) return ! error condition

    if (scan(f(1:1),'0123456789.') > 0) then ! check for begin of a number
        me%immedsize = me%immedsize + 1
        if (allocated(me%immed)) then
            me%immed(me%immedsize) = string_to_real (f)
        end if
        n = cimmed
    else  ! check for a variable
        n = variable_index (f, var)
        if (n > 0) n = varbegin+n-1
    end if

    end function mathitem_index
!*******************************************************************************

!*******************************************************************************
!>
!  Check if function substring F(b:e) is completely enclosed by a pair of parenthesis

    function completely_enclosed (f, b, e) result (res)

    implicit none

    character (len=*), intent(in) :: f    !! function substring
    integer,           intent(in) :: b,e  !! first and last pos. of substring
    logical :: res

    integer :: j,k

    res=.false.
    if (b==0 .or. e==0) return
    if (f(b:b) == '(' .and. f(e:e) == ')') then
        k = 0
        do j=b+1,e-1
            if (f(j:j) == '(') then
                k = k+1
            elseif (f(j:j) == ')') then
                k = k-1
            end if
            if (k < 0) exit
        end do
        if (k == 0) res=.true.    ! all opened parenthesis closed
    end if

    end function completely_enclosed
!*******************************************************************************

!*******************************************************************************
!>
!  Compile i-th function string `f` into bytecode

    recursive subroutine compile_substr (me, f, b, e, var)

    implicit none

    class(fparser),intent(inout)             :: me
    character(len=*),intent(in)              :: f     !! function substring
    integer,intent(in)                       :: b     !! begin position substring
    integer,intent(in)                       :: e     !! end position substring
    character(len=*),dimension(:),intent(in) :: var   !! array with variable names

    integer :: n
    integer :: b2,j,k,io
    character (len=*),parameter :: calpha = 'abcdefghijklmnopqrstuvwxyz'// &
                                            'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    integer :: arg_pos(max_func_args)
    integer :: num_args, iarg

    ! check for special cases of substring
    if (f(b:b) == '+') then                                     ! case 1: f(b:e) = '+...'
        call compile_substr (me, f, b+1, e, var)
        return
    elseif (completely_enclosed (f, b, e)) then                 ! case 2: f(b:e) = '(...)'
        call compile_substr (me, f, b+1, e-1, var)
        return
    elseif (scan(f(b:b),calpha) > 0) then
        n = mathfunction_index (f(b:e), var)
        if (n > 0) then
            b2 = b+index(f(b:e),'(')-1
            if (completely_enclosed(f, b2, e)) then             ! case 3: f(b:e) = 'fcn(...)'

                ! Determine the number of function arguments.
                call find_arg_positions (b2, f, num_args, arg_pos)

                if (num_args > 0) then
                    do iarg = 1, num_args
                        call compile_substr (me, f, b2+1, arg_pos(iarg), var)
                        if (iarg < num_args) then
                            me%stackptr = me%stackptr + 1
                            if (me%stackptr > me%stacksize) me%stacksize = me%stacksize + 1
                            b2 = arg_pos(iarg) + 1
                        end if
                    end do
                else
                    me%stackptr = me%stackptr + 1
                    if (me%stackptr > me%stacksize) me%stacksize = me%stacksize + 1
                end if

                call add_compiled_byte (me, n, num_args)
                return
            end if
        end if
    elseif (f(b:b) == '-') then
        if (completely_enclosed (f, b+1, e)) then               ! case 4: f(b:e) = '-(...)'
            call compile_substr (me, f, b+2, e-1, var)
            call add_compiled_byte (me, cneg)
            return
        elseif (scan(f(b+1:b+1),calpha) > 0) then
            n = mathfunction_index (f(b+1:e), var)
            if (n > 0) then
                b2 = b+index(f(b+1:e),'(')
                if (completely_enclosed(f, b2, e)) then          ! case 5: f(b:e) = '-fcn(...)'

                    ! Determine the number of function arguments.
                    call find_arg_positions (b2, f, num_args, arg_pos)

                    if (num_args > 0) then
                        do iarg = 1, num_args
                            call compile_substr (me, f, b2+1, arg_pos(iarg), var)
                            if (iarg < num_args) then
                                me%stackptr = me%stackptr + 1
                                if (me%stackptr > me%stacksize) me%stacksize = me%stacksize + 1
                                b2 = arg_pos(iarg) + 1
                            end if
                        end do
                    else
                        me%stackptr = me%stackptr + 1
                        if (me%stackptr > me%stacksize) me%stacksize = me%stacksize + 1
                    end if

                    call add_compiled_byte (me, n, num_args)
                    call add_compiled_byte (me, cneg)
                    return
                end if
            end if
        end if
    end if

    ! check for operator in substring: check only base level (k=0), exclude expr. in ()
    do io=cadd,cpow                                          ! increasing priority +-*/^
        k = 0
        do j=e,b,-1
            if (f(j:j) == ')') then
                k = k+1
            elseif (f(j:j) == '(') then
                k = k-1
            end if
        if (k == 0 .and. f(j:j) == operators(io) .and. is_binary_operator (j, f)) then
                if (any(f(j:j) == operators(cmul:cpow)) .and. f(b:b) == '-') then ! case 6: f(b:e) = '-...op...' with op > -
                    call compile_substr (me, f, b+1, e, var)
                    call add_compiled_byte (me, cneg)
                    return
                else                                                        ! case 7: f(b:e) = '...binop...'
                    call compile_substr (me, f, b, j-1, var)
                    call compile_substr (me, f, j+1, e, var)
                    call add_compiled_byte (me, operator_index(operators(io)))
                    me%stackptr = me%stackptr - 1
                    return
                end if
            end if
        end do
    end do

    ! check for remaining items, i.e. variables or explicit numbers
    b2 = b
    if (f(b:b) == '-') b2 = b2+1
    n = mathitem_index(me, f(b2:e), var)
    call add_compiled_byte (me, n)
    me%stackptr = me%stackptr + 1
    if (me%stackptr > me%stacksize) me%stacksize = me%stacksize + 1
    if (b2 > b) call add_compiled_byte (me, cneg)

    end subroutine compile_substr
!*******************************************************************************

!*******************************************************************************
!>
!  Check if operator `f(j:j)` in string `f` is binary operator.
!
!  Special cases already covered elsewhere:              (that is corrected in v1.1)
!  * operator character `f(j:j)` is first character of string (`j=1`)

    function is_binary_operator (j, f) result (res)

    implicit none

    integer,intent(in)           :: j       !! position of operator
    character (len=*),intent(in) :: f       !! string
    logical                      :: res     !! result

    integer :: k
    logical :: dflag,pflag

    res=.true.
    if (f(j:j) == '+' .or. f(j:j) == '-') then              ! plus or minus sign:
        if (j == 1) then                                    ! - leading unary operator ?
            res = .false.
        elseif (scan(f(j-1:j-1),',+-*/^(') > 0) then        ! - other unary operator ?   (or comma from multi-arg functions)
            res = .false.
        elseif (scan(f(j+1:j+1),'0123456789') > 0 .and. &   ! - in exponent of real number ?
                scan(f(j-1:j-1),'eEdD')       > 0) then
            dflag=.false.
            pflag=.false.
            k = j-1
            do while (k > 1)                                !   step to the left in mantissa
                k = k-1
                if (scan(f(k:k),'0123456789') > 0) then
                    dflag=.true.
                elseif (f(k:k) == '.') then
                    if (pflag) then
                        exit                                !   * exit: 2nd appearance of '.'
                    else
                        pflag=.true.                        !   * mark 1st appearance of '.'
                    endif
                else
                    exit                                    !   * all other characters
                end if
            end do
            if (dflag .and. (k == 1 .or. scan(f(k:k),',+-*/^(') > 0)) res = .false.  ! need the comma here too ??
        end if
    end if

    end function is_binary_operator
!*******************************************************************************

!*******************************************************************************
!>
!  Get real number from string.
!
!  Format: `[blanks][+|-][nnn][.nnn][e|E|d|D[+|-]nnn]`

    function string_to_real (str, ibegin, inext, error) result (res)

    implicit none

    real(wp)                       :: res     !! real number
    character (len=*),  intent(in) :: str     !! string
    integer, optional, intent(out) :: ibegin  !! start position of real number
    integer, optional, intent(out) :: inext   !! 1st character after real number
    logical, optional, intent(out) :: error   !! error flag

    integer  :: ib,in,istat
    logical  :: Bflag    !! True at begin of number in str
    logical  :: InMan    !! True in mantissa of number
    logical  :: Pflag    !! True after 1st '.' encountered
    logical  :: Eflag    !! True at exponent identifier 'eEdD'
    logical  :: InExp    !! True in exponent of number
    logical  :: DInMan   !! True if at least 1 digit in mant.
    logical  :: DInExp   !! True if at least 1 digit in exp.
    logical  :: err      !! Local error flag

    Bflag=.true.
    InMan=.false.
    Pflag=.false.
    Eflag=.false.
    InExp=.false.
    DInMan=.false.
    DInExp=.false.
    ib = 1
    in = 1
    do while (in <= len_trim(str))
        select case (str(in:in))
        case (' ')                                   ! Only leading blanks permitted
            ib = ib+1
            if (InMan .or. Eflag .or. InExp) exit
        case ('+','-')                               ! Permitted only
            if (Bflag) then
                InMan=.true.; Bflag=.false.          ! - at beginning of mantissa
            elseif (Eflag) then
                InExp=.true.; Eflag=.false.          ! - at beginning of exponent
            else
                exit                                 ! - otherwise STOP
            endif
        case ('0':'9')                               ! Mark
            if (Bflag) then
                InMan=.true.; Bflag=.false.          ! - beginning of mantissa
            elseif (Eflag) then
                InExp=.true.; Eflag=.false.          ! - beginning of exponent
            endif
            if (InMan) DInMan=.true.                 ! Mantissa contains digit
            if (InExp) DInExp=.true.                 ! Exponent contains digit
        case ('.')
            if (Bflag) then
                Pflag=.true.                         ! - mark 1st appearance of '.'
                InMan=.true.; Bflag=.false.          !   mark beginning of mantissa
            elseif (InMan .and..not.Pflag) then
                Pflag=.true.                         ! - mark 1st appearance of '.'
            else
                exit                                 ! - otherwise STOP
            end if
        case ('e','E','d','D')                       ! Permitted only
            if (InMan) then
                Eflag=.true.; InMan=.false.          ! - following mantissa
            else
                exit                                 ! - otherwise STOP
            endif
        case default
            exit                                     ! STOP at all other characters
        end select
        in = in+1
    end do
    err = (ib > in-1) .or. (.not.DInMan) .or. ((Eflag.or.InExp).and..not.DInExp)
    if (err) then
       res = zero
    else
       read(str(ib:in-1),*,iostat=istat) res
       err = istat /= 0
    end if
    if (present(ibegin)) ibegin = ib
    if (present(inext))  inext  = in
    if (present(error))  error  = err

    end function string_to_real
!*******************************************************************************

!*******************************************************************************
!>
!  Transform upper case letters in `str1` into
!  lower case letters, result is `str2`.

    pure elemental subroutine to_lowercase (str1, str2)

    implicit none

    character (len=*),  intent(in) :: str1
    character (len=*), intent(out) :: str2

    integer :: j,k
    character(len=*), parameter :: lc = 'abcdefghijklmnopqrstuvwxyz'
    character(len=*), parameter :: uc = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

    str2 = str1
    do j=1,len_trim(str1)
        k = index(uc,str1(j:j))
        if (k > 0) str2(j:j) = lc(k:k)
    end do

    end subroutine to_lowercase
!*******************************************************************************

!*******************************************************************************
    end module function_parser
!*******************************************************************************