module cmdline_arguments

  use iso_varying_string
  use string_functions, only: lcase, operator(.lcase.), int, real, double, split, join, count, assignment(=)
  use variable_array, only: shift, unshift, push, splice
  use hash_table, hash_value => value
  use file_functions, only: stderr, stdout
  use precision

  implicit none

  ! private

  ! This module provides an easy to use interface to the iargc and getarg
  ! compatibilty functions provided by the Intel Fortran Compiler (and
  ! supported by most Fortran compilers). These are used to get access to
  ! command line arguments.
  !
  ! Global variables are used to keep track of the argument object for each
  ! particular invocation. This is justifiable as no program has more than 
  ! one set of command line arguments (well I can't think of an instance
  ! where this would be the case, so I'm not supporting it). This makes the
  ! module much simpler and easier to use.

  !! $Log: cmdline_arguments.f90,v $
  !! Revision 1.6  2007/11/12 03:34:31  aidan
  !! Changed internal routine 'current_argument' to 'get_current_argument' to
  !! avoid name conflicts.
  !!
  !! Revision 1.5  2005/06/28 03:51:26  aidan
  !! Added support for double precision real command line options.
  !!
  !! Revision 1.4  2004/09/07 06:30:13  aidan
  !! Removed the 'feature' added in v1.3! The function I altered was always
  !! intended as a stand alone, the one called by the user checked for
  !! initialisation and then called number_of_arguments. Added a global debug
  !! variable and some debugging statements. Fixed a bug in process_options -- the
  !! main loop didn't check to make sure there were still elements in the argument
  !! list before shifting them off. Now the we have a do while construct. Hence I
  !! removed the check for the number of arguments before the loop. This was a
  !! good thing anyway, as we definitely wanted the whole subroutine traversed as
  !! some important status variables are set a the end of it.
  !!
  !! Revision 1.3  2004/09/07 03:48:52  aidan
  !! Added a call to initialise in number_of_arguments so I could call this as
  !! the first cmdline_arguments procedure in a program, i.e. check for ANY
  !! arguments at all, both options and arguments. I had to take out a reference
  !! to nargs in intialise_argument_object, as this created unintentional
  !! recursion.
  !!
  !! Revision 1.2  2004/04/05 23:21:10  aidan
  !! Added command line option parsing. This is a major upgrade to this module.
  !! The option parsing makes heavy use of the hash table module to allow easy
  !! lookup of permitted options, optional values etc. Type coercion of option
  !! values is accomplished through overloading the assignment operator, i.e. an
  !! example would be '--max=8.5', where the value would be stored internally as
  !! a string '8.5' but could be retrieved as a real number by assigning the
  !! result of the 'get_value' function to a real variable. In the same way
  !! comma separated option values, e.g. min=1,2,3,4, can be assigned to an
  !! array and the string will be split on the commas, and the individual values
  !! placed in the array. Care must be taken with whitespace -- the built in
  !! command line processing of the compiler is used to return the command line
  !! arguments, and this usually uses whitespace as the delimiter between
  !! different arguments, so whitespace must be quoted if it is to appear in an
  !! option value. In general whitespace in an option name should be avoided.
  !!
  !! Revision 1.1  2003/09/11 04:30:13  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: cmdline_arguments.f90,v 1.6 2007/11/12 03:34:31 aidan Exp $"

  character(len=1) :: delimiter

  ! We define an argument object which contains an array for the argument
  ! list and two integers, the number of arguments (the length of the 
  ! argument list) and the current argument, which is used to keep track
  ! of our current position when traversing the list
  type argument_object
     integer :: number_of_arguments, current_argument
     type (varying_string), dimension(:), pointer :: argument_list
     type (varying_string), dimension(:), pointer :: bad_option_list => null()
     ! type (varying_string), dimension(:), pointer    :: option_list
     ! type (varying_string), dimension(:)             :: value_list
     ! logical, dimension(:), pointer                  :: value_present
     type (hash_table_object) :: option_list
     type (hash_table_object) :: value_present
     logical :: initialised = .FALSE., options_processed = .FALSE.
  end type argument_object

  ! These types are defined to cope with command line options. One is
  ! an option value, the other the option itself.
  type option_value
     ! integer(kind=1), dimension(:), pointer :: value => null()
     private
     type (varying_string) :: value
  end type option_value

  type single_argument_object
     type (varying_string) :: argument
  end type single_argument_object

  ! We define a number of interfaces to internal routines -- provides 
  ! easy to type and remember versions of longer names

  ! This is not a public procedure -- just a useful little routine which
  ! checks the initialisation of the argument object and calls the
  ! necessary routine if we are not initialised
  interface initialise
     module procedure arg_object_initialised
  end interface

  ! A function which returns a boolean value, true if there are still
  ! arguments remaining on our "stack", false otherwise
  interface have_args
     module procedure arguments_remaining
  end interface

  ! This is the pidgeon-pair to have_args above -- calling this will
  ! "unshift" the next argument from our "stack" of arguments
  interface next_arg
     module procedure next_argument
  end interface

  interface num_args
     module procedure narguments_in_this_object
  end interface

  interface all_args
     module procedure return_all_arguments
  end interface

  ! Return number of unrecognised options
  interface num_bad_options
     module procedure return_num_bad_options
  end interface

  ! If we have an error processing the options this routine allows access
  ! to the bad command line options that were passed to the program
  interface bad_options
     module procedure return_bad_options
  end interface

  interface option_exists
     module procedure option_exists
  end interface

  interface has_value
     module procedure option_has_value
  end interface

  interface get_options
     module procedure process_options
  end interface

  interface get_value
     module procedure get_option_value !, get_option_values, get_option_value_as_character, &
           !get_option_value_as_int, get_option_values_as_ints, get_option_value_as_real
  end interface

  ! interface get_value_array
  !    module procedure get_option_value_array
  ! end interface

!!$  ! Need to provide an explicit interface to the iargc function otherwise
!!$  ! the compiler will complain due to our use of implicit none
!!$  interface 
!!$     integer(kind=4) function iargc()
!!$     end function iargc
!!$  end interface
!!$
!!$  ! Do not *need* an interface to getarg, but it is good practice to 
!!$  ! provide one so the compiler can spot any errors we might make
!!$  interface getarg
!!$     subroutine getarg(argindex, argument)
!!$       integer(4), intent(in) :: argindex
!!$       character(len=*), intent(out) :: argument
!!$     end subroutine getarg
!!$     subroutine getarg_dvf(n, buffer, status)
!!$       integer(2), intent(in) :: n
!!$       character(len=*), intent(out) :: buffer
!!$       integer(2), optional :: status
!!$     end subroutine getarg_dvf
!!$  end interface

  ! These interfaces are not public ...
  interface nargs
     module procedure number_of_arguments
  end interface

  interface current
     module procedure get_current_argument
  end interface

  ! All the magic in converting from the internal formats to character, integer
  ! and real values is done with assignment
  interface assignment(=)
     module procedure assign_option_to_int, assign_option_to_intarray, &
          assign_option_to_real, assign_option_to_realarray, & 
          assign_option_to_double, assign_option_to_doublearray, & 
          assign_option_to_string, assign_option_to_varstring, assign_option_to_vstring_array
  end interface

  ! Public routines
  public :: have_args, next_arg, all_args, num_args
  public :: bad_options, get_options, get_value, option_exists, has_value

  ! Public types
  public :: option_value, argument_object

  ! Overloaded operators
  public :: assignment(=)

  ! Global variables

  ! We have a single command line argument object for each program which uses this module.
  ! We set the save modifier so that we can initialise it once and then change the status 
  ! of the initialised flag
  type (argument_object), save :: this_argument_object

  logical, parameter :: debug = .FALSE.

contains

  subroutine arg_object_initialised()

    ! Check if we need to initialise the global argument object. Do so
    ! if not already done .. otherwise do nothing.
    if (.not. this_argument_object%initialised) then
       call initialise_argument_object(this_argument_object)
    end if

  end subroutine arg_object_initialised

  logical function arguments_remaining()

    ! This is the "main" function in this module. It is the first 
    ! function that is called, so it initialises our global argument
    ! object on the first pass, and then tells the user if there are
    ! any arguments. Subsequent invocations are usually after
    ! the 'next_arg' function, which returns the current argument and
    ! then increments the current argument counter.

    ! This function is used in loops like so
    !
    !         while (have_args()) do
    !            argument = next_arg()
    !         end do
    ! 
    ! or in single statements like
    !
    !         if (have_args()) argument = next_arg()
    ! 

    ! On the first pass through this function we will initialise the
    ! the global argument object (this_argument_object)
    call initialise()

    ! Set default for arguments_remaining
    arguments_remaining = .true.

    ! Set to false if the number of the current object is greater than
    ! or equal to the number of arguments
    if (current(this_argument_object) >= nargs(this_argument_object)) &
         arguments_remaining = .false.

  end function arguments_remaining

  function next_argument()

    ! Return the next argument from the argument array of the globally defined 
    ! argument object 'this_argument_object'

    ! character(len=len(this_argument_object%argument_list(current(this_argument_object)+1))) :: next_argument
    type (varying_string) :: next_argument

    ! On the first pass through this function we will initialise the
    ! the global argument object (this_argument_object)
    call initialise()

    ! This next line (split over several lines to aid readability) does the 
    ! following tasks: 
    !
    !     o gets the current argument number and increments it by one
    !     o passes this incremented value back into the argument object
    !     o uses the return value of setting the incremented value as
    !       an index of the argument array
    !     o sets the next argument to this member of the argument array
    !
    next_argument = this_argument_object%argument_list(                      &
                                 current(this_argument_object,               &
                                         current(this_argument_object) + 1   &
                                        )                                    &
                                 )

  end function next_argument

  function return_all_arguments() result(argument_array)

    type (varying_string), dimension(size(this_argument_object%argument_list)) :: argument_array

    ! On the first pass through this function we will initialise the
    ! the global argument object (this_argument_object)
    call initialise()

    argument_array = this_argument_object%argument_list

  end function return_all_arguments

  integer function narguments_in_this_object()

    ! Simple function to return the number of arguments

    ! On the first pass through this function we will initialise the
    ! the global argument object (this_argument_object)
    call initialise()

    narguments_in_this_object = nargs(this_argument_object)

  end function narguments_in_this_object

  subroutine initialise_argument_object(arg_object)

    ! This is the workhorse routine which grabs the command line
    ! arguments and stores them internally in the argument object
    ! passed to it from 'arguments_remaining'
    
    type (argument_object), intent(out) :: arg_object

    integer    :: error, current_argument, number_of_args
    integer(2) :: i, arglength
    character(len=1000) :: buffer

    ! Initialise the argument list to a null pointer
    nullify(arg_object%argument_list)
    nullify(arg_object%bad_option_list)

    ! Set the current argument index to zero
    current_argument = 0
    current_argument = current(arg_object,current_argument)

    ! print *,"number of arguments = ",iargc()
    number_of_args = iargc()

    ! Grab the number of arguments and assign this internally in
    ! one step (nargs returns the number we put in as a matter of course)
    ! if (nargs(arg_object,iargc()) > 0) then
    if (number_of_args > 0) then

       ! Allocate some memory for the argument array
       allocate(arg_object%argument_list(number_of_args), stat=error)
       if (error /= 0) stop 'Error allocating memory for argument list'

       ! Cycle through the arguments and save in the argument array
       ! do i=1,nargs(arg_object) 
       do i=1,number_of_args
          call getarg(i,buffer)
          ! call getarg(n=i,buffer=buffer,status=arglength)
          if (debug) print '(I0," : ",A)',i,buffer(1:len_trim(buffer))
          arg_object%argument_list(i) = buffer(1:len_trim(buffer))
       enddo

    endif

    ! Set the initialised flag to true
    arg_object%initialised = .true.
    
  end subroutine initialise_argument_object

  integer function number_of_arguments(arg_object)

    ! Simple function to return the number of arguments

    type (argument_object), intent(inout) :: arg_object

    number_of_arguments = size(arg_object%argument_list) 

  end function number_of_arguments

  integer function get_current_argument(arg_object, number)

    ! Allows us to set and return the number of the
    ! current argument -- this is for traversing the argument array

    type (argument_object), intent(inout) :: arg_object
    integer, intent(in), optional         :: number

    ! If we specify a number then set this to tbe the current argument
    if (present(number)) then
       arg_object%current_argument = number
    end if

    ! Always return the value in the current argument
    get_current_argument = arg_object%current_argument 

  end function get_current_argument

  subroutine process_options (options, status)

    ! This routine allows post-processing of the command line arguments
    ! to extract command line options, i.e. arguments that begin with
    ! '--' or '-'.

    ! The list of valid options is passed in as an array of type
    ! varying string. If an option is used that is not in this list
    ! then an error is flagged. If status is passed as an argument
    ! to the subroutine then this is set to a value > 0, otherwise a
    ! fatal error occurs.

    ! Input arguments
    type (varying_string), dimension(:), intent(in) :: options
    integer, intent(out), optional :: status

    ! Local variables
    type (hash_table_object) :: valid_options
    type (varying_string) :: buffer, option, option_value

    integer :: i, option_start, eq_position

    ! Make sure we haven't already processed options
    if (this_argument_object%options_processed) then
       write(stderr,*) 'CMDLINE_ARGUMENTS :: Error! Attempting to process options a second time -- aborted'
       return
    end if

    ! Initialise our argument object if it hasn't been already
    call initialise()

    ! Initialise some hash tables to store information
    call new(valid_options,8)
    call new(this_argument_object%option_list, 6)
    call new(this_argument_object%value_present, 6)

    ! Loop over the list of valid options and add them to our hash table
    do i=1, size(options)
       ! print *,(.lcase.(char(options(i))))
       if (.NOT. add(valid_options,(.lcase.(trim(char(options(i))))),1)) &
       stop "CMDLINE_ARGUMENTS :: Error! Could not add option to table of valid options"
    end do

    ! Initialise the return status to zero, i.e. no error
    if (present(status)) status = 0

    ! In this loop we cycle through all the command line arguments
    ! until we run out options, or encounter a lone double dash '--'
    do while (nargs(this_argument_object) >= 1)

       ! Grab the first argument off the list of arguments
       buffer = shift(this_argument_object%argument_list)

       if (buffer == '--') then 

          ! This is the end of the options, so stop processing them
          exit

       else if (scan(adjustl(buffer),'-') == 1) then 

          ! We have an option, so get rid on any leading space 
          buffer = adjustl(buffer)

          ! Check if we have a single or double dash at the start of
          ! our option
          if (extract(buffer,2,2) == '-') then
             option_start = 3
          else
             option_start = 2
          end if

          ! The command is the string between the leading dashes and
          ! the end of the characters, or an equals sign, whichever
          ! is closer.
          eq_position = index(buffer,'=')
          if (eq_position == 0) eq_position = len_trim(buffer)+1
          option = extract(buffer,option_start,eq_position-1)

          ! Force our option to be lower case
          ! call lcase(option)
          option = (.lcase. char(option))

          ! Check the option is legal by looking for it in the hash table
          ! options passed in by the user

          if (debug) print *,'Is this a valid option: ', char(option)

          if (exists(valid_options,char(option))) then
             ! We have a valid option ... so add this option on to our 
             ! internal list of options passed to the program

             if (debug) print *,'yes!'

             ! But first establish if we were passed a value, if not just 
             ! set it to an empty string
             if (eq_position > len(buffer)) then
                option_value = ''
             else
                option_value = extract(buffer,eq_position+1,len(buffer))
                ! This hash table has an entry for every option that also passed
                ! a value -- the value used is the number of delimiter characters 
                ! that occur in the value, in case we need to split the value into 
                ! an array when accessing it
                if (.NOT. add(this_argument_object%value_present, char(option), count(char(option_value),delimiter))) &
                     stop "CMDLINE_ARGUMENTS :: Error! Could not add option to value present hash"
             end if

             ! Add our option/value key pair to the option hash table
             if (.NOT. add(this_argument_object%option_list, char(option), char(option_value))) &
                  stop "CMDLINE_ARGUMENTS :: Error! Could not add option to hash"
          else
             ! This is an illegal option, so add this to our list of
             ! bad options
             if (push(this_argument_object%bad_option_list, option) < 0) &
                  stop "CMDLINE_ARGUMENTS :: Error! Could not add to bad option list"
          end if

       else

          ! This is not an option .. so put it back on to the front of the
          ! option list
          if (unshift(this_argument_object%argument_list,char(buffer)) < 0) &
               stop "CMDLINE_ARGUMENTS :: Error! Cannot reconstitute the argument list"

          ! And stop looking for command line options
          exit

       end if
       
    end do

    ! Check if we have some bad options
    if (associated(this_argument_object%bad_option_list) .and. size(this_argument_object%bad_option_list) > 0) then
       ! Yes we do .. so return the number of bad options in the status variable
       ! if we have passed this to the subroutine, or barf and stop the program
       if (present(status)) then
          status = size(this_argument_object%bad_option_list)
       else
          write(stderr,*) 'ERROR! Bad options: '
          ! write(stderr,*) 'ERROR! Bad options: ', join(this_argument_object%bad_option_list, " ")
          stop
       end if
    end if

    if (debug) print *,'Number of args left: ', size(this_argument_object%argument_list)

    ! Set flag to indicate we have processed the options
    this_argument_object%options_processed = .TRUE.
    
  end subroutine process_options

  integer function return_num_bad_options ()

    ! This function returns an array of varying strings which were passed as command
    ! line options but did not match any of the specified valid options.

    return_num_bad_options = size(this_argument_object%bad_option_list)

  end function return_num_bad_options

  function return_bad_options () result(bad_options)

    ! This function returns an array of varying strings which were passed as command
    ! line options but did not match any of the specified valid options.

    type (varying_string), dimension(size(this_argument_object%bad_option_list)) :: bad_options

    if (.NOT. this_argument_object%options_processed) then
       write(stderr,*) 'CMDLINE_ARGUMENTS :: Error! Cannot access bad argument list without processing options first!'
       stop
    end if

    bad_options = this_argument_object%bad_option_list

  end function return_bad_options

  logical function option_exists (option)

    ! See if a command line option was specified

    character(len=*), intent(in) :: option
    
    option_exists = exists(this_argument_object%option_list, option)
    
  end function option_exists

  logical function option_has_value (option)

    ! See if a command line option also specified a value, e.g. --max=5.4

    character(len=*), intent(in) :: option
    
    option_has_value = exists(this_argument_object%value_present, option)
    
  end function option_has_value

  function get_option_value (option) result(value)

    ! Return the value specified for a command line option

    ! Input variables
    character(len=*), intent(in)      :: option
    
    ! Return type of function
    type (option_value) :: value

    if (exists(this_argument_object%option_list, option)) then
       if (exists(this_argument_object%value_present, option)) then
          value%value = hash_value(this_argument_object%option_list)
       end if
    end if
    
  end function get_option_value

  ! The magic in the option retrieval part of this module is done with
  ! overloaded assignment. The get_option_value routine returns an 
  ! option value type, which we can then coerce into a number of
  ! different types, depending on what is on the left hand side of the
  ! assignment.

  subroutine assign_option_to_int (intval, value)

    integer, intent(out)            :: intval
    type (option_value), intent(in) :: value

    ! The option value is stored as a varying string. We use the 
    ! overloaded assignment from string functions to convert this
    ! to an integer
    intval = value%value
    
  end subroutine assign_option_to_int

  subroutine assign_option_to_intarray (intval, value)

    type (option_value), intent(in) :: value
    integer, dimension(count(char(value%value),",")+1), intent(out) :: intval
    character(len=len(value%value)) :: buffer
    
    buffer = char(value%value)
    read (buffer,*) intval

    ! Can't use this more concise formulation, compilers don't grok it
    ! intval = split(char(value%value),",")
    
  end subroutine assign_option_to_intarray

  subroutine assign_option_to_real (realval, value)

    real, intent(out)            :: realval
    type (option_value), intent(in) :: value
    
    realval = value%value
    
  end subroutine assign_option_to_real

  subroutine assign_option_to_realarray (realval, value)

    type (option_value), intent(in) :: value
    real, dimension(count(char(value%value),",")+1), intent(out) :: realval
    character(len=len(value%value)) :: buffer
    
    ! realval = split(char(value%value),",")
    buffer = char(value%value)
    read (buffer,*) realval
    
  end subroutine assign_option_to_realarray

  subroutine assign_option_to_double (realval, value)

    real(kind=rd_kind), intent(out) :: realval
    type (option_value), intent(in) :: value
    character(len=len(value%value)) :: buffer
    
    realval = value%value
    
  end subroutine assign_option_to_double

  subroutine assign_option_to_doublearray (realval, value)

    type (option_value), intent(in) :: value
    real(kind=rd_kind), dimension(count(char(value%value),",")+1), intent(out) :: realval
    character(len=len(value%value)) :: buffer
    
    ! realval = split(char(value%value),",")
    buffer = char(value%value)
    read (buffer,*) realval
    
  end subroutine assign_option_to_doublearray

  subroutine assign_option_to_string (stringval, value)

    type (option_value), intent(in) :: value
    character(len=len(value%value)), intent(out) :: stringval
    
    stringval = char(value%value)
    
  end subroutine assign_option_to_string

  subroutine assign_option_to_varstring (stringval, value)

    type(varying_string), intent(out) :: stringval
    type (option_value), intent(in) :: value
    
    stringval = value%value
    
  end subroutine assign_option_to_varstring

  subroutine assign_option_to_vstring_array (varstringval, value)

    type (option_value), intent(in) :: value
    type (varying_string), dimension(count(char(value%value),",")+1), intent(out) :: varstringval

    character(len=len(value%value)), dimension(size(varstringval))    :: stringval
    character(len=len(value%value)) :: buffer
    integer :: i
    
    buffer = char(value%value)
    ! print *,buffer
    stringval = split(buffer,",")
    do i = 1, size(stringval)
       varstringval(i) = trim(stringval(i))
       ! print *,i,varstringval(i)
    end do
    
    ! stringval = split(char(value%value),",")
    
  end subroutine assign_option_to_vstring_array

end module cmdline_arguments
