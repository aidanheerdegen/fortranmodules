module keyword_class

  use iso_varying_string
  use string_functions, only: operator(.ucase.), int, real, double, split, join, count, assignment(=)
  use variable_array, only: shift, unshift, push, splice
  use hash_table, hash_value => value
  use file_functions, only: stderr, stdout, read_buffer, open
  use precision

  implicit none

  private

  ! This module provides a method for parsing input files with keywords
  ! and extracting the data into a hash based lookup table for easy access.

  !! $Log: keyword_class.f90,v $
  !! Revision 1.5  2009/07/02 01:31:53  aidan
  !! Fixed logical error with multiple values. Wasn't checking that
  !! subsequent keywords had a value -- all instances of multiple keywords
  !! must specify a value.
  !!
  !! Revision 1.4  2009/05/21 03:28:49  aidan
  !! Fixed a bug with multiple value keywords. It would only have managed
  !! to store a maximum of two values. Made varying_string wrappers for all
  !! routines which took a character variable. Wrote tests for all new
  !! features.
  !!
  !! Revision 1.3  2009/05/15 01:48:23  aidan
  !! Fixed bug where the intialising routine asked for way too big an
  !! initial hash (used the length of the key array rather than the
  !! log2 of this value). Had to add a log2 function to accomplish this.
  !!
  !! Revision 1.2  2008/11/03 03:55:08  aidan
  !! Working version of keyword class. Implements multiple values using a
  !! kludgey work-around where multiple values are concatenated with a
  !! null character delimiter.
  !!
  !! Revision 1.1  2008/10/09 05:13:50  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: keyword_class.f90,v 1.5 2009/07/02 01:31:53 aidan Exp $"

  character(len=1) :: delimiter = char(38)

  ! We define an argument object which contains an array for the argument
  ! list and two integers, the number of arguments (the length of the 
  ! argument list) and the current argument, which is used to keep track
  ! of our current position when traversing the list
  type keyword_object
     integer :: number_of_keywords, current_argument
     type (varying_string), dimension(:), pointer :: bad_keyword_list
     type (hash_table_object) :: valid_keywords
     type (hash_table_object) :: keyword_list
     type (hash_table_object) :: value_present
     logical :: initialised = .FALSE.
  end type keyword_object

  type keyword_value
     type (varying_string) :: value
     integer :: n
  end type keyword_value

  interface initialise
     module procedure initialise_keyword_object
  end interface

  interface add
     module procedure add_keyword, add_keyword_novalue, add_varstr_keyword
     module procedure add_varstr_key_varstr_word, add_key_varstr_word, add_varstr_keyword_novalue
  end interface

  interface delete
     module procedure delete_keyword, delete_varstr_keyword
  end interface

  interface exists
     module procedure keyword_exists, varstr_keyword_exists
  end interface

  interface has_value
     module procedure keyword_has_value, varstr_keyword_has_value
  end interface

  interface num_value
     module procedure get_num_value, get_num_value_varstr
  end interface

  interface get_value
     module procedure get_keyword_value, get_one_keyword_value, get_keyword_value_varstr, get_one_keyword_value_varstr
  end interface

  interface parse
     module procedure parse_line, parse_varstr_line
  end interface

  interface read
     module procedure open_then_parse, parse_file, open_then_parse_varstr
  end interface

  interface print
     module procedure print_keyword_object, open_then_print, open_then_print_varstr
  end interface

  ! All the magic in converting from the internal formats to character, integer
  ! and real values is done with assignment
  interface assignment(=)
     module procedure assign_keyword_to_int, assign_keyword_to_intarray, &
          assign_keyword_to_real, assign_keyword_to_realarray, & 
          assign_keyword_to_double, assign_keyword_to_doublearray, & 
          assign_keyword_to_string, assign_keyword_to_varstring, assign_keyword_to_vstring_array
  end interface

  ! Public routines
  public :: initialise, add, delete, exists, has_value, get_value, read, parse , print, num_value, split_keyvalue !, log2

  ! Public variables
  public :: keyword_object, keyword_value

  ! Overloaded operators
  public :: assignment(=)

  ! Global variables
  logical, parameter :: debug = .false.

contains

  logical function parse_line (key_object, line) result(parsed_ok)

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: line

    ! Local variables
    integer :: pos
    character(len=len(line)) :: buf

    parsed_ok = .false.

    if (len_trim(line) < 0) return

    buf = adjustl(line)

    ! Get rid of spaces at the beginning of the line
       
    pos = index(buf,' ')

    if (buf(pos:) == " ") then
       parsed_ok = add (key_object, buf(:pos-1))
    else
       parsed_ok = add (key_object, buf(:pos-1), trim(adjustl(buf(pos:))))
    end if
       
  end function parse_line
  
  logical function parse_varstr_line (key_object, line) result(parsed_ok)

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: line

    parsed_ok = parse(key_object, char(line))

  end function parse_varstr_line

  logical function parse_file (key_object, unit, quiet) result(parsed_ok)

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    integer, intent(in)                  :: unit
    logical, intent(in), optional        :: quiet

    ! Local variables
    integer :: pos
    character(len=1000) :: record
    logical :: end_of_file, squawk

    parsed_ok = .false.

    if (present(quiet)) then
       squawk = .not. quiet
    else
       squawk = .true.
    end if

    do
       ! Read through the file until we reach the end 
       call read_buffer(unit, record, end_of_file, comment="!#", removeblanks=.TRUE., inline=.TRUE.)
       if (end_of_file) exit
       if (.not. parse_line(key_object, record)) then
          if (squawk) write(stderr,*) 'KEYWORD_CLASS :: error at this line: ',record
          return
       end if
    end do

    parsed_ok = .true.
       
  end function parse_file
  
  logical function open_then_parse(key_object, file, quiet)

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: file
    logical, intent(in), optional        :: quiet

    ! Internal variables
    logical :: opened_ok
    integer :: unit
    logical :: lquiet


    open_then_parse = .false.

    ! Open the file
    unit = open(file, opened_ok, status='old')

    if (opened_ok) then
       if (present(quiet)) then
          open_then_parse = read(key_object, unit, quiet)
       else
          open_then_parse = read(key_object, unit)
       end if
       close(unit)
    else
       write (stderr,*)'KEYWORD_CLASS :: Unable to find file: ',file
    end if

  end function open_then_parse

  logical function open_then_parse_varstr(key_object, file)

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: file

    open_then_parse_varstr = read(key_object, char(file))

  end function open_then_parse_varstr

  logical function add_keyword (key_object, keyword, value) result(added_ok)

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    character(len=*), intent(in)         :: value

    ! Local variables
    character(len=len_trim(keyword)) :: key
    type (varying_string) :: thisvalue
    integer :: numvalues

    added_ok = .false.

    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    ! Force our option to be lower case
    ! call lcase(option)
    key = (.ucase. trim(keyword))

    if (exists(key_object%valid_keywords,key)) then
       ! We have a valid option ... so add this option on to our 
       ! internal list of keywords passed to the program

       if (debug) print *,'Valid keyword: ',key

       if (exists(key_object, key)) then
          if (.not. has_value(key_object, key)) then
             write (stderr,*) 'KEYWORD_VALUE :: Error! Trying to add another value for keyword that had no value.'
             return
          end if
          numvalues = hash_value(key_object%value_present)
          numvalues = numvalues + 1
          ! Can't use call to interface get_value routine as spits it out as a weird
          ! return value which gets coerced into just the first value ... make sense? 
          ! So .. direct access is the name of the game
          ! value = get_value(key_object, key)
          thisvalue = hash_value(key_object%keyword_list)
          if (debug) print *,'Current value: ',char(thisvalue)
          thisvalue = thisvalue // delimiter // trim(value)
          if (debug) print *,'New value: ',char(thisvalue)
          ! Delete entries for this key/value combination
          if (.NOT. delete(key_object%value_present)) &
               stop "KEYWORD_PARSER :: Error! Could not clear option in value present hash"
          if (.NOT. delete(key_object%keyword_list)) &
               stop "KEYWORD_PARSER :: Error! Could not clear option in keyword/value hash"
       else
          numvalues = 1
          thisvalue = trim(value)
       end if

       if (debug) print *,'value: ',char(thisvalue)
       if (debug) print *,'numvalues: ',numvalues

       ! This hash table has an entry for every option that also passed
       ! a value -- the value used is the number of delimiter characters 
       ! that occur in the value, in case we need to split the value into 
       ! an array when accessing it
       if (.NOT. add(key_object%value_present, key, numvalues)) &
            stop "KEYWORD_PARSER :: Error! Could not add option to value present hash"

       ! Add our keyword/value key pair to the keyword hash table
       if (.NOT. add(key_object%keyword_list, key, char(thisvalue))) &
            stop "KEYWORD_PARSER :: Error! Could not add option to keyword/value hash"

       added_ok = .true.
    else
       ! This is an illegal option, so add this to our list of
       ! bad keywords
       if (push(key_object%bad_keyword_list, key) < 0) &
            stop "CMDLINE_ARGUMENTS :: Error! Could not add to bad option list"
    end if

  end function add_keyword

  logical function add_keyword_novalue (key_object, keyword) result(added_ok)

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword

    ! Local variables
    character(len=len_trim(keyword)) :: key

    added_ok = .false.

    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    ! Force our option to be upper case
    key = (.ucase. trim(keyword))

    if (exists(key_object%valid_keywords,key)) then
       ! We have a valid option ... so add this option on to our 
       ! internal list of keywords passed to the program

       if (debug) print *,'Valid keyword: ',key

       if (exists(key_object, key)) then
          write (stderr,*) 'KEYWORD_CLASS ::: Error!'
          write (stderr,*) 'Encountered keyword with no value but this keyword has already been specified!'
          write (stderr,*) 'Multiple instances of keywords must all have a value.'
          return
       end if

       ! Add our keyword to the keyword hash table (with an empty value)
       if (.NOT. add(key_object%keyword_list, key, '')) &
            stop "KEYWORD_PARSER :: Error! Could not add keyword to keyword/value hash "

       added_ok = .true.
    else
       ! This is an illegal option, so add this to our list of
       ! bad keywords
       if (push(key_object%bad_keyword_list, key) < 0) &
            stop "CMDLINE_ARGUMENTS :: Error! Could not add to bad option list "
    end if

  end function add_keyword_novalue

  logical function add_varstr_keyword (key_object, keyword, value) result(added_ok)

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword
    character(len=*), intent(in)         :: value
    
    added_ok = add(key_object, char(keyword), value)

  end function add_varstr_keyword

  logical function add_varstr_key_varstr_word (key_object, keyword, value) result(added_ok)

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword
    type (varying_string), intent(in)    :: value

    added_ok = add(key_object, char(keyword), char(value))

  end function add_varstr_key_varstr_word

  logical function add_key_varstr_word (key_object, keyword, value) result(added_ok)

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    type (varying_string), intent(in)    :: value
    
    added_ok = add(key_object, keyword, char(value))

  end function add_key_varstr_word

  logical function add_varstr_keyword_novalue (key_object, keyword) result(added_ok)

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword

    added_ok = add(key_object, char(keyword))

  end function add_varstr_keyword_novalue

  subroutine initialise_keyword_object(key_object, keywords)

    ! This is the workhorse routine which grabs the command line
    ! arguments and stores them internally in the argument object
    ! passed to it from 'arguments_remaining'
    
    type (keyword_object), intent(out) :: key_object
    type (varying_string), dimension(:), intent(in) :: keywords

    ! Local variables
    integer :: i, hashsize

    ! Initialise the argument list to a null pointer
    ! nullify(key_object%bad_keyword_list)
    i = splice(key_object%bad_keyword_list, 0)

    ! Set the initialised flag to true
    if (key_object%initialised) then
       if (.not. destroy(key_object%keyword_list)) &
            stop "KEYWORD_PARSER :: Error! Could not free keyword list hash"
       if (.not. destroy(key_object%value_present)) &
            stop "KEYWORD_PARSER :: Error! Could not free value present hash"
       if (.not. destroy(key_object%valid_keywords)) &
            stop "KEYWORD_PARSER :: Error! Could not free valid keywords hash"
    end if

    hashsize = log2(size(keywords)-1) + 1
    
    call new(key_object%keyword_list, hashsize)
    call new(key_object%value_present, hashsize)
    call new(key_object%valid_keywords, hashsize)

    do i = 1, size(keywords)
       if (.NOT. add(key_object%valid_keywords,(.ucase.(trim(char(keywords(i))))),1)) &
       stop "KEYWORD_PARSER :: Error! Could not add keyword to table of valid keywords"
    end do

    ! Set the initialised flag to true
    key_object%initialised = .true.
    
  end subroutine initialise_keyword_object

  logical function keyword_exists (key_object, keyword)

    ! See if a keyword was specified

    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    
    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    keyword_exists = exists(key_object%keyword_list, .ucase. trim(keyword))
    
  end function keyword_exists

  logical function varstr_keyword_exists (key_object, keyword)

    ! Wrapper to exists for varying string

    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword

    varstr_keyword_exists = exists(key_object, char(keyword))

  end function varstr_keyword_exists
    
  logical function keyword_has_value (key_object, keyword)

    ! See if a command line keyword also specified a value, e.g. --max=5.4

    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    
    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    keyword_has_value = exists(key_object%value_present, .ucase. trim(keyword))
    
  end function keyword_has_value

  logical function varstr_keyword_has_value (key_object, keyword)

    ! Wrapper to has_value for varying string

    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword

    varstr_keyword_has_value = has_value(key_object, char(keyword))

  end function varstr_keyword_has_value
    
  integer function get_num_value (key_object, keyword)

    ! See if a command line keyword also specified a value, e.g. --max=5.4

    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword

    type (varying_string) :: buf
    
    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    if (has_value(key_object, keyword)) then
       get_num_value = hash_value(key_object%value_present)
    else
       get_num_value = 0
    end if
    
  end function get_num_value

  integer function get_num_value_varstr (key_object, keyword)

    ! See if a command line keyword also specified a value, e.g. --max=5.4

    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword

    get_num_value_varstr = num_value(key_object, char(keyword))

  end function get_num_value_varstr

  logical function delete_keyword (key_object, keyword) result(deleted_ok)

    ! See if a keyword was specified

    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword

    deleted_ok = .false.
    
    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    if (exists(key_object, keyword)) then
       if (debug) print *,'Deleting ',trim(keyword)
       if (has_value(key_object, keyword)) then
          if (debug) print *,'Deleting value'
          if (.not. delete(key_object%value_present)) return
          ! print *,join(keys(key_object%value_present, " "), ", ")
       end if
       if (.not. delete(key_object%keyword_list)) return
       deleted_ok = .true.
    end if
    
  end function delete_keyword
  
  logical function delete_varstr_keyword (key_object, keyword) result(deleted_ok)

    ! Varying string wrapper for delete

    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword

    deleted_ok = delete (key_object, char(keyword))

  end function delete_varstr_keyword

  function get_keyword_value (key_object, keyword) result(value)

    ! Return the value specified for a command line keyword

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    
    ! Return type of function
    type (keyword_value) :: value

    ! Local variables
    ! character(len=len(keyword)) :: kwd
    integer :: num
    
    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    ! kwd = .ucase. trim(keyword)
    ! kwd = keyword

    value%value = ''
    value%n = 0

    ! if (exists(key_object%keyword_list, kwd)) then
    if (exists(key_object, keyword)) then
       ! if (exists(key_object%value_present, kwd)) then
       if (has_value(key_object, keyword)) then
          value%value = hash_value(key_object%keyword_list)
          value%n = num_value(key_object, keyword)
       else
          print *,'value not found for keyword: ',keyword
          stop
       end if
    else
       print *,'kwd not found: ',keyword
    end if

  end function get_keyword_value

  function get_keyword_value_varstr (key_object, keyword) result(value)

    ! Return the value specified for a command line keyword. Wrapper
    ! to allow variable string access

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword
    
    ! Return type of function
    type (keyword_value) :: value

    value = get_value(key_object, char(keyword))
    
  end function get_keyword_value_varstr

  function get_one_keyword_value (key_object, keyword, valnum) result(value)

    ! Return the value specified for a command line keyword

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    integer, intent(in)                  :: valnum
    
    ! Return type of function
    type (keyword_value) :: value

    if (.not. key_object%initialised) then 
       write(stderr,*) 'KEYWORD_PARSER :: ERROR! Key object not initialised!'
       return
    end if

    ! Grab the entire keyvalue and extract the one we want
    value = extract_keyvalue(get_keyword_value(key_object, keyword), valnum)

  end function get_one_keyword_value

  function get_one_keyword_value_varstr (key_object, keyword, valnum) result(value)

    ! Return the value specified for a command line keyword

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: keyword
    integer, intent(in)                  :: valnum
    
    ! Return type of function
    type (keyword_value) :: value

    ! Grab the entire keyvalue and extract the one we want
    value = extract_keyvalue(get_keyword_value(key_object, char(keyword)), valnum)

  end function get_one_keyword_value_varstr

  subroutine print_keyword_object (key_object, unit)

    type (keyword_object), intent(inout) :: key_object
    integer, intent(in), optional        :: unit

    ! Local variables
    integer :: lunit, i, n
    type (varying_string) :: kwd, val
    
    if (present(unit)) then
       lunit = unit
    else
       lunit = stdout
    end if

    if (first(key_object%keyword_list)) then
       if (debug) print *,'Printing out keyword object'
       do
          kwd = key(key_object%keyword_list," ")
          n = num_value(key_object, char(kwd))
          if (debug) print *,'kwd = ',char(kwd),' n = ',n
          if (n < 1) then
             write(lunit,'(A)') char(kwd)
          else
             do i = 1, n
                if (debug) print *,'i = ',i
                val = get_value(key_object, char(kwd), i)
                write(lunit,'(A)') char(kwd)//" "//char(val)
             end do
          end if
          if (.not. next(key_object%keyword_list)) exit
       end do
    end if

  end subroutine print_keyword_object

  subroutine open_then_print(key_object, file)

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: file

    ! Internal variables
    logical :: opened_ok
    integer :: unit

    ! Open the file
    unit = open(file, opened_ok, status='unknown')

    if (opened_ok) then
       call print_keyword_object(key_object, unit)
       close(unit)
    else
       write (stderr,*)'KEYWORD_CLASS :: Unable to open file: ',file
    end if

  end subroutine open_then_print

  subroutine open_then_print_varstr(key_object, file)

    ! Wrapper to open_then_print with varying string for file name

    ! Input arguments
    type (keyword_object), intent(inout) :: key_object
    type (varying_string), intent(in)    :: file

    call print(key_object, char(file))

  end subroutine open_then_print_varstr

  ! The magic in the keyword retrieval part of this module is done with
  ! overloaded assignment. The get_keyword_value routine returns an 
  ! keyword value type, which we can then coerce into a number of
  ! different types, depending on what is on the left hand side of the
  ! assignment.

  subroutine assign_keyword_to_int (intval, value)

    integer, intent(out)             :: intval
    type (keyword_value), intent(in) :: value

    ! Local variables
    type (keyword_value) :: valbuf

    valbuf = extract_keyvalue(value, 1)

    ! The keyword value is stored as a varying string. We use the 
    ! overloaded assignment from string functions to convert this
    ! to an integer
    intval = valbuf%value
    
  end subroutine assign_keyword_to_int

  subroutine assign_keyword_to_intarray (intval, value)

    type (keyword_value), intent(in) :: value
    integer, dimension(value%n), intent(out) :: intval

    intval = split_keyvalue(value)
    
  end subroutine assign_keyword_to_intarray

  subroutine assign_keyword_to_real (realval, value)

    real, intent(out)            :: realval
    type (keyword_value), intent(in) :: value

    ! Local variables
    type (keyword_value) :: valbuf

    valbuf = extract_keyvalue(value, 1)

    realval = valbuf%value
    
  end subroutine assign_keyword_to_real

  subroutine assign_keyword_to_realarray (realval, value)

    type (keyword_value), intent(in) :: value
    real, dimension(value%n), intent(out) :: realval
    
    realval = split_keyvalue(value)
    
  end subroutine assign_keyword_to_realarray

  subroutine assign_keyword_to_double (realval, value)

    real(kind=rd_kind), intent(out) :: realval
    type (keyword_value), intent(in) :: value

    ! Local variables
    type (keyword_value) :: valbuf

    valbuf = extract_keyvalue(value, 1)

    realval = valbuf%value
    
  end subroutine assign_keyword_to_double

  subroutine assign_keyword_to_doublearray (realval, value)

    type (keyword_value), intent(in) :: value
    real(kind=rd_kind), dimension(value%n), intent(out) :: realval
    
    realval = split_keyvalue(value)
    
  end subroutine assign_keyword_to_doublearray

  subroutine assign_keyword_to_string (stringval, value)

    type (keyword_value), intent(in) :: value
    ! character(len=len(value%value)), intent(out) :: stringval
    character(len=*), intent(out) :: stringval
    
    ! Local variables
    type (keyword_value) :: valbuf

    valbuf = extract_keyvalue(value, 1)

    stringval = valbuf%value
    
  end subroutine assign_keyword_to_string

  subroutine assign_keyword_to_varstring (stringval, value)

    type(varying_string), intent(out) :: stringval
    type (keyword_value), intent(in) :: value
    
    ! Local variables
    type (keyword_value) :: valbuf

    valbuf = extract_keyvalue(value, 1)

    stringval = valbuf%value
    
  end subroutine assign_keyword_to_varstring

  subroutine assign_keyword_to_vstring_array (stringval, value)

    type (keyword_value), intent(in) :: value
    type (varying_string), dimension(value%n), intent(out) :: stringval
    
    stringval = split_keyvalue(value)
    
  end subroutine assign_keyword_to_vstring_array

  function key_value_as_int (key_object, keyword, mold) result(keyvalue)

    ! Return value of keyword as an integer

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    integer, intent(in)                  :: mold

    integer :: keyvalue

    ! Use default assignment to coerce the output of value to integer 
    keyvalue = get_value(key_object, keyword)

  end function key_value_as_int
  
  function key_value_as_vstring (key_object, keyword, mold) result(keyvalue)

    ! Return value of keyword as an integer

    ! Interface variables
    type (keyword_object), intent(inout) :: key_object
    character(len=*), intent(in)         :: keyword
    type (varying_string), intent(in)    :: mold

    type (varying_string) :: keyvalue

    ! Use default assignment to coerce the output of value to integer 
    keyvalue = get_value(key_object, keyword)

  end function key_value_as_vstring
  
  function split_keyvalue (value) result(valarray)

    type (keyword_value), intent(in)          :: value
    type (varying_string), dimension(value%n) :: valarray
    
    valarray = split(char(value%value),delimiter)

  end function split_keyvalue
  
  function extract_keyvalue (value, valnum) result(retvalue)

    type (keyword_value), intent(in) :: value
    integer, intent(in)              :: valnum
    type (keyword_value)             :: retvalue
    
    type (varying_string), pointer :: array(:)
    integer :: i

    retvalue = value

    ! We don't have any values, or only one value (so no need to extract anything)
    if (value%n <= 1) return

    if (valnum > value%n) then 
       write(stderr,'(A,I0,A)') 'KEYWORD_PARSER :: ERROR! Key value does not have ',valnum,' values!'
       return
    end if

    nullify(array)

    ! Make an array out of the values
    i = push(array,split_keyvalue(value))

    ! Pick out the value we want and set the number of values to be 1
    retvalue%value = array(valnum)
    retvalue%n = 1
    
    i = splice(array,0)
    
  end function extract_keyvalue

  function log2(n) result(r)

    integer, intent(in) :: n

    integer :: r

    r = 0

    do while (ishft(n,-r) /= 0)
       r = r + 1
    end do

    r = r-1; ! returns -1 for n==0, floor(log2(n)) otherwise

  end function log2
  
end module keyword_class
