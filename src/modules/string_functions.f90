module string_functions

  ! Some of this nicked from Jay William Ponder's "Tinker" and updated from a usenet
  ! post: http://groups.google.com.au/groups?selm=87wvyk9y1h.fsf%40bglbv.my-dejanews.com

  use iso_varying_string, only: varying_string, len, char, operator(//), var_str, assignment(=)
  use variable_array, only: splice
  use precision, only: rd_kind

  implicit none

  private

  character(len=*), parameter :: version = "$Id: string_functions.f90,v 1.16 2006/05/30 02:21:13 aidan Exp $"

  !! $Log: string_functions.f90,v $
  !! Revision 1.16  2006/05/30 02:21:13  aidan
  !! Changed to splice operation (from push) in expand_range routine. This
  !! is the obvious behaviour, to replace the contents of the array we
  !! return with the expanded range. The previous behaviour pretty much
  !! constitutes an unknown bug.
  !!
  !! Revision 1.15  2006/05/30 01:59:26  aidan
  !! Added a check for a range in expand_range, so that a single integer is
  !! a valid range. This makes it easier to use this function generically
  !! without having to check for a range expression, which is best done
  !! internally for consistency.
  !!
  !! Revision 1.14  2005/06/28 03:40:52  aidan
  !! Added support for conversion of strings to double precision reals. This
  !! includes a new function called 'double', with the same functionality as
  !! the overloaded 'real'.
  !!
  !! Revision 1.13  2005/05/05 02:12:46  aidan
  !! Added a new 'expand' routine. It will expand ranges like this '2-5' to
  !! an integer array like so (/ 2, 3, 4, 5 /). It handles negative increments
  !! as well, i.e. 5-2 -> (/ 5, 4, 3, 2 /). At this point in time only integer
  !! ranges can be expanded like this.
  !!
  !! Revision 1.12  2005/01/20 00:34:30  aidan
  !! Removed adjustc so that this module will compile with ifort v8.1
  !!
  !! Revision 1.11  2004/05/10 06:23:48  aidan
  !! Oops! Forgot to change the next entry in the character maps from 92
  !! to 91 after previous change. Version 1.10 does not compile. Fixed.
  !!
  !! Revision 1.10  2004/05/10 06:20:01  aidan
  !! Corrected an error with the character_maps. The uppercase alpha characters
  !! end at position 90 (0x5A) not 91!
  !!
  !! Revision 1.9  2003/11/26 05:42:38  aidan
  !! Mistake -- forgot to add the varying string split function to the interface.
  !!
  !! Revision 1.8  2003/11/26 05:29:40  aidan
  !! Added support for splitting and counting characters in varying strings.
  !!
  !! Revision 1.7  2003/11/26 04:20:48  aidan
  !! Added support for int and real to accept varying strings as input. This is
  !! a good idea, but was brought on by changes to the return value of the
  !! next_arg() function in the cmdline_arguments module -- this is often then
  !! coerced into an integer or real. Hopefully this update will stop the other
  !! update from breaking too much code.
  !!
  !! Revision 1.6  2003/10/14 05:02:21  aidan
  !! Added a join function, similar to the perl function. This is the opposite of
  !! split, and joins an array of character strings into one string with an
  !! optionally specified delimiter separating the strings.
  !!
  !! Split now uses a varying string as the internal buffer.
  !!
  !! Revision 1.5  2003/10/08 03:00:21  aidan
  !! Added overloaded assignment operator to convert strings and varying_strings
  !! to integers and reals (including rank one arrays).
  !!
  !! Revision 1.4  2003/10/01 06:22:20  aidan
  !! Changed 'split' to return array of varying strings.
  !!
  !! Revision 1.3  2003/09/03 06:14:03  aidan
  !! Made int and real functions elemental--now they can operate on whole
  !! arrays. Added a character count function and a split function. The latter
  !! splits strings given a delimiter character. Crude but works ok.
  !!
  !! Revision 1.2  2003/08/19 04:42:40  aidan
  !! Case changing functions are now called with an operator interface. e.g.
  !! .ucase. rather than ucase(). Added subroutine versions of the case changing
  !! functions, e.g. call ucase(string), which can also handle rank 1 arrays
  !! of strings. Did this because the previous ucase function could not be
  !! made elemental, so we couldn't take array arguments. Solved this in a
  !! slightly ugly way.
  !!
  !! Revision 1.1  2003/08/19 04:03:08  aidan
  !! Initial revision
  !!

  integer :: i

  type character_map
     private
     ! 8 bit character maps only ...
     character(len=1) :: translation (0:255)
  end type character_map

  ! type(character_map), parameter :: ascii_upper_case_map = &
  !      character_map ( (/ (char(i),i=0,96), (char(i-32),i=97,122), &
  !               (char(i),i=123,223), (char(i-32),i=224,246), char(247), &
  !               (char(i-32),i=248,255) /) )

  type(character_map), parameter :: ascii_upper_case_map = &
       character_map ( (/ (char(i),i=0,96), (char(i-32),i=97,122), &
                (char(i),i=123,127), (char(0),i=128,255) /) )

  type(character_map), parameter :: ascii_lower_case_map = &
       character_map ( (/ (char(i),i=0,64), (char(i+32),i=65,90), &
                (char(i),i=91,127), (char(0),i=128,255) /) )

  type(character_map), parameter :: ascii_toggle_case_map = &
       character_map ( (/ (char(i),i=0,64), (char(i+32),i=65,90), (char(i),i=91,96), & 
                (char(i-32),i=97,122), (char(i),i=123,127), (char(0),i=128,255) /) )

  ! This will only map 7 bit characters back on to themselves
  type(character_map), parameter :: ascii_map = character_map ( &
       (/ (char(i),i=0,127), (char(0),i=128,255) /) )
  
  ! type(character_map), parameter, public :: identity_map = character_map ( &
  !      (/ (char(i),i=0,255) /) )
   
  ! All alpha characters + space
  type(character_map), parameter :: alpha_map = &
       character_map ( (/ (char(0),i=0,31), char(32), (char(0),i=33,64), (char(i),i=65,90), &
                (char(0),i=91,96), (char(i),i=97,122), (char(0),i=123,255) /) )
   
  ! All numeric characters + space
  type(character_map), parameter :: numeric_map = &
       character_map ( (/ (char(0),i=0,31), char(32), (char(0),i=33,47), (char(i),i=48,57), &
                (char(0),i=58,255) /) )

  character(len=1), parameter :: blank = ' '

  character(len=26), parameter :: alpha_lower = 'abcdefghijklmnopqrstuvwxyz', &
       alpha_upper = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

  character(len=10), parameter :: numerals = '0123456789'

  ! Convert character strings to integer -- currently not all that useful
  ! as this is not allowed to be part of a write statement -- get a
  ! recursive IO error. Consequently this interface has been removed
  ! interface string
  !    module procedure int_to_str, real_to_str
  ! end interface

  ! Convert character strings to integer
  interface int
     module procedure str_to_int, varstr_to_int
  end interface

  ! Convert character strings to real
  interface real
     module procedure str_to_real, varstr_to_real
  end interface

  ! Convert character strings to real
  interface double
     module procedure str_to_double, varstr_to_double
  end interface

  ! Provide interface to routines which make words uppercase, 
  ! lowercase and toggle (change) case
  interface operator(.ucase.)
     module procedure upper_case
  end interface
  interface operator(.lcase.)
     module procedure lower_case
  end interface
  interface operator(.ccase.)
     module procedure toggle_case
  end interface

  ! Provide interface to routines which make words uppercase, 
  ! lowercase and toggle (change) case
  interface ucase
     module procedure upper_case_array_sub, upper_case_sub
  end interface
  interface lcase
     module procedure lower_case_array_sub, lower_case_sub
  end interface
  interface ccase
     module procedure toggle_case_array_sub, toggle_case_sub
  end interface

  ! Provides a series of operators to determine the composition of
  ! a string

  ! Two complementary operators which test if all alpha characters
  ! in a string are lower or upper case respectively
  interface operator(.islcase.)
     module procedure is_lower_case
  end interface
  interface operator(.isucase.)
     module procedure is_upper_case
  end interface

  ! Two complementary operators which test if all characters
  ! in a string (with leading and trailing whitespace removed)
  ! are alpha or numeric respectively
  interface operator(.isletter.)
     module procedure is_alpha
  end interface
  interface operator(.isnumber.)
     module procedure is_numeric
  end interface

  ! Split a string into array of strings given a delimiter
  interface split
     module procedure split_string ! , split_varstring
  end interface

  ! Counts the number of times a character occurs in a string
  interface count
     module procedure count_character, count_character_vstr
  end interface

  ! Joins array of strings
  interface join
     module procedure join_string, join_string_nodelim, join_varstring, join_varstring_nodelim
  end interface

  ! Expand ranges of integers, i.e. 1-5 is expanded to 1,2,3,4,5
  interface expand
     module procedure expand_range
  end interface

  ! Provides for automatic translation from string and varying_string
  ! to integer and real
  interface assignment(=)
     module procedure assign_str_to_int, assign_str_to_int_array, &
          assign_varstr_to_int, assign_varstr_to_int_array, &
          assign_str_to_real, assign_str_to_real_array, &
          assign_varstr_to_real, assign_varstr_to_real_array, &
          assign_str_to_double, assign_str_to_double_array, &
          assign_varstr_to_double, assign_varstr_to_double_array
          
  end interface

  ! Public variables
  public :: blank, character_map

  ! Public routines
  public :: gettext, getword, ucase, lcase, ccase, is_in_map, split, join !, adjustc
  public :: expand

  ! Overloaded routines
  public :: real, double, int, count

  ! Public operators
  public :: operator(.ucase.), operator(.lcase.), operator(.ccase.)
  public :: operator(.islcase.), operator(.isucase.), operator(.isletter.), operator(.isnumber.)
  public :: assignment(=)

contains

  pure function join_string_nodelim (string) result(joined)
    
    ! Input variable
    character(len=*), dimension(:), intent(in) :: string

    ! Concatenated output -- this is the sum of the length
    ! the all the strings in the input array
    character(len=size(string)*len(string)) :: joined

    joined = join(string,'')

  end function join_string_nodelim

  pure function join_string (string, delimiter) result(joined)
    
    ! Input variable
    character(len=*), dimension(:), intent(in) :: string
    character(len=*), intent(in)               :: delimiter

    ! Concatenated output -- this is the sum of the length
    ! the all the strings in the input array plus the delimiters
    ! we will separate them with
    character(len=size(string)*len(string)+(size(string)-1)*len(delimiter)) :: joined

    ! Local variables
    type (varying_string) :: local_varstring
    integer :: i, lower

    lower = lbound(string,1)
    local_varstring = var_str(string(lower))

    do i = lower + 1, ubound(string,1)
       local_varstring = local_varstring // delimiter // string(i)
    end do

    joined = local_varstring

    ! Reset the local variable so we don't have a memory leak
    local_varstring = ''
    
  end function join_string

  pure function join_varstring_nodelim (varstring) result(joined)
    
    ! Input variable
    type (varying_string), dimension(:), intent(in) :: varstring

    ! Concatenated output -- this is the sum of the length
    ! the all the strings in the input array
    character(len=sum(len(varstring))) :: joined

    joined = join(varstring,'')

  end function join_varstring_nodelim

  pure function join_varstring (varstring, delimiter) result(joined)
    
    ! Input variable
    type (varying_string), dimension(:), intent(in) :: varstring
    character(len=*), intent(in)                    :: delimiter

    ! Concatenated output -- this is the sum of the length
    ! the all the strings in the input array plus the delimiters
    ! we will separate them with
    character(len=sum(len(varstring))+(size(varstring)-1)*len(delimiter)) :: joined

    ! Local variables
    type (varying_string) :: local_varstring
    integer :: i, lower

    lower = lbound(varstring,1)
    local_varstring = char(varstring(lower))

    do i = lower + 1, ubound(varstring,1)
       local_varstring = local_varstring // delimiter // varstring(i)
    end do

    joined = local_varstring

    ! Reset the local variable so we don't have a memory leak
    local_varstring = ''
    
  end function join_varstring
   
   
  elemental function translation (string, map)
    
    ! Translate string with map

    character(len=*), intent(in) :: string
    type(character_map), intent(in) :: map
    character(len=len(string)) :: translation
    forall (i=1:len(string))
       translation(i:i) = map%translation(ichar(string(i:i)))
    end forall

  end function translation

  elemental function is_in_map (string, map)

    ! A string is in a map if it's translation is the same as itself

    character(len=*), intent(in) :: string
    type(character_map), intent(in) :: map
    logical :: is_in_map
    is_in_map = (translation(string,map)==string)

  end function is_in_map

  subroutine set (m, from, to)
    type(character_map), intent(in out) :: m
    character(len=*), intent(in) :: from, to ! No duplicates allowed in FROM
    if (len(to) < len(from)) stop 'error'    ! Exceptions would be nice...
    forall (i=1:len(from))
       m%translation(ichar(from(i:i))) = to(i:i)
    end forall
  end subroutine set

!!$  elemental function adjustc (s)
!!$
!!$    ! Mocks up an adjust-centre function (similar to adjustl and adjustr)
!!$
!!$    character(len=*), intent(in) :: s
!!$    character(len=len(s)) :: adjustc
!!$    
!!$    integer :: ileft, iright, imbalance
!!$    
!!$    ileft = verify (s, ' ', back=.false.)
!!$    if (ileft == 0) then
!!$       adjustc = s
!!$    else
!!$       iright = len (s) - verify (s(ileft:), ' ', back=.true.) - (ileft-1)
!!$       ileft  = ileft - 1
!!$       imbalance = (iright-ileft)/2
!!$       if (imbalance > 0) then
!!$          adjustc = s(len(s)+1-imbalance:) // s
!!$       else
!!$          adjustc = s(1-imbalance:)
!!$       end iF
!!$    end if
!!$  end function adjustc

  Function int_to_str( input )

    ! Return a string given an integer

    character(len=30) :: int_to_str
    integer, intent(in) :: input

    character(len=30) :: tmp

    ! write(int_to_str,'(I)') input
    ! write(int_to_str,'(I5)') input
    write(tmp,*) input
    ! print *,int_to_str
    int_to_str = tmp

  end function int_to_str

  function real_to_str( input )

    ! Return a string given an integer

    character(len=30) :: real_to_str
    real, intent(in) :: input

    write(real_to_str,*) input

  end function real_to_str

  elemental integer function str_to_int( str )

    ! Return an integer given an input string (of integer presumably)

    character(len=*), intent(in) :: str

    ! It would be nice to check the iostatus of this read but
    ! we can't do anything about it as we want this to be an
    ! elemental function which means we can't stop or print
    ! inside here
    read(str,*) str_to_int

  end function str_to_int

  elemental integer function varstr_to_int( varstr )

    ! Return an integer given an input string (of integer presumably)

    type (varying_string), intent(in) :: varstr

    varstr_to_int = int(char(varstr))

  end function varstr_to_int

  ! The following are some assignment rules written using the function
  ! above

  subroutine assign_str_to_int( intval, str )

    ! Assign a string to an integer

    integer, intent(out)         :: intval
    character(len=*), intent(in) :: str
    
    intval = int(str)

  end subroutine assign_str_to_int

  subroutine assign_str_to_int_array( intval, str )

    ! Assign a string array to an integer array

    integer, dimension(:), intent(out)         :: intval
    character(len=*), dimension(:), intent(in) :: str
    
    intval = int(str)

  end subroutine assign_str_to_int_array

  subroutine assign_varstr_to_int( intval, varstr )

    ! Assign a varying string to an integer

    integer, intent(out)             :: intval
    type(varying_string), intent(in) :: varstr
    
    intval = int(char(varstr))

  end subroutine assign_varstr_to_int

  subroutine assign_varstr_to_int_array( intval, varstr )

    ! Assign a varying string array to an integer array

    integer, dimension(:), intent(out)             :: intval
    type(varying_string), dimension(:), intent(in) :: varstr
    
    intval = int(char(varstr))

  end subroutine assign_varstr_to_int_array

  elemental real function str_to_real( str )

    ! Return a real given an input string (of real presumably)

    character(len=*), intent(in) :: str

    ! It would be nice to check the iostatus of this read but
    ! we can't do anything about it as we want this to be an
    ! elemental function which means we can't stop or print
    ! inside here
    read(str,*) str_to_real

  end function str_to_real

  elemental real function varstr_to_real( varstr )

    ! Return a real given an input string (of real presumably)

    type (varying_string), intent(in) :: varstr

    varstr_to_real = real(char(varstr))

  end function varstr_to_real

  elemental real(kind=rd_kind) function str_to_double( str )

    ! Return a real given an input string (of real presumably)

    character(len=*), intent(in) :: str

    ! It would be nice to check the iostatus of this read but
    ! we can't do anything about it as we want this to be an
    ! elemental function which means we can't stop or print
    ! inside here
    read(str,*) str_to_double

  end function str_to_double

  elemental real(kind=rd_kind) function varstr_to_double( varstr )

    ! Return a real given an input string (of real presumably)

    type (varying_string), intent(in) :: varstr

    varstr_to_double = real(char(varstr))

  end function varstr_to_double

  ! The following are some assignment rules written using the function
  ! above

  subroutine assign_str_to_real( realval, str )

    ! Assign a string to a real

    real, intent(out)            :: realval
    character(len=*), intent(in) :: str
    
    realval = real(str)

  end subroutine assign_str_to_real

  subroutine assign_str_to_real_array( realval, str )

    ! Assign a string array to a real array

    real, dimension(:), intent(out)            :: realval
    character(len=*), dimension(:), intent(in) :: str
    
    realval = real(str)

  end subroutine assign_str_to_real_array

  subroutine assign_varstr_to_real( realval, varstr )

    ! Assign a string to a real

    real, intent(out)                :: realval
    type(varying_string), intent(in) :: varstr
    
    realval = real(char(varstr))

  end subroutine assign_varstr_to_real

  subroutine assign_varstr_to_real_array( realval, varstr )

    ! Assign a string array to a real array

    real, dimension(:), intent(out)                :: realval
    type(varying_string), dimension(:), intent(in) :: varstr
    
    realval = real(char(varstr))

  end subroutine assign_varstr_to_real_array

  ! And the same for double precision reals
  
  subroutine assign_str_to_double( doubleval, str )

    ! Assign a string to a double precision real

    real(kind=rd_kind), intent(out) :: doubleval
    character(len=*), intent(in)    :: str
    
    doubleval = double(str)

  end subroutine assign_str_to_double

  subroutine assign_str_to_double_array( doubleval, str )

    ! Assign a string array to a double precision array

    real(kind=rd_kind), dimension(:), intent(out) :: doubleval
    character(len=*), dimension(:), intent(in)    :: str
    
    doubleval = double(str)

  end subroutine assign_str_to_double_array

  subroutine assign_varstr_to_double( doubleval, varstr )

    ! Assign a varying string to a double precision real

    real(kind=rd_kind), intent(out)  :: doubleval
    type(varying_string), intent(in) :: varstr
    
    doubleval = double(char(varstr))

  end subroutine assign_varstr_to_double

  subroutine assign_varstr_to_double_array( doubleval, varstr )

    ! Assign a varying string array to a double precision array

    real(kind=rd_kind), dimension(:), intent(out)  :: doubleval
    type(varying_string), dimension(:), intent(in) :: varstr
    
    doubleval = double(char(varstr))

  end subroutine assign_varstr_to_double_array

  function split_string (string, delimiter)

    ! Split a string into an array of strings based on a supplied
    ! delimiter (single character variable). We don't do anything
    ! fancy with quotes -- just split on a specified character

    ! Input variables
    character(len=*), intent(in) :: string
    character(len=1), intent(in) :: delimiter

    ! Return type of function -- for convenience we return an
    ! array of strings the same length as the input string
    ! with the dimension determined by counting the number of
    ! delimiter characters in the string
    ! character(len=len(string)), dimension(count_character(string,delimiter)+1) :: split_string
    type(varying_string), dimension(count(string,delimiter)+1) :: split_string

    ! Local variables. Buffer has the save attribute so that we
    ! don't have a memory leak -- uses the same varying string
    ! every time
    type (varying_string), save :: buffer
    integer :: split_index, i

    split_index = 1
    buffer = ''

    ! Cycle through the string
    do i = 1, len(string)
       if (string(i:i) == delimiter) then
          split_string(split_index) = buffer
          buffer = ''
          split_index = split_index + 1
       else
          buffer = buffer // string(i:i)
       end if
    end do
    split_string(split_index) = buffer

  end function split_string

!  function split_varstring (varstring, delimiter)
!
!    ! Split a string into an array of strings based on a supplied
!    ! delimiter (single character variable). We don't do anything
!    ! fancy with quotes -- just split on a specified character
!
!    ! Input variables
!    type (varying_string), intent(in) :: varstring
!    character(len=1), intent(in) :: delimiter
!
!    ! Return type of function -- for convenience we return an
!    ! array of strings the same length as the input string
!    ! with the dimension determined by counting the number of
!    ! delimiter characters in the string
!    type(varying_string), dimension(count(varstring,delimiter)+1) :: split_varstring
!    character(len=len(varstring)) :: charstring
!
!    charstring = char(varstring)
!
!    split_varstring = split(charstring,delimiter)
!
!  end function split_varstring

  pure integer function count_character (string, character)

    ! Count the number of occurences of a particular character in
    ! a given string

    character(len=*), intent(in) :: string
    character(len=1), intent(in) :: character

    character(len=1), dimension(len(string)) :: char_array

    char_array = transfer(string,char_array)

    count_character = count(char_array == character)

  end function count_character

  pure integer function count_character_vstr (varstring, character)

    ! Count the number of occurences of a particular character in
    ! a given string

    type (varying_string), intent(in) :: varstring
    character(len=1), intent(in) :: character

    character(len=1), dimension(len(varstring)) :: char_array

    char_array = transfer(char(varstring),char_array)

    count_character_vstr = count(char_array == character)

  end function count_character_vstr


  function gettext (string,next)

    ! "gettext" searchs an input string for the first string of
    ! non-blank characters; the region from a non-blank character
    ! to the first blank space is returned
    
    ! variables and parameters:
    
    !     string    input character string to be searched
    !     next      input with first position of search string;
    !                 output with the position following text

    character(len=*), intent(in)  :: string
    integer, intent(inout)        :: next

    ! Note that the length of gettext is determined by the length of 
    ! the input string and the starting point within it 
    ! character(len=len(string)-next) :: gettext
    ! character(len=len(string)) :: gettext
    character(len=2000) :: gettext
    
    integer :: i, j, first, last, initial, final

    !  move through the string one character at a time,
    !  searching for the first non-blank character
    first = next
    last = 0
    initial = next
    final = next + len(string(next:)) - 1
    do i = initial, final
       if (string(i:i) .gt. blank) then
          first = i
          last = i
          do j = i+1, final
             if (string(j:j) .le. blank) then
                last = j - 1
                next = j
                goto 10
             end if
          end do
       end if
    end do
10  continue
    
    ! transfer the text into the return string
    gettext = string(first:last)

  end function gettext
    
  function getword (string,next)

    ! "getword" searchs an input string for the first alphabetic
    ! character (A-Z or a-z); the region from this first character
    ! to the first blank space or comma is returned
    !
    ! variables and parameters:
    !
    ! string    input character string to be searched
    ! next      input with first position of search string;
    !           output with the position following word

    character(len=*), intent(in)  :: string
    integer, intent(inout)        :: next

    ! Note that the length of getword is determined by the length of 
    ! the input string and the starting point within it 
    character(len=len(string)-next) :: getword
    
    integer   :: i, j, first, last, initial, final
    character :: letter
    
    ! move through the string one character at a time,
    ! searching for the first alphabetic character
    
    first = next
    last = 0
    initial = next
    final = next + len(string(next:)) - 1
    do i = initial, final
       letter = string(i:i)
       if ((letter.ge.'A' .and. letter.le.'Z') .or.                   &
            &       (letter.ge.'a' .and. letter.le.'z')) then
          first = i
          last = i
          do j = i+1, final
             if (string(j:j).le.' ' .or. string(j:j).eq.',') then
                last = j - 1
                next = j
                goto 10
             end if
          end do
       end if
    end do
10  continue
    
    ! transfer the word into the return string
    getword = string(first:last)

    ! skip over the next character when it is a comma
    if (string(next:next) .eq. ',')  next = next + 1

  end function getword
  
  pure function upper_case (string)

    ! Function to return an uppercase version of an input string
    
    character(len=*), intent(in) :: string
    character(len=len(string)) :: upper_case

    upper_case = translation (string, ascii_upper_case_map)
    
  end function upper_case
  
  subroutine upper_case_sub (string)
    
    ! Subroutine wrapper to uppercase function

    character(len=*), intent(inout) :: string
    
    string = upper_case(string)
    
  end subroutine upper_case_sub
  
  subroutine upper_case_array_sub (string)
    
    ! Subroutine wrapper to uppercase function -- allows for arrays
    ! of input strings

    character(len=*), dimension(:), intent(inout) :: string
    
    forall(i=1:size(string))
       string(i) = upper_case (string(i))
    end forall
    
  end subroutine upper_case_array_sub
    
  pure function lower_case(string)

    ! Function to return a lowercase version of an input string
    
    character(len=*), intent(in) :: string
    character(len=len(string)) :: lower_case
    lower_case = translation (string, ascii_lower_case_map)
    
  end function lower_case
  
  subroutine lower_case_sub (string)
    
    ! Subroutine wrapper to lowercase function

    character(len=*), intent(inout) :: string
    
    string = lower_case (string)
    
  end subroutine lower_case_sub
  
  subroutine lower_case_array_sub (string)
    
    ! Subroutine wrapper to lowercase function -- allows for arrays
    ! of input strings

    character(len=*), dimension(:), intent(inout) :: string
    
    forall(i=1:size(string))
       string(i) = lower_case (string(i))
    end forall
    
  end subroutine lower_case_array_sub
    
  pure function toggle_case(string)

    ! Function to return a version of an input string with the case toggled
    ! i.e. lower goes to upper and vice versa

    character(len=*), intent(in) :: string
    character(len=len(string)) :: toggle_case
    toggle_case = translation (string, ascii_toggle_case_map)
    
  end function toggle_case
  
  subroutine toggle_case_sub (string)

    ! Subroutine wrapper to togglecase function
    
    character(len=*), intent(inout) :: string
    
    string = toggle_case (string)
    
  end subroutine toggle_case_sub
  
  subroutine toggle_case_array_sub (string)
    
    ! Subroutine wrapper to togglecase function -- allows for arrays
    ! of input strings

    character(len=*), dimension(:), intent(inout) :: string
    
    forall(i=1:size(string))
       string(i) = toggle_case (string(i))
    end forall
    
  end subroutine toggle_case_array_sub
    
  elemental function is_upper_case(string)

    ! This function returns true if all the characters in the string are
    ! uppercase (ignores non alpha characters)

    character(len=*), intent(in) :: string

    logical                      :: is_upper_case

    is_upper_case = .TRUE.
    if (scan(string,alpha_lower) > 0) is_upper_case = .FALSE.
    
  end function is_upper_case
  
  elemental function is_lower_case(string)

    ! This function returns true if all the characters in the string are
    ! uppercase (ignores non alpha characters)

    character(len=*), intent(in) :: string

    logical                      :: is_lower_case

    is_lower_case = .TRUE.
    if (scan(trim(adjustl(string)),alpha_upper) > 0) is_lower_case = .FALSE.
    
  end function is_lower_case
  
  elemental function is_alpha(string)

    ! This function returns true if all the characters in the string are
    ! alpha (letters)

    character(len=*), intent(in) :: string

    logical                      :: is_alpha

    is_alpha = .TRUE.
    if (verify(trim(adjustl(string)),alpha_upper//alpha_lower) > 0) is_alpha = .FALSE.
    
  end function is_alpha
  
  elemental function is_numeric(string)

    ! This function returns true if all the characters in the string are
    ! numeric (numbers)

    character(len=*), intent(in) :: string

    logical                      :: is_numeric

    is_numeric = .TRUE.
    if (verify(trim(adjustl(string)),numerals) > 0) is_numeric = .FALSE.
    
  end function is_numeric

  elemental function is_ascii (s)

    ! Returns true if all characters in the string are ascii. This
    ! function is not currently available in the public interface,
    ! but is retained to show how is_in_map can be used to verify 
    ! contents of strings -- this method is an order of magnitude
    ! slower for numeric and alpha comparisons compared to the 
    ! functions listed above

    character(len=*), intent(in) :: s
    logical :: is_asciI
    is_ascii = is_in_map (s, ascii_map)

  end function is_ascii

  subroutine expand_range(input,output)

    character(len=*), intent(in) :: input 
    integer, pointer             :: output(:)

    ! Local variables
    integer :: limits(2), junk, increment

    ! Make sure we can deal with a single number, i.e. no
    ! range at all .. this way we can call this generically
    ! without caring if it represents a range or not
    if (scan(input,"-") == 0) then
       junk = splice(output,0,size(output),(/int(input)/))
       return
    end if

    limits = int(split(input,"-"))

    increment = 1
    if (limits(1) > limits(2)) increment = -1

    junk = splice(output,0,size(output),(/ (i, i=limits(1),limits(2),increment) /))

  end subroutine expand_range

end module string_functions

!!$INDEX is not too bad on most modern hardware, which often
!!$have single or special instructions to scan a string for a single
!!$character.  Faster would be the following (not portable if your
!!$program uses characters that aren't in the ASCII set):
!!$
!!$      function upcase(c)
!!$      character*1 upcase, c
!!$      character*1 upper(128)
!!$... DATA to set UPPER to the ASCII sequence, but with
!!$... all lowercase replaced with the corresponding uppercase.
!!$      upcase = upper(iachar(c))
!!$      return
!!$      end
!!$
!!$Hypothetically, a good compiler could use a single instruction
!!$with the lookup table on many hardware platforms (XLAT on
!!$PCs, for example).
  
!!$  elemental function upper_case(string)
!!$
!!$    !    This function reads in a character string, changes the case of letters and 
!!$    !    writes out new string.
!!$    !
!!$    !  Note:
!!$    !
!!$    !    Need to know the difference in the collation sequence of the upper 
!!$    !    and lower case characters - use position of upper and lower case A: 
!!$    !
!!$    !                        iachar('A') - iachar('a')
!!$    
!!$    character(len=*), intent(in) :: string
!!$    character(len=len(string))   :: upper_case
!!$
!!$    INTEGER :: i, lower_to_upper
!!$    
!!$    lower_to_upper = iachar("A") - iachar("a")
!!$    
!!$    do i = 1, len_trim(string)
!!$       select case (string(i:i))
!!$       case ('a':'z')            ! Lower case character found:
!!$          upper_case(i:i) = achar(iachar(string(i:i)) + lower_to_upper)
!!$       case default
!!$          !  No change for any other characters
!!$          upper_case(i:i) = string(i:i)
!!$       end select
!!$    end do
!!$    
!!$  end function upper_case
!!$    
!!$  elemental function lower_case(string)
!!$
!!$    !    This function reads in a character string, changes the case of letters and 
!!$    !    writes out new string.
!!$    !
!!$    !  Note:
!!$    !
!!$    !    Need to know the difference in the collation sequence of the upper 
!!$    !    and lower case characters - use position of upper and lower case A: 
!!$    !
!!$    !                        iachar('A') - iachar('a')
!!$    
!!$    character(len=*), intent(in) :: string
!!$    character(len=len(string))   :: lower_case
!!$
!!$    INTEGER :: i, lower_to_upper
!!$    
!!$    lower_to_upper = iachar("A") - iachar("a")
!!$    
!!$    do i = 1, len_trim(string)
!!$       select case (string(i:i))
!!$       case ('A':'Z')            ! Upper case character found:
!!$          lower_case(i:i) = achar(iachar(string(i:i)) - lower_to_upper)
!!$       case default
!!$          !  No change for any other characters
!!$          lower_case(i:i) = string(i:i)
!!$       end select
!!$    end do
!!$    
!!$  end function lower_case
!!$    
!!$  function toggle_case(string)
!!$
!!$    !    This function reads in a character string, changes the case of letters and 
!!$    !    writes out new string.
!!$    !
!!$    !  Note:
!!$    !
!!$    !    Need to know the difference in the collation sequence of the upper 
!!$    !    and lower case characters - use position of upper and lower case A: 
!!$    !
!!$    !                        iachar('A') - iachar('a')
!!$    
!!$    character(len=*), intent(in) :: string
!!$    character(len=len(string))   :: toggle_case
!!$
!!$    INTEGER :: i, lower_to_upper
!!$    
!!$    lower_to_upper = iachar("A") - iachar("a")
!!$    
!!$    do i = 1, len_trim(string)
!!$       select case (string(i:i))
!!$       case ('a':'z')            ! Lower case character found:
!!$          toggle_case(i:i) = achar(iachar(string(i:i)) + lower_to_upper)
!!$       case ('A':'Z')            ! Upper case character found:
!!$          toggle_case(i:i) = achar(iachar(string(i:i)) - lower_to_upper)
!!$       case default
!!$          !  No change for any other characters
!!$          toggle_case(i:i) = string(i:i)
!!$       end select
!!$    end do
!!$    
!!$  end function toggle_case
!!$  

