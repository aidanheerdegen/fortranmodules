module variable_array

  use iso_varying_string

  implicit none

  private

  ! This module provides functionality that is similar to the perl
  ! array functions shift, unshift, pop, push and splice. These
  ! operators can be used freely on pointer integer, real, double,
  ! character and logical arrays of rank 1 to provide a dynamic array 
  ! functionality. This is particularly useful for applications where 
  ! the length of the given array is likely to grow and shrink in a 
  ! dynamic way. This module has NOT been written with performance in 
  ! mind. Complete copies are done every time an element is added or 
  ! deleted from an array. This is not a particularly expensive operation 
  ! with arrays of less than 1000 or so elements, but don't use these 
  ! operators in the middle of a Monte Carlo loop!

  ! Splice is the general function which is called by the other, specialised 
  ! functions below. This is an interface to integer, real, character and 
  ! varying string splice functions.
  interface splice
     module procedure splice_integer, splice_real, splice_character, splice_vstring, &
          splice_char_into_vstring, splice_logical
  end interface

  ! Place an element or an array of elements on to the end of our array. Returns
  ! the current length of the array.
  interface push
     module procedure push_integer, push_integer_array, push_real, push_real_array, &
          push_character, push_character_array, push_vstring, push_vstring_array, &
          push_char_onto_vstring, push_char_array_onto_vstring, push_logical, push_logical_array
  end interface

  ! Remove an element off the end of our array and return it.
  interface pop
     module procedure pop_integer, pop_real, pop_character, pop_vstring, pop_logical
  end interface

  ! Place an element or an array of elements at the beginning of our array, and shifting
  ! up the index value of all existing elements. Returns the current length of the array.
  interface unshift
     module procedure unshift_integer, unshift_integer_array, unshift_real, unshift_real_array, &
          unshift_character, unshift_character_array, unshift_vstring, unshift_vstring_array, & 
          unshift_char_onto_vstring, unshift_char_array_onto_vstring, unshift_logical, unshift_logical_array
  end interface

  ! Remove an element off the beginning of our array and return it--the index of all 
  ! other elements are shifted down.
  interface shift
     module procedure shift_integer, shift_real, shift_character, shift_vstring, shift_logical
  end interface

  ! Public routines
  public :: splice, push, pop, unshift, shift

  logical, parameter :: debug = .FALSE.

contains

  integer function splice_integer (array, offset, length, list) result(array_length)

    ! This function is meant to emulate the perl function 'splice' but on a
    ! pointer integer array. So, here is the info on the perl function
    ! with some modification:
    !
  
    !   splice ARRAY,OFFSET,LENGTH,LIST
    !   splice ARRAY,OFFSET,LENGTH
    !   splice ARRAY,OFFSET
    !   splice ARRAY
    !           Removes the elements designated by OFFSET (from the lower bound of
    !           the array) and LENGTH from an
    !           array, and replaces them with the elements of LIST, if any.  
    !           Returns the number of elements removed in the array. The array 
    !           grows or shrinks as necessary.
    !           If OFFSET is negative then it starts that far from the end of the
    !           array.  If LENGTH is omitted, removes everything from OFFSET
    !           onward.  If LENGTH is negative, leaves that many elements off the
    !           end of the array.  If both OFFSET and LENGTH are omitted, removes
    !           everything.
    !
    !           The following equivalences hold (assuming the lower bound of array 
    !           is zero):
    !
    !               push(array,x)          splice(array,size(array)+1,0,x)
    !               push(array,(/x,y/))    splice(array,size(array)+1,0,(/x,y/))
    !               pop(array)             splice(array,-1)
    !               shift(array)           splice(array,0,1)
    !               unshift(array,x)       splice(array,0,0,x)
    !               array(x) = y           splice(array,x,1,y)

    ! The main differences between the perl version and this one are the definition
    ! of offset and the return value of the function. Here offset is a length from 
    ! the lower bound of the array (in perl it is an absolute position). This is
    ! so that we can usefully use arrays which start with lower bounds less than zero,
    ! as we want a negative offset to indicate counting from the end of an array.

    ! Some useful equalities to bear in mind when reading this code:
    !
    !          length of array = last element - first element + 1
    !          nth element of array = first element + n - 1

    integer, dimension(:), pointer  :: array
    integer, intent(in), optional   :: offset, length
    integer, dimension(:), optional :: list

    integer, dimension(:), pointer :: buffer_array

    integer :: status, current_length, new_offset
    integer :: replace_length, new_length, lower, tail

    logical ::  array_associated

    array_length = 0

    if (present(offset)) then

       if (associated(array)) then
          ! Grab the number of elements in array
          current_length = size(array)
          ! And the index that the array starts at
          lower = lbound(array,1)
          array_associated = .true.
       else
          ! The array isn't associated, so we set the length to zero
          current_length = 0
          ! Assume the array index starts at one
          lower = 1
          array_associated = .false.
       end if
       if (debug) print *,'Current_length: ',current_length
       if (debug) print *,'Lower bound: ',lower

       ! We don't want to alter the return value of offset so we need to
       ! make a local version of it
       new_offset = offset

       ! A negative offset means we count from the end of the array, and
       ! -1 means the last entry of the array
       if (new_offset < 0) then
          new_offset = lower + current_length + new_offset
          if (new_offset < lower) new_offset = lower
       else
          ! Convert our offset from a relative to an absolute position
          new_offset = lower + new_offset
       end if
       if (debug) print *,'Offset: ',new_offset

       ! We don't want to alter the return value of length so we need to
       ! make a local version of it too
       if (present(length)) then
          new_length = length
          ! Tail is the number of elements left at the end of the array
          ! after offset + length elements have been removed
          tail = max(0,(lower + current_length - 1) - (new_offset + new_length - 1))
       else
          ! No length, so the new_length is to the end of the array, and
          ! there are no tail elements
          new_length = (lower + current_length - 1) - new_offset + 1
          tail = 0
       end if
       if (debug) print *,'Length: ',new_length
       if (debug) print *,'Tail: ',tail

       ! Replace length is the length of list we will put back into
       ! the array to replace the elements we are snipping out. This
       ! is initialised to zero as we are not even required to provide
       ! a replacement list
       replace_length = 0
       if (present(list)) then
          replace_length = size(list)
       end if
       if (debug) print *,'Replace length: ',replace_length

       ! Make a new array that is one less than the offset length plus
       ! the length of the list we are stuffing in plus whatever is left
       ! on the end
       ! allocate(buffer_array(lower:(lower-1)+(new_offset-1)+replace_length+tail), stat=status)
       allocate(buffer_array(lower:lower + &
            (((new_offset-1)-lower+1)+replace_length+tail) & ! length of new array
            - 1 ), stat=status)
       if (status /= 0) then
          print *,'module variable_array :: Error in allocating space for buffer array: ',status
          stop
       end if
       ! buffer_array = -9
       if (debug) print *,'Length of new array: ',size(buffer_array)

       ! Copy over everything up to (not including) offset as long as our
       ! new offset is greater than the start of the array (if not then we
       ! are probably shifting something on to the start of the array)
       if (array_associated .and. (new_offset > lower)) then
          buffer_array(lower:new_offset-1) = array(lower:new_offset-1)
          if (debug) print *,'First bit: ', array(lower:new_offset-1), buffer_array
       end if

       ! If we had a list to insert then place it now. Note that we use
       ! the replace length here ...
       if (present(list)) then
          buffer_array(new_offset:new_offset+replace_length-1) = list
          if (debug) print *,'List bit: ', buffer_array
       end if

       ! If we had elements left on the tail of array then place them on the 
       ! end. Note that we use replace length on the left and new length on
       ! the right. In this way we specify how much space our inserted text
       ! takes -- which can be as low as zero, which implies an insertion and
       ! not a deletion operation
       if (tail > 0) then
          buffer_array(new_offset+replace_length:) = array(new_offset+new_length:)
          if (debug) print *,'Tail bit: ', buffer_array
       end if

       ! Delete the original
       if (array_associated) deallocate(array)

       ! Associate array with the new version
       array => buffer_array

       array_length = size(array)

    else
       nullify(array)
    end if

  end function splice_integer

  integer function push_integer (array, argument) result(array_length)

    integer, dimension(:), pointer :: array
    integer, intent(in)            :: argument

    ! Make an array of length one out of single argument and push
    ! this (calls array push function below)
    array_length = push(array,(/ argument /))

  end function push_integer

  integer function push_integer_array (array, argument) result(array_length)

    integer, dimension(:), pointer :: array
    integer, dimension(:), intent(in) :: argument

    integer :: current_size

    current_size = 0

    ! Grab the number of elements in array
    if (associated(array)) current_size = size(array)

    array_length = splice(array,current_size,0,argument)

  end function push_integer_array

  integer function pop_integer (array, status) result(popped)

    integer, dimension(:), pointer :: array
    logical, optional, intent(out) :: status

    integer :: current_size

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'pop :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       popped = array(ubound(array,1))
       ! Remove the last element from the array
       current_size = splice(array,-1)
    end if

  end function pop_integer

  integer function unshift_integer (array, argument) result(array_length)

    integer, dimension(:), pointer :: array
    integer, intent(in)            :: argument

    ! Make an array of length one out of single argument and unshift
    ! this (calls array unshift function below)
    array_length = unshift(array,(/ argument /))

  end function unshift_integer

  integer function unshift_integer_array (array, argument) result(array_length)

    integer, dimension(:), pointer    :: array
    integer, dimension(:), intent(in) :: argument

    integer :: start

    start = 0

    ! Find the first element (lower bound) in the array
    if (associated(array)) start = lbound(array,1)

    array_length = splice(array,0,0,argument)

  end function unshift_integer_array

  integer function shift_integer (array, status) result(shifted)

    integer, dimension(:), pointer :: array
    logical, optional, intent(out) :: status

    integer :: current_size, lower

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'shift :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       shifted = array(lbound(array,1))
       ! Remove the first element from the array (recall that offset
       ! is relative to the lower bound, so 0 = lower bound of array)
       current_size = splice(array,0,1)
    end if

  end function shift_integer

  integer function splice_real (array, offset, length, list) result(array_length)

    ! This function is meant to emulate the perl function 'splice' but on a
    ! pointer real array. So, here is the info on the perl function
    ! with some modification:
    !
  
    !   splice ARRAY,OFFSET,LENGTH,LIST
    !   splice ARRAY,OFFSET,LENGTH
    !   splice ARRAY,OFFSET
    !   splice ARRAY
    !           Removes the elements designated by OFFSET (from the lower bound of
    !           the array) and LENGTH from an
    !           array, and replaces them with the elements of LIST, if any.  
    !           Returns the number of elements removed in the array. The array 
    !           grows or shrinks as necessary.
    !           If OFFSET is negative then it starts that far from the end of the
    !           array.  If LENGTH is omitted, removes everything from OFFSET
    !           onward.  If LENGTH is negative, leaves that many elements off the
    !           end of the array.  If both OFFSET and LENGTH are omitted, removes
    !           everything.
    !
    !           The following equivalences hold (assuming the lower bound of array 
    !           is zero):
    !
    !               push(array,x)          splice(array,size(array)+1,0,x)
    !               push(array,(/x,y/))    splice(array,size(array)+1,0,(/x,y/))
    !               pop(array)             splice(array,-1)
    !               shift(array)           splice(array,0,1)
    !               unshift(array,x)       splice(array,0,0,x)
    !               array(x) = y           splice(array,x,1,y)

    ! The main differences between the perl version and this one are the definition
    ! of offset and the return value of the function. Here offset is a length from 
    ! the lower bound of the array (in perl it is an absolute position). This is
    ! so that we can usefully use arrays which start with lower bounds less than zero,
    ! as we want a negative offset to indicate counting from the end of an array.

    ! Some useful equalities to bear in mind when reading this code:
    !
    !          length of array = last element - first element + 1
    !          nth element of array = first element + n - 1

    real, dimension(:), pointer   :: array
    integer, intent(in), optional :: offset, length
    real, dimension(:), optional  :: list

    real, dimension(:), pointer :: buffer_array

    integer :: status, current_length, new_offset
    integer :: replace_length, new_length, lower, tail

    logical ::  array_associated

    array_length = 0

    if (present(offset)) then

       if (associated(array)) then
          ! Grab the number of elements in array
          current_length = size(array)
          ! And the index that the array starts at
          lower = lbound(array,1)
          array_associated = .true.
       else
          ! The array isn't associated, so we set the length to zero
          current_length = 0
          ! Assume the array index starts at one
          lower = 1
          array_associated = .false.
       end if
       if (debug) print *,'Current_length: ',current_length
       if (debug) print *,'Lower bound: ',lower

       ! We don't want to alter the return value of offset so we need to
       ! make a local version of it
       new_offset = offset

       ! A negative offset means we count from the end of the array, and
       ! -1 means the last entry of the array
       if (new_offset < 0) then
          new_offset = lower + current_length + new_offset
          if (new_offset < lower) new_offset = lower
       else
          ! Convert our offset from a relative to an absolute position
          new_offset = lower + new_offset
       end if
       if (debug) print *,'Offset: ',new_offset

       ! We don't want to alter the return value of length so we need to
       ! make a local version of it too
       if (present(length)) then
          new_length = length
          ! Tail is the number of elements left at the end of the array
          ! after offset + length elements have been removed
          tail = max(0,(lower + current_length - 1) - (new_offset + new_length - 1))
       else
          ! No length, so the new_length is to the end of the array, and
          ! there are no tail elements
          new_length = (lower + current_length - 1) - new_offset + 1
          tail = 0
       end if
       if (debug) print *,'Length: ',new_length
       if (debug) print *,'Tail: ',tail

       ! Replace length is the length of list we will put back into
       ! the array to replace the elements we are snipping out. This
       ! is initialised to zero as we are not even required to provide
       ! a replacement list
       replace_length = 0
       if (present(list)) then
          replace_length = size(list)
       end if
       if (debug) print *,'Replace length: ',replace_length

       ! Make a new array that is one less than the offset length plus
       ! the length of the list we are stuffing in plus whatever is left
       ! on the end
       ! allocate(buffer_array(lower:(lower-1)+(new_offset-1)+replace_length+tail), stat=status)
       allocate(buffer_array(lower:lower + &
            (((new_offset-1)-lower+1)+replace_length+tail) & ! length of new array
            - 1 ), stat=status)
       if (status /= 0) then
          print *,'module variable_array :: Error in allocating space for buffer array: ',status
          stop
       end if
       ! buffer_array = -9
       if (debug) print *,'Length of new array: ',size(buffer_array)

       ! Copy over everything up to (not including) offset as long as our
       ! new offset is greater than the start of the array (if not then we
       ! are probably shifting something on to the start of the array)
       if (array_associated .and. (new_offset > lower)) then
          buffer_array(lower:new_offset-1) = array(lower:new_offset-1)
          if (debug) print *,'First bit: ', buffer_array
       end if

       ! If we had a list to insert then place it now. Note that we use
       ! the replace length here ...
       if (present(list)) then
          buffer_array(new_offset:new_offset+replace_length-1) = list
          if (debug) print *,'List bit: ', buffer_array
       end if

       ! If we had elements left on the tail of array then place them on the 
       ! end. Note that we use replace length on the left and new length on
       ! the right. In this way we specify how much space our inserted text
       ! takes -- which can be as low as zero, which implies an insertion and
       ! not a deletion operation
       if (tail > 0) then
          buffer_array(new_offset+replace_length:) = array(new_offset+new_length:)
          if (debug) print *,'Tail bit: ', buffer_array
       end if

       ! Delete the original
       if (array_associated) deallocate(array)

       ! Associate array with the new version
       array => buffer_array

       array_length = size(array)

    else
       nullify(array)
    end if

  end function splice_real

  integer function push_real (array, argument) result(array_length)

    real, dimension(:), pointer :: array
    real, intent(in)            :: argument

    ! Make an array of length one out of single argument and push
    ! this (calls array push function below)
    array_length = push(array,(/ argument /))

  end function push_real

  integer function push_real_array (array, argument) result(array_length)

    real, dimension(:), pointer    :: array
    real, dimension(:), intent(in) :: argument

    integer :: current_size

    current_size = 0

    ! Grab the number of elements in array
    if (associated(array)) current_size = size(array)

    array_length = splice(array,current_size,0,argument)

  end function push_real_array

  real function pop_real (array, status) result(popped)

    real, dimension(:), pointer    :: array
    logical, optional, intent(out) :: status

    integer :: current_size

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'pop :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       popped = array(ubound(array,1))
       ! Remove the last element from the array
       current_size = splice(array,-1)
    end if

  end function pop_real

  integer function unshift_real (array, argument) result(array_length)

    real, dimension(:), pointer :: array
    real, intent(in)            :: argument

    ! Make an array of length one out of single argument and unshift
    ! this (calls array unshift function below)
    array_length = unshift(array,(/ argument /))

  end function unshift_real

  integer function unshift_real_array (array, argument) result(array_length)

    real, dimension(:), pointer    :: array
    real, dimension(:), intent(in) :: argument

    integer :: start

    start = 0

    ! Find the first element (lower bound) in the array
    if (associated(array)) start = lbound(array,1)

    array_length = splice(array,0,0,argument)

  end function unshift_real_array

  real function shift_real (array, status) result(shifted)

    real, dimension(:), pointer    :: array
    logical, optional, intent(out) :: status

    integer :: current_size, lower

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'shift :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       shifted = array(lbound(array,1))
       ! Remove the first element from the array (recall that offset
       ! is relative to the lower bound, so 0 = lower bound of array)
       current_size = splice(array,0,1)
    end if

  end function shift_real

  integer function splice_character (array, offset, length, list) result(array_length)

    ! This function is meant to emulate the perl function 'splice' but on a
    ! pointer character array. So, here is the info on the perl function
    ! with some modification:
    !
  
    !   splice ARRAY,OFFSET,LENGTH,LIST
    !   splice ARRAY,OFFSET,LENGTH
    !   splice ARRAY,OFFSET
    !   splice ARRAY
    !           Removes the elements designated by OFFSET (from the lower bound of
    !           the array) and LENGTH from an
    !           array, and replaces them with the elements of LIST, if any.  
    !           Returns the number of elements removed in the array. The array 
    !           grows or shrinks as necessary.
    !           If OFFSET is negative then it starts that far from the end of the
    !           array.  If LENGTH is omitted, removes everything from OFFSET
    !           onward.  If LENGTH is negative, leaves that many elements off the
    !           end of the array.  If both OFFSET and LENGTH are omitted, removes
    !           everything.
    !
    !           The following equivalences hold (assuming the lower bound of array 
    !           is zero):
    !
    !               push(array,x)          splice(array,size(array)+1,0,x)
    !               push(array,(/x,y/))    splice(array,size(array)+1,0,(/x,y/))
    !               pop(array)             splice(array,-1)
    !               shift(array)           splice(array,0,1)
    !               unshift(array,x)       splice(array,0,0,x)
    !               array(x) = y           splice(array,x,1,y)

    ! The main differences between the perl version and this one are the definition
    ! of offset and the return value of the function. Here offset is a length from 
    ! the lower bound of the array (in perl it is an absolute position). This is
    ! so that we can usefully use arrays which start with lower bounds less than zero,
    ! as we want a negative offset to indicate counting from the end of an array.

    ! Some useful equalities to bear in mind when reading this code:
    !
    !          length of array = last element - first element + 1
    !          nth element of array = first element + n - 1

    character(len=*), dimension(:), pointer  :: array
    integer, intent(in), optional            :: offset, length
    character(len=*), dimension(:), optional :: list

    character(len=len(array)), dimension(:), pointer :: buffer_array

    integer :: status, current_length, new_offset
    integer :: replace_length, new_length, lower, tail

    logical ::  array_associated

    array_length = 0

    if (present(offset)) then

       if (associated(array)) then
          ! Grab the number of elements in array
          current_length = size(array)
          ! And the index that the array starts at
          lower = lbound(array,1)
          array_associated = .true.
       else
          ! The array isn't associated, so we set the length to zero
          current_length = 0
          ! Assume the array index starts at one
          lower = 1
          array_associated = .false.
       end if
       if (debug) print *,'Current_length: ',current_length
       if (debug) print *,'Lower bound: ',lower

       ! We don't want to alter the return value of offset so we need to
       ! make a local version of it
       new_offset = offset

       ! A negative offset means we count from the end of the array, and
       ! -1 means the last entry of the array
       if (new_offset < 0) then
          new_offset = lower + current_length + new_offset
          if (new_offset < lower) new_offset = lower
       else
          ! Convert our offset from a relative to an absolute position
          new_offset = lower + new_offset
       end if
       if (debug) print *,'Offset: ',new_offset

       ! We don't want to alter the return value of length so we need to
       ! make a local version of it too
       if (present(length)) then
          new_length = length
          ! Tail is the number of elements left at the end of the array
          ! after offset + length elements have been removed
          tail = max(0,(lower + current_length - 1) - (new_offset + new_length - 1))
       else
          ! No length, so the new_length is to the end of the array, and
          ! there are no tail elements
          new_length = (lower + current_length - 1) - new_offset + 1
          tail = 0
       end if
       if (debug) print *,'Length: ',new_length
       if (debug) print *,'Tail: ',tail

       ! Replace length is the length of list we will put back into
       ! the array to replace the elements we are snipping out. This
       ! is initialised to zero as we are not even required to provide
       ! a replacement list
       replace_length = 0
       if (present(list)) then
          replace_length = size(list)
       end if
       if (debug) print *,'Replace length: ',replace_length

       ! Make a new array that is one less than the offset length plus
       ! the length of the list we are stuffing in plus whatever is left
       ! on the end
       ! allocate(buffer_array(lower:(lower-1)+(new_offset-1)+replace_length+tail), stat=status)
       allocate(buffer_array(lower:lower + &
            (((new_offset-1)-lower+1)+replace_length+tail) & ! length of new array
            - 1 ), stat=status)
       if (status /= 0) then
          print *,'module variable_array :: Error in allocating space for buffer array: ',status
          stop
       end if
       ! buffer_array = -9
       if (debug) print *,'Length of new array: ',size(buffer_array)

       ! Copy over everything up to (not including) offset as long as our
       ! new offset is greater than the start of the array (if not then we
       ! are probably shifting something on to the start of the array)
       if (array_associated .and. (new_offset > lower)) then
          buffer_array(lower:new_offset-1) = array(lower:new_offset-1)
          if (debug) print *,'First bit: ', buffer_array
       end if

       ! If we had a list to insert then place it now. Note that we use
       ! the replace length here ...
       if (present(list)) then
          buffer_array(new_offset:new_offset+replace_length-1) = list
          if (debug) print *,'List bit: ', buffer_array
       end if

       ! If we had elements left on the tail of array then place them on the 
       ! end. Note that we use replace length on the left and new length on
       ! the right. In this way we specify how much space our inserted text
       ! takes -- which can be as low as zero, which implies an insertion and
       ! not a deletion operation
       if (tail > 0) then
          buffer_array(new_offset+replace_length:) = array(new_offset+new_length:)
          if (debug) print *,'Tail bit: ', buffer_array
       end if

       ! Delete the original
       if (array_associated) deallocate(array)

       ! Associate array with the new version
       array => buffer_array

       array_length = size(array)

    else
       nullify(array)
    end if

  end function splice_character

  integer function push_character (array, argument) result(array_length)

    character(len=*), dimension(:), pointer :: array
    character(len=*), intent(in)            :: argument

    ! Make an array of length one out of single argument and push
    ! this (calls array push function below)
    array_length = push(array,(/ argument /))

  end function push_character

  integer function push_character_array (array, argument) result(array_length)

    character(len=*), dimension(:), pointer    :: array
    character(len=*), dimension(:), intent(in) :: argument

    integer :: current_size

    current_size = 0

    ! Grab the number of elements in array
    if (associated(array)) current_size = size(array)

    array_length = splice(array,current_size,0,argument)

  end function push_character_array

  function pop_character (array, status) result(popped)

    character(len=*), dimension(:), pointer :: array
    logical, optional, intent(out)          :: status
    
    ! Return type
    character(len=len(array))               :: popped

    integer :: current_size

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'pop :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       popped = array(ubound(array,1))
       ! Remove the last element from the array
       current_size = splice(array,-1)
    end if

  end function pop_character

  integer function unshift_character (array, argument) result(array_length)

    character(len=*), dimension(:), pointer :: array
    character(len=*), intent(in)            :: argument

    ! Make an array of length one out of single argument and unshift
    ! this (calls array unshift function below)
    array_length = unshift(array,(/ argument /))

  end function unshift_character

  integer function unshift_character_array (array, argument) result(array_length)

    character(len=*), dimension(:), pointer    :: array
    character(len=*), dimension(:), intent(in) :: argument

    integer :: start

    start = 0

    ! Find the first element (lower bound) in the array
    if (associated(array)) start = lbound(array,1)

    array_length = splice(array,0,0,argument)

  end function unshift_character_array

  function shift_character (array, status) result(shifted)

    character(len=*), dimension(:), pointer :: array
    logical, optional, intent(out)          :: status

    ! Return type
    character(len=len(array))               :: shifted

    integer :: current_size, lower

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'shift :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       shifted = array(lbound(array,1))
       ! Remove the first element from the array (recall that offset
       ! is relative to the lower bound, so 0 = lower bound of array)
       current_size = splice(array,0,1)
    end if

  end function shift_character

  integer function splice_vstring (array, offset, length, list) result(array_length)

    ! This function is meant to emulate the perl function 'splice' but on a
    ! pointer character array. So, here is the info on the perl function
    ! with some modification:
    !
  
    !   splice ARRAY,OFFSET,LENGTH,LIST
    !   splice ARRAY,OFFSET,LENGTH
    !   splice ARRAY,OFFSET
    !   splice ARRAY
    !           Removes the elements designated by OFFSET (from the lower bound of
    !           the array) and LENGTH from an
    !           array, and replaces them with the elements of LIST, if any.  
    !           Returns the number of elements removed in the array. The array 
    !           grows or shrinks as necessary.
    !           If OFFSET is negative then it starts that far from the end of the
    !           array.  If LENGTH is omitted, removes everything from OFFSET
    !           onward.  If LENGTH is negative, leaves that many elements off the
    !           end of the array.  If both OFFSET and LENGTH are omitted, removes
    !           everything.
    !
    !           The following equivalences hold (assuming the lower bound of array 
    !           is zero):
    !
    !               push(array,x)          splice(array,size(array)+1,0,x)
    !               push(array,(/x,y/))    splice(array,size(array)+1,0,(/x,y/))
    !               pop(array)             splice(array,-1)
    !               shift(array)           splice(array,0,1)
    !               unshift(array,x)       splice(array,0,0,x)
    !               array(x) = y           splice(array,x,1,y)

    ! The main differences between the perl version and this one are the definition
    ! of offset and the return value of the function. Here offset is a length from 
    ! the lower bound of the array (in perl it is an absolute position). This is
    ! so that we can usefully use arrays which start with lower bounds less than zero,
    ! as we want a negative offset to indicate counting from the end of an array.

    ! Some useful equalities to bear in mind when reading this code:
    !
    !          length of array = last element - first element + 1
    !          nth element of array = first element + n - 1

    type (varying_string), dimension(:), pointer  :: array
    integer, intent(in), optional                 :: offset, length
    type (varying_string), dimension(:), optional :: list

    type (varying_string), dimension(:), pointer :: buffer_array

    integer :: status, current_length, new_offset
    integer :: replace_length, new_length, lower, tail, i

    logical ::  array_associated

    array_length = 0

    if (present(offset)) then

       if (associated(array)) then
          ! Grab the number of elements in array
          current_length = size(array)
          ! And the index that the array starts at
          lower = lbound(array,1)
          array_associated = .true.
       else
          ! The array isn't associated, so we set the length to zero
          current_length = 0
          ! Assume the array index starts at one
          lower = 1
          array_associated = .false.
       end if
       if (debug) print *,'Current_length: ',current_length
       if (debug) print *,'Lower bound: ',lower

       ! We don't want to alter the return value of offset so we need to
       ! make a local version of it
       new_offset = offset

       ! A negative offset means we count from the end of the array, and
       ! -1 means the last entry of the array
       if (new_offset < 0) then
          new_offset = lower + current_length + new_offset
          if (new_offset < lower) new_offset = lower
       else
          ! Convert our offset from a relative to an absolute position
          new_offset = lower + new_offset
       end if
       if (debug) print *,'Offset: ',new_offset

       ! We don't want to alter the return value of length so we need to
       ! make a local version of it too
       if (present(length)) then
          new_length = length
          ! Tail is the number of elements left at the end of the array
          ! after offset + length elements have been removed
          tail = max(0,(lower + current_length - 1) - (new_offset + new_length - 1))
       else
          ! No length, so the new_length is to the end of the array, and
          ! there are no tail elements
          new_length = (lower + current_length - 1) - new_offset + 1
          tail = 0
       end if
       if (debug) print *,'Length: ',new_length
       if (debug) print *,'Tail: ',tail

       ! Replace length is the length of list we will put back into
       ! the array to replace the elements we are snipping out. This
       ! is initialised to zero as we are not even required to provide
       ! a replacement list
       replace_length = 0
       if (present(list)) then
          replace_length = size(list)
       end if
       if (debug) print *,'Replace length: ',replace_length

       ! Make a new array that is one less than the offset length plus
       ! the length of the list we are stuffing in plus whatever is left
       ! on the end
       ! allocate(buffer_array(lower:(lower-1)+(new_offset-1)+replace_length+tail), stat=status)
       allocate(buffer_array(lower:lower + &
            (((new_offset-1)-lower+1)+replace_length+tail) & ! length of new array
            - 1 ), stat=status)
       if (status /= 0) then
          print *,'module variable_array :: Error in allocating space for buffer array: ',status
          stop
       end if
       ! buffer_array = -9
       if (debug) print *,'Length of new array: ',size(buffer_array)

       ! Copy over everything up to (not including) offset as long as our
       ! new offset is greater than the start of the array (if not then we
       ! are probably shifting something on to the start of the array)
       if (array_associated .and. (new_offset > lower)) then
          buffer_array(lower:new_offset-1) = array(lower:new_offset-1)
          if (debug) print *,'First bit: ',char(buffer_array)
       end if

       ! If we had a list to insert then place it now. Note that we use
       ! the replace length here ...
       if (present(list)) then
          buffer_array(new_offset:new_offset+replace_length-1) = list
          ! if (debug) print *,'List bit: ', char(buffer_array)
       end if

       ! If we had elements left on the tail of array then place them on the 
       ! end. Note that we use replace length on the left and new length on
       ! the right. In this way we specify how much space our inserted text
       ! takes -- which can be as low as zero, which implies an insertion and
       ! not a deletion operation
       if (tail > 0) then
          buffer_array(new_offset+replace_length:) = array(new_offset+new_length:)
          ! if (debug) print *,'Tail bit: ', buffer_array
       end if

       ! Delete the original
       if (array_associated) deallocate(array)

       ! Associate array with the new version
       array => buffer_array

       array_length = size(array)

    else
       nullify(array)
    end if

  end function splice_vstring

  integer function splice_char_into_vstring (array, offset, length, list) result(array_length)

    ! This is just a wrapper for the main splice function, which
    ! allows the inserted list to be specified by a character 
    ! array rather than a varying string

    type (varying_string), dimension(:), pointer  :: array
    integer, intent(in), optional                 :: offset, length
    character(len=*), dimension(:)                :: list

    array_length = splice(array, offset, length, var_str(list))

  end function splice_char_into_vstring

  integer function push_vstring (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer :: array
    type (varying_string), intent(in)            :: argument

    ! Make an array of length one out of single argument and push
    ! this (calls array push function below)
    array_length = push(array,(/ argument /))

  end function push_vstring

  integer function push_vstring_array (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer    :: array
    type (varying_string), dimension(:), intent(in) :: argument

    integer :: current_size

    current_size = 0

    ! Grab the number of elements in array
    if (associated(array)) current_size = size(array)

    array_length = splice(array,current_size,0,argument)

  end function push_vstring_array

  integer function push_char_onto_vstring (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer :: array
    character(len=*), intent(in)                 :: argument

    ! Make an array of length one out of single argument and push
    ! this (calls array push function below)
    array_length = push(array,var_str(argument))

  end function push_char_onto_vstring

  integer function push_char_array_onto_vstring (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer :: array
    character(len=*), dimension(:), intent(in)   :: argument

    array_length = push(array, var_str(argument)) 

  end function push_char_array_onto_vstring

  function pop_vstring (array, status) result(popped)

    type (varying_string), dimension(:), pointer :: array
    logical, optional, intent(out)               :: status
    
    ! Return type
    type (varying_string) :: popped

    integer :: current_size

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'pop :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       popped = array(ubound(array,1))
       ! Remove the last element from the array
       current_size = splice(array,-1)
    end if

  end function pop_vstring

  integer function unshift_vstring (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer :: array
    type (varying_string), intent(in)            :: argument

    ! Make an array of length one out of single argument and unshift
    ! this (calls array unshift function below)
    array_length = unshift(array,(/ argument /))

  end function unshift_vstring

  integer function unshift_vstring_array (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer    :: array
    type (varying_string), dimension(:), intent(in) :: argument

    integer :: start

    start = 0

    ! Find the first element (lower bound) in the array
    if (associated(array)) start = lbound(array,1)

    array_length = splice(array,0,0,argument)

  end function unshift_vstring_array

  integer function unshift_char_onto_vstring (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer :: array
    character(len=*), intent(in)                 :: argument

    ! Make the character variable into a variable string and unshift that 
    array_length = unshift(array,var_str(argument))

  end function unshift_char_onto_vstring

  integer function unshift_char_array_onto_vstring (array, argument) result(array_length)

    type (varying_string), dimension(:), pointer :: array
    character(len=*), dimension(:), intent(in)   :: argument

    ! Make the character variable into a an array of variable strings and unshift that 
    array_length = unshift(array, var_str(argument))

  end function unshift_char_array_onto_vstring

  function shift_vstring (array, status) result(shifted)

    type (varying_string), dimension(:), pointer :: array
    logical, optional, intent(out)               :: status

    ! Return type
    type (varying_string) :: shifted

    integer :: current_size, lower

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'shift :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       shifted = array(lbound(array,1))
       ! Remove the first element from the array (recall that offset
       ! is relative to the lower bound, so 0 = lower bound of array)
       current_size = splice(array,0,1)
    end if

  end function shift_vstring

  integer function splice_logical (array, offset, length, list) result(array_length)

    ! This function is meant to emulate the perl function 'splice' but on a
    ! pointer logical array. So, here is the info on the perl function
    ! with some modification:
    !
  
    !   splice ARRAY,OFFSET,LENGTH,LIST
    !   splice ARRAY,OFFSET,LENGTH
    !   splice ARRAY,OFFSET
    !   splice ARRAY
    !           Removes the elements designated by OFFSET (from the lower bound of
    !           the array) and LENGTH from an
    !           array, and replaces them with the elements of LIST, if any.  
    !           Returns the number of elements removed in the array. The array 
    !           grows or shrinks as necessary.
    !           If OFFSET is negative then it starts that far from the end of the
    !           array.  If LENGTH is omitted, removes everything from OFFSET
    !           onward.  If LENGTH is negative, leaves that many elements off the
    !           end of the array.  If both OFFSET and LENGTH are omitted, removes
    !           everything.
    !
    !           The following equivalences hold (assuming the lower bound of array 
    !           is zero):
    !
    !               push(array,x)          splice(array,size(array)+1,0,x)
    !               push(array,(/x,y/))    splice(array,size(array)+1,0,(/x,y/))
    !               pop(array)             splice(array,-1)
    !               shift(array)           splice(array,0,1)
    !               unshift(array,x)       splice(array,0,0,x)
    !               array(x) = y           splice(array,x,1,y)

    ! The main differences between the perl version and this one are the definition
    ! of offset and the return value of the function. Here offset is a length from 
    ! the lower bound of the array (in perl it is an absolute position). This is
    ! so that we can usefully use arrays which start with lower bounds less than zero,
    ! as we want a negative offset to indicate counting from the end of an array.

    ! Some useful equalities to bear in mind when reading this code:
    !
    !          length of array = last element - first element + 1
    !          nth element of array = first element + n - 1

    logical, dimension(:), pointer   :: array
    integer, intent(in), optional :: offset, length
    logical, dimension(:), optional  :: list

    logical, dimension(:), pointer :: buffer_array

    integer :: status, current_length, new_offset
    integer :: replace_length, new_length, lower, tail

    logical ::  array_associated

    array_length = 0

    if (present(offset)) then

       if (associated(array)) then
          ! Grab the number of elements in array
          current_length = size(array)
          ! And the index that the array starts at
          lower = lbound(array,1)
          array_associated = .true.
       else
          ! The array isn't associated, so we set the length to zero
          current_length = 0
          ! Assume the array index starts at one
          lower = 1
          array_associated = .false.
       end if
       if (debug) print *,'Current_length: ',current_length
       if (debug) print *,'Lower bound: ',lower

       ! We don't want to alter the return value of offset so we need to
       ! make a local version of it
       new_offset = offset

       ! A negative offset means we count from the end of the array, and
       ! -1 means the last entry of the array
       if (new_offset < 0) then
          new_offset = lower + current_length + new_offset
          if (new_offset < lower) new_offset = lower
       else
          ! Convert our offset from a relative to an absolute position
          new_offset = lower + new_offset
       end if
       if (debug) print *,'Offset: ',new_offset

       ! We don't want to alter the return value of length so we need to
       ! make a local version of it too
       if (present(length)) then
          new_length = length
          ! Tail is the number of elements left at the end of the array
          ! after offset + length elements have been removed
          tail = max(0,(lower + current_length - 1) - (new_offset + new_length - 1))
       else
          ! No length, so the new_length is to the end of the array, and
          ! there are no tail elements
          new_length = (lower + current_length - 1) - new_offset + 1
          tail = 0
       end if
       if (debug) print *,'Length: ',new_length
       if (debug) print *,'Tail: ',tail

       ! Replace length is the length of list we will put back into
       ! the array to replace the elements we are snipping out. This
       ! is initialised to zero as we are not even required to provide
       ! a replacement list
       replace_length = 0
       if (present(list)) then
          replace_length = size(list)
       end if
       if (debug) print *,'Replace length: ',replace_length

       ! Make a new array that is one less than the offset length plus
       ! the length of the list we are stuffing in plus whatever is left
       ! on the end
       ! allocate(buffer_array(lower:(lower-1)+(new_offset-1)+replace_length+tail), stat=status)
       allocate(buffer_array(lower:lower + &
            (((new_offset-1)-lower+1)+replace_length+tail) & ! length of new array
            - 1 ), stat=status)
       if (status /= 0) then
          print *,'module variable_array :: Error in allocating space for buffer array: ',status
          stop
       end if
       ! buffer_array = -9
       if (debug) print *,'Length of new array: ',size(buffer_array)

       ! Copy over everything up to (not including) offset as long as our
       ! new offset is greater than the start of the array (if not then we
       ! are probably shifting something on to the start of the array)
       if (array_associated .and. (new_offset > lower)) then
          buffer_array(lower:new_offset-1) = array(lower:new_offset-1)
          if (debug) print *,'First bit: ', buffer_array
       end if

       ! If we had a list to insert then place it now. Note that we use
       ! the replace length here ...
       if (present(list)) then
          buffer_array(new_offset:new_offset+replace_length-1) = list
          if (debug) print *,'List bit: ', buffer_array
       end if

       ! If we had elements left on the tail of array then place them on the 
       ! end. Note that we use replace length on the left and new length on
       ! the right. In this way we specify how much space our inserted text
       ! takes -- which can be as low as zero, which implies an insertion and
       ! not a deletion operation
       if (tail > 0) then
          buffer_array(new_offset+replace_length:) = array(new_offset+new_length:)
          if (debug) print *,'Tail bit: ', buffer_array
       end if

       ! Delete the original
       if (array_associated) deallocate(array)

       ! Associate array with the new version
       array => buffer_array

       array_length = size(array)

    else
       nullify(array)
    end if

  end function splice_logical

  integer function push_logical (array, argument) result(array_length)

    logical, dimension(:), pointer :: array
    logical, intent(in)            :: argument

    ! Make an array of length one out of single argument and push
    ! this (calls array push function below)
    array_length = push(array,(/ argument /))

  end function push_logical

  integer function push_logical_array (array, argument) result(array_length)

    logical, dimension(:), pointer    :: array
    logical, dimension(:), intent(in) :: argument

    integer :: current_size

    current_size = 0

    ! Grab the number of elements in array
    if (associated(array)) current_size = size(array)

    array_length = splice(array,current_size,0,argument)

  end function push_logical_array

  logical function pop_logical (array, status) result(popped)

    logical, dimension(:), pointer    :: array
    logical, optional, intent(out) :: status

    integer :: current_size

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'pop :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       popped = array(ubound(array,1))
       ! Remove the last element from the array
       current_size = splice(array,-1)
    end if

  end function pop_logical

  integer function unshift_logical (array, argument) result(array_length)

    logical, dimension(:), pointer :: array
    logical, intent(in)            :: argument

    ! Make an array of length one out of single argument and unshift
    ! this (calls array unshift function below)
    array_length = unshift(array,(/ argument /))

  end function unshift_logical

  integer function unshift_logical_array (array, argument) result(array_length)

    logical, dimension(:), pointer    :: array
    logical, dimension(:), intent(in) :: argument

    integer :: start

    start = 0

    ! Find the first element (lower bound) in the array
    if (associated(array)) start = lbound(array,1)

    array_length = splice(array,0,0,argument)

  end function unshift_logical_array

  logical function shift_logical (array, status) result(shifted)

    logical, dimension(:), pointer    :: array
    logical, optional, intent(out) :: status

    integer :: current_size, lower

    ! Grab the number of elements in array
    if (.not. associated(array)) then
       if (present(status)) then
          status = .false.
       else
          print *,'shift :: Attempted to access zero size array!'
          stop
       end if
    else
       if (present(status)) status = .true.
       shifted = array(lbound(array,1))
       ! Remove the first element from the array (recall that offset
       ! is relative to the lower bound, so 0 = lower bound of array)
       current_size = splice(array,0,1)
    end if

  end function shift_logical


end module variable_array
