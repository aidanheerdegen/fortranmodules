module hash_table

  use file_functions, only: stderr, stdout
  use iso_varying_string

  implicit none

  private

  ! This is almost entirely filched from a public domain C hash table
  ! program by Bob Jenkins (see http://burtleburtle.net/bob/hash/evahash.html).
  ! Translation from C into Fortran by Aidan Heerdegen (2003).
  !
  ! The C version could store anything you could point to with a void pointer,
  ! but we don't have that freedom in Fortran.

  !! $Log: hash_table.f90,v $
  !! Revision 1.6  2006/02/08 04:49:21  aidan
  !! Removed non-intrinsic 'size' function from dimension specification
  !! in two interfaces. Replaced with the intrinsic size function. This
  !! is standard conforming f90.
  !!
  !! Revision 1.5  2004/05/11 00:08:43  aidan
  !! No longer need to use the binary module.
  !!
  !! Revision 1.4  2003/10/17 05:29:06  aidan
  !! Switched to a user defined type (raw_bytes) for the internal storage of
  !! keys and values. This allows all the type coercion on output (from both
  !! value and key access routines) to use the same small number of assignment
  !! routines. Now it is possible to avoid calling the value and key routines
  !! wih a 'mold' parameter, and instead assign the output to an integer/real/
  !! string/varying string etc and have the correct behaviour via assignment.
  !! The hash module can easily be extended by other modules simply by
  !! providing an appropriate hash_add function and assignment function to
  !! coerce the value back to the original type.
  !!
  !! Revision 1.3  2003/09/30 00:31:19  aidan
  !! Added keys and values functions -- these return arrays of all the keys and
  !! values in a hash table. Currently there is no support for returning arrays
  !! of arrays of bytes--the internal format in which all keys and values are
  !! stored. This would require the use of a user defined type to store an array
  !! of bytes--an array of such types would be returned to the user. This is
  !! more than a little ugly and I have decided not to do it at this time. This
  !! is in effect what we do when we return an array of strings--we use the
  !! iso_varying_string module to handle the 'ragged' array of strings. This is
  !! acceptable because the varying string module provides support for easy
  !! conversion from varying string to character for printing etc. If I add in
  !! a full blown byte-array type to a module I might add in support for arrays
  !! of byte value and keys.
  !!
  !! Made the count function pure so it could be used to dimension return
  !! values in the keys and values functions.
  !!
  !! Changed the hash function so that it does what it was supposed to do--namely
  !! return the value of the hash for the currently selected hash item.
  !!
  !! Added some comments to the module interface.
  !!
  !! Revision 1.2  2003/09/26 05:02:36  aidan
  !! Added support for integer and character keys, and integer, character and
  !! real values. Using these to add entries to the hash is transparent, but
  !! retrieving them requires passing a 'mold' variable to the appropriate
  !! routine--this allows the compiler to differentiate between the different
  !! types and return the right one. Internally they are all stored as 1-byte
  !! arrays and type coercion (the 'transfer' function) is used to convert
  !! them to the correct output type.
  !!
  !! Now check for a minimum size in hash_create.
  !!
  !! Also check the sizes of the tables in hash_assign, and increase the
  !! size of the table to be assigned to if it is not large enough. Also
  !! now check if an item in the table being assigned to already exists,
  !! and if so we delete it. Otherwise the assign would fail. The status
  !! of the table being assigned to is now 'intent(inout)', so that it
  !! keeps it's size and any hash items that are not in the table being
  !! assigned from.
  !!
  !! Revision 1.1  2003/09/24 06:46:10  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: hash_table.f90,v 1.6 2006/02/08 04:49:21 aidan Exp $"

  ! This is the internal format for storing keys and values
  type raw_bytes
     integer(kind=1), dimension(:), pointer :: bytes => null() ! Array of raw bytes
  end type raw_bytes

  ! Each key/value pair is stored in a hash_item
  type hash_item
     private
     type (raw_bytes) :: key    ! key that is hashed
     type (raw_bytes) :: value  ! value stored in this item
     integer :: hash            ! hashed key
     type (hash_item), pointer :: next => null() ! next hash item with the same hash key
                                                 ! (chained collision)
  end type hash_item

  ! The hash table is an array of hash items, and a pointer to the
  ! currently selected hash item
  type hash_table_object
     private
     type (hash_item), dimension(:), pointer :: table => null() ! table of hash items
     integer :: logsize                ! log of the size of the table
     integer :: mask                   ! (hash & mask) is the position in the table
     integer :: count                  ! how many items in the hash table so far?
     integer :: apos                   ! current position in the array
     type (hash_item), pointer :: ipos ! pointer to the current item
  end type hash_table_object

  ! Make a new hash or a new item
  interface new
     module procedure hash_create, item_create
  end interface

  ! Destroy a hash or an item
  interface destroy
     module procedure hash_destroy, item_destroy
  end interface

  ! Make a hash table larger
  interface grow
     module procedure hash_grow
  end interface

  ! Add a key/value entry to the hash table
  interface add
     module procedure hash_add, &
          hash_add_string_key_value, hash_add_string_key_real_value, &
          hash_add_string_key_int_value, hash_add_string_key, &
          hash_add_int_key_int_value, hash_add_int_key_string_value, &
          hash_add_int_key_real_value, hash_add_int_key, &
          hash_add_int_value, hash_add_string_value, hash_add_real_value
  end interface

  ! Set the currently selected item to be the next item in the hash
  interface next
     module procedure next_item
  end interface

  ! Set the currently selected item to the first item in the hash table
  interface first
     module procedure first_bucket
  end interface

  ! Tests to see if a key exists in a hash table, if so this key is
  ! made the currently selected item
  interface exists
     module procedure find_item, find_item_int_key, find_item_real_key, find_item_string_key
  end interface

  ! Delete the currently selected key/value entry from the hash table
  interface delete
     module procedure delete_item
  end interface

  ! Mixes 3 32-bit values (used in the hashing algorithm)
  interface mix
     module procedure mix_32_bit
  end interface

  ! Return the key of the currently selected item
  interface key
     module procedure hash_key, hash_key_as_integer, hash_key_as_string, hash_key_as_bytes
  end interface

  ! Return the keys in a hash table as an array
  interface keys
     module procedure hash_keys, hash_keys_as_integers, hash_keys_as_strings
  end interface

  ! Return the value of the currently selected item
  interface value
     module procedure hash_value, hash_value_as_bytes, hash_value_as_int, hash_value_as_real, &
          hash_value_as_string
  end interface

  ! Return the values in a hash table as an array
  interface values
     module procedure hash_values, hash_values_as_ints, hash_values_as_reals, hash_values_as_strings
  end interface

  ! Returns the current number of items in a hash table
  interface count
     module procedure hash_count
  end interface

  ! Returns the current size of a raw bytes array
  interface size
     module procedure size_rawbytes
  end interface

  ! A generic method for printing out the contents of a hash. Not terribly useful as
  ! it has to print out all hashes in the same manner, regardless of what type the
  ! keys and values were input as, and intended to be retrieved as. Useful for 
  ! debugging purposes.
  interface print
     module procedure print_hash_to_stdout, print_hash, print_item_to_stdout, print_item
  end interface

  ! Overloaded assignment so that hash and item assignment makes new copies of the 
  ! data contained in the pointers, rather than just copying their locations. Also
  ! defines assignment from our internal data format to integer, real, string, 
  ! varying string and raw bytes (integer(kind=1) aarray). This is where all the
  ! type coercion for output occurs.
  interface assignment(=)
     module procedure hash_assign, item_assign, ass_rawbytes_to_int, ass_rawbytes_to_real, & 
          ass_rawbytes_to_string, ass_rawbytes_to_varstring, ass_rawbytes_to_rawbytes
  end interface

  ! Public types
  public :: hash_table_object, hash_item, raw_bytes

  ! Public routines
  public :: mix, new, print

  ! Overloaded functions
  public :: count, size

  ! Public functions
  public :: lookup, add, next, first, key, keys, value, values, hash, exists, destroy, delete

  ! Public operators
  public :: assignment(=)

  logical, parameter :: debug =   .FALSE. !.TRUE. !

contains

  elemental integer function size_rawbytes (raw)

    ! Returns the size of the array contained in the raw_bytes object
  
    type (raw_bytes), intent(in) :: raw

    size_rawbytes = size(raw%bytes)

  end function size_rawbytes

  ! Assignment from internal data format to other, more useful, types. These
  ! routines contain all the type coercion for output in this module.

  pure subroutine ass_rawbytes_to_rawbytes (intval, raw)

    ! Assign a raw byte value (probably a key or value returned from a hash)
    ! to an integer(kind=1) array

    type (raw_bytes), intent(in) :: raw
    ! integer(kind=1), dimension(size(raw%bytes)), intent(out)         :: intval
    integer(kind=1), dimension(size(raw%bytes)), intent(out)         :: intval

    intval = transfer(raw%bytes,intval)

  end subroutine ass_rawbytes_to_rawbytes
  
  elemental subroutine ass_rawbytes_to_int (intval, raw)

    ! Assign a raw byte value (probably a key or value returned from a hash)
    ! to an integer

    integer, intent(out)         :: intval
    type (raw_bytes), intent(in) :: raw

    intval = transfer(raw%bytes,intval)

  end subroutine ass_rawbytes_to_int
  
  elemental subroutine ass_rawbytes_to_real (realval, raw)

    ! Assign a raw byte value (probably a key or value returned from a hash)
    ! to a real

    type (raw_bytes), intent(in) :: raw
    real, intent(out)            :: realval

    realval = transfer(raw%bytes,realval)

  end subroutine ass_rawbytes_to_real
  
  pure subroutine ass_rawbytes_to_string (stringval, raw)

    ! Assign a raw byte value (probably a key or value returned from a hash)
    ! to a string

    type (raw_bytes), intent(in)                :: raw
    character(len=size(raw%bytes)), intent(out) :: stringval

    stringval = transfer(raw%bytes,stringval)

  end subroutine ass_rawbytes_to_string
  
  elemental subroutine ass_rawbytes_to_varstring (stringval, raw)

    ! Assign a raw byte value (probably a key or value returned from a hash)
    ! to a string

    type (raw_bytes), intent(in)                :: raw
    ! character(len=size(raw%bytes,1)), intent(out) :: stringval
    type (varying_string), intent(out) :: stringval

    character(len=1) :: mold

    stringval = var_str(transfer(raw%bytes,repeat(mold,size(raw%bytes))))

  end subroutine ass_rawbytes_to_varstring
  
  function hash_key (table)

    ! Return key of currently selected hash table item 

    type (hash_table_object), intent(in) :: table

    ! integer(kind=1), dimension(size(table%ipos%key)) :: hash_key
    type (raw_bytes) :: hash_key
    
    hash_key = table%ipos%key

  end function hash_key
  
  function hash_key_as_bytes (table, mold) result(hashkey)

    ! Return key of currently selected hash table item translated
    ! into a character string

    type (hash_table_object), intent(in) :: table
    integer(kind=1), intent(in)          :: mold

    ! Return type
    integer(kind=1), dimension(size(table%ipos%key%bytes)) :: hashkey

    ! key = table%ipos%key%bytes
    hashkey = key(table)

  end function hash_key_as_bytes
  
  function hash_key_as_integer (table, mold) result(hashkey)

    ! Return key of currently selected hash table item translated
    ! into a character string

    type (hash_table_object), intent(in) :: table
    integer, intent(in)                  :: mold

    ! Return type
    integer :: hashkey

    hashkey = key(table)

  end function hash_key_as_integer
  
  function hash_key_as_string (table, mold) result(hashkey)

    ! Return key of currently selected hash table item translated
    ! into a character string

    type (hash_table_object), intent(in) :: table
    character(len=*), intent(in)         :: mold

    ! Return type
    character(len=size(table%ipos%key%bytes)) :: hashkey

    hashkey = key(table)

  end function hash_key_as_string
  
  function hash_keys (table) result(hashkeys)

    ! Return array of all keys in hash table

    ! Input variables
    type (hash_table_object), intent(inout) :: table

    ! Return type
    type (raw_bytes), dimension(count(table)) :: hashkeys

    integer :: i
    
    i = 0

    ! Traverse the table from the start collecting
    ! all our keys into our array
    if (first(table)) then
       do
          i = i + 1
          hashkeys(i) = table%ipos%key
          if (.NOT. next(table)) exit
       end do
    end if

  end function hash_keys
  
  function hash_keys_as_integers (table, mold) result(hashkeys)

    ! Return array of all keys in hash table translated
    ! into integers

    type (hash_table_object), intent(inout) :: table
    integer, intent(inout)                  :: mold

    ! Return type
    integer, dimension(count(table)) :: hashkeys

    hashkeys = keys(table)

  end function hash_keys_as_integers
  
  function hash_keys_as_strings (table, mold) result(hashkeys)

    ! Return array of all keys in hash table translated
    ! into strings

    type (hash_table_object), intent(inout) :: table
    character(len=*), intent(in)            :: mold

    ! Return type
    type (varying_string), dimension(count(table)) :: hashkeys

    hashkeys = keys(table)

  end function hash_keys_as_strings
  
  function hash (table)

    ! Return value of currently selected hash table item 

    type (hash_table_object), intent(in) :: table

    integer :: hash
    
    hash = table%ipos%hash

  end function hash
  
  function hash_value (table) result(value)

    ! Return value of currently selected hash table item 

    type (hash_table_object), intent(in) :: table

    type(raw_bytes) :: value
    
    value = table%ipos%value !%bytes

  end function hash_value
  
  function hash_value_as_bytes (table, mold) result(hashvalue)

    ! Return value of currently selected hash table item 

    type (hash_table_object), intent(in) :: table
    integer(kind=1), intent(in)          :: mold

    integer(kind=1), dimension(size(table%ipos%value%bytes)) :: hashvalue
    
    ! Use default assignment to coerce the output of value to raw bytes 
    hashvalue = value(table)

  end function hash_value_as_bytes
  
  function hash_value_as_int (table, mold) result(hashvalue)

    ! Return value of currently selected hash table item 
    ! as an integer

    type (hash_table_object), intent(in) :: table
    integer, intent(in)                  :: mold

    integer :: hashvalue

    ! Use default assignment to coerce the output of value to integer 
    hashvalue = value(table)

  end function hash_value_as_int
  
  function hash_value_as_string (table, mold) result(hashvalue)

    ! Return value of currently selected hash table item 
    ! as a character string

    type (hash_table_object), intent(in) :: table
    character(len=*), intent(in)         :: mold

    character(len=size(table%ipos%value%bytes)) :: hashvalue

    ! Use default assignment to coerce the output of value to string 
    hashvalue = value(table)

  end function hash_value_as_string
  
  function hash_value_as_real (table, mold) result(hashvalue)

    ! Return value of currently selected hash table item 
    ! as an integer

    type (hash_table_object), intent(in) :: table
    real, intent(in)                     :: mold

    real :: hashvalue

    ! Use default assignment to coerce the output of value to real 
    hashvalue = value(table)

  end function hash_value_as_real
  
  function hash_values (table) result(hashvalues)

    ! Return array of all values in hash table translated
    ! into integers

    type (hash_table_object), intent(inout) :: table

    ! Return type
    type (raw_bytes), dimension(count(table)) :: hashvalues

    ! Local variables
    integer :: i

    i = 0

    ! Traverse the table from the start collecting
    ! all our keys into our array
    if (first(table)) then
       do
          i = i + 1
          hashvalues(i) = value(table)
          if (.NOT. next(table)) exit
       end do
    end if

  end function hash_values
  
  function hash_values_as_ints (table, mold) result(hashvalues)

    ! Return array of all values in hash table translated
    ! into integers

    type (hash_table_object), intent(inout) :: table
    integer, intent(in)                     :: mold

    ! Return type
    integer, dimension(count(table)) :: hashvalues

    ! Use default assignment to coerce the output of values to integer 
    hashvalues = values(table)

  end function hash_values_as_ints
  
  function hash_values_as_strings (table, mold) result(hashvalues)

    ! Return array of all values in hash table translated
    ! into integers

    type (hash_table_object), intent(inout) :: table
    character(len=*), intent(in)            :: mold

    ! Return type
    type (varying_string), dimension(count(table)) :: hashvalues

    ! Use default assignment to coerce the output of values to varying string 
    hashvalues = values(table)

  end function hash_values_as_strings
  
  function hash_values_as_reals (table, mold) result(hashvalues)

    ! Return array of all values in hash table translated
    ! into integers

    type (hash_table_object), intent(inout) :: table
    real, intent(in)                        :: mold

    ! Return type
    real, dimension(count(table)) :: hashvalues

    ! Use default assignment to coerce the output of values to real 
    hashvalues = values(table)

  end function hash_values_as_reals
  
  pure integer function hash_count (table)

    ! Return the number of items in the hash table

    type (hash_table_object), intent(in) :: table

    hash_count = table%count

  end function hash_count
  
  subroutine hash_create (table, size)

    ! Create a hash table

    ! Input variables
    type (hash_table_object), intent(out) :: table
    integer, intent(in), optional         :: size

    ! Local variables
    integer :: len, error, logsize

    ! The size of our hash table is defined in log units, 
    ! i.e. actual size = 2**(logsize). Default size is 256.
    logsize = 8
    if (present(size)) logsize = size

    if (logsize < 1) then
       write(stderr,*) 'HASH_TABLE :: size must be >=1'
       stop
    end if

    len = ishft(1_4,logsize)

    if (associated(table%table)) then
       if (.NOT. destroy(table)) then
          write(stderr,*) 'HASH_TABLE :: Cannot free memory from previous table when creating new hash table'
          stop
       end if
    end if

    ! The hash table is just an array of hash entries -- allocate
    ! the space for them all now and query the association status
    ! of the key for each entry to determine if it has been filled
    allocate(table%table(0:len), stat=error)
    if (error /=0) then
       write(stderr,*) 'HASH_TABLE :: Error allocating memory for hash table'
       stop
    end if

    ! The mask is used to generate indices into our table from
    ! the 'hash' value using a bitwise AND operation, i.e. the
    ! binary representation of 256 is 100000000 and 255 is 
    ! 11111111 -- a bitwise AND between the latter and a hash
    ! value would pick out the 8 least significant bits of the 
    ! hash.
    table%mask    = len - 1
    table%logsize = logsize
    table%count   = 0
    table%apos    = 0
    nullify(table%ipos)
    
  end subroutine hash_create
  
  subroutine item_create (item, key, value, hash, next)

    ! Create a hash item 

    type (hash_item), intent(out)             :: item
    integer(kind=1), dimension(:), intent(in) :: key, value
    integer, intent(in)                       :: hash
    type (hash_item), pointer, optional       :: next

    ! Local variables
    integer :: error

    ! Allocate enough space to store our key, which is just
    ! an array of bytes 
    allocate(item%key%bytes(size(key)), stat=error)
    if (error /=0) then
       write(stderr,*) 'HASH_TABLE :: Error allocating memory while creating hash item'
       stop
    end if

    ! Allocate enough space to store our value, which is just
    ! an array of bytes 
    allocate(item%value%bytes(size(value)), stat=error)
    if (error /=0) then
       write(stderr,*) 'HASH_TABLE :: Error allocating memory while creating hash item'
       stop
    end if

    ! Assign the values to the new item
    item%key%bytes   = key
    item%value%bytes = value
    item%hash        = hash

    if (present(next)) then
       item%next => next
    else
       ! Allocate the space for our 'chained' item -- no chained item exists yet,
       ! but we need to be able to test the association status of the key array in
       ! the chained item
       allocate(item%next, stat=error)
       if (error /=0) then
          write(stderr,*) 'HASH_TABLE :: Error allocating chained item while creating hash item'
          stop
       end if
    end if
       
  end subroutine item_create
  
  subroutine hash_grow (table)

    ! Double the size of a hash table. 
    
    ! Input variable
    type (hash_table_object), intent(inout) :: table

    ! Local variables
    type (hash_table_object)  :: newtable  ! Temporary table
    type (hash_item), pointer :: hitem     ! Temporary hash item 
    integer :: logsize

    ! Allocate a new, 2x bigger array,
    call new(newtable, table%logsize+1)

    ! Use the assignment operator to copy all the elements in the
    ! hash table into our new, larger version
    newtable = table

    ! Destroy the old table
    if (.NOT. destroy(table)) then
       write(stderr,*) 'HASH_TABLE :: Error while copying elements when growing table'
       stop
    end if

    ! Do a manual 'shallow copy', where we point the 'old' table at 
    ! the new larger version's memory locations and values
    table%table   => newtable%table
    table%logsize =  newtable%logsize
    table%mask    =  newtable%mask 
    table%count   =  newtable%count 
    table%apos    =  newtable%apos
    table%ipos    => newtable%ipos

    ! We don't want to clean up the memory allocated to the new table
    ! as we are pointing to it in our doubled table.

    if (debug) call print(stdout,table)

  end subroutine hash_grow
  
  subroutine hash_assign (tablel, tabler)

    ! Move everything from the right table to the left table

    type (hash_table_object), intent(inout) :: tablel
    type (hash_table_object), intent(in)  :: tabler

    ! Local variables 
    type (hash_item), pointer :: hitem     ! Temporary hash item 
    integer :: i

    if (.NOT. associated(tablel%table)) then
       call new(tablel,tabler%logsize)
    else if (tablel%logsize < tabler%logsize) then
       call new(tablel,tabler%logsize)
    end if

    ! We can't use the nicer first(tabler), next(tabler) syntax as 
    ! this alters the internal state of the hash table object, and 
    ! tabler MUST BE intent(in) when overloading the assignment operator.
    do i = 0, ishft(1_4,tabler%logsize)
       ! Point at the hash item with this index
       hitem => tabler%table(i)
       do while (associated(hitem%key%bytes))
          ! Add this key/value combo to the left table
          if ( exists(tablel, hitem%key%bytes) ) then
             if (.NOT. delete(tablel)) then
                write(stderr,*) 'HASH_TABLE :: Could not overwrite existing entry when copying table'
                stop
             end if
          end if
          if ( .NOT. add(tablel, hitem%key%bytes, hitem%value%bytes)) then
             write(stderr,*) 'HASH_TABLE :: Error while copying elements when copying table '
             stop
          end if
          ! Make the hash item the next in the chain for this key
          hitem => hitem%next
       end do
    end do

  end subroutine hash_assign

  subroutine item_assign (iteml, itemr)

    ! Assign the values of one item to another item

    type (hash_item), intent(out) :: iteml
    type (hash_item), intent(in)  :: itemr

    call new(iteml, itemr%key%bytes, itemr%value%bytes, itemr%hash, itemr%next)

  end subroutine item_assign

  logical recursive function item_destroy (item) result(destroyed)

    ! item_destroy - deallocate all the pointers contained in an item
    type (hash_item), pointer :: item

    ! Local variables
    integer :: i, error1, error2

    destroyed = .TRUE.

    if (associated(item%key%bytes)) then
       ! Try and destroy any items chained to this one before
       ! we kill this item
       if (destroy(item%next)) then
          deallocate(item%value%bytes, stat=error1)
          deallocate(item%key%bytes, stat=error2)
          destroyed = ( (error1 == 0) .AND. (error2 == 0) )
       else
          write(stderr,*) 'HASH_TABLE :: Could not free memory in hash item '
          stop
       end if
    end if

  end function item_destroy
  
  logical function hash_destroy (table) result(destroyed)

    ! hash_destroy - deallocate all the pointers contained in a hash
    type (hash_table_object), intent(inout) :: table

    ! Local variables
    integer :: error

    destroyed = .TRUE.

    ! Cycle through the hash table and destroy each item in turn
    if (first(table)) then
       do 
          ! Destroy the current hash item
          if (.NOT. destroy(table%ipos)) then
             write(stderr,*) 'HASH_TABLE :: Could not destroy hash item'
             destroyed = .FALSE.
             ! stop
          end if
          if (.NOT. next_bucket(table)) exit
       end do
    end if

    ! Free the memory allocated to the hash table
    deallocate(table%table, stat=error)
    if (error /= 0) then
       write(stderr,*) 'HASH_TABLE :: Could not destroy hash table'
       destroyed = .FALSE.
       ! stop
    end if

    ! Reset the values in the hash table
    table%logsize = 0
    table%mask    = 0
    table%count   = 0
    table%apos    = 0
    nullify(table%ipos)

  end function hash_destroy

  ! The following are wrappers to the hash_add function to allow
  ! string and integer valued keys, and integer, string and real
  ! values

  ! String valued keys

  logical function hash_add_string_key_value (table, key, value) result(added_ok)

    ! Add a string key/string value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    character(len=*), intent(in)              :: key, value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_string_key_value

  logical function hash_add_string_key_real_value (table, key, value) result(added_ok)

    ! Add a string key/real value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    character(len=*), intent(in)              :: key
    real, intent(in)                          :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_string_key_real_value

  logical function hash_add_string_key_int_value (table, key, value) result(added_ok)

    ! Add a string key/real value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    character(len=*), intent(in)              :: key
    integer, intent(in)                       :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_string_key_int_value

  logical function hash_add_string_key (table, key, value) result(added_ok)

    ! Add a string key/value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    character(len=*), intent(in)              :: key
    integer(kind=1), dimension(:), intent(in) :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),value)
    
  end function hash_add_string_key

  ! Integer valued keys

  logical function hash_add_int_key_int_value (table, key, value) result(added_ok)

    ! Add a integer key/integer value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer, intent(in)                       :: key, value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_int_key_int_value

  logical function hash_add_int_key_string_value (table, key, value) result(added_ok)

    ! Add a integer key/string value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer, intent(in)                       :: key
    character(len=*), intent(in)              :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_int_key_string_value

  logical function hash_add_int_key_real_value (table, key, value) result(added_ok)

    ! Add a integer key/real value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer, intent(in)                       :: key
    real, intent(in)                          :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_int_key_real_value

  logical function hash_add_int_key (table, key, value) result(added_ok)

    ! Add a integer key/value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer, intent(in)                       :: key
    integer(kind=1), dimension(:), intent(in) :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),value)
    
  end function hash_add_int_key

  ! Integer array keys with the various types of values

  logical function hash_add_int_value (table, key, value) result(added_ok)

    ! Add a integer key/integer value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer(kind=1), dimension(:), intent(in) :: key
    integer, intent(in)                       :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_int_value

  logical function hash_add_string_value (table, key, value) result(added_ok)

    ! Add a integer key/string value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer(kind=1), dimension(:), intent(in) :: key
    character(len=*), intent(in)              :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_string_value

  logical function hash_add_real_value (table, key, value) result(added_ok)

    ! Add a integer key/real value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer(kind=1), dimension(:), intent(in) :: key
    real, intent(in)                          :: value

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,transfer(key,mold),transfer(value,mold))
    
  end function hash_add_real_value

  logical function hash_add_vstring_array_value (table, key, value) result(added_ok)

    ! Add a varying string array value to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer(kind=1), dimension(:), intent(in) :: key
    type (varying_string), intent(in)         :: value(:)

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    added_ok = add(table,key,transfer(value,mold))
    
  end function hash_add_vstring_array_value

  logical function hash_add (table, key, value) result(added_ok)

    ! Add a key/value pair to the hash table

    type (hash_table_object), intent(inout)   :: table
    integer(kind=1), dimension(:), intent(in) :: key, value
    ! integer, intent(in)                       :: value

    ! Local variables
    type (hash_item), pointer :: hitem ! Temporary hash item 

    integer :: hash, index

    ! We don't use the exists function here because we need to
    ! calculate the hash value anyway. We also need to set
    ! hitem to something sensible if we find this key hasn't 
    ! been defined already, and this varies depending if we 
    ! are a chained item or a brand spanking new index
    
    ! Make a 32-bit hash out of the key
    hash = lookup(key)

    ! Convert hash into an index in our hash table by masking it
    index = iand(hash, table%mask)

    ! Point to this item in the hash table
    hitem => table%table(index)

    do while (associated(hitem%key%bytes))

       ! Make sure this key isn't there already ...
       if ( (hash == hitem%hash) .and. all(key == hitem%key%bytes) ) then
          table%apos = index
          table%ipos => hitem
          added_ok = .FALSE.
          return
       end if
       ! Make the hash item the next in the chain for this key
       hitem => hitem%next
    end do

    ! Use the constructor to make a new hash item
    call new(hitem, key, value, hash)

    ! Set our current position appropriately
    table%apos = index
    table%ipos => hitem

    ! Increment the table count
    table%count = table%count + 1

    ! Grow the table if it is getting too big
    if (table%count >= ishft(1_4,table%logsize-1)) then

       ! Double the size of the table
       call grow(table)

       ! Reset the current table position to be the last entry read
       if (.NOT. exists(table,key)) then
          write(stderr,*) 'HASH_TABLE :: Could not reset table position after growth'
          stop
       end if

    end if

    ! Make sure we return a successful value
    added_ok = .TRUE.

  end function hash_add

  logical function delete_item (table) result(deleted)

    ! Delete the currently selected hash item

    type (hash_table_object), intent(inout) :: table

    ! Local variables
    type (hash_item), pointer :: hitem, prev_hitem ! Temporary hash items

    integer :: error1, error2
    logical :: is_root

    ! Point to the root item for the currently selected item
    hitem => table%table(table%apos)
    nullify(prev_hitem)

    is_root = .TRUE.

    error2 = 0
    
    do
       ! Test if this is the item we must delete
       if ( all(table%ipos%key%bytes == hitem%key%bytes) ) then

          ! Free the space used by the key array
          deallocate(hitem%key%bytes, stat=error1)

          ! Is there an item chained to this one?
          if (associated(hitem%next%key%bytes)) then
             if (is_root) then
                ! The deleted item is the root item, so make the 
                ! chained item the root item (this uses the assignment 
                ! operator defined in item_assign)
                table%table(table%apos) = hitem%next
             else
                ! Just skip our deleted item by pointing
                ! the previous item at the next item in the chain
                prev_hitem%next => hitem%next
                ! Free the space being used by the item
                deallocate(hitem, stat=error2)
             end if
          end if

          ! Requires no errors in deallocation 
          deleted = ((error1 == 0) .AND. (error2 == 0))

          ! Reset the current position--we ignore the return status because
          ! we don't care if it had to wrap to find a hash item
          if (next(table)) continue

          ! Decrement the count of the number of items in the table
          if (deleted) table%count = table%count - 1

          return

       end if

       ! Make the hash item the next in the chain for this key
       is_root = .FALSE.
       prev_hitem => hitem
       hitem => hitem%next

    end do

    deleted = .FALSE.

  end function delete_item

  logical function find_item_int_key (table, key) result(found)

    type (hash_table_object), intent(inout) :: table
    integer, intent(in)                     :: key

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    found = exists(table, transfer(key,mold))

  end function find_item_int_key

  logical function find_item_string_key (table, key) result(found)

    type (hash_table_object), intent(inout) :: table
    character(len=*), intent(in)            :: key

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    found = exists(table, transfer(key,mold))

  end function find_item_string_key

  logical function find_item_real_key (table, key) result(found)

    type (hash_table_object), intent(inout) :: table
    real, intent(in)                        :: key

    ! Local variables
    integer(kind=1), dimension(1) :: mold

    found = exists(table, transfer(key,mold))

  end function find_item_real_key

  logical function find_item (table, key) result(found)

    ! Find a hash item with the specified key

    type (hash_table_object), intent(inout)   :: table
    integer(kind=1), dimension(:), intent(in) :: key

    ! Local variables
    type (hash_item), pointer :: hitem ! Temporary hash item 

    integer :: hash, index
    
    ! Will assume it dosen't exist ...
    found = .FALSE.

    ! Get the hash value
    hash = lookup(key)

    ! Make an index out of it
    index = iand(hash, table%mask)

    ! Point at the hash item with this index
    hitem => table%table(index)

    ! While we have a valid key--we can tell if this is the case by 
    ! checking if the key array is associated--keep stepping down
    ! the chain of hash items until we find a match.
    do while (associated(hitem%key%bytes))
       if ( (hash == hitem%hash) .and. all(key == hitem%key%bytes) ) then
          ! Found our entry .. make it the current entry, set the
          ! return value to true and exit the function
          table%apos = index
          table%ipos => hitem
          found = .TRUE.
          return
       end if
       ! Make the hash item the next in the chain for this key
       hitem => hitem%next
    end do

  end function find_item

  logical function next_bucket (table)

    ! Set the current location to be the next available bucket, 
    ! i.e. the root item in the next available index in the hash table 

    type (hash_table_object), intent(inout)   :: table

    ! Local variables
    integer :: i, index, end, oldapos
    
    next_bucket = .FALSE.

    end = ishft(1_4,table%logsize)-1
    oldapos = table%apos

    ! See if the element can be found without wrapping around
    do i = oldapos+1,end
       index = iand(i, table%mask)
       if (associated(table%table(index)%key%bytes)) then
          table%apos = i
          table%ipos => table%table(i)
          next_bucket = .TRUE.
          return
       end if
    end do

    ! Have to wrap around to find the last element
    do i = 0,oldapos
       if (associated(table%table(i)%key%bytes)) then
          table%apos = i
          table%ipos => table%table(i)
          return
       end if
    end do

  end function next_bucket

  logical function next_item (table)

    ! Set the current location to be the next available item

    type (hash_table_object), intent(inout)   :: table

    ! Local variables
    integer :: i, index, end, oldapos
    type (hash_item), pointer :: hitem ! Temporary hash item 

    ! Point at the current hash item
    hitem => table%ipos%next

    ! While we have a valid key--we can tell if this is the case by 
    ! checking if the key array is associated--keep stepping down
    ! the chain of hash items until we find a match.
    if (associated(hitem%key%bytes)) then
       ! if (debug) print *,'next item ',hitem%hash,' ',achar(hitem%key)
       table%ipos => hitem
       next_item = .TRUE.
       return
    else
       next_item = next_bucket(table)
    end if

  end function next_item

  logical function first_bucket (table)

    ! Set the current location to be the first item in the first bucket
    ! in the hash table

    type (hash_table_object), intent(inout)   :: table

    first_bucket = .FALSE.

    ! Set the index of the current position to the maximum
    ! index
    table%apos = table%mask

    ! Now set it to be the first item in the next bucket. We
    ! know the return value of next should be false as it has
    ! to wrap around the end of the table--otherwise it is an error.
    if (next_bucket(table)) then
       write(stderr,*) table%apos
       write(stderr,*) 'HASH_TABLE :: Cannot find first element in hash'
       stop
    end if

    ! We will return true if our current hash item is not empty
    first_bucket = associated(table%ipos%key%bytes)

  end function first_bucket

  subroutine print_hash_to_stdout (table)

    ! Wrapper for print_hash to allow no unit to be specified and then
    ! it is printed to stdout

    type (hash_table_object), intent(inout) :: table

    call print(stdout,table)

  end subroutine print_hash_to_stdout

  subroutine print_hash (unit, table)

    ! Print out an ascii representation of the hash table -- will
    ! default to stderr if no unit is specified

    integer, intent(in)                     :: unit
    type (hash_table_object), intent(inout) :: table

    ! Local variables
    integer :: oldapos, count

    ! Position ourselves at the start of the hash table
    if (.NOT. first(table)) then
       write(stderr,*) 'HASH_TABLE :: Could not find start of table in print'
       stop
    end if

    oldapos = -1
    count = table%count
    do 
       ! We will know when we have changed 'buckets' because the
       ! current hash table 'index' will have changed
       if ( oldapos /= table%apos ) then
          write (unit,fmt='(I4,"  ")',advance='no') table%apos
          oldapos = table%apos
       else
          write (unit,fmt='(A)',advance='no') '    -> '
       end if
       call print(unit,table%ipos)
       count = count - 1
       if (.NOT. next(table)) exit
    end do

    ! Warn and exit if we missed some of the elements
    if (count /= 0) then
       write(stderr,*) 'Did not print all elements! ', count
       stop
    end if

  end subroutine print_hash

  subroutine print_item_to_stdout (hitem)

    ! Wrapper for print_item to allow no unit to be specified and then
    ! it is printed to stdout

    type (hash_item), intent(in)  :: hitem

    call print(stdout,hitem)

  end subroutine print_item_to_stdout

  subroutine print_item (unit, hitem)

    ! Print out an ascii representation of the hash item -- will
    ! default to stderr if no unit is specified

    integer, intent(in)           :: unit
    type (hash_item), intent(in)  :: hitem

    ! Local variables
    ! character(len=size(hitem%key%bytes)) :: mold
    character(len=1) :: mold

    ! write (unit,'(A,A,I4,A,Z8,A)') transfer(achar(hitem%key%bytes),mold), " => ", hitem%value%bytes, " [", hitem%hash, "]"
    write (unit,'(A,A,I4,A,Z8,A)') transfer(achar(hitem%key%bytes),mold), " => ", &
         transfer(hitem%value%bytes,unit), " [", hitem%hash, "]"
    ! write (unit,'(A,A,A,A,Z8,A)') transfer(achar(hitem%key%bytes),mold), " => ", transfer(achar(hitem%value%bytes),mold), " [", hitem%hash, "]"

  end subroutine print_item

  elemental function one2four(onebyte) result(fourbyte)

    ! Convert a one byte integer to a fourbyte and treating the
    ! one byte integer as unsigned--assumes two's complement storage
    integer(kind=1), intent(in) :: onebyte
    integer(kind=4) :: fourbyte

    fourbyte = int(onebyte,4)
    if (onebyte < 0) fourbyte = fourbyte + 256

  end function one2four

  subroutine mix_32_bit (a, b, c) 

    integer, intent(inout) :: a, b, c
    
    ! From the original C program:
    ! --------------------------------------------------------------------
    ! mix -- mix 3 32-bit values reversibly.
    ! For every delta with one or two bit set, and the deltas of all three
    !   high bits or all three low bits, whether the original value of a,b,c
    !   is almost all zero or is uniformly distributed,
    ! * If mix() is run forward or backward, at least 32 bits in a,b,c
    !   have at least 1/4 probability of changing.
    ! * If mix() is run forward, every bit of c will change between 1/3 and
    !   2/3 of the time.  (Well, 22/100 and 78/100 for some 2-bit deltas.)
    ! mix() was built out of 36 single-cycle latency instructions in a 
    !   structure that could supported 2x parallelism, like so:
    !       a -= b; 
    !       a -= c; x = (c>>13);
    !       b -= c; a ^= x;
    !       b -= a; x = (a<<8);
    !       c -= a; b ^= x;
    !       c -= b; x = (b>>13);
    !       ...
    !   Unfortunately, superscalar Pentiums and Sparcs can't take advantage 
    !   of that parallelism.  They've also turned some of those single-cycle
    !   latency instructions into multi-cycle latency instructions.  Still,
    !   this is the fastest good hash I could find.  There were about 2^^68
    !   to choose from.  I only looked at a billion or so.
    ! --------------------------------------------------------------------
    
    a = a - b; a = a - c; a = ieor(a, ishft(c,-13))
    b = b - c; b = b - a; b = ieor(b, ishft(a, +8))
    c = c - a; c = c - b; c = ieor(c, ishft(b,-13))
    a = a - b; a = a - c; a = ieor(a, ishft(c,-12))
    b = b - c; b = b - a; b = ieor(b, ishft(a,+16))
    c = c - a; c = c - b; c = ieor(c, ishft(b, -5))
    a = a - b; a = a - c; a = ieor(a, ishft(c, -3))
    b = b - c; b = b - a; b = ieor(b, ishft(a,+10))
    c = c - a; c = c - b; c = ieor(c, ishft(b,-15))

  end subroutine mix_32_bit

  integer function lookup (onekey, level)

    ! lookup() -- hash a variable-length key into a 32-bit value
    !   key   : the key (the unaligned variable-length array of bytes)
    !   len   : the length of the key, counting by bytes
    !   level : can be any 4-byte value
    ! Returns a 32-bit value.  Every bit of the key affects every bit of
    ! the return value.  Every 1-bit and 2-bit delta achieves avalanche.
    ! About 6len+35 instructions.
    ! 
    ! The best hash table sizes are powers of 2.  There is no need to do
    ! mod a prime (mod is sooo slow!).  If you need less than 32 bits,
    ! use a bitmask.  For example, if you need only 10 bits, do
    !   h = (h & hashmask(10));
    ! In which case, the hash table should have hashsize(10) elements.
    ! 
    ! If you are hashing n strings (ub1 **)k, do it like this:
    !   for (i=0, h=0; i<n; ++i) h = lookup( k[i], len[i], h);
    ! 
    ! By Bob Jenkins, 1996.  bob_jenkins@burtleburtle.net.  You may use this
    ! code any way you wish, private, educational, or commercial.
    ! 
    ! See http://burtleburtle.net/bob/hash/evahash.html
    ! Use for hash table lookup, or anything where one collision in 2^32 is
    ! acceptable.  Do NOT use for cryptographic purposes.

    integer(kind=1), dimension(:) :: onekey

    integer, intent(in), optional :: level

    integer, dimension(size(onekey)) :: key

    integer :: len, a, b, c, base_index

    ! the golden ratio (9e3779b9)--an arbitrary value
    a = -1640531527;
    b = -1640531527;

    len = size(onekey)

    ! If we are supplied a value for level we set c to this 
    ! (it is typically the previously accessed hash value)
    if (present(level)) then
       c = level
    else
       c = 0
    end if

    ! Convert the 'unsigned' one byte key into a four byte 
    ! version. The one byte values are not 'unsigned' as 
    ! such, but the conversion treats them as such, so 'FF'
    ! is converted to 255 and not -1.
    key = one2four(onekey)

    ! write(*,'(A,I4,4000(:,Z4))') "key: ",len,key

    base_index = 0

    ! Handle most of the key
    do while (base_index + 12 <= len)
      ! write(*,'("a,b,c",3(" ",Z8))') a,b,c
      a = a + key(base_index+1) &
           + ishft(key(base_index+2),+8) &
           + ishft(key(base_index+3),+16) &
           + ishft(key(base_index+4),+24)
      b = b + key(base_index+5) &
           + ishft(key(base_index+6),+8) &
           + ishft(key(base_index+7),+16) &
           + ishft(key(base_index+8),+24)
      c = c + key(base_index+9) &
           + ishft(key(base_index+10),+8) &
           + ishft(key(base_index+11),+16) &
           + ishft(key(base_index+12),+24)
      ! write(*,'("a,b,c",3(" ",Z8))') a,b,c
      call mix(a,b,c)
      base_index = base_index + 12
   end do

   ! handle the last 11 bytes
   c = c + size(key)

   ! Go backwards through the last 11 bytes ...
      
   ! This is the fortran bodge-up of a C-style case statement
   ! which dosen't have breaks, i.e. case 11 also executes
   ! cases 10-1.
   select case(len - base_index)
   case(11) ; goto 11
   case(10) ; goto 10
   case(9)  ; goto  9
   case(8)  ; goto  8
   case(7)  ; goto  7
   case(6)  ; goto  6
   case(5)  ; goto  5
   case(4)  ; goto  4
   case(3)  ; goto  3
   case(2)  ; goto  2
   case(1)  ; goto  1
   end select
      
11 c = c + ishft(key(base_index+11),+24)
10 c = c + ishft(key(base_index+10),+16)
9  c = c + ishft(key(base_index+9),+8)
   ! the first byte of c is reserved for the length
8  b = b + ishft(key(base_index+8),+24)
7  b = b + ishft(key(base_index+7),+16)
6  b = b + ishft(key(base_index+6),+8)
5  b = b + key(base_index+5)
4  a = a + ishft(key(base_index+4),+24)
3  a = a + ishft(key(base_index+3),+16)
2  a = a + ishft(key(base_index+2),+8)
1  a = a + key(base_index+1)
      
   ! write(*,'(A,3(" ",Z8))') 'a,b,c ',a,b,c

   call mix(a,b,c)

   ! report the result
   lookup = c

  end function lookup

end module hash_table
