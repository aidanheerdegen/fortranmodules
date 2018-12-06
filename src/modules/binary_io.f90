module binary_io
 
  !  Buffered reading and writing data to binary files.
  !
  !  Version 1.0, February 1998
  !
  !  Written by Jos Bergervoet

  ! Can't use file_functions as it now uses binary_io. Have to
  ! have local versions of these ...
  ! use file_functions,   only: stderr, freeunit
  use string_functions, only: join, operator(.lcase.)

  implicit none                  ! Check all declarations

  private

  ! Define the unit numbers for standard error, input and output
  integer, parameter ::  stderr = 0, stdin = 5, stdout = 6

  character(len=*), parameter :: version = "$Id: binary_io.f90,v 1.7 2007/07/27 04:03:35 aidan Exp $"

  !! $Log: binary_io.f90,v $
  !! Revision 1.7  2007/07/27 04:03:35  aidan
  !! Removed the dependence on file_functions as that module now uses binary_io.
  !! Unfortunately this means there is now a local version of freeunit and
  !! multiple definitions of stdout, stderr and stdin. The latter can be fixed
  !! by absorbing them into a system_dependent module used by both. The former
  !! will disappear when binary_io is no longer strictly necessary.
  !!
  !! Revision 1.6  2007/07/26 04:39:25  aidan
  !! Fixed bug where reading characters from a filehandle would continue beyond
  !! the end of the file. Was not properly setting error codes when switching to
  !! reading the last record.
  !!
  !! Made the read-ahead hack more efficient by reading this into a separate
  !! NextBuf, so as not to read the same record twice.
  !!
  !! Revision 1.5  2007/07/11 00:40:49  aidan
  !! Added a small hack to fix a problem with the intel compiler. It will allow
  !! a record to be read beyond the end of the file, and will only complain if
  !! the beginning of the record to be read is beyond the end of the file. So
  !! now we read the record beyond the one we want and check it's status. This
  !! module should be obselete with stream io so we don't mind a small hack.
  !!
  !! Revision 1.4  2004/05/10 06:57:35  aidan
  !! HUGE revision. Have not changed the basic type that is read/written, namely
  !! character variables, but have changed almost everything else. No longer
  !! specify unit numbers -- use the initialisation to return a file_handle
  !! object. The advantage of this is we can have different buffer lengths for
  !! different jobs, specified at use rather than compile time. To do this I
  !! had to convert the internal buffer into an array of single characters, which
  !! necessitated some changes in converting from character array to string of
  !! characters.
  !!
  !! Also changed the entire interface to the module. Now have nice simple names
  !! for the routines -- read, write, open, close.
  !!
  !! Rationalised some of the routines as well. Now the character reading thingo
  !! actually reads more than one character at a time.
  !!
  !!
  !! revision 1.1 2004/02/09 00:42:38;  aidan
  !! Initial revision
  !! Module written by Jos Bergervoet (v 1.0, Feb 1998) to do binary 
  !! input/output. Includes some buffering and conversion from raw bytes 
  !! to integer.
  !! ----------------------------
  !! revision 1.2 2004/05/06 06:11:09;  aidan
  !! Added skip_chars routine. This allows one to look ahead in the current
  !! buffer for characters that do not match those supplied in the chars
  !! string. Required for the pnm reading routine, where we need to ignore
  !! whitespace but then treat the subsequent chars as part of an integer
  !! value. It was not good enough to keep popping characters off the buffer
  !! until we came to a non-whitespace.
  !! 
  !! Also added support for big and little endian read/write of integers.
  !! 
  !! As a result of adding in skip_chars, I excerpted the bit of code that
  !! incremented the buffer position and placed this in it's own routine.
  !! Called, ironically enough, inc_bufpos.
  !! ----------------------------
  !! Revision 1.3  2004/05/06 06:36:43  aidan
  !! No real changes -- just updated the indenting for the emacs mode so I can
  !! edit this easily.
  !!

  integer, private, parameter :: MinUnit = 60, MaxUnit = 99, BufLen_default = 256

  type binary_filehandle          ! Can be used outside this unit (is public)
     private                      ! but all fields are unknown outside (private)
     character(len=255) :: Name
     integer            :: UnitNum, FilePos, BufLen, BufPos
     logical            :: ForWriting, EOF
     character(len=1), dimension(:), pointer  :: Buf, NextBuf
  end type binary_filehandle

  ! A nice name for the open statement
  interface open
     module procedure open_for_read_or_write
  end interface

  ! Ditto for close
  interface close
     module procedure close_block_io
  end interface

  ! A generic interface to the various read routines
  interface read
     module procedure char_read, integer_read, integer_block_read
  end interface

  ! A generic interface to the various read routines
  interface write
     module procedure char_block_write, integer_write, integer_block_write
  end interface

  ! A nice name for the routine used to skip characters
  interface skip
     module procedure skip_chars
  end interface

  ! Public types
  public :: binary_filehandle

  ! Public operators
  public :: open, close, read, write, skip

contains

  function init_file_block(file, OpenForWrite, buflen) result(fhandle)

    ! This initialises the binary file handle and returns it as a function result

    character(len=*), intent(in)  :: file
    logical, intent(in)           :: OpenForWrite
    integer, intent(in)           :: buflen
    
    ! Return type
    type(binary_filehandle) :: fhandle

    allocate ( fhandle%Buf(buflen), fhandle%NextBuf(buflen) )

    ! Find a free unit number -- this is not thread safe, i.e. if
    ! something else was opening units at the same time it would
    ! not be ok to get the unit now .. but I guess this is a locking
    ! problem ...
    fhandle % UnitNum = freeunit(MinUnit)
    fhandle % FilePos = 0
    fhandle % Name = file
    fhandle % BufLen = buflen
    fhandle % BufPos = 0
    fhandle % ForWriting = OpenForWrite
    fhandle % EOF = .false.

    return
  end function init_file_block

  subroutine open_for_read_or_write(fhandle, file, action, buflen)

    ! Opens a file handle for read or write access, depending on the
    ! action keyword. Note also that it is possible to optionally
    ! define the buffer length -- it is desirable when reading and
    ! writing very large files to have a larger buffer.

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    character(len=*), intent(in)           :: file
    character(len=*), intent(in), optional :: action
    integer, intent(in), optional          :: buflen

    ! Local variables
    character(len=5) :: myaction
    character(len=7) :: mystatus
    logical :: writing
    integer :: i, buffer_length

    ! Use the default buffer length unless we specified one
    buffer_length = BufLen_default 
    if (present(buflen)) buffer_length = buflen

    ! Default to read unless action is specified
    myaction = 'read'
    if (present(action)) myaction = .lcase. action

    select case(myaction)
    case ('read')
       mystatus = 'old'
       writing = .FALSE.
    case ('write')
       mystatus = 'replace'
       writing = .TRUE.
    case default
       write(stderr,*) 'BINARY_IO :: Invalid action ',trim(action),', defaulting to read only'
       do i = 1, len(myaction)
          print '(I2,Z8)',i,ichar(myaction(i:i))
       end do
       mystatus = 'old'
       writing = .FALSE.
    end select
    
    ! Initialise the file handle object
    fhandle = init_file_block(file, OpenForWrite=writing, buflen=buffer_length)

    ! Open the file
    open (unit=fhandle%UnitNum, file=file, form="unformatted", access="direct", &
        recl=fhandle%BufLen, status=mystatus, action=myaction)

    return
  end subroutine open_for_read_or_write

  subroutine flush_last_rec(fhandle)           

    ! Write (last) record < BufLen

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle

    ! Local variables
    integer :: fp, i

    ! close, then re-open with recl=1
    close (unit=fhandle%UnitNum)
    open (unit=fhandle%UnitNum, file=fhandle%Name, form="unformatted", access="direct", &
         recl=1, status="old", action="write")
    
    fp = fhandle % FilePos * fhandle % BufLen
    do i = 1,fhandle % BufPos
       write(unit=fhandle%UnitNum,rec=fp+i) fhandle % buf(i)
    end do
    return
  end subroutine flush_last_rec
  
  subroutine close_block_io(fhandle)

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle

    if (fhandle%ForWriting .and. fhandle%BufPos > 0) then
       call flush_last_rec(fhandle)
    end if
    close(unit=fhandle%UnitNum)
    deallocate ( fhandle%Buf, fhandle%NextBuf )
    return
  end subroutine close_block_io
  
  subroutine char_block_write(fhandle, s)

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    character(len=*), intent(in)           :: s

    ! Local variables
    integer :: ls, ps, fp, bp, slice
    
    fp = fhandle%FilePos
    bp = fhandle%BufPos
    ls = len(s)
    ps = 0
    do  ! write the string in slices of our buffersize BufLen
       if (ls <= ps) then
          exit
       else if (ls-ps > fhandle % BufLen-bp) then
          slice = fhandle % BufLen-bp
       else
          slice = ls-ps
       end if
       ! fhandle % buf(bp+1 : bp+slice) = s(ps+1 : ps+slice)
       call copy_string_to_buffer(fhandle % buf(bp+1 : bp+slice), s(ps+1 : ps+slice))
       bp = bp+slice
       ps = ps+slice
       if (bp == fhandle % BufLen) then
          write(unit=fhandle%UnitNum, rec=fp+1) fhandle%buf
          fp = fp+1
          bp = 0
       end if
    end do
    fhandle%FilePos = fp
    fhandle%BufPos = bp
    return
  end subroutine char_block_write
  
  subroutine read_last_rec(fhandle, NewPtr, ErrCod)  ! Last record < BufLen
    
    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    integer, intent(out)                   :: NewPtr, ErrCod

    ! Local variables
    integer :: fp, bp, Nfound
    
    ! We'll close and re-open with recl=1 to read the last few bytes one by one.
    close (unit=fhandle%UnitNum)
    open (unit=fhandle%UnitNum,file=fhandle%Name,form="unformatted",access="direct", &
         recl=1, status="old", action="read")
    fp = fhandle % FilePos * fhandle % BufLen
    
    Nfound = 0
    do bp = 1,fhandle % BufLen
       read(unit=fhandle%UnitNum,rec=fp+bp,iostat=ErrCod) fhandle % buf(bp)
       if (ErrCod /= 0) then
          Nfound = bp-1
          exit
       end if
       Nfound = bp
    end do
    if (Nfound > 0) then      ! at least we found some valid bytes
       ErrCod = 0
    end if
    ! Now we shift the bytes to high side of buffer
    NewPtr = fhandle % BufLen-Nfound+1
    fhandle%Buf (NewPtr:fhandle % BufLen) = fhandle%Buf (1:Nfound)
    fhandle%EOF = .true.
    return
  end subroutine read_last_rec

  subroutine char_read(fhandle, c, ReadErr)

    ! Read a character string from the file handle

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    character(len=*), intent(out)          :: c
    integer, intent(out), optional         :: ReadErr

    ! Local variables
    integer :: i, fp, bp, ErrCod
    
    do i = 1, len(c)
    
       bp = inc_bufpos(fhandle, ErrCod)

       if (ErrCod /= 0) exit
       
       c(i:i) = fhandle % buf(bp)
       fhandle % BufPos = modulo(bp+1, fhandle % BufLen+1)

    end do
  
    if (present(ReadErr)) then
       ReadErr = ErrCod
    else if (ErrCod > 0) then
       write(unit=*,fmt=*) "Error reading binary file ", trim(fhandle%Name)
       stop
    end if
    
    return

  end subroutine char_read

  function inc_bufpos(fhandle, ErrCod) result(bp)

    ! It is somewhat ironic that inc_bufpos does not in fact increment the
    ! buffer position index. This is due to the fact that the engine of
    ! reading (char_read) increments the buffer after it has finished
    
    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    integer, intent(out)                   :: ErrCod
    
    ! Return value of function
    integer :: bp
  
    ! Local variables
    integer :: fp
    
    ErrCod = 0
    
    bp = fhandle % BufPos
    if (bp == 0) then
       ! Need to check if we've already read the last record and now we've got to the 
       ! end of the file
       if (fhandle % EOF) then
          ErrCod = -1
       else
          fp = fhandle % FilePos
          ! We have already read in the next buffer so assign it to the current buffer
          fhandle % Buf = fhandle % NextBuf
          ! This is a dodgy hack. The intel fortran compiler will let you read a record
          ! even if that takes it beyond the end of the file. It will behave correctly
          ! if the *start* of a record is beyond the end of a file. So, we read the record
          ! beyond the one we are interested in, check that error, and if it is ok we read
          ! in the record we want, if not, we handball it on to read_last_rec, which does
          ! a character by character read, which should be ok. The problem with this is if
          ! we get a read error on the second read, and we're not checking for that. We
          ! shouldn't get the same error as the first time, as that record should exist. 
          read(unit=fhandle%UnitNum,rec=fp+2,iostat=ErrCod) fhandle % NextBuf
          if (ErrCod == 0) then
             ! First time through we have to read in the record
             if (fp == 0) read(unit=fhandle%UnitNum,rec=fp+1,iostat=ErrCod) fhandle % Buf
             bp = 1
          else
             call read_last_rec(fhandle, bp, ErrCod)
          end if
          fhandle % FilePos = fp+1
       end if
    end if
    
  end function inc_bufpos

  subroutine integer_read(fhandle, iread, length, littleendian)
    
    ! Read a single integer from the filehandle

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    integer, intent(out)                   :: iread
    integer, intent(in)                    :: length
    logical, intent(in), optional          :: littleendian

    ! Local variables
    integer :: i, icount, start, end, increment
    character(len=1) :: charbuf(length)
    logical :: bigendian
    
    bigendian = .TRUE.
    
    if (present(littleendian)) bigendian = .NOT. littleendian

    if (bigendian) then
       start = length 
       end = 1
       increment = -1
    else
       start = 1 
       end = length
       increment = 1
    end if
  
    do i=1,length
       call char_read(fhandle, charbuf(i))
    end do
    
    iread = 0
    icount = 0
    do i=start,end,increment
       iread = iread + ichar(charbuf(i))*(256**(icount))
       icount = icount + 1
    end do
  
    return
  end subroutine integer_read

  subroutine integer_block_read(fhandle, iread, length, littleendian)

    ! Read an array of integers from the file handle

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    integer, dimension(:), intent(out)     :: iread
    integer, intent(in)                    :: length
    logical, intent(in), optional          :: littleendian

    ! Local variables
    integer :: i
    logical :: bigendian
    
    bigendian = .TRUE.
    
    if (present(littleendian)) bigendian = .NOT. littleendian

    do i = 1, size(iread)
       call read(fhandle, iread(i), length, .NOT. bigendian)
    end do

    return
  end subroutine integer_block_read

  subroutine integer_write(fhandle, iwrite, length, littleendian)

    ! Write an integer to the filehandle
    
    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    integer, intent(in)             :: iwrite, length
    logical, intent(in), optional   :: littleendian

    integer :: i, icount, start, end, increment
    logical :: bigendian
    character(len=length) :: charbuf
    
    bigendian = .TRUE.

    if (present(littleendian)) bigendian = .NOT. littleendian
    
    if (bigendian) then
       start = length 
       end = 1
       increment = -1
    else
       start = 1 
       end = length
       increment = 1
    end if

    icount = 0
    do i=start,end,increment
       charbuf(i:i) = char(modulo( iwrite/(256**(icount)), 256))
       icount = icount + 1
    end do
    
    call char_block_write(fhandle, charbuf)
    return
  end subroutine integer_write

  subroutine integer_block_write(fhandle, iwrite, length, littleendian)

    ! Write an array of integers to the file handle

    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    integer, intent(in)                    :: iwrite(:), length
    logical, intent(in), optional          :: littleendian

    ! Local variables
    logical :: bigendian
    integer :: i

    bigendian = .TRUE.

    if (present(littleendian)) bigendian = .NOT. littleendian

    do i = 1, size(iwrite)
       call write(fhandle, iwrite(i), length, .NOT. bigendian)
    end do
    
    return
  end subroutine integer_block_write
  
  subroutine skip_chars (fhandle, chars)

    ! Skip over any of the specified characters. Read in a new 
    ! buffer if necessary
    
    ! Interface variables
    type(binary_filehandle), intent(inout) :: fhandle
    character(len=*), intent(in)           :: chars
    
    ! Local variables
    integer :: ErrCod, oldbp
    integer :: bp
    
    ErrCod = 0
    
    do 
       ! Get the current buffer position
       bp = inc_bufpos(fhandle, ErrCod)
       if (ErrCod /=0) then
          ! Exit from the routine if we've reached the end of the file. We
          ! have, by definition, skipped all the space, leave it up to any
          ! subsequent reading routine to flag an error.
          exit
       end if
       ! Scan the buffer for characters that don't match our specified set
       ! and set the buffer pointer to this value
       bp = verify(join(fhandle % buf(bp:)),chars)
       if (bp == 0) then
          ! We haven't found any non-match characters we'll have to read in
          ! a new buffer at this point (the buffer pointer will be zero so
          ! this will force a buffer read)
          fhandle % BufPos = bp
          ! Scan the next buffer
          cycle
       end if
       ! We found non-match characters, and the buffer pointer is now set
       ! to this location, so any subsequent read will get this character
       fhandle % BufPos = fhandle % BufPos + bp - 1
       exit
    end do
    
  end subroutine skip_chars

  subroutine copy_string_to_buffer (buffer, string)

    ! Just a small routine to copy the contents of a character string into
    ! a character buffer (array of characters of length 1)

    ! Interface variables
    character(len=1), dimension(:), intent(out) :: buffer
    character(len=*), intent(in)                :: string
    
    ! Local variables
    integer :: i

    do i = 1, min(len(string),size(buffer))
       buffer(i) = string(i:i)
    end do

  end subroutine copy_string_to_buffer

  ! Included a version of freeunit in this module so we don't have
  ! to rely on file_functions, which uses this module. This results
  ! in an undesirable cross dependence

  integer function freeunit(start)

    ! This was nicked from the 'Tinker' package of Jay William Ponder
    !
    !     "freeunit" finds an unopened Fortran I/O unit and returns
    !     its numerical value from 1 to 99; the units already assigned
    !     to stdin and stdout (usually 5 and 6) are skipped since
    !     they have special meaning as the default I/O units

    ! Have an optional "start" value from which we will start to look for
    ! free units
    integer, intent(in), optional :: start
    
    logical used
      
    ! try each logical unit until an unopened one is found

    if (present(start)) then
       freeunit = start - 1
    else
       freeunit = 0
    end if

    used = .true.
    do while (used)
       freeunit = freeunit + 1
       if (freeunit /= stdin .and. freeunit /= stdout) then
          if (freeunit .gt. 99) then
             write (stderr,*) ' FILE_FUNCTIONS :: freeunit -- no available I/O units'
             stop
          end if
          inquire (unit=freeunit,opened=used)
       end if
    end do
    
  end function freeunit

end module binary_io
