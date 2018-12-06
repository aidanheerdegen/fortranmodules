module file_functions

  use iso_varying_string, only: trim, var_str
  use string_functions, only: operator(.lcase.), join
  use binary_io, only: openbin => open, read, write, close, binary_filehandle

  implicit none

  private

  ! This module provides some routines and parameters that are
  ! useful when using I/O operations in fortran. There is a
  ! procedure to automatically determine the lowest unused unit
  ! number. Functions to determine if a file exists, and another
  ! to count the lines in a file. A wrapper to the open function 
  ! that returns a unit number and also performs some sanity checks.
  ! There is also a routine which reads a buffer from an open unit,
  ! and can filter out commented and blank lines.

  !! $Log: file_functions.f90,v $
  !! Revision 1.14  2007/11/13 00:53:45  aidan
  !! Updated readbuffer routines to not use optional arguments in specification
  !! as well (should have done this with the last commit .. doh!)
  !!
  !! Revision 1.13  2007/11/13 00:40:30  aidan
  !! Added extra routines for commented and uncommented line counting. It is
  !! not standard conforming to have optional dummy arguments in specification
  !! expressions (i.e. commentchar(len(comment)) when comment is optional).
  !!
  !! Revision 1.12  2007/07/26 05:02:42  aidan
  !! Major update. Made the comment character optional in read_buffer. Removed
  !! a couple of interface routines as a result.
  !!
  !! Added 'inline' option to read_buffer to delete inline comments. This is
  !! useful, as inline comments are nice and can really stuff up input file
  !! if the user assumes they are ok (which is fair enough as they are quite
  !! ubiquitous).
  !!
  !! Added full read_buffer options to line_count, so we can count lines in a
  !! file on the same basis we will read them in. This is important for
  !! parsing files which contain comments etc.
  !!
  !! Revision 1.11  2007/07/11 00:36:14  aidan
  !! Added a move_file routine and added an option to force a copy in the
  !! copy_file routine even if the destination file already exists.
  !!
  !! Revision 1.10  2006/07/21 05:48:55  aidan
  !! *** empty log message ***
  !!
  !! Revision 1.9  2006/05/19 05:36:46  aidan
  !! Added a simple logging routine.
  !!
  !! Revision 1.8  2006/04/13 04:43:11  aidan
  !! Added the 'touch' and 'semaphore' functions. The former sort of emulates
  !! the unix touch command but, at the moment, does not update the
  !! modification time. The latter allows for file semaphores,
  !! i.e. appropriately named files in a directory which signal a program. By
  !! default if semaphore finds the file it will delete it and return TRUE,
  !! otherwise false.
  !!
  !! Added the option to 'append' when using the open statement. This is
  !! needed for log files. I should write a 'log' command, to automatically
  !! open with append, write and then close to allow convenient logging.
  !!
  !! Commented out a redundant 'rewind' command from the open statement. At
  !! least I hope it is redundant!
  !!
  !! Revision 1.7  2006/03/30 03:46:40  aidan
  !! Changed behaviour of line_count so that it rewinds at the end rather 
  !! than closing the file. The previous behaviour was surely an unnoticed bug.
  !!
  !! Changed behaviour of file_open. Now when we specify a status but not 
  !! 'opened_ok' the routine will squawk loudly about any problems and STOP the 
  !! program. Previously it would return a unit number of -1 and confuse the hell 
  !! out of the monkey on the other end of the keyboard.
  !!
  !! Revision 1.6  2005/01/20 02:18:32  aidan
  !! Revisited the comment char problem. Multiple characters are very useful as
  !! we can specify more than one valid comment char at once. So I have split
  !! the functionality up further, and made specific routines to cope with no
  !! comment char and made comment a required argument in the main routines.
  !!
  !! Revision 1.5  2005/01/19 23:31:48  aidan
  !! Changed the commentchar in the readbuffer routines to be a fixed length (=1).
  !! This fixed a bug where we used the length of this optional argument in a
  !! specification expression. This is a violation of the standard.
  !!
  !! Revision 1.4  2004/06/10 05:20:40  aidan
  !! Updated delete function to take multiple arguments.
  !!
  !! Revision 1.3  2004/06/04 02:08:04  aidan
  !! Added a routine for deleting files. Surprisingly enough, it is called delete.
  !!
  !! Revision 1.2  2004/05/06 06:46:24  aidan
  !! Added an optional 'start' parameter to freeunit, to allow it to start
  !! searching for free units from a user specified unit number. This is
  !! potentially useful for avoiding conflicts when used in a program with
  !! hard-wired unit numbers known to be below a certain value,
  !!
  !! Revision 1.1  2003/09/10 05:03:41  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: file_functions.f90,v 1.14 2007/11/13 00:53:45 aidan Exp $"

  ! Define the unit numbers for standard error, input and output
  integer, parameter ::  stderr = 0, stdin = 5, stdout = 6

  ! Function to test if a file exists
  interface exists
     module procedure file_exists
  end interface

  ! Generic file opening function (returns a unit number)
  interface open
     module procedure file_open
  end interface

  ! Function to count the number of lines in a file
  interface line_count
     module procedure count_lines_in_file, open_then_count, &
          count_lines_in_file_nocomment, open_then_count_nocomment
  end interface

  ! A standardised routine for reading into a character buffer -- takes care
  ! of filtering out commented and blank lines
  interface read_buffer
     module procedure readbuffer, open_then_readbuffer, readbuffer_nocomment, open_then_readbuffer_nocomment
  end interface

  ! Provide a routine for deleting files
  interface delete
     module procedure delete_file, delete_files
  end interface

  ! Provide a routine for copying files
  interface copy
     module procedure copy_file
  end interface

  ! Provide a routine for moving files
  interface move
     module procedure move_file
  end interface

  ! Provide a routine for making a file if it doesn't exist (currently
  ! does not update modification time of existing file like unix touch command)
  interface touch
     module procedure touch_file
  end interface

  ! Provide a routine for semaphore files. Where you check for the existence
  ! of the file, returning true or false, and by default the semaphore is deleted
  interface semaphore
     module procedure semaphore_file
  end interface

  ! Provide a routine for simple logging to a file
  interface log
     module procedure log_to_file
  end interface

  logical, public :: debug = .FALSE.

  ! Public variables
  public :: stderr, stdin, stdout

  ! Public routines
  public :: freeunit, line_count, read_buffer
  public :: exists, open, delete, touch, semaphore, log, copy, move

contains

  integer function count_lines_in_file(unit, comment, removeblanks) result(nlines)
    
    ! Count the number of lines in a file and return this value

    ! Interface variables
    integer, intent(in)           :: unit
    character(len=*), intent(in)  :: comment
    logical, intent(in), optional :: removeblanks

    ! Local variables
    logical :: finished, noblanks
    character(len=1000) :: buffer

    ! Default to not removing blank lines
    if (present(removeblanks)) then
       noblanks = removeblanks
    else
       noblanks = .false.
    end if

    nlines = 0
    do
       call read_buffer(unit, buffer, finished, comment, noblanks)
       if (finished) exit
       nlines = nlines + 1
    end do

    rewind(unit)
    
  end function count_lines_in_file

  integer function count_lines_in_file_nocomment(unit, removeblanks) result(nlines)
    
    ! Count the number of lines in a file and return this value

    ! Interface variables
    integer, intent(in)           :: unit
    logical, intent(in), optional :: removeblanks

    ! Local variables
    logical :: finished, noblanks, removecomments 
    character(len=1000) :: buffer

    ! Default to not removing blank lines
    if (present(removeblanks)) then
       noblanks = removeblanks
    else
       noblanks = .false.
    end if

    nlines = 0
    do
       call read_buffer(unit, buffer, finished, removeblanks=noblanks)
       if (finished) exit
       nlines = nlines + 1
    end do

    rewind(unit)
    
  end function count_lines_in_file_nocomment

  integer function open_then_count(file, comment, removeblanks) result(nlines)

    ! This is a small wrapper that allows the user to specify a
    ! filename which this routine will open for them, count the
    ! number of lines and then close it again.

    ! filename
    character(len=*), intent(in)  :: file
    character(len=*), intent(in)  :: comment
    logical, intent(in), optional :: removeblanks

    ! Internal variables
    logical :: opened_ok, noblanks
    integer :: unit

    ! Default to not removing blank lines
    if (present(removeblanks)) then
       noblanks = removeblanks
    else
       noblanks = .false.
    end if

    ! Open the file
    unit = open(file, opened_ok, status='old')
    
    if (opened_ok) then
       nlines = line_count(unit, comment, noblanks)
       close(unit)
    else
       write (stderr,*)'FILE_FUNCTIONS :: Unable to find file: ',file
       stop
    end if

  end function open_then_count

  integer function open_then_count_nocomment(file, removeblanks) result(nlines)

    ! This is a small wrapper that allows the user to specify a
    ! filename which this routine will open for them, count the
    ! number of lines and then close it again.

    ! filename
    character(len=*), intent(in) :: file
    logical, intent(in), optional          :: removeblanks

    ! Internal variables
    logical :: opened_ok, noblanks
    integer :: unit

    ! Default to not removing blank lines
    if (present(removeblanks)) then
       noblanks = removeblanks
    else
       noblanks = .false.
    end if

    ! Open the file
    unit = open(file, opened_ok, status='old')
    
    if (opened_ok) then
       nlines = line_count(unit, removeblanks=noblanks)
       close(unit)
    else
       write (stderr,*)'FILE_FUNCTIONS :: Unable to find file: ',file
       stop
    end if

  end function open_then_count_nocomment

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

  integer function file_open(filename, opened_ok, status, append)

    ! A generic interface for opening files. Why? Well we don't have
    ! to worry about units (it is returned to us) and we can examine
    ! the return status to see if the file open worked -- saves 
    ! checking for existence etc

    character(len=*), intent(in)   :: filename
    logical, optional, intent(out) :: opened_ok
    character(len=*), optional, intent(in) :: status
    logical, optional, intent(in) :: append

    ! Local variables
    character(len=7) :: mystatus
    logical :: filename_exists, myappend

    ! Define a default return value for the function should an error
    ! occur
    file_open = -1
    
    ! Define a default status and change it if it has been specified
    ! in the subroutine parameters (make it lower case also)
    mystatus = 'unknown'
    if (present(status)) mystatus = .lcase. status

    ! Check if we're appending to the file
    myappend = .FALSE.
    if (present(append)) myappend = append

    ! We initialise this to false (if it is present as a parameter)
    ! as we will make our tests at the start of the subroutine and
    ! return if we encounter a problem
    if (present(opened_ok)) opened_ok = .false.

    ! We require a filename (ifc would prompt for one in this instance
    ! if we didn't trap it out here)
    if (trim(filename) == '') then
       if (.NOT. present(opened_ok)) then
          write(stderr,*) 'FILE_FUNCTIONS :: No filename specified!'
       else
          return
       end if
    end if

    filename_exists = exists(filename)

    ! We will require that an old file exists, a new file dosen't exist
    if ( (mystatus.eq.'old' .and. .not. filename_exists) .or. &
         (mystatus.eq.'new' .and. filename_exists) ) then
       if (present(opened_ok)) then
          return
       else
          if (filename_exists) then
             write(stderr,*) 'FILE_FUNCTIONS :: ',trim(filename),' already exists!'
          else
             write(stderr,*) 'FILE_FUNCTIONS :: Cannot open "',trim(filename),'". It does not exist!'
          end if
          STOP
       end if
    end if
       
    ! Grab a free unit number
    file_open = freeunit()

    ! Attempt to open the file -- if we get an error we exit without
    ! setting the value of opened_ok to .TRUE.
    if (myappend) then
       open(unit=file_open, file=filename, status=mystatus, err=1, position='append')
    else
       open(unit=file_open, file=filename, status=mystatus, err=1)
       ! rewind(file_open)
    end if

    if (present(opened_ok)) opened_ok = .TRUE.

1   return

  end function file_open

  subroutine delete_files(filenames, status)

    ! Interface variables
    character(len=*), dimension(:), intent(in) :: filenames
    integer, intent(out), optional             :: status

    ! Local variables
    integer :: i, tmp_status, local_status

    local_status = 0

    do i = 1, size(filenames)
       call delete(trim(filenames(i)),tmp_status)
       local_status = local_status + abs(tmp_status)
    end do

    if (present(status)) then
       ! The caller has specified they want the status returned to then, 
       ! so do so ...
       status = local_status
    else
       ! Better warn the caller that there was a problem with deleting the
       ! files, but we won't stop execution.
       if (local_status /= 0) then
          write(stderr,*) 'FILE_FUNCTIONS :: Errors occurred when deleting files: ',join(trim(var_str(filenames))," ")
       end if
    end if

  end subroutine delete_files

  subroutine delete_file(filename, status)

    ! This provides a fortran way of deleting a file, without
    ! having to use platform specific functions. It also 
    ! provides a little sanity checking.

    ! Argument variables 
    character(len=*), intent(in)  :: filename
    integer, intent(out), optional :: status

    ! Local variables 
    integer :: local_status, unit
    logical :: already_open

    local_status = 1

    ! Make sure the file exists before we try and delete it
    if (exists(filename)) then
       inquire (file=filename, opened=already_open)
       if (already_open) then
          ! If it is already open find the unit number associated with it
          inquire (file=filename, number=unit)
          local_status = 0
       else
          unit = freeunit()
          ! Attempt to open the file
          open(unit=unit, file=filename, iostat=local_status)
       end if
       if (local_status == 0) then
          ! We have an open file with a known unit number, so now we
          ! can delete it
          close(unit, status='delete')
       end if
    end if

    if (present(status)) then
       ! The caller has specified they want the status returned to then, 
       ! so do so ...
       status = local_status
    else
       ! Better warn the caller that there was a problem with deleting the
       ! file, but we won't stop execution.
       if (local_status /= 0) then
          write(stderr,*) 'FILE_FUNCTIONS :: Could not delete file: ',filename
       end if
    end if

  end subroutine delete_file

  subroutine copy_file(source, destination, force, status)

    ! This provides a fortran way of copying a file, without
    ! having to use platform specific functions. It also 
    ! provides a little sanity checking.

    ! Argument variables 
    character(len=*), intent(in)   :: source, destination
    logical, intent(in), optional  :: force
    integer, intent(out), optional :: status

    ! Local variables 
    integer :: local_status, error, unit
    logical :: already_open, force_copy
    type (binary_filehandle) :: in, out
    character(len=1) :: c

    local_status = 0

    ! See if we want to force a copy even if the destination file
    ! already exists
    force_copy = .false.
    if (present(force)) then
       force_copy = force
    end if

    ! Make sure the file exists before we try and copy it
    if (.not. exists(source)) then
       write(stderr,*) 'FILE_FUNCTIONS :: Cannot copy file, it does not exist: ',trim(source)
       local_status = 1
       if (present(status)) status = local_status
       return
    end if

    inquire (file=source, opened=already_open)
    if (already_open) then
       ! If it is already open find the unit number associated with it and close it
       inquire (file=source, number=unit)
       close(unit)
    end if

    ! Check if the destination file exists before we try and copy it
    if (exists(destination)) then
       if (force_copy) then
          inquire (file=destination, opened=already_open)
          if (already_open) then
             ! If it is already open find the unit number associated with it and close it
             inquire (file=destination, number=unit)
             close(unit)
          end if
       else
          write(stderr,*) 'FILE_FUNCTIONS :: Cannot copy file, destination file already exists: ',trim(destination)
          local_status = 1
          if (present(status)) status = local_status
          return
       end if
    end if

    ! Open source file
    call openbin(in, source, 'read', 1024)

    ! Open destination file
    call openbin(out, destination, 'write', 1024)

    call read(in, c, error)
    do while (error == 0)
       if (debug) print *,c
       call write(out,c)
       call read(in, c, error)
    end do

    call close(in)
    call close(out)

    if (present(status)) status = local_status

  end subroutine copy_file

  subroutine move_file(source, destination, force, status)

    ! This provides a fortran way of moving a file, without
    ! having to use platform specific functions. It also 
    ! provides a little sanity checking.

    ! Argument variables 
    character(len=*), intent(in)   :: source, destination
    logical, intent(in), optional  :: force
    integer, intent(out), optional :: status

    ! Local variables 
    integer :: local_status, tmp_status
    logical :: force_move

    ! See if we want to force a move even if the destination file
    ! already exists
    force_move = .false.
    if (present(force)) then
       force_move = force
    end if

    ! Make sure the file exists before we try and move it
    if (.not. exists(source)) then
       write(stderr,*) 'FILE_FUNCTIONS :: Cannot move file, it does not exist: ',trim(source)
       local_status = 1
    else
       call copy(source, destination, force_move, tmp_status)
       local_status = tmp_status
       call delete(source, tmp_status)
       local_status = local_status + tmp_status
    end if

    if (present(status)) status = local_status

  end subroutine move_file

  logical function file_exists(filename)

    character(len=*), intent(in) :: filename

    ! We make the decision that an empty filename cannot exist as
    ! the inquire function barfs badly on it
    if (trim(filename) == '') then
       file_exists = .FALSE.
       return
    end if

    inquire(file=filename, exist=file_exists)

  end function file_exists
  
  subroutine readbuffer (unit, buffer, finished, comment, removeblanks, inline)

    ! Reads a line from the given unit into the buffer supplied. Lines that
    ! begin with a comment character (defaults to '#'), and optionally
    ! blank lines, are ignored. When an end of file is reached the logical
    ! variable finished is set to true

    integer, intent(in)           :: unit
    character(len=*), intent(out) :: buffer
    logical, intent(out)          :: finished
    character(len=*), intent(in)  :: comment
    logical, intent(in), optional :: removeblanks, inline

    logical :: removecomments, noblanks, removeinline
    character(len=1) :: firstchar
    integer :: i

    noblanks = .FALSE.
    if (present(removeblanks)) noblanks = removeblanks

    removeinline = .FALSE.
    if (present(inline)) removeinline = inline

    buffer = ''
    ! firstchar = comment(1:1)
    finished = .FALSE.
    
    do 
       read(unit,fmt='(A)',end=1) buffer
       ! Set firstchar to be the first non-blank character in the buffer
       firstchar = adjustl(buffer)
       ! Read another line if we have a line whose first non-blank
       ! character is a quote character or if the entire line is blank 
       ! and we are discarding blank lines
       if ( verify(firstchar, comment) == 0 & 
            .or. (noblanks .and. buffer == ' ') ) cycle
       ! Otherwise we exit this loop
       exit
    end do

    if (removeinline) then
       i = scan(buffer, comment)
       if (i > 0) buffer = buffer(1:i-1)
    end if

    ! Note that we have a different return status if we fall through to
    ! here than if we reach the end of the file, where we skip this
    ! return and set the variables below
    return
    
1   finished = .TRUE.
    buffer = ''

  end subroutine readbuffer

  subroutine readbuffer_nocomment (unit, buffer, finished, removeblanks)

    ! Reads a line from the given unit into the buffer supplied. Optionally 
    ! blank lines, are ignored. When an end of file is reached the logical
    ! variable finished is set to true

    integer, intent(in)                    :: unit
    character(len=*), intent(out)          :: buffer
    logical, intent(out)                   :: finished
    logical, intent(in), optional          :: removeblanks

    logical :: noblanks
    character(len=1) :: firstchar
    integer :: i

    noblanks = .FALSE.
    if (present(removeblanks)) noblanks = removeblanks

    buffer = ''
    ! firstchar = comment(1:1)
    finished = .FALSE.
    
    do 
       read(unit,fmt='(A)',end=1) buffer
       ! Set firstchar to be the first non-blank character in the buffer
       firstchar = adjustl(buffer)
       ! Read another line if entire line is blank and we are discarding blank lines
       if ((noblanks .and. buffer == ' ') ) cycle
       ! Otherwise we exit this loop
       exit
    end do

    ! Note that we have a different return status if we fall through to
    ! here than if we reach the end of the file, where we skip this
    ! return and set the variables below
    return
    
1   finished = .TRUE.
    buffer = ''

  end subroutine readbuffer_nocomment

  subroutine open_then_readbuffer(file, buffer, finished, comment, removeblanks, inline)

    ! This is a small wrapper that allows the user to specify a
    ! filename which this routine will open for them. Each call
    ! to this routine will call read_buffer and and then close 
    ! when it is finished.

    ! filename
    character(len=*), intent(in)  :: file
    character(len=*), intent(out) :: buffer
    logical, intent(out)          :: finished
    character(len=*), intent(in)  :: comment
    logical, intent(in), optional :: removeblanks, inline

    ! Internal variables
    logical :: opened_ok, already_open
    logical :: noblanks, removeinline
    integer :: unit
    character(len=1) :: firstchar

    ! Open the file if we haven't already done so
    inquire (file=file, opened=already_open)
    if (already_open) then
       inquire (file=file, number=unit)
       opened_ok = .TRUE.
    else
       unit = open(file, opened_ok, status='old')
    end if

    if (opened_ok) then

       noblanks = .FALSE.
       if (present(removeblanks)) noblanks = removeblanks
       removeinline = .FALSE.
       if (present(inline)) removeinline = inline

       call read_buffer(unit, buffer, finished, comment, noblanks, removeinline)

       if (finished) close(unit)

    else
       write (stderr,*)'FILE_FUNCTIONS :: Unable to find file: ',file
       stop
    end if

  end subroutine open_then_readbuffer

  subroutine open_then_readbuffer_nocomment(file, buffer, finished, removeblanks)

    ! This is a small wrapper that allows the user to specify a
    ! filename which this routine will open for them. Each call
    ! to this routine will call read_buffer and and then close 
    ! when it is finished.

    ! filename
    character(len=*), intent(in)           :: file
    character(len=*), intent(out)          :: buffer
    logical, intent(out)                   :: finished
    logical, intent(in), optional          :: removeblanks

    ! Internal variables
    logical :: opened_ok, already_open
    logical :: noblanks
    integer :: unit
    character(len=1) :: firstchar

    ! Open the file if we haven't already done so
    inquire (file=file, opened=already_open)
    if (already_open) then
       inquire (file=file, number=unit)
       opened_ok = .TRUE.
    else
       unit = open(file, opened_ok, status='old')
    end if

    if (opened_ok) then

       noblanks = .FALSE.
       if (present(removeblanks)) noblanks = removeblanks

       call read_buffer(unit, buffer, finished, removeblanks=noblanks)

       if (finished) close(unit)

    else
       write (stderr,*)'FILE_FUNCTIONS :: Unable to find file: ',file
       stop
    end if

  end subroutine open_then_readbuffer_nocomment

  subroutine touch_file(file)
    
    character(len=*), intent(in) :: file

    integer :: unit

    unit = open(file)
    close(unit)

  end subroutine touch_file

  logical function semaphore_file(file, deleteafter)

    character(len=*), intent(in)  :: file
    logical, intent(in), optional :: deleteafter

    ! Local variables
    logical :: delete_file

    delete_file = .TRUE.

    if (present(deleteafter)) delete_file = deleteafter
    
    if (exists(file)) then
       semaphore_file = .TRUE.
       if (delete_file) call delete(file)
    else
       semaphore_file = .FALSE.
    end if

  end function semaphore_file

  subroutine log_to_file(file, message)

    ! Convenience routine which opens a file at the end, appends
    ! a message and then closes the file

    character(len=*), intent(in) :: file, message
    
    integer :: unit
    
    unit = open(trim(file), append=.TRUE.)

    write(unit,'(A)') message

    close(unit)

  end subroutine log_to_file

end module file_functions
