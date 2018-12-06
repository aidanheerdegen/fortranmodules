module pnm_class

  use binary_io
  use file_functions, only: stdout, stderr, freeunit, exists
  use string_functions, only: ucase
  use image_transforms, only: rotate_image, rotave_image
  use precision, only: rd_kind, rs_kind

  implicit none

  private

  character(len=*), parameter :: version = "$Id: pnm_class.f90,v 1.7 2005/09/26 02:12:12 aidan Exp $"

  !! $Log: pnm_class.f90,v $
  !! Revision 1.7  2005/09/26 02:12:12  aidan
  !! Fixed memory leak in new_pnm. Should have deallocated the pointer to
  !! the pnm data rather than nullifying it.
  !!
  !! Revision 1.6  2005/09/26 02:05:09  aidan
  !! Alot of changes in this update:
  !!
  !!     * Included support for rotating pnm images in-place
  !!     * Added routines to assign pnm to pnm
  !!     * Added equivalence operators
  !!     * Made a wrapper (pnm_getchar) for our reading routine
  !!       to make the logic of ignoring comments clearer and
  !!       less buggy
  !!
  !! Revision 1.5  2004/06/07 04:36:08  aidan
  !! Cleaned up the logic of the interface to 'new'. We *must* pass data to these
  !! routines, but previously this had been optional. The image type was mandatory
  !! for the more general routine, but in reality we can assume ppm type with
  !! confidence, as the other routine exists to cater for 2D data (pbm and pgm).
  !!
  !! Revision 1.4  2004/05/25 06:15:11  aidan
  !! Fixed bug in assign_3d_data. The array bounds calculation was incorrect, so
  !! rather than fix it I just chose to make it simpler, as it was not necessary
  !! to check for non-unity lower bounds -- any array passed to a subroutine does
  !! not retain this information, apparently.
  !!
  !! Revision 1.3  2004/05/25 05:32:19  aidan
  !! Fixed bug in read_pnm. Was checking the type of the pnm object before it had
  !! been created! Resulted in ppm data being converted into a quasi-greyscale,
  !! as all the red, green and blue pixels were assigned the value of the red.
  !!
  !! Revision 1.2  2004/05/10 23:53:08  aidan
  !! Made module consistent with the new interface of the binary_io module.
  !!
  !! Added some functions for querying the pnm object (size, maxval). Changed
  !! the reading routine to a subroutine -- this makes it consistent with the
  !! read routine in binary_io and in general is a good idea for a heavily
  !! overloaded name like "read". It is easier for the compiler to differentiate
  !! between different versions of read, which avoids ambiguous referencing.
  !!
  !! Also added some data retrieval functions and overloaded the assignment
  !! operator to be able to simply assign the data from the pnm to an array.
  !!
  !! Revision 1.1  2004/05/06 03:15:16  aidan
  !! Initial revision
  !!

  ! An object to encapsulate our pnm information. The imagetype flag is decoded
  ! thus:
  !  1 = ascii pbm
  !  2 = ascii pgm
  !  3 = ascii ppm
  !  4 = binary pbm
  !  5 = binary pgm
  !  6 = binary ppm
  type pnm_object
     private
     ! integer, dimension(:,:,:), pointer :: data
     integer, dimension(:,:,:), allocatable :: data
     integer :: maxval
     character(len=3) :: imagetype
  end type pnm_object

  ! The maximmum values of a pbm, 8 bit pgm/ppm and a 16 bit pgm/ppm respectively
  integer, dimension(0:2), parameter :: MAXIMUM_VALUES = (/ 1, (2**8)-1, (2**16)-1 /)

  ! Whitespace tokens (ASCII collating sequence)
  character(len=1), parameter :: tab = achar(9), lf = achar(10), cr = achar(13), space = achar(32)
  character(len=*), parameter :: whitespace = tab//lf//cr//space

  ! Function to return the type (pbm,pgm,ppm) given a magic number or pnm file as input
  interface type
     module procedure pnm_type, pnm_type_wrapper
  end interface

  ! Function to return the magic number (1-6) given a type (pbm,pgm,pp) as input
  interface magic_number
     module procedure pnm_magic_number
  end interface

  ! Make a new pnm object
  interface new
     module procedure new_pnm, new_pgm_or_pbm 
  end interface

  ! Read a pnm file, return a pnm object
  interface read
     module procedure pnm_read
  end interface

  ! Write a pnm object to a file
  interface write
     module procedure pnm_write
  end interface

  ! These two routines return an array containing the actual image data.
  ! It is probably best to use straight assignment, as this allows the
  ! compiler to choose to return a 3D or 2D array based on what you are
  ! assigning to.
  interface as_array
     module procedure extract_data_as_array
  end interface
  interface as_array_2d
     module procedure extract_data_as_2darray 
  end interface

  ! Print out some simple summary info on a pnm object
  interface info
     module procedure pnm_info
  end interface

  ! Function to return the size of the data array in a pnm object.
  interface size
     module procedure pnm_size
  end interface

  ! Function to return the maxval of the pnm object.
  interface maxval
     module procedure pnm_maxval
  end interface

  ! This allows convenient assignment from a pnm object to a data
  ! array and vice versa -- to allow easy extraction and insertion
  ! of the data in the pnm object. One drawback is the lack of an 
  ! ability to specify the origin of the data -- so by default this 
  ! is assumed to be the 'bottom left' of the pnm image as displayed.
  interface assignment(=)
     module procedure assign_3d_data, assign_2d_data, assign_3d_pnm, assign_2d_pnm
  end interface

  ! Define equivalence operator for pnm objects
  interface operator(==)
     module procedure pnm_eq_pnm
  end interface

  ! Define non-equivalence operator for pnm objects
  interface operator(/=)
     module procedure pnm_neq_pnm
  end interface

  ! 'Overload' rotate_image from image_transforms 
  interface rotate_image
     module procedure rotate_pnm_wrap, rotate_pnm_real_wrap, &
          rotate_pnm_real_wrap_nocov, rotate_pnm_real_wrap_nomiss
  end interface

  ! 'Overload' rotave_image from image_transforms 
  interface rotave_image
     module procedure rotave_pnm_wrap
  end interface

  !!!!!!!!!!!!!!!!
  ! Public stuff !
  !!!!!!!!!!!!!!!!

  ! data types
  public :: pnm_object

  ! public routines
  public :: new, type, write, info, read, as_array, as_array_2d, size, maxval
  public :: rotate_image, rotave_image

  ! public operators
  public :: assignment(=), operator(==), operator(/=)

contains

  subroutine new_pgm_or_pbm ( pnm, data, maximum, origin )

    ! Wrapper routine for the new_pnm subroutine that allows users to
    ! specify a new pgm or pbm using only 2D data

    ! The new pnm object to be returned
    type (pnm_object), intent(inout) :: pnm

    ! Input parameters
    integer, dimension(:,:), intent(in)    :: data
    ! .. rest are optional
    integer, intent(in), optional          :: maximum
    character(len=2), intent(in), optional :: origin

    ! Local parameters
    integer :: local_max
    logical :: bf, rf
    character(len=2) :: local_origin
    character(len=3) :: local_type

    integer, dimension(3,size(data,1),size(data,2)) :: local_data

    if (present(maximum)) then
       local_max = maximum
    else
       local_max = maxval(data)
    end if

    select case(local_max)
       case(0:MAXIMUM_VALUES(0))
          local_type = 'pbm'
       case(MAXIMUM_VALUES(0)+1:MAXIMUM_VALUES(2))
          local_type = 'pgm'
       case default
          local_max = MAXIMUM_VALUES(2)
          local_type = 'pgm'
    end select

    ! Check if an origin has been specified
    if (present(origin)) then
       local_origin = origin
    else
       ! The default case assumes the origin of the array
       ! is at the bottom left
       local_origin = 'bl'
    end if

    local_data = 0
    local_data(1,:,:) = data
    local_data(2,:,:) = data
    local_data(3,:,:) = data

    ! print *, 'Creating a pbm/pgm'

    call new_pnm( pnm, local_data, local_type, maximum = local_max, origin = local_origin )

  end subroutine new_pgm_or_pbm

  subroutine new_pnm ( pnm, data, imagetype, maximum, origin )

    ! General routine to generate a new pnm object. Note that the data parameter
    ! is 3D, i.e. assumes ppm colour "tuples", even if there is no data in the
    ! second and third fields (i.e. bit and greyscale images).

    ! The new pnm object to be returned
    type (pnm_object), intent(inout) :: pnm

    ! Input parameters
    integer, dimension(:,:,:), intent(in)  :: data
    ! .. rest are optional
    character(len=3), intent(in), optional :: imagetype ! 'pbm', 'pgm' or 'ppm'
    integer, intent(in), optional          :: maximum
    character(len=2), intent(in), optional :: origin

    ! Local parameters
    integer :: width, height, depth, xstart, xfinish, xinc, ystart, yfinish, yinc
    logical :: bf, rf
    character(len=2) :: data_origin

    width  = size(data, 2)
    height = size(data, 3)

    if (present(maximum)) then
       pnm%maxval = maximum
    else
       pnm%maxval = maxval(data)
    end if

    if (present(imagetype)) then
       select case(imagetype)
       case('pbm','pgm','ppm')
          pnm%imagetype = imagetype
       case default
          write(stderr,*)'Invalid image type: ',imagetype
          stop
       end select
    else
       ! Default to ppm is no imagetype specified
       pnm%imagetype = 'ppm'
    end if

    ! if (associated(pnm%data)) deallocate(pnm%data)
    if (allocated(pnm%data)) deallocate(pnm%data)
    
    ! Allocate the memory for the data
    allocate(pnm%data(3,width,height))

    ! Check if an origin has been specified
    if (present(origin)) then
       data_origin = origin
    else
       ! The default case assumes the origin of the array
       ! is at the bottom left
       data_origin = 'bl'
    end if

    select case(data_origin)
    case ('bl') ! Bottom left
       xstart  = lbound(data,2)
       xfinish = ubound(data,2)
       xinc    = +1
       ystart  = ubound(data,3)
       yfinish = lbound(data,3)
       yinc    = -1
    case ('br') ! Bottom right
       xstart  = ubound(data,2)
       xfinish = lbound(data,2)
       xinc    = -1
       ystart  = ubound(data,3)
       yfinish = lbound(data,3)
       yinc    = -1
    case ('tl') ! Top left
       xstart  = lbound(data,2)
       xfinish = ubound(data,2)
       xinc    = +1
       ystart  = lbound(data,3)
       yfinish = ubound(data,3)
       yinc    = +1
    case ('tr') ! Top right
       xstart  = ubound(data,2)
       xfinish = lbound(data,2)
       xinc    = -1
       ystart  = lbound(data,3)
       yfinish = ubound(data,3)
       yinc    = +1
    case default
       stop 'PNM_CLASS :: Illegal origin! Should be one of bl, br, tl or tr'
    end select

    pnm%data(:,:,:) = data(:,xstart:xfinish:xinc,ystart:yfinish:yinc)

    ! Clean the data
    where (pnm%data > pnm%maxval) pnm%data = pnm%maxval

  end subroutine new_pnm

  function pnm_type_wrapper ( pnm )

    type (pnm_object), intent(in) :: pnm

    character(len=3) :: pnm_type_wrapper

    pnm_type_wrapper = pnm%imagetype

  end function pnm_type_wrapper

  function pnm_type ( type_flag )

    integer, intent(in) :: type_flag
    character(len=3)    :: pnm_type

    select case (type_flag)
       case (1,4)
          pnm_type = 'pbm'
       case (2,5)
          pnm_type = 'pgm'
       case (3,6)
          pnm_type = 'ppm'
    end select

  end function pnm_type

  integer function pnm_magic_number ( pnm_type, binary )

    character(len=3), intent(in)  :: pnm_type
    logical, intent(in), optional :: binary

    ! Local variables
    logical :: binary_format

    select case (pnm_type)
       case ('pbm')
          pnm_magic_number = 1
       case ('pgm')
          pnm_magic_number = 2
       case ('ppm')
          pnm_magic_number = 3
    end select

    if (present(binary)) then
       binary_format = binary
    else
       binary_format = .TRUE.
    end if

    if (binary_format) pnm_magic_number = pnm_magic_number + 3

  end function pnm_magic_number

  logical function pnm_data_ok ( pnm ) result(data_ok)
    
    type (pnm_object), intent(in) :: pnm

    data_ok = .TRUE.

    select case(type(pnm))

    case('pbm')
       if ( any(pnm%data /= 0 .and. pnm%data /= 1) ) then
          data_ok = .FALSE.
       end if
    case('pgm')
       if ( any(pnm%data <= 0 .and. pnm%data >= pnm%maxval) ) then
          data_ok = .FALSE.
       end if
    case('ppm')
       if ( any(pnm%data <= 0 .and. pnm%data >= pnm%maxval) ) then
          data_ok = .FALSE.
       end if
    end select
    
    return

  end function pnm_data_ok

  function pnm_getchar(fh) result(c)

    ! We need a 'wrapper' around our read method from the filehandle
    ! so we can filter out comments. Comments begin with a '#' and 
    ! end with a line feed or carriage return

    ! Interface variables
    type (binary_filehandle), intent(inout) :: fh

    ! Return value
    character(len=1) :: c

    ! Local variables
    integer :: error

    call read(fh, c, error)
    if (error /= 0) stop "PNM_CLASS :: Unexpected end of file"

    ! Check if we have a comment character
    if (c == "#") then 
       ! Read characters until we find something that is a line 
       ! feed or carriage return
       do while (c /= lf .AND. c /= cr)
          call read(fh, c, error)
          if (error /= 0) stop "PNM_CLASS :: Unexpected end of file"
       end do
    end if

  end function pnm_getchar

  function zap_spaces (fh) result(firstchar)

    type (binary_filehandle) :: fh

    character(len=1) :: firstchar

    ! Local variables
    integer :: error

    ! Read characters from the specified unit until we find something
    ! that is not whitespace (tab, line feed, carriage return and space)
    do 
       firstchar = pnm_getchar(fh) 
       select case(firstchar)
       case (space,tab,lf,cr)
          cycle
       end select
       exit
    end do

  end function zap_spaces

  integer function pnm_getint(fh) result(number)
    
    ! Based this on the function in libpbm4.c

    ! Input parameters
    type (binary_filehandle), intent(inout) :: fh

    ! Local variables
    integer, parameter :: zero = iachar('0')
    integer :: error
    character(len=1) :: charbuf

    ! Read characters from the specified file handle until we find something
    ! that is not whitespace (tab, line feed, carriage return and space)
    charbuf = zap_spaces(fh)

    ! Make sure what we do have is a numeral
    if (charbuf < '0' .OR. charbuf > '9') then
       write(stderr,*) 'PNM_CLASS :: Error! Junk in image file where integer should be ', iachar(charbuf)
       return
    end if
    
    ! Read numerals from the opened unit until we find something that is not
    ! a numeral -- reconstruct the integer from these numerals
    number = 0
    do 
       ! Note that we don't check if the character which is not a numeral is something
       ! that it shouldn't be, i.e. whitespace. This isn't serious, but is a potential
       ! loophole. Note also that we have 'swallowed' this character that is not a numeral.
       ! This is not a problem for pnm files -- this is mostly probably unwanted white space.
       if (charbuf < '0' .OR. charbuf > '9') exit
       number = number * 10 + (iachar(charbuf) - zero)
       charbuf = pnm_getchar(fh) 
    end do

  end function pnm_getint

  subroutine pnm_read(pnm, filename) !, ierror)

    ! Read a binary pbm/pgm/ppm and return a 'pnm object'

    ! Input parameters
    type (pnm_object), intent(out) :: pnm
    character(len=*), intent(in)   :: filename
    ! integer, intent(out), optional :: ierror

    ! Local variables
    integer :: i, j, k, ncolours, type_flag, maximum_value, bit_depth, error 
    integer :: width, height, value, intbuf, row
    character(len=1) :: charbuf
    integer, dimension(:,:,:), allocatable :: data
    type (binary_filehandle) :: fh

    call open(fh, filename)

    ! We use a 'raw' read for the first two characters in the
    ! pnm file, as we don't want any filtering of comments etc
    ! at this stage
    call read(fh, charbuf, error)
    if (error /=0) stop "PNM_CLASS :: Unexpected end of file"

    call ucase(charbuf)

    ! The first byte of a pnm file should be the letter 'P'
    if (charbuf /= 'P') then
       write(stderr,*) 'PNM_CLASS :: Not a pbm/pgm/ppm file: ',filename
       return
    end if

    ! The next byte is a character number from 1-6 (7 is a PAM file) 
    call read(fh, charbuf, error)
    if (error /=0) stop "PNM_CLASS :: Unexpected end of file"
    read(charbuf,'(I1)') type_flag

    if (type_flag < 1 .or. type_flag > 6) stop 'PNM_CLASS :: Not a pbm/pgm/ppm file'
    if (type_flag >= 1 .and. type_flag <= 3) stop 'PNM_CLASS :: ASCII pbm/pgm/pbm not supported!'

    ! Next come the width and the height
    width  = pnm_getint(fh)
    height = pnm_getint(fh)
    
    ! Allocate the memory for the data
    allocate(data(3,width,height))

    ! pbm files contain no 'maxval'
    if (type(type_flag) == 'pbm') then
       maximum_value = 1
    else
       maximum_value = pnm_getint(fh)
    end if

    bit_depth = 1
    if (maximum_value > 256) bit_depth = 2

    if (maximum_value > MAXIMUM_VALUES(bit_depth)) then
       write(stderr,*) 'PNM_CLASS :: Maximum value too large: ',pnm%maxval
       stop
    end if

    ncolours = 1
    if (type(type_flag) == 'ppm') ncolours = 3

    ! There should be a *single* whitespace character between the header info and 
    ! the data. We don't need to skip it -- it has been 'eaten' by the last call
    ! to pnm_getint. This is a bit bodgy, but should work.
    ! call skip(fh, whitespace)

    if (type(type_flag) == 'pbm') then
       do i = 1, height
          ! A pbm image is 'packed', with each bit representing a pixel. Each 'row' 
          ! is padded to be an integer number of bytes. We make sure we read in at 
          ! least this many bytes
          row = 0
          do j = 1, (width+7)/8
             call read(fh, intbuf, 1)
             ! Cycle through the bits in intbuf, most significant bit first, and
             ! extract out each of the bits in turn
             ! write(*,'(A,I,B32)') 'read number as ', intbuf, intbuf
             do k = 7, 0, -1
                row = row + 1
                if (row > width) exit
                data(1,row,i) = ibits(intbuf,k,1) 
                ! write(*,'(I0)',advance='no') data(1,row,i)
             end do
             ! print *
          end do
       end do
    else
       do i = 1, height
          do j = 1, width
             do k = 1, ncolours
                call read(fh, data(k,j,i), bit_depth)
             end do
          end do
       end do
    end if

    ! Copy the values from 'red' into 'green' and 'blue' for pbm and pgm 
    ! images -- these are not assigned values in the reading routine and
    ! we might as well copy them in there I guess?
    if (type(type_flag) /= 'ppm') then
       data(2,:,:) = data(1,:,:)
       data(3,:,:) = data(1,:,:)
    end if

    call close(fh)

    call new( pnm, data, type(type_flag), maximum_value, origin = "tl" )
    
    return

  end subroutine pnm_read

  subroutine pnm_write(pnm, filename, ascii)

    ! Input variables
    type (pnm_object), intent(in) :: pnm
    character(len=*), intent(in)  :: filename
    logical, intent(in), optional :: ascii

    ! Local variables
    logical :: write_ascii
    integer :: ncolours, bit_depth, i, j, k, intbuf, bitpos
    character(len=100) :: buffer

    type (binary_filehandle) :: fh

    if (present(ascii)) then
       ! write_ascii = ascii
       write(stderr,*) 'PNM_CLASS :: Despite appearances, ASCII write of pnm files is currently not supported'
       return
    else
       write_ascii = .FALSE.
    end if

    ncolours = 1
    if (pnm%imagetype == 'ppm') ncolours = 3

    bit_depth = 1
    if (pnm%maxval > 256) bit_depth = 2

    call open(fh, filename, action='write', buflen=8192)

    ! Write the header. We write each of our header values into a 
    ! string first and then pass this to our buffered binary io 
    ! routine with trailing spaces removed and a trailing line feed 
    ! added. Line feed is the way unix files designate a new line. 
    ! In theory any whitespace would have been sufficient but both 
    ! ImageMagick and ImageJ seemed to only accept a new line
    ! separator between the fields in the header. Either they didn't
    ! read the standard too carefully or they didn't care .. oh well.
    ! It might have been slightly more efficient to write the header
    ! into a single string and then write that to the file, but in
    ! reality the saving wouldn't be much and it is way way easier
    ! to get rid of all the clunky spaces and such doing it this way

    ! Magic Number
    write(buffer,'("P",I0)') magic_number(pnm%imagetype, .NOT. write_ascii)
    call write(fh,trim(buffer)//lf)

    ! Width
    write(buffer,'(I0)') size(pnm%data,2)
    call write(fh,trim(buffer)//lf)

    ! Height
    write(buffer,'(I0)') size(pnm%data,3)
    call write(fh,trim(buffer)//lf)

    if (pnm%imagetype /= 'pbm') then
       ! Maxval
       write(buffer,'(I0)') pnm%maxval
       call write(fh,trim(buffer)//lf)
    end if

    ! Raster data
    if (type(pnm) == 'pbm') then
       ! A pbm image is 'packed', with each bit representing a pixel. 
       ! Each 'row' is padded to be an integer number of bytes.
       do i = 1, size(pnm%data,3)
          intbuf = 0
          bitpos = 7
          do j = 1, size(pnm%data,2)
             if (pnm%data(1,j,i) == 1) then
                ! Set the bit in our buffer
                intbuf = ibset(intbuf,bitpos) 
             end if
             bitpos = bitpos - 1
             if (bitpos == -1) then
                ! We have a byte's worth of data, write it out
                call write(fh,intbuf,1,littleendian=.FALSE.)
                ! Reset our buffer and bit position index
                intbuf = 0
                bitpos = 7
             end if
          end do
          ! Write a partially filled byte
          if (bitpos /= 7) then
             call write(fh,intbuf,1,littleendian=.FALSE.)
          end if
       end do
    else
       do i = 1, size(pnm%data,3)
          do j = 1, size(pnm%data,2)
             do k = 1, ncolours
                call write(fh,pnm%data(k,j,i),bit_depth,littleendian=.FALSE.)
             end do
          end do
       end do
    end if

    call close(fh)

  end subroutine pnm_write

  subroutine pnm_info (pnm)

    ! Writes some information about the pnm file to stdout

    ! Interface variables
    type (pnm_object), intent(in) :: pnm

    write(stdout,*) "            Type: ", type(pnm)
    write(stdout,*) "      Size (wxh): ", size(pnm%data,2), ' x ', size(pnm%data,3)
    write(stdout,*) "   Maximum value: ", pnm%maxval
    write(stdout,*) "Actual Max value: ", maxval(pnm%data(1,:,:))
    write(stdout,*) "Actual Min value: ", minval(pnm%data(1,:,:))

  end subroutine pnm_info

  integer function pnm_size (pnm, dim)

    ! Equivalent of the intrinsic size command for pnm objects

    ! Interface variables
    type (pnm_object), intent(in) :: pnm
    integer, intent(in), optional :: dim

    ! Local variables
    integer :: width, height

    width  = size(pnm%data,2)
    height = size(pnm%data,3)

    if (present(dim)) then
       select case (dim)
          case (1)
             pnm_size = width
             return
          case (2)
             pnm_size = height
             return
          case default
             write(stderr,'(A,I0,A)') 'PNM_CLASS :: Illegal dimension: ',dim,'. Must be 1 or 2'
             stop
          end select
    end if
    pnm_size = width * height

  end function pnm_size

  integer function pnm_maxval (pnm)

    ! Sort of equivalent of the maxval intrinsic, but actually returns
    ! the value in the header rather than the actual maximum value

    type (pnm_object), intent(in) :: pnm
    
    pnm_maxval = pnm%maxval

  end function pnm_maxval

  function extract_data_as_array (pnm, origin) result(data)

    ! Extracts the raster data from the pnm object and returns it as an 
    ! integer array

    type (pnm_object), intent(in)           :: pnm
    character(len=2), intent(in), optional  :: origin

    ! Return type
    integer, dimension(3,size(pnm%data,2),size(pnm%data,3)) :: data

    ! Local variables
    character(len=2) :: data_origin
    integer :: xstart, xfinish, xinc, ystart, yfinish, yinc

    ! Check if an origin has been specified
    if (present(origin)) then
       data_origin = origin
    else
       ! The default case assumes the origin of the array
       ! is at the bottom left
       data_origin = 'bl'
    end if

    select case(data_origin)
    case ('bl') ! Bottom left
       xstart  = lbound(data,2)
       xfinish = ubound(data,2)
       xinc    = +1
       ystart  = ubound(data,3)
       yfinish = lbound(data,3)
       yinc    = -1
    case ('br') ! Bottom right
       xstart  = ubound(data,2)
       xfinish = lbound(data,2)
       xinc    = -1
       ystart  = ubound(data,3)
       yfinish = lbound(data,3)
       yinc    = -1
    case ('tl') ! Top left
       xstart  = lbound(data,2)
       xfinish = ubound(data,2)
       xinc    = +1
       ystart  = lbound(data,3)
       yfinish = ubound(data,3)
       yinc    = +1
    case ('tr') ! Top right
       xstart  = ubound(data,2)
       xfinish = lbound(data,2)
       xinc    = -1
       ystart  = lbound(data,3)
       yfinish = ubound(data,3)
       yinc    = +1
    case default
       stop 'PNM_CLASS :: Illegal origin! Should be one of bl, br, tl or tr'
    end select

    data(:,:,:) = pnm%data(:,xstart:xfinish:xinc,ystart:yfinish:yinc)

  end function extract_data_as_array

  function extract_data_as_2darray (pnm, origin) result(data)

    ! Extracts the raster data from the pnm object and returns it as an 
    ! two-dimensional integer array

    type (pnm_object), intent(in)           :: pnm
    character(len=2), intent(in), optional  :: origin

    ! Return type
    integer, dimension(size(pnm%data,2),size(pnm%data,3)) :: data

    ! Local variables
    integer, dimension(3,size(pnm%data,2),size(pnm%data,3)) :: local_data
    character(len=2) :: data_origin

    ! Check if an origin has been specified
    if (present(origin)) then
       data_origin = origin
    else
       ! The default case assumes the origin of the array
       ! is at the bottom left
       data_origin = 'bl'
    end if

    local_data = extract_data_as_array(pnm, data_origin)
    
    if (type(pnm) == 'ppm') then 
       ! Use the quantisation formula given in the ppmtopgm docs, i.e.
       ! .299 red + .587 green + .114 blue
       data = 0.299 * local_data(1,:,:) + 0.587 * local_data(2,:,:) + 0.114 * local_data(3,:,:)
    else
       data = local_data(1,:,:)
    end if

  end function extract_data_as_2darray

  subroutine assign_3d_data (data, pnm)
    
    ! Convenience routine to allow simple assignment to a 3D array
    ! from a pnm object

    ! Interface variables
    integer, dimension(:,:,:), intent(out) :: data
    type (pnm_object), intent(in)          :: pnm

    ! Local variables
    integer :: data_width, data_height, pnm_width, pnm_height, final_width, final_height

    data_width  = size(data,2)
    data_height = size(data,3)

    pnm_width  = size(pnm,1)
    pnm_height = size(pnm,2)

    final_width  = min(data_width,pnm_width)
    final_height = min(data_height,pnm_height)

    if ( (data_width < pnm_width) .OR. (data_height < pnm_height) ) then
       write(stderr,'(A,I0,A,I0,A,I0,A,I0)') 'PNM_CLASS :: Data clipped on assignment: from ',&
            pnm_width,' x ',pnm_height,' to ',final_width,' x ',final_height
    end if

    ! Don't need to account for non-unity lower bounds of data array as these
    ! are reset to one inside a subroutine

    ! Assign the data -- use call to extract function to get the data
    data(:,1:final_width,1:final_height) = extract_data_as_array(pnm)

  end subroutine assign_3d_data

  subroutine assign_2d_data (data, pnm)

    ! Convenience routine to allow simple assignment to a 2D array
    ! from a pnm object

    ! Interface variables
    integer, dimension(:,:), intent(out) :: data
    type (pnm_object), intent(in)        :: pnm

    ! Local variables
    integer :: data_width, data_height, pnm_width, pnm_height, final_width, final_height

    data_width  = size(data,1)
    data_height = size(data,2)

    pnm_width  = size(pnm,1)
    pnm_height = size(pnm,2)

    final_width  = min(data_width,pnm_width)
    final_height = min(data_height,pnm_height)

    if ( (data_width < pnm_width) .OR. (data_height < pnm_height) ) then
       write(stderr,'(A,I0,A,I0,A,I0,A,I0)') 'PNM_CLASS :: Data clipped on assignment: from ',&
            pnm_width,' x ',pnm_height,' to ',final_width,' x ',final_height
    end if

    ! Don't need to account for non-unity lower bounds of data array as these
    ! are reset to one inside a subroutine

    ! Assign the data -- use call to extract function to get the data
    data(1:final_width,1:final_height) = extract_data_as_2darray(pnm)

  end subroutine assign_2d_data

  subroutine assign_3d_pnm (pnm, data)

    ! Convenience routine to allow simple assignment to a pnm object
    ! from a 2D array.

    ! Interface variables
    type (pnm_object), intent(out)        :: pnm
    integer, dimension(:,:,:), intent(in) :: data

    call new(pnm, data)

  end subroutine assign_3d_pnm

  subroutine assign_2d_pnm (pnm, data)

    ! Convenience routine to allow simple assignment to a pnm object
    ! from a 2D array.

    ! Interface variables
    type (pnm_object), intent(out)      :: pnm
    integer, dimension(:,:), intent(in) :: data

    call new(pnm, data)

  end subroutine assign_2d_pnm

  logical function pnm_eq_pnm (pnm_left, pnm_right)

    ! Equivalence function
    
    ! Interface variables
    type (pnm_object), intent(in) :: pnm_left, pnm_right

    pnm_eq_pnm = .FALSE.

    if (pnm_left%maxval /= pnm_right%maxval) return
    if (pnm_left%imagetype /= pnm_right%imagetype) return
    if (any(pnm_left%data /= pnm_right%data)) return
    
    pnm_eq_pnm = .TRUE.

  end function pnm_eq_pnm

  logical function pnm_neq_pnm (pnm_left, pnm_right)

    ! Non-equivalence function
    
    ! Interface variables
    type (pnm_object), intent(in) :: pnm_left, pnm_right

    pnm_neq_pnm = .NOT. (pnm_left == pnm_right)

  end function pnm_neq_pnm

  subroutine rotate_pnm_wrap(data, angle, missingvalue, coverage)

    ! Rotate and/or resize a pnm object with single precision angle and missing value
    ! Interface variables
    type (pnm_object), intent(inout)         :: data
    real(kind=rd_kind), intent(in)           :: angle
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    integer, dimension(size(data%data,2),size(data%data,3)) :: local
    integer, dimension(3,size(data%data,2),size(data%data,3)) :: local_ppm
    integer :: i

    if (type(data) == 'ppm') then
       ! If we have ppm data we'll have to use different input and
       ! output temporary arrays, and loop over the RGB channels in
       ! the ppm file
       local_ppm = data
       if (present(missingvalue)) then
          if (present(coverage)) then
             do i = 1, 3
                call rotate_image(local_ppm(i,:,:),angle,missingvalue,coverage)
             end do
          else
             do i = 1, 3
                call rotate_image(local_ppm(i,:,:),angle,missingvalue)
             end do
          end if
       else
          do i = 1, 3
             call rotate_image(local_ppm(i,:,:),angle)
          end do
       end if
       data = local_ppm
    else
       local = data
       if (present(missingvalue)) then
          if (present(coverage)) then
             call rotate_image(local,angle,missingvalue,coverage)
          else
             call rotate_image(local,angle,missingvalue)
          end if
       else
          call rotate_image(local,angle)
       end if
       data = local
    end if

  end subroutine rotate_pnm_wrap

  subroutine rotate_pnm_real_wrap(data, angle, missingvalue, coverage)

    ! Rotate and/or resize a pnm object with single precision angle, 
    ! missing value and coverage specified

    ! Interface variables
    type (pnm_object), intent(inout) :: data
    real(kind=rs_kind), intent(in)   :: angle
    integer, intent(in)              :: missingvalue
    real(kind=rs_kind), intent(in)   :: coverage

    call rotate_image(data,real(angle,rd_kind),missingvalue,coverage)

  end subroutine rotate_pnm_real_wrap

  subroutine rotate_pnm_real_wrap_nocov(data, angle, missingvalue)

    ! Rotate and/or resize a pnm object with single precision angle, 
    ! missing value specified

    ! Interface variables
    type (pnm_object), intent(inout) :: data
    real(kind=rs_kind), intent(in)   :: angle
    integer, intent(in)              :: missingvalue

    call rotate_image(data,real(angle,rd_kind),missingvalue)

  end subroutine rotate_pnm_real_wrap_nocov

  subroutine rotate_pnm_real_wrap_nomiss(data, angle)

    ! Rotate and/or resize a pnm object with single precision angle

    ! Interface variables
    type (pnm_object), intent(inout) :: data
    real(kind=rs_kind), intent(in)   :: angle

    call rotate_image(data,real(angle,rd_kind))

  end subroutine rotate_pnm_real_wrap_nomiss

  subroutine rotave_pnm_wrap(data, nfold, missingvalue, coverage)

    ! Rotate and/or resize a pnm object with single precision angle and missing value
    ! Interface variables
    type (pnm_object), intent(inout)         :: data
    integer, intent(in)                      :: nfold
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    integer, dimension(size(data%data,2),size(data%data,3)) :: local
    integer, dimension(3,size(data%data,2),size(data%data,3)) :: local_ppm
    integer :: i

    if (type(data) == 'ppm') then
       ! If we have ppm data we'll have to use different input and
       ! output temporary arrays, and loop over the RGB channels in
       ! the ppm file
       local_ppm = data
       if (present(missingvalue)) then
          if (present(coverage)) then
             do i = 1, 3
                call rotave_image(local_ppm(i,:,:),nfold,missingvalue,coverage)
             end do
          else
             do i = 1, 3
                call rotave_image(local_ppm(i,:,:),nfold,missingvalue)
             end do
          end if
       else
          do i = 1, 3
             call rotave_image(local_ppm(i,:,:),nfold)
          end do
       end if
       data = local_ppm
    else
       local = data
       if (present(missingvalue)) then
          if (present(coverage)) then
             call rotave_image(local,nfold,missingvalue,coverage)
          else
             call rotave_image(local,nfold,missingvalue)
          end if
       else
          call rotave_image(local,nfold)
       end if
       data = local
    end if

  end subroutine rotave_pnm_wrap

end module pnm_class
