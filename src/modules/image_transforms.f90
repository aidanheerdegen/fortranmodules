module image_transforms

  use polysample, only: polySamp
  use precision
  use rotmatrix_class
  use fundamental_constants, only: radian, pi
  use variable_array, only: push, splice
  use statistics

  implicit none

  private

  character(len=*), parameter :: version = "$Id: image_transforms.f90,v 1.4 2005/07/01 03:22:50 aidan Exp aidan $"

  !! Routines for rotating and/or scaling images. Currently supports
  !! rotating double, real and integer images. There are subroutine and
  !! functional interfaces. The double precision versions accept only
  !! double precision arguments. The real wrappers will accept double
  !! and single precision angles, but the type of the missing value must
  !! match that of the data that is being manipulated. The integer wrapper 
  !! script will accept either, but if two such arguments are supplied they 
  !! must be the same kind. 
  !! 
  !! The code dealing with direct support for pnm objects has been commented
  !! out (it is probably not in working order either). Rotating pnm's is
  !! problematic -- what if we have a ppm? Really we should rotate all three
  !! fields simultaneously to save computing the complicated overlaps for each
  !! channel. We do not do this currently. More work is required. Possibly
  !! this rotation operation should be moved into the pnm module -- then it
  !! could operate directly on the internal integer data and we wouldn't have
  !! to have a wrapper for it in this module. Hmmmmmm ....

  !! $Log: image_transforms.f90,v $
  !! Revision 1.4  2005/07/01 03:22:50  aidan
  !! Added rotational averaging routines, both subroutine and functional interfaces.
  !!
  !! Revision 1.3  2005/06/29 03:25:03  aidan
  !! Forgot to delete the set_coverage, get_coverage routines from the interface.
  !! Fixed.
  !!
  !! Revision 1.2  2005/06/29 03:20:39  aidan
  !! Added support for specifying minimum coverage.
  !!
  !! Revision 1.1  2004/11/17 02:34:33  aidan
  !! Initial revision
  !!
  !!

  ! Subroutine interface
  interface rotate_image
     module procedure rotate_image_real8_sub, rotate_inplace, &
          rotate_inplace_integer_wrap, rotate_inplace_integerd_wrap, &
          rotate_inplace_real_wrap, rotate_inplace_reald_wrap, &
          rotate_image_integer_wrap, rotate_image_integerd_wrap, &
          rotate_image_real_wrap, rotate_image_reald_wrap
  end interface

  interface rotave_image
     module procedure rotave_real8_sub, rotave_resize_real8_sub
     module procedure rotave_real8_sub_nomiss, rotave_resize_real8_sub_nomiss
     module procedure rotave_integer_sub, rotave_real_sub
     module procedure rotave_real_inplace_sub, rotave_integer_inplace_sub
  end interface

  ! Functional interface
  interface rotimage
     module procedure rotate_resize_image_real8_fn, rotate_image_real8_fn, &
          rotate_resize_image_real_fn, rotate_image_real_fn, &
          rotate_resize_image_reald_fn, rotate_image_reald_fn, &
          rotate_resize_image_integer_fn, rotate_image_integer_fn, &
          rotate_resize_image_integerd_fn, rotate_image_integerd_fn
  end interface

  interface rotave
     module procedure rotave_real8_fn, rotave_real_fn, rotave_integer_fn
  end interface

!!$  interface centre
!!$     module procedure centre_real8
!!$  end interface
!!$
!!$  interface align
!!$     module procedure translational_align_real8
!!$  end interface

  interface translate_image
     module procedure translate_image_real8_sub, translate_image_realdelta
     module procedure translate_image_real_wrap, translate_image_integer_wrap
     module procedure translate_inplace, translate_inplace_realdelta
     module procedure translate_inplace_real_wrap, translate_inplace_integer_wrap
  end interface

  interface translate
     module procedure translate_function, translate_function_realdelta
     module procedure translate_function_real_wrap, translate_function_integer_wrap
  end interface

  public :: rotate_image, rotimage, rotave_image, rotave, translate_image, translate

  logical, parameter :: debug =  .FALSE. !.TRUE.

contains

  ! This is the main workhorse routine for this module. Most of the other
  ! routines are wrappers for this routine.

  subroutine rotate_image_real8_sub(data, angle, output, missingvalue, coverage)

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    real(kind=rd_kind), intent(in)                  :: angle
    real(kind=rd_kind), dimension(:,:), intent(out) :: output
    real(kind=rd_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(3,(size(output,1)+1)*(size(output,2)+1)) :: vertices
    integer :: i, nx, ny, newnx, newny

    type (rotmatrix) :: R

    if (debug) print *,'in rotate_image'

    nx = size(data,1)
    ny = size(data,2)
    newnx = size(output,1)
    newny = size(output,2)

    if (debug) print *,nx,ny,newnx,newny

    ! Reproduce the current vertices (runs from 0 -> nx, 0 -> ny)
    vertices(1,:) = pack(spread(real((/ (i, i = 0, newnx) /),rd_kind) * &
         (real(nx,rd_kind)/real(newnx,rd_kind)), 2, newny+1),.TRUE.)
    vertices(2,:) = pack(spread(real((/ (i, i = 0, newny) /),rd_kind) * &
         (real(ny,rd_kind)/real(newny,rd_kind)), 1, newnx+1),.TRUE.)
    vertices(3,:) = 0.d0

    if (debug) print *,'Original vertices'
    if (debug) print '(2F0.4)',vertices(1:2,1:10)

    ! stop

    vertices(1,:) = vertices(1,:) - real(nx,rd_kind)/2.d0
    vertices(2,:) = vertices(2,:) - real(ny,rd_kind)/2.d0
    
    if (debug) print *,'Change origin dx dy = ',-nx/2,-ny/2
    if (debug) print '(2F0.4)',vertices(1:2,:)

    ! Make a rotation matrix -- rotate by angle about 'z' -- the "third" 
    ! dimension, i.e. out of the plane of the image
    R = as_rotmatrix(3, angle) 

    if (debug) print *,'rotation_matrix:',as_matrix(R)
    call rotate(R,vertices)

    if (debug) print *,'Post rotation'
    if (debug) print '(2F0.4)',vertices(1:2,:)

    vertices(1,:) = vertices(1,:) + real(nx,rd_kind)/2.d0
    vertices(2,:) = vertices(2,:) + real(ny,rd_kind)/2.d0
    
    if (debug) print *,'Origin shifted back'
    if (debug) print '(2F0.4)',vertices(1:2,:)

    if (present(missingvalue)) then
       if (present(coverage)) then
          call polySamp(data, vertices(1,:), vertices(2,:), output, missingvalue, coverage=coverage)
       else
          call polySamp(data, vertices(1,:), vertices(2,:), output, missingvalue)
       end if
    else
       call polySamp(data, vertices(1,:), vertices(2,:), output)
    end if

    ! call reorthogonalise(R)
    ! if (debug) print *,'rotation_matrix:',as_matrix(R)

  end subroutine rotate_image_real8_sub


  ! Subroutine interface


  ! In place rotation

  subroutine rotate_inplace(data, angle, missingvalue, coverage)

    ! Rotate and/or resize a single precision matrix, angle and missing value

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(inout) :: data
    real(kind=rd_kind), intent(in)                    :: angle
    real(kind=rd_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,data,missingvalue,coverage)
       else
          call rotate_image(data,angle,data,missingvalue)
       end if
    else
       call rotate_image(data,angle,data)
    end if

  end subroutine rotate_inplace

  subroutine rotate_inplace_integer_wrap(data, angle, missingvalue, coverage)

    ! Rotate and/or resize an integer matrix, single precision angle

    ! Interface variables
    integer, dimension(:,:), intent(inout)   :: data
    real(kind=rs_kind), intent(in)           :: angle
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output)
    end if

    data = nint(local_output)

  end subroutine rotate_inplace_integer_wrap

  subroutine rotate_inplace_integerd_wrap(data, angle, missingvalue, coverage)

    ! Rotate and/or resize an integer matrix, double precision angle

    ! Interface variables
    integer, dimension(:,:), intent(inout)   :: data
    real(kind=rd_kind), intent(in)           :: angle
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),angle,local_output)
    end if

    data = nint(local_output)

  end subroutine rotate_inplace_integerd_wrap

  subroutine rotate_inplace_real_wrap(data, angle, missingvalue, coverage)

    ! Rotate and/or resize a single precision matrix, angle and missing value

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(inout) :: data
    real(kind=rs_kind), intent(in)                    :: angle
    real(kind=rs_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output)
    end if

    data = real(local_output,rs_kind)

  end subroutine rotate_inplace_real_wrap

  subroutine rotate_inplace_reald_wrap(data, angle, missingvalue, coverage)

    ! Rotate and/or resize a single precision matrix, and missing value, double
    ! precision angle

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(inout) :: data
    real(kind=rd_kind), intent(in)                    :: angle
    real(kind=rs_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),angle,local_output)
    end if

    data = real(local_output,rs_kind)

  end subroutine rotate_inplace_reald_wrap


  ! Output to another variable

  subroutine rotate_image_integer_wrap(data, angle, output, missingvalue, coverage)

    ! Rotate and/or resize an integer matrix, single precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    real(kind=rs_kind), intent(in)           :: angle
    integer, dimension(:,:), intent(out)     :: output
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output)
    end if

    output = nint(local_output)

  end subroutine rotate_image_integer_wrap

  subroutine rotate_image_integerd_wrap(data, angle, output, missingvalue, coverage)

    ! Rotate and/or resize an integer matrix, double precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    real(kind=rd_kind), intent(in)           :: angle
    integer, dimension(:,:), intent(out)     :: output
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind))
       else
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),angle,local_output)
    end if

    output = nint(local_output)

  end subroutine rotate_image_integerd_wrap

  subroutine rotate_image_real_wrap(data, angle, output, missingvalue, coverage)

    ! Rotate and/or resize a single precision matrix, angle and missing value

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in)  :: data
    real(kind=rs_kind), intent(in)                  :: angle
    real(kind=rs_kind), dimension(:,:), intent(out) :: output
    real(kind=rs_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),real(angle,rd_kind),local_output)
    end if

    output = real(local_output,rs_kind)

  end subroutine rotate_image_real_wrap

  subroutine rotate_image_reald_wrap(data, angle, output, missingvalue, coverage)

    ! Rotate and/or resize a single precision matrix, angle and missing value

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in)  :: data
    real(kind=rd_kind), intent(in)                  :: angle
    real(kind=rs_kind), dimension(:,:), intent(out) :: output
    real(kind=rs_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind),coverage)
       else
          call rotate_image(real(data,rd_kind),angle,local_output,real(missingvalue,rd_kind))
       end if
    else
       call rotate_image(real(data,rd_kind),angle,local_output)
    end if

    output = real(local_output,rs_kind)

  end subroutine rotate_image_reald_wrap
  

!!$  subroutine rotate_image_pnm_wrap(data, angle, output, missingvalue)
!!$
!!$    ! Rotate and/or resize a pnm object with single precision angle and missing value
!!$
!!$    ! Interface variables
!!$    type (pnm_object), intent(in)            :: data
!!$    real(kind=rs_kind), intent(in)           :: angle
!!$    type (pnm_object), intent(out)           :: output
!!$    real(kind=rs_kind), intent(in), optional :: missingvalue
!!$
!!$    ! Local variables
!!$    integer, dimension(size(pnm,1),size(pnm,2))                   :: local_input
!!$    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: local_output
!!$
!!$    local_input = data
!!$
!!$    if (present(missingvalue)) then
!!$       call rotate_image(real(local_input,rd_kind),real(angle,rd_kind),local_output,real(missingvalue,rd_kind))
!!$    else
!!$       call rotate_image(real(local_input,rd_kind),real(angle,rd_kind),local_output)
!!$    end if
!!$
!!$    output = nint(local_output)
!!$
!!$  end subroutine rotate_image_pnm_wrap
!!$
!!$  subroutine rotate_image_pnmd_wrap(data, angle, output, missingvalue)
!!$
!!$    ! Rotate and/or resize a pnm object with double precision angle and missing value
!!$
!!$    ! Interface variables
!!$    type (pnm_object), intent(in)            :: data
!!$    real(kind=rs_kind), intent(in)           :: angle
!!$    type (pnm_object), intent(out)           :: output
!!$    real(kind=rs_kind), intent(in), optional :: missingvalue
!!$
!!$    ! Local variables
!!$    integer, dimension(size(pnm,1),size(pnm,2))                   :: local_input
!!$    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: local_output
!!$
!!$    local_input = pnm
!!$
!!$    if (present(missingvalue)) then
!!$       call rotate_image(real(local_input,rd_kind),angle,local_output,missingvalue)
!!$    else
!!$       call rotate_image(real(local_input,rd_kind),angle,local_output)
!!$    end if
!!$
!!$    output = nint(local_output)
!!$
!!$  end subroutine rotate_image_pnmd_wrap


  ! Functional interface

  function rotate_resize_image_real8_fn(data, angle, newsize, missingvalue, coverage) result(rotimage)

    ! Rotate and resize a double precision image

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in) :: data
    real(kind=rd_kind), intent(in)                 :: angle
    integer, dimension(2), intent(in)              :: newsize
    real(kind=rd_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rd_kind), dimension(newsize(1),newsize(2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_resize_image_real8_fn

  function rotate_image_real8_fn(data, angle, missingvalue, coverage) result(rotimage)

    ! Rotate a double precision image

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in) :: data
    real(kind=rd_kind), intent(in)                 :: angle
    real(kind=rd_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_image_real8_fn

  function rotate_resize_image_real_fn(data, angle, newsize, missingvalue, coverage) result(rotimage)

    ! Rotate and resize a real image

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in) :: data
    real(kind=rs_kind), intent(in)                 :: angle
    integer, dimension(2), intent(in)              :: newsize
    real(kind=rs_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rs_kind), dimension(newsize(1),newsize(2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_resize_image_real_fn

  function rotate_image_real_fn(data, angle, missingvalue, coverage) result(rotimage)

    ! Rotate a real image

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in) :: data
    real(kind=rs_kind), intent(in)                 :: angle
    real(kind=rs_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rs_kind), dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_image_real_fn

  function rotate_resize_image_reald_fn(data, angle, newsize, missingvalue, coverage) result(rotimage)

    ! Rotate and resize a real image, double precision angle

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in) :: data
    real(kind=rd_kind), intent(in)                 :: angle
    integer, dimension(2), intent(in)              :: newsize
    real(kind=rs_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rs_kind), dimension(newsize(1),newsize(2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_resize_image_reald_fn

  function rotate_image_reald_fn(data, angle, missingvalue, coverage) result(rotimage)

    ! Rotate a real image, double precision angle

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in) :: data
    real(kind=rd_kind), intent(in)                 :: angle
    real(kind=rs_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rs_kind), dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_image_reald_fn

  function rotate_resize_image_integer_fn(data, angle, newsize, missingvalue, coverage) result(rotimage)

    ! Rotate and resize an image, single precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    real(kind=rs_kind), intent(in)           :: angle
    integer, dimension(2), intent(in)        :: newsize
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Function result
    integer, dimension(newsize(1),newsize(2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_resize_image_integer_fn

  function rotate_image_integer_fn(data, angle, missingvalue, coverage) result(rotimage)

    ! Rotate an image, single precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    real(kind=rs_kind), intent(in)           :: angle
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Function result
    integer, dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_image_integer_fn

  function rotate_resize_image_integerd_fn(data, angle, newsize, missingvalue, coverage) result(rotimage)

    ! Rotate and resize an image, double precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    real(kind=rd_kind), intent(in)           :: angle
    integer, dimension(2), intent(in)        :: newsize
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Function result
    integer, dimension(newsize(1),newsize(2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_resize_image_integerd_fn

  function rotate_image_integerd_fn(data, angle, missingvalue, coverage) result(rotimage)

    ! Rotate an image, double precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    real(kind=rd_kind), intent(in)           :: angle
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Function result
    integer, dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotate_image(data,angle,rotimage,missingvalue,coverage)
       else
          call rotate_image(data,angle,rotimage,missingvalue)
       end if
    else
       call rotate_image(data,angle,rotimage)
    end if

  end function rotate_image_integerd_fn

!!$  function rotate_resize_image_pnm_fn(data, angle, newsize, missingvalue) result(rotimage)
!!$
!!$    ! Rotate and resize a pnm image
!!$
!!$    ! Interface variables
!!$    type (pnm_object), intent(in)     :: data
!!$    real(kind=rs_kind), intent(in)               :: angle
!!$    integer, dimension(2), intent(in) :: newsize
!!$    integer, intent(in), optional     :: missingvalue
!!$
!!$    ! Function result
!!$    type (pnm_object) :: rotimage
!!$
!!$    ! Local variable
!!$    integer, dimension(newsize(1),newsize(2)) :: local_output
!!$
!!$    if (present(missingvalue)) then
!!$       call rotate_image(as_array_2d(data),angle,local_output,missingvalue)
!!$    else
!!$       call rotate_image(as_array_2d(data),angle,local_output)
!!$    end if
!!$
!!$    rotimage = local_output
!!$
!!$  end function rotate_resize_image_pnm_fn
!!$
!!$  function rotate_image_pnm_fn(data, angle, missingvalue) result(rotimage)
!!$
!!$    ! Rotate an image
!!$
!!$    ! Interface variables
!!$    type(pnm_object), intent(in)   :: data
!!$    real(kind=rs_kind), intent(in) :: angle
!!$    integer, intent(in), optional  :: missingvalue
!!$
!!$    ! Function result
!!$    type(pnm_object) :: rotimage
!!$
!!$    if (present(missingvalue)) then
!!$       call rotate_image(data,angle,rotimage,missingvalue)
!!$    else
!!$       call rotate_image(data,angle,rotimage)
!!$    end if
!!$
!!$  end function rotate_image_pnm_fn
  
  subroutine rotave_real8_sub(data, nfold, missingvalue, coverage)

    ! Rotationally average in-place a double precision matrix with 
    ! missing value defined.

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(inout) :: data
    integer, intent(in)                               :: nfold
    real(kind=rd_kind), intent(in)                    :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    ! Local variables
    real(rd_kind), dimension(size(data,1),size(data,2)) :: buffer, totbuffer
    integer, dimension(size(data,1),size(data,2))       :: bufcount, totcount
    real(kind=rd_kind) :: step
    integer :: i, j, k

    ! Initialise the 'total' buffer with the original data
    totbuffer = data

    ! Ditto for the 'total' counts buffer (ignoring missingvalues)
    totcount = 0
    where(totbuffer /= missingvalue) totcount = 1
    
    ! Precompute the rotation corresponding to an nfold rotation
    step = (2.d0/real(nfold,rd_kind))*pi
    
    ! Now step through 
    do i = 1, nfold-1
       buffer = data
       if (present(coverage)) then
          call rotate_image(buffer,real(i,rd_kind)*step,missingvalue,coverage)
       else
          call rotate_image(buffer,real(i,rd_kind)*step,missingvalue)
       end if
       do j = 1, size(data,1)
          do k = 1, size(data,2)
             if (buffer(j,k) /= missingvalue) then
                if (totbuffer(j,k) /= missingvalue) then
                   totbuffer(j,k) = totbuffer(j,k) + buffer(j,k)
                else
                   totbuffer(j,k) = buffer(j,k)
                end if
                totcount(j,k) = totcount(j,k) + 1
             end if
          end do
       end do
    end do

    where(totcount > 0) totbuffer = totbuffer/real(totcount)

    data = totbuffer

  end subroutine rotave_real8_sub
  
  subroutine rotave_resize_real8_sub(data, nfold, output, missingvalue, coverage)

    ! Rotationally average a double precision matrix with missing 
    ! value defined. If output is a different size to data, this 
    ! will be done with re-sizing.

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    integer, intent(in)                             :: nfold
    real(kind=rd_kind), dimension(:,:), intent(out) :: output
    real(kind=rd_kind), intent(in)                  :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(rd_kind), dimension(size(data,1),size(data,2))     :: inbuffer
    real(rd_kind), dimension(size(output,1),size(output,2)) :: outbuffer, totbuffer
    integer, dimension(size(output,1),size(output,2))       :: bufcount, totcount
    real(kind=rd_kind) :: step
    integer :: i, j, k

    totcount = 0
    totbuffer = missingvalue
    
    step = (2.d0/real(nfold,rd_kind))*pi
    
    do i = 0, nfold-1
       inbuffer = data
       outbuffer = missingvalue
       if (present(coverage)) then
          call rotate_image(inbuffer,real(i,rd_kind)*step,outbuffer,missingvalue,coverage)
       else
          call rotate_image(inbuffer,real(i,rd_kind)*step,outbuffer,missingvalue)
       end if
       do j = 1, size(output,1)
          do k = 1, size(output,2)
             if (outbuffer(j,k) /= missingvalue) then
                if (totbuffer(j,k) /= missingvalue) then
                   totbuffer(j,k) = totbuffer(j,k) + outbuffer(j,k)
                else
                   totbuffer(j,k) = outbuffer(j,k)
                end if
                totcount(j,k) = totcount(j,k) + 1
             end if
          end do
       end do
    end do

    where(totcount > 0) totbuffer = totbuffer/real(totcount)

    output = totbuffer  

  end subroutine rotave_resize_real8_sub

  subroutine rotave_real8_sub_nomiss(data, nfold)

    ! Rotationally average a double precision matrix with no missing 
    ! value defined. If output is a different size to data, this 
    ! will be done with re-sizing.

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(inout) :: data
    integer, intent(in)                               :: nfold

    ! Local variables
    real(rd_kind), dimension(size(data,1),size(data,2)) :: buffer, tot
    real(kind=rd_kind) :: step
    integer :: i

    ! Precompute the rotation to create an n-fold average
    step = (2.d0/real(nfold,rd_kind))*pi

    ! Initialise our output to the same as in the input
    tot = data

    ! Now add nfold minus 1 rotated copies on top ...
    do i = 1, nfold-1
       buffer = data
       call rotate_image(buffer,real(i,rd_kind)*step)
       tot = tot + buffer
    end do

    ! ... and divide through by nfold
    data = tot/real(nfold,kind=rd_kind)

  end subroutine rotave_real8_sub_nomiss

  subroutine rotave_resize_real8_sub_nomiss(data, nfold, output)

    ! Rotationally average a double precision matrix with no missing 
    ! value defined.

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    integer, intent(in)                             :: nfold
    real(kind=rd_kind), dimension(:,:), intent(out) :: output

    ! Local variables
    real(rd_kind), dimension(size(data,1),size(data,2)) :: inbuffer
    real(rd_kind), dimension(size(output,1),size(output,2)) :: outbuffer, tot
    real(kind=rd_kind) :: step
    integer :: i

    ! Precompute the rotation to create an n-fold average
    step = (2.d0/real(nfold,rd_kind))*pi

    ! Initialise our output to zero this time
    tot = 0.d0
    
    ! And then create nfold versions of the original, possibly
    ! with resizing ...
    do i = 0, nfold-1
       inbuffer = data
       call rotate_image(inbuffer,real(i,rd_kind)*step,outbuffer)
       tot = tot + outbuffer
    end do

    output = tot/real(nfold,kind=rd_kind)

  end subroutine rotave_resize_real8_sub_nomiss

!!$  subroutine rotave_real8_inplace(data, nfold, missingvalue, coverage)
!!$
!!$    ! Wrapper to allow in-place rotational averaging
!!$
!!$    ! Interface variables
!!$    real(kind=rd_kind), dimension(:,:), intent(inout) :: data
!!$    integer, intent(in)                               :: nfold
!!$    real(kind=rd_kind), intent(in), optional          :: missingvalue
!!$    real(kind=rs_kind), intent(in), optional          :: coverage
!!$
!!$    ! Local variables
!!$    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: output
!!$
!!$    if (present(missingvalue)) then
!!$       if (present(coverage)) then
!!$          call rotave_image(data, nfold, output, missingvalue, coverage)
!!$       else
!!$          call rotave_image(data, nfold, output, missingvalue)
!!$       end if
!!$    else
!!$       call rotave_image(data, nfold, output)
!!$    end if
!!$
!!$    data = output
!!$
!!$  end subroutine rotave_real8_inplace

  subroutine rotave_real_sub(data, nfold, output, missingvalue, coverage)

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in)  :: data
    integer, intent(in)                             :: nfold
    real(kind=rs_kind), dimension(:,:), intent(out) :: output
    real(kind=rs_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: localoutput

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind), coverage)
       else
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind))
       end if
    else
       call rotave_image(real(data,rd_kind), nfold, localoutput)
    end if

    output = real(localoutput,rs_kind)

  end subroutine rotave_real_sub

  subroutine rotave_integer_sub(data, nfold, output, missingvalue, coverage)

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    integer, intent(in)                      :: nfold
    integer, dimension(:,:), intent(out)     :: output
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(output,1),size(output,2)) :: localoutput

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind), coverage)
       else
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind))
       end if
    else
       call rotave_image(real(data,rd_kind), nfold, localoutput)
    end if

    output = nint(localoutput,rs_kind)

  end subroutine rotave_integer_sub

  subroutine rotave_real_inplace_sub(data, nfold, missingvalue, coverage)

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(inout) :: data
    integer, intent(in)                               :: nfold
    real(kind=rs_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: localoutput

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind), coverage)
       else
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind))
       end if
    else
       call rotave_image(real(data,rd_kind), nfold, localoutput)
    end if

    data = real(localoutput,rs_kind)

  end subroutine rotave_real_inplace_sub

  subroutine rotave_integer_inplace_sub(data, nfold, missingvalue, coverage)

    ! Interface variables
    integer, dimension(:,:), intent(inout)   :: data
    integer, intent(in)                      :: nfold
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: localoutput

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind), coverage)
       else
          call rotave_image(real(data,rd_kind), nfold, localoutput, real(missingvalue,rd_kind))
       end if
    else
       call rotave_image(real(data,rd_kind), nfold, localoutput)
    end if

    data = nint(localoutput,rs_kind)

  end subroutine rotave_integer_inplace_sub

  function rotave_real8_fn(data, nfold, missingvalue, coverage) result(rotimage)

    ! Rotate a real image, double precision angle

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in) :: data
    integer, intent(in)                            :: nfold
    real(kind=rd_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(data,nfold,rotimage,missingvalue,coverage)
       else
          call rotave_image(data,nfold,rotimage,missingvalue)
       end if
    else
       call rotave_image(data,nfold,rotimage)
    end if

  end function rotave_real8_fn

  function rotave_real_fn(data, nfold, missingvalue, coverage) result(rotimage)

    ! Rotate a real image, double precision angle

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in) :: data
    integer, intent(in)                            :: nfold
    real(kind=rs_kind), intent(in), optional       :: missingvalue
    real(kind=rs_kind), intent(in), optional       :: coverage

    ! Function result
    real(kind=rs_kind), dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(data,nfold,rotimage,missingvalue,coverage)
       else
          call rotave_image(data,nfold,rotimage,missingvalue)
       end if
    else
       call rotave_image(data,nfold,rotimage)
    end if

  end function rotave_real_fn
  
  function rotave_integer_fn(data, nfold, missingvalue, coverage) result(rotimage)

    ! Rotate a real image, double precision angle

    ! Interface variables
    integer, dimension(:,:), intent(in)      :: data
    integer, intent(in)                      :: nfold
    integer, intent(in), optional            :: missingvalue
    real(kind=rs_kind), intent(in), optional :: coverage

    ! Function result
    integer, dimension(size(data,1),size(data,2)) :: rotimage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call rotave_image(data,nfold,rotimage,missingvalue,coverage)
       else
          call rotave_image(data,nfold,rotimage,missingvalue)
       end if
    else
       call rotave_image(data,nfold,rotimage)
    end if

  end function rotave_integer_fn

!!$  function centre_real8 (image) result(delta)
!!$
!!$    use pnm_class
!!$
!!$    real(kind=rd_kind), dimension(:,:), intent(in) :: image
!!$
!!$    real(rs_kind) :: delta(2)
!!$
!!$    type(pnm_object) :: pnm
!!$
!!$    complex(rd_kind), dimension(size(image,1),size(image,2)) :: c, Q, c2fold
!!$    complex(rd_kind), dimension(1+size(image,1)/2,size(image,2)) :: halfQ
!!$    real(rd_kind), dimension(size(image,1),size(image,2)) :: phase_of_Q
!!$    ! real(rd_kind), dimension(1+size(image,1)/2,size(image,2)) :: phase_of_Q
!!$    complex(rd_kind), dimension(size(image,1),size(image,2)) :: left_singular
!!$    ! complex(rd_kind), dimension(1+size(image,1)/2,size(image,2)) :: left_singular
!!$    complex(rd_kind), dimension(size(image,2),size(image,2)) :: right_singular
!!$    real(rd_kind), dimension(3*size(image)) :: lwork
!!$    real(rd_kind), dimension(6*size(image)) :: rwork
!!$    real(rd_kind), dimension(min(size(image,1),size(image,2))) :: singular_values
!!$    ! real(rd_kind), dimension(min(1+size(image,1)/2,size(image,2))) :: singular_values
!!$    real(rd_kind), dimension(size(image))   :: x, y, sig
!!$    real(rd_kind) :: corr
!!$
!!$    integer :: i, l, m, status, ldim, lensav, limits(2), iwork(8*min(size(image,1),size(image,2)))
!!$
!!$    l = size(image,1)
!!$    m = size(image,2)
!!$
!!$    c = fft(image)
!!$    c2fold = fft(image(l:1:-1,m:1:-1))
!!$
!!$    pnm = nint(atan2(aimag(c),real(c))*(65535./pi))
!!$    call write(pnm,'cphase.pgm')
!!$
!!$    ! Apply a circular shift so that the fourier spectrum look like what 
!!$    ! engineers expect
!!$    c = cshift(cshift(c,l/2,1),m/2,2)
!!$    c2fold = cshift(cshift(c2fold,l/2,1),m/2,2)
!!$
!!$    ! Calculate the phase correlation matrix
!!$    Q = (c*conjg(c2fold))/abs(c2fold*conjg(c2fold))
!!$    halfQ = Q(1:(1+m/2),:)
!!$
!!$    phase_of_Q = atan2(aimag(Q),real(Q))
!!$    pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    call write(pnm,'q.pgm')
!!$
!!$    print *,'about to do SVD ',size(Q,1),size(Q,2)
!!$
!!$    ! Call the svd routine from LAPACK
!!$    ! 'A' = Compute all eigenvalues and eigenvectors [There is currently no option to calculate
!!$    ! a subset of these, only min(M,N), which is bog all good for a symmetric matrix]
!!$    ! call cgesvd('A', 'A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,2), right_singular, size(Q,2), lwork, size(lwork), rwork, size(rwork), status) 
!!$    ! call zgesvd('A', 'A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,2), right_singular, size(Q,2), lwork, size(lwork), rwork, status) 
!!$    call zgesdd('A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,1), right_singular, size(Q,2), lwork, size(lwork), rwork, iwork, status) 
!!$    ! call zgesdd('A', size(halfQ,1), size(halfQ,2), halfQ, size(halfQ,1), singular_values, left_singular, size(halfQ,1), right_singular, size(halfQ,2), lwork, size(lwork), rwork, iwork, status) 
!!$
!!$    print *,'done SVD'
!!$
!!$    phase_of_Q = atan2(aimag(left_singular),real(left_singular))
!!$    ! pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    ! call write(pnm,'left.pgm')
!!$    phase_of_Q(:,1) = unwrap_phase(phase_of_Q(:,1))
!!$    write(1,'(F12.6)') phase_of_Q(:,1)
!!$
!!$    limits = (/ nint(0.2 * size(phase_of_Q,1)), nint(0.8 * size(phase_of_Q,1)) /)
!!$
!!$    print *,limits
!!$
!!$    corr = 0.0
!!$
!!$    ! delta(1) = (sumsquares((/(real(i,rd_kind), i=limits(1), limits(2))/),phase_of_Q(limits(1):limits(2),1))/sumsquares((/(real(i,rd_kind), i=limits(1), limits(2))/)))*(real(l)/(2.*pi))
!!$    ! delta(1) = find_slope((/(real(i,rd_kind), i=1, size(phase_of_Q))/),phase_of_Q(:,1), corr)*(real(l)/(2.*pi))
!!$    ! delta(1) = find_slope((/(real(i,rd_kind), i=limits(1), limits(2))/),phase_of_Q(limits(1):limits(2),1), corr)*(real(l)/(2.*pi))
!!$    ! delta(1) = slope((/(real(i,rd_kind), i=1, size(phase_of_Q,1))/),phase_of_Q(:,1))*(real(l)/(2.*pi))
!!$    ! delta(1) = slope((/(real(i,rd_kind), i=1, size(phase_of_Q,1)/2)/),phase_of_Q(1:size(phase_of_Q,1)/2,1))
!!$    ! delta(1) = (delta(1) + slope((/(real(i,rd_kind), i=1, size(phase_of_Q,1)/2)/),phase_of_Q(size(phase_of_Q,1)/2:,1)))*(real(l)/(4.*pi))
!!$    ! delta(1) = slope( (/(real(i,rd_kind), i=1, int(0.8*size(phase_of_Q,1)/2))/), phase_of_Q(int(0.2*size(phase_of_Q,1)/2):size(phase_of_Q,1),1) )
!!$    ! delta(1) = (delta(1) + slope( (/ (real(i,rd_kind), i=1, int(0.8*size(phase_of_Q,1)/2.)) /), phase_of_Q(size(phase_of_Q,1)/2:int(0.8*size(phase_of_Q,1)),1) ))*(real(l)/(4.*pi))
!!$
!!$    delta(1) = slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,1)/2)/), phase_of_Q(1:size(phase_of_Q,1)/2,1) )*(real(l)/(2.*pi))
!!$    
!!$    limits = (/ nint(0.2 * size(phase_of_Q,2)), nint(0.8 * size(phase_of_Q,2)) /)
!!$
!!$    phase_of_Q = atan2(aimag(right_singular),real(right_singular))
!!$    ! pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    ! call write(pnm,'right.pgm')
!!$    phase_of_Q(1,:) = unwrap_phase(phase_of_Q(1,:))
!!$    write(2,'(F12.6)') phase_of_Q(1,:)
!!$    
!!$    corr = 0.0
!!$
!!$    ! delta(2) = (sumsquares((/(real(i,rd_kind), i=limits(1), limits(2))/),phase_of_Q(1,limits(1):limits(2))) / sumsquares((/(real(i,rd_kind), i=limits(1), limits(2))/)))*(real(m)/(2.*pi))
!!$    ! delta(2) = find_slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,2))/), phase_of_Q(1,:), corr )*(real(m)/(2.*pi))
!!$    ! delta(2) = slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,2))/), phase_of_Q(1,:))*(real(m)/(2.*pi))
!!$    ! delta(2) = slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,2)/2)/), phase_of_Q(1,1:size(phase_of_Q,2)/2) )
!!$    ! delta(2) = ( delta(2) + slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,2)/2)/), phase_of_Q(1,size(phase_of_Q,2)/2:)) )*(real(m)/(4.*pi))
!!$
!!$    delta(2) = slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,2)/2)/), phase_of_Q(1,1:size(phase_of_Q,2)/2) )*(real(m)/(2.*pi))
!!$
!!$    ! We have 2-fold rotated the image and then calculated the 
!!$    ! offset. To get the vector to translate this back to the
!!$    ! centre we need to divide by 2
!!$    delta = delta / 2.
!!$
!!$    print *,'x displacement: ',delta(1)
!!$    print *,'y displacement: ', delta(2)
!!$
!!$    return
!!$
!!$  end function centre_real8
!!$
!!$  function translational_align_real8 (imageone, imagetwo) result(delta)
!!$
!!$    use pnm_class
!!$
!!$    real(kind=rd_kind), dimension(:,:), intent(in) :: imageone, imagetwo
!!$
!!$    real(rs_kind) :: delta(2)
!!$
!!$    type(pnm_object) :: pnm
!!$
!!$    complex(rd_kind), dimension(size(imageone,1),size(imageone,2)) :: c_one, Q, c_two
!!$    complex(rd_kind), dimension(1+size(imageone,1)/2,size(imageone,2)) :: halfQ
!!$    real(rd_kind), dimension(size(imageone,1),size(imageone,2)) :: phase_of_Q
!!$    ! real(rd_kind), dimension(1+size(imageone,1)/2,size(imageone,2)) :: phase_of_Q
!!$    complex(rd_kind), dimension(size(imageone,1),size(imageone,2)) :: left_singular
!!$    ! complex(rd_kind), dimension(1+size(imageone,1)/2,size(imageone,2)) :: left_singular
!!$    complex(rd_kind), dimension(size(imageone,2),size(imageone,2)) :: right_singular
!!$    real(rd_kind), dimension(3*size(imageone)) :: lwork
!!$    real(rd_kind), dimension(6*size(imageone)) :: rwork
!!$    real(rd_kind), dimension(min(size(imageone,1),size(imageone,2))) :: singular_values
!!$    ! real(rd_kind), dimension(min(1+size(imageone,1)/2,size(imageone,2))) :: singular_values
!!$    real(rd_kind), dimension(size(imageone))   :: x, y, sig
!!$    real(rd_kind) :: corr
!!$
!!$    integer :: i, l, m, status, ldim, lensav, limits(2), iwork(8*min(size(imageone,1),size(imageone,2)))
!!$
!!$    l = size(imageone,1)
!!$    m = size(imageone,2)
!!$
!!$    c_one = fft(imageone)
!!$    c_two = fft(imagetwo)
!!$
!!$    pnm = nint(atan2(aimag(c_one),real(c_one))*(65535./pi))
!!$    call write(pnm,'cphase.pgm')
!!$
!!$    ! Apply a circular shift so that the fourier spectrum look like what 
!!$    ! engineers expect
!!$    c_one = cshift(cshift(c_one,l/2,1),m/2,2)
!!$    c_two = cshift(cshift(c_two,l/2,1),m/2,2)
!!$
!!$    ! Calculate the phase correlation matrix
!!$    Q = (c_one*conjg(c_two))/abs(c_two*conjg(c_two))
!!$    halfQ = Q(1:(1+m/2),:)
!!$
!!$    phase_of_Q = atan2(aimag(Q),real(Q))
!!$    pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    call write(pnm,'q.pgm')
!!$
!!$    print *,'about to do SVD ',size(Q,1),size(Q,2)
!!$
!!$    ! Call the svd routine from LAPACK
!!$    ! 'A' = Compute all eigenvalues and eigenvectors [There is currently no option to calculate
!!$    ! a subset of these, only min(M,N), which is bog all good for a symmetric matrix]
!!$    ! call cgesvd('A', 'A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,2), right_singular, size(Q,2), lwork, size(lwork), rwork, size(rwork), status) 
!!$    ! call zgesvd('A', 'A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,2), right_singular, size(Q,2), lwork, size(lwork), rwork, status) 
!!$    call zgesdd('A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,1), right_singular, size(Q,2), lwork, size(lwork), rwork, iwork, status) 
!!$    ! call zgesdd('A', size(halfQ,1), size(halfQ,2), halfQ, size(halfQ,1), singular_values, left_singular, size(halfQ,1), right_singular, size(halfQ,2), lwork, size(lwork), rwork, iwork, status) 
!!$
!!$    print *,'done SVD'
!!$
!!$    phase_of_Q = atan2(aimag(left_singular),real(left_singular))
!!$    ! pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    ! call write(pnm,'left.pgm')
!!$    phase_of_Q(:,1) = unwrap_phase(phase_of_Q(:,1))
!!$    write(1,'(F12.6)') phase_of_Q(:,1)
!!$
!!$    limits = (/ nint(0.2 * size(phase_of_Q,1)), nint(0.8 * size(phase_of_Q,1)) /)
!!$
!!$    print *,limits
!!$
!!$    corr = 0.0
!!$
!!$    delta(1) = slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,1)/2)/), phase_of_Q(1:size(phase_of_Q,1)/2,1) )*(real(l)/(2.*pi))
!!$    
!!$    limits = (/ nint(0.2 * size(phase_of_Q,2)), nint(0.8 * size(phase_of_Q,2)) /)
!!$
!!$    phase_of_Q = atan2(aimag(right_singular),real(right_singular))
!!$    ! pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    ! call write(pnm,'right.pgm')
!!$    phase_of_Q(1,:) = unwrap_phase(phase_of_Q(1,:))
!!$    write(2,'(F12.6)') phase_of_Q(1,:)
!!$    
!!$    corr = 0.0
!!$
!!$    delta(2) = slope( (/(real(i,rd_kind), i=1, size(phase_of_Q,2)/2)/), phase_of_Q(1,1:size(phase_of_Q,2)/2) )*(real(m)/(2.*pi))
!!$
!!$    print *,'x displacement: ',delta(1)
!!$    print *,'y displacement: ', delta(2)
!!$
!!$    return
!!$
!!$  end function translational_align_real8

  function slope(x, y)

    ! Interface variables
    real(rd_kind), dimension(:), intent(in) :: x, y

    ! Return value
    real(rd_kind) :: slope

    ! Local variables
    real(rd_kind) :: ssx, ssy, ssxy, intercept, correlation, new_correlation, residuals(size(y))
    real(rd_kind) :: tolerance, maxcorrelation
    integer :: lb, ub, step

    ! Initialise the totals
    ! tolerance = -0.001
    tolerance = 0.0

    ! step = nint(real(size(x))/10.)
    step = nint(real(size(x))/20.)
    if (step < 1) step = 1

    lb = nint(real(size(x))/2. - step)
    ub = nint(real(size(x))/2. + step)

    ssx = sumsquares(x(lb:ub))
    ssy = sumsquares(y(lb:ub))
    ssxy = sumsquares(x(lb:ub),y(lb:ub))
       
    correlation = (ssxy**2)/(ssx*ssy)
    maxcorrelation = correlation

    slope = ssxy/ssx
    print *,'Size of x ',size(x)
    print *,'lb/ub ',lb,ub
    print *,'Slope ',slope * (real(size(x))/pi)
    print *,'Correlation ',correlation

    ! intercept = mean(y) - slope * mean(x)

    ! residuals = y - (slope*x + intercept)

    ! new_correlation = sqrt((ssxy)**2/(ssx*ssy))

!!$    do 
!!$       ! Add 10% more points in to the fit and recalculate the totals
!!$       ssx = sumsquares(x(lb-step:ub+step))
!!$       ssy = sumsquares(y(lb-step:ub+step))
!!$       ssxy = sumsquares(x(lb-step:ub+step),y(lb-step:ub+step))
!!$       
!!$       new_correlation = (ssxy**2)/(ssx*ssy)
!!$       print *,'New Correlation ',new_correlation
!!$
!!$       maxcorrelation = max(maxcorrelation, new_correlation)
!!$
!!$       ! If we have a significant drop in the 
!!$       ! 'coefficient of determination' (R**2) then
!!$       ! stop what we're doing and use the previous 
!!$       ! estimation of the slope
!!$       print *,(new_correlation - maxcorrelation)
!!$       if ( (new_correlation - maxcorrelation) < tolerance) then
!!$          print *,'Exiting now'
!!$          print *,'Slope ',slope * (real(size(x))/(2.*pi))
!!$          ! slope = zap_outliers(x(lb:ub), y(lb:ub), correlation)
!!$          ! print *,'Slope ',slope * (real(size(x))/(2.*pi))
!!$          exit
!!$       end if
!!$
!!$       ! Update the current upper and lower bounds
!!$       ! and the new slope
!!$       lb = lb - step
!!$       ub = ub + step
!!$       if (ub > size(x)) exit
!!$       if (lb < 1) exit
!!$       slope = ssxy/ssx
!!$       print *,'lb/ub ',lb,ub
!!$       print *,'Slope ',slope * (real(size(x))/(2.*pi))
!!$    end do

    do 
       ! Add 10% more points in to the fit and recalculate the totals
       ssx = sumsquares(x(lb-step:ub))
       ssy = sumsquares(y(lb-step:ub))
       ssxy = sumsquares(x(lb-step:ub),y(lb-step:ub))
       
       new_correlation = (ssxy**2)/(ssx*ssy)
       print *,'New Correlation ',new_correlation

       maxcorrelation = max(maxcorrelation, new_correlation)

       ! If we have a significant drop in the 
       ! 'coefficient of determination' (R**2) then
       ! stop what we're doing and use the previous 
       ! estimation of the slope
       print *,(new_correlation - maxcorrelation)
       if ( (new_correlation - maxcorrelation) < tolerance) then
          print *,'Exiting now'
          print *,'Slope ',slope * (real(size(x))/pi)
          ! slope = zap_outliers(x(lb:ub), y(lb:ub), correlation)
          ! print *,'Slope ',slope * (real(size(x))/(2.*pi))
          exit
       end if

       ! Update the current upper and lower bounds
       ! and the new slope
       lb = lb - step
       if (lb < 1) exit
       slope = ssxy/ssx
       print *,'lb/ub ',lb,ub
       print *,'Slope ',slope * (real(size(x))/pi)
    end do

    do 
       ! Add 10% more points in to the fit and recalculate the totals
       ssx = sumsquares(x(lb:ub+step))
       ssy = sumsquares(y(lb:ub+step))
       ssxy = sumsquares(x(lb:ub+step),y(lb:ub+step))
       
       new_correlation = (ssxy**2)/(ssx*ssy)
       print *,'New Correlation ',new_correlation

       maxcorrelation = max(maxcorrelation, new_correlation)

       ! If we have a significant drop in the 
       ! 'coefficient of determination' (R**2) then
       ! stop what we're doing and use the previous 
       ! estimation of the slope
       print *,(new_correlation - maxcorrelation)
       if ( (new_correlation - maxcorrelation) < tolerance) then
          print *,'Exiting now'
          print *,'Slope ',slope * (real(size(x))/pi)
          ! slope = zap_outliers(x(lb:ub), y(lb:ub), correlation)
          ! print *,'Slope ',slope * (real(size(x))/(2.*pi))
          exit
       end if

       ! Update the current upper and lower bounds
       ! and the new slope
       ub = ub + step
       if (ub > size(x)) exit
       slope = ssxy/ssx
       print *,'lb/ub ',lb,ub
       print *,'Slope ',slope * (real(size(x))/pi)
    end do

  end function slope

  function intercept(x, y, slope)

    ! Interface variables
    real(rd_kind), dimension(:), intent(in) :: x, y
    real(rd_kind), intent(in) :: slope

    ! Return value
    real(rd_kind) :: intercept

    intercept = mean(y) - slope * mean(x)

  end function intercept

  recursive function zap_outliers(x, y, correlation) result(slope)

    ! Interface variables
    real(rd_kind), dimension(:), intent(in) :: x, y
    real(rd_kind), intent(inout) :: correlation

    ! Return value
    real(rd_kind) :: slope

    ! Local variables
    real(rd_kind) :: ssx, ssy, ssxy, intercept, new_correlation, residuals(size(y))
    logical :: kept(size(y))
    real(rd_kind) :: tolerance

    tolerance = 0.10 / real(size(x),rd_kind)

    kept = .TRUE.

    ssx = sumsquares(x)
    ssy = sumsquares(y)
    ssxy = sumsquares(x,y)

    slope = ssxy/ssx

    intercept = mean(y) - slope * mean(x)

    residuals = y - (slope*x + intercept)

    new_correlation = sqrt((ssxy)**2/(ssx*ssy))

    print *,'Slope ',slope * (real(200)/pi)
    print *,'Correlation ',new_correlation
    if ( (new_correlation - correlation) / correlation > tolerance) then
       kept(maxloc(abs(residuals))) = .FALSE.
       print *,'Discarding number ',maxloc(abs(residuals))
       slope = zap_outliers(pack(x,mask=kept), pack(y,mask=kept), new_correlation)
    else
       print *,'Not discarding number ',maxloc(abs(residuals))
    end if

  end function zap_outliers

!!$  function centre_real8_old (input) result(delta)
!!$
!!$    use pnm_class
!!$
!!$    real(kind=rd_kind), dimension(:,:), intent(inout) :: input
!!$
!!$    real(rs_kind) :: delta(2)
!!$
!!$    type(pnm_object) :: pnm
!!$
!!$    complex(rs_kind), dimension(size(input,1),size(input,2)) :: c, Q, c2fold
!!$    real(rs_kind), dimension(size(input,1),size(input,2)) :: phase_of_Q
!!$    complex(rs_kind), dimension(size(input,1),size(input,2)) :: left_singular
!!$    complex(rs_kind), dimension(size(input,2),size(input,2)) :: right_singular
!!$    complex(rs_kind), dimension(3*size(input)) :: lwork, rwork
!!$    real(rs_kind), dimension(size(input,2)) :: singular_values
!!$    real(rs_kind), dimension(size(input))   :: x, y, sig
!!$
!!$    integer :: i, l, m, lenwrk, status, ldim, lensav
!!$    real :: work(size(input)*2)
!!$
!!$    real, allocatable, dimension (:) :: wsave
!!$
!!$    l = size(c,1)
!!$    m = size(c,2)
!!$    lenwrk = size(work)
!!$
!!$    ! We convert the real values to complex ones so we can call the 
!!$    ! appropriate fft routine. We could call one that takes real
!!$    ! values but then I'd have to unpack the stupid half-complex
!!$    ! format. Praps I'll do that in a later version as it would mean
!!$    ! a faster FFT.
!!$
!!$    ! c is the complex version of the input (with zero imaginary values)
!!$    c = cmplx(input(:l,:m))
!!$    ! c2fold is the 2-fold rotated version of c
!!$    c2fold = cmplx(input(l:1:-1,m:1:-1))
!!$
!!$    lensav = 2*(l+m) + int(log(real(l))) + int(log(real(m))) + 8
!!$
!!$    allocate ( wsave(1:lensav) )
!!$
!!$    call cfft2i ( l, m, wsave, lensav, status )
!!$    !
!!$    !  Compute the FFT coefficients.
!!$    !
!!$    ldim = l
!!$
!!$    print *,'about to do FFT'
!!$
!!$    call cfft2f ( ldim, l, m, c, wsave, lensav, work, lenwrk, status )
!!$    call cfft2f ( ldim, l, m, c2fold, wsave, lensav, work, lenwrk, status )
!!$
!!$    ! Apply a circular shift so that the fourier spectrum look like what 
!!$    ! engineers expect
!!$    c = cshift(cshift(c,l/2,1),m/2,2)
!!$    c2fold = cshift(cshift(c2fold,l/2,1),m/2,2)
!!$
!!$    ! Calculate the phase correlation matrix
!!$    Q = (c*conjg(c2fold))/abs(c2fold*conjg(c2fold))
!!$
!!$    phase_of_Q = atan2(aimag(Q),real(Q))
!!$    pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    call write(pnm,'q.pgm')
!!$
!!$    print *,'about to do SVD ',size(Q,1),size(Q,2)
!!$
!!$    ! Call the svd routine from LAPACK
!!$    ! 'A' = Compute all eigenvalues and eigenvectors [There is currently no option to calculate
!!$    ! a subset of these, only min(M,N), which is bog all good for a symmetric matrix]
!!$    call cgesvd('A', 'A', size(Q,1), size(Q,2), Q, size(Q,1), singular_values, left_singular, size(Q,2), right_singular, size(Q,2), lwork, size(lwork), rwork, size(rwork), status) 
!!$
!!$    print *,'done SVD'
!!$
!!$    phase_of_Q = atan2(aimag(left_singular),real(left_singular))
!!$    ! pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    ! call write(pnm,'left.pgm')
!!$ !!   phase_of_Q(:,1) = unwrap_phase(phase_of_Q(:,1))
!!$    ! write(1,*) phase_of_Q(:,1)
!!$
!!$ !!   delta(1) = (sumsquares((/(real(i), i=10, 110)/),phase_of_Q(10:110,1))/sumsquares((/(real(i), i=10, 110)/)))*(l/(2*pi))
!!$    print *,'x displacement: ',delta(1)
!!$
!!$    phase_of_Q = atan2(-aimag(right_singular),real(right_singular))
!!$    ! pnm = nint(phase_of_Q*(65535*maxval(phase_of_Q)))
!!$    ! call write(pnm,'right.pgm')
!!$!!    phase_of_Q(1,:) = unwrap_phase(phase_of_Q(1,:))
!!$    ! write(2,*) phase_of_Q(1,:)
!!$
!!$!!    delta(2) = (sumsquares((/(real(i), i=0, m-1)/),phase_of_Q(1,:)) / sumsquares((/(real(i), i=0, m-1)/)))*(m/(2*pi))
!!$    print *,'y displacement: ', delta(2)
!!$
!!$    delta = 0.
!!$
!!$    deallocate ( wsave )
!!$
!!$    return
!!$
!!$  end function centre_real8_old

!!$  function centre_integer_fn (data, missingvalue) result(shift)
!!$
!!$    integer, dimension(:,:), intent(in) :: data
!!$    integer, intent(in), optional       :: missingvalue
!!$
!!$    ! Function result
!!$    real, dimension(2) :: shift
!!$
!!$    ! Local variables
!!$    integer, dimension(size(data,1),size(data,2)) :: data_twofold
!!$    integer :: width, height, window_width, window_height
!!$    real :: data_median, window_median
!!$
!!$    width = size(data,1)
!!$    height = size(data,1)
!!$
!!$    ! Make a two fold related version of the data
!!$    data_twofold = data(width:1:-1,height:1:-1) 
!!$
!!$    data_median = median(data)
!!$
!!$    ! Now we'll try and match them up
!!$    do
!!$    end do
!!$
!!$  end function centre_integer_fn

  function unwrap_linear_phase(phase) result(pout)

    ! Phase unwrapping where we make the assumption that the
    ! unwrapped phase is linear
    
    ! Interface variables
    real(rd_kind), intent(in) :: phase(:)

    ! Output type
    real(rd_kind)             :: pout(size(phase))

    ! Local variables
    real(rd_kind) :: diffs(size(phase)) 
    real(rd_kind) :: totalphase
    integer :: i
    real(rd_kind) :: ssx, ssy, ssxy

    ! Create an array of the differences between adjacent
    ! elements in phase (eoshifts to the left). The last 
    ! element of diffs is nonsense and never accessed.
    diffs = phase - eoshift(phase,1)

    write(10,'(F12.6)') diffs
    where (abs(diffs) < 0.9*pi) diffs = 0.

    !where (abs(diffs) < 0.9*pi) 
    !   diffs = 0
    !elsewhere
    !where (abs(diffs) > 1.5*pi) diffs = sign(real(2.*pi,rd_kind),diffs)
    !where ((abs(diffs) < 1.5*pi) .and. diffs /= 0.) diffs = sign(real(pi,rd_kind),diffs)

    where (abs(diffs) < pi) 
       diffs = 0
    elsewhere
       diffs = sign(real(2.*pi,rd_kind),diffs)
    end where

    write(11,'(F12.6)') diffs
    
    totalphase = 0.
    pout(1) = phase(1)

    ! ssx = 1
    ! ssy = ssy + sumsquares(y(lb:ub))
    ! ssxy = ssxy + sumsquares(x(lb:ub),y(lb:ub))
       

    do i = 2, size(phase)
       totalphase = totalphase + diffs(i-1)
       pout(i) = phase(i) + totalphase
    end do
    
  end function unwrap_linear_phase

  function unwrap_phase(phase) result(pout)

    use variable_array
    
    ! Interface variables
    real(rd_kind), intent(in) :: phase(:)

    ! Output type
    real(rd_kind)             :: pout(size(phase))

    ! Local variables
    real(rd_kind) :: diffs(size(phase)) 
    real(rd_kind) :: totalphase, meandiff, sddiff
    integer :: i, index
    integer, dimension(:), pointer :: indices => null()

    ! Create an array of the differences between adjacent
    ! elements in phase (eoshifts to the left). The last 
    ! element of diffs is nonsense and never accessed.
    diffs = phase - eoshift(phase,1)

    write(10,'(F12.6)') diffs
    ! where (abs(diffs) < 0.9*pi) diffs = 0.

    !where (abs(diffs) < 0.9*pi) 
    !   diffs = 0
    !elsewhere
    !where (abs(diffs) > 1.5*pi) diffs = sign(real(2.*pi,rd_kind),diffs)
    !where ((abs(diffs) < 1.5*pi) .and. diffs /= 0.) diffs = sign(real(pi,rd_kind),diffs)

    where (abs(diffs) < pi) 
       diffs = 0
    elsewhere
       diffs = sign(real(2.*pi,rd_kind),diffs)
    end where

    ! We think a discontinuity occurs where there is a difference that is more than
    ! 3 sigma different from the mean

    meandiff = mean(diffs)
    sddiff = stddev(diffs)

    ! Make an array of indices where there are potential discontinuities
    i = push(indices, pack((/(i,i=1,size(diffs))/), mask=abs(diffs) > 3. * sddiff))

    ! if (size(indices) > 0) then
    if (.FALSE.) then
       do i = 1, size(indices)
          index = indices(i)
          print *,'index ',index
          do 
             ! Keep looping while adding multiples of pi to our
             ! difference brings it closer to the mean of differences
             ! if (abs(diffs(indices(i))) = sign(pi,diffs(indices(i)))
             if ( abs(diffs(index) - meandiff) < abs((diffs(index) - sign(pi,diffs(index))) - meandiff) ) exit
             diffs(index) = diffs(index) - sign(pi,diffs(index))
          end do
       end do

       ! Free up any memory associated with indices
       i =  splice(indices,0)

    end if

    write(11,'(F12.6)') diffs
    
    totalphase = 0.
    pout(1) = phase(1)

    do i = 2, size(phase)
       totalphase = totalphase + diffs(i-1)
       pout(i) = phase(i) + totalphase
    end do

  end function unwrap_phase


  real function sumsquares(x,y)

    use statistics, only: mean

    real(rd_kind), intent(in) :: x(:)
    real(rd_kind), intent(in), optional :: y(:)

    if (present(y)) then
       sumsquares = sum((x - mean(x))*(y - mean(y)))
    else
       sumsquares = sum((x - mean(x))**2)
    end if

  end function sumsquares

  
  ! Image translation

  subroutine translate_image_real8_sub(data, delta, output, missingvalue, coverage)

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    real(kind=rd_kind), dimension(2), intent(in)    :: delta
    real(kind=rd_kind), dimension(:,:), intent(out) :: output
    real(kind=rd_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(3,(size(output,1)+1)*(size(output,2)+1)) :: vertices
    integer :: i, nx, ny, newnx, newny

    if (debug) print *,'in translate_image'

    nx = size(data,1)
    ny = size(data,2)
    newnx = size(output,1)
    newny = size(output,2)

    if (debug) print *,nx,ny,newnx,newny

    ! Reproduce the current vertices (runs from 0 -> nx, 0 -> ny)
    vertices(1,:) = pack(spread(real((/ (i, i = 0, newnx) /),rd_kind) * &
         (real(nx,rd_kind)/real(newnx,rd_kind)), 2, newny+1),.TRUE.)
    vertices(2,:) = pack(spread(real((/ (i, i = 0, newny) /),rd_kind) * &
         (real(ny,rd_kind)/real(newny,rd_kind)), 1, newnx+1),.TRUE.)
    vertices(3,:) = 0.d0

    if (debug) print *,'Original vertices'
    if (debug) print '(2F0.4)',vertices(1:2,1:10)

    ! Translate the vertices by delta
    vertices(1,:) = vertices(1,:) + delta(1)
    vertices(2,:) = vertices(2,:) + delta(2)
    
    if (debug) print '(2F0.4)',vertices(1:2,:)

    if (present(missingvalue)) then
       if (present(coverage)) then
          call polySamp(data, vertices(1,:), vertices(2,:), output, missingvalue, coverage=coverage)
       else
          call polySamp(data, vertices(1,:), vertices(2,:), output, missingvalue)
       end if
    else
       call polySamp(data, vertices(1,:), vertices(2,:), output)
    end if

  end subroutine translate_image_real8_sub

  subroutine translate_image_realdelta(data, delta, output, missingvalue, coverage)

    ! Translate a double precision matrix, but single precision delta

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    real(kind=rs_kind), dimension(2), intent(in)    :: delta
    real(kind=rd_kind), dimension(:,:), intent(out) :: output
    real(kind=rd_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,real(delta,rd_kind),output,missingvalue,coverage)
       else
          call translate_image(data,real(delta,rd_kind),output,missingvalue)
       end if
    else
       call translate_image(data,real(delta,rd_kind),output)
    end if

  end subroutine translate_image_realdelta
  
  subroutine translate_image_real_wrap(data, delta, output, missingvalue, coverage)

    ! Translate a single precision matrix

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in)  :: data
    real(kind=rs_kind), dimension(2), intent(in)    :: delta
    real(kind=rs_kind), dimension(:,:), intent(out) :: output
    real(kind=rs_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(real(data,rd_kind),real(delta,rd_kind),local_output,real(missingvalue,rd_kind),coverage)
       else
          call translate_image(real(data,rd_kind),real(delta,rd_kind),local_output,real(missingvalue,rd_kind))
       end if
    else
       call translate_image(real(data,rd_kind),real(delta,rd_kind),local_output)
    end if

    output = real(local_output,rs_kind)

  end subroutine translate_image_real_wrap

  subroutine translate_image_integer_wrap(data, delta, output, missingvalue, coverage)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    integer, dimension(:,:), intent(in)          :: data
    real(kind=rd_kind), dimension(2), intent(in) :: delta
    integer, dimension(:,:), intent(out)         :: output
    integer, intent(in), optional                :: missingvalue
    real(kind=rs_kind), intent(in), optional     :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(real(data,rd_kind),delta,local_output,real(missingvalue,rd_kind),coverage)
       else
          call translate_image(real(data,rd_kind),delta,local_output,real(missingvalue,rd_kind))
       end if
    else
       call translate_image(real(data,rd_kind),delta,local_output)
    end if

    output = nint(local_output)

  end subroutine translate_image_integer_wrap


  ! In place translation

  subroutine translate_inplace(data, delta, missingvalue, coverage)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(inout) :: data
    real(kind=rd_kind), dimension(2), intent(in)      :: delta
    real(kind=rd_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,delta,data,missingvalue,coverage)
       else
          call translate_image(data,delta,data,missingvalue)
       end if
    else
       call translate_image(data,delta,data)
    end if

  end subroutine translate_inplace

  subroutine translate_inplace_realdelta(data, delta, missingvalue, coverage)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(inout) :: data
    real(kind=rs_kind), dimension(2), intent(in)      :: delta
    real(kind=rd_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,real(delta,rd_kind),data,missingvalue,coverage)
       else
          call translate_image(data,real(delta,rd_kind),data,missingvalue)
       end if
    else
       call translate_image(data,real(delta,rd_kind),data)
    end if

  end subroutine translate_inplace_realdelta

  subroutine translate_inplace_real_wrap(data, delta, missingvalue, coverage)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(inout) :: data
    real(kind=rs_kind), dimension(2), intent(in)      :: delta
    real(kind=rs_kind), intent(in), optional          :: missingvalue
    real(kind=rs_kind), intent(in), optional          :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(real(data,rd_kind),delta,local_output,real(missingvalue,rd_kind),coverage)
       else
          call translate_image(real(data,rd_kind),delta,local_output,real(missingvalue,rd_kind))
       end if
    else
       call translate_image(real(data,rd_kind),delta,local_output)
    end if

    data = real(local_output,rs_kind)

  end subroutine translate_inplace_real_wrap

  subroutine translate_inplace_integer_wrap(data, delta, missingvalue, coverage)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    integer, dimension(:,:), intent(inout)       :: data
    real(kind=rd_kind), dimension(2), intent(in) :: delta
    integer, intent(in), optional                :: missingvalue
    real(kind=rs_kind), intent(in), optional     :: coverage

    ! Local variables
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: local_output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(real(data,rd_kind),delta,local_output,real(missingvalue,rd_kind),coverage)
       else
          call translate_image(real(data,rd_kind),delta,local_output,real(missingvalue,rd_kind))
       end if
    else
       call translate_image(real(data,rd_kind),delta,local_output)
    end if

    data = nint(local_output)

  end subroutine translate_inplace_integer_wrap

  
  ! Functional interface

  function translate_function(data, delta, missingvalue, coverage) result(output)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    real(kind=rd_kind), dimension(2), intent(in)    :: delta
    real(kind=rd_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Function return value
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,delta,output,missingvalue,coverage)
       else
          call translate_image(data,delta,output,missingvalue)
       end if
    else
       call translate_image(data,delta,output)
    end if

  end function translate_function

  function translate_function_realdelta(data, delta, missingvalue, coverage) result(output)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    real(kind=rd_kind), dimension(:,:), intent(in)  :: data
    real(kind=rs_kind), dimension(2), intent(in)    :: delta
    real(kind=rd_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Function return value
    real(kind=rd_kind), dimension(size(data,1),size(data,2)) :: output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,real(delta,rd_kind),output,missingvalue,coverage)
       else
          call translate_image(data,real(delta,rd_kind),output,missingvalue)
       end if
    else
       call translate_image(data,real(delta,rd_kind),output)
    end if

  end function translate_function_realdelta

  function translate_function_real_wrap(data, delta, missingvalue, coverage) result(output)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    real(kind=rs_kind), dimension(:,:), intent(in)  :: data
    real(kind=rs_kind), dimension(2), intent(in)    :: delta
    real(kind=rs_kind), intent(in), optional        :: missingvalue
    real(kind=rs_kind), intent(in), optional        :: coverage

    ! Local variables
    real(kind=rs_kind), dimension(size(data,1),size(data,2)) :: output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,delta,output,missingvalue,coverage)
       else
          call translate_image(data,delta,output,missingvalue)
       end if
    else
       call translate_image(data,delta,output)
    end if

  end function translate_function_real_wrap

  function translate_function_integer_wrap(data, delta, missingvalue, coverage) result(output)

    ! Translate in-place a double precision  precision matrix

    ! Interface variables
    integer, dimension(:,:), intent(in)          :: data
    real(kind=rd_kind), dimension(2), intent(in) :: delta
    integer, intent(in), optional                :: missingvalue
    real(kind=rs_kind), intent(in), optional     :: coverage

    ! Local variables
    integer, dimension(size(data,1),size(data,2)) :: output

    if (present(missingvalue)) then
       if (present(coverage)) then
          call translate_image(data,delta,output,missingvalue,coverage)
       else
          call translate_image(data,delta,output,missingvalue)
       end if
    else
       call translate_image(data,delta,output)
    end if

  end function translate_function_integer_wrap

end module image_transforms
