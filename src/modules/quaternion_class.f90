module quaternion_class

  use rotmatrix_class
  use vector_class
  use fundamental_constants
  use file_functions, only: stderr
  use precision

  implicit none

  private

  character(len=*), parameter :: version = "$Id: quaternion_class.f90,v 1.2 2004/08/09 02:46:49 aidan Exp aidan $"

  !! This module provides a 'class' for quaternions. Quaternions form a group
  !! that is isomorphic with the group SO(3), which is all (3x3) orthogonal
  !! matrices with determinant +1. This means quaternions can also represent
  !! proper rotations (see http://mathworld.wolfram.com/Quaternion.html for
  !! more details on quaternions).
  
  !! Because quaternions can only represent proper rotations and we wish to
  !! use quaternions in real crystallographic systems which also have improper
  !! (or roto-inversion) symmetry, we define a quaternion type that stores a
  !! flag to indicate if the quaternion in question also represents an improper
  !! rotation. In general this is just added functionality, but it does mean we
  !! cannot add two quaternions which differ in their improper values (note that
  !! multiplication is fine).
  
  !! This module provides a number of convenient features for making new
  !! quaternions, converting them to and from rotation matrices, a complete
  !! algebra for adding and subtracting them, as well as equalities for
  !! comparing quaternions.

  !! $Log: quaternion_class.f90,v $
  !! Revision 1.2  2004/08/09 02:46:49  aidan
  !! Added some options to two of the as_quaternion routines, to enable the unit
  !! and improper flags to specified at the time of creating the quaternion. This
  !! means the user doesn't have to use the improper function once the quaternion
  !! is created.
  !!
  !! Altered one of the routine names so that it fitted under the 31 character
  !! limit.
  !!
  !! Revision 1.1  2003/08/15 03:26:04  aidan
  !! Initial revision
  !!

  type quaternion
     private
     ! q(4) = (w, x, y, z) where w is the scalar part
     real :: q(4)
     ! Flags to tell us if we want to keep the quaternion unitised, and
     ! also if the quaternion represents an improper rotation (rotation + inversion)
     logical :: unit
     logical :: improper
  end type quaternion

  ! Define a tolerance for equivalence and "effectively zero"
  real :: tolerance = 1e-6 

  ! Define quaternion conjugate operator -- only monadic
  interface operator (.conjugate.)
     module procedure conjq
  end interface

  ! Define quaternion norm operator -- only monadic
  interface operator (.norm.)
     module procedure normq
  end interface

  ! Define quaternion inverse operator -- only monadic
  interface operator (.inverse.)
     module procedure invq
  end interface

  ! Define quaternion rotation routine (on a vector, or matrix of vectors)
  interface rotate
     module procedure qrotv
  end interface

  ! Define quaternion multiplication
  interface operator (*)
     module procedure qmultq, smultq, qmults, imultq, qmulti
  end interface

  ! Define quaternion division (only division by a scalar and 
  ! vice versa) -- quaternions in general don't commute, and division
  ! dosen't tell us what the order of multiplication of the inverses
  ! should be.
  interface operator (/)
     module procedure qdivs, sdivq
  end interface

  ! Define quaternion addition
  interface operator (+)
     module procedure qplusq, iplusq, qplusi, rplusq, qplusr, plusq, qplusarray, arrayplusq
  end interface

  ! Define quaternion subtraction
  interface operator (-)
     module procedure qminusq, iminusq, qminusi, rminusq, qminusr, minusq
  end interface

  ! Define quaternion assignment
  interface assignment (=)
     module procedure assign_array, assign_rotmatrix, assign_rmat_from_quaternion !, assign_xyz
  end interface

  ! Define quaternion equivalence operator (also defines the .eq. operator by default)
  interface operator (==)
     module procedure qeqq
  end interface

  ! Define quaternion non equivalence operator (also defines the .ne. operator by default)
  interface operator (/=)
     module procedure qneq
  end interface

  ! Define a quaternion constructor
  interface new
     module procedure new_quaternion
  end interface

  ! Define procedures to quaternion
  interface as_quaternion
     module procedure coerce_array_to_quaternion, coerce_rotmatrix_to_quaternion, &
                      coerce_axis_angle_to_quaternion !, coerce_xyz_to_quaternion
  end interface

  ! Define a procedure to coerce quaternion to array
  interface as_array
     module procedure as_array_quaternion
  end interface

  ! Define a procedure to coerce quaternion to rotation matrix
  interface as_rotmatrix
     module procedure coerce_quaternion_to_rotmatrix
  end interface

  ! Define a procedure to set and/or return the status of the "improper"
  ! flag, i.e. if the rotation being represented is an improper one
  interface improper
     module procedure improper_quaternion
  end interface

  interface normalise
     module procedure normalise_quaternion
  end interface

  interface axis_angle
     module procedure recover_axis_angle_from_quat
  end interface

  !!!!!!!!!!!!!!!!
  ! Public stuff !
  !!!!!!!!!!!!!!!!

  ! data types
  public :: quaternion

  ! routines
  public :: new, as_array, as_quaternion, as_rotmatrix, rotate, improper, axis_angle

  ! new operators
  public :: operator(.norm.), operator(.conjugate.), operator(.inverse.)

  ! overloaded operators
  public :: operator(/), operator(*), operator(+), operator(-), assignment(=)
  public :: operator(==), operator(/=)

contains

  pure subroutine new_quaternion(q, array, unit, improper)

    ! Initialise a new quaternion

    type (quaternion), intent(out) :: q
    real, optional, intent(in)     :: array(4)
    logical, optional, intent(in)  :: unit, improper

    ! Check if we have passed a unit flag -- if not we will assume
    ! that we want the quaternion to remain a 'unit quaternion' under
    ! all operations -- implies automatic normalising
    if ( present(unit) ) then
       q%unit = unit
    else
       q%unit = .TRUE.
    end if

    ! Check if we have passed an improper flag -- if not we will assume
    ! that the quaternion represents a proper rotation
    if ( present(improper) ) then
       q%improper = improper
    else
       q%improper = .FALSE.
    end if

    if ( present(array) ) then
       ! ******** deleted feature
       ! If the first element of the quaternion is negative then multiply
       ! the whole damn lot by -1 -- this is an exactly equivalent rotation,
       ! but it enforces a standard output from this routine
       ! if (array(1) < 0) array = -1.*array
       ! ******** deleted feature
       ! No longer do this sign change -- led to problems in qrotv (quaternion
       ! rotation) as the 'point' we were rotating was somewhat randomly changed
       ! sign depending on a the sign of the nonsense first value. We could
       ! reinstate this if we have a dedicated routine for quaternion rotation.
       ! Otherwise we will have to make equivalence insensitive to sign ..

       ! Don't use default assignment to initialise q
       q%q = array
    else
       ! Initialise quaternion to the equivalent of a 0 deg rotation
       ! if not passed an array.
       q%q = (/ 1, 0, 0, 0 /)
    end if

    ! Normalise quaternion if we have set this flag
    if (q%unit) call normalise_quaternion(q)

  end subroutine new_quaternion

  logical function improper_quaternion(q,switch)

    ! Improper rotation  switch

    type (quaternion), intent(inout) :: q
    logical, optional, intent(in)    :: switch
    if ( present(switch) ) then
       q%improper = switch
    end if
    improper_quaternion = q%improper

  end function improper_quaternion

  elemental subroutine normalise_quaternion(q)

    ! Normalise the quaternion by dividing the quaternion by it's
    ! 'reduced norm'

    type (quaternion), intent(inout) :: q

    q%q = q%q / (sqrt(sum(q%q**2)))
    ! q%q = q%q / (.norm. q)

  end subroutine normalise_quaternion

  elemental function normq(q)

    ! The 'reduced norm' of a quaternion is just the square root of the
    ! arithmetic sum of the squares of it's elements

    real :: normq
    type (quaternion), intent(in) :: q

    normq = sqrt(sum(q%q**2))

  end function normq
  
  function as_array_quaternion(q)

    ! Return quaternion as a real array -- useful for printing
    ! the order of elements is w, x, y, z

    real :: as_array_quaternion(4)
    type (quaternion) :: q

    as_array_quaternion = q%q

  end function as_array_quaternion
  
  function coerce_array_to_quaternion(array, unit, improper)

    ! Coerce an array with 4 elements to a quaternion

    ! Input variables
    real, dimension(:), intent(in) :: array
    logical, intent(in), optional  :: unit, improper

    ! Return type of function
    type (quaternion) :: coerce_array_to_quaternion
    
    ! Local variables
    integer :: dimlength
    logical :: improper_flag, unit_flag

    improper_flag = .FALSE.
    if (present(improper)) improper_flag = improper
    unit_flag = .TRUE.
    if (present(unit)) unit_flag = unit

    dimlength = size(array)

    if (dimlength == 4) then
       ! Call quaternion constructor
       call new(coerce_array_to_quaternion, array, unit_flag, improper_flag)
    else if (dimlength == 3) then
       call new(coerce_array_to_quaternion, (/ 0.0, array /), unit_flag, improper_flag)
    else
       stop "Attempted to initialise quaternion with invalid number of elements: "
    end if

  end function coerce_array_to_quaternion
  
  function coerce_xyz_to_quaternion(xyz)

    ! Coerce a array of 3 elements (assume it is a point in cartesian
    ! space) to a quaternion

    real, dimension(3), intent(in) :: xyz
    type (quaternion) :: coerce_xyz_to_quaternion

    ! Call quaternion constructor
    call new(coerce_xyz_to_quaternion, (/ 0.0, xyz /))

  end function coerce_xyz_to_quaternion
  
  function coerce_axis_angle_to_quaternion(array, angle, improper)

    ! Coerce an axis plus an angle (in radians) to a quaternion

    ! Input variables
    real, intent(in)              :: array(3), angle
    logical, intent(in), optional :: improper

    ! Return type of function
    type (quaternion) :: coerce_axis_angle_to_quaternion

    ! Local variables
    real    :: halfangle
    logical :: improper_flag

    improper_flag = .FALSE.
    if (present(improper)) improper_flag = improper

    halfangle = angle/2.

    ! Call the quaternion constructor with the appropriate arguments
    call new(coerce_axis_angle_to_quaternion, (/ cos(halfangle), sin(halfangle)*array /), improper=improper_flag)

  end function coerce_axis_angle_to_quaternion
  
  subroutine recover_axis_angle_from_quat(q, axis, angle)

    ! Recover an axis plus an angle (in radians) from a quaternion

    type (quaternion), intent(in) :: q
    real, intent(out) :: axis(3), angle

    if ( sum(q%q(2:4)) == 0. ) then
       !  If we do not have a rotation axis choose the z-axis (arbitrary choice)
       angle = 0.
       axis  = (/ 0., 0., 1. /)
    else if ( q%q(1) == 0.0 ) then
       !  If we do not have a rotation angle choose the 180 degrees (pi)
       angle   = pi
       axis = q%q(2:4)
    else
       ! Ensure this was a unit vector
       axis = q%q(2:4)
       call normalise(axis)
       angle = 2.0 * atan2 ( (.norm. q%q(2:4)), real(q%q(1),rd_kind) )
    end if

  end subroutine recover_axis_angle_from_quat

  function coerce_quaternion_to_rotmatrix(q)

    ! Make a 3x3 rotation matrix from a quaternion -- stolen 
    ! from the JPL 'spicelib' routine q2m

    type (quaternion) :: q
    type (rotmatrix) :: coerce_quaternion_to_rotmatrix
    ! real, dimension(3,3) :: coerce_quaternion_to_rotmatrix
    real(kind=rd_kind) :: R(3,3)
    real(kind=rd_kind) :: Q01, Q02, Q03, Q12, Q13, Q23, Q1S, Q2S, Q3S
    real(kind=rd_kind) :: L2, SHARPN
 
    Q01  =  q%q(1) * q%q(2)
    Q02  =  q%q(1) * q%q(3)
    Q03  =  q%q(1) * q%q(4)
    
    Q12  =  q%q(2) * q%q(3)
    Q13  =  q%q(2) * q%q(4)
    
    Q23  =  q%q(3) * q%q(4)
 
    Q1S  =  q%q(2) * q%q(2)
    Q2S  =  q%q(3) * q%q(3)
    Q3S  =  q%q(4) * q%q(4)

    !  We sharpen the computation by effectively converting Q to
    !  a unit quaternion if it isn't one already.

    L2   =  q%q(1) * q%q(1) + Q1S + Q2S + Q3S
 
    
    IF ( L2 .NE. 1.0D0 .AND. L2 .NE. 0.0D0 ) THEN
 
       SHARPN = 1.0D0 / L2
 
       Q01    = Q01 * SHARPN
       Q02    = Q02 * SHARPN
       Q03    = Q03 * SHARPN

       Q12    = Q12 * SHARPN
       Q13    = Q13 * SHARPN
 
       Q23    = Q23 * SHARPN
 
       Q1S    = Q1S * SHARPN
       Q2S    = Q2S * SHARPN
       Q3S    = Q3S * SHARPN

    END IF

    R(1,1) =  1.D0  -  2.D0 * ( Q2S + Q3S )
    R(2,1) =           2.D0 * ( Q12 + Q03 )
    R(3,1) =           2.D0 * ( Q13 - Q02 )
    
    R(1,2) =           2.D0 * ( Q12 - Q03 )
    R(2,2) =  1.D0  -  2.D0 * ( Q1S + Q3S )
    R(3,2) =           2.D0 * ( Q23 + Q01 )
 
    R(1,3) =           2.D0 * ( Q13 + Q02 )
    R(2,3) =           2.D0 * ( Q23 - Q01 )
    R(3,3) =  1.D0  -  2.D0 * ( Q1S + Q2S )
    
    ! If the quaternion is improper it represents an improper rotation
    ! i.e. a rotoinversion, so we have to invert the rotation matrix
    if (q%improper) R = -1.d0 * R

    ! Return a single precision result ... also use default assignment
    ! (as defined in the rotmatrix_class module
    coerce_quaternion_to_rotmatrix = real(R)

  end function coerce_quaternion_to_rotmatrix

  function coerce_rotmatrix_to_quaternion(R)

    ! Make a quaternion from a rotation matrix

    type (quaternion) :: coerce_rotmatrix_to_quaternion
    type (rotmatrix), intent(in) :: R
    real :: matrix(3,3)
    real(kind=rd_kind) :: c, cc4, factor, mtrace, s(3), s114, s224, s334, trace
    logical :: improper

    ! Extract 3x3 rotation matrix from rotmatrix
    matrix = as_matrix(R)

    ! Check if the rotation matrix is improper
    if ((.det. R) < 0.0 ) then
       matrix = -1 * matrix
       improper = .true.
    else
       improper = .false.
    end if

    !     Q              is a unit quaternion representing R.  Q is a
    !                    4-dimensional vector.  If R rotates vectors by an
    !                    angle of r radians about a unit vector A, where
    !                    r is in [0, pi], then if h = r/2,
    !
    !                       Q = ( cos(h), sin(h)A ,  sin(h)A ,  sin(h)A ).
    !                                            1          2          3
    !
    !                    The restriction that r must be in the range [0, pi]
    !                    determines the output quaternion Q uniquely
    !                    except when r = pi; in this special case, both of
    !                    the quaternions
    !
    !                       Q = ( 0,  A ,  A ,  A  )
    !                                  1    2    3
    !                    and
    !
    !                       Q = ( 0, -A , -A , -A  )
    !                                  1    2    3
    !
    !                   are possible outputs, if A is a choice of rotation
    !                   axis for R.

    !     If our quaternion is C, S1, S2, S3 (the S's being the imaginary
    !     part) and we let
    !
    !        CSi = C  * Si
    !        Sij = Si * Sj
    !
    !     then the rotation matrix corresponding to our quaternion is:
    !
    !        R(1,1)      = 1.0D0 - 2*S22 - 2*S33
    !        R(2,1)      =         2*S12 + 2*CS3
    !        R(3,1)      =         2*S13 - 2*CS2
    !
    !        R(1,2)      =         2*S12 - 2*CS3
    !        R(2,2)      = 1.0D0 - 2*S11 - 2*S33
    !        R(3,2)      =         2*S23 + 2*CS1
    !
    !        R(1,3)      =         2*S13 + 2*CS2
    !        R(2,3)      =         2*S23 - 2*CS1
    !        R(3,3)      = 1.0D0 - 2*S11 - 2*S22
    !
    !        From the above we can see that
    !
    !           TRACE = 3 - 4*(S11 + S22 + S33)
    !
    !        so that
    !
    !
    !           1.0D0 + TRACE = 4 - 4*(S11 + S22 + S33)
    !                         = 4*(CC + S11 + S22 + S33)
    !                         - 4*(S11 + S22 + S33)
    !                         = 4*CC
    !
    !        Thus up to sign
    !
    !          C = 0.5D0 * DSQRT( 1.0D0 + TRACE )
    !
    !        But we also have
    !
    !          1.0D0 + TRACE - 2.0D0*R(i,i) = 4.0D0 - 4.0D0(Sii + Sjj + Skk)
    !                                       - 2.0D0 + 4.0D0(Sjj + Skk )
    !
    !                                       = 2.0D0 - 4.0D0*Sii
    !
    !        So that
    !
    !           1.0D0 - TRACE + 2.0D0*R(i,i) = 4.0D0*Sii
    !
    !        and so up to sign
    !
    !           Si = 0.5D0*DSQRT( 1.0D0 - TRACE + 2.0D0*R(i,i) )
    !
    !        in addition to this observation, we note that all of the
    !        product pairs can easily be computed
    !
    !         CS1 = (R(3,2) - R(2,3))/4.0D0
    !         CS2 = (R(1,3) - R(3,1))/4.0D0
    !         CS3 = (R(2,1) - R(1,2))/4.0D0
    !         S12 = (R(2,1) + R(1,2))/4.0D0
    !         S13 = (R(3,1) + R(1,3))/4.0D0
    !         S23 = (R(2,3) + R(3,2))/4.0D0
    !
    !     But taking sums or differences of numbers that are nearly equal
    !     or nearly opposite results in a loss of precision. As a result
    !     we should take some care in which terms to select when computing
    !     C, S1, S2, S3.  However, by simply starting with one of the
    !     large quantities cc, S11, S22, or S33 we can make sure that we
    !     use the best of the 6 quantities above when computing the
    !     remaining components of the quaternion.
    !
    
    trace  = matrix(1,1) + matrix(2,2) + matrix(3,3)
    mtrace = 1.0d0 - trace
 
    cc4    = 1.0d0  + trace
    s114   = mtrace + 2.0d0*matrix(1,1)
    s224   = mtrace + 2.0d0*matrix(2,2)
    s334   = mtrace + 2.0d0*matrix(3,3)
 
    !
    !     Note that if you simply add CC4 + S114 + S224 + S334
    !     you get four. Thus at least one of the 4 terms is greater than 1.
    !
    if ( cc4 >= 1.0d0 ) then
       c      =  sqrt  ( cc4 * 0.25d0 )
       factor =  1.0d0 /( c   * 4.0d0  )
       
       s(1)   = ( matrix(3,2) - matrix(2,3) )*factor
       s(2)   = ( matrix(1,3) - matrix(3,1) )*factor
       s(3)   = ( matrix(2,1) - matrix(1,2) )*factor
    else if ( s114 >= 1.0d0 ) then
       s(1)   = sqrt  ( s114 * 0.25d0 )
       factor = 1.0d0 /( s(1) * 4.0d0  )
       
       c      = ( matrix(3,2) - matrix(2,3) ) * factor
       s(2)   = ( matrix(1,2) + matrix(2,1) ) * factor
       s(3)   = ( matrix(1,3) + matrix(3,1) ) * factor
    else if ( s224 >= 1.0 ) then
       s(2)   = sqrt  ( s224 * 0.25d0 )
       factor = 1.0d0 /( s(2) * 4.0d0  )
       
       C      = ( matrix(1,3) - matrix(3,1) ) * factor
       S(1)   = ( matrix(1,2) + matrix(2,1) ) * factor
       S(3)   = ( matrix(2,3) + matrix(3,2) ) * factor
    else
       s(3)   = sqrt  ( s334 * 0.25d0 )
       factor = 1.0d0 /( s(3) * 4.0d0  )
       
       c      = ( matrix(2,1) - matrix(1,2) ) * factor
       s(1)   = ( matrix(1,3) + matrix(3,1) ) * factor
       s(2)   = ( matrix(2,3) + matrix(3,2) ) * factor
    end if

    ! Initialise new quaternion
    call new(coerce_rotmatrix_to_quaternion,array=real((/ c, s /)),improper=improper)

  end function coerce_rotmatrix_to_quaternion

  !
  ! ASSIGNMENT
  !

  subroutine assign_array(q,array)

    ! Assign to quaternion the elements in array

    type (quaternion), intent(out)  :: q
    real, intent(in) :: array(4)

    ! Use the quaternion constructor to assign values to q
    call new(q,array=array)

  end subroutine assign_array

  subroutine assign_rotmatrix(q,R)

    ! Allow simple assignment from a rotation matrix to a quaternion

    type (quaternion), intent(out)  :: q
    type (rotmatrix), intent(in)    :: R

    ! Call the coercion function ...
    q =  as_quaternion(R)

  end subroutine assign_rotmatrix

  subroutine assign_rmat_from_quaternion(R,q)

    ! Allow simple assignment from quaternion to a rotation matrix

    type (rotmatrix), intent(out)  :: R
    type (quaternion), intent(in)  :: q

    ! Call the coercion function ...
    R =  as_rotmatrix(q)

  end subroutine assign_rmat_from_quaternion

  subroutine assign_xyz(q,xyz)

    ! A point in cartesian space (x,y,z) can be represented 
    ! as a quaternion -- the x, y and z values become the 
    ! 'vector' part, with the scalar part set to zero. 
    ! 
    ! This subroutine allows simple assignment of a point
    ! to a quaternion -- so an array of size 3 is assumed to
    ! be a point in cartesian space.
    !
    ! Currently not used -- might 'trap' someone who inadvertently
    ! puts only three elements in their assigned array. I think it
    ! is better to force the programmer to coerce the xyz values 
    ! to a quaternion use coerce_xyz_to_quaternion.

    type (quaternion), intent(out)  :: q
    real, intent(in) :: xyz(3)

    ! Use default assignment to initialise q
    q = (/ 0.0, xyz /)
    
  end subroutine assign_xyz

  !
  ! CONJUGATION
  !

  elemental function conjq(q)

    ! Return conjugate of quaternion 

    type (quaternion)  :: conjq
    type (quaternion), intent(in) :: q

    conjq = q

    ! The conjugate of q is defined to be the negative of the vector
    ! part of the quaternion, i.e.
    conjq%q(2:4) = conjq%q(2:4)*(-1)

  end function conjq

  !
  ! INVERSION
  !

  elemental function invq(q)

    ! Return inverse of a quaternion 

    type (quaternion)  :: invq
    type (quaternion), intent(in) :: q

    ! The inverse of q is defined to be the conjugate divided by the square of the
    ! norm. Had to explicitly call smultq -- complained about the multiplication 
    ! being ambiguous ..
    invq = smultq((1/(.norm. q)**2),(.conjugate. q))

  end function invq


  !
  ! DIVISION
  !

  elemental function qdivs(q, s)

    ! Return quaternion divided by a scalar -- this is just
    ! translated to multiplying by the inverse of the scalar

    type (quaternion)  :: qdivs
    type (quaternion), intent(in) :: q
    real, intent(in) :: s
    real :: news

    ! if (s == 0) stop "Attempted to divide quaternion by zero!"

    news = s**(-1)

    qdivs = smultq(news, q)

  end function qdivs

  elemental function sdivq(s, q)

    ! Return scalar divided by a quaternion -- this is just
    ! translated to the inverse of the quaternion multiplied 
    ! by the scalar

    type (quaternion)  :: sdivq
    real, intent(in) :: s
    type (quaternion), intent(in) :: q

    sdivq = .inverse. q

    ! Had to explicitly call smultq -- complained about the multiplication 
    ! being ambiguous ..
    sdivq = smultq(s,sdivq)

  end function sdivq

  !
  ! MULTIPLICATION
  !

  elemental function qmultq(qL,qR)

    ! Return quaternion product qL * qR

    type (quaternion)  :: qmultq
    type (quaternion), intent(in) :: qL, qR
    real(kind=rd_kind) :: Rmat(4,4), dbl_array(4)
    real :: w, v(3)
    ! w = w1w2-v1.v2, v = v1 X v2 + w1v2 + w2v1

    ! Currently not correct ...
    !
    ! w = qL%q(1)*qR%q(1) - qL%q(2:4) .dot. qR%q(2:4) 
    ! v = qL%q(2:4) .cross. qR%q(2:4)
    ! v = v + qL%q(1)*qR%q(2:4)
    ! v = v + qR%q(1)*qL%q(2:4) 

    ! call new(qmultq, array=(/ w, v /), unit=(qL%unit .and. qR%unit), &
    !      improper=((qL%improper .neqv. qR%improper)))
    
    ! Make a 4x4 matrix out of the coefficients of the
    ! RHS quaternion ... works ok
    ! Rmat = reshape( &
    !      (/  qR%q(1),  qR%q(2),  qR%q(3),  qR%q(4),  &
    !      -qR%q(2),  qR%q(1), -qR%q(4),  qR%q(3),  &
    !      -qR%q(3),  qR%q(4),  qR%q(1), -qR%q(2),  &
    !      -qR%q(4), -qR%q(3),  qR%q(2),  qR%q(1) /), &
    !      (/4,4/), order=(/2,1/) )

    ! dbl_array = matmul( qL%q, Rmat )

    dbl_array(1) = qL%q(1)*qR%q(1) - qL%q(2)*qR%q(2) - qL%q(3)*qR%q(3) - qL%q(4)*qR%q(4) 
    dbl_array(2) = qL%q(1)*qR%q(2) + qL%q(2)*qR%q(1) + qL%q(3)*qR%q(4) - qL%q(4)*qR%q(3) 
    dbl_array(3) = qL%q(1)*qR%q(3) - qL%q(2)*qR%q(4) + qL%q(3)*qR%q(1) + qL%q(4)*qR%q(2) 
    dbl_array(4) = qL%q(1)*qR%q(4) + qL%q(2)*qR%q(3) - qL%q(3)*qR%q(2) + qL%q(4)*qR%q(1) 

    ! Don't use default assignment for result of matrix multiplication.
    ! Call new directly, so we can determine if the resulting quaternion
    ! should be normalised -- only normalise if *both* quaternions have the
    ! unit flag set. Also perform exclusive or (XOR) operation on the
    ! improper flags, i.e. resulting quaternion is only improper if either
    ! qL or qR is improper, but not both.
    call new(qmultq, array=real(dbl_array), unit=(qL%unit .and. qR%unit), &
         improper=((qL%improper .neqv. qR%improper)))

    ! qmultq%q = matmul( qL%q, Rmat )
    ! qmultq = qmultq / sqrt(sum(qmultq**2))

    ! Or ...
    ! 
    ! w = w1w2-v1.v2, v = v1 X v2 + w1v2 + w2v1

    ! Or ...
    !
    ! P*Q = <(P:h*Q:h - P:i*Q:i - P:j*Q:j - P:k*Q:k), 
    !         P:h*Q:i + P:i*Q:h + P:j*Q:k - P:k*Q:j), 
    !         P:h*Q:j - P:i*Q:k + P:j*Q:h + P:k*Q:i), 
    !         P:h*Q:k + P:i*Q:j - P:j*Q:i + P:k*Q:h)>
    !
    ! Or ..
    !
    ! w = w1w2 - x1x2 - y1y2 - z1z2
    ! x = w1x2 + x1w2 + y1z2 - z1y2
    ! y = w1y2 + y1w2 + z1x2 - x1z2
    ! z = w1z2 + z1w2 + x1y2 - y1x2

  end function qmultq

  elemental function smultq(s, q)

    ! Return quaternion multiplied by a scalar

    type (quaternion)  :: smultq
    real, intent(in) :: s
    type (quaternion), intent(in) :: q

    ! We assign the quaternion directly (don't use default assignment)
    ! as we don't want to normalise the resulting quaternion -- why 
    ! would we!? We have just multiplied it by a scalar!
    ! smultq%q = s * q%q
    call new(smultq, array=q%q * s, unit=q%unit, improper=q%improper)

  end function smultq

  elemental function qmults(q, s)

    ! Return quaternion multiplied by a scalar
    ! (same as above, but swapped order of quat & scal)

    type (quaternion)  :: qmults
    real, intent(in) :: s
    type (quaternion), intent(in) :: q

    qmults = smultq(s, q)
    ! qmults%q = q%q * s

  end function qmults

  elemental function qmulti(q, i)

    ! Return quaternion multiplied by an integer

    type (quaternion)  :: qmulti
    integer, intent(in) :: i
    type (quaternion), intent(in) :: q

    qmulti = smultq(real(i),q)

  end function qmulti

  elemental function imultq(i, q)

    ! Return quaternion multiplied by an integer

    type (quaternion)  :: imultq
    integer, intent(in) :: i
    type (quaternion), intent(in) :: q

    imultq = smultq(real(i),q)

  end function imultq

  !
  ! ADDITION
  !

  function qplusq(qL, qR)

    ! Return quaternion added to quaternion. This is not an elemental function
    ! because we need to check the status of the improper flags of the two
    ! quaternions -- we don't know how to handle the situation where they are
    ! not the same

    type (quaternion)  :: qplusq
    type (quaternion), intent(in) :: qL, qR

    ! We stop catastrophically if either quaternion is improper and 
    ! the other is not -- we don't know how to deal with that situation!
    if (qL%improper .neqv. qR%improper) then
       print *,'One of the quaternions you are attempting to add together represents ',&
       'an improper rotation, but the other does not. We cannot cater for this situation!'
       stop
    end if

    ! Don't use default assignment. The result of the addition of
    ! the two arrays will be an array, not a quaternion, which we
    ! pass to the quaternion constructor 'new'. We normalise the
    ! resulting quaternion if either quaternions have the 'unit' flag
    ! set. We make the resulting quaternion improper if both of the 
    ! quaternions are improper, and vice versa (need only use the 
    ! improper flag from one of the quaternions, as we know if we 
    ! have progressed this far that both will have the same value). 
    call new(qplusq, array=(qL%q + qR%q), unit=(qL%unit.or.qR%unit), improper=qL%improper)

  end function qplusq

  elemental function qplusr(q, r)

    ! Return quaternion added to real

    ! Return type of function
    type (quaternion)  :: qplusr

    ! Input variables
    type (quaternion), intent(in) :: q
    real, intent(in) :: r

    ! Add the real value to the quaternion as an array
    call new(qplusr, array=(q%q + (/ r, 0.0, 0.0, 0.0 /)), unit=q%unit, improper=q%improper)

  end function qplusr

  elemental function rplusq(r, q)

    ! Return real added to quaternion 

    type (quaternion)  :: rplusq
    type (quaternion), intent(in) :: q
    real, intent(in) :: r

    ! Add quaternion to real -- will call qplusr above
    rplusq = q + r

  end function rplusq

  elemental function qplusi(q, i)

    ! Return quaternion added to integer

    type (quaternion)  :: qplusi
    type (quaternion), intent(in) :: q
    integer, intent(in) :: i

    ! Add quaternion to integer (coerced to real -- will call qplusr above)
    qplusi = q + real(i)

  end function qplusi

  elemental function iplusq(i, q)

    ! Return integer added to quaternion

    type (quaternion)  :: iplusq
    type (quaternion), intent(in) :: q
    integer, intent(in) :: i

    ! Add quaternion to integer (coerced to real -- will call qplusr above)
    iplusq = q + real(i)

  end function iplusq

  elemental function plusq(q)

    ! Return unary addition of quaternion

    type (quaternion)  :: plusq
    type (quaternion), intent(in) :: q

    plusq = q

  end function plusq

  function qplusarray(q, array)

    ! Return quaternion added to array of 4 real numbers

    ! Return type
    type (quaternion)  :: qplusarray

    ! Input variables
    type (quaternion), intent(in) :: q
    real, dimension(:), intent(in) :: array

    ! Local variables
    real, dimension(4) :: dummy_array
    character(len=40) :: sizestring

    select case(size(array))
       case (3)
          dummy_array = (/ 0., array /)
       case (4)
          dummy_array = array
       case default
          write(sizestring,*) size(array)
          write(stderr,*)"QUATERNION_CLASS :: Do not know how add an array of length ",trim(adjustl(sizestring))," to a quaternion!"
          stop
    end select

    ! We pass to the quaternion constructor 'new'. We use the flags set
    ! in the quaternion to determine if we normalise. We make the 
    ! resulting quaternion improper if the quaternion is improper.
    call new(qplusarray, array=(q%q + dummy_array), unit=q%unit, improper=q%improper)

  end function qplusarray

  function arrayplusq(array, q)

    ! Return array of 4 real numbers added to quaternion

    type (quaternion)  :: arrayplusq
    type (quaternion), intent(in) :: q

    real, dimension(:), intent(in) :: array

    ! Swap the order of the arguments and pass to qplusarray
    arrayplusq = qplusarray(q, array)
    
  end function arrayplusq

  !
  ! SUBTRACTION
  !

  function qminusq(qL, qR)

    ! Return quaternion subtracted from quaternion

    type (quaternion)  :: qminusq
    type (quaternion), intent(in) :: qL, qR

    ! Use default assignment, i.e. the result of the subtraction of
    ! the two arrays will be an array, not a quaternion, so the
    ! rule for assigning arrays will be invoked -- takes care of
    ! normalising the resulting quaternion
    qminusq = qL + (-qR)

  end function qminusq

  elemental function qminusr(q, r)

    ! Return real subtracted from quaternion

    type (quaternion)  :: qminusr
    type (quaternion), intent(in) :: q
    real, intent(in) :: r

    ! Add the quaternion to minus the real value -- will call
    ! qplusr above
    qminusr = q + (-r)

  end function qminusr

  elemental function rminusq(r, q)

    ! Return quaternion subtracted from real 

    type (quaternion)  :: rminusq
    type (quaternion), intent(in) :: q
    real, intent(in) :: r

    ! Add the real to minus the quaternion -- will call
    ! rplusq above
    rminusq = r + (-q)

  end function rminusq

  elemental function qminusi(q, i)

    ! Return quaternion subtracted to integer

    type (quaternion)  :: qminusi
    type (quaternion), intent(in) :: q
    integer, intent(in) :: i

    ! Will call rminusq above
    qminusi = q - real(i)

  end function qminusi

  elemental function iminusq(i, q)

    ! Return integer subtracted to quaternion

    type (quaternion)  :: iminusq
    type (quaternion), intent(in) :: q
    integer, intent(in) :: i

    ! Will call rminusq above
    iminusq = real(i) - q

  end function iminusq

  elemental function minusq(q)

    ! Return unary subtraction of quaternion

    type (quaternion)  :: minusq
    type (quaternion), intent(in) :: q

    ! Copy all the attributes of q to minusq
    minusq = q

    ! ... and then 'negate' the quaternion array
    minusq%q = -1 * q%q

  end function minusq

  function arrayminusq(array, q)

    ! Return array of 4 real numbers added to quaternion

    type (quaternion)  :: arrayminusq
    type (quaternion), intent(in) :: q

    real, dimension(3), intent(in) :: array

    ! Swap the order of the arguments and pass to qplusarray
    arrayminusq = qplusarray(-q, array)
    
  end function arrayminusq

  function qminusarray(array, q)

    ! Return array of 4 real numbers added to quaternion

    type (quaternion)  :: qminusarray
    type (quaternion), intent(in) :: q

    real, dimension(3), intent(in) :: array

    ! Swap the order of the arguments and pass to qplusarray
    qminusarray = qplusarray(q, -array)
    
  end function qminusarray

  !
  ! EQUIVALENCE
  !

  function qeqq(qL, qR)

    ! Determines if two quaternions are equivalent (within
    ! an arbitrary tolerance)

    logical :: qeqq
    type (quaternion), intent(in) :: qL, qR

    real :: positive, negative

    qeqq = .false.

    positive = sum(abs(qL%q + qR%q))
    negative = sum(abs(qL%q - qR%q))

    ! See if the difference betwixt the quaternions is greater than the tolerance
    if ((positive < tolerance) .or. (negative < tolerance)) qeqq = .true.

  end function qeqq

  function qneq(qL, qR)

    ! Determines if two quaternions are not equivalent (within
    ! an arbitrary tolerance)

    logical :: qneq
    type (quaternion), intent(in) :: qL, qR

    qneq = .true.

    ! Just call the equivalence routine and return the opposite value
    if (qL == qR) qneq = .false.

  end function qneq

  !
  ! ROTATION
  !

  subroutine qrotv(qR,v)

    ! Rotate a real space coordinate by a quaternion. 
    !
    ! To rotate a real space coordinate by a quaternion the coordinate
    ! is converted into it's quaternion representation and the rotated
    ! quaternionis is defined by
    !
    !     qR * q * qR^-1
    !
    ! where qR is the quaternion we are rotating by, and q is the 
    ! quaternion being rotated. Recall that a point (x,y,z) can be
    ! represented as a quaternion, i.e. (0,x,y,z).
    !
    ! After rotation real space coordinates are extracted out of the
    ! quaternion.
 
    real, dimension(3), intent(inout) :: v
    ! qR is "inout" even though we don't intend on altering it -- we call
    ! the improper function on qR later on, and this has the ability to alter
    ! qR (though it dosen't in this case)
    type (quaternion), intent(inout)  :: qR

    ! Temporary variables we will use to hold intermediate results
    type (quaternion)  :: q
    real, dimension(4) :: array

    ! Make a quaternion version of the vector v -- don't use default
    ! constructor so we can specify *not* to automatically make q a 
    ! unit quaternion
    ! q = as_quaternion(v)
    call new(q, array=(/ 0.0, v /), unit=.false.)
    
    ! Perform quaternion rotation
    q = (qR * q) * (.inverse. qR)

    ! Extract out the values of the quaternion as a 4-member array
    array = as_array(q)

    ! The resulting quaternion should have a zero in the first (scalar, w)
    ! position, i.e. it represents a position in real space.
    if (array(1) > tolerance) stop "Something is rotten in the state of Denmark!"

    ! If qR is an improper rotation then invert the coordinates
    if (improper(qR)) then
       array = -array
    end if

    ! The second, third and fourth elements of the quaternion are the 
    ! rotated real space coordinates
    v = array(2:4)

  end subroutine qrotv

  subroutine qrotq(qR,q)

    ! Rotate a quaternion by another quaternion. This is not used.
    ! Not sure it is valid for general quaternions (i.e. with w /= 0)
 
    type (quaternion), intent(in) :: qR
    type (quaternion), intent(inout) :: q

    q = qR * q * (.inverse. qR)

  end subroutine qrotq

end module quaternion_class
