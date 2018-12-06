module rotmatrix_class

  use file_functions, only: stderr
  use precision
  use vector_class

  implicit none

  private

  character(len=*), parameter :: version = "$Id: rotmatrix_class.f90,v 1.1 2004/08/04 04:35:48 aidan Exp aidan $"

  !! Implements rotation matrices. Provides methods for comparing, querying,
  !! creating rotation matrices.

  !! $Log: rotmatrix_class.f90,v $
  !! Revision 1.1  2004/08/04 04:35:48  aidan
  !! Initial revision
  !!

  type rotmatrix
     private
     real(kind=rd_kind) :: R(3,3)
     logical :: inverse
  end type rotmatrix

  ! Define rotmatrix determinant operator -- only monadic
  interface operator (.det.)
     module procedure detr, detm
  end interface

  ! Define operator for determining if something is a valid rotation matrix
  interface operator (.isrot.)
     module procedure isrotm, isrotrm
  end interface

  ! Define rotate routine (on a vector, or matrix of vectors)
  interface rotate
     module procedure rrotv, rrotm, rrotsinglem, rrotsinglev
  end interface

  ! Define rotmatrix multiplication
  interface operator (*)
     module procedure rmultr, smultr, rmults, imultr, rmulti
  end interface

  ! Define rotmatrix addition
  interface operator (+)
     module procedure rplusr, iplusr, rplusi, splusr, ssingleplusr, rpluss, rplussingles, plusr
  end interface

  ! Define rotmatrix subtraction
  interface operator (-)
     module procedure rminusr, iminusr, rminusi, sminusr, ssingleminusr, rminuss, rminussingles, minusr
  end interface

  ! Define rotmatrix assignment
  interface assignment (=)
     module procedure assign_matrix, assign_doublematrix !, assign_xyz
  end interface

  ! Define rotmatrix equivalence operator (also defines the .eq. operator by default)
  interface operator (==)
     module procedure reqr
  end interface

  ! Define rotmatrix inequivalence operator (also defines the .ne. operator by default)
  interface operator (/=)
     module procedure rner
  end interface

  ! Define a rotmatrix constructor
  interface new
     module procedure new_rotmatrix, new_rotmatrix_real_wrapper
  end interface

  ! Define an interface to coerce various types of data to a rotation matrix
  interface as_rotmatrix
     module procedure coerce_matrix_to_rotmatrix, &
          coerce_caxis_angle_to_rotmatrix, coerce_caxis_singang_to_rotmat !, coerce_aaxis_angle_to_rotmatrix
  end interface

  ! Define a function to coerce rotmatrix to matrix
  interface as_matrix
     module procedure as_matrix_rotmatrix
  end interface

  interface reorthogonalise
     module procedure reorthogonalise_matrix, reorthogonalise_rotmatrix
  end interface

  ! interface normalise
  !    module procedure normalise_rotmatrix
  ! end interface

  !!!!!!!!!!!!!!!!
  ! Public stuff !
  !!!!!!!!!!!!!!!!

  ! data types
  public :: rotmatrix

  ! routines
  public :: new, as_matrix, as_rotmatrix, rotate, reorthogonalise

  ! new operators
  public :: operator(.det.), operator(.isrot.)

  ! overloaded operators
  public :: operator(*), operator(+), operator(-), assignment(=)
  public :: operator(==), operator(/=)

contains

  subroutine new_rotmatrix(R, matrix, inverse)

    ! Initialise a new rotmatrix. This is _the_ definitive routine for
    ! assigning rotation matrices. All other routines call this one in
    ! some way.

    type (rotmatrix), intent(out) :: R
    real(kind=rd_kind), optional  :: matrix(3,3)
    logical, optional             :: inverse

    integer :: i
    logical :: dummy_inverse

    ! Don't use default assignment to initialise R -- this would be
    ! circular as the default assignment calls this routine.
    if ( present(matrix) ) then
       ! Check we have a rotation matrix
       if (.not. .isrot. matrix) then
          call reorthogonalise(matrix)
          if (.not. .isrot. matrix) then
             write(stderr,*) "Not a rotation matrix: "
             do i=1,3
                write(*,'(4X,3F10.4)') matrix(i,1:3)
             end do
             stop
          end if
       end if
       R%R = matrix
    else
       ! Initialise rotmatrix to the identity matrix if not passed a matrix
       R%R = reshape( (/ 1., 0., 0., 0., 1., 0., 0., 0., 1. /), (/3, 3/) )
    end if

    ! Initialise inverse to false if not passed a flag for this
    if ( present(inverse) ) then
       R%inverse = inverse
    else
       R%inverse = .false.
    end if

  end subroutine new_rotmatrix

  subroutine new_rotmatrix_real_wrapper(R, matrix, inverse)

    ! Initialise a new rotmatrix with a single precision input

    ! Interface variables
    type (rotmatrix), intent(out) :: R
    real(kind=rs_kind)            :: matrix(3,3)
    logical, optional             :: inverse

    if ( present(inverse) ) then
       call new(R, real(matrix,rd_kind), inverse)
    else
       call new(R, real(matrix,rd_kind))
    end if

  end subroutine new_rotmatrix_real_wrapper

  subroutine gram_schmidt_orthogonalisation (a, r)

    real(kind=rd_kind), dimension(:,:), intent(inout) :: a, r
    real(kind=rd_kind)  :: sum
    integer           :: n, m
    integer           :: i, j, k

    do k = 1, m
       sum = 0
       do i = 1, n
          sum = sum + a(i, k)*a(i, k)
       end do
       r(k, k) = sqrt(sum)
       do i = 1, n
          a(i, k) = a(i, k) / r(k, k)
       end do
       do j = k + 1, m
          sum = 0
          do i = 1, n
             sum = sum + a(i, k) * a(i, j)
          end do
          r(k, j) = sum
          do i = 1, n
             a(i, j) = a(i, j) - a(i, k) * r(k, j)
          end do
       end do
    end do
    
  end subroutine gram_schmidt_orthogonalisation

  ! subroutine normalise_rotmatrix(q)

    ! Normalise the rotmatrix by dividing the rotmatrix by it's
    ! 'reduced norm'

    ! type (rotmatrix) :: q

    ! q = q / (.norm. q)

  ! end subroutine normalise_rotmatrix

  ! function normq(q)

    ! The 'reduced norm' of a rotmatrix is just the arithmetic sum
    ! of the squares of the it's elements

    ! real :: normq
    ! type (rotmatrix), intent(in) :: q

    ! normq = sum(q%q**2)

  ! end function normq
  
  function as_matrix_rotmatrix(R)

    ! Return rotmatrix as a real array -- useful for printing
    ! the order of elements is w, x, y, z

    real :: as_matrix_rotmatrix(3,3)
    type (rotmatrix) :: R

    as_matrix_rotmatrix = real(R%R,rs_kind)

  end function as_matrix_rotmatrix
  
  function coerce_matrix_to_rotmatrix(matrix)

    ! Coerce a matrix to a rotmatrix

    real, dimension(3,3) :: matrix
    type (rotmatrix)     :: coerce_matrix_to_rotmatrix
    ! integer :: i
    ! logical :: inverse = .false.

    ! Check if the matrix contains an inversion
    ! if ((.det. matrix) < 0.0 ) then
       ! matrix = -1 * matrix
    !    inverse = .true.
    ! end if

    ! Use default assignment (defined below)
    ! coerce_matrix_to_rotmatrix = matrix
    call new(coerce_matrix_to_rotmatrix, matrix=matrix)

  end function coerce_matrix_to_rotmatrix
  
  function coerce_caxis_singang_to_rotmat(axis, angle) result(rmat)

    ! Coerce a coordinate axis (x=1,y=2,z=3) and a angle to a rotmatrix. 
    ! Nicked from spicelib.

    ! Interface variables
    integer, intent(in)            :: axis
    real(kind=rs_kind), intent(in) :: angle
    
    ! Return type
    type (rotmatrix) :: rmat

    rmat = as_rotmatrix(axis, real(angle,rd_kind))

  end function coerce_caxis_singang_to_rotmat
  
  function coerce_caxis_angle_to_rotmatrix(axis, angle) result(rmat)

    ! Coerce a coordinate axis (x=1,y=2,z=3) and a angle to a rotmatrix. 
    ! Nicked from spicelib.

    ! Interface variables
    integer, intent(in)            :: axis
    real(kind=rd_kind), intent(in) :: angle
    
    ! Return type
    type (rotmatrix) :: rmat

    integer :: temp, indexs(5), i1, i2, i3
    real(kind=rd_kind) :: s, c, R(3,3)

    save indexs

    indexs = (/ 3, 1, 2, 3, 1 /)

    ! Get the sine and cosine of ANGLE

    s = sin(angle)
    c = cos(angle)

    ! Get indices for axes. The first index is for the axis of rotation.
    ! The next two axes follow in right hand order (XYZ).  First get the
    ! non-negative value of IAXIS mod 3 .
    
    temp = mod ( mod(axis,3) + 3, 3 )
 
    i1 = indexs( temp + 1 )
    i2 = indexs( temp + 2 )
    i3 = indexs( temp + 3 )
 
    ! Construct the rotation matrix

    ! R(i1,i1) = 1.0; R(i1,i2) = 0.0; R(i1,i3) = 0.0
    ! R(i2,i1) = 0.0; R(i2,i2) = c;   R(i2,i3) = s
    ! R(i3,i1) = 0.0; R(i3,i2) = -s;  R(i3,i3) = c

    R(i1,i1) = 1.0; R(i1,i2) = 0.0; R(i1,i3) = 0.0
    R(i2,i1) = 0.0; R(i2,i2) = c;   R(i2,i3) = -s
    R(i3,i1) = 0.0; R(i3,i2) = s;   R(i3,i3) = c

    ! Call constructor
    call new(rmat, R)

  end function coerce_caxis_angle_to_rotmatrix

  function coerce_aaxis_angle_to_rotmatrix() !(axis, angle)

    type (rotmatrix)     :: coerce_aaxis_angle_to_rotmatrix

    call new(coerce_aaxis_angle_to_rotmatrix)

    !      This routine computes the result of rotating (in a right handed
    !      sense) the vector V about the axis represented by AXIS through
    !      an angle of THETA radians.
    !
    !      If W is a unit vector parallel to AXIS, then R is given by:
    !
    !          R = V + ( 1 - cos(THETA) ) Wx(WxV) + sin(THETA) (WxV)
    !
    !      where "x" above denotes the vector cross product.
    !
    !      NOTE: If the input AXIS is the zero vector R will be returned
    !      as V.
    !
    !$ Examples
    !
    !      If AXIS = ( 0, 0, 1 ) and THETA = PI/2 then the following results
    !      for R will be obtained
    !
    !           V                           R
    !      -------------             ----------------
    !      ( 1, 2, 3 )                ( -2, 1, 3 )
    !      ( 1, 0, 0 )                (  0, 1, 0 )
    !      ( 0, 1, 0 )                ( -1, 0, 0 )
    !
    !
    !      If AXIS = ( 0, 1, 0 ) and THETA = PI/2 then the following results
    !      for R will be obtained
    !
    !           V                           R
    !      -------------             ----------------
    !      ( 1, 2, 3 )                (  3, 2, -1 )
    !      ( 1, 0, 0 )                (  0, 0, -1 )
    !      ( 0, 1, 0 )                (  0, 1,  0 )
    !
    !
    !      If AXIS = ( 1, 1, 1 ) and THETA = PI/2 then the following results
    !      for R will be obtained
    !
    !           V                                     R
    !      -----------------------------      -----------------------------
    !      ( 1.0,     2.0,     3.0     )      ( 1.422.., 3.154.., 1.422.. )
    !      ( 1.422.., 3.154.., 1.422.. )      ( 3.0      2.0,     1.0     )
    !      ( 3.0      2.0,     1.0     )      ( 2.577.., 0.845.., 2.577.. )
    !      ( 2.577.., 0.845.., 2.577.. )      ( 1.0      2.0,     3.0     )
    
    !     Local Variables
    !
    ! DOUBLE PRECISION C
    ! DOUBLE PRECISION S
    ! DOUBLE PRECISION RPLANE (3)
    ! DOUBLE PRECISION P      (3)
    ! DOUBLE PRECISION V1     (3)
    ! DOUBLE PRECISION V2     (3)
    ! DOUBLE PRECISION X      (3)
    
    ! We don't need to check the bona fides of the axes we are rotating about
    ! as we know they are legit basis vectors ...
    
    ! R = V + ( 1 - cos(THETA) ) Wx(WxV) + sin(THETA) (WxV)
    ! 
    !     Just in case the user tries to rotate about the zero vector -
    !     check, and if so return the input vector
    !
    !      IF ( VNORM( AXIS ) .EQ. 0.D0 ) THEN
    !         CALL MOVED (V,3,R)
    !         RETURN
    !      END IF
    ! 
    !
    !     Compute the unit vector that lies in the direction of the
    !     AXIS.  Call it X.
    !
    !      CALL VHAT ( AXIS, X )
    ! 
    !
    !     Compute the projection of V onto AXIS.  Call it P.
    !
    !      CALL VPROJ ( V, X, P )
    ! 
    !
    !     Compute the component of V orthogonal to the AXIS.  Call it V1.
    !
    !      CALL VSUB  ( V, P,  V1 )
    ! 
    !
    !     Rotate V1 by 90 degrees about the AXIS and call the result V2.
    !
    !      CALL VCRSS  ( X, V1, V2)
    ! 
    !
    !     Compute COS(THETA)*V1 + SIN(THETA)*V2. This is V1 rotated about
    !     the AXIS in the plane normal to the axis, call the result RPLANE
    !
    !      C = DCOS (THETA)
    !      S = DSIN (THETA)
    ! 
    !      CALL VLCOM ( C, V1, S, V2, RPLANE )
    ! 
    !
    !     Add the rotated component in the normal plane to AXIS to the
    !     projection of V onto AXIS (P) to obtain R.
    !
    !      CALL VADD ( RPLANE, P, R )
    ! 
    !
    !      RETURN
    !      END
    
  end function coerce_aaxis_angle_to_rotmatrix
  
  
  !
  ! ASSIGNMENT
  !

  subroutine assign_doublematrix(R,matrix)

    ! Assign to rotmatrix the elements in array

    type (rotmatrix), intent(out)  :: R
    real(kind=rd_kind), intent(in) :: matrix(3,3)

    ! Use the rotmatrix constructor to assign values to q
    call new(R,matrix)

  end subroutine assign_doublematrix

  subroutine assign_matrix(R,matrix)

    ! Assign to rotmatrix the elements in array

    type (rotmatrix), intent(out)  :: R
    real(kind=rs_kind), intent(in) :: matrix(3,3)

    ! Use the rotmatrix constructor to assign values to q
    call new(R,matrix)

  end subroutine assign_matrix

  !
  ! TRANSPOSITION
  !

  function transposer(R)

    ! Return conjugate of rotmatrix 

    type (rotmatrix)  :: transposer
    type (rotmatrix), intent(in) :: R

    transposer = transpose(R%R)

  end function transposer

  !
  ! INVERSION
  !

  function invr(R)

    ! Return inverse of a rotmatrix 

    type (rotmatrix)  :: invr
    type (rotmatrix), intent(in) :: R

    ! The inverse of an orthogonal matrix (i.e. a rotation matrix) is
    ! the transpose
    invr = transposer(R)

  end function invr

  !
  ! DETERMINANT
  !

  real function detr(R)

    ! Determinant of a rotation matrix
 
    type (rotmatrix), intent(in) :: R

    detr = detm(R%R)

  end function detr

  real function detm(matrix)

    ! Determinant of a 3x3 matrix
 
    real(kind=rd_kind), intent(in) :: matrix(3,3)

    detm = matrix(1,1) * ( matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2) ) & 
         - matrix(1,2) * ( matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1) ) & 
         + matrix(1,3) * ( matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1) )

  end function detm

  !
  ! IS A ROTATION?
  !

  function isrotrm(R)

    ! Determinant of a rotation matrix
 
    logical :: isrotrm
    type (rotmatrix), intent(in) :: R

    isrotrm = isrotm(R%R)

  end function isrotrm

  function isrotm(matrix)

    ! Determine if matrix is a valid rotation matrix
 
    ! Interface variable
    real(kind=rd_kind), intent(in)   :: matrix(3,3)
    
    ! Local variables
    logical :: isrotm, normok, detok
    real(kind=rd_kind) :: unit(3,3), mag(3)
    real(kind=rd_kind), parameter :: dtol = 0.01, ntol=0.01

    !     The function returns the value .TRUE. if and only if M is found
    !     to be a rotation matrix.  The criteria that M must meet are:
    !
    !
    !        1) The norm of each column of M must satisfy the relation
    !
    !              1.0 - NTOL  <   || column ||   <  1.0 + NTOL.
    !                           -                  -
    !
    !        2) The determinant of the matrix whose columns are the
    !           unitized columns of M must satisfy
    !
    !              1.0 - DTOL  <   determinant   <  1.0 + DTOL.
    !                           -                 -
    !     This routine is an error checking `filter'; its purpose is to
    !     detect gross errors, such as uninitialized matrices.  Matrices
    !     that do not pass the tests used by this routine hardly qualify as
    !     rotation matrices.  The test criteria can be adjusted by varying
    !     the parameters NTOL and DTOL.
    !
    !     A property of rotation matrices is that their columns form a
    !     right-handed, orthonormal basis in 3-dimensional space.  The
    !     converse is true:  all 3x3 matrices with this property are
    !     rotation matrices.
    !
    !     An ordered set of three vectors V1, V2, V3 forms a right-handed,
    !     orthonormal basis if and only if
    !
    !        1)   || V1 ||  =  || V2 ||  =  || V3 ||  =  1
    !
    !        2)   V3 = V1 x V2.  Since V1, V2, and V3 are unit vectors,
    !             we also have
    !
    !             < V3, V1 x V2 > = 1.
    !
    !             This quantity is the determinant of the matrix whose
    !             colums are V1, V2 and V3.
    !
    !     When finite precision numbers are used, rotation matrices will
    !     usually fail to satisfy these criteria exactly.  We must use
    !     criteria that indicate approximate conformance to the criteria
    !     listed above.  We choose
    !
    !        1)   |   || Vi ||  -  1   |   <   NTOL,  i = 1, 2, 3.
    !                                      -
    !
    !        2)   Let
    !
    !                       Vi
    !                Ui = ------ ,   i = 1, 2, 3.
    !                     ||Vi||
    !
    !             Then we require
    !
    !                | < U3, U1 x U2 > - 1 |  <  DTOL;
    !                                         -
    !
    !             equivalently, letting U be the matrix whose columns
    !             are U1, U2, and U3, we insist on
    !
    !                | det(U) - 1 |  <  DTOL.
    !                                _

    !     Tolerances must be non-negative.
    !
    !
    !     The columns of M must resemble unit vectors.  If the norms are
    !     outside of the allowed range, M is not a rotation matrix.
    !
    !     Also, the columns of M are required to be pretty nearly
    !     orthogonal.  The discrepancy is gauged by taking the determinant
    !     of the matrix UNIT, computed below, whose columns are the
    !     unitized columns of M.

    ! Initialise return to false
    isrotm = .false.

    ! Generate an array whose values are the rms of the columns of the
    ! matrix
    mag = sqrt(sum(matrix**2,dim=1))

    ! Check that all the columns of the matrix have a magnitude of 1.0 +/- ntol
    normok = ( maxval(abs(mag - 1.0)) < ntol )

    ! Exit from the routine without setting isrotm to true
    if (.not. normok) return

    ! Make a "unitized" version of matrix by dividing the columns by their 
    ! rms magnitude
    unit(1:3,1) = matrix(1:3,1)/mag(1)
    unit(1:3,2) = matrix(1:3,2)/mag(2)
    unit(1:3,3) = matrix(1:3,3)/mag(3)

    ! Check that the determinant of this unitized matrix is +/- 1.0 
    detok  = ( abs( abs(.det. unit) - 1.0 ) < DTOL )
    
    ! Exit from the routine without setting isrotm to true
    if ( .not. detok ) return

    ! We have passed all the tests -- this must be a value rotation matrix
    isrotm = .true.
    
  end function isrotm

  !
  ! ROTATION
  !

  subroutine rrotv(R,v)

    ! Rotate a vector by a rotation matrix
 
    type (rotmatrix), intent(in)      :: R
    real(kind=rd_kind), intent(inout) :: v(3)

    v = matmul(R%R, v)

  end subroutine rrotv

  subroutine rrotm(R,m)

    ! Rotate a matrix (of vectors) by a rotation matrix
 
    type (rotmatrix), intent(in)                      :: R
    real(kind=rd_kind), dimension(:,:), intent(inout) :: m

    m = matmul(R%R, m)

  end subroutine rrotm

  subroutine rrotsinglev(R,v)

    ! Rotate a vector by a rotation matrix
 
    type (rotmatrix), intent(in)      :: R
    real(kind=rs_kind), intent(inout) :: v(3)

    v = matmul(R%R, v)

  end subroutine rrotsinglev

  subroutine rrotsinglem(R,m)

    ! Rotate a matrix (of vectors) by a rotation matrix
 
    type (rotmatrix), intent(in)                      :: R
    real(kind=rs_kind), dimension(:,:), intent(inout) :: m

    m = matmul(R%R, m)

  end subroutine rrotsinglem



  !
  ! MULTIPLICATION
  !

  function rmultr(RL,RR)

    ! Return rotmatrix product RL * RR

    type (rotmatrix)  :: rmultr
    type (rotmatrix), intent(in) :: RL, RR
    
    ! Use default assignment for result of matrix multiplication.
    rmultr = matmul( RL%R, RR%R )

  end function rmultr

  function smultr(s, R)

    ! Return rotmatrix multiplied by a scalar

    type (rotmatrix)  :: smultr
    real, intent(in) :: s
    type (rotmatrix), intent(in) :: R

    smultr = R
    smultr%R = s * smultr%R

  end function smultr

  function rmults(R, s)

    ! Return rotmatrix multiplied by a scalar
    ! (same as above, but swapped order of rmat & scal)

    type (rotmatrix)  :: rmults
    real, intent(in) :: s
    type (rotmatrix), intent(in) :: R

    rmults = smultr(s, R)

  end function rmults

  function rmulti(R, i)

    ! Return rotmatrix multiplied by an integer

    type (rotmatrix)  :: rmulti
    integer, intent(in) :: i
    type (rotmatrix), intent(in) :: R

    rmulti = smultr(real(i),R)

  end function rmulti

  function imultr(i, R)

    ! Return rotmatrix multiplied by an integer

    type (rotmatrix)  :: imultr
    integer, intent(in) :: i
    type (rotmatrix), intent(in) :: R

    imultr = smultr(real(i),R)

  end function imultr

  !
  ! ADDITION
  !

  function rplusr(RL, RR)

    ! Return rotmatrix added to rotmatrix

    type (rotmatrix)  :: rplusr
    type (rotmatrix), intent(in) :: RL, RR

    ! Use default assignment, i.e. the result of the addition of
    ! the two matrices will be an matrix, not a rotmatrix, so the
    ! rule for assigning matrices will be invoked
    rplusr = RL%R + RR%R

  end function rplusr

  function rpluss(R, s)

    ! Return rotmatrix added to double precision scalar

    type (rotmatrix), intent(in)   :: R
    real(kind=rd_kind), intent(in) :: s
    
    ! Return type
    type (rotmatrix)  :: rpluss

    ! Add together and use default assignment
    rpluss = R%R + s

  end function rpluss

  function rplussingles(R, s)

    ! Return rotmatrix added to single precision scalar

    type (rotmatrix), intent(in)   :: R
    real(kind=rs_kind), intent(in) :: s
    
    ! Return type
    type (rotmatrix)  :: rplussingles

    ! Add together and use default assignment
    rplussingles = R%R + s

  end function rplussingles

  function splusr(s, R)

    ! Return double precision scalar added to rotmatrix 

    type (rotmatrix), intent(in)   :: R
    real(kind=rd_kind), intent(in) :: s

    ! Return type
    type (rotmatrix) :: splusr

    ! Add rotmatrix to real
    splusr = R + s

  end function splusr

  function ssingleplusr(s, R)

    ! Return single precision scalar added to rotmatrix 

    type (rotmatrix), intent(in)   :: R
    real(kind=rs_kind), intent(in) :: s

    ! Return type
    type (rotmatrix) :: ssingleplusr

    ! Add rotmatrix to real
    ssingleplusr = R + s

  end function ssingleplusr

  function rplusi(R, i)

    ! Return rotmatrix added to integer

    type (rotmatrix)  :: rplusi
    type (rotmatrix), intent(in) :: R
    integer, intent(in) :: i

    ! Add rotmatrix to integer (coerced to real -- will call qplusr above)
    rplusi = R + real(i,rd_kind)

  end function rplusi

  function iplusr(i, R)

    ! Return integer added to rotmatrix

    type (rotmatrix)  :: iplusr
    type (rotmatrix), intent(in) :: R
    integer, intent(in) :: i

    ! Add rotmatrix to integer (coerced to real -- will call rpluss above)
    iplusr = R + real(i,rd_kind)

  end function iplusr

  function plusr(R)

    ! Return unary addition of rotmatrix

    type (rotmatrix)  :: plusr
    type (rotmatrix), intent(in) :: R

    plusr = R

  end function plusr

  !
  ! SUBTRACTION
  !

  function rminusr(RL, RR)

    ! Return rotmatrix subtracted from rotmatrix

    type (rotmatrix)  :: rminusr
    type (rotmatrix), intent(in) :: RL, RR

    ! Use default assignment, i.e. the result of the addition of
    ! the two matrices will be an matrix, not a rotmatrix, so the
    ! rule for assigning matrices will be invoked
    rminusr = RL%R - RR%R

  end function rminusr

  function rminuss(R, s)

    ! Return scalar subtracted from rotmatrix

    type (rotmatrix)  :: rminuss
    type (rotmatrix), intent(in)   :: R
    real(kind=rd_kind), intent(in) :: s

    ! Assign using default assignment
    rminuss = R%R - s

  end function rminuss

  function rminussingles(R, s)

    ! Return scalar subtracted from rotmatrix

    type (rotmatrix)  :: rminussingles
    type (rotmatrix), intent(in)   :: R
    real(kind=rs_kind), intent(in) :: s

    ! Assign using default assignment
    rminussingles = R%R - real(s,rd_kind)

  end function rminussingles

  function sminusr(s, R)

    ! Return rotmatrix subtracted from scalar

    type (rotmatrix)  :: sminusr
    type (rotmatrix), intent(in)   :: R
    real(kind=rd_kind), intent(in) :: s

    ! Add minus rotmatrix to scalar
    sminusr = -R + s

  end function sminusr

  function ssingleminusr(s, R)

    ! Return rotmatrix subtracted from scalar

    type (rotmatrix)  :: ssingleminusr
    type (rotmatrix), intent(in)   :: R
    real(kind=rs_kind), intent(in) :: s

    ! Add minus rotmatrix to scalar
    ssingleminusr = -R + s

  end function ssingleminusr

  function rminusi(R, i)

    ! Return rotmatrix added to integer

    type (rotmatrix)  :: rminusi
    type (rotmatrix), intent(in) :: R
    integer, intent(in) :: i

    ! Subtract integer from rotmatrix
    rminusi = R - real(i,rd_kind)

  end function rminusi

  function iminusr(i, R)

    ! Return integer added to rotmatrix

    type (rotmatrix)  :: iminusr
    type (rotmatrix), intent(in) :: R
    integer, intent(in) :: i

    ! Add minus rotmatrix to integer (coerced to real -- will call sminusr above)
    iminusr = real(i,rd_kind) - R

  end function iminusr

  function minusr(R)

    ! Return unary negation of rotmatrix

    type (rotmatrix)  :: minusr
    type (rotmatrix), intent(in) :: R

    minusr = R
    minusr%R = -1*minusr%R

  end function minusr

  !
  ! EQUIVALENCE
  !

  function reqr(RL, RR)

    ! Determines if two rotation matrices are equivalent (within
    ! an arbitrary tolerance)

    logical :: reqr
    type (rotmatrix), intent(in) :: RL, RR
    real(kind=rd_kind) :: tolerance = 0.01

    reqr = .false.

    ! See if the RMS difference twixt the matrices is greater than the tolerance
    if (sum(sqrt( (RL%R - RR%R)**2 )) < tolerance) reqr = .true.

  end function reqr

  function rner(RL, RR)

    ! Determines if two rotation matrices are not equivalent (within
    ! an arbitrary tolerance)

    logical :: rner
    type (rotmatrix), intent(in) :: RL, RR

    rner = .true.

    ! Just call the equivalence routine and return the opposite value
    if (RL == RR) rner = .false.

  end function rner

  subroutine reorthogonalise_matrix (matrix)
    
    ! From:
    !
    ! http://groups.google.com.au/groups?hl=en&lr=&ie=UTF-8&safe=off&selm=0Fhl6.217%24kC2.9621%40sea-read.news.verio.net
    !
      ! Treating the first three columns as vectors X, Y and Z, you need to do
    ! something like this:
    ! 
    ! Y = Z cross X
    ! X = Y cross Z
    ! Normalize(X)
    ! Normalize(Y)
    ! Normalize(Z)
    ! 
    ! This is, of course, assuming that things haven't gotten so bad that two of
    ! the column vectors are congruent.
    
    ! And from spicelib sharpr.f:
    
    ! Given a matrix that is "nearly" a rotation, adjust the columns
    ! (from left to right in the usual printed presentation of a matrix)
    ! so that the columns are numerically unit length and orthogonal.
    
    ! This routine "sharpens" the orthogonality of a potential
    ! rotation matrix.  It is intended for use in those situations
    ! in which you have a rotation matrix that may be derived
    ! from single precision inputs or that may have experienced
    ! round off errors in its construction.
      
    real(kind=rd_kind), intent(inout) :: matrix(3,3)

    ! Normalise 'X'
    matrix(:,1) = unit(matrix(:,1))
    ! 'Z' = norm(X cross Y)
    matrix(:,3) = unit(matrix(:,1) .cross. matrix(:,2))
    ! 'Y' = norm(Z cross X)
    matrix(:,2) = unit(matrix(:,3) .cross. matrix(:,1))

  end subroutine reorthogonalise_matrix

  subroutine reorthogonalise_rotmatrix (R)

    ! Wrapper to routine above
    type (rotmatrix), intent(inout) :: R

    call reorthogonalise(R%R)

  end subroutine reorthogonalise_rotmatrix
  
  
end module rotmatrix_class
