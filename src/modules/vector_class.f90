module vector_class

  ! A great deal (almost all!) of this is adapted from the spicelib 
  ! routines (NASA), and so work really really well (how many NASA
  ! things have fallen out of the sky ....)

  use fundamental_constants
  use file_functions, only: stderr
  use precision

  implicit none

  private

  character(len=*), parameter :: version = "$Id: vector_class.f90,v 1.2 2006/09/29 04:37:59 aidan Exp $"

  !! Implements vector maths functions, as well as making a new plane type.

  !! $Log: vector_class.f90,v $
  !! Revision 1.2  2006/09/29 04:37:59  aidan
  !! Added support for the precision module. Made double and single precision
  !! versions of all routines.
  !!
  !! Fixed a bug with the normal and constant routines. They an ambiguity in
  !! the module procedure specification because they did not differ in their
  !! non-optional arguments.
  !!
  !! Revision 1.1  2004/08/05 05:16:46  aidan
  !! Initial revision
  !!

  real(kind=rd_kind) :: tolerance = 1e-4

  ! A plane is defined by a unit normal and the distance from the plane
  ! to the origin (constant). Every point X in the plane satisifies
  !
  !           X .dot. normal =  constant
  !
  ! The normal vector, scaled by the constant, is the closest point in 
  ! the plane to the origin.
  type plane_object 
     ! private
     real(kind=rd_kind) :: normal(3), constant
  end type plane_object

  ! Define assignment from a filename
  ! interface assignment (=)
  !    module procedure psv2pl
  ! end interface

  ! There are three ways of defining a plane 
  ! 1. With a normal and a constant
  ! 2. With a normal and a point on the plane
  ! 3. With two spanning vectors and a point on the plane
  ! There are routines to account for any combination of single and double
  ! precision with the above 3 methods
  interface new
     module procedure new_plane, new_plane_single, new_plane_single_double, new_plane_double_single, &
          normal_point_to_plane, normal_point_to_plane_single, normal_point_to_plane_sd, &
          normal_point_to_plane_ds, & 
          point_span_to_plane_sss, point_span_to_plane_ssd, point_span_to_plane_sdd, & 
          point_span_to_plane_sds, point_span_to_plane_dds, point_span_to_plane_dss, & 
          point_span_to_plane_dsd, point_span_to_plane
  end interface

  ! Define vector norm operator -- only monadic
  interface operator (.norm.)
     module procedure normv, normv_wrap
  end interface

  ! Define cross product operator -- only diadic
  interface operator (.cross.)
     module procedure cross_product_double, cross_product_single, & 
          cross_product_single_double, cross_product_double_single
  end interface

  ! Define dot product operator -- only diadic
  interface operator (.dot.)
     module procedure dot_product_double, dot_product_single, dot_product_single_double, dot_product_double_single
  end interface

  ! Define angle operator -- only diadic
  interface operator (.angle.)
     module procedure vsep, vsep_single, vsep_double_single, vsep_single_double
  end interface

  ! Define vector equivalence operator. Can't overload the == operator
  ! as this is functionally equivalent to the == operator apparently
  ! (even though trying an equivalence statement on two vectors dosen't work!)
  interface operator (.veqv.)
     module procedure veqv, veqv_single, veqv_single_double, veqv_double_single
  end interface

  ! Define vector non-equivalence operator
  interface operator (.vnev.)
     module procedure vnev, vnev_single, vnev_single_double, vnev_double_single
  end interface

  ! Define plane equivalence operator
  interface operator (==)
     module procedure plane_eq_plane
  end interface

  ! Define plane equivalence operator
  interface operator (/=)
     module procedure plane_neq_plane
  end interface

  ! Define a subroutine interface to normalise a vector
  interface normalise
     module procedure normalise_vector, normalise_vector_wrap
  end interface

  ! Define a functional interface to normalise a vector
  interface unit
     module procedure normalise_vector_function, normalise_vector_function_wrap
  end interface

  ! Define angle procedure
  interface angle
     module procedure vsep, vsep_single, vsep_double_single, vsep_single_double
  end interface

  ! Interface to project function -- projects a vector on to a plane
  interface project
     module procedure project_vector_onto_plane, project_vector_onto_plane_wrap
  end interface

  ! Access routine for normal array of plane
  interface normal
     module procedure get_normal, set_normal, set_normal_single
  end interface

  ! Access routine for constant value of plane
  interface constant
     module procedure get_constant, set_constant, set_constant_single
  end interface

  !!!!!!!!!!!!!!!!
  ! Public stuff !
  !!!!!!!!!!!!!!!!

  ! types
  public :: plane_object

  ! routines
  public :: normalise, angle, unit, new, project, normal, constant

  ! new operators
  public :: operator(.norm.), operator(.cross.), operator(.dot.), operator(.angle.)
  public :: operator(.veqv.), operator(.vnev.), operator(==), operator(/=)

contains

  subroutine new_plane ( plane, normal, const )

    ! Make a plane from a normal vector and a constant.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: normal
    real(kind=rd_kind), intent(in)                 :: const

    ! normal, const are, respectively, a normal vector and constant
    ! defining a plane.  normal need not be a unit vector.  Let the 
    ! The geometric plane is the set of vectors X in three-dimensional 
    ! space that satisfy
    !
    !                       X .dot.  normal =  const

    real(kind=rd_kind) :: mag

    mag = .norm. normal
    
    ! The normal vector must be non-zero.
    if (mag == 0.0) THEN
       write(stderr,*) 'VECTOR_CLASS :: Normal must be non-zero.'
       stop
    end if

    plane%normal  =  normal / mag

    ! To find the plane constant corresponding to the unitized normal
    ! vector, we observe that
    !
    !      X .dot. normal = const,
    !
    !     so
    !
    !     X .dot. (normal / || normal || )   =   const / || normal ||
    !
    !
    plane%constant  =  const / mag
    
    ! The constant should be the distance of the plane from the
    ! origin.  If the constant is negative, negate both it and the
    ! normal vector.
    if ( plane%constant < 0.0 ) then
       plane%constant = -plane%constant
       plane%normal   = -plane%normal
    end if
    return
  end subroutine new_plane

  subroutine new_plane_single ( plane, normal, const )

    ! Make a plane from a normal vector and a constant.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: normal
    real(kind=rs_kind), intent(in)                 :: const

    call new(plane, real(normal,rd_kind), real(const,rd_kind))

    return
  end subroutine new_plane_single

  subroutine new_plane_single_double ( plane, normal, const )

    ! Make a plane from a normal vector and a constant.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: normal
    real(kind=rd_kind), intent(in)                 :: const

    call new(plane, real(normal,rd_kind), const)

    return
  end subroutine new_plane_single_double

  subroutine new_plane_double_single ( plane, normal, const )

    ! Make a plane from a normal vector and a constant.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: normal
    real(kind=rs_kind), intent(in)                 :: const

    call new(plane, normal, real(const,rd_kind))

    return
  end subroutine new_plane_double_single

  subroutine point_span_to_plane ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: point, span1, span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, point) 

    return
  end subroutine point_span_to_plane

  subroutine point_span_to_plane_sss ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: point, span1, span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, real(point,rd_kind)) 

    return
  end subroutine point_span_to_plane_sss

  subroutine point_span_to_plane_ssd ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: point
    real(kind=rs_kind), dimension(3), intent(in)   :: span1
    real(kind=rd_kind), dimension(3), intent(in)   :: span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, point) 

    return
  end subroutine point_span_to_plane_ssd

  subroutine point_span_to_plane_sdd ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: point
    real(kind=rd_kind), dimension(3), intent(in)   :: span1
    real(kind=rd_kind), dimension(3), intent(in)   :: span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, point) 

    return
  end subroutine point_span_to_plane_sdd

  subroutine point_span_to_plane_sds ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: point
    real(kind=rd_kind), dimension(3), intent(in)   :: span1
    real(kind=rs_kind), dimension(3), intent(in)   :: span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, real(point,rd_kind)) 

    return
  end subroutine point_span_to_plane_sds

  subroutine point_span_to_plane_dds ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: point
    real(kind=rd_kind), dimension(3), intent(in)   :: span1
    real(kind=rs_kind), dimension(3), intent(in)   :: span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, real(point,rd_kind)) 

    return
  end subroutine point_span_to_plane_dds

  subroutine point_span_to_plane_dss ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: point
    real(kind=rs_kind), dimension(3), intent(in)   :: span1
    real(kind=rs_kind), dimension(3), intent(in)   :: span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, real(point,rd_kind)) 

    return
  end subroutine point_span_to_plane_dss

  subroutine point_span_to_plane_dsd ( plane, span1, span2, point )

    ! Make a plane from a point and two spanning vectors. Makes a
    ! normal from the two span vectors and calls new_plane

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: point
    real(kind=rs_kind), dimension(3), intent(in)   :: span1
    real(kind=rd_kind), dimension(3), intent(in)   :: span2

    ! The spanning vectors span1 and span2 must be linearly independent, but 
    ! they need not be orthogonal or unitized.

    ! The cross product of SPAN1 and SPAN2 is normal vector, or possibly its 
    ! inverse.
    call new(plane, span1 .cross. span2, point) 

    return
  end subroutine point_span_to_plane_dsd

  subroutine normal_point_to_plane ( plane, normal, point )

    ! Make a plane from a normal vector and a point.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: normal, point

    ! normal and point are, respectively, a normal vector and point that
    ! define a plane in three-dimensional space.  normal need not be a 
    ! unit vector. Let the symbol .dot. indicate the inner product of 
    ! vectors a and b; then the geometric plane is the set of vectors X
    ! in three-dimensional space that satisfy
    !
    !                       (X - POINT) .dot. NORMAL =  0.

    call new(plane, normal, unit(normal) .dot. point)
    return

  end subroutine normal_point_to_plane

  subroutine normal_point_to_plane_single ( plane, normal, point )

    ! Make a plane from a normal vector and a point.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: normal, point

    ! normal and point are, respectively, a normal vector and point that
    ! define a plane in three-dimensional space.  normal need not be a 
    ! unit vector. Let the symbol .dot. indicate the inner product of 
    ! vectors a and b; then the geometric plane is the set of vectors X
    ! in three-dimensional space that satisfy
    !
    !                       (X - POINT) .dot. NORMAL =  0.

    call new(plane, real(normal,rd_kind), unit(normal) .dot. point)
    return

  end subroutine normal_point_to_plane_single

  subroutine normal_point_to_plane_sd ( plane, normal, point )

    ! Make a plane from a normal vector and a point.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rs_kind), dimension(3), intent(in)   :: normal
    real(kind=rd_kind), dimension(3), intent(in)   :: point

    ! normal and point are, respectively, a normal vector and point that
    ! define a plane in three-dimensional space.  normal need not be a 
    ! unit vector. Let the symbol .dot. indicate the inner product of 
    ! vectors a and b; then the geometric plane is the set of vectors X
    ! in three-dimensional space that satisfy
    !
    !                       (X - POINT) .dot. NORMAL =  0.

    call new(plane, real(normal,rd_kind), unit(normal) .dot. point)
    return

  end subroutine normal_point_to_plane_sd

  subroutine normal_point_to_plane_ds ( plane, normal, point )

    ! Make a plane from a normal vector and a point.

    ! Output variables
    type (plane_object), intent(out) :: plane

    ! Input variables
    real(kind=rd_kind), dimension(3), intent(in)   :: normal
    real(kind=rs_kind), dimension(3), intent(in)   :: point

    ! normal and point are, respectively, a normal vector and point that
    ! define a plane in three-dimensional space.  normal need not be a 
    ! unit vector. Let the symbol .dot. indicate the inner product of 
    ! vectors a and b; then the geometric plane is the set of vectors X
    ! in three-dimensional space that satisfy
    !
    !                       (X - POINT) .dot. NORMAL =  0.

    call new(plane, normal, unit(normal) .dot. point)
    return

  end subroutine normal_point_to_plane_ds

  function project_vector_onto_plane ( vector, plane ) result(vproj)

    ! Project a vector onto a specified plane, orthogonally.

    ! Output variable
    real(kind=rd_kind), dimension(3) :: vproj

    real(kind=rd_kind), dimension(3), intent(in) :: vector
    type (plane_object), intent(in)              :: plane

    ! Projecting a vector V orthogonally onto a plane can be thought of
    ! as finding the closest vector in the plane to V.  This `closest
    ! vector' always exists; it may be coincident with the original
    ! vector.
    !
    ! vector differs from its projection onto plane by some multiple of
    ! normal.  That multiple is
    !
    !    ((vector - vproj) .dot. normal) *  normal
    !
    !        =   ((vector .dot. normal) - (vproj .dot. normal))  *  normal
    !
    !        =   ((vector .dot. normal) - const)  *  normal
    !
    ! Subtracting this multiple of normal from vector yields vproj.

    vproj = vector + (plane%constant - (vector .dot. plane%normal)) * plane%normal

    return
  end function project_vector_onto_plane

  function project_vector_onto_plane_wrap ( vector, plane ) result(vproj)

    ! Project a vector onto a specified plane, orthogonally.

    ! Output variable
    real(kind=rd_kind), dimension(3) :: vproj

    real(kind=rs_kind), dimension(3), intent(in) :: vector
    type (plane_object), intent(in)              :: plane

    vproj = project( real(vector,rd_kind), plane ) 

  end function project_vector_onto_plane_wrap

  subroutine normalise_vector(v)

    ! Normalise the vector by invoking the normalise function

    real(kind=rd_kind), dimension(3), intent(inout) :: v

    v = unit(v)

  end subroutine normalise_vector

  subroutine normalise_vector_wrap(v)

    ! Normalise the vector by invoking the normalise function

    real(kind=rs_kind), dimension(3), intent(inout) :: v

    v = real(unit(real(v,rd_kind)),rs_kind)

  end subroutine normalise_vector_wrap

  function normalise_vector_function(v) result(vnorm)

    ! Normalise the vector by dividing the vector by it's
    ! norm

    real(kind=rd_kind), dimension(3) :: vnorm
    real(kind=rd_kind), dimension(3), intent(in) :: v

    real(kind=rd_kind) :: vmag

    vmag = .norm. v

    if (vmag > 0.d0) then
       vnorm = v / vmag
    else
       vnorm = 0.d0
    end if

  end function normalise_vector_function

  function normalise_vector_function_wrap(v) result(vnorm)

    ! Normalise the vector by invoking the double precision
    ! functional version
    real(kind=rs_kind), dimension(3), intent(in) :: v

    ! Return type
    real(kind=rd_kind), dimension(3) :: vnorm

    vnorm = unit(real(v,rd_kind))

  end function normalise_vector_function_wrap

  function normv(v)

    ! The norm of a vector is just the square root of the
    ! arithmetic sum of the squares of it's elements
    real(kind=rd_kind), dimension(3), intent(in) :: v

    ! Return value
    real(kind=rd_kind) :: normv
    real(kind=rd_kind) :: dummy

    real(kind=rd_kind) :: vmax

    vmax = maxval(abs(v))

    ! Use an intermediate variable so we don't lose precision
    ! dummy = sqrt(sum(real(v,8)**2))
    ! normv = real(dummy)
    ! normv = sqrt(sum(v**2))
    if (vmax == 0.0) then
       normv = 0.0
    else
       normv = vmax * sqrt( sum((v/vmax)**2) )
    end if

  end function normv

  function normv_wrap(v)

    ! The norm of a vector is just the square root of the
    ! arithmetic sum of the squares of it's elements
    real(kind=rs_kind), dimension(3), intent(in) :: v

    ! Return value
    real(kind=rd_kind) :: normv_wrap

    normv_wrap = .norm. real(v,rd_kind)

  end function normv_wrap

  function cross_product_double(v1,v2) result(xproduct)

    ! Return cross product of two 3d vectors

    real(kind=rd_kind), dimension(3)  :: xproduct
    real(kind=rd_kind), dimension(3), intent(in) :: v1, v2

    xproduct = (/ v1(2)*v2(3)-v1(3)*v2(2),       &
         v1(3)*v2(1)-v1(1)*v2(3),       &
         v1(1)*v2(2)-v1(2)*v2(1) /)

  end function cross_product_double

  function cross_product_single(v1,v2) result(xproduct)

    ! Return cross product of two 3d vectors

    real(kind=rd_kind), dimension(3)  :: xproduct
    real(kind=rs_kind), dimension(3), intent(in) :: v1, v2

    xproduct = real(v1,rd_kind) .cross. real(v2,rd_kind) 

  end function cross_product_single

  function cross_product_single_double(v1,v2) result(xproduct)

    ! Return cross product of two 3d vectors

    real(kind=rd_kind), dimension(3)  :: xproduct
    real(kind=rs_kind), dimension(3), intent(in) :: v1
    real(kind=rd_kind), dimension(3), intent(in) :: v2

    xproduct = real(v1,rd_kind) .cross. v2

  end function cross_product_single_double

  function cross_product_double_single(v1,v2) result(xproduct)

    ! Return cross product of two 3d vectors

    real(kind=rd_kind), dimension(3)  :: xproduct
    real(kind=rd_kind), dimension(3), intent(in) :: v1
    real(kind=rs_kind), dimension(3), intent(in) :: v2

    xproduct = v1 .cross. real(v2,rd_kind)

  end function cross_product_double_single
  
  function dot_product_double(v1,v2) result(dotp)

    ! Return dot product of two 3d vectors

    real(kind=rd_kind) :: dotp
    real(kind=rd_kind), dimension(3), intent(in) :: v1, v2

    dotp = dot_product(v1,v2) 

  end function dot_product_double
  
  function dot_product_single(v1,v2) result(dotp)

    ! Return dot product of two 3d vectors

    real(kind=rs_kind) :: dotp
    real(kind=rs_kind), dimension(3), intent(in) :: v1, v2

    dotp = dot_product(v1,v2) 

  end function dot_product_single
  
  function dot_product_double_single(v1,v2) result(dotp)

    ! Return dot product of two 3d vectors

    real(kind=rd_kind) :: dotp
    real(kind=rd_kind), dimension(3), intent(in) :: v1
    real(kind=rs_kind), dimension(3), intent(in) :: v2

    dotp = dot_product(v1,v2) 

  end function dot_product_double_single
  
  function dot_product_single_double(v1,v2) result(dotp)

    ! Return dot product of two 3d vectors

    real(kind=rd_kind) :: dotp
    real(kind=rs_kind), dimension(3), intent(in) :: v1
    real(kind=rd_kind), dimension(3), intent(in) :: v2

    dotp = dot_product(v1,v2) 

  end function dot_product_single_double

  function vsep ( v1, v2 )
    
    ! Find the separation angle in radians between two 
    ! 3-dimensional vectors.  This angle is defined as zero
    ! if either vector is zero.

    real(kind=rd_kind), dimension(3), intent(in) :: v1, v2
    real(kind=rd_kind) :: vsep
    
    ! VSEP    is the angle between V1 and V2 expressed in radians.
    !         VSEP is strictly non-negative.  If either V1 or V2 is
    !         the zero vector, then VSEP is defined to be 0 radians.

    !      In the plane, it is a simple matter to calculate the angle
    !      between two vectors once the two vectors have been made to be
    !      unit length.  Then, since the two vectors form the two equal
    !      sides of an isosceles triangle, the length of the third side
    !      is given by the expression
    !
    !            LENGTH = 2.0 * SINE ( VSEP/2.0 )
    !
    !      The length is given by the magnitude of the difference of the
    !      two unit vectors
    !
    !            LENGTH = NORM ( U1 - U2 )
    !
    !      Once the length is found, the value of VSEP may be calculated
    !      by inverting the first expression given above as
    !
    !            VSEP = 2.0 * ARCSINE ( LENGTH/2.0 )
    !
    !      This expression becomes increasingly unstable when VSEP gets
    !      larger than PI/2 or 90 degrees.  In this situation (which is
    !      easily detected by determining the sign of the dot product of
    !      V1 and V2) the supplementary angle is calculated first and
    !      then VSEP is given by
    !
    !            VSEP = PI - SUPPLEMENTARY_ANGLE

    !     The following declarations represent, respectively:
    !        Magnitudes of V1, V2
    !

    real(kind=rd_kind) :: dmag1, dmag2, u1(3), u2(3)
    real(kind=rd_kind) :: vtemp(3)
    
    !  Calculate the magnitudes of v1 and v2; if either is 0, vsep = 0

    dmag1 = .norm. v1
    if ( dmag1 .eq. 0.0 ) then
       vsep = 0.0
       return
    else
       u1 = v1 / dmag1
    end if

    dmag2 = .norm. v2
    if ( dmag2 .eq. 0.0 ) then
       vsep = 0.0
       return
    else
       u2 = v2 / dmag2
    end if

    ! Use the most numerically stable form of asin ...
    if ( (u1 .dot. u2) > 0 ) then
       vtemp = u1 - u2
       vsep = 2.00 * asin (0.50 * (.norm. vtemp))
    else if ( (u1 .dot. u2) < 0 ) then
       vtemp = u1 + u2
       vsep = pi - 2.00 * asin (0.50 * (.norm. vtemp))
    else
       vsep = pi / 2.00
    end if

    return

  end function vsep

  function vsep_single ( v1, v2 ) result(vsep)
    
    ! Find the separation angle in radians between two 
    ! 3-dimensional vectors.  This angle is defined as zero
    ! if either vector is zero.

    real(kind=rs_kind), dimension(3), intent(in) :: v1, v2
    real(kind=rd_kind) :: vsep

    vsep = angle(real(v1,rd_kind), real(v2,rd_kind))

  end function vsep_single

  function vsep_single_double ( v1, v2 ) result(vsep)
    
    ! Find the separation angle in radians between two 
    ! 3-dimensional vectors.  This angle is defined as zero
    ! if either vector is zero.

    real(kind=rs_kind), dimension(3), intent(in) :: v1
    real(kind=rd_kind), dimension(3), intent(in) :: v2
    real(kind=rd_kind) :: vsep

    vsep = angle(real(v1,rd_kind), v2)

  end function vsep_single_double

  function vsep_double_single ( v1, v2 ) result(vsep)
    
    ! Find the separation angle in radians between two 
    ! 3-dimensional vectors.  This angle is defined as zero
    ! if either vector is zero.

    real(kind=rd_kind), dimension(3), intent(in) :: v1
    real(kind=rs_kind), dimension(3), intent(in) :: v2
    real(kind=rd_kind) :: vsep

    vsep = angle(v1, real(v2,rd_kind))

  end function vsep_double_single
    
  !
  ! EQUIVALENCE
  !

  function veqv(vL, vR)

    ! Determines if two vectors are equivalent (within
    ! an arbitrary tolerance)

    logical :: veqv
    real(kind=rd_kind), dimension(3), intent(in) :: vL, vR

    ! See if the RMS difference twixt the vectors is greater than the tolerance
    veqv = (all(abs(vL - vR) < tolerance))

  end function veqv

  function veqv_single(vL, vR) result(veqv)

    ! Determines if two vectors are equivalent (within
    ! an arbitrary tolerance)

    logical :: veqv
    real(kind=rs_kind), dimension(3), intent(in) :: vL, vR

    ! See if the RMS difference twixt the vectors is greater than the tolerance
    veqv = (all(abs(vL - vR) < tolerance))

  end function veqv_single

  function veqv_single_double(vL, vR) result(veqv)

    ! Determines if two vectors are equivalent (within
    ! an arbitrary tolerance)

    logical :: veqv
    real(kind=rs_kind), dimension(3), intent(in) :: vL
    real(kind=rd_kind), dimension(3), intent(in) :: vR

    ! See if the RMS difference twixt the vectors is greater than the tolerance
    veqv = (all(abs(real(vL,rd_kind) - vR) < tolerance))

  end function veqv_single_double

  function veqv_double_single(vL, vR) result(veqv)

    ! Determines if two vectors are equivalent (within
    ! an arbitrary tolerance)

    logical :: veqv
    real(kind=rd_kind), dimension(3), intent(in) :: vL
    real(kind=rs_kind), dimension(3), intent(in) :: vR

    ! See if the RMS difference twixt the vectors is greater than the tolerance
    veqv = (all(abs(vL - real(vR,rd_kind)) < tolerance))

  end function veqv_double_single

  function vnev(vL, vR)

    ! Determines if two vectors are not equivalent (within
    ! an arbitrary tolerance)

    logical :: vnev
    real(kind=rd_kind), dimension(3), intent(in) :: vL, vR

    ! Just call the equivalence routine and return the opposite value
    vnev = (.not. (vL .veqv. vR))

  end function vnev

  function vnev_single(vL, vR) result(vnev)

    ! Determines if two vectors are not equivalent (within
    ! an arbitrary tolerance)

    logical :: vnev
    real(kind=rs_kind), dimension(3), intent(in) :: vL, vR

    ! Just call the equivalence routine and return the opposite value
    vnev = (.not. (vL .veqv. vR))

  end function vnev_single

  function vnev_single_double(vL, vR) result(vnev)

    ! Determines if two vectors are not equivalent (within
    ! an arbitrary tolerance)

    logical :: vnev
    real(kind=rs_kind), dimension(3), intent(in) :: vL
    real(kind=rd_kind), dimension(3), intent(in) :: vR

    ! Just call the equivalence routine and return the opposite value
    vnev = (.not. (real(vL,rd_kind) .veqv. vR))

  end function vnev_single_double

  function vnev_double_single(vL, vR) result(vnev)

    ! Determines if two vectors are not equivalent (within
    ! an arbitrary tolerance)

    logical :: vnev
    real(kind=rd_kind), dimension(3), intent(in) :: vL
    real(kind=rs_kind), dimension(3), intent(in) :: vR

    ! Just call the equivalence routine and return the opposite value
    vnev = (.not. (vL .veqv. real(vR,rd_kind)))

  end function vnev_double_single

  logical function plane_eq_plane(planeL, planeR)

    type (plane_object), intent(in) :: planeL, planeR

    ! Planes are equal if their normals and constants are equal
    ! within a given tolerance
    plane_eq_plane = ( (planeL%normal .veqv. planeR%normal) .and. & 
         (abs(planeL%constant - planeR%constant) < tolerance) )

  end function plane_eq_plane

  logical function plane_neq_plane(planeL, planeR)

    type (plane_object), intent(in) :: planeL, planeR

    ! Planes are equal if their normals and constants are equal
    ! within a given tolerance
    plane_neq_plane = ( .not. planeL == planeR )

  end function plane_neq_plane

  ! Access routines

  function get_normal ( plane )

    ! Access the normal vector in a plane object
    
    ! Return value
    real(kind=rd_kind), dimension(3) :: get_normal

    type (plane_object), intent(inout) :: plane

    get_normal = plane%normal

    return
  end function get_normal

  function set_normal ( plane, vector )

    ! Access the normal vector in a plane object
    
    ! Return value
    real(kind=rd_kind), dimension(3) :: set_normal

    type (plane_object), intent(inout)           :: plane
    real(kind=rd_kind), dimension(:), intent(in) :: vector

    call new(plane, vector, plane%constant)

    set_normal = plane%normal

    return
  end function set_normal

  function set_normal_single ( plane, vector ) result(set_normal)

    ! Access the normal vector in a plane object
    
    ! Return value
    real(kind=rd_kind), dimension(3) :: set_normal

    type (plane_object), intent(inout)           :: plane
    real(kind=rs_kind), dimension(:), intent(in) :: vector

    call new(plane, vector, plane%constant)

    set_normal = plane%normal

    return
  end function set_normal_single

  function get_constant ( plane )

    ! Access the constant in a plane object
    
    ! Return value
    real(kind=rd_kind) :: get_constant

    type (plane_object), intent(inout) :: plane

    get_constant = plane%constant

    return
  end function get_constant

  function set_constant ( plane, constant )

    ! Access the constant in a plane object
    
    ! Return value
    real(kind=rd_kind) :: set_constant

    type (plane_object), intent(inout) :: plane
    real(kind=rd_kind), intent(in)     :: constant

    call new(plane, plane%normal, constant)

    set_constant = plane%constant

    return
  end function set_constant

  function set_constant_single ( plane, constant ) result(set_constant)

    ! Access the constant in a plane object
    
    ! Return value
    real(kind=rd_kind) :: set_constant

    type (plane_object), intent(inout) :: plane
    real(kind=rs_kind), intent(in)     :: constant

    call new(plane, plane%normal, constant)

    set_constant = plane%constant

    return
  end function set_constant_single

end module vector_class
