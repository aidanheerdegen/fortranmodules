module array_functions

  implicit none

  private

  real :: tolerance = 1e-3 

  ! Function which returns the reverse an array
  interface reverse
     module procedure reverse_int, reverse_real, reverse_double, reverse_char, reverse_logical
  end interface

  ! Overload comparison operator so that it will work with arrays
  interface operator (.alleq.)
     module procedure int_eq_int, real_eq_real, double_eq_double, char_eq_char, logical_eq_logical
  end interface

  ! Overload comparison operator so that it will work with arrays
  interface operator (.allneq.)
     module procedure int_neq_int, real_neq_real, double_neq_double, logical_neq_logical, char_neq_char
  end interface

  ! Public routines
  public :: reverse

  ! Overloaded operators
  public :: operator(.alleq.), operator(.allneq.)

contains

  ! Reverse functions

  function reverse_int( array )

    ! Return the reverse of an array

    ! Input variables
    integer, dimension(:), intent(in) :: array

    ! Return value
    integer, dimension(size(array)) :: reverse_int

    reverse_int = array(ubound(array,1):lbound(array,1):-1)

  end function reverse_int

  function reverse_real( array )

    ! Return the reverse of an array

    ! Input variables
    real, dimension(:), intent(in) :: array

    ! Return value
    real, dimension(size(array)) :: reverse_real

    reverse_real = array(ubound(array,1):lbound(array,1):-1)

  end function reverse_real

  function reverse_double( array )

    ! Return the reverse of an array

    ! Input variables
    real(kind=8), dimension(:), intent(in) :: array

    ! Return value
    real(kind=8), dimension(size(array)) :: reverse_double

    reverse_double = array(ubound(array,1):lbound(array,1):-1)

  end function reverse_double

  function reverse_char( array )

    ! Return the reverse of an array

    ! Input variables
    character, dimension(:), intent(in) :: array

    ! Return value
    character, dimension(size(array)) :: reverse_char

    reverse_char = array(ubound(array,1):lbound(array,1):-1)
    
  end function reverse_char

  function reverse_logical( array )

    ! Return the reverse of an array

    ! Input variables
    logical, dimension(:), intent(in) :: array

    ! Return value
    logical, dimension(size(array)) :: reverse_logical

    reverse_logical = array(ubound(array,1):lbound(array,1):-1)
    
  end function reverse_logical

  ! Equivalence
  
  logical function int_eq_int(arrayL,arrayR)

    integer, dimension(:), intent(in) :: arrayL, arrayR

    int_eq_int = all(arrayL == arrayR)

  end function int_eq_int
  
  logical function int_neq_int(arrayL,arrayR)

    integer, dimension(:), intent(in) :: arrayL, arrayR

    int_neq_int = (.not. (arrayL .alleq. arrayR)) 

  end function int_neq_int
  
  logical function real_eq_real(arrayL,arrayR)

    real, dimension(:), intent(in) :: arrayL, arrayR

    real_eq_real = all(arrayL == arrayR)

  end function real_eq_real
  
  logical function real_neq_real(arrayL,arrayR)

    real, dimension(:), intent(in) :: arrayL, arrayR

    real_neq_real = (.not. (arrayL .alleq. arrayR)) 

  end function real_neq_real
  
  logical function double_eq_double(arrayL,arrayR)

    real(kind=8), dimension(:), intent(in) :: arrayL, arrayR

    double_eq_double = all(arrayL == arrayR)

  end function double_eq_double
  
  logical function double_neq_double(arrayL,arrayR)

    real(kind=8), dimension(:), intent(in) :: arrayL, arrayR

    double_neq_double = (.not.(arrayL .alleq. arrayR)) 

  end function double_neq_double
  
  
  logical function char_eq_char(arrayL,arrayR)

    character, dimension(:), intent(in) :: arrayL, arrayR

    char_eq_char = all(arrayL == arrayR)

  end function char_eq_char
  
  logical function char_neq_char(arrayL,arrayR)

    character, dimension(:), intent(in) :: arrayL, arrayR

    char_neq_char = (.not.(arrayL .alleq. arrayR)) 

  end function char_neq_char

  logical function logical_eq_logical(arrayL,arrayR)

    logical, dimension(:), intent(in) :: arrayL, arrayR

    logical_eq_logical = all(arrayL .eqv. arrayR)

  end function logical_eq_logical

  logical function logical_neq_logical(arrayL,arrayR)

    logical, dimension(:), intent(in) :: arrayL, arrayR
    integer, dimension(size(arrayL)) :: intmaskL, intmaskR

    logical_neq_logical = (.not.(arrayL .alleq. arrayR)) 

  end function logical_neq_logical

end module array_functions
