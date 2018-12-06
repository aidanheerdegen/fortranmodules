module zmatrix_class

  use vector_class, only: operator(.angle.), operator(.cross.), operator(.veqv.), operator(.dot.), normalise
  use variable_array, only: push
  use cartesian_class, only: connectivity_entry
  use file_functions, only: stderr, stdout, open
  use string_functions, only: getword, gettext, ucase
  use fundamental_constants, only: radian
  use mol2_class, only: mol2_object, coords, order, labels, name, comment, connectivity, atom_num
  use sort_functions, only: sort, operator(.sort.)

  implicit none

  private

  character(len=*), parameter :: version = "$Id: zmatrix_class.f90,v 1.11 2005/02/11 05:10:54 aidan Exp aidan $" 

  !! $Log: zmatrix_class.f90,v $
  !! Revision 1.11  2005/02/11 05:10:54  aidan
  !! About to make a major change to allocatable components, so am saving this
  !! intermediate version that is showing a fairdeal of difference between itself
  !! and v.1.10.
  !!
  !! Revision 1.10  2004/08/16 02:29:28  aidan
  !! Added the ability to construct a z-matrix from a template zmatrix object.
  !! Minimal error checking and options are available at this stage -- the
  !! cartesian coordinates supplied with the template must be in the same order
  !! as they appear in the template. A wrapper for a mol2 object is also supplied.
  !! The workhorse for these functions is a new routine called calculate_internals
  !! which does what it's name suggests. I toyed with the idea of using this
  !! routine to calculate the values of the internal degrees of freedom for the
  !! makeint routine, but decided against it. In determining the best atoms with
  !! which to define a bond angle makeint has to determine the bond angle, so it
  !! seemed silly not to do the trivial calculations for bond length and dihedral
  !! angle as well.
  !!
  !! Re-wrote the reorder_zmatrix routine -- it was stupidly inefficient.
  !!
  !! Fixed a bug which was tripped in the delete_zmatrix (which itself was called
  !! from assign_from_zmatrix) which queried the association status of an undefined
  !! pointer (the zmatrices component of a z-matrix). This component is now
  !! initialised to null.
  !!
  !! Added another variable to the zmatrix_object -- "ordered". This is a logical
  !! which tells us if the z-matrix has been reordered. Not that useful currently
  !! as the only place this needs to be done is in makeint. May be useful down the
  !! track.
  !!
  !! Added a global debug flag. It is a parameter which means all statements that
  !! look like this:
  !!
  !!     if (debug) print *,"something"
  !!
  !! are optimised away in the compilation. So almost as good as preprocessing in
  !! that debugging flags are not included in the executable. I confirmed this by
  !! producing assembly code from the compiler with the flag set to .TRUE. and
  !! .FALSE..
  !!
  !! Revision 1.9  2004/07/16 03:52:33  aidan
  !! A *major* revision because I neglected to check it in before now. There is a
  !! big improvement in the z-matrix construction algorithm (makeint) which will
  !! hopefully remedy alot of the "duplication of torsions" that was occurring in
  !! the benzil type structures. This new algorithm attempts to make more local
  !! "pair of triangles" type dihedral angles in preference to the classic linear
  !! "torsion angle" sort of dihedrals.
  !!
  !! Have added a "set_zmatrix_single_param" function to the generic parameter
  !! interface. This allows a single parameter to set. It is just a wrapper to
  !! the general interface, but is necessary to avoid mucky confusing user code,
  !! and is pitifully easy to implement.
  !!
  !! There is a new function "zmatrix_from_template" which was added to the
  !! as_zmatrix interface, but has been removed. This does not do what it patently
  !! should, i.e. use the connectivity from an existing zmatrix to create a new
  !! z-matrix using, possibly different, cartesian coordinates. Instead the
  !! connectivity was collected together and blindly passed to the z-matrix
  !! creation algorithm. This did not produce the same result. Probably not
  !! surprising. This function should be re-written as it is potentially useful.
  !!
  !! Revision 1.8  2004/01/11 22:51:08  aidan
  !! Removed automatic assignment between z-matrix objects -- the constructor
  !! function was subject to this assignment also -- causing a memory leak! Doh!
  !! Made a new subroutine interface called 'copy' to accomplish a deep copy.
  !! Also had to check and ignore the return value of the deallocation as the
  !! pointer was (apparently) not created with 'allocate'. Odd.
  !!
  !! Revision 1.7  2004/01/09 05:20:41  aidan
  !! Added a routine to free the memory associated with a z-matrix object so as
  !! to avoid memory leaks.
  !!
  !! Revision 1.6  2004/01/09 03:40:02  aidan
  !! Added a deep copy  routine for assignment between two z-matrix objects. This
  !! fixed a bug where altering the contents of a local copy of the z-matrix also
  !! altered the contents of the original z-matrix. This was not a good thing in
  !! the context of monte carlo type operations.
  !!
  !! Revision 1.5  2003/08/21 05:26:19  aidan
  !! Updated the mol2->zmatrix conversion routine to take an optional substructure
  !! argument and reflect the changes to the mol2 access routines, i.e. no longer
  !! make temporary arrays and fill them with coords etc, just call the z-matrix
  !! creation routine with the function calls as parameters.
  !!
  !! Revision 1.4  2003/08/15 07:01:00  aidan
  !! Added a method for accessing the connectivity in a z-matrix (i.e. the
  !! atoms that define the bonds and angles within a z-matrix).
  !!
  !! Added RCS identifying strings
  !!
  !! Cleaned out old code in the parameter setting subroutine
  !!
  !! Deallocated memory in the mol2_to_zmatrix routine
  !!

  type zmatrix_entry
     private
     integer :: conn(3), number
     real :: bond_length, bond_angle, dihedral_angle
     character(len=6) :: label
  end type zmatrix_entry

  type zmatrix_object
     private
     type(zmatrix_entry), dimension(:), allocatable :: zmatrices
     character(len=180) :: info
     logical :: ordered = .FALSE., allocated = .FALSE.
  end type zmatrix_object

  ! Define a tolerance for equivalence and "effectively zero"
  real :: tolerance = 5e-5 

  ! Make a debugging flag
  logical, parameter :: debug = .FALSE.

  ! Procedure to load a zmatrix object
  interface load_zmatrix
     module procedure open_then_read, readzmatrix
  end interface

  ! Procedure to load a zmatrix object
  interface delete
     module procedure delete_zmatrix
  end interface

  ! Make a generic interface for bond distance
  interface distance
     module procedure bond_distance
  end interface

  ! Make a generic interface for angle calculations
  interface angle
     module procedure bond_angle, dihedral
  end interface

  ! Make a generic interface to zmatrix construction function
  interface as_zmatrix
     module procedure makeint, mol2_to_zmatrix, zmatrix_from_template, zmatrix_from_template_mol2
  end interface

  ! Make a generic interface to xyz construction function
  interface as_xyz
     module procedure makexyz_rasmol
  end interface

  ! Make a generic interface to xyz construction function
  interface as_cartesian
     module procedure makexyz_rasmol
  end interface

  ! Make a generic interface to zmatrix print routine
  interface print
     module procedure print_zmatrix, open_then_print, print_zmatrix_to_stdout
  end interface

  ! Make a generic interface to setting zmatrix parameters
  interface parameter
     module procedure set_zmatrix_param_general, set_zmatrix_single_param
  end interface

  ! Make a generic interface to setting zmatrix connectivity
  interface connectivity
     module procedure set_zmatrix_conn_general
  end interface

  ! Make a generic interface to setting zmatrix information
  interface info
     module procedure set_zmatrix_info
  end interface

  ! Make a generic interface to the inquiry function for finding out the number of atoms
  interface num
     module procedure get_num_entries
  end interface

  ! Function to set and get the atoms labels from a zmatrix object
  interface labels
     module procedure set_atom_labels, set_one_atom_label, get_all_atom_labels
  end interface

  interface copy
     module procedure assign_from_zmatrix
  end interface

  ! Define zmatrix equivalence operator (also defines the .eq. operator by default)
  interface operator (==)
     module procedure zmat_eq_zmat
  end interface

  ! Define zmatrix equivalence operator (also defines the .ne. operator by default)
  interface operator (/=)
     module procedure zmat_neq_zmat
  end interface

  ! Define assignment from a filename
  interface assignment (=)
     module procedure assign_from_filename, assign_from_mol2
  end interface

  !!!!!!!!!!!!!!!!
  ! Public stuff !
  !!!!!!!!!!!!!!!!

  ! Public data types
  public :: zmatrix_object, zmatrix_entry

  ! Public routines
  public :: angle, distance, as_zmatrix, print, as_xyz, parameter, delete, copy
  public :: as_cartesian, load_zmatrix, info, num, labels, connectivity

  ! overloaded operators
  public :: operator(==), operator(/=), assignment(=)

contains

  subroutine delete_zmatrix(zmatrix)

    ! Free the memory in a z-matrix object
  
    type (zmatrix_object), intent(inout) :: zmatrix

    integer :: error

    ! print *,'in delete'
    ! deallocate(zmatrix%zmatrices, stat=error)
    ! if (associated(zmatrix%zmatrices)) deallocate(zmatrix%zmatrices, stat=error)
    ! if (zmatrix%allocated) deallocate(zmatrix%zmatrices, stat=error)
    ! deallocate(zmatrix%zmatrices,stat=error)
    if (allocated(zmatrix%zmatrices)) deallocate(zmatrix%zmatrices)
    ! print *,'here'
    ! nullify(zmatrix%zmatrices)

    ! Ignore status of deallocation ... must fix this
    ! print *,'delete ok'

  end subroutine delete_zmatrix

  ! Accessing zmatrix files ...

  subroutine assign_from_mol2(zmatrix,mol2)

    ! Another wrapper to allow direct assignment of a zmatrix 
    ! object from a filename
  
    ! LHS of assignment
    type (zmatrix_object), intent(inout) :: zmatrix
    ! RHS of assignment
    type (mol2_object), intent(in) :: mol2

    call delete(zmatrix)

    zmatrix = as_zmatrix(mol2)

  end subroutine assign_from_mol2

  subroutine assign_from_zmatrix(zmatrixl,zmatrixr)

    ! Do a deep copy of data from one zmatrix object to another
  
    ! LHS of assignment
    type (zmatrix_object), intent(inout) :: zmatrixl
    ! RHS of assignment
    type (zmatrix_object), intent(in)  :: zmatrixr

    ! Local variables
    integer :: num_entries, i

    call delete(zmatrixl)

    num_entries = size(zmatrixr%zmatrices)
    allocate(zmatrixl%zmatrices(num_entries))
    zmatrixl%allocated = .TRUE.

    do i = 1, num_entries
       zmatrixl%zmatrices(i) = zmatrixr%zmatrices(i)
    end do

    zmatrixl%info = zmatrixr%info
    zmatrixl%ordered = zmatrixr%ordered

  end subroutine assign_from_zmatrix

  subroutine assign_from_filename(zmatrix,file)

    ! Another wrapper to allow direct assignment of a zmatrix 
    ! object from a filename
  
    type (zmatrix_object), intent(out) :: zmatrix
    ! zmatrix filename
    character(len=*), intent(in) :: file

    zmatrix = open_then_read(file)

  end subroutine assign_from_filename

  function open_then_read(file)

    ! This is a small wrapper that allows the user to specify the
    ! name of a zmatrix file which this routine will open for them, 
    ! call the zmatrix reading routine and then close the file again.

    ! Return type of function
    type (zmatrix_object) :: open_then_read
    
    ! zmatrix filename
    character(len=*), intent(in) :: file

    ! Internal variables
    logical :: opened_ok
    integer :: unit

    ! Open the zmatrix file
    unit = open(file, opened_ok, status='old')

    if (opened_ok) then
       open_then_read = readzmatrix(unit)
       close(unit)
    else
       write (stderr,*)'ZMATRIX_CLASS :: Unable to find file: ',file
       stop
    end if

  end function open_then_read

  function readzmatrix (unit)

    ! Main zmatrix parsing function. It will read from an already opened
    ! unit number and store the results in a zmatrix object which it returns

    ! Return type of function
    type(zmatrix_object) :: readzmatrix

    ! Unit number
    integer, intent(in) :: unit

    ! Local variables
    integer :: line_number, entry_number, num_entries, error, next, id(3)
    real(kind=8) :: value(3)
    character(len=200) :: record
    character(len=6) :: label

    logical :: reading_entries
    
10  format(a200)

    line_number = 0
    entry_number = 0
    id = 0
    value = 0
    reading_entries = .FALSE.

    do
       ! Read through the zmatrix file until we reach the end 
       read (unit,10,end=1)  record

       if (debug) print *,record
       record = adjustl(record)

       ! Ignore comment lines
       if(record(1:1) .eq. '#') cycle 

       line_number = line_number + 1

       ! Check if we read the last entry last time through the loop
       if (entry_number /= 0 .and. (entry_number == num_entries)) exit 

       ! Use the (non-comment) line number to choose what action
       ! we will take to parse the line
       select case (line_number)
       case (1) 

          ! the first line is the comment string
          readzmatrix%info = trim(record(:len(readzmatrix%info))) 

       case (2) 

          ! The second line is the number of entries in the zmatrix
          read (record,*)  num_entries
          ! Allocate memory for the z-matrix entries
          ! nullify(readzmatrix%zmatrices)
          allocate(readzmatrix%zmatrices(num_entries), stat=error)
          readzmatrix%allocated = .TRUE.
          if (error /= 0) then
             write(stderr,*) 'ZMATRIX_CLASS :: readzmatrix :: Cannot allocate memory for zmatrix'
             stop
          end if

       case (3) 

          reading_entries = .TRUE.

          ! This is the first entry in the z-matrix
          entry_number = 1

          ! Grab the atom label -- may or may not contain more than just an atom 
          ! label, but we will ignore it all anyway ....
          next = 1
          label = trim(gettext(record,next))

       case (4) 

          ! This is the second entry in the z-matrix
          entry_number = 2

          ! Grab the atom label ...
          next = 1
          label = trim(gettext(record,next))

          ! ... and the id and bond length
          read (record(next:),*) id(1), value(1)

       case (5) 

          ! This is the third entry in the z-matrix
          entry_number = 3

          ! Grab the atom label ...
          next = 1
          label = trim(gettext(record,next))

          ! ... and the id and bond length and another id and bond angle
          read (record(next:),*) id(1), value(1), id(2), value(2)

       case default

          ! All subsequent entries
          entry_number = entry_number + 1

          ! Grab the atom label ...
          next = 1
          label = trim(gettext(record,next))

          ! ... and the id and bond length and another id and bond angle, and another id 
          ! and dihedral angle
          read (record(next:),*) id(1), value(1), id(2), value(2), id(3), value(3)

       end select

       ! Save z-matrix entry to object if we are reading entries
       if (reading_entries) &
       readzmatrix%zmatrices(entry_number) = zmatrix_entry( id, entry_number, value(1), value(2)/radian, value(3)/radian, label )
    
    end do
1   return

  end function readzmatrix

  real function bond_distance(p1,p2)

    ! Two points in real space which we will use to determine the 
    ! euclidian distance between
    real, dimension(3), intent(in) :: p1, p2

    bond_distance = sqrt(sum((p2 - p1)**2))

  end function bond_distance

  function bond_angle(p1,p2,p3)

    ! Three points in real space which we will use to determine the bond
    ! angle. 
    real, dimension(3), intent(in) :: p1, p2, p3

    real :: bond_angle

    ! Internal variables
    real, dimension(3) :: a, b

    ! Make vectors out of the three points .. both pointing from
    ! the central point (p2) to the other points (p1 and p3)
    a = p1 - p2
    call normalise(a)
    b = p3 - p2
    call normalise(b)

    bond_angle = a .angle. b

  end function bond_angle

  function dihedral(p1,p2,p3,p4)

    ! Given four points in cartesian space, the dihedral angle is 
    ! defined as the angle between the plane formed by p1, p2 and
    ! p3 and the plane formed by p2, p3 and p4. The dihedral angle
    ! is also defined as the angle between the unit normals of these
    ! two planes. The normal to a plane is the cross product of two 
    ! vectors in the plane. 
    !
    ! So our strategy is to form three unit vectors (a,b,c) traversing 
    ! the points from p1 to p4. Then we take the cross products of 
    ! these vectors in turn. This produces two unit normal vectors.
    ! The dihedral angle is the angular separation of these two unit
    ! normals.

    ! Four points in real space 
    real, dimension(3), intent(in) :: p1, p2, p3, p4

    real :: dihedral

    ! Internal variables
    real, dimension(3) :: a, b, c, norm1, norm2
    real :: sign

    ! Make the three unit vectors. Note that they are the vectors
    ! *from* p1 *to* p2 and *from* p2 *to* p3 and *from* p3 *to* p4.
    a = p2 - p1
    call normalise(a)
    b = p3 - p2
    call normalise(b)
    c = p4 - p3
    call normalise(c)

    ! The two unit normals are the cross products of the unit vectors
    norm1 = a .cross. b
    norm2 = b .cross. c

    ! If either of our cross products produces the zero vector then we
    ! will return zero
    if ( (norm1 .veqv. (/ 0., 0., 0./)) .or.   &
         (norm2 .veqv. (/ 0., 0., 0./)) ) then
       dihedral = 0.0
       return
    end if

    ! The dihedral angle is the angular separation of the two unit normals
    dihedral = norm1 .angle. norm2

    ! The sign of the dihedral angle is given by:
    sign = a .dot. (b .cross. c)

    ! Reverse the sign of the dihedral angle if needed
    if (sign < 0) then
       ! Sign is positive, should be negative
       if (dihedral > 0) dihedral = -dihedral
    else
       ! Sign is negative, should be positive
       if (dihedral < 0) dihedral = -dihedral
    end if

  end function dihedral

  function dihedral_array(array_of_points) result(dihedral)

    ! This is just a wrapper to the dihedral function but allows input of
    ! the points in a single 3X4 matrix. This is a convenient capability
    ! as we can just throw it a subset of a cartesian coordinate array

    ! Four points in real space 
    real, dimension(3,4), intent(in) :: array_of_points

    real :: dihedral

    dihedral = angle(array_of_points(:,1),array_of_points(:,2),array_of_points(:,3),array_of_points(:,4))

  end function dihedral_array

  function mol2_to_zmatrix(mol2, substructure)

    ! Wrapper function to makeint -- extracts all the necessaries out
    ! of the mol2 object and makes a z-matrix using makeint. Convenience
    ! routine in a sense -- but also saves someone from having to understand
    ! those complicated variable declarations (see local variables below)
    
    ! Return type of function
    type(zmatrix_object) :: mol2_to_zmatrix

    ! Input variable
    type (mol2_object), intent(in) :: mol2
    integer, intent(in), optional  :: substructure

    ! Local variables
    character(len=80)  :: mol2name, mol2comment
    character(len=180) :: info_string
    integer :: subst_id

    subst_id = 1

    if (present(substructure)) subst_id = substructure

    ! This was required to make the test program compile. I know .. seriously odd.
    ! write(*,'(A)',advance='no') ''

    ! if (associated(mol2_to_zmatrix%zmatrices)) print *,'Houston, we have a problem'

    ! call delete(mol2_to_zmatrix)

    ! Feed z-matrix making routine the coordinates, the connectivity, the ordering 
    ! and the labels for the requested substructure ...
    mol2_to_zmatrix = as_zmatrix(coords(mol2,subst_id), &
                                 connectivity(mol2,subst_id), &
                                 order(mol2,subst_id), &
                                 labels(mol2,subst_id))

    ! See if we can grab some useful info from the mol2 file
    mol2name = name(mol2)
    mol2comment = comment(mol2)

    if (trim(mol2name) .ne. '') then
       if (trim(mol2comment) .ne. '') then
          info_string = trim(mol2name)//' : '//trim(mol2comment)
       else
          info_string = mol2name
       end if
    else
       info_string = mol2comment
    end if
    
    ! Save our information string into our z-matrix object
    info_string = info(mol2_to_zmatrix,info_string) 

  end function mol2_to_zmatrix

  function zmatrix_from_template_mol2(template, mol2, substructure) result(zmatrix)

    ! Wrapper function to zmatrix_from_template
    
    ! Return type of function
    type(zmatrix_object) :: zmatrix

    ! Input variable
    type (mol2_object), intent(in)   :: mol2
    type(zmatrix_object), intent(in) :: template
    integer, intent(in), optional    :: substructure

    ! Local variables
    character(len=80)  :: mol2name, mol2comment
    character(len=180) :: info_string
    integer :: subst_id

    subst_id = 1
    if (present(substructure)) subst_id = substructure

    ! Check we have enough coordinates in the specified substructure
    if (atom_num(mol2,subst_id) /= num(template)) then
       write(stderr,&
            "('ZMATRIX_CLASS :: Number of atoms in zmatrix template (',I0,') does not match number those in mol2 (',I0,')')")&
            num(template),atom_num(mol2,subst_id)
    end if

    ! Invoke the generic routine
    zmatrix = as_zmatrix(template, coords(mol2,subst_id))

    ! See if we can grab some useful info from the mol2 file
    mol2name = name(mol2)
    mol2comment = comment(mol2)

    if (trim(mol2name) .ne. '') then
       if (trim(mol2comment) .ne. '') then
          info_string = trim(mol2name)//' : '//trim(mol2comment)
       else
          info_string = mol2name
       end if
    else
       info_string = mol2comment
    end if
    
    ! Save our information string into our z-matrix object
    info_string = info(zmatrix,info_string) 

  end function zmatrix_from_template_mol2

  function zmatrix_from_template(template, coordinates) result(zmatrix)

    ! Use an existing z-matrix as a template for a new z-matrix -- the
    ! bond lengths, angles and dihedrals are recalculated based on the
    ! coordinates supplied and the connectivity in the template
    
    ! Return type of function
    type(zmatrix_object) :: zmatrix

    ! Input variable
    type(zmatrix_object), intent(in) :: template
    real, dimension(:,:), intent(in) :: coordinates

    if (size(coordinates,2) /= num(template)) then
       write(stderr,&
    "('ZMATRIX_CLASS :: Number of atoms in zmatrix template (',I0,') does not match number those in coordinate array (',I0,')')")&
       num(template),size(coordinates,2)
       ! stop
    end if

    ! Copy contents of template in zmatrix
    call copy(zmatrix,template)

    ! Calculate all the bond lengths, angles and dihedrals
    call calculate_internals(zmatrix, coordinates)

  end function zmatrix_from_template

  function makeint(cartesian, connectivity, order, labels)
    
    ! "makeint" converts Cartesian to internal coordinates and returns
    ! a pointer to a zmatrix "object" (actually an array of zmatrix_entries)
    !
    ! Note: The entries in the z-matrix are indexed by their position in
    ! the list of cartesian coordinates. This is so we can rely on the
    ! connectivity listing to be accurate. We note down the actual position
    ! in which each entry is placed in the z-matrix (which can be altered
    ! by the optional parameter 'order') and then sort the z-matrix into 
    ! this order once we have constructed it.

    ! Note also that we don't use the 'calculate internals' routine to 
    ! determine the bond lengths, angles and dihedrals as we use some feedback
    ! from the bond angle in choosing suitable candidates to connect to. It
    ! seemed sort of pointless not to calculate the other two values as well.

    ! Return type of function (is a pointer, so need to use pointer
    ! assignment, e.g. zmatrix <= makeint( ... )
    type(zmatrix_object) :: makeint

    ! Input variables

    ! Cartesian is a 3 x num_atoms matrix of cartesian coordinates
    real, dimension(:,:), intent(in)                   :: cartesian 

    ! Connectivity is a 'ragged array' of connections between the atoms
    ! specified in cartesian (should be length num_atoms)
    type(connectivity_entry), dimension(:), intent(in) :: connectivity

    ! An optional ordering array that specifies a particular ordering of the atoms
    integer, dimension(:), optional, intent(in) :: order

    ! An optional label array that specifies labels for the atoms
    character(len=*), dimension(:), optional, intent(in) :: labels

    ! Local variables
    integer :: i, error, num_atoms, base_atom, second_atom, third_atom, fourth_atom, dummy
    logical :: have_order, have_labels
    real :: bond_angle

    integer, dimension(:), pointer :: skip_list
    
    nullify(skip_list)

    second_atom = 0; third_atom = 0; fourth_atom = 0

    num_atoms = size(cartesian,2)

    ! if (associated(makeint%zmatrices)) print *,'Houston, we have a problem'

    ! call delete(makeint)

    ! Allocate memory for the z-matrix entries
    ! nullify(makeint%zmatrices)
    allocate(makeint%zmatrices(num_atoms), stat=error)
    makeint%allocated = .TRUE.

    if (error /= 0) then
       write(stderr,*) 'ZMATRIX_CLASS :: Cannot allocate memory for zmatrix'
       stop
    end if

    ! Initialise the info to something useful
    makeint%info = ''

    have_labels = present(labels)
    have_order = present(order)

    ! Intialise the zmatrix entries
    do i = 1, num_atoms
       makeint%zmatrices(i) = zmatrix_entry( (/0,0,0/), 0, 0.0, 0.0, 0.0, '' )
    end do

    ! Cycle through the array of atoms
    do i = 1, num_atoms

       ! Default ordering goes from 1 .. num_atoms
       base_atom = i

       ! However, if we have provided an ordering array, then use the
       ! i'th entry from this
       if (have_order) base_atom = order(i)

       if (have_labels) makeint%zmatrices(base_atom)%label = labels(base_atom)

       if (debug) print *,i,base_atom

       if (i >= 2) then
          ! Define a bond length for the base atom
          
          if (debug) print *,'Looking for bond length'
          ! Find a second atom to bond to
          second_atom = adjacent (base_atom, (/ 0 /), makeint, connectivity)

          ! If we fail to find any suitable atoms connected to the base
          ! atom then tell the punter about it and stop
          if (second_atom .eq. 0) then 
             call squeal(base_atom, i)
             write(stderr,"('Error occurred attempting to find second atom to define bond length')")
             write(stderr,'("Current zmatrix : ")')
             call print(stderr,makeint)
             stop
          end if

          ! Calculate the bond length
          makeint%zmatrices(base_atom)%bond_length =    &
               distance(cartesian(:, base_atom),cartesian(:, second_atom))
          
          if (i >= 3) then
             ! define a bond angle for the base atom
             
             nullify(skip_list)
             dummy = push(skip_list, base_atom)

             ! Cycle through all the eligible adjacent atoms until we find one that
             ! is not linear
             do
             
                if (debug) print *,'Looking for an atom to form a bond angle'
                ! Find a third atom so we can generate a bond angle
                third_atom = adjacent (second_atom, skip_list, makeint, connectivity)
                
                ! If we fail to find any suitable atoms connected to the second
                ! atom then tell the punter about it and stop
                if (third_atom .eq. 0) then 
                   call squeal(base_atom, i)
                   write(stderr,"('Error occurred attempting to find third atom to define bond angle')")
                   write(stderr,'("second atom: ",I4)') second_atom
                   write(stderr,'("Current zmatrix : ")')
                   call print(stderr,makeint)
                   stop
                end if
                
                ! Calculate the bond angle
                bond_angle = angle(cartesian(:,third_atom), &
                     cartesian(:,second_atom), cartesian(:,base_atom))

                ! Accept this bond angle if it is not close to 180 deg
                if (abs(bond_angle*radian - 180.) > 5.) exit

                ! Push this atom on to the skip list and get another contender
                dummy = push(skip_list, third_atom)

             end do

             makeint%zmatrices(base_atom)%bond_angle = bond_angle
             
             if (i >= 4) then
                ! define the dihedral angle for the current atom

                if (debug) print *,'Looking for an atom to form a dihedral angle'
                ! Find a fourth atom so we can generate a bond angle
                ! We will look for a candidate atom bonded to the second atom,
                ! as this tends to produce a better outcome with bridged molecules
                ! like benzil -- producing more 'local' z-matrices.
                fourth_atom = adjacent (second_atom, (/ base_atom, third_atom /), makeint, connectivity)

                ! We failed to find any suitable atoms connected to the second atom
                if (fourth_atom .eq. 0) then 

                   ! Try and find a fourth atom connected to the third atom
                   fourth_atom = adjacent (third_atom, (/ second_atom /), makeint, connectivity)
                
                   ! If we fail to find any suitable atoms connected to the third
                   ! or second atom then tell the punter about it and stop
                   if (fourth_atom .eq. 0) then 
                      call squeal(base_atom, i)
                      write(stderr,"('Error occurred attempting to find fourth atom to define dihedral angle')")
                      write(stderr,'("second atom: ",I4)') second_atom
                      write(stderr,'("third atom: ",I4)') third_atom
                      write(stderr,'("Current zmatrix : ")')
                      call print(stderr,makeint)
                      stop
                   end if

                end if

                ! Calculate the dihedral angle
                makeint%zmatrices(base_atom)%dihedral_angle = angle(cartesian(:,base_atom),  &
                     cartesian(:,second_atom), cartesian(:,third_atom),  &
                     cartesian(:,fourth_atom))
             end if
          end if
       end if

       ! Initialise the connections array
       makeint%zmatrices(base_atom)%conn = 0

       if (debug) print *,i,base_atom, second_atom, third_atom, fourth_atom

       ! transfer defining atoms to permanent array if they are > 0
       if (second_atom > 0) makeint%zmatrices(base_atom)%conn(1) = makeint%zmatrices(second_atom)%number
       if (third_atom  > 0) makeint%zmatrices(base_atom)%conn(2) = makeint%zmatrices(third_atom)%number
       if (fourth_atom > 0) makeint%zmatrices(base_atom)%conn(3) = makeint%zmatrices(fourth_atom)%number
         
       ! mark the current base atom as finished by adding a number entry for it
       makeint%zmatrices(base_atom)%number = i

    end do

    ! Put the z-matrix into the correct order
    call reorder_zmatrix(makeint)

    return
    
  end function makeint

  subroutine calculate_internals(zmatrix, cartesian)
    
    ! Given the skeleton of a z-matrix, i.e. all the connections
    ! definining the atoms, this routine "fills in the blanks" of
    ! the bond lengths, angles and dihedral angles.

    ! Interface variables

    ! This is the z-matrix "template" that we will fill in
    type(zmatrix_object), intent(inout) :: zmatrix

    ! Cartesian is a 3 x num_atoms matrix of cartesian coordinates
    ! in the same order as the z-matrix entries
    real, dimension(:,:), intent(in)                   :: cartesian 

    ! Local variables
    integer :: i, conn(3)
    integer :: base_atom, second_atom, third_atom, fourth_atom

    ! Intialise the first three zmatrix entries
    do i = 1, 3
       zmatrix%zmatrices(i)%bond_length = 0.0
       zmatrix%zmatrices(i)%bond_angle = 0.0
       zmatrix%zmatrices(i)%dihedral_angle = 0.0
    end do

    ! Cycle through the entries in zmatrix
    do base_atom = 2, num(zmatrix)

       ! conn = connectivity(zmatrix,spread(base_atom,3,1),(/1,2,3/))
       conn = connectivity(zmatrix,(/base_atom,base_atom,base_atom/),(/1,2,3/))

       if (debug) print *,base_atom,zmatrix%zmatrices(base_atom)%number,conn

       ! Grab the second atom as defined in the z-matrix
       second_atom = conn(1)

       ! Calculate the bond length
       zmatrix%zmatrices(base_atom)%bond_length =    &
               distance(cartesian(:, base_atom),cartesian(:, second_atom))
          
       if (base_atom < 3) cycle

       ! define a bond angle for the base atom
       
       ! Grab the third atom so we can generate a bond angle
       third_atom = conn(2)
                
       ! Calculate the bond angle
       zmatrix%zmatrices(base_atom)%bond_angle = angle(cartesian(:,third_atom), &
            cartesian(:,second_atom), cartesian(:,base_atom))
          
       if (base_atom < 4) cycle

       ! define the dihedral angle for the current atom
       
       ! Grab the fourth atom so we can generate a dihedral angle
       fourth_atom = conn(3)

       ! Calculate the dihedral angle
       zmatrix%zmatrices(base_atom)%dihedral_angle = angle(cartesian(:,base_atom),  &
                     cartesian(:,second_atom), cartesian(:,third_atom),  &
                     cartesian(:,fourth_atom))

    end do

    return
    
  end subroutine calculate_internals
    
  function adjacent (base_atom_id, skip, zmatrix, connectivity)
      
    ! "adjacent" finds an atom connected to our "base_atom_id" 
    ! other than any atoms specified in the array "skip". If no 
    ! such atom exists, then a zero is returned
    
    integer :: adjacent
    
    integer, intent(in) :: base_atom_id, skip(:)
    type (zmatrix_object), intent(in) :: zmatrix
    type (connectivity_entry), dimension(:), intent(in) :: connectivity
    
    integer atom_id, j, k, num_valid
    integer, dimension(:), pointer ::  valid_connections
    double precision ::  dist, short
    
    logical :: defining_adjacent

    ! Initialise the list of eligible atoms bonded to the atom of interest
    nullify(valid_connections)
    num_valid = 0
    
    if (debug) print *,'Base Atom id:',base_atom_id
    ! Cycle through all the 1-2 bonds for atom in reverse order.
    ! We reverse the order so that the z-matrix is as 'local'
    ! as possible, i.e. defined by other atoms that are as close
    ! as possible in the z-matrix order. In this way changing a 
    ! torision angle will alter the whole structure, whereas 
    ! connections which 'reproduce' previous torsion angles would 
    ! stop this from happening.
    ! Geddit? Good, because this was a terrible explanation!
    do j = size(connectivity(base_atom_id)%connectivity), 1, -1

       atom_id = connectivity(base_atom_id)%connectivity(j)
       if (debug) print *,'Atom id:',atom_id

       ! Check if atom i has already been defined in the zmatrix and 
       ! is not atom skip
       if (zmatrix%zmatrices(atom_id)%number /= 0 .and. all(atom_id /= skip)) then
          ! if (skip == 0) then
             ! Push the atom id onto the array keeping track of valid connections 
             num_valid = push(valid_connections, atom_id)
          ! else
             ! If atom id is already defined as being bonded to atom
             ! or vice-versa then we accept it
             ! NB: This might be an excessively stringent criterion for
             ! acceptance -- we might only care if the two atoms are 
             ! connected at all, not just connected in the z-matrix sense?
             ! if (zmatrix%zmatrices(atom_id)%conn(1) == base_atom_id .or.    &
             !     zmatrix%zmatrices(base_atom_id)%conn(1) == atom_id .or.    &
             !     any(connectivity(base_atom_id)%connectivity == atom_id) &
             !    ) then
             !    num_valid = push(valid_connections, atom_id)
             ! end if
          ! end if
       end if

    end do
    
    ! if no bonded atom is eligible, use the nearest neighbor? Prolly best
    ! to return zero and search for nearest neighbour in routine above
    
    if (num_valid .eq. 0) then
       adjacent = 0
       ! short = 1000000.0d0
       ! do i = 1, n
       !    if (iz0(i).ne.0 .and. i.ne.atom .and. i.ne.skip) then
       !       dist = sqrt((x(i)-x(atom))**2 + (y(i)-y(atom))**2         &
       !            &                              + (z(i)-z(atom))**2)
       !       if (dist .lt. short) then
       !          short = dist
       !          adjacent = i
       !       end if
       !    end if
       ! end do
       ! if (skip .eq. 0) then
       !    ndel = ndel + 1
       !    idel(1,ndel) = adjacent
       !    idel(2,ndel) = atom
       ! end if
    else

       ! use an adjacent atom bonded to undefined atoms
       adjacent = valid_connections(1)

       return

       ! cycle through all the valid connections
       loop: do k = 1, num_valid
          ! cycle through all the atoms connected to each of
          ! the valid connections
          do j = 1, size(connectivity(valid_connections(k))%connectivity)
             atom_id = connectivity(valid_connections(k))%connectivity(j)
             ! if we find one that is in the z-matrix and not
             ! our base atom, then we'll accept it and stop
             ! searching
             if (zmatrix%zmatrices(atom_id)%number /= 0 .and. atom_id /= base_atom_id) then
                adjacent = valid_connections(k)
                exit loop
             end if
          end do
       end do loop

    end if
    
    return
    
  end function adjacent

  subroutine squeal(atom_id, atom_num)

    ! A small semi-standard squeal to stderr

    integer, intent(in) :: atom_id, atom_num

1   format ("zmatrix_class :: Connectivity error -- atom id: ",I4,"   atom number: ",I4)
    write(stderr,1) atom_id, atom_num

  end subroutine squeal

  subroutine open_then_print(file, zmatrix)

    ! This is a small wrapper that allows the user to specify a filename
    ! to print to rather than a unit number. This routine will open the
    ! file for them, print the zmatrix and then close it again.

    type (zmatrix_object), intent(in) :: zmatrix

    ! zmatrix filename
    character(len=*), intent(in) :: file

    ! Internal variables
    logical :: opened_ok
    integer :: unit

    ! Open the zmatrix file
    unit = open(file, opened_ok)

    if (opened_ok) then
       call print(unit,zmatrix)
       close(unit)
    else
       write (stderr,*)'ZMATRIX_CLASS :: Unable to open file: ',file
       stop
    end if

  end subroutine open_then_print

  subroutine print_zmatrix (unit, zmatrix)

    ! Output a zmatrix in a standard format

    type (zmatrix_object), intent(in) :: zmatrix

    integer, intent(in) :: unit

    ! Local variables
    integer :: i
    character(len=6) :: label
    character(len=180) :: info_string

    ! Make sure we have a properly ordered zmatrix before we print it out
    ! call reorder_zmatrix(zmatrix)

    info_string = zmatrix%info
    if (trim(info_string) .eq. '') info_string = 'Written from ZMATRIX_CLASS'

    write(unit,'(A)') trim(info_string)
    write(unit,'(I0)') size(zmatrix%zmatrices)
    do i=1,size(zmatrix%zmatrices)
       if (zmatrix%zmatrices(i)%label .eq. '') then
          write(label,'(I6)') zmatrix%zmatrices(i)%number
       else
          label = zmatrix%zmatrices(i)%label 
       end if
       if (zmatrix%zmatrices(i)%number > 0) &
            write(unit,'(A,I5,F10.5,2(I5,F12.5))') label, &
            zmatrix%zmatrices(i)%conn(1),                 &
            zmatrix%zmatrices(i)%bond_length,             & 
            zmatrix%zmatrices(i)%conn(2),                 &
            zmatrix%zmatrices(i)%bond_angle*radian,       & 
            zmatrix%zmatrices(i)%conn(3),                 &
            zmatrix%zmatrices(i)%dihedral_angle*radian
    end do

    print *,'Printed ok'

  end subroutine print_zmatrix

  subroutine print_zmatrix_to_stdout (zmatrix)

    ! Wrapper for print_zmatrix when we do not specify a unit
    ! number (previously had unit as optional in above subroutine,
    ! but wanted the ordering to remain, which requires the user to
    ! utilise a keyword for the zmatrix argument -- simpler to just
    ! make this wrapper)

    type (zmatrix_object), intent(in) :: zmatrix

    ! Print the zmatrix to stdout
    call print(stdout,zmatrix)

  end subroutine print_zmatrix_to_stdout

  subroutine reorder_zmatrix (zmatrix)

    ! Wah? Why have this routine? Well for a damn good reason! It was
    ! deemed easier to index the z-matrix with 'atom ids' whilst 
    ! constructing it (so we could find the correct connections etc) 
    ! and then put it in the correct order after we finished

    type (zmatrix_object), intent(inout) :: zmatrix

    ! Local variables
    ! type (zmatrix_object) :: new_zmatrix
    integer, dimension(size(zmatrix%zmatrices)) :: new_order
    integer :: i, j

    ! Don't bother sorting this zmatrix if it has already been done
    if (zmatrix%ordered) return

    ! Allocate memory for the temporary zmatrix
    ! allocate(new_zmatrix%zmatrices(size(zmatrix%zmatrices)))

    ! Copy all the entries into the temporary z-matrix
    do i=1,size(zmatrix%zmatrices)
       ! new_zmatrix%zmatrices(zmatrix%zmatrices(i)%number) = zmatrix%zmatrices(i)
       new_order(i) = zmatrix%zmatrices(i)%number
    end do

    ! print *,new_order
    ! print *,.sort. new_order

    ! Now copy them back into the original zmatrix
    ! zmatrix%zmatrices(.sort. new_order) = zmatrix%zmatrices
    zmatrix%zmatrices = zmatrix%zmatrices(.sort. new_order) 

    ! Set the flag to say we have ordered this z-matrix
    zmatrix%ordered = .TRUE.

    ! Free the memory for the temporary z-matrix
    ! deallocate(new_zmatrix%zmatrices)

  end subroutine reorder_zmatrix
  
  ! Access functions

  function set_zmatrix_param_general (zmatrix, num_array, param_array, value_array)

    ! This the general routine for accessing the parameters in a 
    ! z-matrix object. 

    ! num_array contains in the indices of the the zmatrix entries 
    ! we wish to change 

    ! param_array contains numbers from 1 to 3 
    !     * 1 = bond length
    !     * 2 = bond angle
    !     * 3 = dihedral angle

    ! value_array contains the actual values we want to change the 
    ! parameters to -- if this is absent the parameters are unchanged

    ! Note: This function ALWAYS returns an array of the values 
    ! referred to by num_array and param_array -- not providing a
    ! value_array allows access to current values

    ! Input variables
    type (zmatrix_object) :: zmatrix
    integer, intent(in), dimension(:)          :: num_array, param_array
    real, intent(in), dimension(:), optional   :: value_array

    ! Return type
    real, dimension(size(num_array))           :: set_zmatrix_param_general

    ! Internal variables
    integer :: num_changes, i, zmat_index, error
    integer, dimension(size(num_array))        :: param_internal
    real, dimension(size(num_array))           :: value_internal

    ! Grab the total number of parameters we are accessing
    num_changes = size(num_array,1) 

    ! Ensure our parameter array is the same length as our number
    ! array -- we fill it with itself repeatedly if it is too short.
    ! This means a user can specify a single number in an array to 
    ! specify that they want all bonds, or all dihedral angles, for example
    param_internal = reshape(source=param_array, shape=(/ num_changes /), pad=param_array)

    set_zmatrix_param_general = 0.

    ! Can't do this prettily *and* efficiently, so choose the latter, as 
    ! we may be doing this in the inner loop of a Monte Carlo simulation!
    ! So .. check to see if we are setting values or just returning them
    if (present(value_array)) then
       ! Ensure our value array is the same length as our number
       ! array -- we fill it with itself repeatedly if it is too short.
       ! This means a user can specify a single number in an array to 
       ! specify that they want the same value for all the specified
       ! parameters
       value_internal = reshape(source=value_array, shape=(/ num_changes /), pad=value_array)
       do i=1,num_changes
          zmat_index = num_array(i)
          ! Use case statement to choose between the three parameters
          select case (param_internal(i))
          case (1)
             zmatrix%zmatrices(zmat_index)%bond_length    = value_internal(i)
          case (2)
             zmatrix%zmatrices(zmat_index)%bond_angle     = value_internal(i)
          case (3)
             zmatrix%zmatrices(zmat_index)%dihedral_angle = value_internal(i)
          end select
          ! Make the return array the same as the new values we took as
          ! input
          set_zmatrix_param_general = value_internal 
       end do
    else
       do i=1,num_changes
          zmat_index = num_array(i)
          ! Use case statement to choose between the three parameters
          select case (param_internal(i))
          case (1)
             set_zmatrix_param_general(i) = zmatrix%zmatrices(zmat_index)%bond_length
          case (2)
             set_zmatrix_param_general(i) = zmatrix%zmatrices(zmat_index)%bond_angle
          case (3)
             set_zmatrix_param_general(i) = zmatrix%zmatrices(zmat_index)%dihedral_angle
          end select
       end do
    end if

  end function set_zmatrix_param_general

  function set_zmatrix_single_param (zmatrix, num, param, value)

    ! This the specific routine for accessing a single parameter in a 
    ! z-matrix object. Just a wrapper for the general routine 

    ! num contains in the index of the the zmatrix entry 
    ! we wish to change 

    ! param contains a number from 1 to 3 
    !     * 1 = bond length
    !     * 2 = bond angle
    !     * 3 = dihedral angle

    ! value contains the actual value we want to change the 
    ! parameter to -- if this is absent the parameter are unchanged

    ! Note: This function ALWAYS returns a value 
    ! referred to by num and param -- not providing a
    ! value allows access to current values

    ! Input variables
    type (zmatrix_object)      :: zmatrix
    integer, intent(in)        :: num, param
    real, intent(in), optional :: value

    ! Return type
    real :: set_zmatrix_single_param

    ! Internal variables
    real, dimension(1) :: temp

    if (present(value)) then
       temp = parameter(zmatrix, (/num/), (/param/), (/value/))
    else
       temp = parameter(zmatrix, (/num/), (/param/))
    end if

    set_zmatrix_single_param = temp(1)

  end function set_zmatrix_single_param

  function set_zmatrix_conn_general (zmatrix, num_array, conn_array, value_array)

    ! This the general routine for accessing the connectivity in a 
    ! z-matrix object. 

    ! num_array contains in the indices of the the zmatrix entries 
    ! we wish to access 

    ! conn_array contains numbers from 1 to 3, signifying the connection
    ! defines one of:
    !     * 1 = bond length
    !     * 2 = bond angle
    !     * 3 = dihedral angle

    ! value_array contains the actual values we want to change the 
    ! connections to -- if this is absent the connections are unchanged

    ! Note: This function ALWAYS returns an array of the values 
    ! referred to by num_array and param_array -- not providing a
    ! value_array allows access to current values

    ! Input variables
    type (zmatrix_object) :: zmatrix
    integer, intent(in), dimension(:)           :: num_array, conn_array
    integer, intent(in), dimension(:), optional :: value_array

    ! Return type
    integer, dimension(size(num_array))         :: set_zmatrix_conn_general

    ! Internal variables
    integer :: num_changes, i, zmat_index, error
    integer, dimension(size(num_array))         :: conn_internal
    integer, dimension(size(num_array))         :: value_internal

    ! Grab the total number of parameters we are accessing
    num_changes = size(num_array,1) 

    ! Ensure our connectivity array is the same length as our number
    ! array -- we fill it with itself repeatedly if it is too short.
    ! This means a user can specify a single number in an array to 
    ! specify that they want all bond contacts, or all dihedral  contacts, 
    ! for example
    conn_internal = reshape(source=conn_array, shape=(/ num_changes /), pad=conn_array)

    set_zmatrix_conn_general = 0

    ! Check to see if we are setting values or just returning them
    if (present(value_array)) then
       ! Ensure our value array is the same length as our number
       ! array -- we fill it with itself repeatedly if it is too short.
       ! This means a user can specify a single number in an array to 
       ! specify that they want the same value for all the specified
       ! connections
       value_internal = reshape(source=value_array, shape=(/ num_changes /), pad=value_array)
       do i=1,num_changes
          zmat_index = num_array(i)
          zmatrix%zmatrices(zmat_index)%conn(conn_internal(i)) = value_internal(i)
          ! Make the return array the same as the new values we took as input
          set_zmatrix_conn_general = value_internal 
       end do
    else
       do i=1,num_changes
          zmat_index = num_array(i)
          set_zmatrix_conn_general(i) = zmatrix%zmatrices(zmat_index)%conn(conn_internal(i))
       end do
    end if

  end function set_zmatrix_conn_general

  pure integer function get_num_entries (zmatrix)

    ! This routine returns the number of entries in a z-matrix object

    ! Input variables
    type (zmatrix_object), intent(in) :: zmatrix
    
    get_num_entries = 0
    if (allocated(zmatrix%zmatrices)) get_num_entries = size(zmatrix%zmatrices)

  end function get_num_entries

  function set_atom_labels (zmatrix, num_array, label_array) result(labels)

    ! This the general routine for accessing the labels in a 
    ! z-matrix object. 

    ! num_array contains in the indices of the the zmatrix entries 
    ! we wish to access 

    ! value_array contains the new labels -- if this is absent the 
    ! atom labels are unchanged

    ! Note: This function ALWAYS returns an array of the labels 
    ! referred to by num_array -- not providing a label_array 
    ! allows access to current values

    ! Input variables
    type (zmatrix_object), intent(inout)                 :: zmatrix
    integer, intent(in), dimension(:)                    :: num_array
    character(len=*), intent(in), dimension(:), optional :: label_array

    ! Return type
    character(len=6), dimension(size(num_array))         :: labels

    ! Internal variables
    integer :: num_changes, i, num_labels, label_array_index
    ! character(len=6), dimension(size(num_array))         :: label_internal

    ! Grab the total number of labels we are accessing
    num_changes = size(num_array,1) 

    labels = ''

    ! So .. check to see if we are setting values or just returning them
    if (present(label_array)) then

       num_labels = size(label_array)
       label_array_index = 0

       ! Ensure our value array is the same length as our number
       ! array -- we fill it with itself repeatedly if it is too short.
       ! This means a user can specify a single number in an array to 
       ! specify that they want the same value for all the specified
       ! parameters
       ! label_internal = reshape(source=label_array, shape=(/ num_changes /), pad=label_array)
       do i=1,num_changes
          label_array_index = label_array_index + 1
          if (label_array_index > num_labels) label_array_index = 1
          zmatrix%zmatrices(num_array(i))%label = label_array(label_array_index)(:min(len(label_array),len(labels)))
          ! Make the return array the same as the new labels we took as input
          ! labels = label_internal 
       end do
    end if
    do i=1,num_changes
       labels(i) = zmatrix%zmatrices(num_array(i))%label
    end do

  end function set_atom_labels

  function get_all_atom_labels (zmatrix) result(labels)

    ! This a specific routine for accessing all the labels in a 
    ! z-matrix object. 

    ! Input variables
    type (zmatrix_object), intent(in)                    :: zmatrix

    ! Return type
    character(len=6), dimension(size(zmatrix%zmatrices)) :: labels

    labels = zmatrix%zmatrices(:)%label

  end function get_all_atom_labels

  function set_one_atom_label (zmatrix, index, newlabel) result(label)

    ! This a specific routine for accessing just one label in a 
    ! z-matrix object. 

    ! Input variables
    type (zmatrix_object), intent(inout)   :: zmatrix
    integer, intent(in)                    :: index
    character(len=*), intent(in), optional :: newlabel

    ! Return type
    character(len=6) :: label

    if (present(newlabel)) zmatrix%zmatrices(index)%label = newlabel(:min(len(newlabel),len(label)))

    label = zmatrix%zmatrices(index)%label

  end function set_one_atom_label
  
  function set_zmatrix_info (zmatrix, info)

    ! This the general routine for accessing the information string
    ! in a z-matrix object. 

    ! info contains the information to save inside the object. This
    ! is optional
    
    ! Note: This function ALWAYS returns the information string, 
    ! so not specifying info has the side-effect of just returning
    ! the current value

    ! Return type
    ! type (varying_string) :: set_zmatrix_info
    character(len=180) :: set_zmatrix_info

    ! Input variables
    type (zmatrix_object)      :: zmatrix
    character(len=*), optional :: info
    
    if(present(info)) zmatrix%info = info  

    set_zmatrix_info = zmatrix%info

  end function set_zmatrix_info

  function makexyz_rasmol(zmatrix)

    ! Input variables
    type (zmatrix_object), intent(in) :: zmatrix

    ! Return type of function
    real, dimension(3,size(zmatrix%zmatrices)) :: makexyz_rasmol

    ! This code was converted from the babel 1.6 C routine intcart:

    ! /*-----------------------------------------------------
    ! Notice of blatant theft.  This code was totally lifted
    ! from Roger Sayle's Internal2Cartesian routine in RasMol 2.6.
    ! Thanks Roger.
    
    ! Roger's code was a whole lot cleaner than the earlier version
    ! which was stolen from MOPAC 5.  However, by looking at 
    ! Roger's varialbe names it looks like his code was derived
    ! from MOPAC also.  Oh well, what comes around goes around :-)
    
    ! ------------------------------------------------------*/

    real(kind=8) :: cosph,sinph,costh,sinth,coskh,sinkh
    real(kind=8) :: cosa,sina,cosd,sind
    real(kind=8) :: dist,angle,dihed
    
    real(kind=8) :: xpd,ypd,zpd,xqd,yqd,zqd
    real(kind=8) :: xa,ya,za,xb,yb,zb
    real(kind=8) :: rbc,xyb,yza,temp
    real(kind=8) :: xpa,ypa,zqa
    real(kind=8) :: xd,yd,zd
    logical :: flag
    integer :: i, na, nb, nc, num_atoms

    num_atoms = size(zmatrix%zmatrices)
    
    ! allocate(makexyz_rasmol(3,num_atoms))

    ! Atom #1
    makexyz_rasmol(:,1) = (/ 0.0, 0.0, 0.0 /)
    
    if (num_atoms == 1) return
  
    ! Atom #2 - by convention this atom is placed along the z-axis
    makexyz_rasmol(:,2) = (/ 0.0, 0.0, zmatrix%zmatrices(2)%bond_length /)
  
    if (num_atoms == 2) return
  
    ! Atom #3 - by convention this atom is placed in the x-z plane
    dist = zmatrix%zmatrices(3)%bond_length
    angle = zmatrix%zmatrices(3)%bond_angle
    cosa = cos(angle);
    sina = sin(angle);
  
    ! Check to see if atom 3 is bonded to atom 1 or atom 2 and
    ! invoke the appropriate simple trigonometric relation to
    ! determine the z position
    if( zmatrix%zmatrices(3)%conn(1) == 1 ) then
       makexyz_rasmol(3,3) = makexyz_rasmol(3,1) + cosa*dist;
    else ! must be equal to 2
       makexyz_rasmol(3,3) = makexyz_rasmol(3,2) - cosa*dist;
    end if
    makexyz_rasmol(1,3) = sina*dist;
    makexyz_rasmol(2,3) = 0.0;
  
    ! Cycle through the rest of the atoms
    do i = 4, num_atoms

       dist = zmatrix%zmatrices(i)%bond_length
       angle = zmatrix%zmatrices(i)%bond_angle
       dihed = zmatrix%zmatrices(i)%dihedral_angle
       
       na = zmatrix%zmatrices(i)%conn(1)
       nb = zmatrix%zmatrices(i)%conn(2)
       nc = zmatrix%zmatrices(i)%conn(3)
    
       xb = makexyz_rasmol(1,nb) - makexyz_rasmol(1,na)
       yb = makexyz_rasmol(2,nb) - makexyz_rasmol(2,na)
       zb = makexyz_rasmol(3,nb) - makexyz_rasmol(3,na)
    
       rbc = xb*xb + yb*yb + zb*zb;
       if( rbc < 0.0001 ) then
          ! Atoms are coincident. Fatal error
          write(stderr,*)'Error! Atoms ',na,' and ',nb,' are coincident'
          stop
       end if
       rbc = 1.0d0/sqrt(rbc)
    
       cosa = cos(angle)
       sina = sin(angle)
    
       if( abs(cosa) >= 0.999999 ) then
          ! Colinear
          temp = dist*rbc*cosa;
          makexyz_rasmol(1,i) = makexyz_rasmol(1,na) + temp*xb
          makexyz_rasmol(2,i) = makexyz_rasmol(2,na) + temp*yb
          makexyz_rasmol(3,i) = makexyz_rasmol(3,na) + temp*zb
       else
          xa = makexyz_rasmol(1,nc) - makexyz_rasmol(1,na)
          ya = makexyz_rasmol(2,nc) - makexyz_rasmol(2,na)
          za = makexyz_rasmol(3,nc) - makexyz_rasmol(3,na)
          
          sind = -sin(dihed)
          cosd = cos(dihed)
      
          xd = dist*cosa
          yd = dist*sina*cosd
          zd = dist*sina*sind
      
          xyb = sqrt(xb*xb + yb*yb)

          if( xyb < 0.1 ) then
             ! Rotate about y-axis!
             temp = za; za = -xa; xa = temp;
             temp = zb; zb = -xb; xb = temp;
             xyb = sqrt(xb*xb + yb*yb)
             flag = .TRUE.
          else 
             flag = .FALSE.
          end if
      
          costh = xb/xyb
          sinth = yb/xyb
          xpa = costh*xa + sinth*ya
          ypa = costh*ya - sinth*xa
          
          sinph = zb*rbc
          cosph = sqrt(1.0 - sinph*sinph)
          zqa = cosph*za  - sinph*xpa
          
          yza = sqrt(ypa*ypa + zqa*zqa)
          
          if( yza > 1.0E-10 ) then
             coskh = ypa/yza
             sinkh = zqa/yza
             ypd = coskh*yd - sinkh*zd
             zpd = coskh*zd + sinkh*yd
          else
             ! Angle so small as to be unimportant
             ! i.e. coskh = 1.0 & sinkh = 0.0
             ypd = yd
             zpd = zd
          end if

       end if
      
       xpd = cosph*xd  - sinph*zpd
       zqd = cosph*zpd + sinph*xd
       xqd = costh*xpd - sinth*ypd
       yqd = costh*ypd + sinth*xpd
      
       if( flag ) then
          ! Rotate about y-axis!
          makexyz_rasmol(1,i) = makexyz_rasmol(1,na) - zqd
          makexyz_rasmol(2,i) = makexyz_rasmol(2,na) + yqd
          makexyz_rasmol(3,i) = makexyz_rasmol(3,na) + xqd
       else
          makexyz_rasmol(1,i) = makexyz_rasmol(1,na) + xqd
          makexyz_rasmol(2,i) = makexyz_rasmol(2,na) + yqd
          makexyz_rasmol(3,i) = makexyz_rasmol(3,na) + zqd
       end if
    end do
  end function makexyz_rasmol

  ! Equivalence functions

  function zmat_eq_zmat(zmatrixl, zmatrixr)

    ! This the general routine for deciding if two z-matrices are equivalent

    ! Return type
    logical :: zmat_eq_zmat

    ! Input variables
    type (zmatrix_object), intent(in) :: zmatrixl, zmatrixr

    ! Internal variables
    integer :: numl, numr, i

    ! Initialise the return to false, so that any return below will
    ! indicate failure
    zmat_eq_zmat = .FALSE.

    ! Grab the total number of entries in each zmatrix
    numl = size(zmatrixl%zmatrices) 
    numr = size(zmatrixr%zmatrices) 

    if (numl /= numr) return

    do i=1,numl
       if (sum(zmatrixl%zmatrices(i)%conn - zmatrixr%zmatrices(i)%conn) /= 0) return
       if (abs(zmatrixl%zmatrices(i)%bond_length - zmatrixr%zmatrices(i)%bond_length) > tolerance) return
       if (abs(zmatrixl%zmatrices(i)%bond_angle - zmatrixr%zmatrices(i)%bond_angle) > tolerance) return
       if (abs(zmatrixl%zmatrices(i)%dihedral_angle - zmatrixr%zmatrices(i)%dihedral_angle) > tolerance) return
    end do

    ! If we checked all the entries in the zmatrix and had no problems
    zmat_eq_zmat = .TRUE.

  end function zmat_eq_zmat

  function zmat_neq_zmat(zmatrixl, zmatrixr)

    ! This the general routine for deciding if two z-matrices are *not*
    ! equivalent (just calls equivalency routine and negates it)

    ! Return type
    logical :: zmat_neq_zmat

    ! Input variables
    type (zmatrix_object), intent(in) :: zmatrixl, zmatrixr

    zmat_neq_zmat = .TRUE.

    ! Just call the equivalence routine and return the opposite value
    if (zmatrixl == zmatrixr) zmat_neq_zmat = .FALSE.

  end function zmat_neq_zmat

end module zmatrix_class
