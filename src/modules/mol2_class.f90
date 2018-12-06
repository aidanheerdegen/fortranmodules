module mol2_class


  use file_functions, only: stdout, stderr, open, read_buffer
  use string_functions, only: getword, gettext, operator(.ucase.)
  use variable_array, only: push, splice
  use sort_functions, only: sort
  use cartesian_class, only: connectivity_entry, operator(==), operator(/=)

  implicit none

  private

  ! Using this module you can read a Sybyl mol2 format file into
  ! a mol2 object. Routines are provided to access the information
  ! contained in the mol2 object.

  ! Much of this is nicked from Jay William Ponder (author of Tinker)
  ! but has been extensively reworked.  

  !! $Log: mol2_class.f90,v $
  !! Revision 1.12  2007/11/05 01:28:00  aidan
  !! No longer use the bond_id field. This caused a bug when bonds were deleted.
  !! Now just use the loop counter which goes from 1 .. num_bonds.
  !!
  !! Revision 1.11  2006/06/14 04:08:40  aidan
  !! Added a type/types interface to access the atom types (identical
  !! to the label/labels interface). Also added a function 'as_element'
  !! to convert from a Sybyl atom_type to an element symbol. It is an
  !! elemental function, so will operate on single atom types, or whole
  !! arrays.
  !!
  !! Revision 1.10  2006/06/14 03:22:55  aidan
  !! Added an extra interface to label so that we can get a label
  !! by referring to an atom position relative to a substructure.
  !! This is actually really important, as we can retrieve atom
  !! lists by substructure and then wouldn't have been able to find
  !! out the atom label.
  !!
  !! Also replaced a bunch of interface variables called subst_id
  !! with subst. This is cleaner and clearer to a user I reckon.
  !!
  !! Revision 1.9  2006/05/30 06:15:27  aidan
  !! Found a bug in readmol2 with regards to initialising the root_atom
  !! of the substructures to zero before parsing. This has been incorporated
  !! into the init_mol2 routine.
  !!
  !! Revision 1.8  2006/05/29 05:15:15  aidan
  !! Added a 'label' function which returns the label of a single atom in the
  !! mol2 object.
  !!
  !! Revision 1.7  2006/05/25 05:06:12  aidan
  !! Added a 'bonds' enquiry function. Took the private attribute off all the
  !! minor defined types. Silly to make them private if there are no functions
  !! defined to enquire and set their values.
  !!
  !! Revision 1.6  2006/05/18 01:13:30  aidan
  !! Removed the redundant ordering information in the substructure type.
  !! Initialised all the pointer components of a mol12 object to null().
  !! I hope this doesn't lead to memory leaks.
  !! Added explicit init and destroy routines for mol2 objects.
  !! Added a function to add one mol2 object to another.
  !! Added a bond number enquiry function.
  !! Memory allocation removed from readmol2 and replaced with call to
  !! init().
  !! The last update commented out the code in readmol2 that figured
  !! out the last_atom numbers for substructures. This was a mistake and
  !! has been undone.
  !! Added a test to the equality operator to test that the connectivity
  !! arrays in the mol2 objects were the same.
  !!
  !! Revision 1.5  2006/05/08 04:28:09  aidan
  !! Added a bond type access function.
  !!
  !! Added an explicit mol2 to mol2 assignment routine (does a deep copy)
  !!
  !! Made a clone (effectively an alias) of the print interface, called write.
  !!
  !! Revision 1.4  2003/09/10 06:51:36  aidan
  !! Changed the mol2 reading function so that now it uses a consistent
  !! interface for reading from the file--a read_buffer routine from
  !! the file_functions module ensures that all comment lines are correctly
  !! parsed (i.e. ignored).
  !!
  !! Revision 1.3  2003/09/10 00:02:31  aidan
  !! Added a subroutine to print mol2 objects. Added equality (==) and
  !! inequality (/=) operators which compare two mol2 objects. Made a small
  !! change to the mol2 reading function--a comment is only read if there
  !! is a place holder for the 'status bits' field (in the MOLECULE section).
  !!
  !! Revision 1.2  2003/08/21 03:51:12  aidan
  !! Re-wrote most of the object access functions, converting them from a single
  !! function with an optional substructure argument to two separate functions.
  !! Did this so I can dynamically size the return arrays/strings with calls to
  !! pure functions or simple reference to mol object parameters. Could not do
  !! this with an optional argument approach. This makes for a cleaner interface
  !! which no longer returns pointers (so code that assumes this will break) and
  !! does no direct memory allocation. This should prevent potential memory
  !! leaks.
  !!
  !! As a consequence of the rewrite I decided to change the way substructure
  !! access functions returned their data. Previously the module had honoured
  !! the absolute id numbers in the mol2 object so that ordering and connectivity
  !! arrays would still be valid. This doesn't really work with these auto-sized
  !! arrays, and if someone had allocated enough space for an array off their own
  !! bat and just assigned the results of a call to coords to their array the
  !! special indexing would be lost anyway. So now the module adjusts ordering
  !! and connectivity arrays so that they are relative to their own substructure,
  !! assuming that the function was called with a substructure argument. So
  !! atom ids are renumbered from 1 and the appropriate offsets are subtracted
  !! from ordering and connectivity arrays. Inter-substructure bonding is
  !! problematic in this case, but this easily overcome by returning ALL the
  !! coordinates, ordering and connectivity entries and then just choosing the
  !! substructure required.
  !!
  !! Removed all the _mol2_ qualifiers in the routine names. Made them long,
  !! ugly and hard to read and they were totally redundant.
  !!
  !! Revision 1.1  2003/08/20 00:25:44  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: mol2_class.f90,v 1.12 2007/11/05 01:28:00 aidan Exp $"

  type atom_entry

     ! Elements defined in an @<TRIPOS>ATOM entry

     ! Coordinates (x,y,z) and the charge on the atom
     real             :: coords(3), charge
     ! Atom name and type (type is the SYBYL atom type, e.g. C.2 = sp2 
     ! hybridised Carbon
     character(len=6) :: atom_name, atom_type
     ! Substructure id
     integer          :: subst_id

  end type atom_entry

  type subst_entry

     ! Elements defined in  an @<TRIPOS>SUBSTRUCTURE entry

     ! Substructure name
     character(len=6) :: subst_name
     ! Root atom of the substructure (i.e. where it starts)
     integer          :: root_atom, last_atom

     !  There is other info associated with substructure (namely the
     ! type, dictionary, chain to which it belonds, subtype, number of
     ! inter substructure bonds, the ubiquituous status which we are
     ! ignoring completely, and a comment) but we will not support that
     ! currently as it is related to biological molecules.

  end type subst_entry

  type bond_entry

     ! Elements defined in an @<TRIPOS>BOND entry

     ! Connectivity array, the two integers are the atom_ids of the 
     ! origin and target atoms respectively
     integer :: conn(2)

     ! Bond type - integer value from 1 to 8, see bond_type_lookup for
     ! what these values mean
     integer :: bond_type

  end type bond_entry

  ! Bond type lookup
  character(len=2), dimension(8), parameter  :: bond_type_lookup = (/ &
       '1 ',   & ! = single
       '2 ',   & ! = double
       '3 ',   & ! = triple
       'am',   & ! = amide
       'ar',   & ! = aromatic
       'du',   & ! = dummy
       'un',   & ! = unknown
       'nc'  /)  ! = not connected

  type mol2_object

     private

     ! For more information on the mol2 file format, see
     !
     !        http://www.tripos.com/custResources/mol2Files/
     !
     
     ! Elements defined in @<TRIPOS>MOLECULE

     character(len=80) :: mol_name, mol_comment
     ! type (varying_string) :: mol_name, mol_comment

     ! Molecule type can be SMALL, PROTEIN, NUCLEIC_ACID or SACCHARIDE
     character(len=20) :: mol_type

     ! Charge type can be one of:
     !
     !     NO_CHARGES, DEL_RE, GASTEIGER, GAST_HUCK, HUCKEL, PULLMAN, 
     !     GAUSS80_CHARGES, AMPAC_CHARGES, MULLIKEN_CHARGES, DICT_CHARGES, 
     !     MMFF94_CHARGES or USER_CHARGES
     character(len=20) :: charge_type

     integer :: num_atoms, num_bonds, num_subst, num_feat, num_sets

     ! Elements defined in @<TRIPOS>ATOM 
     !
     !     o The yet to be determined dimension
     !       corresponds to the atom_id and is of
     !       length num_atoms

     type (atom_entry), dimension(:), pointer :: atoms => null()

     ! This is an optional extra -- we take a note of the order in which
     ! the atoms are specified in the mol2 file -- this is to allow an
     ! easy method for reordering the atoms, e.g. when later converting to
     ! a z-matrix
     integer, dimension(:), pointer :: order => null()

     ! Elements defined in @<TRIPOS>BOND
     !
     !       The yet to be determined dimension 
     !       corresponds to the bond_id and is of
     !       length num_bonds

     type (bond_entry), dimension(:), pointer :: bonds => null()
     
     ! This is 'derived information' from the bond section
     ! above -- an array of the same length as the number of atoms, 
     ! each entry being an array of atom ids to which that atom is 
     ! connected. We use the variable_array module to dynamically update
     ! these arrays
     type (connectivity_entry), dimension(:), pointer :: connections => null()

     ! Elements defined in @<TRIPOS>SUBSTRUCTURE
     !
     !     o The yet to be determined dimension
     !       corresponds to the subst_id and is of
     !       length num_subst

     type (subst_entry), dimension(:), pointer :: substructures => null()

     ! Elements defined in @<TRIPOS>CRYSIN

     ! The cell parameters a, b, c, alpha, beta, gamma
     real, dimension(6) :: cell
     ! The space group and setting (integer values as defined in the
     ! Intl Tables).
     integer            :: space_grp, setting

  end type mol2_object

  ! Initialise a mol2 object
  interface init
     module procedure init_mol2
  end interface

  ! Destroy a mol2 object
  interface destroy
     module procedure destroy_mol2
  end interface

  ! Procedure to load a mol2 object
  interface load_mol2
     module procedure open_then_read, readmol2
  end interface

  ! Make a generic interface to mol2 print routine
  interface print
     module procedure print_mol2, open_then_print, print_mol2_to_stdout
  end interface

  ! Make write a copy of print
  interface write
     module procedure print_mol2, open_then_print, print_mol2_to_stdout
  end interface

  ! Function to grab the coordinates from a mol2 object
  interface coords
     module procedure get_all_coords, get_coords
  end interface

  ! Function to grab the bonds from a mol2 object
  interface bonds
     module procedure get_all_bonds
  end interface

  ! Returns an ordering array which are used to index coordinate arrays 
  interface order
     module procedure get_all_order, get_order
  end interface

  ! Function to grab the connectivity from a mol2 object
  interface connectivity
     module procedure get_all_connectivity, get_connectivity
  end interface

  ! Function to grab the atoms labels from a mol2 object
  interface labels
     module procedure get_all_atom_labels, get_atom_labels
  end interface

  ! Function to grab the atom types from a mol2 object
  interface types
     module procedure get_all_atom_types, get_atom_types
  end interface

  ! Function to convert a Sybyl atom type to an element symbol
  interface as_element
     module procedure get_element_symbol_from_type
  end interface

  ! Function to grab the label for a single atom
  interface label
     module procedure get_atom_label, get_relative_atom_label
  end interface

  ! Function to grab the label for a single atom
  interface type
     module procedure get_atom_type, get_relative_atom_type
  end interface

  ! Function to get the number of atoms from a mol2 object
  interface atom_num
     module procedure get_num_atoms
  end interface

  ! Function to get the number of atoms from a mol2 object
  interface bond_num
     module procedure get_num_bonds
  end interface

  ! Function to get the cell parameters
  interface cell
     module procedure get_cell_parameters
  end interface

  ! Function to get the space group
  interface group
     module procedure get_space_group
  end interface

  ! Function to get the setting info
  interface setting
     module procedure get_setting
  end interface

  ! Function to get the number of substructures from a mol2 object
  interface sub_num
     module procedure get_num_substructures, get_substructure_root 
  end interface

  ! Function to get the name of the mol2 object
  interface name
     module procedure get_name, set_name
  end interface

  ! Function to get the comment field out of the mol2 object
  interface comment
     module procedure get_comment, set_comment
  end interface

  ! Define assignment from a filename and a deep copy from another mol2 object
  interface assignment (=)
     module procedure assign_from_filename, assign_from_mol2
  end interface

  ! Define a routine to add atoms from one mol2 file to another
  interface add
     module procedure add_mol2
  end interface

  ! Define equality
  interface operator (==)
     module procedure mol2_eq_mol2
  end interface

  ! Define inequality
  interface operator (/=)
     module procedure mol2_neq_mol2
  end interface

  ! Generic interface to bond type routines (private)
  interface btype
     module procedure bond_type_from_string, bond_type_from_integer
  end interface

  ! Public data types
  public :: atom_entry, bond_entry, subst_entry, mol2_object

  ! Public routines
  public :: init, destroy, load_mol2, add
  public :: coords, atom_num, bonds, bond_num, sub_num, connectivity, order, cell
  public :: labels, label, types, type, as_element
  public :: name, comment, group, setting, print, write, btype

  ! Overloaded operators
  public :: assignment(=), operator(==), operator(/=)
       
contains


  ! Access routines .. 


  ! There are generally two separate routines for accessing data
  ! in the mol2_object. One for accessing all the relevant data
  ! and another for accessing just that from one substructure.
  ! This is so that we can automatically determine the size of
  ! the function result using a call to atom_num -- otherwise we
  ! would have to use an optional substructure id argument which
  ! wouldn't be able to be automatically passed to atom_num

  function get_all_coords(mol2) result(coords)

    ! Return a 3 x n real array, where n is the number of atoms
    ! in the mol2 object

    ! Input variable
    type (mol2_object), intent(in)     :: mol2

    ! Return value -- the size is determined by the number of atom
    ! entries in the mol2 object
    real, dimension(3,size(mol2%atoms)) :: coords

    integer :: i

    do i = 1, mol2%num_atoms
       coords(:,i) = mol2%atoms(i)%coords(:)
    end do
  end function get_all_coords

  function get_coords(mol2, subst) result(coords)

    ! Return a 3 x n real array, where n is the number of atoms
    ! in the specified substructure of the mol2 object

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst

    ! Return value -- call to atom_num determines the size
    real, dimension(3,atom_num(mol2,subst)) :: coords

    integer :: i, index

    if (subst > mol2%num_subst) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    index = 0
    do i = mol2%substructures(subst)%root_atom, mol2%substructures(subst)%last_atom
       index = index+1
       coords(:,index) = mol2%atoms(i)%coords(:)
    end do

  end function get_coords

  function get_all_bonds(mol2) result(bonds)

    ! Input variable
    type (mol2_object), intent(in)     :: mol2

    ! Return value
    type (bond_entry), dimension(size(mol2%bonds)) :: bonds

    bonds = mol2%bonds

  end function get_all_bonds

  function get_all_order(mol2) result(ordering)

    ! Returns an integer array of atom ids in the order that they
    ! appeared in the mol2 file. 

    ! Arguments
    type (mol2_object), intent(in) :: mol2

    ! Return type of function -- the size is determined by the number of atom
    ! entries in the mol2 object
    integer, dimension(size(mol2%order)) :: ordering
    
    ordering = mol2%order

  end function get_all_order

  function get_order(mol2, subst) result(ordering)

    ! Returns an integer array of atom ids in the order that they
    ! appeared in the mol2 file. Optionally a substructure id can
    ! be specified which will limit the array to only that substructure.

    ! Arguments
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst

    ! Return type of function -- the size is determined by a call to atom_num
    integer, dimension(atom_num(mol2,subst)) :: ordering
    
    ! Internal variables
    integer :: i, index, root_atom, last_atom, atom_id

    if (subst > mol2%num_subst) then
       ! Attempting to get the ordering for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    root_atom = mol2%substructures(subst)%root_atom
    last_atom = mol2%substructures(subst)%last_atom

    ! Cycle through the ordering array pushing atom ids onto our array
    ! if they are within the range dictated by the substructure limits
    index = 0
    do i = lbound(mol2%order,1), ubound(mol2%order,1)
       atom_id = mol2%order(i)
       if ( root_atom <= atom_id .and. atom_id <= last_atom )  then
          index = index + 1
          ! Note that we don't pass the *actual* id, rather one that
          ! is relative to its root atom -- this is because it is
          ! intended to index the arrays returned from the routines above,
          ! and their indices go from 1 .. num_atoms_in_substructure
          ordering(index) = atom_id - root_atom + 1
       end if
    end do

  end function get_order

  function get_all_atom_labels(mol2) result(atom_labels)

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    ! Return type
    character(len=6), dimension(size(mol2%atoms))  :: atom_labels

    integer :: i

    do i = 1, mol2%num_atoms
       atom_labels(i) = mol2%atoms(i)%atom_name
    end do

  end function get_all_atom_labels

  function get_atom_labels(mol2, subst) result(atom_labels)

    ! Return a character array of atom labels

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst

    ! Return value -- call to atom_num determines the size
    character(len=6), dimension(atom_num(mol2,subst)) :: atom_labels

    integer :: i, index

    if (subst > mol2%num_subst) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    index = 0
    do i = mol2%substructures(subst)%root_atom, mol2%substructures(subst)%last_atom
       index = index+1
       atom_labels(index) = mol2%atoms(i)%atom_name
    end do

  end function get_atom_labels

  function get_atom_label(mol2, atom) result(label)

    ! Return a character array of atom labels

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: atom

    ! Return value
    character(len=6) :: label

    integer :: i, index

    if (atom > mol2%num_atoms) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Atom id does not exist: ', atom
       stop
    end if

    label = mol2%atoms(atom)%atom_name

  end function get_atom_label

  function get_relative_atom_label(mol2, subst, atom) result(label)

    ! Return the atom label specified relative to a substructure

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst, atom

    ! Return value
    character(len=6) :: label

    integer :: i, index

    if (subst > mol2%num_subst) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    if (atom > atom_num(mol2,subst)) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,'(A,I0,A,I0)')'MOL2_CLASS :: Atom id ',atom,' does not exist in substructure ',subst
       stop
    end if

    label = mol2%atoms(atom+sub_num(mol2,subst)-1)%atom_name

  end function get_relative_atom_label

  function get_all_atom_types(mol2) result(atom_types)

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    ! Return type
    character(len=6), dimension(size(mol2%atoms))  :: atom_types

    integer :: i

    do i = 1, mol2%num_atoms
       atom_types(i) = mol2%atoms(i)%atom_type
    end do

  end function get_all_atom_types

  function get_atom_types(mol2, subst) result(atom_types)

    ! Return a character array of atom types

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst

    ! Return value -- call to atom_num determines the size
    character(len=6), dimension(atom_num(mol2,subst)) :: atom_types

    integer :: i, index

    if (subst > mol2%num_subst) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    index = 0
    do i = mol2%substructures(subst)%root_atom, mol2%substructures(subst)%last_atom
       index = index+1
       atom_types(index) = mol2%atoms(i)%atom_type
    end do

  end function get_atom_types

  function get_atom_type(mol2, atom) result(type)

    ! Return a character array of atom types

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: atom

    ! Return value
    character(len=6) :: type

    integer :: i, index

    if (atom > mol2%num_atoms) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Atom id does not exist: ', atom
       stop
    end if

    type = mol2%atoms(atom)%atom_type

  end function get_atom_type

  function get_relative_atom_type(mol2, subst, atom) result(type)

    ! Return the atom type specified relative to a substructure

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst, atom

    ! Return value
    character(len=6) :: type

    integer :: i, index

    if (subst > mol2%num_subst) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    if (atom > atom_num(mol2,subst)) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,'(A,I0,A,I0)')'MOL2_CLASS :: Atom id ',atom,' does not exist in substructure ',subst
       stop
    end if

    type = mol2%atoms(atom+sub_num(mol2,subst)-1)%atom_type

  end function get_relative_atom_type

  elemental function get_element_symbol_from_type(type) result(symbol)

    character(len=6), intent(in) ::  type

    character(len=2) ::  symbol

    symbol = type(1:2)

    if (symbol(2:2) == '.') symbol(2:2) = ' '

  end function get_element_symbol_from_type

  function get_all_connectivity(mol2) result(connectivity)

    ! Returns a ragged array of connectivities. In this case the length
    ! of this array should equal the number of atoms in the mol2 object.
    ! The array index is atom ids, and each entry is a variable length 
    ! integer array (a connectivity entry) of atom ids to which that 
    ! atom id is connected.

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    ! Return type
    type (connectivity_entry), dimension(size(mol2%connections)) :: connectivity

    integer :: i

    connectivity = mol2%connections

  end function get_all_connectivity

  function get_connectivity(mol2, subst) result(connectivity)

    ! Returns a ragged array of connectivities. In this case the length
    ! of this array should equal the number of atoms in the specified
    ! substructure. The array index is atom ids, and each entry is a 
    ! variable length integer array (a connectivity entry) of atom ids 
    ! to which that atom id is connected. Note that the indices and ids
    ! in the arrays are relative to the atoms within the substructure,
    ! and do not correspond to the numbers in the mol2 file itself. They
    ! are internally consistent with the other routines which return
    ! coordinates and ordering arrays of individual substructures.

    ! Input variable
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst

    ! Return type
    type (connectivity_entry), dimension(atom_num(mol2,subst)) :: connectivity

    integer :: i, index, dummy

    if (subst > mol2%num_subst) then
       ! Attempting to get coordinates for a substructure that doesn't exist
       write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
       stop
    end if

    index = 0
    do i = mol2%substructures(subst)%root_atom, mol2%substructures(subst)%last_atom
       index = index+1
       nullify(connectivity(index)%connectivity)
       ! Push the connectivity array for this atom on to the array
       ! we are returning. The id of the root atom for this substructure
       ! is subracted from the connectivity array--this connectivity array will
       ! be used with the coordinate arrays above, which are indexed from
       ! 1 .. num_atom_in_substructure. No provision is made for bonds between 
       ! substructures.
       dummy = push(connectivity(index)%connectivity, &
            mol2%connections(i)%connectivity - mol2%substructures(subst)%root_atom + 1)
    end do

  end function get_connectivity

  elemental integer function get_num_atoms(mol2, subst) result(num_atoms)

    ! This function returns the number of atoms in a particular
    ! substructure, or the total number of atoms if subst is
    ! omitted. It is important that this function is "pure" (an
    ! elemental function is, by definition, pure) as we use it
    ! to dimension the results of a number of the access routines.
    ! Being pure is a requirement for such functionality.

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, optional, intent(in)  :: subst

    num_atoms = 0

    if (present(subst)) then
       if (subst > mol2%num_subst) then
          ! Attempting to access a substructure that doesn't exist
          ! Can't output any error messages because we are inside a
          ! pure function
          ! write(stderr,*)'MOL2_CLASS :: Substructure id does not exist: ', subst
          return
       end if
       num_atoms = mol2%substructures(subst)%last_atom -  &
            mol2%substructures(subst)%root_atom + 1
    else
       num_atoms = mol2%num_atoms
    end if

  end function get_num_atoms

  elemental integer function get_num_bonds(mol2) result(num_bonds)

    ! This function returns the number of bonds in a mol2 object

    ! Input variables
    type (mol2_object), intent(in) :: mol2

    num_bonds = mol2%num_bonds

  end function get_num_bonds

  elemental integer function get_num_substructures(mol2)

    ! Returns the number of substructures in a mol2 object

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    get_num_substructures = mol2%num_subst

  end function get_num_substructures

  elemental integer function get_substructure_root(mol2, subst)

    ! Returns the actual atom id of the root atom of a given
    ! substructure, i.e. the lowest atom id for a substructure

    ! Input variables
    type (mol2_object), intent(in) :: mol2
    integer, intent(in)            :: subst

    if (subst > mol2%num_subst) then
       ! Return zero if there is no such substructure
       get_substructure_root = 0
    else
       get_substructure_root = mol2%substructures(subst)%root_atom
    end if

  end function get_substructure_root

  function get_cell_parameters(mol2)

    ! Returns the cell parameters of a mol2 object

    ! Return type of function -- cell parameters: a, b, c, alpha, beta, gamma
    real, dimension(6) :: get_cell_parameters

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    get_cell_parameters = mol2%cell

  end function get_cell_parameters

  elemental function get_space_group(mol2)

    ! Returns the space group of a mol2 object

    ! Return type of function -- space group number as per the international tables
    integer :: get_space_group

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    get_space_group = mol2%space_grp

  end function get_space_group

  elemental function get_setting(mol2)

    ! Returns the setting number of a mol2 object

    ! Return type of function -- setting number as per the international tables
    integer :: get_setting

    ! Input variable
    type (mol2_object), intent(in) :: mol2

    get_setting = mol2%setting

  end function get_setting
  
  function get_name(mol2)

    ! Returns the molecule name in a mol2_object. This functionality
    ! used to be incorporated into set_name, but by making dedicated
    ! functions we can also return character strings of exactly the
    ! right length. Nice feature.

    ! Input variable
    type (mol2_object) :: mol2

    ! Return type of function -- string
    character(len=len_trim(mol2%mol_name)) :: get_name

    get_name = mol2%mol_name(1:len(get_name))

  end function get_name
  
  function set_name(mol2, name)

    ! Sets (and returns) the molecule name in a mol2_object. Why a
    ! functional interface? We want to have just one command to
    ! return and set the name, and f90 dosen't allow functions and
    ! subroutines in the same generic interface
    
    ! Input variables
    type (mol2_object) :: mol2
    character(len=*)   :: name

    ! Return type of function -- string. Note that the length is 
    ! determined so that it is the same as that returned by get_name.
    ! The input name will be truncated to the length of the character
    ! variable that is used to store these things internally, and then
    ! it is trimmed of blank space when returned.
    character(len=len_trim(name(:len(mol2%mol_name)))) :: set_name

    mol2%mol_name = name

    set_name = mol2%mol_name

  end function set_name
  
  function get_comment(mol2)

    ! Input variable
    type (mol2_object) :: mol2

    ! Return type 
    character(len=len_trim(mol2%mol_comment)) :: get_comment 

    get_comment = mol2%mol_comment

  end function get_comment
  
  function set_comment(mol2, comment)

    ! Input variables
    type (mol2_object) :: mol2
    character(len=*)   :: comment

    ! Return type of function -- string. Note that the length is 
    ! determined so that it is the same as that returned by get_comment.
    ! The input comment will be truncated to the length of the character
    ! variable that is used to store these things internally, and then
    ! it is trimmed of blank space when returned.
    character(len=len_trim(comment(:len(mol2%mol_comment)))) :: set_comment

    mol2%mol_comment = comment

    set_comment = mol2%mol_comment

  end function set_comment


  subroutine destroy_mol2(mol2)
    
    ! Free the memory associated with a mol2 object
    
    type (mol2_object), intent(inout) :: mol2

    ! Local variables
    integer :: i, j

    if (associated(mol2%atoms)) deallocate(mol2%atoms)
    if (associated(mol2%order)) deallocate(mol2%order)
    if (associated(mol2%connections)) then
       do i = 1, mol2%num_atoms
          ! deallocate(mol2%connections(i)%connectivity)
          j  = splice(mol2%connections(i)%connectivity,0)
       end do
       deallocate(mol2%connections)
    end if

    if (associated(mol2%bonds)) deallocate(mol2%bonds)
    if (associated(mol2%substructures)) deallocate(mol2%substructures)
    
  end subroutine destroy_mol2

  subroutine init_mol2(mol2, natoms, nbonds, nsub, error)

    ! Allocate memory all the objects inside a mol2 object

    type (mol2_object), intent(inout) :: mol2
    integer, intent(in)               :: natoms, nbonds, nsub
    integer, intent(out), optional    :: error

    ! Local variables
    integer :: lerror, subst_id

    call destroy(mol2)

    mol2%num_atoms = natoms
    mol2%num_bonds = nbonds
    mol2%num_subst = nsub

    allocate(mol2%atoms(natoms), mol2%order(natoms), mol2%connections(natoms), &
         mol2%bonds(nbonds), mol2%substructures(nsub), stat=lerror)
    
    do subst_id = 1, nsub
       mol2%substructures(subst_id)%root_atom = 0
    end do

    if (present(error)) error = lerror

  end subroutine init_mol2


  ! Assignment from one mol2 object to another

  subroutine assign_from_mol2(mol2l,mol2r)

    ! This allows the user to copy the data from one mol2 object to another 

    ! Left hand side of the assignment
    type (mol2_object), intent(inout) :: mol2l
    ! Right hand side of the assignment
    type (mol2_object), intent(in)  :: mol2r

    integer :: i, j, error

    call init(mol2l, mol2r%num_atoms, mol2r%num_bonds, mol2r%num_subst, error)

    mol2l%mol_name    = mol2r%mol_name
    mol2l%num_atoms   = mol2r%num_atoms
    mol2l%num_bonds   = mol2r%num_bonds
    mol2l%num_subst   = mol2r%num_subst
    mol2l%mol_type    = mol2r%mol_type
    mol2l%charge_type = mol2r%charge_type
    mol2l%mol_comment = mol2r%mol_comment

    if (error /= 0) then
       write(stderr,*) 'MOL2_CLASS :: Error allocating memory for mol2 object'
       stop
    end if

    mol2l%order = mol2r%order
    do i = 1, mol2l%num_atoms
       mol2l%atoms(i) = mol2r%atoms(i)
       j = push(mol2l%connections(i)%connectivity,mol2r%connections(i)%connectivity)
    end do

    do i = 1, mol2l%num_bonds
       mol2l%bonds(i) = mol2r%bonds(i)
    end do
    
    do i = 1, mol2l%num_subst
       mol2l%substructures(i) = mol2r%substructures(i)
    end do

    mol2l%cell      = mol2r%cell
    mol2l%space_grp = mol2r%space_grp
    mol2l%setting   = mol2r%setting

  end subroutine assign_from_mol2


  ! Internal routines

  integer function bond_type_from_string(string)

    character(len=2), intent(in) :: string

    ! Local variables
    integer :: i, result(1)

    result = pack((/(i, i=1,size(bond_type_lookup))/),bond_type_lookup==string)

    bond_type_from_string = result(1)
    
  end function bond_type_from_string

  function bond_type_from_integer(int)

    integer, intent(in) :: int

    ! Function result
    character(len=2) :: bond_type_from_integer

    if ((int < 0) .OR. (int > size(bond_type_lookup))) then
       bond_type_from_integer = 'un'
    else
       bond_type_from_integer = bond_type_lookup(int)
    end if

  end function bond_type_from_integer

  
  ! Reading mol2 files ...


  subroutine assign_from_filename(mol2,file)

    ! This allows the user to specify a mol2 filename in an assignment
    ! statement 

    ! This is the object on the left hand side of the assignment
    type (mol2_object), intent(out) :: mol2
    
    ! mol2 filename
    character(len=*), intent(in) :: file

    mol2 = open_then_read(file)

  end subroutine assign_from_filename

  function open_then_read(file)

    ! This is a small wrapper that allows the user to specify a mol2 
    ! filename which this routine will open for them, read the mol2 
    ! file and then close it again.

    ! Return type of function
    type (mol2_object) :: open_then_read
    
    ! mol2 filename
    character(len=*), intent(in) :: file

    ! Internal variables
    logical :: opened_ok
    integer :: unit

    ! Open the mol2 file
    unit = open(file, opened_ok, status='old')

    if (opened_ok) then
       open_then_read = readmol2(unit)
       close(unit)
    else
       write (stderr,*)'MOL2_CLASS :: Unable to find file: ',file
       stop
    end if

  end function open_then_read

  function readmol2 (unit)

    ! Main mol2 parsing function. It will read from an already opened
    ! unit number and store the results in a mol2 object which it returns

    ! Return type of function
    type(mol2_object) :: readmol2

    ! Unit number
    integer, intent(in) :: unit

    ! Local variables
    integer :: i, dummy, atom_id, bond_id, subst_id, next, error, atom_id_2, bond_type, root_atom
    character(len=80) :: record, string
    character(len=2)  :: bond_type_string
    logical :: end_of_file
    
    do
       ! Read through the mol2 file until we reach the end 
       call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)
       if (end_of_file) exit

       ! Get the first string in the record
       next = 1
       string = gettext(record,next)

       ! Compare this string to tripos commands we are interested in
       select case (.ucase.(trim(string)))
       case ('@<TRIPOS>MOLECULE') 

          ! the title line and the number of atoms and bonds
          call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)
          readmol2%mol_name = trim(adjustl(record))

          ! The number of atoms, bonds and substructures
          call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)
          read (record,*) readmol2%num_atoms, readmol2%num_bonds, readmol2%num_subst

          ! The molecule type
          call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)
          read (record,*)  readmol2%mol_type

          ! The charge type
          call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)
          read (record,*) readmol2%charge_type

          ! Status bits .. ignore but check if we have a place holder
          ! string (four asterisks, i.e. "****"), and if so, we try
          ! and continue on and read a comment
          call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)

          if (record == '****') then

             ! The comment string
             call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)
             readmol2%mol_comment = trim(adjustl(record))

          end if

          call init(readmol2, readmol2%num_atoms, readmol2%num_bonds, readmol2%num_subst, error)

          if (error /= 0) then
             write(stderr,*) 'MOL2_CLASS :: Error allocating memory for mol2 object'
             stop
          end if

!!$          ! Allocate the memory required to store the atom positions
!!$          allocate(readmol2%atoms(readmol2%num_atoms),stat=error)
!!$
!!$          if (error /= 0) then
!!$             write(stderr,*) 'MOL2_CLASS :: Cannot allocate memory for atom records'
!!$             stop
!!$          end if
!!$    
!!$          ! Allocate the memory required to ordering information
!!$          allocate(readmol2%order(readmol2%num_atoms),stat=error)
!!$
!!$          if (error /= 0) then
!!$             write(stderr,*) 'MOL2_CLASS :: Cannot allocate memory for order array'
!!$             stop
!!$          end if
!!$    
!!$          ! Allocate the memory required to store the substructure data
!!$          allocate(readmol2%substructures(readmol2%num_subst),stat=error)
!!$          do subst_id = 1,readmol2%num_subst
!!$             nullify(readmol2%substructures(subst_id)%order)
!!$             readmol2%substructures(subst_id)%root_atom = 0
!!$             readmol2%substructures(subst_id)%last_atom = 0
!!$          end do
!!$
!!$          if (error /= 0) then
!!$             write(stderr,*) 'MOL2_CLASS :: Cannot allocate memory for substructures'
!!$             stop
!!$          end if
    
       case ('@<TRIPOS>ATOM') ! the atom names and coordinates

          do i = 1, readmol2%num_atoms

             call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)

             ! Grab the atom id -- should be the same as i if they are numbered
             ! sequentially, but just in case ..
             read (record,*) atom_id

             ! Make a note of the order that the atoms appear in
             readmol2%order(i) = atom_id

             ! Grab the atom name
             next = 1
             readmol2%atoms(atom_id)%atom_name = getword (record,next)

             ! The next three numbers after the atom name are the coordinates
             read (record(next:),*)  readmol2%atoms(atom_id)%coords

             ! The next word is the atom type
             readmol2%atoms(atom_id)%atom_type = getword (record,next)

             ! After that the next integer is the substructure id
             ! read (record(next:),*)  readmol2%atoms(atom_id)%subst_id
             read (record(next:),*)  subst_id

             readmol2%atoms(atom_id)%subst_id = subst_id

             ! Here we fill in the root atom information so that we can
             ! do a consistency check later on when we red in the substructure
             ! records
             if ((atom_id < readmol2%substructures(subst_id)%root_atom) .or. &
                  (readmol2%substructures(subst_id)%root_atom == 0)) &
                  readmol2%substructures(subst_id)%root_atom = atom_id

          end do

       case ('@<TRIPOS>BOND') ! bond list

!!$          ! Allocate the memory required to store the bond information
!!$          allocate(readmol2%bonds(readmol2%num_bonds),stat=error)
!!$
!!$          if (error /= 0) then
!!$             write(stderr,*) 'We haz an herror! Check it!'
!!$             stop
!!$          end if
!!$
!!$          ! Allocate the memory required to store the connections information
!!$          allocate(readmol2%connections(readmol2%num_bonds),stat=error)
!!$
!!$          if (error /= 0) then
!!$             write(stderr,*) 'We haz an herror! Check it!'
!!$             stop
!!$          end if
    
          do i = 1, readmol2%num_bonds

             call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)

             ! Read the the bond id, the two atom ids and the bond_type
             read (record,*)  bond_id, atom_id, atom_id_2, bond_type_string

             bond_type = btype(bond_type_string)

             ! Save these values into the bonds section of the mol2 object
             readmol2%bonds(i)%conn = (/ atom_id, atom_id_2 /)
             readmol2%bonds(i)%bond_type = bond_type

             ! Make some 'derived data' from the bond information -- basically
             ! a list of atoms connected to each atom

             ! Push the second atom id onto the list of bonded atoms for 
             ! the first atom (we are ignoring the return value of the function)
             dummy = push (readmol2%connections(atom_id)%connectivity, atom_id_2)

             ! As above but vice-versa (pushing the atom id of the first atom on
             ! to the list of bonded atoms for the second atom)
             dummy = push (readmol2%connections(atom_id_2)%connectivity, atom_id)

          end do

          ! for each atom, sort its list of attached atoms
          do i = 1, readmol2%num_atoms
             if (associated(readmol2%connections(i)%connectivity)) then
                call sort(readmol2%connections(i)%connectivity)
             else
                !  print *,'Atom: ',i,' has no connections!'
             end if
          end do

       case ('@<TRIPOS>SUBSTRUCTURE') ! substructure info

          do i = 1, readmol2%num_subst

             call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)

             ! Grab the substructure id -- should be the same as i if they 
             ! are numbered sequentially, but just in case ..
             read (record,*) subst_id

             ! Grab the substructure name
             next = 1
             readmol2%substructures(subst_id)%subst_name = getword (record,next)

             ! The next number is the root atom of the substructure
             ! read (record(next:),*)  readmol2%substructures(subst_id)%root_atom
             read (record(next:),*)  root_atom
             if (root_atom .ne. readmol2%substructures(subst_id)%root_atom) then
                write(stderr,'(A,I0,A,I0,A,I0,A)') 'Error! The root atom (',root_atom, &
                     ') of substructure ', subst_id,' is inconsistent with the atom list (',&
                     readmol2%substructures(subst_id)%root_atom,')' 
                stop
             end if
             
          end do

          ! We want to know where a substructure finishes as well as where it starts,
          ! so we cycle through the all the substructures except for the last, taking 
          ! one from the root atom of the next substructure. The last substructure is
          ! treated separately after the loop has exited
          do i = 1, readmol2%num_subst - 1
             readmol2%substructures(i)%last_atom = readmol2%substructures(i+1)%root_atom - 1
          end do
          readmol2%substructures(readmol2%num_subst)%last_atom= readmol2%num_atoms

       case ('@<TRIPOS>CRYSIN') ! crystal information

          call read_buffer(unit, record, end_of_file, comment="#", removeblanks=.TRUE.)

          read (record,*) readmol2%cell, readmol2%space_grp, readmol2%setting

       end select
    
    end do
    
    ! do i = 1, size(readmol2%substructures)
    !    print *,readmol2%substructures(i)%root_atom,readmol2%substructures(i)%last_atom
    ! end do

  end function readmol2

  subroutine add_mol2(orig, addition)

    ! Add one mol2 object to another

    type (mol2_object), intent(inout) :: orig
    type (mol2_object), intent(in)    :: addition

    ! Local variables
    type (mol2_object) :: tmp
    integer :: error, i, j, index
    
    ! Make a copy of the original
    tmp = orig
    
    ! Now expand the original to incorporate the space needed
    call init(orig, atom_num(tmp) + atom_num(addition), &
         bond_num(tmp) + bond_num(addition), &
         sub_num(tmp) + sub_num(addition), error)

    if (error /=0) then
       write(stderr,*) 'MOL2_CLASS :: Error allocating memory in add_mol2'
       stop
    end if

    ! Copy back the atoms and connections from tmp
    index = 0
    do i = 1, tmp%num_atoms
       index = index + 1
       orig%atoms(index) = tmp%atoms(i)
       j = push(orig%connections(index)%connectivity,tmp%connections(i)%connectivity)
    end do

    ! Add the new atoms and connections from addition
    do i = 1, addition%num_atoms
       index = index + 1
       orig%atoms(index) = addition%atoms(i)
       orig%atoms(index)%subst_id = addition%atoms(i)%subst_id + tmp%num_subst
       j = push(orig%connections(index)%connectivity,addition%connections(i)%connectivity+tmp%num_atoms)
    end do

    ! Copy back the bonds from tmp
    index = 0
    do i = 1, tmp%num_bonds
       index = index + 1
       orig%bonds(index) = tmp%bonds(i)
    end do
    
    ! Add the new bonds from addition
    do i = 1, addition%num_bonds
       index = index + 1
       orig%bonds(index) = addition%bonds(i)
       orig%bonds(index)%conn = orig%bonds(index)%conn + tmp%num_atoms
    end do

    ! Copy back the substructure info from tmp
    index = 0
    do i = 1, tmp%num_subst
       index = index + 1
       orig%substructures(index) = tmp%substructures(i)
    end do

    ! Add the new substructure info from addition
    do i = 1, addition%num_subst
       index = index + 1
       orig%substructures(index) = addition%substructures(i)
       orig%substructures(index)%root_atom = orig%substructures(index)%root_atom + tmp%num_atoms
       orig%substructures(index)%last_atom = orig%substructures(index)%last_atom + tmp%num_atoms
    end do
    
    orig%order = (/ tmp%order, addition%order + tmp%num_atoms /)

    call destroy(tmp)

  end subroutine add_mol2


  ! Writing mol2 files ...


  subroutine open_then_print(file, mol2)

    ! This is a small wrapper that allows the user to specify a filename
    ! to print to rather than a unit number. This routine will open the
    ! file for them, print the mol2 object and then close it again.

    type (mol2_object) :: mol2

    ! mol2 filename
    character(len=*), intent(in) :: file

    ! Internal variables
    logical :: opened_ok
    integer :: unit

    ! Open the zmatrix file
    unit = open(file, opened_ok)

    if (opened_ok) then
       call print(unit,mol2)
       close(unit)
    else
       write (stderr,*)'MOL2_CLASS :: Unable to open file: ',file
       stop
    end if

  end subroutine open_then_print

  subroutine print_mol2 (unit, mol2)

    ! Output a mol2 object in a standard format

    type (mol2_object) :: mol2

    integer, intent(in) :: unit

    ! Local variables
    integer :: i, atom_id

    write(unit,'(A)') '@<TRIPOS>MOLECULE' 
    write(unit,'(A)') mol2%mol_name
    write(unit,'(3I5)') mol2%num_atoms, mol2%num_bonds, mol2%num_subst
    write(unit,'(A)') mol2%mol_type
    write(unit,'(A)') mol2%charge_type
    if (mol2%mol_comment /= '') then 
       ! We only write out the comment string if we have to, and so we
       ! need a 'blank' value holder for the mol2 status bits, which we
       ! know nothing about
       write(unit,'(A)') '****'
       write(unit,'(A)') mol2%mol_comment
    end if
    write(unit,*) 

    write(unit,'(A)') '@<TRIPOS>ATOM' 
    do i = 1, mol2%num_atoms
       atom_id = mol2%order(i)
       write(unit,'(I5,1X,A,3F12.6,2X,A,I4,2X,A)') atom_id, &
            mol2%atoms(atom_id)%atom_name, & 
            mol2%atoms(atom_id)%coords, &
            mol2%atoms(atom_id)%atom_type, &
            mol2%atoms(atom_id)%subst_id, &
            mol2%substructures(mol2%atoms(atom_id)%subst_id)%subst_name
    end do

    write(unit,'(A)') '@<TRIPOS>BOND'
    do i = 1, mol2%num_bonds
       write(unit,'(3I5,A6)') i, mol2%bonds(i)%conn, btype(mol2%bonds(i)%bond_type)
    end do

    write(unit,'(A)') '@<TRIPOS>SUBSTRUCTURE'
    do i = 1, mol2%num_subst
       write(unit,'(I5,2X,A,I5)') i, mol2%substructures(i)%subst_name, mol2%substructures(i)%root_atom
    end do

    ! Write out crystallography information if it is available
    if (any(mol2%cell /= 0.0)) then 
       write(unit,'(A)') '@<TRIPOS>CRYSIN'
       write(unit,'(3F10.5,3F10.5,2X,2I3)') mol2%cell, mol2%space_grp, mol2%setting
    end if

  end subroutine print_mol2

  subroutine print_mol2_to_stdout (mol2)

    ! Wrapper for print_mol2 when we do not specify a unit
    ! number (previously had unit as optional in above subroutine,
    ! but wanted the ordering to remain, which requires the user to
    ! utilise a keyword for the mol2 argument -- simpler to just
    ! make this wrapper)

    type (mol2_object), intent(in) :: mol2

    ! Print the mol2 to stdout
    call print(stdout,mol2)

  end subroutine print_mol2_to_stdout


  ! function coerce_mol2_to_cartesian (mol2, status)

    ! type (mol2_object), intent(in) :: mol2
    ! type (cartesian_entry), pointer, dimension(:) :: coerce_mol2_to_cartesian 

    ! integer, optional, intent(in)  :: subst_id
    ! logical, optional, intent(out) :: status

    ! integer :: error

    ! Grab the necessaries out of a mol2 object to make a generic cartesian
    ! object -- could do this automatically, i.e. have the internal mol2 rep
    ! be a cartesian object .. we'll see

    ! allocate(coerce_mol2_to_cartesian(), stat=error)

    ! if (error /= 0) then
    !    write(stderr,*) 'mol2_class :: coerce_mol2_to_cartesian :: Error allocating memory!'
    !    stop
    ! end if

  ! end function coerce_mol2_to_cartesian


  ! Relational operators

  ! Equals

  logical function mol2_eq_mol2 (mol2l, mol2r) result(equal)

    ! Input variables
    type (mol2_object), intent(in) :: mol2l, mol2r

    ! Local variables
    integer :: i, atom_idl, atom_idr

    ! Initialise the result of the function to false so that we can
    ! simply return at any time if we find something that dosen't
    ! agree between the two mol2 objects
    equal = .FALSE.

    if (mol2l%mol_name /= mol2r%mol_name) return
    if (mol2l%num_atoms /= mol2r%num_atoms) return
    if (mol2l%num_bonds /= mol2r%num_bonds) return
    if (mol2l%num_subst /= mol2r%num_subst) return
    if (mol2l%mol_type /= mol2r%mol_type) return
    if (mol2l%charge_type /= mol2r%charge_type) return
    if (mol2l%mol_comment /= mol2r%mol_comment) return

    do i = 1, mol2l%num_atoms
       atom_idl = mol2l%order(i)
       atom_idr = mol2r%order(i)
       if (atom_idl /= atom_idr) return
       if (mol2l%atoms(atom_idl)%atom_name /= mol2r%atoms(atom_idr)%atom_name) return
       if (any(mol2l%atoms(atom_idl)%coords /= mol2r%atoms(atom_idr)%coords)) return
       if (mol2l%atoms(atom_idl)%atom_type /= mol2r%atoms(atom_idr)%atom_type) return
       if (mol2l%atoms(atom_idl)%subst_id /= mol2r%atoms(atom_idr)%subst_id) return
       ! if (mol2l%connections(atom_idl) /= mol2r%connections(atom_idr)) return
    end do

    do i = 1, mol2l%num_bonds
       if (any(mol2l%bonds(i)%conn /= mol2r%bonds(i)%conn)) return
       if (mol2l%bonds(i)%bond_type /= mol2r%bonds(i)%bond_type) return
    end do
    
    do i = 1, mol2l%num_subst
       if (mol2l%substructures(i)%subst_name /= mol2r%substructures(i)%subst_name) return
       if (mol2l%substructures(i)%root_atom /= mol2r%substructures(i)%root_atom) return
    end do

    if (any(mol2l%cell /= mol2r%cell)) return
    if (mol2l%space_grp /= mol2r%space_grp) return
    if (mol2l%setting /= mol2r%setting) return

    ! We have passed the entire battery of tests, so the two
    ! mol2 objects are equal
    equal = .TRUE.

  end function mol2_eq_mol2

  logical function mol2_neq_mol2 (mol2l, mol2r) result(equal)

    ! Input variables
    type (mol2_object), intent(in) :: mol2l, mol2r

    equal = (.not. (mol2l == mol2r))

  end function mol2_neq_mol2

end module mol2_class
