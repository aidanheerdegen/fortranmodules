module cartesian_class

  implicit none

  private

  !! $Log: cartesian_class.f90,v $
  !! Revision 1.1  2006/05/18 00:30:15  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: cartesian_class.f90,v 1.1 2006/05/18 00:30:15 aidan Exp $"

  type connectivity_entry

     ! private

     ! A connectivity array containing a unique list of atoms that this
     ! atom is bonded to -- variable length array
     integer, dimension(:), pointer :: connectivity => null()

  end type connectivity_entry

  type cartesian_entry

     ! private

     ! A description of an atom in cartesian coordinates

     ! x, y, z coordinates in real space (typically angstroms though no
     ! unit is specified)
     real :: coordinates(3)

     ! type (connectivity_entry) :: connectivity

  end type cartesian_entry

  interface operator(==)
     module procedure conn_eq_conn
  end interface

  interface operator(/=)
     module procedure conn_neq_conn
  end interface

  ! Public data types
  public :: connectivity_entry, cartesian_entry

  ! Public operators
  public :: operator(==), operator(/=)

contains

  logical function conn_eq_conn(conn1, conn2) result(same)

     type (connectivity_entry), intent(in) :: conn1, conn2

     same = .false.

     if (associated(conn1%connectivity) .neqv. associated(conn2%connectivity)) return
     if (size(conn1%connectivity) /= size(conn2%connectivity)) return
     if (any(conn1%connectivity /= conn2%connectivity)) return

     same = .true.

  end function conn_eq_conn

  logical function conn_neq_conn(conn1, conn2) result(different)

     type (connectivity_entry), intent(in) :: conn1, conn2

     different = .not. (conn1 == conn2)

   end function conn_neq_conn

end module cartesian_class
