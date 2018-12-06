module fundamental_constants

  use precision

  implicit none

  public

  ! Pi -- reference http://www.lacim.uqam.ca/piDATA/pi.html
  ! Ratio of the circumference of a circle to it's diameter
  ! Can be computed as acos( -1.0d0 )
  real(kind=rd_kind), parameter :: pi = 3.14159265358979_rd_kind

  ! Radian (in degrees) = 180/pi degrees
  real(kind=rd_kind), parameter :: radian = 57.29577951308232088_rd_kind

  ! From here on these definitions were from http://physics.nist.gov/cuu/index.html

  ! Electron volt in J
  real(kind=rd_kind), parameter :: eV = 1.60217653e-19_rd_kind

  ! Planck's constant in J s
  real(kind=rd_kind), parameter :: h = 6.6260693e-34_rd_kind

  ! Planck's constant over 2 pi in J s
  real(kind=rd_kind), parameter :: hbar = 1.05457168e-34_rd_kind

  ! Planck's constant in eV s
  real(kind=rd_kind), parameter :: h_eV = 4.13566743e-15_rd_kind

  ! Planck's constant over 2 pi in eV s
  real(kind=rd_kind), parameter :: hbar_eV = 6.58211915e-16_rd_kind

  ! Speed of light in a vacuum m s-1
  real(kind=rd_kind), parameter :: c = 299792458_rd_kind

end module fundamental_constants 
