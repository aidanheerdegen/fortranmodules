module statistics

  use sort_functions
  use precision
  use variable_array

  implicit none

  private

  !! $Log: statistics.f90,v $
  !! Revision 1.2  2006/03/14 03:45:13  aidan
  !! Added .isodd. and .iseven. operators.
  !!
  !! Revision 1.1  2005/10/21 04:51:37  aidan
  !! Initial revision
  !!

  ! Revision control software updates this character parameter.
  ! The 'ident' command can extract this version string from an
  ! object file or executable, which means one can identify which
  ! version of the module was used to compile it.
  character(len=*), parameter :: version = "$Id: statistics.f90,v 1.2 2006/03/14 03:45:13 aidan Exp aidan $"

  interface median
     module procedure median_qsort_integer, median_qsort_double, median_qsort_real ! iterativemedian
  end interface

  interface mean
     module procedure mean_integer, mean_real, mean_double
  end interface

  interface mode
     module procedure mode_integer, mode_real ! , mode_double
  end interface

  interface variance
     module procedure variance_integer, variance_real, variance_double
  end interface

  interface stddev
     module procedure stddev_integer, stddev_real, stddev_double
  end interface

  interface standard_deviation
     module procedure stddev_integer, stddev_real, stddev_double
  end interface

  interface operator(.isodd.)
     module procedure is_odd
  end interface

  interface operator(.iseven.)
     module procedure is_even
  end interface

  public :: mean, median, mode, variance, stddev, standard_deviation
  public :: is_even, is_odd, operator(.isodd.), operator(.iseven.)

contains

  logical function is_odd (number)
    
    integer, intent(in) :: number

    is_odd = (modulo(number, 2) == 1)

  end function is_odd
    
  logical function is_even (number)
    
    integer, intent(in) :: number

    is_even = .not. is_odd(number)

  end function is_even
    
  real(rd_kind) function mean_integer (array) result(mean)
      
    integer, intent(in) :: array(:)

    ! Local variables

    ! mean = real(sum(int(array,i8_kind)),rd_kind)/real(size(array),rd_kind)
    mean = real(sum(real(array,rd_kind)/real(size(array),rd_kind)))
    
  end function mean_integer
    
  real(rd_kind) function mean_real (array) result(mean)
      
    real, intent(in) :: array(:)

    ! Local variables

    ! mean = real(sum(int(array,i8_kind)),rd_kind)/real(size(array),rd_kind)
    mean = real(sum(real(array,rd_kind)/real(size(array),rd_kind)))
    
  end function mean_real

  real(rd_kind) function mean_double (array) result(mean)
      
    real(rd_kind), intent(in) :: array(:)

    ! Local variables

    ! mean = real(sum(int(array,i8_kind)),rd_kind)/real(size(array),rd_kind)
    mean = real(sum(real(array,rd_kind)/real(size(array),rd_kind)))
    
  end function mean_double

  real(rd_kind) function variance_integer (array, meanval) result(variance)
      
    integer, dimension(:), intent(in)   :: array
    real(rd_kind), intent(in), optional :: meanval
    
    real(rd_kind) :: arraymean

    if (present(meanval)) then
       arraymean = meanval
    else
       arraymean = mean(array)
    end if
    
    variance = sum((real(array,rd_kind)-arraymean)**2)/real(size(array)-1,rd_kind)
    
  end function variance_integer

  real(rd_kind) function variance_real (array, meanval) result(variance)
      
    real, dimension(:), intent(in)      :: array
    real(rd_kind), intent(in), optional :: meanval
    
    real(rd_kind) :: arraymean

    if (present(meanval)) then
       arraymean = meanval
    else
       arraymean = mean(array)
    end if
    
    variance = sum((real(array,rd_kind)-arraymean)**2)/real(size(array)-1,rd_kind)
    
  end function variance_real

  real(rd_kind) function variance_double (array, meanval) result(variance)
      
    real(rd_kind), dimension(:), intent(in) :: array
    real(rd_kind), intent(in), optional     :: meanval
    
    real(rd_kind) :: arraymean

    if (present(meanval)) then
       arraymean = meanval
    else
       arraymean = mean(array)
    end if
    
    variance = sum((real(array,rd_kind)-arraymean)**2)/real(size(array)-1,rd_kind)
    
  end function variance_double

  real(rd_kind) function stddev_integer (array) result(stddev)
      
    integer, dimension(:), intent(in) :: array
    
    stddev = sqrt(variance(array))
    
  end function stddev_integer

  real(rd_kind) function stddev_real (array) result(stddev)
      
    real, dimension(:), intent(in) :: array
    
    stddev = sqrt(variance(array))
    
  end function stddev_real

  real(rd_kind) function stddev_double (array) result(stddev)
      
    real(rd_kind), dimension(:), intent(in) :: array
    
    stddev = sqrt(variance(array))
    
  end function stddev_double

  real function iterativemedian(array,initial) result(median)

    integer, dimension(:), intent(in) :: array
    real, intent(in), optional        :: initial

    real :: upperbound, lowerbound

    integer :: asize, numover
    integer :: maximum, minimum, step, position

    logical, dimension(size(array)) :: overmask, undermask, maskbuffer
    
    asize = size(array)

    maximum = maxval(array)
    minimum = minval(array)

    upperbound = maximum
    lowerbound = minimum

    if (present(initial)) then
       median = initial
    else
       ! median = (maximum + minimum)/2
       median = mean(array)
    end if

    ! print *,median

    overmask = .false.
    undermask = .false.
    maskbuffer = .false.
    numover = 0

    do
       ! print *,median,step,numover,upperbound,lowerbound
       maskbuffer = overmask
       where (.not. (overmask .or. undermask)) maskbuffer = array > median
       ! print *,count(.NOT. (overmask .OR. undermask))
       ! print *,count(.NOT. (maskbuffer .OR. undermask))
       numover = count(maskbuffer)
       
       if (abs(numover*2 - asize) <= 1) then
          overmask = maskbuffer
          if ((numover*2 - asize) == 0) then 
             ! print *,'even'
             median = (maxval(array,mask=(.not. overmask)) + minval(array,mask=overmask))/2.
          else
             ! print *,'odd'
             median = maxval(array,mask=(.not. overmask))
          end if
          exit
       end if
       if (numover < asize/2) then
          overmask = maskbuffer
          if (median < upperbound) upperbound = median
          step = (median - lowerbound)/2
          median = median - step
       else
          undermask = .not. maskbuffer
          if (median > lowerbound) lowerbound = median
          step = (upperbound - median)/2
          median = median + step
       end if

       ! if (step == 0 .OR. abs(numover*2 - asize) <= 1000) then
       if (step == 0) then
          ! print *,abs(numover*2 - asize)

          ! We have reached a point where we must sort 
          ! to establish a median value.

          overmask = (array > upperbound)
          ! maskbuffer = (array < lowerbound)
          ! position = asize/2  - count(maskbuffer)
          position = asize/2 - count(undermask)

          if (modulo(asize,2) == 0) then
             ! print *,'even'
             ! Even
             ! print *,pack(array,mask=(.NOT. (overmask .OR. maskbuffer)))
             ! print *,position
             median = (select_value_by_position(pack(array,mask=(.not. (overmask .or. undermask))),position) &
                   + select_value_by_position(pack(array,mask=(.not. (overmask .or. undermask))),position+1)) / 2.
          else
             ! print *,'odd'
             ! Odd
             median = select_value_by_position(pack(array,mask=(.not. (overmask .or. undermask))),position)
          end if
          exit
          
       end if

    end do

  end function iterativemedian

  real function select_value_by_position(array, position) result(value)

    integer, intent(in) :: array(:), position

    integer :: local_copy(size(array)), asize

    local_copy = array
    asize = size(array)

    if (position > asize) then
       value = array(asize)
    end if
    
    call sort(local_copy)
    
    value = local_copy(position)

  end function select_value_by_position

  real function mediansel(inarray)

    integer, dimension(:), intent(in) :: inarray

    integer, dimension(size(inarray)) :: array

    logical :: lefttoright
    integer :: i, j, left, step, asize, threshold, remainder

    left = 1; step = 1; threshold = sqrt(real(size(inarray))) ! 100000

    lefttoright = .false.

    array = inarray

    asize = size(array)

    do while (asize > threshold)
       
       remainder = mod(asize,3)

       print *,asize, step, lefttoright, remainder

       lefttoright = .not. lefttoright

       if (lefttoright) then
          i = left
       else
          i = left + (3 + remainder)*step
       end if

       do j = 1, (asize/3 - 1)
          call triplet_adjust(array, i, step)
          i = i + 3 * step
       end do
       
       if (lefttoright) then
          i = left + step
       else
          i = left
          left = left + (1 + remainder)*step
       end if

       call selection_sort(array, i, 3 + remainder, step)

       if (remainder == 2) then
          if (lefttoright) then
             call swap(array,i+step,i+2*step)
          else
             call swap(array,i+2*step,i+3*step)
          end if
       end if

       step = 3 * step
       asize = asize / 3

    end do

    print *,left, asize, step
    call selection_sort(array, left, asize, step)
    print *,left, step, asize, left + step*((asize-1)/2)

    mediansel = array(left + step*((asize-1)/2))

  end function mediansel

  subroutine selection_sort(array, left, asize, step)
    
    ! This procedure sorts asize elements of array located at
    ! positions left, left + step, left + 2*step, ....

    integer, intent(inout) :: array(:)
    integer, intent(in)    :: left, asize, step 

    integer :: i, j, min

    do i = left, left + (asize - 1)*step, step
       min = i
       ! print *,'array(i) = ',array(i)
       do j = i + step, left + asize*step-1, step
          ! print *,i,j
          if (array(j) < array(min)) min = j
       end do
       call swap(array,i,min)
    end do

  end subroutine selection_sort

  subroutine triplet_adjust(array, i, step)
    
    integer, intent(inout) :: array(:)
    integer, intent(in)    :: i, step

    integer :: j, k

    j = i + step
    k = i + 2*step

    if (array(i) < array(j)) then
       if (array(k) < array(i)) then
          call swap(array,i,j)
       else if (array(k) < array(j)) then
          call swap(array,j,k)
       end if
    else
       if (array(i) < array(k)) then
          call swap(array,i,j)
       else if (array(k) > array(j)) then
          call swap(array,j,k)
       end if
    end if

  end subroutine triplet_adjust

  subroutine swap(array, i, j)
    
    integer, intent(inout) :: array(:)
    integer, intent(in)    :: i, j
    
    integer :: temp

    temp = array(i)
    array(i) = array(j)
    array(j) = temp

  end subroutine swap

  real(rd_kind) function median_qsort_double(array) result(xmed)

    ! Find the median of X(1), ... , X(N), using as much of the quicksort
    ! algorithm as is needed to isolate it.
    ! N.B. On exit, the array X is partially ordered.

    !     Latest revision - 26 November 1996
    implicit none

    real(rd_kind), intent(IN), dimension(:) :: array

    ! Local variables

    real(rd_kind), dimension(size(array)) :: x

    integer :: n
    real(rd_kind)    :: temp, xhi, xlo, xmax, xmin
    logical :: odd
    integer :: hi, lo, nby2, nby2p1, mid, i, j, k

    x = array
    n = size(array)

    nby2 = n / 2
    nby2p1 = nby2 + 1
    odd = .true.

    !     HI & LO are position limits encompassing the median.

    if (n == 2 * nby2) odd = .false.
    lo = 1
    hi = n
    if (n < 3) then
       if (n < 1) then
          xmed = 0.0
          return
       end if
       xmed = x(1)
       if (n == 1) return
       xmed = 0.5*(xmed + x(2))
       return
    end if

    !     Find median of 1st, middle & last values.

10  mid = (lo + hi)/2
    xmed = x(mid)
    xlo = x(lo)
    xhi = x(hi)
    if (xhi < xlo) then          ! Swap xhi & xlo
       temp = xhi
       xhi = xlo
       xlo = temp
    end if
    if (xmed > xhi) then
       xmed = xhi
    else if (xmed < xlo) then
       xmed = xlo
    end if

    ! The basic quicksort algorithm to move all values <= the sort key (XMED)
    ! to the left-hand end, and all higher values to the other end.

    i = lo
    j = hi
50  do
       if (x(i) >= xmed) exit
       i = i + 1
    end do
    do
       if (x(j) <= xmed) exit
       j = j - 1
    end do
    if (i < j) then
       temp = x(i)
       x(i) = x(j)
       x(j) = temp
       i = i + 1
       j = j - 1

       !     Decide which half the median is in.

       if (i <= j) GO TO 50
    end if

    if (.not. odd) then
       if (j == nby2 .and. i == nby2p1) GO TO 130
       if (j < nby2) lo = i
       if (i > nby2p1) hi = j
       if (i /= j) GO TO 100
       if (i == nby2) lo = nby2
       if (j == nby2p1) hi = nby2p1
    else
       if (j < nby2p1) lo = i
       if (i > nby2p1) hi = j
       if (i /= j) GO TO 100

       ! Test whether median has been isolated.

       if (i == nby2p1) return
    end if
100 if (lo < hi - 1) GO TO 10

    if (.not. odd) then
       xmed = 0.5*(x(nby2) + x(nby2p1))
       return
    end if
    temp = x(lo)
    if (temp > x(hi)) then
       x(lo) = x(hi)
       x(hi) = temp
    end if
    xmed = x(nby2p1)
    return

    ! Special case, N even, J = N/2 & I = J + 1, so the median is
    ! between the two halves of the series.   Find max. of the first
    ! half & min. of the second half, then average.

130 xmax = x(1)
    do k = lo, j
       xmax = max(xmax, x(k))
    end do
    xmin = x(n)
    do k = i, hi
       xmin = min(xmin, x(k))
    end do
    xmed = 0.5*(xmin + xmax)

    return
  end function median_qsort_double

  real function median_qsort_real(array) result(xmed)

    ! Find the median of X(1), ... , X(N), using as much of the quicksort
    ! algorithm as is needed to isolate it.
    ! N.B. On exit, the array X is partially ordered.

    !     Latest revision - 26 November 1996
    implicit none

    real, intent(IN), dimension(:) :: array

    ! Local variables

    real, dimension(size(array)) :: x

    integer :: n
    real    :: temp, xhi, xlo, xmax, xmin
    logical :: odd
    integer :: hi, lo, nby2, nby2p1, mid, i, j, k

    x = array
    n = size(array)

    nby2 = n / 2
    nby2p1 = nby2 + 1
    odd = .true.

    !     HI & LO are position limits encompassing the median.

    if (n == 2 * nby2) odd = .false.
    lo = 1
    hi = n
    if (n < 3) then
       if (n < 1) then
          xmed = 0.0
          return
       end if
       xmed = x(1)
       if (n == 1) return
       xmed = 0.5*(xmed + x(2))
       return
    end if

    !     Find median of 1st, middle & last values.

10  mid = (lo + hi)/2
    xmed = x(mid)
    xlo = x(lo)
    xhi = x(hi)
    if (xhi < xlo) then          ! Swap xhi & xlo
       temp = xhi
       xhi = xlo
       xlo = temp
    end if
    if (xmed > xhi) then
       xmed = xhi
    else if (xmed < xlo) then
       xmed = xlo
    end if

    ! The basic quicksort algorithm to move all values <= the sort key (XMED)
    ! to the left-hand end, and all higher values to the other end.

    i = lo
    j = hi
50  do
       if (x(i) >= xmed) exit
       i = i + 1
    end do
    do
       if (x(j) <= xmed) exit
       j = j - 1
    end do
    if (i < j) then
       temp = x(i)
       x(i) = x(j)
       x(j) = temp
       i = i + 1
       j = j - 1

       !     Decide which half the median is in.

       if (i <= j) GO TO 50
    end if

    if (.not. odd) then
       if (j == nby2 .and. i == nby2p1) GO TO 130
       if (j < nby2) lo = i
       if (i > nby2p1) hi = j
       if (i /= j) GO TO 100
       if (i == nby2) lo = nby2
       if (j == nby2p1) hi = nby2p1
    else
       if (j < nby2p1) lo = i
       if (i > nby2p1) hi = j
       if (i /= j) GO TO 100

       ! Test whether median has been isolated.

       if (i == nby2p1) return
    end if
100 if (lo < hi - 1) GO TO 10

    if (.not. odd) then
       xmed = 0.5*(x(nby2) + x(nby2p1))
       return
    end if
    temp = x(lo)
    if (temp > x(hi)) then
       x(lo) = x(hi)
       x(hi) = temp
    end if
    xmed = x(nby2p1)
    return

    ! Special case, N even, J = N/2 & I = J + 1, so the median is
    ! between the two halves of the series.   Find max. of the first
    ! half & min. of the second half, then average.

130 xmax = x(1)
    do k = lo, j
       xmax = max(xmax, x(k))
    end do
    xmin = x(n)
    do k = i, hi
       xmin = min(xmin, x(k))
    end do
    xmed = 0.5*(xmin + xmax)

    return
  end function median_qsort_real
  
  real function median_qsort_integer(array) result(xmed)

    ! Find the median of X(1), ... , X(N), using as much of the quicksort
    ! algorithm as is needed to isolate it.
    ! N.B. On exit, the array X is partially ordered.

    !     Latest revision - 26 November 1996
    implicit none

    integer, intent(IN), dimension(:) :: array

    ! Local variables

    integer, dimension(size(array)) :: x

    integer :: temp, xhi, xlo, xmax, xmin, n
    logical :: odd
    integer :: hi, lo, nby2, nby2p1, mid, i, j, k

    x = array

    n = size(array)

    nby2 = n / 2
    nby2p1 = nby2 + 1
    odd = .true.

    !     HI & LO are position limits encompassing the median.

    if (n == 2 * nby2) odd = .false.
    lo = 1
    hi = n
    if (n < 3) then
       if (n < 1) then
          xmed = 0.0
          return
       end if
       xmed = x(1)
       if (n == 1) return
       xmed = 0.5*(xmed + x(2))
       return
    end if

    !     Find median of 1st, middle & last values.

10  mid = (lo + hi)/2
    xmed = x(mid)
    xlo = x(lo)
    xhi = x(hi)
    if (xhi < xlo) then          ! Swap xhi & xlo
       temp = xhi
       xhi = xlo
       xlo = temp
    end if
    if (xmed > xhi) then
       xmed = xhi
    else if (xmed < xlo) then
       xmed = xlo
    end if

    ! The basic quicksort algorithm to move all values <= the sort key (XMED)
    ! to the left-hand end, and all higher values to the other end.

    i = lo
    j = hi
50  do
       if (x(i) >= xmed) exit
       i = i + 1
    end do
    do
       if (x(j) <= xmed) exit
       j = j - 1
    end do
    if (i < j) then
       temp = x(i)
       x(i) = x(j)
       x(j) = temp
       i = i + 1
       j = j - 1

       !     Decide which half the median is in.

       if (i <= j) GO TO 50
    end if

    if (.not. odd) then
       if (j == nby2 .and. i == nby2p1) GO TO 130
       if (j < nby2) lo = i
       if (i > nby2p1) hi = j
       if (i /= j) GO TO 100
       if (i == nby2) lo = nby2
       if (j == nby2p1) hi = nby2p1
    else
       if (j < nby2p1) lo = i
       if (i > nby2p1) hi = j
       if (i /= j) GO TO 100

       ! Test whether median has been isolated.

       if (i == nby2p1) return
    end if
100 if (lo < hi - 1) GO TO 10

    if (.not. odd) then
       xmed = 0.5*(x(nby2) + x(nby2p1))
       return
    end if
    temp = x(lo)
    if (temp > x(hi)) then
       x(lo) = x(hi)
       x(hi) = temp
    end if
    xmed = x(nby2p1)
    return

    ! Special case, N even, J = N/2 & I = J + 1, so the median is
    ! between the two halves of the series.   Find max. of the first
    ! half & min. of the second half, then average.

130 xmax = x(1)
    do k = lo, j
       xmax = max(xmax, x(k))
    end do
    xmin = x(n)
    do k = i, hi
       xmin = min(xmin, x(k))
    end do
    xmed = 0.5*(xmin + xmax)

    return
  end function median_qsort_integer

  integer function mode_integer(array) result(mode)
    
    integer, intent(in) :: array(:)

    ! Internal variables
    integer, dimension(1+maxval(array)-minval(array)) :: bins

    integer :: i, minarray, maxlocarray(1)

    bins = 0
    minarray = minval(array)
    
    ! print *,minarray,size(bins),size(array)

    do i = 1, size(array)
       bins(1+array(i)-minarray) = bins(1+array(i)-minarray) + 1
    end do

    maxlocarray = maxloc(bins)
    mode = maxlocarray(1) + minarray - 1

  end function mode_integer

  real function mode_real(array) result(mode)
    
    real, intent(in) :: array(:)

    ! Internal variables
    real, dimension(size(array)) :: local_array
    real, pointer :: unique_values(:)
    integer, allocatable :: mode_counts(:)

    integer :: i, nunique, ivalue, maxlocarray(1)

    local_array = array

    call sort(local_array)

    ! Make a list of unique values
    nunique = push(unique_values,unique(local_array))

    allocate(mode_counts(nunique))
    mode_counts = 0

    ivalue = 1
    do i = 1, size(array)
       if (local_array(i) /= unique_values(ivalue)) ivalue = ivalue + 1 
       mode_counts(ivalue) = mode_counts(ivalue) + 1 
    end do
    
    maxlocarray = maxloc(mode_counts)
    mode = unique_values(maxlocarray(1))

    deallocate(mode_counts)
    
    nunique = splice(unique_values,0)

  end function mode_real

!!$  real(rd_kind) function mode_double(array) result(mode)
!!$    
!!$    real(rd_kind), intent(in) :: array(:)
!!$
!!$    ! Internal variables
!!$    real(rd_kind), dimension(size(array)) :: local_array
!!$    real(rd_kind), pointer :: unique_values(:)
!!$    integer, allocatable :: mode_counts(:)
!!$
!!$    integer :: i, nunique, ivalue, maxlocarray(1)
!!$
!!$    call sort(local_array)
!!$
!!$    ! Make a list of unique values
!!$    nunique = push(unique_values,unique(local_array))
!!$
!!$    allocate(mode_counts(nunique))
!!$    mode_counts = 0
!!$
!!$    ivalue = 1
!!$    do i = 1, size(array)
!!$       if (local_array(i) /= unique_values(ivalue)) ivalue = ivalue + 1 
!!$       mode_counts(ivalue) = mode_counts(ivalue) + 1 
!!$    end do
!!$    
!!$    maxlocarray = maxloc(mode_counts)
!!$    mode = unique_values(maxlocarray(1))
!!$
!!$    deallocate(mode_counts)
!!$    
!!$    nunique = splice(unique_values,0)
!!$
!!$  end function mode_double

end module statistics
