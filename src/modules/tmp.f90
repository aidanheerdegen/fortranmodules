program tmp

    character(len=20) :: b
    integer :: array(4)
    character(len=10) :: carray(4)

    b = '1,3,5,8'

    read (b,*) array
    print *,array

    read (b,*) carray
    print *,carray

    print *,count(transfer(b, 'a', len(b)) == ",")

end program tmp
