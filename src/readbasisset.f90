module readbasisset

    use parsing
    use global
    use kinds

    contains
    subroutine readbasissetparams(basisset, atom, angularmomentum, numberofpgf, exponents, coefficients)
        implicit none
        character(len=*), intent(in) :: basisset
        character(len=*), intent(in) :: atom
        integer, dimension(:), allocatable, intent(out) :: angularmomentum
        integer, dimension(:), allocatable, intent(out) :: numberofpgf
        real(r8), dimension(:), allocatable, intent(out) :: exponents
        real(r8), dimension(:), allocatable, intent(out) :: coefficients

        integer :: ios, numberofangularmomentums, n, numpgf, lenpath, numOfSP, temp, i, sumcgfsSP
        integer, dimension(6) :: numOfAngMomArray
        logical :: is_atom = .false.
        logical :: pIsCor = .false.
        character(len=200) :: val1, val2, val3
        character(len=500) :: pathtobasisset
        real(r8), dimension(:), allocatable :: tempCoeArray
        real(r8) :: exp, coe, tempCoe

        lenpath = len(trim(path))
        pathtobasisset = path(1:lenpath) // '/basissets/' // trim(basisset) // '.basis'
        numberofangularmomentums = 0
        numOfAngMomArray = 0
        numOfSP = 0
        open(unit = 1, file=pathtobasisset, iostat=ios)

        !get the number of angularmomentums
        is_atom = .false.
        do while (.not. is_atom)
            read(1,*) val1
            if(val1=='****') then
                read(1,*) val1
                if(val1==atom) then
                    is_atom = .true.
                    do while (.true.)
                        read(1,*) val1, val2
                        if(val1 == 'S' .or. val1 == 'P' .or. val1 == 'D' .or. val1 == 'G' .or. val1 == 'F') then
                            numberofangularmomentums = numberofangularmomentums + 1
                        end if
                        if(val1 == 'SP') then
                            numberofangularmomentums = numberofangularmomentums + 2
                            numOfSP = numOfSP + 1
                        end if
                        if(val1 == '****') exit
                    end do
                end if
            end if
        end do


        !allocate the arrays angularmomentum and numberofpgf#
        allocate(angularmomentum(numberofangularmomentums))
        allocate(numberofpgf(numberofangularmomentums))
        angularmomentum = 0
        numberofpgf = 0

        !get the angularmomentums and the number of primitive gauss functions (pgf)
        rewind(unit = 1)
        is_atom=.false.
        n = 1
        do while (.not. is_atom)
            read(1,*) val1
            if(val1=='****') then
                read(1,*) val1
                if(val1==atom) then
                    is_atom = .true.
                    do while (.true.)
                        read(1,*) val1, val2
                        if(val1 == 'S' .or. val1 == 'P' .or. val1 == 'D' .or. val1 == 'G' .or. val1 == 'F'&
                                & .or. val1 == 'SP') then
                            select case (val1)
                                case ('SP')
                                angularmomentum(n) = 0
                                angularmomentum(n + numOfSP) = 1
                                numOfAngMomArray(1) = numOfAngMomArray(1) + 1
                                call strtoint(val2, numpgf)
                                numberofpgf(n) = numpgf
                                numberofpgf(n + numOfSP) = numpgf
                                n = n + 1
                                case ('S')
                                angularmomentum(n) = 0
                                numOfAngMomArray(2) = numOfAngMomArray(2) + 1
                                call strtoint(val2, numpgf)
                                numberofpgf(n) = numpgf
                                n = n + 1
                                case ('P')
                                angularmomentum(n) = 1
                                numOfAngMomArray(3) = numOfAngMomArray(3) + 1
                                call strtoint(val2, numpgf)
                                numberofpgf(n) = numpgf
                                n = n + 1
                                case ('D')
                                angularmomentum(n + numOfSP) = 2
                                numOfAngMomArray(4) = numOfAngMomArray(4) + 1
                                call strtoint(val2, numpgf)
                                numberofpgf(n + numOfSP) = numpgf
                                n = n + 1
                                case ('F')
                                angularmomentum(n + numOfSP) = 3
                                numOfAngMomArray(5) = numOfAngMomArray(5) + 1
                                call strtoint(val2, numpgf)
                                numberofpgf(n + numOfSP) = numpgf
                                n = n + 1
                                case ('G')
                                angularmomentum(n + numOfSP) = 4
                                numOfAngMomArray(6) = numOfAngMomArray(6) + 1
                                call strtoint(val2, numpgf)
                                numberofpgf(n + numOfSP) = numpgf
                                n = n + 1
                            end select
                        end if
                        if(val1 == '****') then
                            exit
                        end if
                    end do
                end if
            end if
        end do
!        write(*,*) numOfAngMomArray
!        write(*,*) numberofpgf

        !allocating the array exponents and coefficients
        allocate(exponents(SUM(numberofpgf)))
        allocate(coefficients(SUM(numberofpgf)))
        exponents = 0._r8
        coefficients = 0._r8
        sumcgfsSP = SUM(numberofpgf(1 + numOfAngMomArray(1):numOfAngMomArray(1) * 2))

        !get the exponents and the coefficients
        rewind(unit = 1)
        is_atom=.false.
        pIsCor=.false.
        n = 1
        do while (.not. is_atom)
            read(1,*) val1
            if(val1=='****') then
                read(1,*) val1
                if(val1==atom) then
                    is_atom = .true.
                    if (numOfAngMomArray(1) == 0) then !check if SP is in Basisset
                        do while (.true.)
                            read(1,*) val1, val2
                            if(val1 == 'S' .or. val1 == 'P' .or. val1 == 'D' .or. val1 == 'G' .or. val1 == 'F') then
                                read(1,*) val1, val2
                                call strtoreal8(val1, exp)
                                call strtoreal8(val2, coe)
                                exponents(n) = exp
                                coefficients(n) = coe
                                n = n + 1

                            elseif(val1 == '****') then
                                exit

                            else
                                call strtoreal8(val1, exp)
                                call strtoreal8(val2, coe)
                                exponents(n) = exp
                                coefficients(n) = coe
                                n = n + 1

                            end if
                        end do
                    else
                        read(1,*) val1, val2
                        do while (.true.)
                            if(val1 == 'S') then
                                do while (.true.)
                                    read(1,*) val1, val2
                                    if(val1 == 'S' .or. val1 == 'P' .or. val1 == 'D' .or. val1 == 'G' .or. val1 == 'F'&
                                            & .or. val1 == 'SP' .or. val1 == '****') then
                                        exit
                                    end if
                                    call strtoreal8(val1, exp)
                                    call strtoreal8(val2, coe)
                                    exponents(n) = exp
                                    coefficients(n) = coe
                                    n = n + 1
                                end do
                            else if (val1 == 'SP') then
                                allocate(tempCoeArray(numberofpgf(2) * 2))
                                tempCoeArray = 0._r8
                                temp = 1
                                do while (.true.)
                                    read(1,*) val1, val2, val3
                                    if(val1 == 'S' .or. val1 == 'P' .or. val1 == 'D' .or. val1 == 'G' .or. val1 == 'F'&
                                            & .or. val1 == 'SP' .or. val1 == '****') then
                                        exit
                                    end if
                                    call strtoreal8(val1, exp)
                                    call strtoreal8(val2, coe)
                                    call strtoreal8(val3, tempCoe)
                                    exponents(n) = exp
                                    coefficients(n) = coe
                                    tempCoeArray(temp) = tempCoe
                                    temp = temp + 1
                                    n = n + 1
                                end do
                                do i = 1, temp - 1
!                                    write(*,*) n, temp - 1, n + i + sumcgfsSP - temp, n - temp + i
!                                    write(*,*) tempCoeArray
                                    exponents(n + sumcgfsSP - temp + i) = exponents(n - temp + i)
                                    coefficients(n + sumcgfsSP - temp + i) = tempCoeArray(i)
                                end do
                                deallocate(tempCoeArray)
                            else if (val1 == 'D' .or. val1 == 'G' .or. val1 == 'F') then
                                if (.not. pIsCor) then
                                    n = n + sumcgfsSP
                                    pIsCor = .true.
                                end if

                                do while (.true.)
                                    read(1,*) val1, val2
                                    if(val1 == 'S' .or. val1 == 'P' .or. val1 == 'D' .or. val1 == 'G' .or. val1 == 'F'&
                                            & .or. val1 == 'SP' .or. val1 == '****') then
                                        exit
                                    end if
                                    call strtoreal8(val1, exp)
                                    call strtoreal8(val2, coe)
                                    exponents(n) = exp
                                    coefficients(n) = coe
                                    n = n + 1
                                end do
                            end if
                            if (val1 == '****') then
                                exit
                            end if
                        end do
                    end if
                end if
            end if
        end do
        close(1)
    end subroutine readbasissetparams

end module readbasisset