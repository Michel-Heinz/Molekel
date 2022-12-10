! Created by michel on 08.12.22.

module outputParsing

    use parsing
    use kinds
    use global
    use atomparams

    contains

    function createOrbLabels(atoms) result(orbLabels)
        implicit none
        class(atom), dimension(:), intent(in) :: atoms
        character(len=9), dimension(:), allocatable :: orbLabels
        character(len=40) :: ps, ds, fs, num
        character(len=2):: num2
        integer :: i, j, s, p, d, f, l, m, n, numberOfOrbs, orb

        s = 1
        p = 1
        d = 1
        f = 1
        ps = 'PxPyPz'
        ds = 'Dx2Dy2Dz2DxyDxzDyz'
        fs = 'Fx3 Fy3 Fz3 Fxy2Fx2yFx2zFxz2Fyz2Fy2zFxyz'
        numberOfOrbs = 0

        do i = 1, SIZE(atoms)
            numberOfOrbs = numberOfOrbs + atoms(i)%numberofGTO
        end do
        allocate(orbLabels(numberOfOrbs))
! IT WAS LATE OKAY... SORRY FOR HARDCODE FIESTA
        orb = 1
        do i = 1, SIZE(atoms)
            write(num, '(I2)') i
            num2 = trim(num)
            do j = 1, atoms(i)%numberofGTO
                if (atoms(i)%contractedGaussFunctions(j)%angularmomentum == 0) then
                    orbLabels(orb) = num2//atoms(i)%atomName//' '//'S'
                    orb = orb + 1
                elseif (atoms(i)%contractedGaussFunctions(j)%angularmomentum == 1) then
                    if (p == 1) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ps(1:2)
                        p = p + 1
                        orb = orb + 1
                    elseif (p == 2) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ps(3:4)
                        p = p + 1
                        orb = orb + 1
                    elseif (p == 3) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ps(5:6)
                        p = 1
                        orb = orb + 1
                    end if
                elseif (atoms(i)%contractedGaussFunctions(j)%angularmomentum == 2) then
                    if (d == 1) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ds(1:3)
                        d = d + 1
                        orb = orb + 1
                    elseif (d == 2) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ds(4:6)
                        d = d + 1
                        orb = orb + 1
                    elseif (d == 3) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ds(7:9)
                        d = d + 1
                        orb = orb + 1
                    elseif (d == 4) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ds(10:12)
                        d = d + 1
                        orb = orb + 1
                    elseif (d == 5) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ds(13:15)
                        d = d + 1
                        orb = orb + 1
                    elseif (d == 6) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//ds(16:18)
                        d = 1
                        orb = orb + 1
                    end if
                elseif (atoms(i)%contractedGaussFunctions(j)%angularmomentum == 3) then
                    if (f == 1) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(1:4)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 2) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(5:8)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 3) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(9:12)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 4) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(13:16)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 5) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(17:20)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 6) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(21:24)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 7) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(25:28)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 8) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(29:32)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 9) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(33:36)
                        f = f + 1
                        orb = orb + 1
                    elseif (f == 10) then
                        orbLabels(orb) = num2//atoms(i)%atomName//' '//fs(37:40)
                        f = 1
                        orb = orb + 1
                    end if
                end if
            end do
        end do

    end function createOrbLabels

    subroutine printMatrix(matrix, columnLabels, rowLabels)
        implicit none
        real(r8), dimension(:,:), intent(in) :: matrix
        character(len=9), dimension(:), intent(in), optional :: rowLabels
        character(len=10), dimension(:), intent(in), optional :: columnLabels
        character(len=40) :: num
        character(len=9) :: numrl
        character(len=5) :: numcl
        integer :: i, j, l, cl, rl, cols, rows, colsidx, rowsidx

        rl = 0
        cl = 0
        cols = SIZE(matrix(1,:))
        rows = SIZE(matrix(:,1))
        if (PRESENT(rowLabels)) rl = 1
        if (PRESENT(columnLabels)) cl = 1

        colsidx = cols
        rowsidx = rows
        l = 1
        do while (colsidx > 0)
            if (colsidx >= 7) then
                if (cl == 0) then
                    write(*, '(A10)', Advance = 'No') ''
                    do i = l, 5 + l
                        write(num, '(I2)') i
                        numcl = trim(num)
                        write(*, '(A10)', Advance = 'No') numcl
                    end do
                    write(num, '(I2)') i
                    numcl = trim(num)
                    write(*, '(A10)') numcl
                else
                    write(*, '(A10)', Advance = 'No') ''
                    do i = l, 5 + l
                        write(*, '(A10)', Advance = 'No') columnLabels(i)
                    end do
                    write(*, '(A10)') columnLabels(i)
                end if
                do i = 1, rows
                    if (rl == 1) then
                        write(*, '(A9)', Advance = 'No') rowLabels(i)
                    else
                        write(num, '(I2)') i
                        numrl = trim(num)
                        write(*, '(A9)', Advance = 'No') numrl
                    end if
                    do j = l, 5 + l
                        write(*, '(F10.4,i5)', Advance = 'No') matrix(i, j)
                    end do
                    write(*, '(F10.4,i5)') matrix(i, j)
                end do
                write(*,*) ''
                write(*,*) ''
                l = l + 7
            else
                if (cl == 1) then
                    write(*, '(A10)', Advance = 'No') ''
                    do i = l, cols - 1
                        write(*, '(A10)', Advance = 'No') columnLabels(i)
                    end do
                    write(*, '(A6)') columnLabels(cols)
                else
                    write(*, '(A10)', Advance = 'No') ''
                    do i = l, cols - 1
                        write(num, '(I2)') i
                        numcl = trim(num)
                        write(*, '(A10)', Advance = 'No') numcl
                    end do
                    write(num, '(I2)') i
                    numcl = trim(num)
                    write(*, '(A10)') numcl
                end if
                do i = 1, rows
                    if (rl == 1) then
                        write(*, '(A9)', Advance = 'No') rowLabels(i)
                    else
                        write(num, '(I2)') i
                        numrl = trim(num)
                        write(*, '(A9)', Advance = 'No') numrl
                    end if
                    do j = l, cols - 1
                        write(*, '(F10.4)', Advance = 'No') matrix(i, j)
                    end do
                    write(*, '(F10.4)') matrix(i, cols)
                end do
            end if
            colsidx = colsidx - 7

        end do
    end subroutine printMatrix

    subroutine printTitle(title)
        character(len=*), intent(in) :: title
        character(len=40) :: left
        integer :: titleLen, l, i
        logical :: even

        write(*,*) ' '
        write(*,*) ' '
        titleLen = LEN(title)
        even = (MODULO(titleLen, 2) == 0)
        if (.not. even) then
            l = (MAXCHAR - titleLen) / 2
            left = repeat('*',l-4)
            write(*,*) repeat(' ',l-5)//repeat('*',titleLen+10)//repeat(' ',l-5)
            write(*,*) TRIM(left)//repeat(' ',4)//title//repeat(' ',4)//Trim(left)
            write(*,*) repeat(' ',l-5)//repeat('*',titleLen+10)//repeat(' ',l-5)
        else
            l = (MAXCHAR - titleLen-1) / 2
            left = repeat('*',l-4)
            write(*,*) repeat(' ',l-5)//repeat('*',titleLen+11)//repeat(' ',l-5)
            write(*,*) TRIM(left)//repeat(' ',4)//title//repeat(' ',5)//Trim(left)
            write(*,*) repeat(' ',l-5)//repeat('*',titleLen+11)//repeat(' ',l-5)
        end if
        write(*,*) ' '
    end subroutine printTitle

    subroutine printCenter(title)
        character(len=*), intent(in) :: title
        character(len=40) :: left
        integer :: titleLen, l, i
        logical :: even

        titleLen = LEN(title)
        even = (MODULO(titleLen, 2) == 0)
        if (.not. even) then
            l = (MAXCHAR - titleLen) / 2
            left = repeat(' ',l-4)
            write(*,*) repeat(' ',l-4)//repeat(' ',4)//title//repeat(' ',4)//repeat(' ',l-4)
        else
            l = (MAXCHAR - titleLen-1) / 2
            left = repeat(' ',l-4)
            write(*,*) repeat(' ',l-4)//repeat(' ',4)//title//repeat(' ',5)//repeat(' ',l-4)
        end if
    end subroutine printCenter

    subroutine printInt(text, int)
        implicit none
        character(len=*), intent(in) :: text
        integer, intent(in) :: int
        integer :: spaces
        character(len=40) :: num
        character(len=16) :: num2

        write(num, '(I5)') int
        num2 = trim(num)
        write(*,*) text//repeat(' ', MAXCHAR - len(text) - 16 - 34)//adjustr(num2)
    end subroutine printInt

    subroutine printFloat(text, float, unit)
        implicit none
        character(len=*), intent(in) :: text
        character(len=*), intent(in), optional :: unit
        real(r8), intent(in) :: float
        integer :: spaces
        character(len=40) :: num
        character(len=16) :: num2

        if (PRESENT(unit)) then
            write(num, '(F16.8)') float
            num2 = trim(num)
            write(*,*) text//repeat(' ', MAXCHAR - len(text) - 16 - 34)//num2//' '//unit
        else
            write(num, '(F16.8)') float
            num2 = trim(num)
            write(*,*) text//repeat(' ', MAXCHAR - len(text) - 16 - 34)//num2
        end if
    end subroutine printFloat

end module outputParsing