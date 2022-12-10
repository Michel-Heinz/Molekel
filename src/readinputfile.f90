module readinputfile
    use machine
    use parsing
    use kinds
    use global
    contains

    subroutine readinputparams(method, basisset, atoms, coordinates)

        character(len=*), intent(out) :: method, basisset
        character(len=*), dimension(:), allocatable, intent(out) :: atoms
        real(r8), dimension(:,:), allocatable, intent(out) :: coordinates

        character(len=500) :: filename, token, xs, ys, zs
        character(len=2) :: atom
        integer :: ierr, numberofarg, lenparam, numberofatoms, i, j
        real(r8) :: x, y, z

        !getting the file name from the commandline
        numberofarg = 1
        call getarguments(numberofarg, filename, ierr)
        if (ierr/=0) STOP 'Molekel needs an inputfile (.in)'
        !write(*,*) filename

        !getting the parameters from the different keywords
        open(unit = 1, file=filename, iostat=ios)
        do while(.true.)
            read(1,*) token
            if(token(1:6) == 'method') then
                lenparam = len(trim(token))
                method = token(8:lenparam)
            elseif(token(1:8) == 'basisset') then
                lenparam = len(trim(token))
                basisset = token(10:lenparam)
            elseif(token(1:4) == '****') then
                numberofatoms = 0
                !how many atoms are there
                do while(.true.)
                    read(1,*) token
                    if(token(1:4) == '****') exit
                    numberofatoms = numberofatoms + 1
                end do

                !allocating atoms and numberofatoms
                allocate(atoms(numberofatoms))
                allocate(coordinates(numberofatoms,3))

                !getting the atoms and the coordinates
                rewind(1)
                do while(.true.)
                    read(1,*) token
                    if(token(1:4) == '****') exit
                end do
                do i = 1, numberofatoms
                    read(1,*) atom, xs, ys, zs
                    atoms(i) = atom
                    call strtoreal8(xs, x)
                    call strtoreal8(ys, y)
                    call strtoreal8(zs, z)
                    coordinates(i,1) = angstoau*x
                    coordinates(i,2) = angstoau*y
                    coordinates(i,3) = angstoau*z
                end do
                exit
            end if
        end do
        close(1)
    end subroutine readinputparams

end module readinputfile