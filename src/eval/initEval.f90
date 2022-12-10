! Created by michel on 07.12.22.

module initEval

    use parsing
    use machine
    use global
    use readinputfile, only: readinputparams
    use elements, only: getElementInfo
    use atomparams, only: create_atom, atom
    use geometry, only: setGeom, geom
    use errorAssert
    use RHF
    use outputParsing

    contains

    subroutine energyEvaluation()
        integer :: numberofarg, ierr, nLines, idx, nbl, iflag, method, maxiter, guess, numberofatoms, verbose, i
        !method: 1 = RHF
        !guess:  1 = HCore, 2 = GWH
        real(r8) :: precision
        class(atom), dimension(:), allocatable :: atoms
        class(geom), allocatable :: geo
        character(len=120) :: filename, token, basisset, guessStr
        character(len=MAXLEN), dimension(MAXLINES) :: allLines, blockLines
        character(len=2), dimension(:), allocatable :: atomsNames
        real(r8), dimension(:,:), allocatable :: coordinates

        numberofarg = 1
        call getarguments(numberofarg, filename, ierr)
        call assert(ierr==0, 'Molekel needs an inputfile (.in)')


        call readFile(filename, allLines, nLines)

        call printTitle('Input')
        do i = 1, nLines
            write(*,*) TRIM(allLines(i))
        end do

        verbose = 1

        idx = 1

        do
            call getNextBlock(allLines,nLines,idx,'$',')','!',MAXLINES,blockLines,nbl)
            if (nbl == 0) exit

            token = getToken(blockLines(1),'$','(')
            select case(token)
            case('general')
                call getinta(blockLines, nbl, 'verbose=', verbose, iflag)
                if (iflag == 1) then
                    verbose = 1
                else
                    call assert(verbose >0, 'Please input verbosity > 0!')
                end if
            case('method')
                if (finda(blockLines,nbl,'RHF')) then
                    method = 1
                else
                    call assert(.false., 'You need to specify an method (currently only RHF)')
                end if
                call getstra(blockLines, nbl, 'basisset=', basisset, iflag)
                call assert(iflag == 0, 'You need to specifiy a basisset.')
                call getinta(blockLines, nbl, 'maxiter=', maxiter, iflag)
                if (iflag == 1) maxiter = 100
                call getdbla(blockLines, nbl, 'precision=', precision, iflag)
                if (iflag == 1) precision = 1e-8_r8
                call getstra(blockLines, nbl, 'guess=', guessStr, iflag)
                if (iflag == 1 .or. guessStr == 'HCore') then
                    guess = 1
                elseif (guessStr == 'GWH') then
                    guess = 2
                else
                    call assert(.false., 'Please specifiy a valid guess (either HCore or GWH))!')
                end if
            case('geometry')
                call assert(TRIM(blockLines(nbl)) == ')', 'The geometry section must end with a ")" in the next &
                        &line after all the atoms')
                call readGeometry(atomsNames, coordinates, blockLines(2:nbl-1), nbl)
                numberofatoms = SIZE(atomsNames)
            end select
        end do

        allocate(atoms(numberofatoms))
        allocate(geo)
        call setGeom(geo, atomsNames, coordinates)

        do i = 1, numberofatoms
            call create_atom(atoms(i), basisset, atomsNames(i))
        end do

        Energy = doRHF(atoms, geo, 0, 100, 1e-8_r8, guess, verbose)


    end subroutine energyEvaluation

    subroutine readGeometry(atoms, coordinates, lines, nbl)
        implicit none
        character(len=2), dimension(:), allocatable, intent(out) :: atoms
        integer, intent(in) :: nbl
        real(r8), dimension(:,:), allocatable, intent(out) :: coordinates
        character(len=MAXLEN), dimension(nbl - 1), intent(in) :: lines
        character(len=120), dimension(:), allocatable :: atomAndCoords
        integer :: i, numberofatoms
        real(r8) :: x, y, z

        numberofatoms = nbl-2

        !allocating atoms and numberofatoms
        allocate(atoms(numberofatoms))
        allocate(coordinates(numberofatoms,3))

        !getting the atoms and the coordinates
        do i = 1, numberofatoms
            atomAndCoords = Py_split(lines(i))
            atoms(i) = TRIM(atomAndCoords(1))
            call strtoreal8(TRIM(atomAndCoords(2)), x)
            call strtoreal8(TRIM(atomAndCoords(3)), y)
            call strtoreal8(TRIM(atomAndCoords(4)), z)
            coordinates(i,1) = angstoau*x
            coordinates(i,2) = angstoau*y
            coordinates(i,3) = angstoau*z
        end do

    end subroutine readGeometry

end module initEval