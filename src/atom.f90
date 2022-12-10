module atomparams

    use global
    use kinds
    use PGF
    use readinputfile
    use readbasisset
    use elements
    use geometry

    type atom
        integer :: numberofGTO
        character(len=2) :: atomName
        class(cgfs), dimension(:), allocatable :: contractedGaussFunctions

        contains!TODO setter for atom and getter for geometry
!        procedure :: set_atom
        procedure :: get_contractedGaussFunctions
        procedure :: get_numberofGTO

    end type atom

    contains

!    subroutine set_atom(this, coordinates, atomname, atomicnumber, numberofGTO, contractedGaussFunctions)
!        class(atom), intent(inout) :: this
!        real(r8), dimension(3), intent(in) :: coordinates
!        character(len=2), intent(in) :: atomname
!        integer, intent(in) :: atomicnumber, numberofGTO
!        class(cgfs), dimension(:), intent(in) :: contractedGaussFunctions
!
!        integer :: numberofcgfs
!
!        numberofcgfs = SIZE(contractedGaussFunctions)
!        allocate(this%contractedGaussFunctions(numberofcgfs))
!
!        this%coordinates = coordinates
!        this%atomname = atomname
!        this%numberofGTO = numberofGTO
!        this%contractedGaussFunctions = contractedGaussFunctions
!    end subroutine set_atom

    function get_contractedGaussFunctions(this)
        implicit none
        class(atom), intent(in) :: this
        class(cgfs), dimension(:), allocatable :: get_contractedGaussFunctions
        allocate(get_contractedGaussFunctions(SIZE(this%contractedGaussFunctions)))
        get_contractedGaussFunctions = this%contractedGaussFunctions
    end function get_contractedGaussFunctions

    subroutine create_contractedGaussFunctions(contractedGaussFunctions, atomname)
        implicit none
        class(cgfs), dimension(:), allocatable, intent(out) :: contractedGaussFunctions
        character(len=2) :: atomname

    end subroutine create_contractedGaussFunctions

    function get_numberofGTO(this)
        implicit none
        class(atom), intent(in) :: this
        integer :: get_numberofGTO
        get_numberofGTO = this%numberofGTO

    end function get_numberofGTO

    subroutine create_atom(atom_, basisset, atomname)
        implicit none
        class(atom), intent(out) :: atom_
        character(len=2), intent(in) :: atomname
        character(len=*), intent(in) :: basisset

        integer, dimension(:), allocatable :: angularmomentum
        integer, dimension(:), allocatable :: numberofpgf
        real(r8), dimension(:), allocatable :: exponents
        real(r8), dimension(:), allocatable :: coefficients
        integer :: numberofprotons, numberofGTO
        real(r8) :: mass
        integer :: a, b, c, index, indexincrement, indexnext
        integer, dimension(:,:), allocatable :: ijk

        !getting all the data
        call readbasissetparams(basisset, atomname, angularmomentum, numberofpgf, exponents, coefficients)
        call getElementInfo(atomname, numberofprotons, mass)

        !ordering the data and creating GTOs
        numberofGTO = 0
        do a = 1, SIZE(angularmomentum)
            select case(angularmomentum(a))
                case(0)
                numberofGTO = numberofGTO + 1
                case(1)
                numberofGTO = numberofGTO + 3
                case(2)
                numberofGTO = numberofGTO + 6
                case(3)
                numberofGTO = numberofGTO + 10
                case default
                stop 'no implementation of G orbitals!'
            end select
        end do

        allocate(atom_%contractedGaussFunctions(numberofGTO))
        index = 1
        c = 1
        indexincrement = 1

        do a = 1, SIZE(angularmomentum)
            indexincrement = numberofpgf(a) - 1
            indexnext = index + indexincrement
            call create_ijk(angularmomentum(a), ijk)
            select case(angularmomentum(a))
                case(0)
                call atom_%contractedGaussFunctions(c)%set_cgfs(angularmomentum(a), ijk, exponents(index:indexnext),&
                        coefficients(index:indexnext))
                c = c + 1
                case(1)
                do b = 1, 3
                    call atom_%contractedGaussFunctions(c)%set_cgfs(angularmomentum(a), [ijk(b,:)],&
                            exponents(index:indexnext), coefficients(index:indexnext))
                    c = c + 1
                end do
                case(2)
                do b = 1, 6
                    call atom_%contractedGaussFunctions(c)%set_cgfs(angularmomentum(a), [ijk(b,:)],&
                            exponents(index:indexnext), coefficients(index:indexnext))
                    c = c + 1
                end do
                case(3)
                do b = 1, 10
                    call atom_%contractedGaussFunctions(c)%set_cgfs(angularmomentum(a), [ijk(b,:)],&
                            exponents(index:indexnext), coefficients(index:indexnext))
                    c = c + 1
                end do
            end select
            index = index + indexincrement + 1

        end do
!        write(*,*) 'gogogo'
!        do a = 1, numberofGTO
!            write(*,*) atom_%contractedGaussFunctions(a)%get_angularmomentum()
!            write(*,*) atom_%contractedGaussFunctions(a)%get_ijk()
!            write(*,*) atom_%contractedGaussFunctions(a)%get_exponents()
!            write(*,*) atom_%contractedGaussFunctions(a)%get_coefficients()
!        end do
!        write(*,*) 'gogogo'
        atom_%numberofGTO = numberofGTO
        atom_%atomName = atomname

    end subroutine create_atom

end module atomparams