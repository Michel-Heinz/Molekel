! Created by michel on 07.12.22.

module geometry

    use kinds
    use global
    use elements

    type geom
        real(r8), dimension(:,:), allocatable :: coordinates
        real(r8), dimension(:), allocatable :: centerOfMass, centerOfCharge
        integer, dimension(:), allocatable :: atomicNumber
        real(r8), dimension(:), allocatable :: mass

    contains
!TODO setter and getters
        procedure :: setGeom
        procedure :: calcCenterOfMass
        procedure :: calcCenterOfCharge

    end type geom

    contains
        subroutine setGeom(this, atomName, coordinates)
            class(geom), intent(inout) :: this
            real(r8), dimension(:,:), intent(in) :: coordinates
            character(len=2), dimension(:), intent(in) :: atomName
            integer :: numberOfAtoms, i

            numberOfAtoms = SIZE(atomName)
            allocate(this%coordinates(numberOfAtoms, 3), this%atomicNumber(numberOfAtoms), this%mass(numberOfAtoms))
            this%coordinates = coordinates
            do i = 1, numberOfAtoms
                call getElementInfo(atomName(i), this%atomicNumber(i), this%mass(i))
            end do

            this%centerOfMass = this%calcCenterOfMass()
            this%centerOfCharge = this%calcCenterOfCharge()
        end subroutine setGeom

        function calcCenterOfMass(this) result(com)
            class(geom), intent(in) :: this
            integer :: i, moleculeMass
            real(r8), dimension(3) :: com

            com = 0._r8
            moleculeMass = 0._r8
            do i = 1, SIZE(this%atomicNumber)
                com = com + this%mass(i) * this%coordinates(i,:)
                moleculeMass = moleculeMass + this%mass(i)
            end do
            com = com / moleculeMass
        end function calcCenterOfMass

        function calcCenterOfCharge(this) result(coc)
            class(geom), intent(in) :: this
            integer :: i, chargeOfAllNuclei
            real(r8), dimension(3) :: coc

            coc = 0._r8
            chargeOfAllNuclei = 0._r8
            do i = 1, SIZE(this%atomicNumber)
                coc = coc + this%atomicNumber(i) * this%coordinates(i,:)
                chargeOfAllNuclei = chargeOfAllNuclei + this%atomicNumber(i)
            end do
            coc = coc / chargeOfAllNuclei
        end function calcCenterOfCharge

end module geometry