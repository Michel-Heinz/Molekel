! Created by michel on 04/12/22.

module integralsMain

    use nuclearAttractionIntegralJoshuaGoings
    use kinds
    use global
    use PGF
    use atomparams
    use overlapJoshuaGoings
    use nuclearRepulsion
    use kineticEJoshuaGoings
    use electronRepulsionIntegralsJoshuaGoings
    use dipoleIntegrals

    contains

    subroutine getIntegrals(atoms, geo, S_, T_, V_, int2e, nucRep)
        implicit none
        class(atom), dimension(:), intent(in) :: atoms
        class(geom), intent(in) :: geo
        real(r8), dimension(:,:), intent(inout), allocatable :: S_, T_, V_
        real(r8), dimension(:), intent(inout), allocatable :: int2e
        real(r8), intent(inout) :: nucRep
        integer :: i, numberOfOrbs, order
        real(r8), dimension(3) :: CCoords

        CCoords = 0._8
        order = 1

        numberOfOrbs = 0
        do i = 1, SIZE(atoms)
            numberOfOrbs = numberOfOrbs + SIZE(atoms(i)%contractedGaussFunctions)
        end do

        allocate(S_(numberOfOrbs,numberOfOrbs), T_(numberOfOrbs,numberOfOrbs), V_(numberOfOrbs,numberOfOrbs),&
                &int2e(numberOfUniqueInts(numberOfOrbs)))
        S_ = 0._r8
        T_ = 0._r8
        V_ = 0._r8
        int2e = 0._r8
        S_ =  overlapMatrix(atoms, geo, SIZE(atoms), numberOfOrbs)
        T_ =  kineticEMatrix(atoms, geo, SIZE(atoms), numberOfOrbs)
        V_ =  nuclearAMatrix(atoms, geo, SIZE(atoms), numberOfOrbs)
        int2e =  eriArray(atoms, geo, SIZE(atoms), numberOfOrbs)
        nucRep = vNuc(atoms, geo)
        end subroutine getIntegrals

end module integralsMain