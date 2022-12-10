! Created by michel on 09.12.22.

module dipolMoment

    use global
    use kinds
    use dipoleIntegrals

    contains

    function getDipolMoments(atoms, geo, numberOfOrbs, D) result(dipolMoments)
        integer, intent(in) :: numberOfOrbs
        real(r8), dimension(:,:), intent(in) :: D
        class(atom), dimension(:), intent(in) :: atoms
        class(geom), intent(in) :: geo
        real(r8), dimension(4) :: dipolMoments
        real(r8), dimension(numberOfOrbs,numberOfOrbs,3) :: Dip
        integer :: i, A

        Dip = dipoleMatrices(atoms, geo, SIZE(Atoms), numberOfOrbs, 1)
        dipolMoments = 0._r8
        do A = 1, SIZE(geo%atomicnumber)
            dipolMoments(:3) = dipolMoments(:3) + geo%atomicnumber(A) * (geo%coordinates(A,:) - geo%centerofcharge)
        end do

        do i = 1, numberOfOrbs
            do j = 1, numberOfOrbs
                dipolMoments(:3) = dipolMoments(:3) - D(i, j) * Dip(i, j, :)
            end do
        end do

        dipolMoments = dipolMoments * auToDebye
        dipolMoments(4) = NORM2(dipolMoments(:3))

        do i = 1, 4
            if (ABS(DipolMoments(i)) < 1.e-6_r8) then
                DipolMoments(i) = 0._r8
            end if
        end do
    end function getDipolMoments
end module dipolMoment