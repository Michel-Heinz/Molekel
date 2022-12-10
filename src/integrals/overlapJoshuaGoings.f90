! Created by michel on 19/05/22.
! Source https://joshuagoings.com/2017/04/28/integrals/ 21.5.22 20:46

module overlapJoshuaGoings

    use kinds
    use global
    use PGF
    use atomparams

    contains
    recursive function E(i,j,t,Qx,a,b) result(EResult)
!        Recursive definition of Hermite Gaussian coefficients.
!        Returns a float.
!        a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
!        b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
!        i,j: orbital angular momentum number on Gaussian 'a' and 'b'
!        t: number nodes in Hermite (depends on type of integral,
!                e.g. always zero for overlap integrals)
!        Qx: distance between origins of Gaussian 'a' and 'b'
!        Copied from above mentioned source
        implicit none
        integer, value :: i, j, t !pass by value with fortran 2003 faster for recusrive functions
        real(r8), value :: Qx, a, b
        real(r8) :: EResult, p, q

        p = a + b
        q = a * b / p
        EResult = 0._r8

        if ((t < 0) .or. (t > (i + j))) then
            !out of bounds for t
            EResult = 0.0_r8
        else if (i == j .and. i == t .and. i == 0) then
            !base case
            EResult = exp(-q * Qx * Qx)
        else if (j == 0) then
            !decrement index i
            EResult = (1 / (2 * p)) * E(i - 1, j, t-1, Qx, a, b) - (q * Qx / a) * E(i - 1, j, t, Qx, a, b) &
                    & + (t + 1) * E(i - 1, j, t+1, Qx, a, b)
        else
            !decrement index j
            EResult = (1 / (2 * p)) * E(i, j - 1, t-1, Qx, a, b) + (q * Qx / b) * E(i, j - 1, t, Qx, a, b) &
                    & + (t + 1) * E(i, j - 1, t+1, Qx, a, b)
        end if

    end function E

    function overlap(a, ijk1, Axyz, b, ijk2, Bxyz)
!        Evaluates overlap integral between two Gaussians
!        Returns a float.
!        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
!        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
!        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
!        for Gaussian 'a'
!        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
!        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
!        B:    list containing origin of Gaussian 'b'
!        Copied from above mentioned source
        integer, dimension(3), intent(in) :: ijk1, ijk2
        real(r8), dimension(3), intent(in) :: Axyz, Bxyz
        real(r8), intent(in) :: a, b
        real(r8) :: S1, S2, S3, overlap

        S1 = E(ijk1(1), ijk2(1), 0, (Axyz(1) - Bxyz(1)), a, b) !x
        S2 = E(ijk1(2), ijk2(2), 0, (Axyz(2) - Bxyz(2)), a, b) !y
        S3 = E(ijk1(3), ijk2(3), 0, (Axyz(3) - Bxyz(3)), a, b) !z
        overlap = S1 * S2 * S3 * (pi / (a + b)) ** (1.5_r8)

        end function overlap

    function S(a, aCoords, b, bCoords, ijk1, ijk2)
!        Evaluates overlap between two contracted Gaussians
!        Returns float.
!        Arguments:
!        a: contracted Gaussian 'a', BasisFunction object
!        b: contracted Gaussian 'b', BasisFunction object
!        Copied from above mentioned source
        class(cgfs), intent(in) :: a, b
        real(r8), dimension(3), intent(in) :: aCoords, bCoords
        integer, dimension(3), optional, intent(in) :: ijk1, ijk2
        integer(i8) :: ia, ib
        real(r8) :: primitiveOverlap, S
        integer, dimension(3) :: ijka, ijkb

        if (PRESENT(ijk1)) then
            ijka = ijk1
            ijkb = ijk2
        else
            ijka = a%ijk
            ijkb = b%ijk
        end if
        primitiveOverlap = 1.0_r8
        S = 0.0_r8
        do ia = 1, SIZE(a%coefficients)
            do ib = 1, SIZE(b%coefficients)
                !write(*,*) overlap(a%exponents(ia), a%ijk, aCoords, b%exponents(ib), b%ijk, bCoords)
                primitiveOverlap = overlap(a%exponents(ia), ijka, aCoords, b%exponents(ib), ijkb, bCoords)
                S = S + a%norms(ia) * b%norms(ib) * a%coefficients(ia) * b%coefficients(ib) * primitiveOverlap * &
                        & 1!a%normCGTO * b%normCGTO
            end do
        end do
        !write(*,*) a%normCGTO
        end function S

        function overlapMatrix(atoms, geo, sizeAtoms, sizeSMatrix) result(sMatrix)

            integer, intent(in)                                         :: sizeAtoms, sizeSMatrix
            class(atom), dimension(sizeAtoms), intent(in)               :: atoms
            class(geom), intent(in)                                     :: geo
            real(r8), dimension(sizeSMatrix,sizeSMatrix)                :: sMatrix
            integer                                                     :: atomAIndex, atomBIndex, aoAIndex, aoBIndex
            integer                                                     :: aoAIMax, aoBIMax, a, b

            sMatrix = 0._r8
            t = 0
            a = 1
            b = 1

            do atomAIndex = 1, sizeAtoms
                aoAIMax = SIZE(atoms(atomAIndex)%contractedGaussFunctions)
                do aoAIndex = 1, aoAIMax
                    do atomBIndex = 1, sizeAtoms
                        aoBIMax = SIZE(atoms(atomBIndex)%contractedGaussFunctions)
                        do aoBIndex = 1, aoBIMax
                            sMatrix(a, b) = S(atoms(atomAIndex)%contractedGaussFunctions(aoAIndex), &
                                    &geo%coordinates(atomAIndex,:), atoms(atomBIndex)%contractedGaussFunctions(aoBIndex),&
                                    &geo%coordinates(atomBIndex,:))
                            b = b + 1
                        end do
                    end do
                    b = 1
                    a = a + 1
                end do
            end do
        end function overlapMatrix


end module overlapJoshuaGoings