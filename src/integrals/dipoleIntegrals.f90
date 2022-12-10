! Created by michel on 07.12.22.

module dipoleIntegrals

    use kinds
    use global
    use PGF
    use atomparams
    use math
    use overlapJoshuaGoings

    contains

    recursive function M(e, t, PCx, p) result(Mval)
        !Helgaker p356-357
        integer, value :: e, t
        real(r8), value :: PCx, p
        real(r8) :: Mval
        Mval = 0._r8

        if (t > e .or. t < 0) then
            Mval = 0._r8
        elseif (t == e .and. e == 0) then
            Mval = Mval + SQRT(pi / p)
        else
            Mval = t * M(e-1, t-1, PCx, p) + PCx * M(e-1, t, PCx, p) + 1 / (2 * p) * M(e-1, t+1, PCx, p)
        end if
    end function M

    function SDipIntxyz(a, ijk1, ACoord, b, ijk2, BCoord, CCoord, order) result(SDipIntvalxyz)
        !        Evaluates multipole integral between two primitiv Gaussians
        !        Returns a float.
        !        Cxyz: origin of cartesian multipole moments
        !        order: order of multipole (1 for dipole...)
        !        a:    orbital exponent on Gaussian 'a' (e.g. alpha in the text)
        !        b:    orbital exponent on Gaussian 'b' (e.g. beta in the text)
        !        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
        !        for Gaussian 'a'
        !        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
        !        A:    list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
        !        B:    list containing origin of Gaussian 'b'
        !        Copied from above mentioned source
        integer, dimension(3), intent(in) :: ijk1, ijk2
        integer, intent(in) :: order
        real(r8), dimension(3), intent(in) :: ACoord, BCoord, CCoord
        real(r8), intent(in) :: a, b
        real(r8) :: p, S1, S2, S3
        real(r8), dimension(3) :: Se, PCx, SDipIntvalxyz
        integer :: i, t

        PCx = gaussianProductCenter(a,ACoord,b,BCoord)
        p = a + b
        Se = 0._r8
        SDipIntvalxyz = 1._r8
        do i = 1, 3
            select case(i)
            case(1)
                do t = 0, MIN(ijk1(1) + ijk2(1), order)
                    Se(1) = Se(1) + E(ijk1(1),ijk2(1),t,(ACoord(1) - BCoord(1)),a,b) * M(order, t,&
                    &(PCx(1) - CCoord(1)), p)
                end do
                !S1 = E(ijk1(1), ijk2(1), 0, (Axyz(1) - Bxyz(1)), a, b) !x
                S2 = E(ijk1(2), ijk2(2), 0, (ACoord(2) - BCoord(2)), a, b) !y
                S3 = E(ijk1(3), ijk2(3), 0, (ACoord(3) - BCoord(3)), a, b) !z
                SDipIntvalxyz(1) = Se(1) * S2 * S3 * pi/p
            case(2)
                do t = 0, MIN(ijk1(2) + ijk2(2), order)
                    Se(2) = Se(2) + E(ijk1(2),ijk2(2),t,(ACoord(2) - BCoord(2)),a,b) * M(order, t,&
                            &(PCx(2) - CCoord(2)), p)
                end do
                S1 = E(ijk1(1), ijk2(1), 0, (ACoord(1) - BCoord(1)), a, b) !x
                !S2 = E(ijk1(2), ijk2(2), 0, (ACoord(2) - BCoord(2)), a, b) !y
                S3 = E(ijk1(3), ijk2(3), 0, (ACoord(3) - BCoord(3)), a, b) !z
                SDipIntvalxyz(2) = Se(2) * S1 * S3 * pi/p
            case(3)
                do t = 0, MIN(ijk1(3) + ijk2(3), order)
                    Se(3) = Se(3) + E(ijk1(3),ijk2(3),t,(ACoord(3) - BCoord(3)),a,b) * M(order, t,&
                            &(PCx(3) - CCoord(3)), p)
                end do
                S1 = E(ijk1(1), ijk2(1), 0, (ACoord(1) - BCoord(1)), a, b) !x
                S2 = E(ijk1(2), ijk2(2), 0, (ACoord(2) - BCoord(2)), a, b) !y
                !S3 = E(ijk1(3), ijk2(3), 0, (ACoord(3) - BCoord(3)), a, b) !z
                SDipIntvalxyz(3) = Se(3) * S1 * S2 * pi/p
            end select
        end do

    end function SDipIntxyz

    function dipoleIntxyz(a, aCoords, b, bCoords, CCoords, order, ijk1, ijk2) result(DipIntxyz)
        !        Evaluates the multipole integral of nth order two contracted Gaussians
        !        Returns float.
        !        Arguments:
        !        Cxyz: origin of cartesian multipole moments
        !        order: order of multipole (1 for dipole...)
        !        a: contracted Gaussian 'a', BasisFunction object
        !        b: contracted Gaussian 'b', BasisFunction object
        !        Copied from above mentioned source
        class(cgfs), intent(in) :: a, b
        real(r8), dimension(3), intent(in) :: aCoords, bCoords, CCoords
        integer, dimension(3), optional, intent(in) :: ijk1, ijk2
        integer, intent(in) :: order
        integer(i8) :: ia, ib
        real(r8), dimension(3) :: DipIntxyz
        integer, dimension(3) :: ijka, ijkb

        if (PRESENT(ijk1)) then
            ijka = ijk1
            ijkb = ijk2
        else
            ijka = a%ijk
            ijkb = b%ijk
        end if
        DipIntxyz = 0._r8

        do ia = 1, SIZE(a%coefficients)
            do ib = 1, SIZE(b%coefficients)
                DipIntxyz = DipIntxyz + a%norms(ia) * b%norms(ib) * a%coefficients(ia) * b%coefficients(ib) *&
                        & SDipIntxyz(a%exponents(ia), ijka, aCoords, b%exponents(ib), ijkb, bCoords, CCoords, order)
            end do
        end do

        end function dipoleIntxyz

        function dipoleMatrices(atoms, geo, sizeAtoms, sizeSMatrix, order) result(dipoleMatrixxyz)

            integer, intent(in)                                         :: sizeAtoms, sizeSMatrix
            class(atom), dimension(sizeAtoms), intent(in)               :: atoms
            class(geom), intent(in)                                     :: geo
            real(r8), dimension(sizeSMatrix,sizeSMatrix,3)              :: dipoleMatrixxyz
            integer, intent(in)                                         :: order
            integer                                                     :: atomAIndex, atomBIndex, aoAIndex, aoBIndex
            integer                                                     :: aoAIMax, aoBIMax, a, b

            MPMatrix = 0._r8
            t = 0
            a = 1
            b = 1

            do atomAIndex = 1, sizeAtoms
                aoAIMax = SIZE(atoms(atomAIndex)%contractedGaussFunctions)
                do aoAIndex = 1, aoAIMax
                    do atomBIndex = 1, sizeAtoms
                        aoBIMax = SIZE(atoms(atomBIndex)%contractedGaussFunctions)
                        do aoBIndex = 1, aoBIMax
                            dipoleMatrixxyz(a, b, :) = dipoleIntxyz(atoms(atomAIndex)%contractedGaussFunctions(aoAIndex),&
                                    &geo%coordinates(atomAIndex,:),atoms(atomBIndex)%contractedGaussFunctions(aoBIndex),&
                                    &geo%coordinates(atomBIndex,:), geo%centerOfCharge, order)
                            b = b + 1
                        end do
                    end do
                    b = 1
                    a = a + 1
                end do
            end do
        end function dipoleMatrices



end module dipoleIntegrals