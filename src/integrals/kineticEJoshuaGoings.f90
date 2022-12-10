! Created by michel on 02/12/22.
! Source https://joshuagoings.com/2017/04/28/integrals/ 21.5.22 20:46

module kineticEJoshuaGoings

    use kinds
    use overlapJoshuaGoings
    use global
    use PGF
    use atomparams

    contains

    function kinetic(a, ijk1, Axyz, b, ijk2, Bxyz) result(eKin)
!        Evaluates kinetic energy integral between two Gaussians
!        Returns a float.
!        a:        orbital exponent on Gaussian 'a' (e.g. alpha in the text)
!        b:        orbital exponent on Gaussian 'b' (e.g. beta in the text)
!        lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
!        for Gaussian 'a'
!        lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
!        A:        list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
!        B:        list containing origin of Gaussian 'b'
        real(r8), intent(in) :: a, b
        real(r8), intent(in), dimension(3) :: Axyz, Bxyz
        integer, intent(in), dimension(3) :: ijk1, ijk2
        real(r8) :: eKin, term0, term1, term2
        integer, dimension(3) :: temp1, temp2, temp3

        eKin = 0
        term0 = b*(2*(ijk2(1)+ijk2(2)+ijk2(3))+3)*overlap(a,ijk1,Axyz,b,ijk2,Bxyz)

        temp1 = ijk2
        temp1(1) = temp1(1) + 2
        temp2 = ijk2
        temp2(2) = temp2(2) + 2
        temp3 = ijk2
        temp3(3) = temp3(3) + 2

        term1 = -2*(b**2)*(overlap(a,ijk1,Axyz,b,temp1,Bxyz) +&
                &overlap(a,ijk1,Axyz,b,temp2,Bxyz) +&
                &overlap(a,ijk1,Axyz,b,temp3,Bxyz))

        temp1 = ijk2
        temp1(1) = temp1(1) - 2
        temp2 = ijk2
        temp2(2) = temp2(2) - 2
        temp3 = ijk2
        temp3(3) = temp3(3) - 2

        term2 = -0.5*(ijk2(1)*(ijk2(1)-1)*overlap(a,ijk1,Axyz,b,temp1,Bxyz) +&
                &ijk2(2)*(ijk2(2)-1)*overlap(a,ijk1,Axyz,b,temp2,Bxyz) +&
                &ijk2(3)*(ijk2(3)-1)*overlap(a,ijk1,Axyz,b,temp3,Bxyz))
        eKin = term0+term1+term2
    end function kinetic

    function T(a, aCoords, b, bCoords) result(eKin)
!        Evaluates kinetic energy between two contracted Gaussians
!        Returns float.
!        Arguments:
!        a: contracted Gaussian 'a', BasisFunction object
!        b: contracted Gaussian 'b', BasisFunction object
        class(cgfs), intent(in) :: a, b
        real(r8), dimension(3), intent(in) :: aCoords, bCoords
        real(r8) :: eKin
        integer :: ia, ib

        eKin = 0._r8
        do ia = 1, SIZE(a%coefficients)
            do ib = 1, SIZE(b%coefficients)
                eKin = eKin + a%norms(ia) * b%norms(ib) * a%coefficients(ia) * b%coefficients(ib) *&
                        &kinetic(a%exponents(ia), a%ijk, aCoords, b%exponents(ib), b%ijk, bCoords)
            end do
        end do

    end function T

    function kineticEMatrix(atoms, geo, sizeAtoms, sizeMatrix) result(kinEMatrix)
        integer, intent(in) :: sizeAtoms, sizeMatrix
        type(atom), dimension(sizeAtoms), intent(in) :: atoms
        class(geom), intent(in) :: geo
        real(r8), dimension(sizeMatrix, sizeMatrix) :: kinEMatrix
        integer :: atomAIndex, atomBIndex, aoAIndex, aoBIndex, aoAIMax, aoBIMax, a, b, i

        kinEMatrix = 0._r8
        a = 1
        b = 1

        do atomAIndex = 1, sizeAtoms
            aoAIMax = SIZE(atoms(atomAIndex)%contractedGaussFunctions)
            do aoAIndex = 1, aoAIMax
                do atomBIndex = 1, sizeAtoms
                    aoBIMax = SIZE(atoms(atomBIndex)%contractedGaussFunctions)
                    do aoBIndex = 1, aoBIMax
                        kinEMatrix(a, b) = T(atoms(atomAIndex)%contractedGaussFunctions(aoAIndex), &
                                &geo%coordinates(atomAIndex,:), atoms(atomBIndex)%contractedGaussFunctions(aoBIndex),&
                                &geo%coordinates(atomBIndex,:))
                        b = b + 1
                    end do
                end do
                b = 1
                a = a + 1
            end do
        end do

    end function kineticEMatrix


end module kineticEJoshuaGoings