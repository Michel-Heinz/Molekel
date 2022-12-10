! Created by michel on 04/12/22.

module nuclearAttractionIntegralJoshuaGoings

    use kinds
    use global
    use PGF
    use atomparams
    use math
    use overlapJoshuaGoings

    contains

        recursive function R(t,u,v,n,p,PCx,PCy,PCz,RPC) result(val)
!            Returns the Coulomb auxiliary Hermite integrals
!            Returns a float.
!            Arguments:
!            t,u,v: order of Coulomb Hermite derivative in x,y,z
!            (see defs in Helgaker and Taylor)
!            n: order of Boys function
!            PCx,y,z: Cartesian vector distance between Gaussian
!            composite center P and nuclear center C
!            RPC: Distance between P and C

            implicit none
            integer, value :: t, u, v, n
            real(r8), value :: PCx, PCy, PCz, RPC, p
            real(r8) :: val, Tb
            Tb = p * RPC * RPC
            val = 0.0_r8
            if (t == u .and. t == v .and. t == 0) then
                val = val + ((-2*p)**n) * Boys_func(n,Tb)
            elseif (t == u .and. t == 0) then
                if (v > 1) then
                    val = val + (v-1)*R(t,u,v-2,n+1,p,PCx,PCy,PCz,RPC)
                end if
                val = val + PCz * R(t,u,v-1,n+1,p,PCx,PCy,PCz,RPC)
            elseif (t == 0) then
                if (u > 1) then
                    val = val + (u-1) * R(t,u-2,v,n+1,p,PCx,PCy,PCz,RPC)
                end if
                val = val + PCy * R(t,u-1,v,n+1,p,PCx,PCy,PCz,RPC)
            else
                if (t > 1) then
                    val = val + (t-1) * R(t-2,u,v,n+1,p,PCx,PCy,PCz,RPC)
                end if
                val = val + PCx * R(t-1,u,v,n+1,p,PCx,PCy,PCz,RPC)
            end if

        end function R

        function nuclearAttraction(a,lmn1,ACoord,b,lmn2,BCoord,CCoord) result(nucA)
!            Evaluates kinetic energy integral between two Gaussians
!            Returns a float.
!            a: orbital exponent on Gaussian 'a' (e.g. alpha in the text)
!            b: orbital exponent on Gaussian 'b' (e.g. beta in the text)
!            lmn1: int tuple containing orbital angular momentum (e.g. (1,0,0))
!            for Gaussian 'a'
!            lmn2: int tuple containing orbital angular momentum for Gaussian 'b'
!            ACoord: list containing origin of Gaussian 'a', e.g. [1.0, 2.0, 0.0]
!            BCoord: list containing origin of Gaussian 'b'
!            CCoord: list containing origin of nuclear center 'C'

            integer, dimension(3), intent(in) :: lmn1, lmn2
            real(r8), intent(in) :: a, b
            real(r8), dimension(3), intent(in) :: ACoord, BCoord, CCoord
            real(r8), dimension(3) :: gpc
            real(r8) :: nucA, p, RPC
            integer :: t, u, v, l1, l2, m1, m2, n1, n2

            l1 = lmn1(1)
            m1 = lmn1(2)
            n1 = lmn1(3)
            l2 = lmn2(1)
            m2 = lmn2(2)
            n2 = lmn2(3)

            p = a + b
            gpc = gaussianProductCenter(a,ACoord,b,BCoord)
            RPC = 0._r8
            RPC = NORM2(gpc - CCoord)

            nucA = 0.0_r8
            do t = 0, l1 + l2
                do u = 0, m1 + m2
                    do v = 0, n1 + n2
                        nucA = nucA + E(l1,l2,t,(ACoord(1)-BCoord(1)),a,b) * E(m1,m2,u,(ACoord(2)-BCoord(2)),a,b) *&
                                &E(n1,n2,v,(ACoord(3)-BCoord(3)),a,b) *&
                                &R(t,u,v,0,p,(gpc(1)-CCoord(1)),(gpc(2)-CCoord(2)),(gpc(3)-CCoord(3)),RPC)
                    end do
                end do
            end do
            nucA = nucA * 2 * pi / p
        end function nuclearAttraction

        function V(a, ACoords, b, BCoords, CCoords) result(nucA)
!            Evaluates overlap between two contracted Gaussians
!            Returns float.
!            Arguments:
!            a: contracted Gaussian 'a', BasisFunction object
!            b: contracted Gaussian 'b', BasisFunction object
!            C: center of nucleus
            implicit none
            class(cgfs), intent(in) :: a, b
            real(r8), dimension(3), intent(in) :: ACoords, BCoords, CCoords
            real(r8) :: nucA
            integer :: ia, ib

            nucA = 0._r8
            do ia = 1, SIZE(a%coefficients)
                do ib = 1, SIZE(b%coefficients)
                    nucA = nucA + a%norms(ia) * b%norms(ib) * a%coefficients(ia) * b%coefficients(ib) *&
                            &nuclearAttraction(a%exponents(ia), a%ijk, ACoords, b%exponents(ib), b%ijk, BCoords,&
                            &CCoords)
                end do
            end do
        end function V

        function nuclearAMatrix(atoms, geo, sizeAtoms, sizeMatrix) result(nucAMatrix)
            integer, intent(in) :: sizeAtoms, sizeMatrix
            type(atom), dimension(sizeAtoms), intent(in) :: atoms
            class(geom), intent(in) :: geo
            real(r8), dimension(sizeMatrix, sizeMatrix) :: nucAMatrix
            integer :: atomAIndex, atomBIndex, atomCIndex, aoAIndex, aoBIndex, aoAIMax, aoBIMax, a, b, i

            a = 1
            b = 1
            nucAMatrix = 0._r8

            do atomAIndex = 1, sizeAtoms
                aoAIMax = SIZE(atoms(atomAIndex)%contractedGaussFunctions)
                do aoAIndex = 1, aoAIMax
                    do atomBIndex = 1, sizeAtoms
                        aoBIMax = SIZE(atoms(atomBIndex)%contractedGaussFunctions)
                        do aoBIndex = 1, aoBIMax
                            do atomCIndex = 1, sizeAtoms
                                nucAMatrix(a, b) = nucAMatrix(a, b) - V(atoms(atomAIndex)%contractedGaussFunctions(aoAIndex), &
                                        &geo%coordinates(atomAIndex,:), &
                                        &atoms(atomBIndex)%contractedGaussFunctions(aoBIndex),&
                                        &geo%coordinates(atomBIndex,:), geo%coordinates(atomCIndex,:)) *&
                                        &geo%atomicNumber(atomCIndex)
                            end do
                            b = b + 1
                        end do
                    end do
                    b = 1
                    a = a + 1
                end do
            end do

        end function nuclearAMatrix

end module nuclearAttractionIntegralJoshuaGoings