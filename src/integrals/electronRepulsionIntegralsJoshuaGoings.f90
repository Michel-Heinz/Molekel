! Created by michel on 04/12/22.

module electronRepulsionIntegralsJoshuaGoings

    use nuclearAttractionIntegralJoshuaGoings
    use kinds
    use global
    use PGF
    use atomparams
    use math
    use overlapJoshuaGoings

    contains

    function electronRepulsion(a,lmn1,ACoord,b,lmn2,BCoord,c,lmn3,CCoord,d,lmn4,DCoord) result(eri)
!    Evaluates kinetic energy integral between two Gaussians
!    Returns a float.
!    a,b,c,d: orbital exponent on Gaussian 'a','b','c','d'
!    lmn1,lmn2
!    lmn3,lmn4: int tuple containing orbital angular momentum
!    for Gaussian 'a','b','c','d', respectively
!    A,B,C,D: list containing origin of Gaussian 'a','b','c','d'

        implicit none
        real(r8), intent(in) :: a, b, c, d
        real(r8), dimension(3), intent(in) :: ACoord, BCoord, CCoord, DCoord
        integer, dimension(3), intent(in) :: lmn1, lmn2, lmn3, lmn4
        integer :: t, u, v, tau, nu, phi
        real(r8) :: p, q, alpha, RPQ, eri
        real(r8), dimension(3) :: pgpc, qgpc

        p = 0._r8
        q = 0._r8
        alpha = 0._r8
        pgpc = 0._r8
        qgpc = 0._r8
        RPQ = 0._r8

        p = a + b
        q = c + d
        alpha = (p * q) / (p + q)
        pgpc = gaussianProductCenter(a, ACoord, b, BCoord)
        qgpc = gaussianProductCenter(c, CCoord, d, DCoord)
        RPQ = NORM2(pgpc - qgpc)
        eri = 0._r8

        do t = 0, lmn1(1) + lmn2(1)
            do u = 0, lmn1(2) + lmn2(2)
                do v = 0, lmn1(3) + lmn2(3)
                    do tau = 0, lmn3(1) + lmn4(1)
                        do nu = 0, lmn3(2) + lmn4(2)
                            do phi = 0, lmn3(3) + lmn4(3)
                                eri = eri + E(lmn1(1),lmn2(1),t,  ACoord(1)-BCoord(1),a,b) *&
                                           &E(lmn1(2),lmn2(2),u,  ACoord(2)-BCoord(2),a,b) *&
                                           &E(lmn1(3),lmn2(3),v,  ACoord(3)-BCoord(3),a,b) *&
                                           &E(lmn3(1),lmn4(1),tau,CCoord(1)-DCoord(1),c,d) *&
                                           &E(lmn3(2),lmn4(2),nu, CCoord(2)-DCoord(2),c,d) *&
                                           &E(lmn3(3),lmn4(3),phi,CCoord(3)-DCoord(3),c,d) *&
                                           &(-1)**(tau + nu + phi) *&
                                           &R(t+tau,u+nu,v+phi,0,alpha,(pgpc(1)-qgpc(1)),(pgpc(2)-qgpc(2)),&
                                           &(pgpc(3)-qgpc(3)),RPQ)
                            end do
                        end do
                    end do
                end do
            end do
        end do
        eri = eri * 2 * (pi**(2.5_r8)) / (p * q * SQRT(p+q))
        end function electronRepulsion

        function eri(a, ACoord, b, BCoord, c, CCoord, d, DCoord) result(eri_)
!        Evaluates the electon repulsion intergral
!        Returns float.
!        Arguments:
!        a: contracted Gaussian 'a', BasisFunction object
!        b: contracted Gaussian 'b', BasisFunction object
!        c: contracted Gaussian 'b', BasisFunction object
!        d: contracted Gaussian 'b', BasisFunction object

            class(cgfs), intent(in) :: a, b, c, d
            real(r8), dimension(3), intent(in) :: ACoord, BCoord, CCoord, DCoord
            real(r8) :: eri_
            integer :: ia, ib, ic, id, ll

            eri_ = 0._r8
            ll = 1

            do ia = 1, SIZE(a%coefficients)
                do ib = 1, SIZE(b%coefficients)
                    do ic = 1, SIZE(c%coefficients)
                        do id = 1, SIZE(d%coefficients)
                            eri_ = eri_ + a%norms(ia) * b%norms(ib) * c%norms(ic) * d%norms(id) *&
                                    &a%coefficients(ia) * b%coefficients(ib) * c%coefficients(ic) * d%coefficients(id)&
                                    &* electronRepulsion(a%exponents(ia),a%ijk,ACoord,b%exponents(ib),b%ijk,BCoord,&
                                    &c%exponents(ic),c%ijk,CCoord,d%exponents(id),d%ijk,DCoord)
                        end do
                    end do
                end do
            end do
        end function eri

        function numberOfUniqueInts(nOrbs) result(n)
            implicit none
            integer, intent(in) :: nOrbs
            integer :: n
            n = 0
            n = (nOrbs * (nOrbs + 1) * (nOrbs**2 + nOrbs + 2))/8
            end function numberOfUniqueInts

        function eriArray(atoms, geo, sizeAtoms, numberOfOrbs) result(eriList)
            implicit none
            integer, intent(in) :: sizeAtoms, numberOfOrbs
            type(atom), dimension(sizeAtoms), intent(in) :: atoms
            class(geom), intent(in) :: geo
            real(r8), dimension(numberOfOrbs, 3) :: Coords
            real(r8), dimension(:), allocatable :: eriList
            integer :: atomAIndex, atomBIndex, atomCIndex, atomDIndex, aoAIndex, aoBIndex, aoCIndex, aoDIndex
            integer :: aoAIMax, aoBIMax, aoCIMax, aoDIMax, i, j, k, l, x, y, nUniqueInts, orb, idx, ij, kl
            integer, dimension(numberOfOrbs) :: atomNumber, orbIDX

            nUniqueInts = numberOfUniqueInts(numberOfOrbs)
            ALLOCATE(eriList(nUniqueInts))
            eriList = 0._r8
            atomNumber = 0
            orbIDX = 0
            Coords = 0._r8

            orb = 1
            do x = 1, sizeAtoms
                do y = 1, SIZE(atoms(x)%contractedGaussFunctions)
                    atomNumber(orb) = x
                    orbIDX(orb) = y
                    Coords(orb,:) = geo%coordinates(x,:)
                    orb = orb + 1
                end do
            end do
            do i = 1, numberOfOrbs
                do j = i, numberOfOrbs
                    ij = (i * (i + 1) / 2 + j)
                    do k = 1, numberOfOrbs
                        do l = k, numberOfOrbs
                            kl = (k * (k + 1) / 2 + l)
                            if (ij >= kl) then
                                idx = eriIndex(i-1, j-1, k-1, l-1)
                                eriList(idx+1) = eri(atoms(atomNumber(i))%contractedGaussFunctions(orbIDX(i)), Coords(i&
                                        &,:), atoms(atomNumber(j))%contractedGaussFunctions(orbIDX(j)), Coords(j,:), &
                                        &atoms(atomNumber(k))%contractedGaussFunctions(orbIDX(k)), Coords(k,:), &
                                        &atoms(atomNumber(l))%contractedGaussFunctions(orbIDX(l)), Coords(l,:))
                            end if
                        end do
                    end do
                end do
            end do
        end function eriArray

end module electronRepulsionIntegralsJoshuaGoings