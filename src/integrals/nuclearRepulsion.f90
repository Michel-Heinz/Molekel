! Created by michel on 02/12/22.

module nuclearRepulsion

    use kinds
    use global
    use atomparams

    contains

    function vNuc(atoms, geo) result(vNuc_)

        type(atom), intent(in), dimension(:) :: atoms
        class(geom), intent(in) :: geo
        real(r8) :: vNuc_, dist
        integer :: atomIndex1, atomIndex2, numAtoms, n1, n2, luaf

        vNuc_ = 0._r8
        numAtoms = SIZE(atoms)
        luaf = 1

        do atomIndex1=1, numAtoms
            n1 = geo%atomicNumber(atomIndex1)
            do atomIndex2=atomIndex1 + 1, numAtoms
                n2 = geo%atomicNumber(atomIndex2)
                dist = NORM2(geo%coordinates(atomIndex1,:) - geo%coordinates(atomIndex2,:))
                vNuc_ = vNuc_ + (n1*n2)/dist
            end do
        end do

    end function vnuc
end module nuclearRepulsion