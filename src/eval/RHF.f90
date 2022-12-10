! Created by michel on 04/12/22.

module RHF

    use kinds
    use global
    use math
    use integralsMain
    use atomparams
    use outputParsing
    use dipolMoment

    contains

    function getHCore(T, V) result(HCore)
        implicit none
        real(r8), dimension(:,:), intent(in) :: T, V
        real(r8), dimension(SIZE(T(1,:)),SIZE(T(1,:))) :: HCore
        integer :: mu, nu

        HCore = T + V
    end function getHCore

    function getDensityMatrix(C, N) result(D)
        implicit none
        real(r8), dimension(:,:), intent(in) :: C
        integer, intent(in) :: N
        real(r8), dimension(SIZE(C(1,:)),SIZE(C(1,:))) :: D
        integer :: mu, nu, i, numberOfOrbs, NOcc

        numberOfOrbs = SIZE(C(1,:))
        NOcc = N / 2
        D = 0._r8
        do mu = 1, numberOfOrbs
            do nu = 1, numberOfOrbs
                do i = 1, NOcc
                    D(mu, nu) = D(mu, nu) + 2 * C(mu, i) * C(nu, i)
                end do
            end do
        end do

    end function getDensityMatrix

    function getTwoeContrib(int2e, D, numberOfOrbs) result(G)
        real(r8), dimension(:,:), intent(in) :: D
        integer, intent(in) :: numberOfOrbs
        real(r8), dimension(:), intent(in) :: int2e
        integer :: lambda, sigma, mu, nu
        real(r8), dimension(numberOfOrbs,numberOfOrbs) :: G

        G = 0._r8
        do mu = 1, numberOfOrbs
            do nu = 1, numberOfOrbs
                do lambda = 1, numberOfOrbs
                    do sigma = 1, numberOfOrbs
                        G(mu, nu) = G(mu, nu) + D(lambda, sigma) * (int2e(eriIndex(mu-1, nu-1, lambda-1, sigma-1)+1)&
                                &- 0.5_r8 * int2e(eriIndex(mu-1, sigma-1, lambda-1, nu-1)+1))
                    end do
                end do
            end do
        end do

    end function getTwoeContrib

    function getFockMatrix(HCore, D, int2e, N) result(F)
        implicit none
        real(r8), dimension(:,:), intent(in) :: HCore, D
        real(r8), dimension(:), intent(in) :: int2e
        integer, intent(in) :: N
        real(r8), dimension(SIZE(HCore(1,:)),SIZE(HCore(1,:))) :: F
        integer :: NOcc, mu, nu, lambda, sigma, numberOfOrbs, i

        numberOfOrbs = SIZE(HCore(1,:))
        NOcc = N / 2
        F = HCore

        do mu = 1, numberOfOrbs
           do nu = 1, numberOfOrbs
               do lambda = 1, numberOfOrbs
                   do sigma = 1, numberOfOrbs
                       F(mu, nu) = F(mu, nu) + D(lambda, sigma) * (int2e(eriIndex(mu-1, nu-1, lambda-1, sigma-1)+1)&
                               &- 0.5_r8 * int2e(eriIndex(mu-1, sigma-1, lambda-1, nu-1)+1))
                   end do
               end do
           end do
        end do

    end function getFockMatrix

    function getRHFEnergy(F, D, HCore, Vnuc, numberOfOrbs) result(ERHF)
        implicit none
        real(r8), dimension(:,:), intent(in) :: F, D, HCore
        real(r8), intent(in) :: Vnuc
        integer, intent(in) :: numberOfOrbs
        real(r8) :: ERHF
        integer :: mu ,nu

        ERHF = 0._r8

        do mu = 1, numberOfOrbs
            do nu = 1, numberOfOrbs
                ERHF = ERHF + (D(mu ,nu) * (HCore(mu, nu) + F(mu ,nu)))
            end do
        end do
        ERHF = 0.5_r8 * ERHF + Vnuc
    end function getRHFEnergy

    !M. Wolfsberg and L. Helmholtz, J. Chem. Phys. 20, 837 (1952).
    function GWHGuess(HCore, S, numberOfOrbs, c_) result(guess)
        implicit none
        real(r8), dimension(:,:), intent(in) :: HCore, S
        real(r8), intent(in), optional :: c_
        integer, intent(in) :: numberOfOrbs
        real(r8), dimension(numberOfOrbs, numberOfOrbs) :: guess
        real(r8) :: c
        integer :: i, j

        if (.not. PRESENT(c_)) then
            c = 2.3_r8
        else
            c = c_
        end if


        do i = 1, numberOfOrbs
            do j = 1, numberOfOrbs
                if (i == j) then
                    guess(i, j) = HCore(i ,j)
                else
                    guess(i, j) = c * S(i ,j) * (HCore(i, i) + HCore(j, j)) / 2
                end if
            end do
        end do


    end function GWHGuess

    function getKinEnergy(D, T, numberOfOrbs) result(eKin)
        integer, intent(in) :: numberOfOrbs
        real(r8), dimension(numberOfOrbs,numberOfOrbs), intent(in) :: D, T
        integer :: mu, nu
        real(r8) :: eKin

        eKin = 0._r8

        do nu = 1, numberOfOrbs
            do mu = 1, numberOfOrbs
                eKin = eKin + D(mu, nu) * T(mu, nu)
            end do
        end do

    end function getKinEnergy

    function getTwoEEnergy(int2e , D, numberOfOrbs) result(twoEE)
        real(r8), dimension(:,:), intent(in) :: D
        integer, intent(in) :: numberOfOrbs
        real(r8), dimension(:), intent(in) :: int2e
        integer :: mu, nu
        real(r8), dimension(numberOfOrbs,numberOfOrbs) :: G
        real(r8) :: twoEE

        G = getTwoeContrib(int2e , D, numberOfOrbs)
        twoEE = 0._r8
        do nu = 1, numberOfOrbs
            do mu = 1, numberOfOrbs
                twoEE = twoEE + D(mu, nu) * G(mu, nu) * 1 / 2
            end do
        end do

    end function getTwoEEnergy

    function getPotEnergy(D, V, numberOfOrbs) result(ePot)
        integer, intent(in) :: numberOfOrbs
        real(r8), dimension(numberOfOrbs,numberOfOrbs), intent(in) :: D, V
        integer :: mu, nu
        real(r8) :: ePot

        ePot = 0._r8
        do nu = 1, numberOfOrbs
            do mu = 1, numberOfOrbs
                ePot = ePot + D(mu, nu) * V(mu, nu)
            end do
        end do

    end function getPotEnergy


    function doRHF(atoms, geo, charge, maxIter, precision, guess, verbose) result(ERHF)
        implicit none
        class(atom), dimension(:), intent(in) :: atoms
        class(geom), intent(in) :: geo
        real(r8), dimension(:,:), allocatable :: S, T, V, HCore, D, SLapack, F, guessMat, C
        real(r8), dimension(:), allocatable :: int2, epsilon, work
        integer, intent(in) :: charge, maxIter, guess, verbose
        real(r8), intent(in) :: precision
        real(r8), dimension(4) :: DipMom
        character(len=9), dimension(:), allocatable :: rowLabels
        character(len=10), dimension(:), allocatable :: columnLabels
        real(r8) :: ERHF, Vnuc, oldERHF, mux, suma
        integer :: i, j, N, numberOfOrbs, lwork, liwork, info, iter, A
        integer, dimension(:), allocatable :: iwork
        logical :: converged

        N = 0
        do i = 1, SIZE(atoms)
            N = N + geo%atomicNumber(i)
        end do
        N = N - charge

        call getIntegrals(atoms, geo, S, T, V, int2, Vnuc)
        numberOfOrbs = SIZE(S(1,:))
        allocate(D(numberOfOrbs, numberOfOrbs), SLapack(numberOfOrbs, numberOfOrbs), epsilon(numberOfOrbs),&
                &guessMat(numberOfOrbs, numberOfOrbs), rowLabels(numberOfOrbs), columnLabels(numberOfOrbs), &
                &C(numberOfOrbs, numberOfOrbs))
        if (verbose >= 3) then
            rowlabels = createOrbLabels(atoms)
            columnLabels = createOrbLabels(atoms)
            call printTitle('Overlap Matrix')
            call printMatrix(S, columnLabels, rowlabels)
        end if
        C = 0._r8
        SLapack = S
        lwork = 1 + 6 * numberOfOrbs + 2 * numberOfOrbs**2
        liwork = 3 + 5 * numberOfOrbs
        allocate(work(lwork), iwork(liwork))
        work = 0._r8
        iwork = 0

        Hcore = getHCore(T, V)
        guessMat = GWHGuess(HCore, S, numberOfOrbs)

        call DSYGV(1, 'V', 'U', numberOfOrbs, guessMat, numberOfOrbs, SLapack, numberOfOrbs, epsilon, work, lwork, info)

        D = getDensityMatrix(guessMat, N)
        F = getFockMatrix(HCore, D, int2, N)
        oldERHF = 0._r8

        converged = .false.

        if (verbose >= 3) call printTitle('SCF Started')
        if (verbose >= 3) then
            write(*,'(A9, A20, A30)') 'Iteration', 'Energy (Eh)', 'Energy Difference (Eh)'
        end if
        do iter = 1, maxIter
            SLapack = S

            call DSYGV(1, 'V', 'U', numberOfOrbs, F, numberOfOrbs, SLapack, numberOfOrbs, epsilon, work, lwork, info)

            C = F
            D = getDensityMatrix(F, N)
            F = getFockMatrix(HCore, D, int2, N)

            ERHF = getRHFEnergy(F, D, HCore, Vnuc, numberOfOrbs)
            if (verbose >= 3) then
                write(*,'(I9, F20.10, F30.10)') iter, ERHF, ERHF-oldERHf
            end if
            if (ABS(oldERHF - ERHF) < precision) then
                converged = .true.
                exit
            end if
            oldERHF = ERHF
        end do

        DipMom = getDipolMoments(atoms, geo, numberOfOrbs, D)

        if (converged .and. verbose >= 1) then
            call printTitle('SCF converged')
            call printInt('Iterations:', iter)
            call printFloat('Final RHF Energy:', ERHF, 'Eh')
            if (verbose >= 2) then
                write(*,*) ' '
                call printFloat('Nuclear Repulsion Energy:', Vnuc, 'Eh')
                call printFloat('1 Electron Energy:', getKinEnergy(D, T, numberOfOrbs) + getPotEnergy(D, V,&
                        &numberOfOrbs), 'Eh')
                call printFloat('2 Electron Energy:', getTwoEEnergy(int2 , D, numberOfOrbs), 'Eh')
                write(*,*) ' '
                call printFloat('Kinetic Energy:', getKinEnergy(D, T, numberOfOrbs), 'Eh')
                call printFloat('Potental Energy:', getPotEnergy(D, V, numberOfOrbs) + &
                        &Vnuc + getTwoEEnergy(int2 , D, numberOfOrbs), 'Eh')
                call printFloat('Virial Ratio:', -(getPotEnergy(D, V, numberOfOrbs) + &
                        &Vnuc + getTwoEEnergy(int2 , D, numberOfOrbs)) / getKinEnergy(D, T, numberOfOrbs))
                if (verbose >= 3) then
                    call printTitle('Orbital Energies and Occupation')
                    write(*,'(A18, A15, A15)') 'Molecular Orbital', 'Occupation', 'Energy (Eh)'
                    do i = 1, numberOfOrbs
                        if (i * 2 <= N) write(*,'(I18, F15.4, F15.4)') i, 2._r8, epsilon(i)
                        if (i * 2 > N) write(*,'(I18, F15.4, F15.4)') i, 0._r8, epsilon(i)
                    end do
                    call printTitle('Orbital Coefficients')
                    call printMatrix(C, rowlabels=rowlabels)
                end if
                call printTitle('Dipolemoment')
                call printFloat('Dipolemoment (x):', DipMom(1), 'D')
                call printFloat('Dipolemoment (y):', DipMom(2), 'D')
                call printFloat('Dipolemoment (z):', DipMom(3), 'D')
                call printFloat('Dipolemoment (norm):', DipMom(4), 'D')
            end if
        end if



        do i = 1, numberOfOrbs
            rowLabels(i) = 'test'
            columnLabels(i) = 'test2'
        end do
        rowlabels = createOrbLabels(atoms)
        columnLabels = createOrbLabels(atoms)

        !call printMatrix(S, columnLabels=columnLabels, rowlabels=rowLabels)
        !call printMatrix(C, rowlabels=rowLabels)
        !print*, epsilon
        !call printMatrix
    end function doRHF
end module RHF