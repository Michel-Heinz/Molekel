module PGF

    use global
    use kinds
    use math

    type cgfs

        real(r8), dimension(:), allocatable :: exponents
        real(r8), dimension(:), allocatable :: coefficients
        real(r8), dimension(:), allocatable :: norms
        integer :: angularmomentum
        integer, dimension(3) :: ijk

        contains
        procedure :: set_cgfs
        procedure :: get_exponents
        procedure :: get_coefficients
        procedure :: get_angularmomentum
        procedure :: get_ijk
        procedure :: get_norms

    end type cgfs

    contains

    subroutine set_cgfs(this, angularmomentum, ijk, exponents, coefficients)
        class(cgfs), intent(inout) :: this
        integer, intent(in) :: angularmomentum
        integer, dimension(3), intent(in) :: ijk
        real(r8), dimension(:), intent(in) :: exponents
        real(r8), dimension(:), intent(in) :: coefficients


        integer :: i, numberofcgfs

        numberofcgfs = SIZE(exponents)
        allocate(this%exponents(numberofcgfs))
        allocate(this%coefficients(numberofcgfs))
        allocate(this%norms(numberofcgfs))
        this%angularmomentum = angularmomentum
        this%ijk = ijk

        do i = 1, numberofcgfs
            this%exponents(i) = exponents(i)
            this%coefficients(i) = coefficients(i)
        end do

        call normalize_PGF(this%norms, this%exponents, this%ijk)

        normCGTO = 0.0_r8
        call normalizeCGTO(this%exponents, this%coefficients, this%norms, this%ijk)

    end subroutine set_cgfs

    function get_exponents(this)
        implicit none
        class(cgfs), intent(in) :: this
        real(r8), dimension(:), allocatable :: get_exponents
        get_exponents = this%exponents

    end function get_exponents

    function get_coefficients(this)
        implicit none
        class(cgfs), intent(in) :: this
        real(r8), dimension(:), allocatable :: get_coefficients
        get_coefficients = this%coefficients

    end function get_coefficients

    function get_angularmomentum(this)
        implicit none
        class(cgfs), intent(in) :: this
        integer :: get_angularmomentum
        get_angularmomentum = this%angularmomentum

    end function get_angularmomentum

    function get_ijk(this)
        implicit none
        class(cgfs), intent(in) :: this
        integer, dimension(3) :: get_ijk
        get_ijk = this%ijk

    end function get_ijk

    function get_norms(this)
        implicit none
        class(cgfs), intent(in) :: this
        integer, dimension(:), allocatable :: get_norms
        get_norms = this%norms

    end function get_norms

    subroutine create_ijk(angularmomentum, ijk_a)
        implicit none
        integer, intent(in) :: angularmomentum
        integer, dimension(:,:), allocatable, intent(out) :: ijk_a
        integer :: a, b, c, d

        select case(angularmomentum)
            case(0)
            allocate(ijk_a(1,3))
            ijk_a = 0
            case(1)
            allocate(ijk_a(3,3))
            ijk_a = 0
            do a = 1, 3
                ijk_a(a,a) = 1
            end do
            case(2)
            allocate(ijk_a(6,3))
            ijk_a = 0
            do a = 1, 3
                ijk_a(a,a) = 2
            end do
            ijk_a(4:6,:) = 1
            do a = 1, 3
                ijk_a(a+3,4-a) = 0
            end do
            case(3)
            allocate(ijk_a(10,3))
            ijk_a = 0
            do a = 1, 3
                ijk_a(a,a) = 3
            end do
            ijk_a(10,:) = 1
            ijk_a(4,1) = 1
            ijk_a(4,2) = 2
            ijk_a(5,1) = 2
            ijk_a(5,2) = 1
            ijk_a(6,1) = 2
            ijk_a(6,3) = 1
            ijk_a(7,1) = 1
            ijk_a(7,3) = 2
            ijk_a(8,2) = 1
            ijk_a(8,3) = 2
            ijk_a(9,2) = 2
            ijk_a(9,3) = 1
!            ijk_a = 0
!            ijk_a(4,:) = 1
!            ijk_a(4:10,:) = 2
!            do a = 1, 3
!                ijk_a(a,a) = 3
!            end do
!            do a = 5, 7
!                ijk_a(a,a-4) = 0
!                if(a == 7) then
!                    ijk_a(a,2) = 1
!                else
!                    ijk_a(a,3) = 1
!                end if
!            end do
!            do a = 8, 10
!                ijk_a(a,a-7) = 0
!                if(a == 8) then
!                    ijk_a(a,2) = 1
!                else
!                    ijk_a(a,1) = 1
!                end if
!            end do
        end select

    end subroutine create_ijk

    subroutine normalize_PGF(norms, exps, ijk)
        implicit none
        real(r8), dimension(:), intent(inout) :: norms
        real(r8), dimension(:), intent(in) :: exps
        integer, dimension(3), intent(in) :: ijk
        integer :: i, j, k, L, o, fi, fj, fk
        real(r8), dimension(3) :: coord


        i = ijk(1)
        j = ijk(2)
        k = ijk(3)
        L = i + j + k

        fi = fact2(2 * i - 1)
        fj = fact2(2 * j - 1)
        fk = fact2(2 * k - 1)

        norms = DSQRT((2 ** (2 * L + 1.5_r8) * exps ** (L + 1.5_r8)) / (fi * fj * fk * pi ** 1.5_r8))

        end subroutine normalize_PGF

    subroutine normalizeCGTO(exps, coeffs, norms, ijk)
        real(r8), dimension(:), intent(in) :: exps, norms
        real(r8), dimension(:), intent(inout) :: coeffs
        integer, dimension(3), intent(in) :: ijk
        real(r8) :: prefactor, N
        integer :: ia, ib

        i = ijk(1)
        j = ijk(2)
        k = ijk(3)
        L = i + j + k

        fi = fact2(2 * i - 1)
        fj = fact2(2 * j - 1)
        fk = fact2(2 * k - 1)

        prefactor = (pi ** (1.5_r8) * fi * fj * fk) / (2 ** (L))

        N = 0.0_r8

        do ia = 1, SIZE(exps)
            do ib = 1, SIZE(exps)
                !write(*,*) ((exps(ia) * exps(ib)) ** (1))
                N = N + ((norms(ia) * norms(ib) * coeffs(ia) * coeffs(ib)) / ((exps(ia) + exps(ib)) ** (L + 1.5_r8)))
            end do
        end do

        N = N * prefactor
        coeffs = coeffs * N ** (-0.5_r8)

        end subroutine normalizeCGTO

end module PGF