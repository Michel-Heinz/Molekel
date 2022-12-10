! Copyright (C) 2012, 2018 luechow
!
! SPDX-License-Identifier: GPL-3.0-or-later

module errorAssert
    use kinds, only: r8
    implicit none

    ! error handling and assert
    ! program termination handled only by error(): adapt if necessary
    ! in time critical routines use:
    ! if (asserts) call assert(condition)
    ! to allow compiler to remove asserts with assert = .false.

    logical, parameter :: asserts = .true.
    logical, parameter :: debug = .true.

contains

    subroutine abortp(s)
        ! general abnormal thread termination
        character(*), intent(in) :: s

        write(*,*) ' *** ABNORMAL TERMINATION *** '
        write(*,*) s

        error stop
    end subroutine abortp

    subroutine error(string)
        character(len=*), optional, intent(in) :: string

        if (present(string)) then
            call abortp(string)
        else
            call abortp("ERROR")
        end if
    end subroutine error


    subroutine assert(logicalExp, string)
        logical, intent(in)                    :: logicalExp
        character(len=*), optional, intent(in) :: string
        character(len=120)                      :: errstring

        if (.not. logicalExp) then
            if (present(string)) then
                errstring = 'assertion violated: '//string
            else
                errstring = 'assertion violated!'
            endif
            call error(errstring)
        end if
    end subroutine assert


    subroutine assertEqualAbsolute(t1, t2, tol, string)
        real(r8), intent(in)                     :: t1, t2, tol
        character(len=*), optional, intent(in) :: string
        character(len=120)                      :: errstring

        if (abs(t1-t2) > tol) then
            if (present(string)) then
                errstring = 'assertion violated: '//string
            else
                errstring = 'assertion violated!'
            endif
            call error(errstring)
        end if
    end subroutine assertEqualAbsolute


    subroutine assertEqualRelative(t1, t2, tol, string)
        real(r8), intent(in)                     :: t1, t2, tol
        character(len=*), optional, intent(in) :: string
        character(len=120)                      :: errstring

        if (abs((t1-t2)/t1) > tol) then
            if (present(string)) then
                errstring = 'assertion violated: '//string
            else
                errstring = 'assertion violated!'
            endif
            call error(errstring)
        end if
    end subroutine assertEqualRelative


end module errorAssert