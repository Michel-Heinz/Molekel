module machine

    contains

    subroutine getarguments(n,s,ierr)
        !copied from Amolqc code ;)

        integer, intent(in)              :: n    ! get n-th command line argument
        character(len=*), intent(out)    :: s    ! return argument
        integer, intent(out)             :: ierr ! error code: 1=no argument

        ierr = 0
        if (iargc() < n) then
            ierr = 1
        else
            call getarg(n,s)
        endif
    end subroutine getarguments


end module machine