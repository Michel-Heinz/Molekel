module init

    contains

    subroutine getenvvar(name, path)

        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: path

        call getenv(name, path)
        path = trim(path)

        if(path == '') STOP 'Environment variable MOLEKEL not set!'

    end subroutine getenvvar

end module init