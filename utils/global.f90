module global

    use kinds

    !global parameters
    character(len=*), parameter :: envvar = 'MOLEKEL'
    real(r8), parameter :: angstoau = 1.889726134_r8
    real(r8), parameter :: pi = 4.0_r8 * DATAN(1.0_r8)
    real(r8), parameter :: auToDebye= 2.541765
    integer, parameter :: MAXLINES=1000
    integer, parameter :: MAXLEN=200
    integer, parameter :: MAXCHAR=79


    !global variables(r8)
    character(len=500) :: path

end module global