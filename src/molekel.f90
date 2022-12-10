program molekel

    use readbasisset
    use readinputfile
    use init
    use global
    use kinds
    use PGF
    use RHF
    use elements
    use atomparams
    use integralsMain
    use geometry
    use parsing
    use outputParsing
    use initEval

    implicit none
    character(len=13) :: basisset
    character(len=2), dimension(:), allocatable :: atoms

    character(len=500) :: method, basis, filename
    real(r8), dimension(:,:), allocatable :: coordinates

    integer :: i, numberofarg, ierr, nLines, char, space
    character(len=MAXLEN), dimension(MAXLINES)       :: allLines
    class(atom), dimension(:), allocatable :: ohoh
    class(geom), allocatable :: geo
    real(r8), dimension(:,:), allocatable :: S_, T_, V_
    real(r8), dimension(:), allocatable :: int2
    real(r8) :: Energy

    real(r8) :: testVar, test, nucRep
    char = len('    __  _______  __    ________ __ ________  ')
    space = (MAXCHAR - CHAR) / 2

    write(*,*) repeat(' ', space)//'    __  _______  __    ________ __ ________  '
    write(*,*) repeat(' ', space)//'   /  |/  / __ \/ /   / ____/ //_// ____/ /  '
    write(*,*) repeat(' ', space)//'  / /|_/ / / / / /   / __/ / ,<  / __/ / /   '
    write(*,*) repeat(' ', space)//' / /  / / /_/ / /___/ /___/ /| |/ /___/ /___ '
    write(*,*) repeat(' ', space)//'/_/  /_/\____/_____/_____/_/ |_/_____/_____/ '
    write(*,*) repeat(' ', space)//' '
    call printCenter('Michel V. Heinz')


    call getenvvar(envvar, path)

    call energyEvaluation()




end program molekel