module elements
    use kinds
!!!! __________________________SOURCE: https://github.com/smillerc/fortran_periodic_table
    type :: element_t
        !< A simple type to represent an element from the periodic table
        character(len=:), allocatable :: name
        character(len=:), allocatable :: symbol
        real(r8) :: mass = 0.0_r8 !< Mass of the element in amu
        integer :: n_protons = 0
    end type

    contains
    subroutine getElementInfo(atomname, atomicnumber, mass)
        implicit none
        character(len=2), intent(in) :: atomname
        real(r8), intent(out) :: mass
        class(element_t), allocatable :: elem
        integer, intent(out) :: atomicnumber

        select case(atomname)
        case("hydrogen", "protium", "hydrogen-1", "H")
            elem = element_t(name="hydrogen", symbol="H", mass=1.00794_r8, n_protons=1)
        case("deuterium", "hydrogen-2", "2H", "D")
            elem = element_t(name="deuterium", symbol="D", mass=2.014101778_r8, n_protons=1)
        case("tritium", "hydrogen-3", "3H", "T")
            elem = element_t(name="tritium", symbol="T", mass=3.0160492675_r8, n_protons=1)
        case("helium", "He")
            elem = element_t(name="helium", symbol="He", mass=4.002602_r8, n_protons=2)
        case("helium-3", "3He", "He-3")
            elem = element_t(name="helium", symbol="He", mass=3.0160293097_r8, n_protons=2)
        case("lithium", "Li")
            elem = element_t(name="lithium", symbol="Li", mass=6.941_r8, n_protons=3)
        case("beryllium", "Be")
            elem = element_t(name="beryllium", symbol="Be", mass=9.012182_r8, n_protons=4)
        case("boron", "B")
            elem = element_t(name="boron", symbol="B", mass=10.811_r8, n_protons=5)
        case("carbon", "C")
            elem = element_t(name="carbon", symbol="C", mass=12.0107_r8, n_protons=6)
        case("nitrogen", "N")
            elem = element_t(name="nitrogen", symbol="N", mass=14.0067_r8, n_protons=7)
        case("oxygen", "O")
            elem = element_t(name="oxygen", symbol="O", mass=15.9994_r8, n_protons=8)
        case("fluorine", "F")
            elem = element_t(name="fluorine", symbol="F", mass=18.9984032_r8, n_protons=9)
        case("neon", "Ne")
            elem = element_t(name="neon", symbol="Ne", mass=20.1797_r8, n_protons=10)
        case("sodium", "Na")
            elem = element_t(name="sodium", symbol="Na", mass=22.98977_r8, n_protons=11)
        case("magnesium", "Mg")
            elem = element_t(name="magnesium", symbol="Mg", mass=24.305_r8, n_protons=12)
        case("aluminum", "Al")
            elem = element_t(name="aluminum", symbol="Al", mass=26.981538_r8, n_protons=13)
        case("silicon", "Si")
            elem = element_t(name="silicon", symbol="Si", mass=28.0855_r8, n_protons=14)
        case("phosphorus", "P")
            elem = element_t(name="phosphorus", symbol="P", mass=30.973761_r8, n_protons=15)
        case("sulfur", "S")
            elem = element_t(name="sulfur", symbol="S", mass=32.065_r8, n_protons=16)
        case("chlorine", "Cl")
            elem = element_t(name="chlorine", symbol="Cl", mass=35.453_r8, n_protons=17)
        case("argon", "Ar")
            elem = element_t(name="argon", symbol="Ar", mass=39.948_r8, n_protons=18)
        case("potassium", "K")
            elem = element_t(name="potassium", symbol="K", mass=39.0983_r8, n_protons=19)
        case("calcium", "Ca")
            elem = element_t(name="calcium", symbol="Ca", mass=40.078_r8, n_protons=20)
        case("scandium", "Sc")
            elem = element_t(name="scandium", symbol="Sc", mass=44.95591_r8, n_protons=21)
        case("titanium", "Ti")
            elem = element_t(name="titanium", symbol="Ti", mass=47.867_r8, n_protons=22)
        case("vanadium", "V")
            elem = element_t(name="vanadium", symbol="V", mass=50.9415_r8, n_protons=23)
        case("chromium", "Cr")
            elem = element_t(name="chromium", symbol="Cr", mass=51.9961_r8, n_protons=24)
        case("manganese", "Mn")
            elem = element_t(name="manganese", symbol="Mn", mass=54.938049_r8, n_protons=25)
        case("iron", "Fe")
            elem = element_t(name="iron", symbol="Fe", mass=55.845_r8, n_protons=26)
        case("cobalt", "Co")
            elem = element_t(name="cobalt", symbol="Co", mass=58.9332_r8, n_protons=27)
        case("nickel", "Ni")
            elem = element_t(name="nickel", symbol="Ni", mass=58.6934_r8, n_protons=28)
        case("copper", "Cu")
            elem = element_t(name="copper", symbol="Cu", mass=63.546_r8, n_protons=29)
        case("zinc", "Zn")
            elem = element_t(name="zinc", symbol="Zn", mass=65.409_r8, n_protons=30)
        case("gallium", "Ga")
            elem = element_t(name="gallium", symbol="Ga", mass=69.723_r8, n_protons=31)
        case("germanium", "Ge")
            elem = element_t(name="germanium", symbol="Ge", mass=72.64_r8, n_protons=32)
        case("arsenic", "As")
            elem = element_t(name="arsenic", symbol="As", mass=74.9216_r8, n_protons=33)
        case("selenium", "Se")
            elem = element_t(name="selenium", symbol="Se", mass=78.96_r8, n_protons=34)
        case("bromine", "Br")
            elem = element_t(name="bromine", symbol="Br", mass=79.904_r8, n_protons=35)
        case("krypton", "Kr")
            elem = element_t(name="krypton", symbol="Kr", mass=83.798_r8, n_protons=36)
        case("rubidium", "Rb")
            elem = element_t(name="rubidium", symbol="Rb", mass=85.4678_r8, n_protons=37)
        case("strontium", "Sr")
            elem = element_t(name="strontium", symbol="Sr", mass=87.62_r8, n_protons=38)
        case("yttrium", "Y")
            elem = element_t(name="yttrium", symbol="Y", mass=88.90585_r8, n_protons=39)
        case("zirconium", "Zr")
            elem = element_t(name="zirconium", symbol="Zr", mass=91.224_r8, n_protons=40)
        case("niobium", "Nb")
            elem = element_t(name="niobium", symbol="Nb", mass=92.90638_r8, n_protons=41)
        case("molybdenum", "Mo")
            elem = element_t(name="molybdenum", symbol="Mo", mass=95.94_r8, n_protons=42)
        case("technetium", "Tc")
            elem = element_t(name="technetium", symbol="Tc", mass=98_r8, n_protons=43)
        case("ruthenium", "Ru")
            elem = element_t(name="ruthenium", symbol="Ru", mass=101.07_r8, n_protons=44)
        case("rhodium", "Rh")
            elem = element_t(name="rhodium", symbol="Rh", mass=102.9055_r8, n_protons=45)
        case("palladium", "Pd")
            elem = element_t(name="palladium", symbol="Pd", mass=106.42_r8, n_protons=46)
        case("silver", "Ag")
            elem = element_t(name="silver", symbol="Ag", mass=107.8682_r8, n_protons=47)
        case("cadmium", "Cd")
            elem = element_t(name="cadmium", symbol="Cd", mass=112.411_r8, n_protons=48)
        case("indium", "In")
            elem = element_t(name="indium", symbol="In", mass=114.818_r8, n_protons=49)
        case("tin", "Sn")
            elem = element_t(name="tin", symbol="Sn", mass=118.71_r8, n_protons=50)
        case("antimony", "Sb")
            elem = element_t(name="antimony", symbol="Sb", mass=121.76_r8, n_protons=51)
        case("tellurium", "Te")
            elem = element_t(name="tellurium", symbol="Te", mass=127.6_r8, n_protons=52)
        case("iodine", "I")
            elem = element_t(name="iodine", symbol="I", mass=126.90447_r8, n_protons=53)
        case("xenon", "Xe")
            elem = element_t(name="xenon", symbol="Xe", mass=131.293_r8, n_protons=54)
        case("cesium", "Cs")
            elem = element_t(name="cesium", symbol="Cs", mass=132.90545_r8, n_protons=55)
        case("barium", "Ba")
            elem = element_t(name="barium", symbol="Ba", mass=137.327_r8, n_protons=56)
        case("lanthanum", "La")
            elem = element_t(name="lanthanum", symbol="La", mass=138.9055_r8, n_protons=57)
        case("cerium", "Ce")
            elem = element_t(name="cerium", symbol="Ce", mass=140.116_r8, n_protons=58)
        case("praseodymium", "Pr")
            elem = element_t(name="praseodymium", symbol="Pr", mass=140.90765_r8, n_protons=59)
        case("neodymium", "Nd")
            elem = element_t(name="neodymium", symbol="Nd", mass=144.24_r8, n_protons=60)
        case("promethium", "Pm")
            elem = element_t(name="promethium", symbol="Pm", mass=145_r8, n_protons=61)
        case("samarium", "Sm")
            elem = element_t(name="samarium", symbol="Sm", mass=150.36_r8, n_protons=62)
        case("europium", "Eu")
            elem = element_t(name="europium", symbol="Eu", mass=151.964_r8, n_protons=63)
        case("gadolinium", "Gd")
            elem = element_t(name="gadolinium", symbol="Gd", mass=157.25_r8, n_protons=64)
        case("terbium", "Tb")
            elem = element_t(name="terbium", symbol="Tb", mass=158.92534_r8, n_protons=65)
        case("dysprosium", "Dy")
            elem = element_t(name="dysprosium", symbol="Dy", mass=162.5_r8, n_protons=66)
        case("holmium", "Ho")
            elem = element_t(name="holmium", symbol="Ho", mass=164.93032_r8, n_protons=67)
        case("erbium", "Er")
            elem = element_t(name="erbium", symbol="Er", mass=167.259_r8, n_protons=68)
        case("thulium", "Tm")
            elem = element_t(name="thulium", symbol="Tm", mass=168.93421_r8, n_protons=69)
        case("ytterbium", "Yb")
            elem = element_t(name="ytterbium", symbol="Yb", mass=173.04_r8, n_protons=70)
        case("lutetium", "Lu")
            elem = element_t(name="lutetium", symbol="Lu", mass=174.967_r8, n_protons=71)
        case("hafnium", "Hf")
            elem = element_t(name="hafnium", symbol="Hf", mass=178.49_r8, n_protons=72)
        case("tantalum", "Ta")
            elem = element_t(name="tantalum", symbol="Ta", mass=180.9479_r8, n_protons=73)
        case("tungsten", "W")
            elem = element_t(name="tungsten", symbol="W", mass=183.84_r8, n_protons=74)
        case("rhenium", "Re")
            elem = element_t(name="rhenium", symbol="Re", mass=186.207_r8, n_protons=75)
        case("osmium", "Os")
            elem = element_t(name="osmium", symbol="Os", mass=190.23_r8, n_protons=76)
        case("iridium", "Ir")
            elem = element_t(name="iridium", symbol="Ir", mass=192.217_r8, n_protons=77)
        case("platinum", "Pt")
            elem = element_t(name="platinum", symbol="Pt", mass=195.078_r8, n_protons=78)
        case("gold", "Au")
            elem = element_t(name="gold", symbol="Au", mass=196.96655_r8, n_protons=79)
        case("mercury", "Hg")
            elem = element_t(name="mercury", symbol="Hg", mass=200.59_r8, n_protons=80)
        case("thallium", "Tl")
            elem = element_t(name="thallium", symbol="Tl", mass=204.3833_r8, n_protons=81)
        case("lead", "Pb")
            elem = element_t(name="lead", symbol="Pb", mass=207.2_r8, n_protons=82)
        case("bismuth", "Bi")
            elem = element_t(name="bismuth", symbol="Bi", mass=208.98038_r8, n_protons=83)
        case("polonium", "Po")
            elem = element_t(name="polonium", symbol="Po", mass=209_r8, n_protons=84)
        case("astatine", "At")
            elem = element_t(name="astatine", symbol="At", mass=210_r8, n_protons=85)
        case("radon", "Rn")
            elem = element_t(name="radon", symbol="Rn", mass=222_r8, n_protons=86)
        case("francium", "Fr")
            elem = element_t(name="francium", symbol="Fr", mass=223_r8, n_protons=87)
        case("radium", "Ra")
            elem = element_t(name="radium", symbol="Ra", mass=226_r8, n_protons=88)
        case("actinium", "Ac")
            elem = element_t(name="actinium", symbol="Ac", mass=227_r8, n_protons=89)
        case("thorium", "Th")
            elem = element_t(name="thorium", symbol="Th", mass=232.0381_r8, n_protons=90)
        case("protactinium", "Pa")
            elem = element_t(name="protactinium", symbol="Pa", mass=231.03588_r8, n_protons=91)
        case("uranium", "U")
            elem = element_t(name="uranium", symbol="U", mass=238.02891_r8, n_protons=92)
        case("neptunium", "Np")
            elem = element_t(name="neptunium", symbol="Np", mass=237_r8, n_protons=93)
        case("plutonium", "Pu")
            elem = element_t(name="plutonium", symbol="Pu", mass=244_r8, n_protons=94)
        case("americium", "Am")
            elem = element_t(name="americium", symbol="Am", mass=243_r8, n_protons=95)
        case("curium", "Cm")
            elem = element_t(name="curium", symbol="Cm", mass=247_r8, n_protons=96)
        case("berkelium", "Bk")
            elem = element_t(name="berkelium", symbol="Bk", mass=247_r8, n_protons=97)
        case("californium", "Cf")
            elem = element_t(name="californium", symbol="Cf", mass=251_r8, n_protons=98)
        case("einsteinium", "Es")
            elem = element_t(name="einsteinium", symbol="Es", mass=252_r8, n_protons=99)
        case("fermium", "Fm")
            elem = element_t(name="fermium", symbol="Fm", mass=257_r8, n_protons=100)
        case("mendelevium", "Md")
            elem = element_t(name="mendelevium", symbol="Md", mass=258_r8, n_protons=101)
        case("nobelium", "No")
            elem = element_t(name="nobelium", symbol="No", mass=259_r8, n_protons=102)
        case("lawrencium", "Lr")
            elem = element_t(name="lawrencium", symbol="Lr", mass=262_r8, n_protons=103)
        case("rutherfordium", "Rf")
            elem = element_t(name="rutherfordium", symbol="Rf", mass=261_r8, n_protons=104)
        case("dubnium", "Db")
            elem = element_t(name="dubnium", symbol="Db", mass=262_r8, n_protons=105)
        case("seaborgium", "Sg")
            elem = element_t(name="seaborgium", symbol="Sg", mass=266_r8, n_protons=106)
        case("bohrium", "Bh")
            elem = element_t(name="bohrium", symbol="Bh", mass=264_r8, n_protons=107)
        case("hassium", "Hs")
            elem = element_t(name="hassium", symbol="Hs", mass=277_r8, n_protons=108)
        case("meitnerium", "Mt")
            elem = element_t(name="meitnerium", symbol="Mt", mass=268_r8, n_protons=109)
        case("darmstadtium", "Ds")
            elem = element_t(name="darmstadtium", symbol="Ds", mass=281_r8, n_protons=110)
        case("roentgenium", "Rg")
            elem = element_t(name="roentgenium", symbol="Rg", mass=272_r8, n_protons=111)
        case("copernicium", "Cn")
            elem = element_t(name="copernicium", symbol="Cn", mass=285_r8, n_protons=112)
        case("nihonium", "Nh")
            elem = element_t(name="nihonium", symbol="Nh", mass=286_r8, n_protons=113)
        case("flerovium", "Fl")
            elem = element_t(name="flerovium", symbol="Fl", mass=289_r8, n_protons=114)
        case("moscovium", "Mc")
            elem = element_t(name="moscovium", symbol="Mc", mass=289_r8, n_protons=115)
        case("livermorium", "Lv")
            elem = element_t(name="livermorium", symbol="Lv", mass=293_r8, n_protons=116)
        case("tennessine", "Ts")
            elem = element_t(name="tennessine", symbol="Ts", mass=294_r8, n_protons=117)
        case("oganesson", "Og")
            elem = element_t(name="oganesson", symbol="Og", mass=294_r8, n_protons=118)
        case default
            error stop "Unknown element or symbol"
        end select

        atomicnumber = elem%n_protons
        mass = elem%mass

    end subroutine getElementInfo
end module elements