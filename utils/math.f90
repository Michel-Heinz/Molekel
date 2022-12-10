module math

    use kinds
    use global
    use libBoys_data

    contains
        function eriIndex(i, j, k, l) result(idx)
            implicit none
            integer, intent(in) :: i, j, k, l
            integer :: ij, kl, idx

            ij = 0
            kl = 0
            idx = 0
            if (i < j)then
                ij = j * (j + 1) / 2 + i
            else
                ij = i * (i + 1) / 2 + j
            end if
            if (k < l) then
                kl = l * (l + 1) / 2 + k
            else
                kl = k * (k + 1) / 2 + l
            end if
            if (ij < kl) then
                idx = (kl * (kl + 1) / 2 + ij)
            else
                idx = (ij * (ij + 1) / 2 + kl)
            end if
        end function eriIndex

        function gaussianProductCenter(a, ACoord, b, BCoord) result(gpc)
            implicit none
            real(r8), intent(in) :: a, b
            real(r8), dimension(3), intent(in) :: ACoord, BCoord
            real(r8), dimension(3) :: gpc

            gpc =  (a*ACoord+b*BCoord)/(a+b)
        end function gaussianProductCenter

        function Boys_func(n, x) result(res)
            !
            ! libboys - a FORTRAN library to numerically calculate the Boys function
            ! Copyright (C) 2014-2016 Michael Böhme <boehme.mic@gmail.com>
            !
            ! This program is free software; you can redistribute it and/or
            ! modify it under the terms of the GNU Lesser General Public
            ! License as published by the Free Software Foundation; either
            ! version 2.1 of the License, or (at your option) any later version.
            !
            ! This program is distributed in the hope that it will be useful,
            ! but WITHOUT ANY WARRANTY; without even the implied warranty of
            ! MDRCHANTABILITY or FITNDSS FOR A PARTICULAR PURPOSD.  See the GNU
            ! Lesser General Public License for more details.
            !
            ! You should have received a copy of the GNU Lesser General Public
            ! License along with this program; if not, write to the Free Software
            ! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
            !

            ! ***************************************************************************************
            ! * Boys function F_n(x)
            ! ***************************************************************************************
            use libboys_data
            implicit none
            !
            integer, intent(in) :: n
            double precision, intent(in) :: x
            double precision :: res, sfac, dx, dxi
            integer :: i, j
            double precision, dimension(6) :: n_fac_dble = (/ 1d0, 2d0, 6d0, 24d0, 120d0, 720d0 /)
            double precision, dimension(31) :: n_fac2_dble = (/ &
                    1d0, 1d0, 1d0, 2d0, 3d0, 8d0, 15d0, 48d0, 105d0, 384d0, 945d0, 3840d0, &
                            10395d0,46080d0,135135d0,645120d0,2027025d0,10321920d0,34459425d0,185794560d0,654729075d0,3715891200d0,13749310575d0,81749606400d0, &
                            316234143225d0,1961990553600d0,7905853580625d0,51011754393600d0,213458046676875d0,1428329123020800d0,6190283353629375d0 /)
            double precision, parameter :: Pi = 3.1415926535897932d0
            double precision, parameter :: eps = 1d-14
            integer, parameter :: MAX_RECURSION = 6
            double precision :: epsrel
            !
            res = 0d0

            if ( n .eq. 0 ) then
                if ( x .lt. eps ) then
                    res = 1d0
                else
                    res = dsqrt( Pi / (4d0*x) ) * derf(dsqrt(x))
                end if
            else
                if ( n .gt. libBoysMaxN ) then
                    write(*,*) "libboys error: not implemented!"
                    return
                end if
                if ( n .lt. 0 ) then
                    write(*,*) "libboys error: n < 0!"
                    return
                end if

                if ( x .lt. eps ) then
                    res = 1d0 / ( 2d0*dble(n) + 1d0 )
                else if ( x .gt. 50d0 ) then
                    res = n_fac2_dble(2*n-1 +2) / 2d0**(n+1) * dsqrt(Pi/x**(2*n+1))
                else

                    if ( x .ge. 10d0 ) then
                        j = int((x-9.95d0)*10d0) + 1
                        dx = BoysFuncValuesL(j, 1) - x
                        dxi = dx
                        res = BoysFuncValuesL(j, n + 2)
                        epsrel = res * eps
                        do i = 1, MAX_RECURSION
                            sfac = BoysFuncValuesL(j, n + 2 + i) * dxi / n_fac_dble(i)
                            res = res + sfac
                            if ( abs(sfac) .lt. epsrel ) then
                                return
                            end if
                            dxi = dxi * dx
                        end do

                    else if ( x .ge. 5d0 ) then
                        j = int((x-4.975d0)*20d0) + 1
                        dx = BoysFuncValuesM(j, 1) - x
                        dxi = dx
                        res = BoysFuncValuesM(j, n + 2)
                        epsrel = res * eps
                        do i = 1, MAX_RECURSION
                            sfac = BoysFuncValuesM(j, n + 2 + i) * dxi / n_fac_dble(i)
                            res = res + sfac
                            if ( abs(sfac) .lt. epsrel ) then
                                return
                            end if
                            dxi = dxi * dx
                        end do

                    else
                        j = int(x*40d0+0.5d0) + 1
                        dx = BoysFuncValuesS(j, 1) - x
                        dxi = dx
                        res = BoysFuncValuesS(j, n + 2)
                        epsrel = res * eps
                        do i = 1, MAX_RECURSION
                            sfac = BoysFuncValuesS(j, n + 2 + i) * dxi / n_fac_dble(i)
                            res = res + sfac
                            if ( abs(sfac) .lt. epsrel ) then
                                return
                            end if
                            dxi = dxi * dx
                        end do

                    end if
                end if
            end if

        end function


        recursive function fact2(n) result(fact2Result)

        integer, value :: n
        integer :: fact2Result

        if (n == 0 .or. n == 1 .or. n == -1) then
            fact2Result = 1
            return
        end if
        fact2result = n * fact2(n - 2)

        end function fact2
end module math