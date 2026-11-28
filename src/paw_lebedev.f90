!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2026 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief   Generation of the spherical Lebedev grids.
!>          The original 14 grids were generated with a precision of at least
!>          33 digits (Mathematica). Thus the Lebedev abscissas and weights
!>          are accurate to 32 digits (quadruple precision).
!> \version 1.0
!>
!> \par Literature
!>          - V. I. Lebedev, Zh. Vychisl. Mat. Mat. Fiz. 15, 48 (1975)
!>          - V. I. Lebedev, Zh. Vychisl. Mat. Mat. Fiz. 16, 293 (1976)
!>          - V. I. Lebedev, Sibirsk Mat. Zh. 18, 132 (1977)
!>          - V. I. Lebedev and A. L. Skorokhodov, Russ. Acad. Sci. Dokl.
!>            Math. 45, 587 (1992)
!>          - V. I. Lebedev, Russ. Acad. Sci. Dokl. Math. 50, 283 (1995)
!> \par Notes
!>          Implemented Lebedev grids:
!>            No.    l    n(nonred.)    n(total)
!>             1     3         1             6
!>             2     5         2            14
!>             3     7         3            26
!>             4     9         3            38
!>             5    11         4            50
!>             6    15         5            86
!>             7    17         6           110
!>             8    19         7           146
!>             9    23         9           194
!>            10    29        12           302
!>            11    35        16           434
!>            12    41        20           590
!>            13    47        25           770
!>            14    53        30           974
!>            15    59        36          1202
!>            16    65        42          1454
!>            17    95        81          3074
!>            18   131       144          5810
!>          The rules above order 53 use 16-digit Lebedev-Laikov coefficients;
!>          the 32-digit precision statement above applies to rules 1-14.
!>          Source/attribution: tests/unittests/skala_reconstruction/LEBEDEV.md.
!> \par Variables:
!>        - lebedev_grid: Lebedev grids.
!>                        l: Angular momentum quantum number l.
!>                        n: Number of grid points.
!>                        w: Grid point weights.
!>                        r: Abscissas (grid point coordinates)
!>        - max_np      : Maximum number of nonredundant grid points.
!>        - na1         : Number of point for subsystem A1.
!>        - na2         : Number of point for subsystem A2.
!>        - na3         : Number of point for subsystem A3.
!>        - nb          : Number of point for subsystem B.
!>        - nc          : Number of point for subsystem C.
!>        - nd          : Number of point for subsystem D.
!>        - nlg         : Number of implemented Lebedev grids.
!>        - nlgp        : Grid point counter.
!>        - r           : Working storage (grid point coordinates).
!>        - w           : Working storage (grid point weights).
! **************************************************************************************************
MODULE PAW_LEBEDEV_MODULE
   IMPLICIT NONE

   SAVE

   PRIVATE

   INTEGER, PARAMETER  :: na1 = 6, &
                          na2 = 12, &
                          na3 = 8, &
                          nb = 24, &
                          nc = 24, &
                          nd = 48, &
                          nlg = 18, &
                          max_np = 36

! **************************************************************************************************
   TYPE oh_grid
      INTEGER                           :: l = -1, n = -1
      REAL(8), DIMENSION(:), POINTER   :: w => NULL()
      REAL(8), DIMENSION(:, :), POINTER :: r => NULL()
   END TYPE oh_grid

   TYPE(oh_grid), DIMENSION(nlg), TARGET :: lebedev_grid = oh_grid()

   REAL(8), PARAMETER :: one = 1.0_8, rs2 = one/SQRT(2.0_8), rs3 = one/SQRT(3.0_8), zero = 0.0_8
   INTEGER  :: nlgp = -1
   LOGICAL  :: init_lebedev_grids_done = .FALSE.

   REAL(8), DIMENSION(max_np/3) :: w = 0.0_8
   REAL(8), DIMENSION(max_np)   :: r = 0.0_8

   PUBLIC :: paw_lebedev_get, paw_lebedev_size, paw_lebedev_test
   INTEGER, PARAMETER, PUBLIC :: PAW_LEBEDEV_MAX_EXACTNESS = 131

CONTAINS

! **************************************************************************************************
!> \brief Return the size of the smallest grid exact through lmin.
! **************************************************************************************************
   SUBROUTINE paw_lebedev_size(lmin, npoint, lexact)
      INTEGER, INTENT(IN)  :: lmin
      INTEGER, INTENT(OUT) :: npoint
      INTEGER, INTENT(OUT), OPTIONAL :: lexact
      INTEGER              :: ilg

      ilg = get_number_of_lebedev_grid(l=lmin)
      npoint = lebedev_grid(ilg)%n
      IF (PRESENT(lexact)) lexact = lebedev_grid(ilg)%l
   END SUBROUTINE paw_lebedev_size

! **************************************************************************************************
!> \brief Copy the smallest grid exact through lmin.
!> \note The angular weights sum to one. Multiply them by 4*pi for a solid-angle integral.
! **************************************************************************************************
   SUBROUTINE paw_lebedev_get(lmin, npoint, direction, weight)
      INTEGER, INTENT(IN) :: lmin, npoint
      REAL(8), INTENT(OUT) :: direction(3, npoint), weight(npoint)
      INTEGER              :: ilg

      ilg = get_number_of_lebedev_grid(l=lmin)
      IF (npoint /= lebedev_grid(ilg)%n) THEN
         CALL paw_lebedev_abort("Inconsistent Lebedev grid size")
      END IF
      direction = lebedev_grid(ilg)%r
      weight = lebedev_grid(ilg)%w
   END SUBROUTINE paw_lebedev_get

! **************************************************************************************************
!> \brief Check normalization, unit radii, and second angular moments.
! **************************************************************************************************
   SUBROUTINE paw_lebedev_test(maxerror)
      REAL(8), INTENT(OUT) :: maxerror
      INTEGER              :: i, ilg, ip
      REAL(8)              :: moment

      IF (.NOT. init_lebedev_grids_done) CALL init_lebedev_grids
      maxerror = 0.0_8
      DO ilg = 1, nlg
         maxerror = MAX(maxerror, ABS(SUM(lebedev_grid(ilg)%w) - 1.0_8))
         DO ip = 1, lebedev_grid(ilg)%n
            maxerror = MAX(maxerror, &
               ABS(SUM(lebedev_grid(ilg)%r(:, ip)**2) - 1.0_8))
         END DO
         DO i = 1, 3
            moment = SUM(lebedev_grid(ilg)%w*lebedev_grid(ilg)%r(i, :)**2)
            maxerror = MAX(maxerror, ABS(moment - 1.0_8/3.0_8))
         END DO
      END DO
   END SUBROUTINE paw_lebedev_test

! **************************************************************************************************
   SUBROUTINE paw_lebedev_abort(message)
      CHARACTER(*), INTENT(IN) :: message

      CALL ERROR$MSG(TRIM(message))
      CALL ERROR$STOP('PAW_LEBEDEV')
   END SUBROUTINE paw_lebedev_abort

! **************************************************************************************************
!> \brief  Get the number of the Lebedev grid, which has the requested
!>          angular momentum quantnum number l or size n.
!> \param l ...
!> \param n ...
!> \return ...
!> \date   05.05.99
!> \par Variables
!> \author Matthias Krack
!> \version 1.0
!> \note         - l : Minimal angular momentum quantum number of the requested Lebedev grid.
!> \note         - n : Minimal size of the requested Lebedev grid.
! **************************************************************************************************
   FUNCTION get_number_of_lebedev_grid(l, n) RESULT(number_of_lebedev_grid)

      INTEGER, INTENT(IN), OPTIONAL                      :: l, n
      INTEGER                                            :: number_of_lebedev_grid

      INTEGER                                            :: i
      INTEGER, ALLOCATABLE, DIMENSION(:)                 :: ll, nn
      INTEGER, DIMENSION(1)                              :: lgnum

      IF (.NOT. init_lebedev_grids_done) CALL init_lebedev_grids

      lgnum(1) = 0
      ALLOCATE (nn(SIZE(lebedev_grid, 1)), ll(SIZE(lebedev_grid, 1)))
      DO i = 1, SIZE(nn, 1)
         nn(i) = lebedev_grid(i)%n
         ll(i) = lebedev_grid(i)%l
      END DO
      IF (PRESENT(l)) THEN
         IF (l < 1 .OR. l > PAW_LEBEDEV_MAX_EXACTNESS) THEN
            CALL paw_lebedev_abort("Lebedev minimum exactness must lie in 1..131")
         END IF
!      lgnum(:) = MINLOC(lebedev_grid(:)%n,MASK=(lebedev_grid(:)%l >= l))
         lgnum(:) = MINLOC(nn(:), MASK=(ll(:) >= l))
      ELSE IF (PRESENT(n)) THEN
!      lgnum(:) = MINLOC(lebedev_grid(:)%n,MASK=(lebedev_grid(:)%n >= n))
         lgnum(:) = MINLOC(nn(:), MASK=(nn(:) >= n))
      ELSE
         CALL paw_lebedev_abort("*** ERROR in FUNCTION get_number_of_lebedev_grid ***")
      END IF

      IF (lgnum(1) == 0) THEN
         CALL paw_lebedev_abort("*** ERROR in FUNCTION get_number_of_lebedev_grid ***")
      ELSE
         number_of_lebedev_grid = lgnum(1)
      END IF
      DEALLOCATE (nn, ll)
   END FUNCTION get_number_of_lebedev_grid

! **************************************************************************************************
!> \brief  Load the Lebedev grid points of a subsystem.
!> \param subsystem ...
!> \param lgnum ...
!> \param np ...
!> \date    06.05.99
!> \par Variables
!>       - subsystem: String, which specifies the selected subsystem.
!>       - lgnum    : Number of the current Lebedev grid.
!>       - np       : Number of the nonredundant grid points.
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE load_sub_grid(subsystem, lgnum, np)
      CHARACTER(*), INTENT(IN)                           :: subsystem
      INTEGER, INTENT(IN)                                :: lgnum, np

      INTEGER                                            :: i, j
      REAL(8)                                      :: x, y, z

!   *** Check argument values ***

      IF ((lgnum < 1) .OR. (lgnum > nlg)) THEN
         CALL paw_lebedev_abort("Argument value #2 (lgnum) is out of range")
      END IF

      IF ((np < 0) .OR. (np > max_np)) THEN
         CALL paw_lebedev_abort("Argument value #3 (np) is out of range")
      END IF

      SELECT CASE (subsystem)
      CASE ("A1")
         lebedev_grid(lgnum)%w(1:na1) = w(1)
         lebedev_grid(lgnum)%r(1:3, 1:na1) = &
            RESHAPE([zero, zero, one, &
                     zero, zero, -one, &
                     zero, one, zero, &
                     zero, -one, zero, &
                     one, zero, zero, &
                     -one, zero, zero], [3, na1])
         nlgp = na1
      CASE ("A2")
         lebedev_grid(lgnum)%w(nlgp + 1:nlgp + na2) = w(1)
         lebedev_grid(lgnum)%r(1:3, nlgp + 1:nlgp + na2) = &
            RESHAPE([zero, rs2, rs2, &
                     zero, rs2, -rs2, &
                     zero, -rs2, rs2, &
                     zero, -rs2, -rs2, &
                     rs2, zero, rs2, &
                     rs2, zero, -rs2, &
                     -rs2, zero, rs2, &
                     -rs2, zero, -rs2, &
                     rs2, rs2, zero, &
                     rs2, -rs2, zero, &
                     -rs2, rs2, zero, &
                     -rs2, -rs2, zero], [3, na2])
         nlgp = nlgp + na2
      CASE ("A3")
         lebedev_grid(lgnum)%w(nlgp + 1:nlgp + na3) = w(1)
         lebedev_grid(lgnum)%r(1:3, nlgp + 1:nlgp + na3) = &
            RESHAPE([rs3, rs3, rs3, &
                     rs3, rs3, -rs3, &
                     rs3, -rs3, rs3, &
                     -rs3, rs3, rs3, &
                     rs3, -rs3, -rs3, &
                     -rs3, rs3, -rs3, &
                     -rs3, -rs3, rs3, &
                     -rs3, -rs3, -rs3], [3, na3])
         nlgp = nlgp + na3
      CASE ("B")
         DO i = 1, np
            x = r(i)
            y = rs2*SQRT(one - x**2)
            lebedev_grid(lgnum)%w(nlgp + nb*(i - 1) + 1:nlgp + nb*i) = w(i)
            lebedev_grid(lgnum)%r(1:3, nlgp + nb*(i - 1) + 1:nlgp + nb*i) = &
               RESHAPE([x, y, y, &
                        x, y, -y, &
                        x, -y, y, &
                        x, -y, -y, &
                        -x, y, y, &
                        -x, y, -y, &
                        -x, -y, y, &
                        -x, -y, -y, &
                        y, x, y, &
                        y, x, -y, &
                        -y, x, y, &
                        -y, x, -y, &
                        y, -x, y, &
                        y, -x, -y, &
                        -y, -x, y, &
                        -y, -x, -y, &
                        y, y, x, &
                        y, -y, x, &
                        -y, y, x, &
                        -y, -y, x, &
                        y, y, -x, &
                        y, -y, -x, &
                        -y, y, -x, &
                        -y, -y, -x], [3, nb])
         END DO
         nlgp = nlgp + nb*np
      CASE ("C")
         DO i = 1, np
            x = r(i)
            y = SQRT(one - x**2)
            lebedev_grid(lgnum)%w(nlgp + nc*(i - 1) + 1:nlgp + nc*i) = w(i)
            lebedev_grid(lgnum)%r(1:3, nlgp + nc*(i - 1) + 1:nlgp + nc*i) = &
               RESHAPE([x, y, zero, &
                        x, -y, zero, &
                        -x, y, zero, &
                        -x, -y, zero, &
                        y, x, zero, &
                        y, -x, zero, &
                        -y, x, zero, &
                        -y, -x, zero, &
                        x, zero, y, &
                        x, zero, -y, &
                        -x, zero, y, &
                        -x, zero, -y, &
                        y, zero, x, &
                        y, zero, -x, &
                        -y, zero, x, &
                        -y, zero, -x, &
                        zero, x, y, &
                        zero, x, -y, &
                        zero, -x, y, &
                        zero, -x, -y, &
                        zero, y, x, &
                        zero, y, -x, &
                        zero, -y, x, &
                        zero, -y, -x], [3, nc])
         END DO
         nlgp = nlgp + nc*np
      CASE ("D")
         IF (MODULO(np, 3) == 0) THEN
            DO i = 1, np, 3
               j = (i + 2)/3
               x = r(i)
               y = r(i + 1)
               z = r(i + 2)
               lebedev_grid(lgnum)%w(nlgp + nd*(j - 1) + 1:nlgp + nd*j) = w(j)
               lebedev_grid(lgnum)%r(1:3, nlgp + nd*(j - 1) + 1:nlgp + nd*j) = &
                  RESHAPE([x, y, z, x, y, -z, x, -y, z, -x, y, z, &
                           x, -y, -z, -x, y, -z, -x, -y, z, -x, -y, -z, &
                           x, z, y, x, z, -y, x, -z, y, -x, z, y, &
                           x, -z, -y, -x, z, -y, -x, -z, y, -x, -z, -y, &
                           y, x, z, y, x, -z, y, -x, z, -y, x, z, &
                           y, -x, -z, -y, x, -z, -y, -x, z, -y, -x, -z, &
                           y, z, x, y, z, -x, y, -z, x, -y, z, x, &
                           y, -z, -x, -y, z, -x, -y, -z, x, -y, -z, -x, &
                           z, x, y, z, x, -y, z, -x, y, -z, x, y, &
                           z, -x, -y, -z, x, -y, -z, -x, y, -z, -x, -y, &
                           z, y, x, z, y, -x, z, -y, x, -z, y, x, &
                           z, -y, -x, -z, y, -x, -z, -y, x, -z, -y, -x], [3, nd])
            END DO
         ELSE
            CALL paw_lebedev_abort("Subsytem D: np is not modulo 3 (check argument #3)")
         END IF
         nlgp = nlgp + nd*np/3
      CASE DEFAULT
         CALL paw_lebedev_abort("The invalid subsystem <"//TRIM(subsystem)// &
                                "> was specified (check argument #1)")
      END SELECT

!   *** Reinitialize working storage ***

      w(:) = zero
      r(:) = zero

   END SUBROUTINE load_sub_grid

! **************************************************************************************************
!> \brief ...
! **************************************************************************************************
   SUBROUTINE deallocate_lebedev_grids()

      INTEGER                                            :: ilg

      IF (.NOT. init_lebedev_grids_done) RETURN
      init_lebedev_grids_done = .FALSE.

      DO ilg = 1, nlg
         DEALLOCATE (lebedev_grid(ilg)%r, lebedev_grid(ilg)%w)
      END DO

   END SUBROUTINE deallocate_lebedev_grids

   ! Map the published orbit codes onto the existing symmetry generator.
   SUBROUTINE load_laikov_orbit(lgnum, code, a, b, weight)
      INTEGER, INTENT(IN) :: lgnum, code
      REAL(8), INTENT(IN) :: a, b, weight
      w(1) = weight
      SELECT CASE (code)
      CASE (0)
         CALL load_sub_grid("A1", lgnum, 0)
      CASE (1)
         CALL load_sub_grid("A2", lgnum, 0)
      CASE (2)
         CALL load_sub_grid("A3", lgnum, 0)
      CASE (3)
         ! CP2K's B orbit takes the unequal coordinate of (a,a,c).
         r(1) = SQRT(1.0_8 - 2.0_8*a*a)
         CALL load_sub_grid("B", lgnum, 1)
      CASE (4)
         r(1) = a
         CALL load_sub_grid("C", lgnum, 1)
      CASE (5)
         r(1:3) = [a, b, SQRT(1.0_8 - a*a - b*b)]
         CALL load_sub_grid("D", lgnum, 3)
      CASE DEFAULT
         CALL paw_lebedev_abort("Invalid Lebedev-Laikov orbit type")
      END SELECT
   END SUBROUTINE load_laikov_orbit

! **************************************************************************************************
!> \brief Load the coordinates and weights of the nonredundant Lebedev grid
!>         points.
!> \date    04.05.99
!> \author  Matthias Krack
!> \version 1.0
! **************************************************************************************************
   SUBROUTINE init_lebedev_grids()
      INTEGER                                            :: ilg

      DO ilg = 1, nlg
         IF (ASSOCIATED(lebedev_grid(ilg)%r)) DEALLOCATE (lebedev_grid(ilg)%r)
         IF (ASSOCIATED(lebedev_grid(ilg)%w)) DEALLOCATE (lebedev_grid(ilg)%w)
      END DO

!   *** Load the angular momentum quantum numbers l of the Lebedev grids ***

      lebedev_grid(1:nlg)%l = [3, 5, 7, 9, 11, 15, 17, 19, 23, 29, 35, 41, 47, 53, 59, 65, 95, 131]

!   *** Load the total number of grid points for each Lebedev grid ***

      lebedev_grid(1:nlg)%n = [6, 14, 26, 38, 50, 86, 110, 146, 194, 302, 434, 590, 770, &
                               974, 1202, 1454, 3074, 5810]

!   *** Allocate storage for the Lebedev grids ***

      DO ilg = 1, nlg
         ALLOCATE (lebedev_grid(ilg)%r(3, lebedev_grid(ilg)%n), &
                   lebedev_grid(ilg)%w(lebedev_grid(ilg)%n))
      END DO

!   *** Load the coordinates and weights of the Lebedev grid points ***

!   *** 1. l = 3 (6 points, octahedron) ***

      w(1) = 1.0_8/6.0_8
      CALL load_sub_grid("A1", 1, 0)

!   *** 2. l = 5 (14 points, capped octahedron) ***

      w(1) = 1.0_8/15.0_8
      CALL load_sub_grid("A1", 2, 0)

      w(1) = 3.0_8/40.0_8
      CALL load_sub_grid("A3", 2, 0)

!   *** 3. l = 7 (26 points) ***

      w(1) = 1.0_8/21.0_8
      CALL load_sub_grid("A1", 3, 0)

      w(1) = 4.0_8/105.0_8
      CALL load_sub_grid("A2", 3, 0)

      w(1) = 9.0_8/280.0_8
      CALL load_sub_grid("A3", 3, 0)

!   *** 4. l = 9 (38 points) ***

      w(1) = 1.0_8/105.0_8
      CALL load_sub_grid("A1", 4, 0)

      w(1) = 9.0_8/280.0_8
      CALL load_sub_grid("A3", 4, 0)

      w(1) = 1.0_8/35.0_8
      r(1) = rs2*SQRT(1.0_8 + rs3)
      CALL load_sub_grid("C", 4, 1)

!   *** 5. l = 11 (50 points) ***

      w(1) = 4.0_8/315.0_8
      CALL load_sub_grid("A1", 5, 0)

      w(1) = 64.0_8/2835.0_8
      CALL load_sub_grid("A2", 5, 0)

      w(1) = 27.0_8/1280.0_8
      CALL load_sub_grid("A3", 5, 0)

      w(1) = 14641.0_8/725760.0_8
      r(1) = 3.0_8/SQRT(11.0_8)
      CALL load_sub_grid("B", 5, 1)

!   *** 6. l = 15 (86 points) ***

      w(1) = 1.15440115440115440115440115440115E-2_8
      CALL load_sub_grid("A1", 6, 0)

      w(1) = 1.19439090858562823236989259736470E-2_8
      CALL load_sub_grid("A3", 6, 0)

      w(1) = 1.11105557106034025109468482160140E-2_8
      r(1) = 8.52518311701267605338736780155357E-1_8
      w(2) = 1.18765012945371420137882805994025E-2_8
      r(2) = 1.89063552885395482707075847005288E-1_8
      CALL load_sub_grid("B", 6, 2)

      w(1) = 1.18123037469044753644792263073650E-2_8
      r(1) = 9.27330657151172465678969739310097E-1_8
      CALL load_sub_grid("C", 6, 1)

!   *** 7. l = 17 (110 points) ***

      w(1) = 3.82827049493716160382827049493716E-3_8
      CALL load_sub_grid("A1", 7, 0)

      w(1) = 9.79373751248751248751248751248751E-3_8
      CALL load_sub_grid("A3", 7, 0)

      w(1) = 8.21173728319111097598993405227308E-3_8
      r(1) = 9.65124035086594105655529546158531E-1_8
      w(2) = 9.59547133607096284945318117290259E-3_8
      r(2) = 8.28769981252592210694031500711749E-1_8
      w(3) = 9.94281489117810328140065828526450E-3_8
      r(3) = 2.15957291845848832354930328945946E-1_8
      CALL load_sub_grid("B", 7, 3)

      w(1) = 9.69499636166302832969499636166303E-3_8
      r(1) = 8.78158910604066133449110592671618E-1_8
      CALL load_sub_grid("C", 7, 1)

!   *** 8. l = 19 (146 points) ***

      w(1) = 5.99631368862138092907323676554446E-4_8
      CALL load_sub_grid("A1", 8, 0)

      w(1) = 7.37299971862075642305743268410561E-3_8
      CALL load_sub_grid("A2", 8, 0)

      w(1) = 7.21051536014448777763305996834282E-3_8
      CALL load_sub_grid("A3", 8, 0)

      w(1) = 7.57439415905403372268748574713806E-3_8
      r(1) = 9.74888643677173235480043489928723E-1_8
      w(2) = 6.75382948631447744073541732486459E-3_8
      r(2) = 8.07089818359582501629574989264033E-1_8
      w(3) = 7.11635549311755538760089284953968E-3_8
      r(3) = 2.91298882209526746288335878313128E-1_8
      CALL load_sub_grid("B", 8, 3)

      w(1) = 6.99108735330326239417148508057599E-3_8
      r(1) = 1.40355381171318328571556780746292E-1_8
      r(2) = 4.49332832326955734884695382705101E-1_8
      r(3) = 8.82270011260322631916665753046583E-1_8
      CALL load_sub_grid("D", 8, 3)

!   *** 9. l = 23 (194 points) ***

      w(1) = 1.78234044724461115736727104869868E-3_8
      CALL load_sub_grid("A1", 9, 0)

      w(1) = 5.71690594997710189299212838832099E-3_8
      CALL load_sub_grid("A2", 9, 0)

      w(1) = 5.57338317884873796836784958446647E-3_8
      CALL load_sub_grid("A3", 9, 0)

      w(1) = 5.51877146727361369172768460119380E-3_8
      r(1) = 7.77493219314767127213777704028801E-1_8
      w(2) = 5.15823771180538310324916154718792E-3_8
      r(2) = 9.12509096867473692992168387555434E-1_8
      w(3) = 5.60870408258799684374936673855184E-3_8
      r(3) = 3.14196994182586079225390955831946E-1_8
      w(4) = 4.10677702816939409072861128564582E-3_8
      r(4) = 9.82972302707253296863729593965148E-1_8
      CALL load_sub_grid("B", 9, 4)

      w(1) = 5.05184606461480847598931196006390E-3_8
      r(1) = 9.38319218137591520905616389195671E-1_8
      CALL load_sub_grid("C", 9, 1)

      w(1) = 5.53024891623309370129768269143303E-3_8
      r(1) = 1.59041710538352952424263581362096E-1_8
      r(2) = 5.25118572443642024905268207753317E-1_8
      r(3) = 8.36036015482458885943746377793309E-1_8
      CALL load_sub_grid("D", 9, 3)

!   *** 10. l = 29 (302 points) ***

      w(1) = 8.54591172512814813423121032618880E-4_8
      CALL load_sub_grid("A1", 10, 0)

      w(1) = 3.59911928502557145886397858961119E-3_8
      CALL load_sub_grid("A3", 10, 0)

      w(1) = 3.65004580767725542865433220112651E-3_8
      r(1) = 1.29238672710514925339493976600550E-1_8
      w(2) = 3.60482260141988171131480913104353E-3_8
      r(2) = 3.71034178384821189405344300671593E-1_8
      w(3) = 3.57672966174336707556208137560857E-3_8
      r(3) = 7.43452042987555751133204566760209E-1_8
      w(4) = 3.44978842430588331001302771048283E-3_8
      r(4) = 8.67643624544083327473890295837421E-1_8
      w(5) = 3.10895312241367525484587698082774E-3_8
      r(5) = 9.49454317226443084214869821724903E-1_8
      w(6) = 2.35210141368916437879217118337424E-3_8
      r(6) = 9.90705621379408123821774729602922E-1_8
      CALL load_sub_grid("B", 10, 6)

      w(1) = 3.60082093221646027279920634177099E-3_8
      r(1) = 8.20326419827759303328870367500864E-1_8
      w(2) = 2.98234496317180385195111046924520E-3_8
      r(2) = 9.64408914879206014987053264857406E-1_8
      CALL load_sub_grid("C", 10, 2)

      w(1) = 3.57154055427338708123297920312395E-3_8
      r(1) = 2.51003475177046506904110494820777E-1_8
      r(2) = 5.44867737258077380269062623720797E-1_8
      r(3) = 8.00072749407395172686232883863933E-1_8
      w(2) = 3.39231220500617018197882653945702E-3_8
      r(4) = 1.23354853258332742165467555277499E-1_8
      r(5) = 4.12772408316853095996383989497138E-1_8
      r(6) = 9.02442529533000401095678759353334E-1_8
      CALL load_sub_grid("D", 10, 6)

!   *** 11. l = 35 (434 points) ***

      w(1) = 5.26589796822443623921598397444095E-4_8
      CALL load_sub_grid("A1", 11, 0)

      w(1) = 2.54821997200260718024899528063707E-3_8
      CALL load_sub_grid("A2", 11, 0)

      w(1) = 2.51231741892730716751285677008112E-3_8
      CALL load_sub_grid("A3", 11, 0)

      w(1) = 1.46249562159461384222611198788353E-3_8
      r(1) = 9.94255912631277812357804047239411E-1_8
      w(2) = 2.01427902091852819554717351102614E-3_8
      r(2) = 9.67987158791472796844440561853288E-1_8
      w(3) = 2.30269478222741577118243755624539E-3_8
      r(3) = 9.14472801120872429366691177286487E-1_8
      w(4) = 2.44537343731298000188170923067793E-3_8
      r(4) = 8.31584400419232294499745594669204E-1_8
      w(5) = 2.50172516840293614549327311336461E-3_8
      r(5) = 7.19016501040843432905788800294849E-1_8
      w(6) = 2.51326717459756436118956600426018E-3_8
      r(6) = 4.07712664897769512483055658274718E-1_8
      w(7) = 2.53040380118635500209853222454729E-3_8
      r(7) = 2.12646824707552073524018894972799E-1_8
      CALL load_sub_grid("B", 11, 7)

      w(1) = 1.91095128217953227363116198076372E-3_8
      r(1) = 9.77642811118264871364737333213224E-1_8
      w(2) = 2.41744237563898077608784837021936E-3_8
      r(2) = 8.81813287779428812899447039158355E-1_8
      CALL load_sub_grid("C", 11, 2)

      w(1) = 2.23660776043784866397235325901501E-3_8
      r(1) = 9.92176963642923726861469072682696E-2_8
      r(2) = 3.34436314534345491435541212141297E-1_8
      r(3) = 9.37180985855372235789161206437430E-1_8
      w(2) = 2.41693004432477530950382860647420E-3_8
      r(4) = 2.05482369640304370826435780923978E-1_8
      r(5) = 4.50233038258262537898138529993494E-1_8
      r(6) = 8.68946032287241197946412945028781E-1_8
      w(3) = 2.51223685456349510228853571762152E-3_8
      r(7) = 1.06801826075804827808519981079378E-1_8
      r(8) = 5.90515704892527108597639202779953E-1_8
      r(9) = 7.99927854385728570494616568630393E-1_8
      w(4) = 2.49664405455308596359511734587523E-3_8
      r(10) = 3.10428403516654146818973681387686E-1_8
      r(11) = 5.55015236107680716116395510647568E-1_8
      r(12) = 7.71746262691590088133470564383333E-1_8
      CALL load_sub_grid("D", 11, 12)

!   *** 12. l = 41 (590 points) ***

      w(1) = 3.09512129530618734224885916984167E-4_8
      CALL load_sub_grid("A1", 12, 0)

      w(1) = 1.85237969859748902097779927120955E-3_8
      CALL load_sub_grid("A3", 12, 0)

      w(1) = 9.76433116505105003063946798281256E-4_8
      r(1) = 9.96278129754016372294779397439010E-1_8
      w(2) = 1.38473723485169190049558481046633E-3_8
      r(2) = 9.78480583762693839428807199465825E-1_8
      w(3) = 1.61721064725441119207240564027320E-3_8
      r(3) = 9.41414158220402542954812123169147E-1_8
      w(4) = 1.74956465728115411858337657547711E-3_8
      r(4) = 8.83078727934132546639472279693012E-1_8
      w(5) = 1.81847177816276877817839645350131E-3_8
      r(5) = 8.02836877335273789544708232847842E-1_8
      w(6) = 1.84671595615124182496431643927999E-3_8
      r(6) = 7.00768575373572920755474238992223E-1_8
      w(7) = 1.85202882829621309713188574624172E-3_8
      r(7) = 4.33373868777154278470427918169993E-1_8
      w(8) = 1.85881258543831701688713685607976E-3_8
      r(8) = 2.70356088359165054450644840302674E-1_8
      w(9) = 1.87179063927774375088198703268059E-3_8
      r(9) = 9.21904070768989460773325829762845E-2_8
      CALL load_sub_grid("B", 12, 9)

      w(1) = 1.30032168588604773254689500223665E-3_8
      r(1) = 9.85013335028001910431569760301413E-1_8
      w(2) = 1.70515399639586401999426890257133E-3_8
      r(2) = 9.18045287711453949033261544688478E-1_8
      w(3) = 1.85716119677407798248150549753942E-3_8
      r(3) = 7.91101929626901988791072435283520E-1_8
      CALL load_sub_grid("C", 12, 3)

      w(1) = 1.55521360339680849695800543399046E-3_8
      r(1) = 8.21302158193251139256011433978433E-2_8
      r(2) = 2.77867319058624428710339115104831E-1_8
      r(3) = 9.57102074310072578548507660672158E-1_8
      w(2) = 1.80223912800852549116687840059028E-3_8
      r(4) = 8.99920584207487492736600381248490E-2_8
      r(5) = 5.03356427107511721798695016575784E-1_8
      r(6) = 8.59379855890721226547417613471625E-1_8
      w(3) = 1.84983056044366016523678105066982E-3_8
      r(7) = 1.81664084036020946133369276322890E-1_8
      r(8) = 5.98412649788537962548021778158698E-1_8
      r(9) = 7.80320742479920327279647082439668E-1_8
      w(4) = 1.71390450710670866804340083221112E-3_8
      r(10) = 1.72079522565687812382959749101837E-1_8
      r(11) = 3.79103540769556328153581995776882E-1_8
      r(12) = 9.09213475092373607255560905376387E-1_8
      w(5) = 1.80265893437745115894760427403985E-3_8
      r(13) = 2.63471665593794963182134911777324E-1_8
      r(14) = 4.74239284255198024276965888711084E-1_8
      r(15) = 8.40047488359050429417090817376105E-1_8
      w(6) = 1.84286647290528563239873317969286E-3_8
      r(16) = 3.51828092773351899746468433493901E-1_8
      r(17) = 5.61026380862206018993247924039616E-1_8
      r(18) = 7.49310611904115932021924350883623E-1_8
      CALL load_sub_grid("D", 12, 18)

!   *** 13. l = 47 (770 points) ***

      w(1) = 2.19294208818118413191132531546907E-4_8
      CALL load_sub_grid("A1", 13, 0)

      w(1) = 1.43643361731907982089311290725143E-3_8
      CALL load_sub_grid("A2", 13, 0)

      w(1) = 1.42194034433587736471699229343586E-3_8
      CALL load_sub_grid("A3", 13, 0)

      w(1) = 6.79812351105050201620332392505941E-4_8
      r(1) = 9.97408677652823043744886496130338E-1_8
      w(2) = 9.91318423529491220896771270468706E-4_8
      r(2) = 9.84799753572301179668683394920869E-1_8
      w(3) = 1.18020783323894879580204544917243E-3_8
      r(3) = 9.58036675983391439732458619060098E-1_8
      w(4) = 1.29659960208092067423594734151509E-3_8
      r(4) = 9.15317950483154831654282453968846E-1_8
      w(5) = 1.36587142742831640181200652406724E-3_8
      r(5) = 8.55901928697886435270454067917142E-1_8
      w(6) = 1.40298860477532526267417750582642E-3_8
      r(6) = 7.79621319527635171330859492071650E-1_8
      w(7) = 1.41864556359560936142032713234525E-3_8
      r(7) = 6.86644447264154195285405155618666E-1_8
      w(8) = 1.42137674185166176103362636565150E-3_8
      r(8) = 4.52311920313658441442123883818590E-1_8
      w(9) = 1.42399647549096162940277119054908E-3_8
      r(9) = 3.12521305001653125265256223110153E-1_8
      w(10) = 1.43155404217856675285551659613462E-3_8
      r(10) = 1.60155803498828974610387120985563E-1_8
      CALL load_sub_grid("B", 13, 10)

      w(1) = 9.25440149986536789398468419487000E-4_8
      r(1) = 9.89477537495598496934466234590453E-1_8
      w(2) = 1.25023999505350931588903916766368E-3_8
      r(2) = 9.40776878793758755393692765525837E-1_8
      w(3) = 1.39436584332923012313782491233921E-3_8
      r(3) = 8.45749305193653306776973233373898E-1_8
      CALL load_sub_grid("C", 13, 3)

      w(1) = 1.12708909467174883469499826293069E-3_8
      r(1) = 6.94402439334941301856132689108761E-2_8
      r(2) = 2.35518789424232641745294975361896E-1_8
      r(3) = 9.69385863498432080572651262482050E-1_8
      w(2) = 1.34575376091067007284708159178371E-3_8
      r(4) = 2.26900410952945985812144550510454E-1_8
      r(5) = 4.10218247404573022502091360294078E-1_8
      r(6) = 8.83310360522112719095559130461882E-1_8
      w(3) = 1.42495728331678280442457616193820E-3_8
      r(7) = 8.02557460777533890494192154117455E-2_8
      r(8) = 6.21430241748160458082227723944198E-1_8
      r(9) = 7.79348105702660973885784203335648E-1_8
      w(4) = 1.26152334123774998512114871120604E-3_8
      r(10) = 1.46799952789657197013820141783396E-1_8
      r(11) = 3.24528434571739443574137218155171E-1_8
      r(12) = 9.34414827052402165196138643907493E-1_8
      w(5) = 1.39254710605269594012443620424449E-3_8
      r(13) = 1.57150776982472710564392240262888E-1_8
      r(14) = 5.22448218969662967424959711591943E-1_8
      r(15) = 8.38064133458312497889043403162912E-1_8
      w(6) = 1.41876167787765638610216672064623E-3_8
      r(16) = 2.36570299315724563080851921006450E-1_8
      r(17) = 6.01754663408955808458505175650313E-1_8
      r(18) = 7.62840624604669826410820436887813E-1_8
      w(7) = 1.33836668447955413932118842030907E-3_8
      r(19) = 7.71481586676573204146393256189282E-2_8
      r(20) = 4.34657551614116278920169057936696E-1_8
      r(21) = 8.97285336132833396499244991766618E-1_8
      w(8) = 1.39370086267613137196129133534721E-3_8
      r(22) = 3.06293666621073020830889612421600E-1_8
      r(23) = 4.90882658903761621546950183900134E-1_8
      r(24) = 8.15609223203975422051589658409773E-1_8
      w(9) = 1.41591475746693200265535028223575E-3_8
      r(25) = 3.82247737952478700050265491776632E-1_8
      r(26) = 5.64876814909950046651375249432877E-1_8
      r(27) = 7.31300793659765701947278263710234E-1_8
      CALL load_sub_grid("D", 13, 27)

!   *** 14. l = 53 (974 points) ***

      w(1) = 1.43829419052743111472634384792800E-4_8
      CALL load_sub_grid("A1", 14, 0)

      w(1) = 1.12577228828700411922446003897444E-3_8
      CALL load_sub_grid("A3", 14, 0)

      w(1) = 4.94802934194924095056628951683154E-4_8
      r(1) = 9.98155345023846501191494527575073E-1_8
      w(2) = 7.35799010912547049057334785150092E-4_8
      r(2) = 9.88883224354685545466876075112182E-1_8
      w(3) = 8.88913277130438430936840534873175E-4_8
      r(3) = 9.68890220434707404463853311788372E-1_8
      w(4) = 9.88834783892143487490379018339493E-4_8
      r(4) = 9.36602730407163187902332491169310E-1_8
      w(5) = 1.05329968170947064968554154238842E-3_8
      r(5) = 8.91267942647606054508357503868536E-1_8
      w(6) = 1.09277880701457848752501645357106E-3_8
      r(6) = 8.32596723702351872809065001851482E-1_8
      w(7) = 1.11438939406322716311544826926779E-3_8
      r(7) = 7.60582905315251495676259855686005E-1_8
      w(8) = 1.12372478805155528926212353749792E-3_8
      r(8) = 6.75400969108414309734349145075511E-1_8
      w(9) = 1.12523932524381359890224350348300E-3_8
      r(9) = 4.66858905695743283152769352096062E-1_8
      w(10) = 1.12615327181590500171763113615557E-3_8
      r(10) = 3.44613654237438223690005816795002E-1_8
      w(11) = 1.13028693112384079609639852411886E-3_8
      r(11) = 2.11954151850184645694140401473565E-1_8
      w(12) = 1.13498653436395488643714589963293E-3_8
      r(12) = 7.16244014499556615845369803814026E-2_8
      CALL load_sub_grid("B", 14, 12)

      w(1) = 6.82336792710993096637138363442284E-4_8
      r(1) = 9.92323565431490196070394244122238E-1_8
      w(2) = 9.45415816044709583567874030695969E-4_8
      r(2) = 9.55781512496548488574596747552425E-1_8
      w(3) = 1.07442997538567914061819358161104E-3_8
      r(3) = 8.82785980701181710576015026411856E-1_8
      w(4) = 1.12930008656913165822354345410170E-3_8
      r(4) = 7.73778447257374736807517446423091E-1_8
      CALL load_sub_grid("C", 14, 4)

      w(1) = 8.43688450090195438498026316671055E-4_8
      r(1) = 5.97404861418134181366741615801962E-2_8
      r(2) = 2.02912875277752280821562682821390E-1_8
      r(3) = 9.77372722845309994986723236843860E-1_8
      w(2) = 1.07525572044888463296011392358171E-3_8
      r(4) = 1.37576040847363647998707229269443E-1_8
      r(5) = 4.60262194248405393969672054729315E-1_8
      r(6) = 8.77058461865802689033970747916718E-1_8
      w(3) = 1.10857723686446203483870096411918E-3_8
      r(7) = 3.39101652633628571879780724877562E-1_8
      r(8) = 5.03067399966203571325278413235681E-1_8
      r(9) = 7.94942299964208514029117644252761E-1_8
      w(4) = 9.56647532378335729987022997768218E-4_8
      r(10) = 1.27167519143981950887694563391455E-1_8
      r(11) = 2.81760642244213432084049569237994E-1_8
      r(12) = 9.51020169374389952612600596224493E-1_8
      w(5) = 1.08066325071739071464463316273343E-3_8
      r(13) = 2.69312074041351249441688496163892E-1_8
      r(14) = 4.33156129172015739916409763707103E-1_8
      r(15) = 8.60143461601761963548074094490711E-1_8
      w(6) = 1.12679713119629459182839576027292E-3_8
      r(16) = 1.41978645260191825726693154370171E-1_8
      r(17) = 6.25616735858081415467234398457286E-1_8
      r(18) = 7.67102186220558360145428572382623E-1_8
      w(7) = 1.02256871535806121135514702022316E-3_8
      r(19) = 6.70928460073825495111114896527700E-2_8
      r(20) = 3.79839521685915669705644712365302E-1_8
      r(21) = 9.22616110730809020138255560099254E-1_8
      w(8) = 1.10896026771310753764243318312336E-3_8
      r(22) = 7.05773818325617225225155364651637E-2_8
      r(23) = 5.51750542142352049132113965320893E-1_8
      r(24) = 8.31017552413474225750006298786621E-1_8
      w(9) = 1.12279065343576582686387381292198E-3_8
      r(25) = 2.78388847788215460217462722823059E-1_8
      r(26) = 6.02961915615918685584436745852250E-1_8
      r(27) = 7.47620610834085792324588316104014E-1_8
      w(10) = 1.03240184711745983700951101978720E-3_8
      r(28) = 1.97957893891740690574322079335028E-1_8
      r(29) = 3.58960632958909578333772677316348E-1_8
      r(30) = 9.12118378409121483660044157731511E-1_8
      w(11) = 1.10724938228385393683388780275631E-3_8
      r(31) = 2.08730706110327399856945819943227E-1_8
      r(32) = 5.34866643813547649971323458627733E-1_8
      r(33) = 8.18748536281021806431067610000242E-1_8
      w(12) = 1.12178004851997205823535727177375E-3_8
      r(34) = 4.05512213787283588596898056051174E-1_8
      r(35) = 5.67499754607437348401393912269671E-1_8
      r(36) = 7.16591845467023718833743633176626E-1_8
      CALL load_sub_grid("D", 14, 36)

      ! Lebedev-Laikov numerical tables, as distributed in PySCF (Apache-2.0).
      ! See LEBEDEV.md for pinned provenance and the original literature.

      ! l = 59, 1202 points.
      CALL load_laikov_orbit(15,0,0.0_8,0.0_8, &
                             0.1105189233267572E-3_8)
      CALL load_laikov_orbit(15,1,0.0_8,0.0_8, &
                             0.9205232738090741E-3_8)
      CALL load_laikov_orbit(15,2,0.0_8,0.0_8, &
                             0.9133159786443561E-3_8)
      CALL load_laikov_orbit(15,3,0.3712636449657089E-1_8,0.0_8, &
                             0.3690421898017899E-3_8)
      CALL load_laikov_orbit(15,3,0.9140060412262223E-1_8,0.0_8, &
                             0.5603990928680660E-3_8)
      CALL load_laikov_orbit(15,3,0.1531077852469906E+0_8,0.0_8, &
                             0.6865297629282609E-3_8)
      CALL load_laikov_orbit(15,3,0.2180928891660612E+0_8,0.0_8, &
                             0.7720338551145630E-3_8)
      CALL load_laikov_orbit(15,3,0.2839874532200175E+0_8,0.0_8, &
                             0.8301545958894795E-3_8)
      CALL load_laikov_orbit(15,3,0.3491177600963764E+0_8,0.0_8, &
                             0.8686692550179628E-3_8)
      CALL load_laikov_orbit(15,3,0.4121431461444309E+0_8,0.0_8, &
                             0.8927076285846890E-3_8)
      CALL load_laikov_orbit(15,3,0.4718993627149127E+0_8,0.0_8, &
                             0.9060820238568219E-3_8)
      CALL load_laikov_orbit(15,3,0.5273145452842337E+0_8,0.0_8, &
                             0.9119777254940867E-3_8)
      CALL load_laikov_orbit(15,3,0.6209475332444019E+0_8,0.0_8, &
                             0.9128720138604181E-3_8)
      CALL load_laikov_orbit(15,3,0.6569722711857291E+0_8,0.0_8, &
                             0.9130714935691735E-3_8)
      CALL load_laikov_orbit(15,3,0.6841788309070143E+0_8,0.0_8, &
                             0.9152873784554116E-3_8)
      CALL load_laikov_orbit(15,3,0.7012604330123631E+0_8,0.0_8, &
                             0.9187436274321654E-3_8)
      CALL load_laikov_orbit(15,4,0.1072382215478166E+0_8,0.0_8, &
                             0.5176977312965694E-3_8)
      CALL load_laikov_orbit(15,4,0.2582068959496968E+0_8,0.0_8, &
                             0.7331143682101417E-3_8)
      CALL load_laikov_orbit(15,4,0.4172752955306717E+0_8,0.0_8, &
                             0.8463232836379928E-3_8)
      CALL load_laikov_orbit(15,4,0.5700366911792503E+0_8,0.0_8, &
                             0.9031122694253992E-3_8)
      CALL load_laikov_orbit(15,5,0.9827986018263947E+0_8,0.1771774022615325E+0_8, &
                             0.6485778453163257E-3_8)
      CALL load_laikov_orbit(15,5,0.9624249230326228E+0_8,0.2475716463426288E+0_8, &
                             0.7435030910982369E-3_8)
      CALL load_laikov_orbit(15,5,0.9402007994128811E+0_8,0.3354616289066489E+0_8, &
                             0.7998527891839054E-3_8)
      CALL load_laikov_orbit(15,5,0.9320822040143202E+0_8,0.3173615246611977E+0_8, &
                             0.8101731497468018E-3_8)
      CALL load_laikov_orbit(15,5,0.9043674199393299E+0_8,0.4090268427085357E+0_8, &
                             0.8483389574594331E-3_8)
      CALL load_laikov_orbit(15,5,0.8912407560074747E+0_8,0.3854291150669224E+0_8, &
                             0.8556299257311812E-3_8)
      CALL load_laikov_orbit(15,5,0.8676435628462708E+0_8,0.4932221184851285E+0_8, &
                             0.8803208679738260E-3_8)
      CALL load_laikov_orbit(15,5,0.8581979986041619E+0_8,0.4785320675922435E+0_8, &
                             0.8811048182425720E-3_8)
      CALL load_laikov_orbit(15,5,0.8396753624049856E+0_8,0.4507422593157064E+0_8, &
                             0.8850282341265444E-3_8)
      CALL load_laikov_orbit(15,5,0.8165288564022188E+0_8,0.5632123020762100E+0_8, &
                             0.9021342299040653E-3_8)
      CALL load_laikov_orbit(15,5,0.8015469370783529E+0_8,0.5434303569693900E+0_8, &
                             0.9010091677105086E-3_8)
      CALL load_laikov_orbit(15,5,0.7773563069070351E+0_8,0.5123518486419871E+0_8, &
                             0.9022692938426915E-3_8)
      CALL load_laikov_orbit(15,5,0.7661621213900394E+0_8,0.6394279634749102E+0_8, &
                             0.9158016174693465E-3_8)
      CALL load_laikov_orbit(15,5,0.7553584143533510E+0_8,0.6269805509024392E+0_8, &
                             0.9131578003189435E-3_8)
      CALL load_laikov_orbit(15,5,0.7344305757559503E+0_8,0.6031161693096310E+0_8, &
                             0.9107813579482705E-3_8)
      CALL load_laikov_orbit(15,5,0.7043837184021765E+0_8,0.5693702498468441E+0_8, &
                             0.9105760258970126E-3_8)
      IF (nlgp /= 1202) CALL paw_lebedev_abort("Invalid Lebedev 1202-point table")

      ! l = 65, 1454 points.
      CALL load_laikov_orbit(16,0,0.0_8,0.0_8, &
                             0.7777160743261247E-4_8)
      CALL load_laikov_orbit(16,2,0.0_8,0.0_8, &
                             0.7557646413004701E-3_8)
      CALL load_laikov_orbit(16,3,0.3229290663413854E-1_8,0.0_8, &
                             0.2841633806090617E-3_8)
      CALL load_laikov_orbit(16,3,0.8036733271462222E-1_8,0.0_8, &
                             0.4374419127053555E-3_8)
      CALL load_laikov_orbit(16,3,0.1354289960531653E+0_8,0.0_8, &
                             0.5417174740872172E-3_8)
      CALL load_laikov_orbit(16,3,0.1938963861114426E+0_8,0.0_8, &
                             0.6148000891358593E-3_8)
      CALL load_laikov_orbit(16,3,0.2537343715011275E+0_8,0.0_8, &
                             0.6664394485800705E-3_8)
      CALL load_laikov_orbit(16,3,0.3135251434752570E+0_8,0.0_8, &
                             0.7025039356923220E-3_8)
      CALL load_laikov_orbit(16,3,0.3721558339375338E+0_8,0.0_8, &
                             0.7268511789249627E-3_8)
      CALL load_laikov_orbit(16,3,0.4286809575195696E+0_8,0.0_8, &
                             0.7422637534208629E-3_8)
      CALL load_laikov_orbit(16,3,0.4822510128282994E+0_8,0.0_8, &
                             0.7509545035841214E-3_8)
      CALL load_laikov_orbit(16,3,0.5320679333566263E+0_8,0.0_8, &
                             0.7548535057718401E-3_8)
      CALL load_laikov_orbit(16,3,0.6172998195394274E+0_8,0.0_8, &
                             0.7554088969774001E-3_8)
      CALL load_laikov_orbit(16,3,0.6510679849127481E+0_8,0.0_8, &
                             0.7553147174442808E-3_8)
      CALL load_laikov_orbit(16,3,0.6777315251687360E+0_8,0.0_8, &
                             0.7564767653292297E-3_8)
      CALL load_laikov_orbit(16,3,0.6963109410648741E+0_8,0.0_8, &
                             0.7587991808518730E-3_8)
      CALL load_laikov_orbit(16,3,0.7058935009831749E+0_8,0.0_8, &
                             0.7608261832033027E-3_8)
      CALL load_laikov_orbit(16,4,0.9955546194091857E+0_8,0.0_8, &
                             0.4021680447874916E-3_8)
      CALL load_laikov_orbit(16,4,0.9734115901794209E+0_8,0.0_8, &
                             0.5804871793945964E-3_8)
      CALL load_laikov_orbit(16,4,0.9275693732388626E+0_8,0.0_8, &
                             0.6792151955945159E-3_8)
      CALL load_laikov_orbit(16,4,0.8568022422795103E+0_8,0.0_8, &
                             0.7336741211286294E-3_8)
      CALL load_laikov_orbit(16,4,0.7623495553719372E+0_8,0.0_8, &
                             0.7581866300989608E-3_8)
      CALL load_laikov_orbit(16,5,0.5707522908892223E+0_8,0.4387028039889501E+0_8, &
                             0.7538257859800743E-3_8)
      CALL load_laikov_orbit(16,5,0.5196463388403083E+0_8,0.3858908414762617E+0_8, &
                             0.7483517247053123E-3_8)
      CALL load_laikov_orbit(16,5,0.4646337531215351E+0_8,0.3301937372343854E+0_8, &
                             0.7371763661112059E-3_8)
      CALL load_laikov_orbit(16,5,0.4063901697557691E+0_8,0.2725423573563777E+0_8, &
                             0.7183448895756934E-3_8)
      CALL load_laikov_orbit(16,5,0.3456329466643087E+0_8,0.2139510237495250E+0_8, &
                             0.6895815529822191E-3_8)
      CALL load_laikov_orbit(16,5,0.2831395121050332E+0_8,0.1555922309786647E+0_8, &
                             0.6480105801792886E-3_8)
      CALL load_laikov_orbit(16,5,0.2197682022925330E+0_8,0.9892878979686097E-1_8, &
                             0.5897558896594636E-3_8)
      CALL load_laikov_orbit(16,5,0.1564696098650355E+0_8,0.4598642910675510E-1_8, &
                             0.5095708849247346E-3_8)
      CALL load_laikov_orbit(16,5,0.6027356673721295E+0_8,0.3376625140173426E+0_8, &
                             0.7536906428909755E-3_8)
      CALL load_laikov_orbit(16,5,0.5496032320255096E+0_8,0.2822301309727988E+0_8, &
                             0.7472505965575118E-3_8)
      CALL load_laikov_orbit(16,5,0.4921707755234567E+0_8,0.2248632342592540E+0_8, &
                             0.7343017132279698E-3_8)
      CALL load_laikov_orbit(16,5,0.4309422998598483E+0_8,0.1666224723456479E+0_8, &
                             0.7130871582177445E-3_8)
      CALL load_laikov_orbit(16,5,0.3664108182313672E+0_8,0.1086964901822169E+0_8, &
                             0.6817022032112776E-3_8)
      CALL load_laikov_orbit(16,5,0.2990189057758436E+0_8,0.5251989784120085E-1_8, &
                             0.6380941145604121E-3_8)
      CALL load_laikov_orbit(16,5,0.6268724013144998E+0_8,0.2297523657550023E+0_8, &
                             0.7550381377920310E-3_8)
      CALL load_laikov_orbit(16,5,0.5707324144834607E+0_8,0.1723080607093800E+0_8, &
                             0.7478646640144802E-3_8)
      CALL load_laikov_orbit(16,5,0.5096360901960365E+0_8,0.1140238465390513E+0_8, &
                             0.7335918720601220E-3_8)
      CALL load_laikov_orbit(16,5,0.4438729938312456E+0_8,0.5611522095882537E-1_8, &
                             0.7110120527658118E-3_8)
      CALL load_laikov_orbit(16,5,0.6419978471082389E+0_8,0.1164174423140873E+0_8, &
                             0.7571363978689501E-3_8)
      CALL load_laikov_orbit(16,5,0.5817218061802611E+0_8,0.5797589531445219E-1_8, &
                             0.7489908329079234E-3_8)
      IF (nlgp /= 1454) CALL paw_lebedev_abort("Invalid Lebedev 1454-point table")

      ! l = 95, 3074 points.
      CALL load_laikov_orbit(17,0,0.0_8,0.0_8, &
                             0.2599095953754734E-4_8)
      CALL load_laikov_orbit(17,1,0.0_8,0.0_8, &
                             0.3603134089687541E-3_8)
      CALL load_laikov_orbit(17,2,0.0_8,0.0_8, &
                             0.3586067974412447E-3_8)
      CALL load_laikov_orbit(17,3,0.1886108518723392E-1_8,0.0_8, &
                             0.9831528474385880E-4_8)
      CALL load_laikov_orbit(17,3,0.4800217244625303E-1_8,0.0_8, &
                             0.1605023107954450E-3_8)
      CALL load_laikov_orbit(17,3,0.8244922058397242E-1_8,0.0_8, &
                             0.2072200131464099E-3_8)
      CALL load_laikov_orbit(17,3,0.1200408362484023E+0_8,0.0_8, &
                             0.2431297618814187E-3_8)
      CALL load_laikov_orbit(17,3,0.1595773530809965E+0_8,0.0_8, &
                             0.2711819064496707E-3_8)
      CALL load_laikov_orbit(17,3,0.2002635973434064E+0_8,0.0_8, &
                             0.2932762038321116E-3_8)
      CALL load_laikov_orbit(17,3,0.2415127590139982E+0_8,0.0_8, &
                             0.3107032514197368E-3_8)
      CALL load_laikov_orbit(17,3,0.2828584158458477E+0_8,0.0_8, &
                             0.3243808058921213E-3_8)
      CALL load_laikov_orbit(17,3,0.3239091015338138E+0_8,0.0_8, &
                             0.3349899091374030E-3_8)
      CALL load_laikov_orbit(17,3,0.3643225097962194E+0_8,0.0_8, &
                             0.3430580688505218E-3_8)
      CALL load_laikov_orbit(17,3,0.4037897083691802E+0_8,0.0_8, &
                             0.3490124109290343E-3_8)
      CALL load_laikov_orbit(17,3,0.4420247515194127E+0_8,0.0_8, &
                             0.3532148948561955E-3_8)
      CALL load_laikov_orbit(17,3,0.4787572538464938E+0_8,0.0_8, &
                             0.3559862669062833E-3_8)
      CALL load_laikov_orbit(17,3,0.5137265251275234E+0_8,0.0_8, &
                             0.3576224317551411E-3_8)
      CALL load_laikov_orbit(17,3,0.5466764056654611E+0_8,0.0_8, &
                             0.3584050533086076E-3_8)
      CALL load_laikov_orbit(17,3,0.6054859420813535E+0_8,0.0_8, &
                             0.3584903581373224E-3_8)
      CALL load_laikov_orbit(17,3,0.6308106701764562E+0_8,0.0_8, &
                             0.3582991879040586E-3_8)
      CALL load_laikov_orbit(17,3,0.6530369230179584E+0_8,0.0_8, &
                             0.3582371187963125E-3_8)
      CALL load_laikov_orbit(17,3,0.6718609524611158E+0_8,0.0_8, &
                             0.3584353631122350E-3_8)
      CALL load_laikov_orbit(17,3,0.6869676499894013E+0_8,0.0_8, &
                             0.3589120166517785E-3_8)
      CALL load_laikov_orbit(17,3,0.6980467077240748E+0_8,0.0_8, &
                             0.3595445704531601E-3_8)
      CALL load_laikov_orbit(17,3,0.7048241721250522E+0_8,0.0_8, &
                             0.3600943557111074E-3_8)
      CALL load_laikov_orbit(17,4,0.5591105222058232E-1_8,0.0_8, &
                             0.1456447096742039E-3_8)
      CALL load_laikov_orbit(17,4,0.1407384078513916E+0_8,0.0_8, &
                             0.2252370188283782E-3_8)
      CALL load_laikov_orbit(17,4,0.2364035438976309E+0_8,0.0_8, &
                             0.2766135443474897E-3_8)
      CALL load_laikov_orbit(17,4,0.3360602737818170E+0_8,0.0_8, &
                             0.3110729491500851E-3_8)
      CALL load_laikov_orbit(17,4,0.4356292630054665E+0_8,0.0_8, &
                             0.3342506712303391E-3_8)
      CALL load_laikov_orbit(17,4,0.5321569415256174E+0_8,0.0_8, &
                             0.3491981834026860E-3_8)
      CALL load_laikov_orbit(17,4,0.6232956305040554E+0_8,0.0_8, &
                             0.3576003604348932E-3_8)
      CALL load_laikov_orbit(17,5,0.9469870086838469E-1_8,0.2778748387309470E-1_8, &
                             0.1921921305788564E-3_8)
      CALL load_laikov_orbit(17,5,0.1353170300568141E+0_8,0.6076569878628364E-1_8, &
                             0.2301458216495632E-3_8)
      CALL load_laikov_orbit(17,5,0.1771679481726077E+0_8,0.9703072762711040E-1_8, &
                             0.2604248549522893E-3_8)
      CALL load_laikov_orbit(17,5,0.2197066664231751E+0_8,0.1354112458524762E+0_8, &
                             0.2845275425870697E-3_8)
      CALL load_laikov_orbit(17,5,0.2624783557374927E+0_8,0.1750996479744100E+0_8, &
                             0.3036870897974840E-3_8)
      CALL load_laikov_orbit(17,5,0.3050969521214442E+0_8,0.2154896907449802E+0_8, &
                             0.3188414832298066E-3_8)
      CALL load_laikov_orbit(17,5,0.3472252637196021E+0_8,0.2560954625740152E+0_8, &
                             0.3307046414722089E-3_8)
      CALL load_laikov_orbit(17,5,0.3885610219026360E+0_8,0.2965070050624096E+0_8, &
                             0.3398330969031360E-3_8)
      CALL load_laikov_orbit(17,5,0.4288273776062765E+0_8,0.3363641488734497E+0_8, &
                             0.3466757899705373E-3_8)
      CALL load_laikov_orbit(17,5,0.4677662471302948E+0_8,0.3753400029836788E+0_8, &
                             0.3516095923230054E-3_8)
      CALL load_laikov_orbit(17,5,0.5051333589553359E+0_8,0.4131297522144286E+0_8, &
                             0.3549645184048486E-3_8)
      CALL load_laikov_orbit(17,5,0.5406942145810492E+0_8,0.4494423776081795E+0_8, &
                             0.3570415969441392E-3_8)
      CALL load_laikov_orbit(17,5,0.5742204122576457E+0_8,0.4839938958841502E+0_8, &
                             0.3581251798496118E-3_8)
      CALL load_laikov_orbit(17,5,0.1865407027225188E+0_8,0.3259144851070796E-1_8, &
                             0.2543491329913348E-3_8)
      CALL load_laikov_orbit(17,5,0.2321186453689432E+0_8,0.6835679505297343E-1_8, &
                             0.2786711051330776E-3_8)
      CALL load_laikov_orbit(17,5,0.2773159142523882E+0_8,0.1062284864451989E+0_8, &
                             0.2985552361083679E-3_8)
      CALL load_laikov_orbit(17,5,0.3219200192237254E+0_8,0.1454404409323047E+0_8, &
                             0.3145867929154039E-3_8)
      CALL load_laikov_orbit(17,5,0.3657032593944029E+0_8,0.1854018282582510E+0_8, &
                             0.3273290662067609E-3_8)
      CALL load_laikov_orbit(17,5,0.4084376778363622E+0_8,0.2256297412014750E+0_8, &
                             0.3372705511943501E-3_8)
      CALL load_laikov_orbit(17,5,0.4499004945751427E+0_8,0.2657104425000896E+0_8, &
                             0.3448274437851510E-3_8)
      CALL load_laikov_orbit(17,5,0.4898758141326335E+0_8,0.3052755487631557E+0_8, &
                             0.3503592783048583E-3_8)
      CALL load_laikov_orbit(17,5,0.5281547442266309E+0_8,0.3439863920645423E+0_8, &
                             0.3541854792663162E-3_8)
      CALL load_laikov_orbit(17,5,0.5645346989813992E+0_8,0.3815229456121914E+0_8, &
                             0.3565995517909428E-3_8)
      CALL load_laikov_orbit(17,5,0.5988181252159848E+0_8,0.4175752420966734E+0_8, &
                             0.3578802078302898E-3_8)
      CALL load_laikov_orbit(17,5,0.2850425424471603E+0_8,0.3562149509862536E-1_8, &
                             0.2958644592860982E-3_8)
      CALL load_laikov_orbit(17,5,0.3324619433027876E+0_8,0.7330318886871096E-1_8, &
                             0.3119548129116835E-3_8)
      CALL load_laikov_orbit(17,5,0.3785848333076282E+0_8,0.1123226296008472E+0_8, &
                             0.3250745225005984E-3_8)
      CALL load_laikov_orbit(17,5,0.4232891028562115E+0_8,0.1521084193337708E+0_8, &
                             0.3355153415935208E-3_8)
      CALL load_laikov_orbit(17,5,0.4664287050829722E+0_8,0.1921844459223610E+0_8, &
                             0.3435847568549328E-3_8)
      CALL load_laikov_orbit(17,5,0.5078458493735726E+0_8,0.2321360989678303E+0_8, &
                             0.3495786831622488E-3_8)
      CALL load_laikov_orbit(17,5,0.5473779816204180E+0_8,0.2715886486360520E+0_8, &
                             0.3537767805534621E-3_8)
      CALL load_laikov_orbit(17,5,0.5848617133811376E+0_8,0.3101924707571355E+0_8, &
                             0.3564459815421428E-3_8)
      CALL load_laikov_orbit(17,5,0.6201348281584888E+0_8,0.3476121052890973E+0_8, &
                             0.3578464061225468E-3_8)
      CALL load_laikov_orbit(17,5,0.3852191185387871E+0_8,0.3763224880035108E-1_8, &
                             0.3239748762836212E-3_8)
      CALL load_laikov_orbit(17,5,0.4325025061073423E+0_8,0.7659581935637135E-1_8, &
                             0.3345491784174287E-3_8)
      CALL load_laikov_orbit(17,5,0.4778486229734490E+0_8,0.1163381306083900E+0_8, &
                             0.3429126177301782E-3_8)
      CALL load_laikov_orbit(17,5,0.5211663693009000E+0_8,0.1563890598752899E+0_8, &
                             0.3492420343097421E-3_8)
      CALL load_laikov_orbit(17,5,0.5623469504853703E+0_8,0.1963320810149200E+0_8, &
                             0.3537399050235257E-3_8)
      CALL load_laikov_orbit(17,5,0.6012718188659246E+0_8,0.2357847407258738E+0_8, &
                             0.3566209152659172E-3_8)
      CALL load_laikov_orbit(17,5,0.6378179206390117E+0_8,0.2743846121244060E+0_8, &
                             0.3581084321919782E-3_8)
      CALL load_laikov_orbit(17,5,0.4836936460214534E+0_8,0.3895902610739024E-1_8, &
                             0.3426522117591512E-3_8)
      CALL load_laikov_orbit(17,5,0.5293792562683797E+0_8,0.7871246819312640E-1_8, &
                             0.3491848770121379E-3_8)
      CALL load_laikov_orbit(17,5,0.5726281253100033E+0_8,0.1187963808202981E+0_8, &
                             0.3539318235231476E-3_8)
      CALL load_laikov_orbit(17,5,0.6133658776169068E+0_8,0.1587914708061787E+0_8, &
                             0.3570231438458694E-3_8)
      CALL load_laikov_orbit(17,5,0.6515085491865307E+0_8,0.1983058575227646E+0_8, &
                             0.3586207335051714E-3_8)
      CALL load_laikov_orbit(17,5,0.5778692716064976E+0_8,0.3977209689791542E-1_8, &
                             0.3541196205164025E-3_8)
      CALL load_laikov_orbit(17,5,0.6207904288086192E+0_8,0.7990157592981152E-1_8, &
                             0.3574296911573953E-3_8)
      CALL load_laikov_orbit(17,5,0.6608688171046802E+0_8,0.1199671308754309E+0_8, &
                             0.3591993279818963E-3_8)
      CALL load_laikov_orbit(17,5,0.6656263089489130E+0_8,0.4015955957805969E-1_8, &
                             0.3595855034661997E-3_8)
      IF (nlgp /= 3074) CALL paw_lebedev_abort("Invalid Lebedev 3074-point table")

      ! l = 131, 5810 points.
      CALL load_laikov_orbit(18,0,0.0_8,0.0_8, &
                             0.9735347946175486E-5_8)
      CALL load_laikov_orbit(18,1,0.0_8,0.0_8, &
                             0.1907581241803167E-3_8)
      CALL load_laikov_orbit(18,2,0.0_8,0.0_8, &
                             0.1901059546737578E-3_8)
      CALL load_laikov_orbit(18,3,0.1182361662400277E-1_8,0.0_8, &
                             0.3926424538919212E-4_8)
      CALL load_laikov_orbit(18,3,0.3062145009138958E-1_8,0.0_8, &
                             0.6667905467294382E-4_8)
      CALL load_laikov_orbit(18,3,0.5329794036834243E-1_8,0.0_8, &
                             0.8868891315019135E-4_8)
      CALL load_laikov_orbit(18,3,0.7848165532862220E-1_8,0.0_8, &
                             0.1066306000958872E-3_8)
      CALL load_laikov_orbit(18,3,0.1054038157636201E+0_8,0.0_8, &
                             0.1214506743336128E-3_8)
      CALL load_laikov_orbit(18,3,0.1335577797766211E+0_8,0.0_8, &
                             0.1338054681640871E-3_8)
      CALL load_laikov_orbit(18,3,0.1625769955502252E+0_8,0.0_8, &
                             0.1441677023628504E-3_8)
      CALL load_laikov_orbit(18,3,0.1921787193412792E+0_8,0.0_8, &
                             0.1528880200826557E-3_8)
      CALL load_laikov_orbit(18,3,0.2221340534690548E+0_8,0.0_8, &
                             0.1602330623773609E-3_8)
      CALL load_laikov_orbit(18,3,0.2522504912791132E+0_8,0.0_8, &
                             0.1664102653445244E-3_8)
      CALL load_laikov_orbit(18,3,0.2823610860679697E+0_8,0.0_8, &
                             0.1715845854011323E-3_8)
      CALL load_laikov_orbit(18,3,0.3123173966267560E+0_8,0.0_8, &
                             0.1758901000133069E-3_8)
      CALL load_laikov_orbit(18,3,0.3419847036953789E+0_8,0.0_8, &
                             0.1794382485256736E-3_8)
      CALL load_laikov_orbit(18,3,0.3712386456999758E+0_8,0.0_8, &
                             0.1823238106757407E-3_8)
      CALL load_laikov_orbit(18,3,0.3999627649876828E+0_8,0.0_8, &
                             0.1846293252959976E-3_8)
      CALL load_laikov_orbit(18,3,0.4280466458648093E+0_8,0.0_8, &
                             0.1864284079323098E-3_8)
      CALL load_laikov_orbit(18,3,0.4553844360185711E+0_8,0.0_8, &
                             0.1877882694626914E-3_8)
      CALL load_laikov_orbit(18,3,0.4818736094437834E+0_8,0.0_8, &
                             0.1887716321852025E-3_8)
      CALL load_laikov_orbit(18,3,0.5074138709260629E+0_8,0.0_8, &
                             0.1894381638175673E-3_8)
      CALL load_laikov_orbit(18,3,0.5319061304570707E+0_8,0.0_8, &
                             0.1898454899533629E-3_8)
      CALL load_laikov_orbit(18,3,0.5552514978677286E+0_8,0.0_8, &
                             0.1900497929577815E-3_8)
      CALL load_laikov_orbit(18,3,0.5981009025246183E+0_8,0.0_8, &
                             0.1900671501924092E-3_8)
      CALL load_laikov_orbit(18,3,0.6173990192228116E+0_8,0.0_8, &
                             0.1899837555533510E-3_8)
      CALL load_laikov_orbit(18,3,0.6351365239411131E+0_8,0.0_8, &
                             0.1899014113156229E-3_8)
      CALL load_laikov_orbit(18,3,0.6512010228227200E+0_8,0.0_8, &
                             0.1898581257705106E-3_8)
      CALL load_laikov_orbit(18,3,0.6654758363948120E+0_8,0.0_8, &
                             0.1898804756095753E-3_8)
      CALL load_laikov_orbit(18,3,0.6778410414853370E+0_8,0.0_8, &
                             0.1899793610426402E-3_8)
      CALL load_laikov_orbit(18,3,0.6881760887484110E+0_8,0.0_8, &
                             0.1901464554844117E-3_8)
      CALL load_laikov_orbit(18,3,0.6963645267094598E+0_8,0.0_8, &
                             0.1903533246259542E-3_8)
      CALL load_laikov_orbit(18,3,0.7023010617153579E+0_8,0.0_8, &
                             0.1905556158463228E-3_8)
      CALL load_laikov_orbit(18,3,0.7059004636628753E+0_8,0.0_8, &
                             0.1907037155663528E-3_8)
      CALL load_laikov_orbit(18,4,0.3552470312472575E-1_8,0.0_8, &
                             0.5992997844249967E-4_8)
      CALL load_laikov_orbit(18,4,0.9151176620841283E-1_8,0.0_8, &
                             0.9749059382456978E-4_8)
      CALL load_laikov_orbit(18,4,0.1566197930068980E+0_8,0.0_8, &
                             0.1241680804599158E-3_8)
      CALL load_laikov_orbit(18,4,0.2265467599271907E+0_8,0.0_8, &
                             0.1437626154299360E-3_8)
      CALL load_laikov_orbit(18,4,0.2988242318581361E+0_8,0.0_8, &
                             0.1584200054793902E-3_8)
      CALL load_laikov_orbit(18,4,0.3717482419703886E+0_8,0.0_8, &
                             0.1694436550982744E-3_8)
      CALL load_laikov_orbit(18,4,0.4440094491758889E+0_8,0.0_8, &
                             0.1776617014018108E-3_8)
      CALL load_laikov_orbit(18,4,0.5145337096756642E+0_8,0.0_8, &
                             0.1836132434440077E-3_8)
      CALL load_laikov_orbit(18,4,0.5824053672860230E+0_8,0.0_8, &
                             0.1876494727075983E-3_8)
      CALL load_laikov_orbit(18,4,0.6468283961043370E+0_8,0.0_8, &
                             0.1899906535336482E-3_8)
      CALL load_laikov_orbit(18,5,0.6095964259104373E-1_8,0.1787828275342931E-1_8, &
                             0.8143252820767350E-4_8)
      CALL load_laikov_orbit(18,5,0.8811962270959388E-1_8,0.3953888740792096E-1_8, &
                             0.9998859890887728E-4_8)
      CALL load_laikov_orbit(18,5,0.1165936722428831E+0_8,0.6378121797722990E-1_8, &
                             0.1156199403068359E-3_8)
      CALL load_laikov_orbit(18,5,0.1460232857031785E+0_8,0.8985890813745037E-1_8, &
                             0.1287632092635513E-3_8)
      CALL load_laikov_orbit(18,5,0.1761197110181755E+0_8,0.1172606510576162E+0_8, &
                             0.1398378643365139E-3_8)
      CALL load_laikov_orbit(18,5,0.2066471190463718E+0_8,0.1456102876970995E+0_8, &
                             0.1491876468417391E-3_8)
      CALL load_laikov_orbit(18,5,0.2374076026328152E+0_8,0.1746153823011775E+0_8, &
                             0.1570855679175456E-3_8)
      CALL load_laikov_orbit(18,5,0.2682305474337051E+0_8,0.2040383070295584E+0_8, &
                             0.1637483948103775E-3_8)
      CALL load_laikov_orbit(18,5,0.2989653312142369E+0_8,0.2336788634003698E+0_8, &
                             0.1693500566632843E-3_8)
      CALL load_laikov_orbit(18,5,0.3294762752772209E+0_8,0.2633632752654219E+0_8, &
                             0.1740322769393633E-3_8)
      CALL load_laikov_orbit(18,5,0.3596390887276086E+0_8,0.2929369098051601E+0_8, &
                             0.1779126637278296E-3_8)
      CALL load_laikov_orbit(18,5,0.3893383046398812E+0_8,0.3222592785275512E+0_8, &
                             0.1810908108835412E-3_8)
      CALL load_laikov_orbit(18,5,0.4184653789358347E+0_8,0.3512004791195743E+0_8, &
                             0.1836529132600190E-3_8)
      CALL load_laikov_orbit(18,5,0.4469172319076166E+0_8,0.3796385677684537E+0_8, &
                             0.1856752841777379E-3_8)
      CALL load_laikov_orbit(18,5,0.4745950813276976E+0_8,0.4074575378263879E+0_8, &
                             0.1872270566606832E-3_8)
      CALL load_laikov_orbit(18,5,0.5014034601410262E+0_8,0.4345456906027828E+0_8, &
                             0.1883722645591307E-3_8)
      CALL load_laikov_orbit(18,5,0.5272493404551239E+0_8,0.4607942515205134E+0_8, &
                             0.1891714324525297E-3_8)
      CALL load_laikov_orbit(18,5,0.5520413051846366E+0_8,0.4860961284181720E+0_8, &
                             0.1896827480450146E-3_8)
      CALL load_laikov_orbit(18,5,0.5756887237503077E+0_8,0.5103447395342790E+0_8, &
                             0.1899628417059528E-3_8)
      CALL load_laikov_orbit(18,5,0.1225039430588352E+0_8,0.2136455922655793E-1_8, &
                             0.1123301829001669E-3_8)
      CALL load_laikov_orbit(18,5,0.1539113217321372E+0_8,0.4520926166137188E-1_8, &
                             0.1253698826711277E-3_8)
      CALL load_laikov_orbit(18,5,0.1856213098637712E+0_8,0.7086468177864818E-1_8, &
                             0.1366266117678531E-3_8)
      CALL load_laikov_orbit(18,5,0.2174998728035131E+0_8,0.9785239488772918E-1_8, &
                             0.1462736856106918E-3_8)
      CALL load_laikov_orbit(18,5,0.2494128336938330E+0_8,0.1258106396267210E+0_8, &
                             0.1545076466685412E-3_8)
      CALL load_laikov_orbit(18,5,0.2812321562143480E+0_8,0.1544529125047001E+0_8, &
                             0.1615096280814007E-3_8)
      CALL load_laikov_orbit(18,5,0.3128372276456111E+0_8,0.1835433512202753E+0_8, &
                             0.1674366639741759E-3_8)
      CALL load_laikov_orbit(18,5,0.3441145160177973E+0_8,0.2128813258619585E+0_8, &
                             0.1724225002437900E-3_8)
      CALL load_laikov_orbit(18,5,0.3749567714853510E+0_8,0.2422913734880829E+0_8, &
                             0.1765810822987288E-3_8)
      CALL load_laikov_orbit(18,5,0.4052621732015610E+0_8,0.2716163748391453E+0_8, &
                             0.1800104126010751E-3_8)
      CALL load_laikov_orbit(18,5,0.4349335453522385E+0_8,0.3007127671240280E+0_8, &
                             0.1827960437331284E-3_8)
      CALL load_laikov_orbit(18,5,0.4638776641524965E+0_8,0.3294470677216479E+0_8, &
                             0.1850140300716308E-3_8)
      CALL load_laikov_orbit(18,5,0.4920046410462687E+0_8,0.3576932543699155E+0_8, &
                             0.1867333507394938E-3_8)
      CALL load_laikov_orbit(18,5,0.5192273554861704E+0_8,0.3853307059757764E+0_8, &
                             0.1880178688638289E-3_8)
      CALL load_laikov_orbit(18,5,0.5454609081136522E+0_8,0.4122425044452694E+0_8, &
                             0.1889278925654758E-3_8)
      CALL load_laikov_orbit(18,5,0.5706220661424140E+0_8,0.4383139587781027E+0_8, &
                             0.1895213832507346E-3_8)
      CALL load_laikov_orbit(18,5,0.5946286755181518E+0_8,0.4634312536300553E+0_8, &
                             0.1898548277397420E-3_8)
      CALL load_laikov_orbit(18,5,0.1905370790924295E+0_8,0.2371311537781979E-1_8, &
                             0.1349105935937341E-3_8)
      CALL load_laikov_orbit(18,5,0.2242518717748009E+0_8,0.4917878059254806E-1_8, &
                             0.1444060068369326E-3_8)
      CALL load_laikov_orbit(18,5,0.2577190808025936E+0_8,0.7595498960495142E-1_8, &
                             0.1526797390930008E-3_8)
      CALL load_laikov_orbit(18,5,0.2908724534927187E+0_8,0.1036991083191100E+0_8, &
                             0.1598208771406474E-3_8)
      CALL load_laikov_orbit(18,5,0.3236354020056219E+0_8,0.1321348584450234E+0_8, &
                             0.1659354368615331E-3_8)
      CALL load_laikov_orbit(18,5,0.3559267359304543E+0_8,0.1610316571314789E+0_8, &
                             0.1711279910946440E-3_8)
      CALL load_laikov_orbit(18,5,0.3876637123676956E+0_8,0.1901912080395707E+0_8, &
                             0.1754952725601440E-3_8)
      CALL load_laikov_orbit(18,5,0.4187636705218842E+0_8,0.2194384950137950E+0_8, &
                             0.1791247850802529E-3_8)
      CALL load_laikov_orbit(18,5,0.4491449019883107E+0_8,0.2486155334763858E+0_8, &
                             0.1820954300877716E-3_8)
      CALL load_laikov_orbit(18,5,0.4787270932425445E+0_8,0.2775768931812335E+0_8, &
                             0.1844788524548449E-3_8)
      CALL load_laikov_orbit(18,5,0.5074315153055574E+0_8,0.3061863786591120E+0_8, &
                             0.1863409481706220E-3_8)
      CALL load_laikov_orbit(18,5,0.5351810507738336E+0_8,0.3343144718152556E+0_8, &
                             0.1877433008795068E-3_8)
      CALL load_laikov_orbit(18,5,0.5619001025975381E+0_8,0.3618362729028427E+0_8, &
                             0.1887444543705232E-3_8)
      CALL load_laikov_orbit(18,5,0.5875144035268046E+0_8,0.3886297583620408E+0_8, &
                             0.1894009829375006E-3_8)
      CALL load_laikov_orbit(18,5,0.6119507308734495E+0_8,0.4145742277792031E+0_8, &
                             0.1897683345035198E-3_8)
      CALL load_laikov_orbit(18,5,0.2619733870119463E+0_8,0.2540047186389353E-1_8, &
                             0.1517327037467653E-3_8)
      CALL load_laikov_orbit(18,5,0.2968149743237949E+0_8,0.5208107018543989E-1_8, &
                             0.1587740557483543E-3_8)
      CALL load_laikov_orbit(18,5,0.3310451504860488E+0_8,0.7971828470885599E-1_8, &
                             0.1649093382274097E-3_8)
      CALL load_laikov_orbit(18,5,0.3646215567376676E+0_8,0.1080465999177927E+0_8, &
                             0.1701915216193265E-3_8)
      CALL load_laikov_orbit(18,5,0.3974916785279360E+0_8,0.1368413849366629E+0_8, &
                             0.1746847753144065E-3_8)
      CALL load_laikov_orbit(18,5,0.4295967403772029E+0_8,0.1659073184763559E+0_8, &
                             0.1784555512007570E-3_8)
      CALL load_laikov_orbit(18,5,0.4608742854473447E+0_8,0.1950703730454614E+0_8, &
                             0.1815687562112174E-3_8)
      CALL load_laikov_orbit(18,5,0.4912598858949903E+0_8,0.2241721144376724E+0_8, &
                             0.1840864370663302E-3_8)
      CALL load_laikov_orbit(18,5,0.5206882758945558E+0_8,0.2530655255406489E+0_8, &
                             0.1860676785390006E-3_8)
      CALL load_laikov_orbit(18,5,0.5490940914019819E+0_8,0.2816118409731066E+0_8, &
                             0.1875690583743703E-3_8)
      CALL load_laikov_orbit(18,5,0.5764123302025542E+0_8,0.3096780504593238E+0_8, &
                             0.1886453236347225E-3_8)
      CALL load_laikov_orbit(18,5,0.6025786004213506E+0_8,0.3371348366394987E+0_8, &
                             0.1893501123329645E-3_8)
      CALL load_laikov_orbit(18,5,0.6275291964794956E+0_8,0.3638547827694396E+0_8, &
                             0.1897366184519868E-3_8)
      CALL load_laikov_orbit(18,5,0.3348189479861771E+0_8,0.2664841935537443E-1_8, &
                             0.1643908815152736E-3_8)
      CALL load_laikov_orbit(18,5,0.3699515545855295E+0_8,0.5424000066843495E-1_8, &
                             0.1696300350907768E-3_8)
      CALL load_laikov_orbit(18,5,0.4042003071474669E+0_8,0.8251992715430854E-1_8, &
                             0.1741553103844483E-3_8)
      CALL load_laikov_orbit(18,5,0.4375320100182624E+0_8,0.1112695182483710E+0_8, &
                             0.1780015282386092E-3_8)
      CALL load_laikov_orbit(18,5,0.4699054490335947E+0_8,0.1402964116467816E+0_8, &
                             0.1812116787077125E-3_8)
      CALL load_laikov_orbit(18,5,0.5012739879431952E+0_8,0.1694275117584291E+0_8, &
                             0.1838323158085421E-3_8)
      CALL load_laikov_orbit(18,5,0.5315874883754966E+0_8,0.1985038235312689E+0_8, &
                             0.1859113119837737E-3_8)
      CALL load_laikov_orbit(18,5,0.5607937109622117E+0_8,0.2273765660020893E+0_8, &
                             0.1874969220221698E-3_8)
      CALL load_laikov_orbit(18,5,0.5888393223495521E+0_8,0.2559041492849764E+0_8, &
                             0.1886375612681076E-3_8)
      CALL load_laikov_orbit(18,5,0.6156705979160163E+0_8,0.2839497251976899E+0_8, &
                             0.1893819575809276E-3_8)
      CALL load_laikov_orbit(18,5,0.6412338809078123E+0_8,0.3113791060500690E+0_8, &
                             0.1897794748256767E-3_8)
      CALL load_laikov_orbit(18,5,0.4076051259257167E+0_8,0.2757792290858463E-1_8, &
                             0.1738963926584846E-3_8)
      CALL load_laikov_orbit(18,5,0.4423788125791520E+0_8,0.5584136834984293E-1_8, &
                             0.1777442359873466E-3_8)
      CALL load_laikov_orbit(18,5,0.4760480917328258E+0_8,0.8457772087727143E-1_8, &
                             0.1810010815068719E-3_8)
      CALL load_laikov_orbit(18,5,0.5085838725946297E+0_8,0.1135975846359248E+0_8, &
                             0.1836920318248129E-3_8)
      CALL load_laikov_orbit(18,5,0.5399513637391218E+0_8,0.1427286904765053E+0_8, &
                             0.1858489473214328E-3_8)
      CALL load_laikov_orbit(18,5,0.5701118433636380E+0_8,0.1718112740057635E+0_8, &
                             0.1875079342496592E-3_8)
      CALL load_laikov_orbit(18,5,0.5990240530606021E+0_8,0.2006944855985351E+0_8, &
                             0.1887080239102310E-3_8)
      CALL load_laikov_orbit(18,5,0.6266452685139695E+0_8,0.2292335090598907E+0_8, &
                             0.1894905752176822E-3_8)
      CALL load_laikov_orbit(18,5,0.6529320971415942E+0_8,0.2572871512353714E+0_8, &
                             0.1898991061200695E-3_8)
      CALL load_laikov_orbit(18,5,0.4791583834610126E+0_8,0.2826094197735932E-1_8, &
                             0.1809065016458791E-3_8)
      CALL load_laikov_orbit(18,5,0.5130373952796940E+0_8,0.5699871359683649E-1_8, &
                             0.1836297121596799E-3_8)
      CALL load_laikov_orbit(18,5,0.5456252429628476E+0_8,0.8602712528554394E-1_8, &
                             0.1858426916241869E-3_8)
      CALL load_laikov_orbit(18,5,0.5768956329682385E+0_8,0.1151748137221281E+0_8, &
                             0.1875654101134641E-3_8)
      CALL load_laikov_orbit(18,5,0.6068186944699046E+0_8,0.1442811654136362E+0_8, &
                             0.1888240751833503E-3_8)
      CALL load_laikov_orbit(18,5,0.6353622248024907E+0_8,0.1731930321657680E+0_8, &
                             0.1896497383866979E-3_8)
      CALL load_laikov_orbit(18,5,0.6624927035731797E+0_8,0.2017619958756061E+0_8, &
                             0.1900775530219121E-3_8)
      CALL load_laikov_orbit(18,5,0.5484933508028488E+0_8,0.2874219755907391E-1_8, &
                             0.1858525041478814E-3_8)
      CALL load_laikov_orbit(18,5,0.5810207682142106E+0_8,0.5778312123713695E-1_8, &
                             0.1876248690077947E-3_8)
      CALL load_laikov_orbit(18,5,0.6120955197181352E+0_8,0.8695262371439526E-1_8, &
                             0.1889404439064607E-3_8)
      CALL load_laikov_orbit(18,5,0.6416944284294319E+0_8,0.1160893767057166E+0_8, &
                             0.1898168539265290E-3_8)
      CALL load_laikov_orbit(18,5,0.6697926391731260E+0_8,0.1450378826743251E+0_8, &
                             0.1902779940661772E-3_8)
      CALL load_laikov_orbit(18,5,0.6147594390585488E+0_8,0.2904957622341456E-1_8, &
                             0.1890125641731815E-3_8)
      CALL load_laikov_orbit(18,5,0.6455390026356783E+0_8,0.5823809152617197E-1_8, &
                             0.1899434637795751E-3_8)
      CALL load_laikov_orbit(18,5,0.6747258588365477E+0_8,0.8740384899884715E-1_8, &
                             0.1904520856831751E-3_8)
      CALL load_laikov_orbit(18,5,0.6772135750395347E+0_8,0.2919946135808105E-1_8, &
                             0.1905534498734563E-3_8)
      IF (nlgp /= 5810) CALL paw_lebedev_abort("Invalid Lebedev 5810-point table")

      init_lebedev_grids_done = .TRUE.

   END SUBROUTINE init_lebedev_grids

END MODULE PAW_LEBEDEV_MODULE
