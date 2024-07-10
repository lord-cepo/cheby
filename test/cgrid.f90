MODULE test_cgrid
   USE testdrive,     ONLY : new_unittest, unittest_type, error_type, check
   USE q_grids,       ONLY : q_grid, setup_simple_grid
   USE coeff_grid,    ONLY : cgrid
   USE kinds,         ONLY : DP
   USE ph_system,     ONLY : ph_system_info
   USE input_fc,      ONLY : forceconst2_grid, read_fc2, aux_system, div_mass_fc2
   USE asr2_module,   ONLY : impose_asr2

   implicit none
   private
   TYPE(cgrid) :: unit_cgrid
   TYPE(ph_system_info) :: unit_S

   TYPE(ph_system_info) :: Si_S
   TYPE(forceconst2_grid) :: Si_fc2
   TYPE(cgrid) :: Si_cgrid
   REAL(DP), PARAMETER :: THRESHOLD = 1d-8
   public :: collect_cgrid

CONTAINS
   !> Collect all exported unit tests
   subroutine collect_cgrid(testsuite)
      !> Collection of tests
      type(unittest_type), allocatable, intent(out) :: testsuite(:)
      CALL setup_test()

      testsuite = [ &
         new_unittest("testing cgrid_coeffs", test_cgrid_coeffs), &
         new_unittest("testing cgrid_coeffs1", test_cgrid_coeffs1), &
         new_unittest("setting Si example", set_Si_example), &
         new_unittest("testing cgrid%coefficients", test_coeff_value) &
         ]
   end subroutine collect_cgrid

   FUNCTION unit_matrix() result(M)
      REAL(DP) :: M(3,3)
      M = 0._dp
      M(1,1) = 1._dp
      M(2,2) = 1._dp
      M(3,3) = 1._dp
   END FUNCTION unit_matrix

   !> sets a q_grid with 8 points: 000, 001, ...
   SUBROUTINE setup_test()
      TYPE(cgrid) :: grid

      ! local variables
      TYPE(q_grid) :: qgrid
      TYPE(ph_system_info) :: S
      REAL(DP) :: bg(3,3)
      INTEGER :: iq

      bg = unit_matrix()

      CALL setup_simple_grid(bg, 3,5,1, qgrid)
      grid%qgrid = qgrid
      ALLOCATE(grid%coefficients(1,1,qgrid%nqtot))
      DO iq = 1, qgrid%nqtot
         grid%coefficients(1,1,iq) = CMPLX(REAL(iq, DP), 0._dp)
      END DO
      grid%nR = 1
      unit_cgrid = grid

      !> sets ph_system_info
      S%at = unit_matrix()
      S%bg = unit_matrix()
      S%nat3 = 1
      unit_S = S

   END SUBROUTINE setup_test

   SUBROUTINE test_cgrid_coeffs(error)

      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      ! local variables
      REAL(DP) :: xq_test(3), cgrid_coeffs(1,1), delta_x(3)
      INTEGER :: iq
      TYPE(q_grid) :: g
      g = unit_cgrid%qgrid
      !> the grid is [0.0 0.0 0.0], [0.0 0.0 0.5], ...
      !> the coeffs are          1,             2, ...
      !> We want points in [0.25,0.75] to be identified in the 0.5 point

      DO iq = 1, g%nq
         xq_test = g%xq(:,iq)
         delta_x = (/MOD(iq, 2), -MOD(iq,3), MOD(iq,4)/)
         cgrid_coeffs = unit_cgrid%coeffs(xq_test + delta_x, unit_S)
         CALL check(error, REAL(cgrid_coeffs(1,1), DP), REAL(unit_cgrid%coefficients(1,1,iq), DP))
         if (allocated(error)) return
      END DO

      ! DO iq = 1, g%nq
      !    SGN = -1
      !    IF(MOD(iq, 2) == 1) SGN = 1
      !    xq_test = g%xq(:,iq)
      !    delta_x = (/0._dp, 0._dp, SGN*1._dp/g%nq/)
      !    cgrid_coeffs = unit_cgrid%coeffs(xq_test + delta_x, unit_S)
      !    CALL check(error, cgrid_coeffs(1,1), &
      !     unit_cgrid%coefficients(1,1,iq) * (1 - 1._dp/g%nq) + unit_cgrid%coefficients(1,1,iq + SGN)/g%nq)
      !    IF (allocated(error)) THEN
      !       WRITE(*,*) "xq", g%xq(:,iq)
      !       WRITE(*,*) "iq", iq, "SGN", SGN
      !       ! WRITE(*,*)
      !    END IF
      ! END DO

   END SUBROUTINE test_cgrid_coeffs

   SUBROUTINE test_cgrid_coeffs1(error)

      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      ! local variables
      TYPE(q_grid) g
      REAL(DP) :: xq_test(3), cgrid_coeffs(1,1), delta_x(3)
      INTEGER :: iq, SGN, iq_
      g = unit_cgrid%qgrid


      !> the grid is [0.0 0.0 0.0], [0.0 0.0 0.5], ...
      !> the coeffs are          1,             2, ...
      !> we expect a linear interpolation in each direction

      ! xq_test = (/ 0.0, 0.0, 0.0/)
      ! cgrid_coeffs = unit_cgrid%coeffs1(xq_test, unit_S)
      ! CALL check(error, cgrid_coeffs(1,1), 1._dp)

      ! xq_test = (/ 0.0, 0.0, 0.25/)
      ! cgrid_coeffs = unit_cgrid%coeffs1(xq_test, unit_S)
      ! CALL check(error, cgrid_coeffs(1,1), 1.5_dp)

      ! xq_test = (/ 0.0, 0.5, 0.0/)
      ! cgrid_coeffs = unit_cgrid%coeffs1(xq_test, unit_S)
      ! CALL check(error, cgrid_coeffs(1,1), 3._dp)

      ! xq_test = (/ -4.1, 0.0, 0.0/)
      ! cgrid_coeffs = unit_cgrid%coeffs1(xq_test, unit_S)
      ! CALL check(error, ABS(cgrid_coeffs(1,1) - 1.8_dp) < 1d-4)

      ! xq_test = (/ 0.75, 0.0, 0.0/)
      ! cgrid_coeffs = unit_cgrid%coeffs1(xq_test, unit_S)
      ! CALL check(error, cgrid_coeffs(1,1), 3._dp)

      DO iq = 1, g%nq
         SGN = -1
         IF(MOD(iq, 2) == 1) SGN = 1
         xq_test = g%xq(:,iq)
         delta_x = (/0._dp, 0._dp, REAL(SGN, DP)/g%nq/)
         cgrid_coeffs = unit_cgrid%coeffs1(xq_test + delta_x, unit_S)
         iq_ = iq + SGN
         IF( iq_ <= 0 .OR. (iq-1)/g%n(3) /= (iq_-1)/g%n(3) ) iq_ = iq_ - SGN*g%n(3)
         CALL check(error, REAL(cgrid_coeffs(1,1), DP), &
            REAL(unit_cgrid%coefficients(1,1,iq) * (1 - REAL(g%n(3), DP)/g%nq) + unit_cgrid%coefficients(1,1,iq_)*g%n(3)/g%nq, DP), &
            "testing 3rd coordinate", thr=THRESHOLD)
         IF (allocated(error)) THEN
            WRITE(*,'(A, 3(F5.2, 1X))') " >>>>>>>>>> xq + delta_x     ",  g%xq(:,iq) + delta_x
            WRITE(*,*) ">>>>>>>>>> iq", iq, "SGN", SGN
            WRITE(*,*) ">>>>>>>>>> iq_", iq_
            WRITE(*,*) unit_cgrid%coefficients(1,1,iq), unit_cgrid%coefficients(1,1,iq_)
         END IF
      END DO

      DO iq = 1, g%nq
         SGN = -1
         IF(MOD(iq, 2) == 1) SGN = 1
         xq_test = g%xq(:,iq)
         delta_x = (/-4._dp, REAL(SGN, DP)/g%nq, 3._dp/)
         cgrid_coeffs = unit_cgrid%coeffs1(xq_test + delta_x, unit_S)
         iq_ = iq + SGN*g%n(3)
         IF( iq_ <= 0 .OR. (iq-1)/(g%n(3)*g%n(2)) /= (iq_-1)/(g%n(3)*g%n(2)) ) iq_ = iq_ - SGN*g%n(3)*g%n(2)
         CALL check(error, REAL(cgrid_coeffs(1,1), DP), &
            REAL(unit_cgrid%coefficients(1,1,iq) * (1 - REAL(g%n(2), DP)/g%nq) + unit_cgrid%coefficients(1,1,iq_)*g%n(2)/g%nq, DP), &
            "testing 2nd coordinate", thr=THRESHOLD)
         IF (allocated(error)) THEN
            WRITE(*,*) ">>>>>>>>>> testing 2nd coordinate"
            WRITE(*,'(A, 3(F5.2, 1X))') " >>>>>>>>>> xq + delta_x     ",  g%xq(:,iq) + delta_x
            WRITE(*,*) ">>>>>>>>>> iq", iq, "SGN", SGN
            WRITE(*,*) ">>>>>>>>>> iq_", iq_
         END IF
      END DO
   END SUBROUTINE test_cgrid_coeffs1

   SUBROUTINE set_Si_example(error)
      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      CHARACTER(len=256), PARAMETER :: path = '../Examples/Silicon/reference/mat2R'

      CALL read_fc2(path, Si_S,  Si_fc2)
      CALL aux_system(Si_S)
      CALL impose_asr2('simple', Si_S%nat, Si_fc2, Si_S%zeu)
      CALL div_mass_fc2(Si_S, Si_fc2)

      CALL check(error, Si_S%nat3, 6)
      IF(ALLOCATED(error)) WRITE(*,*) "Something wrong in reading", path
      CALL Si_cgrid%setup(Si_S, Si_fc2, 1)
   END SUBROUTINE set_Si_example

   SUBROUTINE test_coeff_value(error)
      USE fc2_interpolate,  ONLY : freq_phq_safe, add_rgd_blk_d3, fftinterp_mat2_safe, mat2_diag_pure, freq_phq
      USE constants,        ONLY : tpi
      USE functions,        ONLY : rotate_d2

      TYPE(error_type), ALLOCATABLE, INTENT(OUT) :: error
      ! local variables
      INTEGER, PARAMETER :: iq = 2
      !> mean of coefficients : 7E-8
      REAL(DP) :: freq(Si_S%nat3), xq(3), arg
      COMPLEX(DP) :: U(Si_S%nat3,Si_S%nat3), D(Si_S%nat3,Si_S%nat3), UDU(Si_S%nat3,Si_S%nat3), freq_coeff(Si_S%nat3)
      INTEGER :: iR, ibnd, i, j

      D = (0._dp, 0._dp)
      ! DO iR = 1, 5
      !    WRITE(*,*) Si_fc2%yR(:,iR), Si_fc2%FC(1,1,iR)/0.3906429607*1E4
      ! END DO
      xq = Si_cgrid%qgrid%xq(:,iq)
      CALL freq_phq_safe(xq, Si_S, Si_fc2, freq, U)
      DO iR = 1, Si_cgrid%nR
         arg = tpi * SUM(xq*Si_fc2%xR(:,iR))
         freq_coeff = freq_coeff + Si_cgrid%coefficients(:,iR,iq) * CMPLX(Cos(arg), -Sin(arg), DP)
         ! freq_coeff = freq_coeff + REAL(Si_cgrid%coefficients(:,iR,iq), DP) * Cos(arg) &
         !    +  AIMAG(Si_cgrid%coefficients(:,iR,iq)) * Sin(arg)
      END DO
      ! IF(Si_S%lrigid) THEN
      !    CALL add_rgd_blk_d3(xq, Si_S, Si_fc2, D)
      !    UDU = MATMUL(TRANSPOSE(U), D)
      !    UDU = MATMUL(UDU, U)
      ! END IF
      DO ibnd = 1, Si_S%nat3
         ! IF(Si_S%lrigid) freq_coeff(ibnd) = freq_coeff(ibnd) + REAL(UDU(ibnd,ibnd), DP)
         CALL check(error, SQRT(REAL(freq_coeff(ibnd), DP)), freq(ibnd), thr=THRESHOLD)
         IF(ALLOCATED(error)) THEN
            WRITE(*,*) "---------------------------------------"
            ! WRITE(*,'(A,(3(F5.2)))') "xq =", xq
            ! WRITE(*,*) "calculated", SQRT(REAL(freq_coeff(ibnd), DP))
            WRITE(*,*) "---------------------------------------"
         END IF
      END DO
   END SUBROUTINE

END MODULE test_cgrid
