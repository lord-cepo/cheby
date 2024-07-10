
MODULE coeff_grid
   USE kinds,            ONLY: DP
   USE ph_system,        ONLY : ph_system_info
   USE q_grids,          ONLY : q_grid
   !
   EXTERNAL cryst_to_cart, errore
   !> contains all the coefficients of the approximated polynomial, in terms like <eqs|D(R)|eqs>.
   !> You only need ph_system_info and a quality factor to initiate it.
   TYPE cgrid
      INTEGER :: nR = 0
      !! equal to n(1)*n(2)*n(3) in ph/el calculation
      TYPE(q_grid) :: qgrid
      !! final q grid in which coefficients are defined
      COMPLEX(DP), ALLOCATABLE :: coefficients(:,:,:)
      !! <eqs|D(R)|eqs>, size: nat3,nR,nout
   CONTAINS
      PROCEDURE :: setup => cgrid_setup
      PROCEDURE :: destroy => cgrid_destroy
      PROCEDURE :: coeffs => cgrid_coeffs
      PROCEDURE :: coeffs1 => cgrid_coeffs1
   END TYPE cgrid
   !
CONTAINS

   !> Safe implementation of indexing (i,j,k) => index, with Periodic Boundary Conditions and 
   !> negative (i,j,k) allowed
   FUNCTION index3(grid, ijk) result(i)
      TYPE(q_grid), INTENT(IN) :: grid
      INTEGER, INTENT(IN) :: ijk(3)
      INTEGER :: i, ijk_(3)

      ijk_ = MOD(ijk, grid%n)
      WHERE(ijk_ < 0)
         ijk_ = ijk_ + grid%n
      END WHERE
      i = MOD(grid%n(2)*grid%n(3)*ijk_(1) + grid%n(3)*ijk_(2) + ijk_(3), grid%nqtot)
      i = i + 1
   END FUNCTION index3

   !> trivial sign function
   function signum(x) result(res)
      implicit none
      real(dp), intent(in) :: x
      integer :: res

      IF (x < 0.0_dp) THEN
         res = -1
      ELSE
         res = 1
      END IF
   end function signum

   !> Sets nR, the q_grid and the coefficients <eqs|D(R)|eqs>
   SUBROUTINE cgrid_setup(grid, S, fc2, quality)
      USE input_fc,         ONLY : forceconst2_grid
      USE fc2_interpolate,  ONLY : freq_phq_safe
      USE q_grids,          ONLY : q_grid, setup_simple_grid

      IMPLICIT NONE

      TYPE(ph_system_info), INTENT(IN) :: S
      !! ph_system_info from read_fc2
      TYPE(forceconst2_grid), INTENT(IN) :: fc2
      !! force constants from read_fc2
      INTEGER, INTENT(IN) :: quality
      !! the q-grid will have dimensions multiplied by this number
      CLASS(cgrid), INTENT(INOUT) :: grid
      ! local variables
      TYPE(q_grid) :: fine_grid
      INTEGER :: iq, iR, ibnd, j, k
      REAL(DP) :: freq(S%nat3)
      COMPLEX(DP) :: U(S%nat3,S%nat3), UDU(S%nat3,S%nat3)

      grid%nR = fc2%n_R

      CALL setup_simple_grid(S%bg, fc2%nq(1)*quality, fc2%nq(2)*quality, fc2%nq(3)*quality, fine_grid)
      grid%qgrid = fine_grid
      IF(ALLOCATED(grid%coefficients)) CALL errore('cgrid_setup', 'routine called twice', 1)
      ALLOCATE(grid%coefficients(S%nat3,grid%nR,fine_grid%nq))

      DO iq = 1, fine_grid%nq
         CALL freq_phq_safe(fine_grid%xq(:,iq), S, fc2, freq, U)
         DO iR = 1, grid%nR
            UDU = MATMUL(TRANSPOSE(CONJG(U)), MATMUL(fc2%FC(:,:,iR),U))
            DO ibnd = 1, S%nat3
               grid%coefficients(ibnd,iR,iq) = UDU(ibnd, ibnd)
               ! WRITE(*,*) grid%coefficients(ibnd,iR,iq)
            !    DO k = 1, S%nat3
            !       DO j = 1, S%nat3
            !          grid%coefficients(ibnd,iR,iq) = grid%coefficients(ibnd,iR,iq) &
            !             + CONJG(U(j, ibnd)) * fc2%FC(j,k,iR) * U(k, ibnd)
            !       END DO
            !    END DO
            END DO
         END DO
      END DO
   END SUBROUTINE

   SUBROUTINE cgrid_destroy(grid)
      CLASS(cgrid), INTENT(INOUT) :: grid
      IF(ALLOCATED(grid%coefficients)) DEALLOCATE(grid%coefficients)
      grid%nR = 0
      CALL grid%qgrid%destroy()
   END SUBROUTINE

   !> EVALUATION METHOD: takes the point in the grid the nearest wrt xq_cart
   FUNCTION cgrid_coeffs(grid, xq_cart, S)
      CLASS(cgrid), INTENT(IN) :: grid
      REAL(DP), INTENT(IN) :: xq_cart(3)
      !! q-point in cartesian coordinates
      TYPE(ph_system_info) :: S
      !! ph_system_info

      ! local variables
      REAL(DP) cgrid_coeffs(S%nat3,grid%nR)
      REAL(DP) :: xq_cryst(3)
      INTEGER :: iq, ijk(3)

      xq_cryst = xq_cart
      CALL cryst_to_cart(1, xq_cryst, S%at, -1)
      ijk = NINT(xq_cryst * grid%qgrid%n)
      iq = index3(grid%qgrid, ijk)

      cgrid_coeffs = grid%coefficients(:,:,iq)
   END FUNCTION

   !> EVALUATION METHOD: does a linear (tetrahedral) interpolation between the nearest point and the
   !> three in the quadrant
   FUNCTION cgrid_coeffs1(grid, xq_cart, S)
      CLASS(cgrid), INTENT(IN) :: grid
      REAL(DP), INTENT(IN) :: xq_cart(3)
      !! q-point in cartesian coordinates
      TYPE(ph_system_info) :: S
      !! ph_system_info

      ! local variables
      REAL(DP) cgrid_coeffs1(S%nat3,grid%nR)
      REAL(DP) :: xq_cryst(3), xyz(3)
      INTEGER :: iq0, ijk(3), ijk_(3), iq(3), i

      xq_cryst = xq_cart
      CALL cryst_to_cart(1, xq_cryst, S%at, -1)
      ijk = NINT(xq_cryst * grid%qgrid%n)
      iq0 = index3(grid%qgrid, ijk)

      xyz =  xq_cryst * grid%qgrid%n - ijk

      DO i = 1, 3
         ijk_ = ijk
         ijk_(i) = ijk(i) + signum(xyz(i))
         iq(i) = index3(grid%qgrid, ijk_)
      END DO
      xyz = ABS(xyz)
      cgrid_coeffs1 = grid%coefficients(:,:,iq0) * (1 - SUM(xyz)) + &
         grid%coefficients(:,:,iq(1)) * xyz(1) + &
         grid%coefficients(:,:,iq(2)) * xyz(2) + &
         grid%coefficients(:,:,iq(3)) * xyz(3)
   END FUNCTION

END MODULE coeff_grid
