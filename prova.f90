PROGRAM prova
   integer :: a(2,2), b(2,2), i, j, k, c(2,2), d(2)

   DO j = 1, 2
      DO i = 1, 2
         a(i,j) = i+j
         b(i,j) = i-j
      ENDDO
   ENDDO

   DO i = 1, 2
      DO j = 1, 2
         DO k = 1, 2
            C(i, j) = C(i, j) + A(i, k) * B(k, j)
         END DO
      END DO
   END DO

   WRITE(*,*) c(1,:)
   c = MATMUL(a,b)
   WRITE(*,*) c(1,:)
   d = (/1,-1/)
   WRITE(*,*) MATMUL(b,d)


END PROGRAM prova
