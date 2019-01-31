!//////////////////////////////////////////////////////////////////
!//             Module for hand-coded block diagonalizer 
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

MODULE hcbd
  IMPLICIT NONE

  CONTAINS

!---------------------------------------------------------------------
!	block_diag
!		James H. Thorpe
!	- uses Jacobi Transforms to block-diagonalize a sym. real matrix 
!	  in the form specified by `blocks`
!	- very poorly coded, but it will do for small matrices
!	- based off algorithms described in Numerical Recipies
!	- uses upper triangular form
!---------------------------------------------------------------------
! Variables
! n		: int, size of matrix
! blocks	: 2D int, block structure of matrix
! A		: 2D real*8, real, symmetric matrix to be block-diag
! V		: 2D real*8, transform matrix 
! info		: int, information about how calculation went
! w1,w2		: 1D real*8, working vectors

  SUBROUTINE block_diag(n,blocks,A,V,info)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: A,V
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: blocks
    INTEGER, INTENT(INOUT) :: info
    INTEGER, INTENT(IN) :: n
    !internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: w1,w2
    REAL(KIND=8) :: conv,tol
    INTEGER :: i,j

    info = 0
    tol = 1.0D-16

    !initialize
    V = 0
    DO i=0,n-1
      V(i,i) = 1
    END DO 

    DO i=0,49 !max of 50 sweeps
    
      !check convergence
      CALL get_conv(conv,A,blocks,n)    
      IF (conv .LT. tol) EXIT

      !otherwise, do Jacobi sweep
      CALL jacobi_sweep(n,blocks,A,V,i)

    END DO

    !Mirror
    DO i=0,n-1
      DO j=i+1,n-1
        A(j,i) = A(i,j)
      END DO
    END DO

    IF (conv .LT. tol) THEN
      info = 0
    ELSE
      info = 1
    END IF 

  END SUBROUTINE block_diag

!---------------------------------------------------------------------
!	get_conv
!		James H. Thorpe
!	- checks convergence of block-diagonal matrix
!---------------------------------------------------------------------
! conv		: real*8, convergence
! A		: 2D real*8, matrix
! blocks	: 2D int, block info
! n		: int, rank

  SUBROUTINE get_conv(conv,A,blocks,n)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: A 
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: blocks
    REAL(KIND=8), INTENT(INOUT) :: conv
    INTEGER, INTENT(IN) :: n
    !internal
    INTEGER :: i,j

    conv = 0.0D0

    DO i=0,n-1
      DO j=i+1,n-1
        conv = conv + 2*(A(i,j)**2.0D0)*blocks(i,j)
      END DO 
    END DO
    
  END SUBROUTINE get_conv

!---------------------------------------------------------------------
!	jacobi_sweep
!		James H. Thorpe
!	-performs jacobi sweep of matrix to diagonalize elements in 
!	 blocks, in upper triangular form
!	- uses code from NR 77
!---------------------------------------------------------------------
! n		: int, rank of matrix
! A		: 2D real*8, matrix to be reduced
! V		: 2D real*8, eigenvectors 
! w1,w2		: 1D real*8, working vectors
! blocks	: 2D int, block info
! p,q		: int, indicies
! c,s,t		: real*8, same as in Num. Rec.
! th		: real*8, theta(θ) in Num. Rec. 
! ta		: real*8, tau(τ) in Num. Rec.
! h		: real*8, differences between A(q,q) and A(p,p)
! itr		: int, iteration of sweep

  SUBROUTINE jacobi_sweep(n,blocks,A,V,itr)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: A,V 
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: blocks
    INTEGER, INTENT(IN) :: n,itr
    !internal
    REAL(KIND=8) :: c,s,t,o,th,ta,tol,h,g,f
    INTEGER :: p,q,i,j,k,r

    tol = 1.0D-16

    !sweep through indices
    DO p=0,n-2
      DO q=p+1,n-1
        IF (blocks(p,q) .EQ. 0) CYCLE !skip these elements
        g = 100.0D0*ABS(A(p,q))       !numerical stuff

        IF ( (itr .GT. 4) .AND. (ABS(A(p,p))+g .EQ. ABS(A(p,p))) &
             .AND. (ABS(A(q,q))+g .EQ. ABS(A(q,q))) ) THEN
          A(p,q) = 0.0D0
        ELSE IF (ABS(A(p,q)) .GT. tol) THEN !element is not zero, eliminate
          h = A(q,q) - A(p,p)
          IF ((ABS(h)+g) .EQ. ABS(h)) THEN !machine overflow
            t = A(p,q)/h
          ELSE !normal case
            th = 0.5D0*h/A(p,q)
            t = 1.0D0/(ABS(th)+SQRT(1.0D0 + th**2.0D0))
            IF (th .LT. 0) t = -t
          END IF !machine overflow

          !get parameters
          c = 1.0D0/SQRT(1.0D0 + t**2.0D0)
          s = t*c
          ta = s/(1.0D0 + c)
          f = t*A(p,q) 

          !change A matrix 
          DO r=0,p-1 !1 <= r < p         
            g = A(r,p)  
            h = A(r,q)
            A(r,p) =  g - s*(h + ta*g)
            A(r,q) =  h + s*(g - ta*h)
          END DO

          A(p,p) = A(p,p) - f !r = p, points p,p and p,q 
          A(p,q) = 0.0D0 

          DO r=p+1,q-1 !p < r < q
            g = A(p,r)  
            h = A(r,q)
            A(p,r) =  g - s*(h + ta*g)
            A(r,q) =  h + s*(g - ta*h)
          END DO

          !A(q,p) = 0.0D0 !points q,q and q,p
          A(q,q) = A(q,q) + f 

          DO r=q+1,n-1 !q < j <= n
            g = A(p,r)
            h = A(q,r)
            A(p,r) =  g - s*(h + ta*g)
            A(q,r) =  h + s*(g - ta*h)
          END DO

          !change V matrix
          DO r=0,n-1
            g = V(r,p)
            h = V(r,q)
            V(r,p) = g - s*(h + ta*g) 
            V(r,q) = h + s*(g - ta*h)
          END DO

        END IF !nonzero elements

      END DO
    END DO

  END SUBROUTINE jacobi_sweep
!---------------------------------------------------------------------
END MODULE hcbd

