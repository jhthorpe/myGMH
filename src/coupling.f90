!//////////////////////////////////////////////////////////////////
!//             Module for finding coupling of GMH diabats 
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

MODULE coupling
  USE hcbd
  IMPLICIT NONE

  CONTAINS

!---------------------------------------------------------------------
!	calc_Hab
!		James H. Thorpe
!	-control for calculation Hab couplings between GMH diabats
!---------------------------------------------------------------------
! Variables
! Hab		: 2D real*8, diabatic hamiltonian 
! n	: int, number of states 
! u_mat		: 2D real*8, projected dipole matrix 
! E_mat		: 1D real*8, adiabatic hamiltonian
! diab_mat	: int, coupling connections matrix 
! flag		: bool, flag
! S_mat		: 2D real*8, unitary transform matrix from ad to di

  SUBROUTINE calc_Hab(Hab,n,u_mat,E_mat,diab_mat,flag)
    IMPLICIT NONE

    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Hab
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: u_mat
    REAL(KIND=8), DIMENSION(0:), INTENT(IN) :: E_mat
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: diab_mat
    INTEGER, INTENT(IN) :: n
    LOGICAL, INTENT(INOUT) :: flag

    !internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: S_mat,int_mat,dummy,A
    INTEGER, ALLOCATABLE, DIMENSION(:,:) :: one
    INTEGER :: info
    INTEGER :: i,j

    flag = .FALSE.

    WRITE(*,*) 
    WRITE(*,*) "----------------------------------------------"
    WRITE(*,*) "                  DIABATS                     "
    WRITE(*,*) 
    WRITE(*,*) "Calculating diabatic Hamiltonian"

    ALLOCATE(S_mat(0:n-1,0:n-1))
    ALLOCATE(int_mat(0:n-1,0:n-1))
    ALLOCATE(dummy(0:n-1,0:n-1))
    ALLOCATE(A(0:n-1,0:n-1))
    ALLOCATE(one(0:n-1,0:n-1))
    S_mat = 0
    dummy = 0
    one = 1
    DO i=0,n-1
      dummy(i,i) = E_mat(i)
    END DO
    A = u_mat

    !1) Diagonalize the dipoles
    WRITE(*,*)
    WRITE(*,*) "Diagonalizing dipole matrix to μab = 0"
    !WARNING : don't use the Lapack routines, they change the ordering
    !CALL diag_mat(n,u_mat,S_mat) !LAPACK full diagonalization
    !CALL diag_mat(n,A,S_mat) !LAPACK full diagonalization
    CALL block_diag(n,one,A,S_mat,info)
    IF (info .NE. 0) THEN
      flag = .TRUE.
      RETURN
    END IF
    WRITE(*,*) 
    WRITE(*,*) "Fully Diagonalized μ matrix:"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') A(i,0:n-1)
    END DO
    WRITE(*,*) 
    WRITE(*,*) "Transformation matrix:"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') S_mat(i,0:n-1)
    END DO

    !2) Use transform matrix to get fully diabatic hamiltonian
    CALL linal_xy_2Dreal8(n,n,n,dummy,TRANSPOSE(S_mat),int_mat) 
    CALL linal_xy_2Dreal8(n,n,n,S_mat,int_mat,Hab) 
    WRITE(*,*) 
    WRITE(*,*) "Fully Diabatic Hamiltonian"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') Hab(i,0:n-1)
    END DO

    !3) Get locally adiabatic GMH diabats
    WRITE(*,*) "----------------------------------------------"
    WRITE(*,*) "               GMH DIABATS                    "

    !block diagonalize  diabatic Hamiltonian
    CALL block_diag(n,ABS(one-diab_mat),Hab,S_mat,info)
    IF (info .NE. 0) THEN
      flag = .TRUE.
      RETURN
    END IF
    WRITE(*,*) 
    WRITE(*,*) "GMH (locally-adiabatic) Diabatic Hamiltonian"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') Hab(i,0:n-1)
    END DO 
    WRITE(*,*) 
    WRITE(*,*) "Original Adiabatic Hamiltonian"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') dummy(i,0:n-1)
    END DO

    !transform dipole matrix
    CALL linal_xy_2Dreal8(n,n,n,A,TRANSPOSE(S_mat),int_mat) 
    CALL linal_xy_2Dreal8(n,n,n,S_mat,int_mat,dummy) 
    WRITE(*,*) 
    WRITE(*,*) "GMH Dipole Matrix"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') dummy(i,0:n-1)
    END DO
    WRITE(*,*)
    WRITE(*,*) "Original Dipole Matrix"
    DO i=0,n-1
      WRITE(*,'(999(F15.10))') u_mat(i,0:n-1)
    END DO

  END SUBROUTINE calc_Hab
!---------------------------------------------------------------------
!	diag_umat
!		James H. Thorpe
!	-diagonalize an nxn matrix, return transform matrix
!---------------------------------------------------------------------
! n		: int, rank of matrix
! mat		: 2D real*8, matrix to diagonalize
! S		: 2D real*8, transform matrix

  SUBROUTINE diag_mat(n,mat,S)
    IMPLICIT NONE

    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: S
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: mat 
    INTEGER, INTENT(IN) :: n

    !internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A,dummy!,B,C
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: W,WORK
    INTEGER :: LWORK,INFO,i

    LWORK = 4*n+1
    ALLOCATE(A(0:n-1,0:n-1)) 
    !ALLOCATE(B(0:n-1,0:n-1)) 
    !ALLOCATE(C(0:n-1,0:n-1)) 
    ALLOCATE(dummy(0:n-1,0:n-1)) 
    ALLOCATE(WORK(0:LWORK-1))
    ALLOCATE(W(0:n-1))

    WORK = 0
    W = 0
    INFO = 0
    A = mat

    WRITE(*,*) "WARNING, DIAGONALIZATION HAS BEEN SILENCED"
    !diagonalize
    !CALL DSYEV('V','U',n,A(0:n-1,0:n-1),n,W,WORK(0:LWORK-1),LWORK,INFO)
    !CALL DSYEV('V','L',n,A(0:n-1,0:n-1),n,W,WORK(0:LWORK-1),LWORK,INFO)
    dummy = 0
    WRITE(*,*)
    WRITE(*,*) "Eigenvalues matrix:"
    DO i=0,n-1
      dummy(i,i) = W(i)
      WRITE(*,*) dummy(i,0:n-1)
    END DO 
    WRITE(*,*)
    WRITE(*,*) "Unitary Transform Matrix"
    DO i=0,n-1
      WRITE(*,*) A(i,0:n-1)
    END DO

    S = A
    mat = 0
    DO i=0,n-1
      mat(i,i) = W(i)
    END DO
    !B = 0
    !C = 0
    !
    !WRITE(*,*) "TESTING TESTING TESTING" 
    !CALL linal_xy_2Dreal8(n,n,n,mat,S,B)
    !CALL linal_xTy_2Dreal8(n,n,n,S,B,C)
    !WRITE(*,*) C
    !WRITE(*,*) ""
    !B = 0
    !C = 0
    !CALL linal_xy_2Dreal8(n,n,n,dummy,S,B)
    !CALL linal_xTy_2Dreal8(n,n,n,S,B,C)
    !WRITE(*,*) C

    DEALLOCATE(A)
    !DEALLOCATE(B)
    !DEALLOCATE(C)
    DEALLOCATE(dummy)
    DEALLOCATE(WORK)
    DEALLOCATE(W) 

  END SUBROUTINE diag_mat
!---------------------------------------------------------------------
!	inv_mat
!		James H. Thorpe
!	-invert a matrix via LU factorization and then inversion
!	- not needed here, but maybe in future work
!---------------------------------------------------------------------
! n		: int, rank of matrix
! mat		: 2D real*8, matrix to be inverted
! imat		: 2D real*8, inverted matrix
  SUBROUTINE inv_mat(n,mat,imat)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: imat
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: mat 
    INTEGER, INTENT(IN) :: n
    !internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: A
    INTEGER, ALLOCATABLE, DIMENSION(:) :: IPIV,WORK
    INTEGER :: INFO,LWORK

    LWORK = 2*n
    ALLOCATE(A(0:n-1,0:n-1))
    ALLOCATE(IPIV(0:n-1))
    ALLOCATE(WORK(0:LWORK-1)) 

    A = mat
    WORK = 0
   
    WRITE(*,*) "WARNING - INVERSION HAS BEEN SILENCED"

    !1) get LU decomposition of matrix
    !CALL DGETRF(n,n,A(0:n-1,0:n-1),n,IPIV(0:n-1),INFO)
    !2) invert 
    !CALL DGETRI(n,A(0:n-1,0:n-1),n,IPIV(0:n-1),WORK,LWORK,INFO) 
 
    imat = A
   
    DEALLOCATE(A)
    DEALLOCATE(IPIV) 
    DEALLOCATE(WORK)

  END SUBROUTINE inv_mat

!---------------------------------------------------------------------
!       linal_xTy_2Dreal8
!               James H. Thorpe
!               Dec 9, 2018
!       -matrix multiplication using BLAS
!       -performs Z = X**T.Y
!---------------------------------------------------------------------
  ! n1          : int, number of rows of X**T 
  ! n2          : int, number of cols of Y
  ! n3          : int, number of cols/rows of X**T/Y
  ! X           : 2D real*8, matrix to be transposed, size [n2,n1]
  ! Y           : 2D real*8, matrix size [n3,n2]
  ! Z           : 2D real*8, matrix size [n1,n3]
  SUBROUTINE linal_xTy_2Dreal8(n1,n2,n3,X,Y,Z)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Z
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: X,Y
    INTEGER, INTENT(IN) :: n1,n2,n3

    Z = 0.0D0
!    CALL DGEMM('T','N',n1,n2,n3,1.0D0,X(0:n3-1,0:n1-1),n3,&
!            Y(0:n3-1,0:n2-1),n3,0.0D0,Z(0:n1-1,0:n2-1),n1)

  END SUBROUTINE linal_xTy_2Dreal8

!---------------------------------------------------------------------
!       linal_xy_2Dreal8
!               James H. Thorpe
!               Dec 9, 2018
!       -matrix multiplication using BLAS  
!       -performs Z = X.Y
!	-hand-coded, now, so we don't need BLAS/LAPACK for anything
!---------------------------------------------------------------------
  ! n1          : int, number of rows of X
  ! n2          : int, number of cols of Y
  ! n3          : int, number of cols/rows of X/Y
  ! X           : 2D real*8, matrix size [n1,n2]
  ! Y           : 2D real*8, matrix size [n2,n3]
  ! Z           : 2D real*8, matrix size [n1,n3]
  SUBROUTINE linal_xy_2Dreal8(n1,n2,n3,X,Y,Z)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: Z
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(IN) :: X,Y
    INTEGER, INTENT(IN) :: n1,n2,n3
    INTEGER :: i,j,k

    Z = 0.0D0
!    CALL DGEMM('N','N',n1,n2,n3,1.0D0,X(0:n1-1,0:n3-1),n1,&
!            Y(0:n3-1,0:n2-1),n3,0.0D0,Z(0:n1-1,0:n2-1),n1)

    DO i=0,n1-1
      DO j=0,n2-1
        DO k=0,n3-1
          Z(i,j) = Z(i,j) + X(i,k)*Y(k,j) 
        END DO
      END DO
    END DO

  END SUBROUTINE linal_xy_2Dreal8
!---------------------------------------------------------------------

END MODULE coupling

