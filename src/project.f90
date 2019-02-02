!//////////////////////////////////////////////////////////////////
!//             Module for finding projections of dipole moments
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

MODULE project
  IMPLICIT NONE
  
  CONTAINS

!---------------------------------------------------------------------
!	project_dipoles
!		James H. Thorpe
!	-subroutine that projects dipoles onto average of dipole
!	  difference vectors
!---------------------------------------------------------------------
! Variables
! dipoles	: 3D real*8, unprojected dipole moments
! u_mat		: 2D real*8, projected dipole components
! nstates	: int, number of matrix elements
! diab_mat	: 2D int, matrix of corresponding elements
! ntrans	: int, number of adiabat transitions
! udiff		: 1D real*8, averaged dipole difference

  SUBROUTINE project_dipoles(dipoles,u_mat,nstates,diab_mat,flag)
    IMPLICIT NONE
    !inout
    REAL(KIND=8), DIMENSION(0:,0:), INTENT(INOUT) :: u_mat
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(IN) :: dipoles
    INTEGER, DIMENSION(0:,0:), INTENT(IN) :: diab_mat
    LOGICAL, INTENT(INOUT) :: flag 
    INTEGER, INTENT(IN) :: nstates
    !internal
    REAL(KIND=8), DIMENSION(0:2) :: udiff
    REAL(KIND=8) :: norm
    INTEGER :: ntrans,ctv
    INTEGER :: i,j,k

    flag = .FALSE.

    CALL check_ctv(ctv)

    !build average dipole difference
    ntrans = 0
    udiff = 0
    DO i=0,nstates-2
      DO j=i+1,nstates-1
        ntrans = ntrans + 1
        udiff(0) = udiff(0) + dipoles(i,i,0) - dipoles(j,j,0) 
        udiff(1) = udiff(1) + dipoles(i,i,1) - dipoles(j,j,1) 
        udiff(2) = udiff(2) + dipoles(i,i,2) - dipoles(j,j,2) 
      END DO
    END DO 

    !check for special symmetric case
    IF (MAXVAL(ABS(udiff)) .LT. 1.0D-16 .AND. ctv .NE. 1) THEN
      WRITE(*,*)
      WRITE(*,*) "You have a truely symmetric system..."
      WRITE(*,*) "Using transition moments as the effective vector"
      udiff = 0.0D0
      ntrans = 0
      DO i=0,nstates-1
        DO j=0,i-1
          ntrans = ntrans + 1
          udiff(0) = udiff(0) + dipoles(i,j,0)
          udiff(1) = udiff(1) + dipoles(i,j,1)
          udiff(2) = udiff(2) + dipoles(i,j,2)
        END DO
      END DO 
    !check for my alternate scheme
    ELSE IF (ctv .EQ. 1) THEN
      DO i=0,nstates-1
        DO j=0,i-1
          ntrans = ntrans + 1
          udiff(0) = udiff(0) + dipoles(i,j,0)
          udiff(1) = udiff(1) + dipoles(i,j,1)
          udiff(2) = udiff(2) + dipoles(i,j,2)
        END DO
      END DO 
    END IF

    udiff  = udiff/(1.0D0*ntrans)

    !print average dipole differences vector
    WRITE(*,*)
    WRITE(*,*) "Average adiabat dipole differences vector:"
    WRITE(*,'(2x,A2,999(F15.10))') "X:", udiff(0)
    WRITE(*,'(2x,A2,999(F15.10))') "Y:", udiff(1)
    WRITE(*,'(2x,A2,999(F15.10))') "Z:", udiff(2)

    DO i=0,nstates-1
      DO j=0,nstates-1
        u_mat(i,j) = scalar_proj(dipoles(i,j,0:2),udiff(0:2))
      END DO
    END DO
    
    WRITE(*,*) 
    WRITE(*,*) "Dipoles projected on <Δμ12> :" 
    DO i=0,nstates-1 
      WRITE(*,'(999(F15.10))') u_mat(i,0:nstates-1)
    END DO

  END SUBROUTINE project_dipoles

!---------------------------------------------------------------------
! scalar_proj
!		James H. Thorpe
!	-calcuates scaler projection of vector a onto b	
!---------------------------------------------------------------------
  REAL(KIND=8) FUNCTION scalar_proj(a,b)
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(0:2), INTENT(IN) :: a,b
    REAL(KIND=8) :: bnorm,temp

    bnorm = SQRT(b(0)**2.0+b(1)**2.0+b(2)**2.0)
    temp = a(0)*b(0)+a(1)*b(1)+a(2)*b(2)
    temp = temp / bnorm
    scalar_proj = temp

  END FUNCTION scalar_proj
!---------------------------------------------------------------------
!	check_ctv
!		James H. Thorpe
!	-checks for the %ctv keyword, which can be used to control
!	the scheme for the projection of the charge transfer vector
!---------------------------------------------------------------------
  SUBROUTINE check_ctv(ctv)
    IMPLICIT NONE
    !inout
    INTEGER, INTENT(INOUT) :: ctv
    !internal
    INTEGER :: fline, exitval
    LOGICAL :: flag

    CALL EXECUTE_COMMAND_LINE('grep -q \%cvt gmh.data',.TRUE.,exitval) 

    IF (exitval .EQ. 0) THEN
      WRITE(*,*) 
      WRITE(*,*) "WARNING"
      WRITE(*,*) "Using average of dipole differences and"
      WRITE(*,*) "transition dipoles as charge transfer vector." 
      WRITE(*,*)
      WRITE(*,*) "This has not been tested or published..."
      WRITE(*,*) "be very very careful."
      ctv = 1
    ELSE
      ctv = 0
    END IF 
    
  END SUBROUTINE check_ctv
!------------------------------------------------------------------
!               Check if file exists and get number of lines    
!------------------------------------------------------------------
  SUBROUTINE getfline(flag,fline)
    IMPLICIT NONE
    !Inout
    INTEGER, INTENT(INOUT) :: fline
    LOGICAL, INTENT(INOUT) :: flag
    !Internal
    INTEGER :: io
    LOGICAL :: ex
    flag = .FALSE.
    INQUIRE(file='gmh.data',EXIST=ex)
    IF (.NOT. ex) THEN
      flag = .TRUE.
      fline = -1
      RETURN
    END IF
    fline = 0
    io = 0
    OPEN(unit=1,file='gmh.data',status='old',access='sequential')
    DO WHILE (io .EQ. 0)
      READ(1,*,iostat=io)
      IF (io .EQ. 0) fline = fline + 1
    END DO
    CLOSE(unit=1)
  END SUBROUTINE getfline
!------------------------------------------------------------------


END MODULE project
