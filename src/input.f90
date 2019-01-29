!//////////////////////////////////////////////////////////////////
!//             Parses myGMH input
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

MODULE input
  IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------
! read_input
!       James H. Thorpe
!       -reads input files for GMH calculation
!---------------------------------------------------------------------
!Variables
! nstates             : int, number of states
! u_mat               : 2D real*8, matrix of adiabatic dipole moments...
!                                  projected onto average of dipole ...
!                                  difference vectors
! E_mat               : 1D real*8, matrix of adiabatic (vertical?) energies
! diab_mat            : 2D int, matrix that connects the various diabatic
!                               states

  SUBROUTINE read_input(nstates,u_mat,E_mat,diab_mat,flag)
    IMPLICIT NONE

    !Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: u_mat
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: E_mat
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: diab_mat
    INTEGER, INTENT(INOUT) :: nstates
    LOGICAL, INTENT(INOUT) :: flag

    !Internal
    INTEGER :: id
    WRITE(*,*) "Reading input file"

    flag = .FALSE.
    CALL check_files(flag)
    IF (flag) RETURN

    id = 100
    OPEN(file='gmh.data',unit=id,status='old')
    READ(id,*) nstates
    WRITE(*,*) "Number of states: ", nstates
    ALLOCATE(u_mat(0:nstates-1,0:2))
    ALLOCATE(E_mat(0:nstates-1))
    ALLOCATE(diab_mat(0:nstates-1,0:nstates-1))
    READ(id,*) 
    CALL read_diab_mat(id,diab_mat,nstates) 
    CLOSE(unit=100)

  END SUBROUTINE read_input  
!---------------------------------------------------------------------
! check_files
!       James H. Thorpe
!       -checks all files that are needed are present
!---------------------------------------------------------------------
  SUBROUTINE check_files(flag)
    IMPLICIT NONE
    LOGICAL, INTENT(INOUT) :: flag

    INQUIRE(file='gmh.data',EXIST=flag)    
    IF (.NOT. flag) THEN
      WRITE(*,*) "You are missing the gmh.data file"
      flag = .TRUE.
      RETURN 
    END IF

    INQUIRE(file='dipole.data',EXIST=flag)    
    IF (.NOT. flag) THEN
      WRITE(*,*) "You are missing the dipole.data file"
      flag = .TRUE.
      RETURN 
    END IF

    INQUIRE(file='adiabat.data',EXIST=flag)    
    IF (.NOT. flag) THEN
      WRITE(*,*) "You are missing the adiabat.data file"
      flag = .TRUE.
      RETURN 
    END IF

    flag = .FALSE.
  
  END SUBROUTINE check_files
!---------------------------------------------------------------------
! read_diab_mat
!       James H. Thorpe
!       -read which diabatic states connect to which
!       -information is stored in lower-triangular form
!---------------------------------------------------------------------
!Variables
! id            : int, read file unit
! diab_mat      : 2D int, diabatic connection matrix
! nstates       : int, number of states

  SUBROUTINE read_diab_mat(id,diab_mat,nstates)
    IMPLICIT NONE
     
    !inout
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: diab_mat
    INTEGER, INTENT(IN) :: nstates,id
   
    !internal 
    INTEGER :: i,j,k 
   
    WRITE(*,*) 
    WRITE(*,*) "Coupling matrix:" 
    diab_mat = 0

    !read states
    DO i=0,nstates-1
      READ(id,*) diab_mat(i,0:i)
    END DO

    !mirror it 
    DO i=0,nstates-1
      DO j=i,nstates-1
        diab_mat(i,j) = diab_mat(j,i)
      END DO
    END DO
   
    !print states 
    DO i=0,nstates-1
      WRITE(*,*) diab_mat(i,0:nstates-1) 
    END DO 

  END SUBROUTINE read_diab_mat
!---------------------------------------------------------------------

END MODULE input

