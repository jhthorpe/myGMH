!//////////////////////////////////////////////////////////////////
!//             Parses myGMH input
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

MODULE input
  USE project
  IMPLICIT NONE

CONTAINS
!---------------------------------------------------------------------
! read_input
!       James H. Thorpe
!       -reads input files for GMH calculation
!---------------------------------------------------------------------
!Variables
! nstates               : int, number of states
! u_mat                 : 2D real*8, matrix of adiabatic dipole moments...
!                                  projected onto average of dipole ...
!                                  difference vectors
! E_mat                 : 1D real*8, matrix of adiabatic (vertical?) energies
! diab_mat              : 2D int, matrix that connects the various diabatic
!                               states
! dipoles		: 3D real*8, unprojected dipole moments
! Hab			: 2D real*8, diabatic matrix

  SUBROUTINE read_input(nstates,u_mat,E_mat,diab_mat,flag)
    IMPLICIT NONE

    !Inout
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: u_mat
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: E_mat
    INTEGER, ALLOCATABLE, DIMENSION(:,:), INTENT(INOUT) :: diab_mat
    INTEGER, INTENT(INOUT) :: nstates
    LOGICAL, INTENT(INOUT) :: flag

    !Internal
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: dipoles
    INTEGER :: id
    WRITE(*,*) "Reading input file"

    flag = .FALSE.
    CALL check_files(flag)
    IF (flag) RETURN

    !get states and coupling information
    id = 100
    OPEN(file='gmh.data',unit=id,status='old')
    READ(id,*) nstates
    WRITE(*,*) "Number of states: ", nstates
    ALLOCATE(u_mat(0:nstates-1,0:nstates-1))
    ALLOCATE(E_mat(0:nstates-1))
    ALLOCATE(diab_mat(0:nstates-1,0:nstates-1))
    ALLOCATE(dipoles(0:nstates-1,0:nstates-1,0:2))
    READ(id,*) 
    CALL read_diab_mat(id,diab_mat,nstates) 
    CLOSE(unit=id)

    !get adiabatic energies
    id = 102
    OPEN(file='adiabat.data',unit=id,status='old')
    CALL read_adiabat(id,E_mat,nstates)
    CLOSE(unit=id)

    !get dipole moments
    id = 101
    OPEN(file='dipole.data',unit=id,status='old')
    CALL read_dipole(id,dipoles,nstates)
    CLOSE(unit=id)
    CALL project_dipoles(dipoles,u_mat,nstates,diab_mat,flag)
    IF (flag) RETURN


    DEALLOCATE(dipoles)

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
    INTEGER, DIMENSION(0:,0:), INTENT(INOUT) :: diab_mat
    INTEGER, INTENT(IN) :: nstates,id
   
    !internal 
    INTEGER :: i,j,k 
   
    WRITE(*,*)
    WRITE(*,*) "----------------------------------------------"
    WRITE(*,*) "                  COUPLINGS                   "
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
! read_dipole
!	James H. Thorpe
!	-reads in unprojected dipole moments in lower-triangular form 
!---------------------------------------------------------------------
! Variables
! id		: int, read file unit
! dipoles	: 3D real*8, unprojected dipoles
! nstates	: int, number of states

  SUBROUTINE read_dipole(id,dipoles,nstates)
    IMPLICIT NONE

    !inout
    REAL(KIND=8), DIMENSION(0:,0:,0:), INTENT(INOUT) :: dipoles
    INTEGER, INTENT(IN) :: nstates,id
    
    !internal
    INTEGER :: i,j,k
    WRITE(*,*) 
    WRITE(*,*) "----------------------------------------------"
    WRITE(*,*) "                  DIPOLES                     "

    DO k=0,2

      !read each one
      DO i=0,nstates-1
        READ(id,*) dipoles(i,0:i,k)
      END DO

      !flip
      DO i=0,nstates-1
        DO j=i,nstates-1
          dipoles(i,j,k) = dipoles(j,i,k)
        END DO
      END DO

      !blank line
      IF (k .LT. 2) READ(id,*)
    END DO 

    !print
    WRITE(*,*) 
    WRITE(*,*) "X component of dipole moments:"
    DO i=0,nstates-1
      WRITE(*,*) dipoles(i,0:nstates-1,0) 
    END DO 
    WRITE(*,*) 
    WRITE(*,*) "Y component of dipole moments:"
    DO i=0,nstates-1
      WRITE(*,*) dipoles(i,0:nstates-1,1) 
    END DO 
    WRITE(*,*) 
    WRITE(*,*) "Z component of dipole moments:"
    DO i=0,nstates-1
      WRITE(*,*) dipoles(i,0:nstates-1,2) 
    END DO 

  END SUBROUTINE read_dipole

!---------------------------------------------------------------------
!	read_adiabat
!		James H. Thorpe
!	-reads in adiabatic energies from adiabat.data
!---------------------------------------------------------------------
  SUBROUTINE read_adiabat(id,E_mat,nstates)
    IMPLICIT NONE

    !INOUT
    REAL(KIND=8), DIMENSION(0:), INTENT(INOUT) :: E_mat
    INTEGER, INTENT(IN) :: nstates,id

    !INTERNAL
    REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: dummy
    INTEGER :: i

    WRITE(*,*) 
    WRITE(*,*) "----------------------------------------------"
    WRITE(*,*) "                  ADIABATS                    "
    WRITE(*,*) 
    WRITE(*,*) "Adiabatic energies:"

    ALLOCATE(dummy(0:nstates-1,0:nstates-1))
    dummy = 0
 
    DO i=0,nstates-1
      READ(id,*) E_mat(i)
      dummy(i,i) = E_mat(i)
    END DO

    DO i=0,nstates-1
      WRITE(*,*) dummy(i,0:nstates-1)
    END DO

    DEALLOCATE(dummy)

  END SUBROUTINE read_adiabat
!---------------------------------------------------------------------

END MODULE input

