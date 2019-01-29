!//////////////////////////////////////////////////////////////////
!//             Control program for myGMH
!//
!//             James H. Thorpe, in the Group of John Stanton
!//             The Univeristy of Florida
!//
!//////////////////////////////////////////////////////////////////

!-------------------------------------------------------------------------
! input
!       James H. Thorpe
!       -control program for GMH
!-------------------------------------------------------------------------
!Variables
! nstates             : int, number of states
! u_mat               : 2D real*8, matrix of adiabatic dipole moments...
!                                  projected onto average of dipole ...
!                                  difference vectors
! E_mat               : 1D real*8, matrix of adiabatic energies
! diab_mat            : 2D int, matrix that connects the various diabatic
!                               states

PROGRAM gmh
  USE input
  !USE linal
  !USE calc
  IMPLICIT NONE

  !Internal
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:) :: u_mat
  REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: E_mat
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: diab_mat 
  INTEGER :: nstates
  LOGICAL :: flag

  CALL print_start()
  CALL read_input(nstates,u_mat,E_mat,diab_mat,flag)
  IF (flag) STOP 1

CONTAINS
  
!-------------------------------------------------------------------------
! print_start
!       James H. Thorpe
!       -prints begining of program
!-------------------------------------------------------------------------
  SUBROUTINE print_start()
    IMPLICIT NONE
    WRITE(*,*) 
    WRITE(*,*) "Beging GMH calculation"

  END SUBROUTINE print_start
  
!-------------------------------------------------------------------------

END PROGRAM gmh
