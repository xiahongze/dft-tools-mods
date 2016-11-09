PROGRAM MATDYNREADER
    IMPLICIT NONE
    ! index
    INTEGER :: i,j
    ! variables for the input list
    INTEGER :: nks,nbnd
    REAL(8), ALLOCATABLE :: freq(:,:),qvec(:,:)
    COMPLEX(8), ALLOCATABLE :: z(:,:,:)
    CHARACTER(LEN=256) :: flfreq,flmode

    WRITE(*,*) "your freq file?"
    READ(5,*) flfreq
    ! read from QE output
    CALL readband(flfreq,freq,qvec,nbnd,nks)
    ! 
    WRITE(*,*) "your mode file?"
    READ(5,*) flmode
    ! read from QE output
    CALL readmode(flmode,z,nbnd,nks)
    ! End
    DEALLOCATE( freq,qvec,z )

CONTAINS
    subroutine readband(flfreq,eval,qvec,nbnd,nks)
        CHARACTER(LEN=256), intent(in) :: flfreq
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: eval, qvec
        integer :: nbnd,nks
        ! local
        integer :: i,j
        NAMELIST /plot/ nbnd,nks
        open(unit=1001,file=TRIM(flfreq),action='read')
        read(1001,nml=plot)
        ALLOCATE( eval(nbnd,nks),qvec(3,nks) )
        do i = 1,nks
            read(1001,*) (qvec(j,i),j=1,3)
            read(1001,*) (eval(j,i),j=1,nbnd)
        end do
        close(1001)
    end subroutine readband
    
    subroutine readmode(flmode,z,nbnd,nks)
        CHARACTER(LEN=256), intent(in) :: flmode
        COMPLEX(8), ALLOCATABLE, intent(inout) :: z(:,:,:)
        integer, intent(in) :: nbnd,nks
        ! local
        integer :: i,j,k,nat,ipol
        !
        nat = nbnd/3
        ! allocate mode array
        ALLOCATE(z(nbnd,nbnd,nks))
        !
        open(unit=1001,file=TRIM(flmode),action='read')
        do i = 1,nks
            READ(1001,*) ! skip reading diagonalizing
            READ(1001,*) ! skip reading q
            READ(1001,*) ! skip blank line
            READ(1001,*) ! skip reading ****
            do j = 1,nbnd
!                 WRITE(*,*) "kpoint",i,"lambda",j ! test
                READ(1001,*) ! skip reading freq
                do k =1,nat
                     READ(1001,9020) (z((k-1)*3+ipol,j,i),ipol=1,3)
!                      WRITE(*,9010) (z((k-1)*3+ipol,j,i),ipol=1,3)
                end do
            end do
            READ(1001,*) ! skip reading ****
        end do
    9010 format (1x,'(',3 (f10.6,1x,f10.6,3x),')') ! test
    9020 format (2x,3 (f10.6,1x,f10.6,3x),1x)
    end subroutine readmode
            
END PROGRAM MATDYNREADER