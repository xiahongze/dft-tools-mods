! --------------------------------------------------------------------------------------------------
! This fortran program runs with outputs from Matdyn.x. It basically
! needs two inputs, i.e., the k-mesh, the frequencies.
! The branch of LTO should be set manually, i.g., for Silicone
! TO and LO degenerate at zone centre.
! It finds all possible decay paths for the selected phonon mode. The
! accuracy/tolerance control can be set manually or use default.
! It eventually dumps complete decay information:
! q-point, j-mode and its frequency.
! 
! Hongze Xia, 2014-03-23
! 2014-04-03: read from QE output band file directly
! 2014-04-04: read q-weight from stdin and append it to the output
! 2014-04-06: real q-weight is used
! --------------------------------------------------------------------------------------------------
! Input cards: namelist &input
!   w0  : the target LTO frequency in cm^-1. ***MUST BE SET***
!   errtol : error tolerance for calc. Default value is 0.1d-2.
!   flfreq : the filename of freq ouput by QE. default name is "freq.txt".
PROGRAM DECAYSCAN
    IMPLICIT NONE
    ! index, q-point, j-mode, no. of paths 
    INTEGER :: q,j,j1,n=1
    ! variables for the input list
    INTEGER :: nks,nbnd
    REAL(8), ALLOCATABLE :: freq(:,:),qvec(:,:),qwght(:)
    REAL(8) :: errtol, w0
    CHARACTER(LEN=256) :: flfreq
    NAMELIST /input/ w0, errtol, flfreq
    !
    ! set default values
    errtol = 0.1d-2
    flfreq = "freq.txt"
    !
    ! read the namelist from file
    READ(5,nml=input)
    ! read from QE output
    CALL readband(flfreq,freq,qvec,nbnd,nks)
    ! read q-weight
    ALLOCATE( qwght(nks) )
    READ(5,*) (qwght(q),q=1,nks)
    qwght = qwght/sum(qwght)
    ! End of initiation

    ! Scan the whole grid for possible decay paths
    do q = 1,nks
        do j = 1,nbnd
            do j1 = j,nbnd
                ! skip the target mode itself & accoustic modes
                if ( ABS(freq(j,q)-w0)<errtol*w0 .or. & 
                    & ABS(freq(j1,q)-w0)<errtol*w0) CYCLE
                ! w0 = wj+wj1 ?
                if ( ABS( freq(j,q)+freq(j1,q)-w0 )<errtol*w0 ) then
                    write(*,9000) n,qvec(:,q)
                    write(*,9001) j,j1,q,qwght(q)
                    write(*,9002) freq(j,q),freq(j1,q)
                    write(*,*)
                    n = n+1
                endif
            end do
        end do
    end do

    ! End
    DEALLOCATE( freq,qvec,qwght )
    
9000 format(5x,"Decay path No. ",i0.0,1x,"at ",'(',3 (f8.4,1x),')')
9001 format(5x,"The two modes are",i3.0,2x,"and",i3.0,"  at nq =",i5.0,"  qw = ",e10.4)
9002 format(5x,"Frequencies are",1x,f10.4,3x,"and",1x,f10.4,4x,"cm-1")

CONTAINS
    subroutine readband(flband,eval,qvec,nbnd,nks)
        CHARACTER(LEN=256), intent(in) :: flband
        REAL(8), DIMENSION(:,:), ALLOCATABLE :: eval, qvec
        integer :: nbnd,nks
        ! local
        integer :: i,j
        NAMELIST /plot/ nbnd,nks
        open(unit=1001,file=TRIM(flband),action='read')
        read(1001,nml=plot)
        ALLOCATE( eval(nbnd,nks),qvec(3,nks) )
        do i = 1,nks
            read(1001,*) (qvec(j,i),j=1,3)
            read(1001,*) (eval(j,i),j=1,nbnd)
        end do
        close(1001)
    end subroutine readband
END PROGRAM DECAYSCAN