      SUBROUTINE MKMPI_GROUPS(MASTER,MYID,NGROUPS,NPROCS,MYCOMM,
     ;                        IROW,MYNID,MYHD_COMM,MYSLV_COMM)
! 
!     Makes the mpi groups, here is the idea pictorally, say we begin
!     with 12 processors, with id's 0-11; 3 groups, and 4 rows per 
!     group, i.e. 4 rows. 
!     
!     myid:  0  1  2  3  4  5  6  7  8  9 10 11 
!     is split into 3 groups
!     irow:  0  0  0  0  1  1  1  1  2  2  2  2 
!     with new ids
!     jcol:  0  1  2  3  0  1  2  3  0  1  2  3 
!
!     So my head process can talk to the new master id processes, that
!     is processes 4 and 8 we need to create a new group communicator 
!     handle: myhd_comm 
! 
!     Likewise, each newly promoted master node needs to be able to 
!     communicate with it's subprocesses after 
! 
!     master node id:   0     | 4     | 8
!     subprocess id :   1 2 3 | 5 6 7 | 9 10 11
! 
!     So we define a new communicator myslv_comm for this.  Lastly, 
!     for mpi reasons we need to define the communicator ranks, this 
!     is why we not only create the groups, but call mpi_comm_rank.    
! 
!     - BB March 2011
!
!     INPUT      MEANING 
!     -----      ------- 
!     MASTER     head process from driver routine 
!     MYID       global id numbers 
!     NGROUPS    number of groups to break processes into 
!     NPROCS     number of processes 
!  
!     OUTPUT     MEANING 
!     ------     ------- 
!     IROW       process group - 1  
!     MASCOM     mass communication for mumps 
!     MYNID      new id (equivalent to jcol) 
!     MYHD_COMM  head node communication handle
!     MYSLV_COMM slave communication handle 
!
      INCLUDE 'mpif.h'
      ALLOCATABLE IHDGRP_ID(:), ISLVGRP_ID(:) 
      INTEGER*4 IHDGRP_ID, ISLVGRP_ID
! 
!----------------------------------------------------------------------!
! 
!.... create teams
      CALL MPI_BCAST(NGROUPS,1,MPI_INTEGER,MASTER,MYCOMM,MPIERR)
      IF (NPROCS.EQ.1) THEN
         IROW = 0
         MYNID = MYID 
         MYHD_COMM = MYCOMM
         MYSLV_COMM = MYCOMM 
         RETURN
      ENDIF 
      NROWS = NPROCS/NGROUPS
      IROW = MYID/NROWS
      JCOL = MOD(MYID,NROWS)
      MYNID = JCOL 
! 
!.... create the head node group
      ALLOCATE(IHDGRP_ID(NGROUPS))
      DO 1 I=1,NGROUPS 
         IHDGRP_ID(I) = (I - 1)*NROWS 
    1 CONTINUE
      NGP = NGROUPS 
      CALL MPI_COMM_GROUP(MYCOMM,MYHD_GRP,MPIERR)
      CALL MPI_GROUP_INCL(MYHD_GRP,NGP,IHDGRP_ID,MYHD_ROW,MPIERR)
      CALL MPI_COMM_CREATE(MYCOMM,MYHD_ROW,MYHD_COMM,MPIERR)
      IF (MYNID.EQ.MASTER) THEN
         CALL MPI_COMM_RANK(MYHD_COMM,MYHD_RANK,MPIERR)
         CALL MPI_COMM_SIZE(MYHD_COMM,NHEADS,MPIERR) 
      ENDIF 
! 
!.... now create the slave groups
      ALLOCATE(ISLVGRP_ID(NROWS))
      DO 2 I=1,NROWS
         ISLVGRP_ID(I) = IROW*NROWS + I - 1 
    2 CONTINUE 
      NRP = NROWS
      CALL MPI_COMM_GROUP(MYCOMM,MYSLV_GRP,MPIERR)
      CALL MPI_GROUP_INCL(MYSLV_GRP,NRP,ISLVGRP_ID,MYSLV_ROW,MPIERR)
      CALL MPI_COMM_CREATE(MYCOMM,MYSLV_ROW,MYSLV_COMM,MPIERR)
! 
!.... and the group info 
      DO 3 I=1,NROWS
         IF (I-1.EQ.IROW) THEN
            CALL MPI_COMM_RANK(MYSLV_COMM,MYSLV_RANK,MPIERR)
            CALL MPI_COMM_SIZE(MYSLV_COMM,NSUBS,MPIERR) 
         ENDIF
    3 CONTINUE
      DEALLOCATE(IHDGRP_ID)
      DEALLOCATE(ISLVGRP_ID) 
! 
!.... block until we are done
      CALL MPI_BARRIER(MYCOMM,MPIERR) 
      RETURN 
      END 
!                                                                      !
!======================================================================!
!                                                                      !
      SUBROUTINE mkMATgroup(MYCOMM,MYID,NPARTS, NDMAT,MYMAT_COMM) 
! 
!     Makes the distributed matrix MPI groups.  The approach is simple
!     we choose ny/2 as an upper bound, so if the number of processes 
!     in a group is greater than ny/2 then we exclude group ids 
!     greater than ny/2 + 1 and simply leave those extra processes 
!     for MUMPS to use as it sees fit.  
!  
      INCLUDE 'mpif.h'
      INTEGER*4 MYCOMM,MYID,NPARTS
      INTEGER*4 NDMAT,MYMAT_COMM
      INTEGER*4, ALLOCATABLE :: IMAT_ID(:)  
      INTEGER*4 MYMAT_GRP,I,NPR,MPIERR
! 
!----------------------------------------------------------------------!
!
      ALLOCATE(IMAT_ID(NPARTS)) 
      DO 1 I=0,NPARTS-1
        IMAT_ID(I+1) = I 
    1 CONTINUE 
      NPR = NPARTS
      CALL MPI_COMM_GROUP(MYCOMM,MYMAT_GRP,MPIERR) 
      CALL MPI_GROUP_INCL(MYMAT_GRP,NPR,IMAT_ID,MYMAT_ROW,MPIERR)
      CALL MPI_COMM_CREATE(MYCOMM,MYMAT_ROW,MYMAT_COMM,MPIERR)
      IF (MYID.LE.NPARTS-1) THEN 
         CALL MPI_COMM_RANK(MYMAT_COMM,MYMAT_RANK,MPIERR)
         CALL MPI_COMM_SIZE(MYMAT_COMM,NDMAT,MPIERR)
      ENDIF
      CALL MPI_BARRIER(MYCOMM,MPIERR)
      RETURN
      END
