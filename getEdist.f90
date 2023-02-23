
MODULE VARS
   IMPLICIT NONE
   INTEGER, PARAMETER  :: REAL64 = SELECTED_REAL_KIND(15, 307)  
   INTEGER :: NNODES = 0
   INTEGER :: NEDGES = 0
   INTEGER, ALLOCATABLE :: EDGES(:,:)
   REAL(KIND=REAL64), ALLOCATABLE :: WEIGHTS(:)
END MODULE VARS

PROGRAM EDIST
   USE VARS
   IMPLICIT NONE
   INTEGER :: NMIN
   INTEGER :: NTS
   REAL(KIND = REAL64), ALLOCATABLE :: EMIN(:), ETS(:), DISTS(:)
   INTEGER, ALLOCATABLE :: CONNECTIVITY(:,:)
   INTEGER :: DUNIT
   INTEGER :: I, STARTMIN
   CHARACTER(LEN=8) :: A

   CALL GETARG(1,A)
   READ(A, *) NMIN
   CALL GETARG(2,A)
   READ(A, *) NTS
   ! try whether an additional argument is provided - if we do start the iteration on a later minimum
   IF (IARGC().GT.2) THEN
      CALL GETARG(3,A) 
      READ(A, *) STARTMIN
   ELSE
      STARTMIN = 1
   END IF
   ! get number of minima and ts from file length
  ! NMIN = FILE_LENGTH("min.data")
  ! NTS = FILE_LENGTH("ts.data")

   ! allocate data arrays
   ALLOCATE(EMIN(NMIN))
   ALLOCATE(ETS(NTS))
   ALLOCATE(CONNECTIVITY(NTS,2))

   ! parse data files
   CALL PARSE_MIN_DATA("min.data",NMIN,EMIN)
   CALL PARSE_TS_DATA("ts.data",NTS,ETS,CONNECTIVITY)

   WRITE(*,*) " Parsed min.data and ts.data files"
   WRITE(*,*) " Found ", NMIN, " minima and ", NTS, " transition states"
   WRITE(*,*) " "   

   ! create graph from the data
   CALL CREATE_GRAPH(NMIN,NTS,EMIN,ETS,CONNECTIVITY)

   WRITE(*,*) " Created graph with ", NNODES, " nodes and ", NEDGES, " edges"
   WRITE(*,*) " " 

   ! check connectivity of graph
!   CALL CHECK_CONNECTIVITY

   ! get paths
   ALLOCATE(DISTS(NNODES))
   DO I=STARTMIN,NNODES
      DISTS(1:NNODES) = 0.0
      WRITE(*,*) " Starting Dijkstra path searches from node ", I     
      CALL DIJKSTRA(I,DISTS)
      CALL FILE_OPEN("distances.dat",DUNIT,.TRUE.)
      WRITE(DUNIT,*) DISTS
      CLOSE(DUNIT)
      WRITE(*,*) " Completed Dijkstra calulations for this minimum"
      WRITE(*,*) " "  
   END DO

   ! deallocate arrays in main program
   DEALLOCATE(DISTS, CONNECTIVITY, EMIN, ETS)
   ! deallocate arrays from VARS
   DEALLOCATE(EDGES, WEIGHTS) 


   CONTAINS
      SUBROUTINE CHECK_CONNECTIVITY()
         USE VARS
         IMPLICIT NONE
         LOGICAL :: CONNECTED(NNODES) 
         INTEGER :: CURRSET(NNODES)              ! which set does each minimum belong to?
         INTEGER :: NEXTSET                      ! what id has the next new set
         INTEGER :: NINSET(NNODES)               ! how many minima are in each set
         INTEGER :: ALLSETS(NNODES,NNODES)       ! All sets, first index is set id, second indes is index for minimum, final entry is given by NINSET(set id)
         INTEGER :: I, K, IDX, JDX, STATI, STATJ
         INTEGER :: NMINMERGE, MINTOADD, MININSET
         INTEGER :: NUNCONNECTED

         WRITE(*,*) "Hello"
         NINSET(:) = 0
         NINSET(1:NNODES) = 0
         ALLSETS(1:NNODES,1:NNODES) = 0
         CONNECTED(1:NNODES) = .FALSE.
         CURRSET(1:NNODES) = 0
         NEXTSET = 1
         
         DO I=1,NEDGES-1,2
            IDX = EDGES(I,1)
            JDX = EDGES(I,2)
            STATI = CURRSET(IDX)
            STATJ = CURRSET(JDX)
            ! base case a new set
            IF ((STATI.EQ.0).AND.(STATJ.EQ.0)) THEN
               CURRSET(IDX) = NEXTSET
               CURRSET(JDX) = NEXTSET
               NINSET(NEXTSET) = 2
               ALLSETS(NEXTSET,1) = IDX
               ALLSETS(NEXTSET,2) = JDX
               NEXTSET = NEXTSET + 1
               IF (CURRSET(IDX).EQ.1) THEN
                  CONNECTED(IDX) = .TRUE.
                  CONNECTED(JDX) = .TRUE.
               END IF
            ! second case - i is in a set, j isn't
            ELSE IF ((STATI.NE.0).AND.(STATJ.EQ.0)) THEN
               NINSET(STATI) = NINSET(STATI) + 1
               ALLSETS(STATI,NINSET(STATI)) = JDX
               CURRSET(JDX) = STATI
               IF (STATI.EQ.1) CONNECTED(JDX) = .TRUE.
            ! third case - j is in a set, i isn't  
            ELSE IF ((STATI.EQ.0).AND.(STATJ.NE.0)) THEN
               NINSET(STATJ) = NINSET(STATJ) + 1
               ALLSETS(STATJ,NINSET(STATJ)) = IDX
               CURRSET(IDX) = STATJ
               IF (STATJ.EQ.1) CONNECTED(IDX) = .TRUE.
            ! final case - both are in sets, if it's the same set we skip, otherwise we merge sets
            ELSE
               IF (STATI.EQ.STATJ) CYCLE
               ! first case I is lower than J (set id that is)
               IF (STATI.LT.STATJ) THEN
                  NMINMERGE = NINSET(STATJ)
                  MININSET = NINSET(STATI)
                  DO K=1,NMINMERGE
                     MINTOADD = ALLSETS(STATJ,K)
                     ALLSETS(STATI,MININSET+K) = MINTOADD
                     CURRSET(MINTOADD) = STATI
                     IF (STATI.EQ.1) CONNECTED(MINTOADD) = .TRUE.
                  END DO
                  NINSET(STATJ) = 0
                  NINSET(STATI) = MININSET + NMINMERGE
               ELSE
                  NMINMERGE = NINSET(STATI)
                  MININSET = NINSET(STATJ)
                  DO K=1,NMINMERGE
                     MINTOADD = ALLSETS(STATI,K)
                     ALLSETS(STATJ,MININSET+K) = MINTOADD
                     CURRSET(MINTOADD) = STATJ
                     IF (STATJ.EQ.1) CONNECTED(MINTOADD) = .TRUE.                     
                  END DO
                  NINSET(STATI) = 0
                  NINSET(STATJ) = MININSET + NMINMERGE
               END IF
            END IF
         END DO

         NUNCONNECTED = 0 
         DO I=1,NNODES
            IF (.NOT.CONNECTED(I)) THEN
               NUNCONNECTED = NUNCONNECTED + 1
               WRITE(*,*) " unconnected minimum: ", I
            END IF 
         END DO
      END SUBROUTINE CHECK_CONNECTIVITY

      SUBROUTINE CREATE_GRAPH(NMIN,NTS,EMIN,ETS,CONNECTIVITY)
         USE VARS
         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NMIN
         INTEGER, INTENT(IN) :: NTS
         REAL(KIND = REAL64), INTENT(IN) :: EMIN(NMIN)
         REAL(KIND = REAL64), INTENT(IN) :: ETS(NTS)
         INTEGER, INTENT(IN) :: CONNECTIVITY(NTS,2)

         REAL(KIND = REAL64), PARAMETER :: EPS = 1.0D-6
         REAL(KIND = REAL64) :: EMIN1, EMIN2, ETHISTS

         LOGICAL :: USED(NMIN)
         REAL(KIND = REAL64) :: TSFORTHISMIN(NMIN)

         INTEGER :: DUMMYEDGES(2*NTS,2)
         REAL(KIND = REAL64) :: DUMMYWEIGHTS(2*NTS)

         INTEGER :: IDX, JDX, I

         !first we can set the number of nodes
         NNODES = NMIN
         NEDGES = 0
         ! initialise edges and weights
         DUMMYEDGES(1:2*NTS,1:2) = 0
         DUMMYWEIGHTS(1:2*NTS) = 0.0D0
         !then we can parse the edges, for this we need to account for three special cases
         ! 1. TS connecting the itself
         ! 2. multiple edges between the same pair
         ! 3. TS energies rougly equal or lower than the minima 
         ! (1) we ignore, (2) we select the lowest energy one, (3) we use an epsilon parameter
         DO IDX = 1,NMIN
            ! (re)set used list and TS energies
            USED(1:NMIN) = .FALSE.
            TSFORTHISMIN(1:NMIN) = 0.0D0

            EMIN1 = EMIN(IDX)
            DO I=1,NTS
               IF (CONNECTIVITY(I,1).EQ.IDX) THEN
                  JDX = CONNECTIVITY(I,2)
                  IF (IDX.EQ.JDX) CYCLE
                  ETHISTS = ETS(I)
                  IF (USED(JDX)) THEN
                     IF (TSFORTHISMIN(JDX).GT.ETHISTS) THEN
                        TSFORTHISMIN(JDX) = ETHISTS
                     END IF
                  ELSE
                     USED(JDX) = .TRUE.
                     TSFORTHISMIN(JDX) = ETHISTS
                  END IF                  
               END IF
            END DO

            DO JDX=1,NMIN
               IF (.NOT.USED(JDX)) CYCLE
               EMIN2 = EMIN(JDX)
               ETHISTS = TSFORTHISMIN(JDX)
               IF (EMIN1.GT.ETHISTS) ETHISTS = EMIN1+EPS
               IF (EMIN2.GT.ETHISTS) ETHISTS = EMIN2+EPS
               NEDGES = NEDGES + 1
               DUMMYEDGES(NEDGES,1) = IDX
               DUMMYEDGES(NEDGES,2) = JDX
               DUMMYWEIGHTS(NEDGES) = ETHISTS-EMIN1
               NEDGES = NEDGES + 1
               DUMMYEDGES(NEDGES,2) = IDX
               DUMMYEDGES(NEDGES,1) = JDX
               DUMMYWEIGHTS(NEDGES) = ETHISTS-EMIN2              
            END DO       
         END DO
         ALLOCATE(EDGES(NEDGES,2))
         ALLOCATE(WEIGHTS(NEDGES))
         WEIGHTS(1:NEDGES) = DUMMYWEIGHTS(1:NEDGES)
         EDGES(1:NEDGES,:) = DUMMYEDGES(1:NEDGES,:)
            
      END SUBROUTINE CREATE_GRAPH

      INTEGER FUNCTION FILE_LENGTH(FILE_NAME)
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         INTEGER                       :: FILE_UNIT
         INTEGER                       :: IO_STATUS

         CALL FILE_OPEN(FILE_NAME, FILE_UNIT, .FALSE.)
         FILE_LENGTH = 0
         DO
            READ(FILE_UNIT, *, IOSTAT=IO_STATUS)
            IF (IO_STATUS /= 0) EXIT
            FILE_LENGTH = FILE_LENGTH + 1
         ENDDO
         CLOSE(FILE_UNIT)
      END FUNCTION FILE_LENGTH 

      SUBROUTINE FILE_OPEN(FILE_NAME, FILE_UNIT, APPEND)     
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN)  :: FILE_NAME
         LOGICAL, INTENT(IN)           :: APPEND
         INTEGER, INTENT(OUT)          :: FILE_UNIT
            
         IF (APPEND) THEN
            OPEN(NEWUNIT=FILE_UNIT, FILE=FILE_NAME, POSITION='APPEND')
         ELSE
            OPEN(NEWUNIT=FILE_UNIT, FILE=FILE_NAME)
         END IF   
      END SUBROUTINE FILE_OPEN

      INTEGER FUNCTION GETUNIT()
         IMPLICIT NONE
         LOGICAL :: INUSE
         INTEGER :: UNITNUM

         INUSE=.TRUE.
         UNITNUM=200

         DO WHILE (INUSE)
            INQUIRE(UNIT=UNITNUM,OPENED=INUSE)
            IF (.NOT.INUSE) THEN
               GETUNIT=UNITNUM 
            ELSE     
               UNITNUM=UNITNUM+1
            ENDIF
         ENDDO
      END FUNCTION GETUNIT

      SUBROUTINE READ_LINE(LINE,NWORDS,WORDSOUT)
         CHARACTER(*), INTENT(IN) :: LINE
         INTEGER, INTENT(IN) :: NWORDS
         CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
         INTEGER:: J1,START_IND,END_IND,J2
         CHARACTER(25) :: WORD
         START_IND=0
         END_IND=0
         J1=1
         J2=0
         DO WHILE(J1.LE.LEN(LINE))
             IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
                START_IND=J1
             ENDIF
             IF (START_IND.GT.0) THEN
                IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
                IF (J1.EQ.LEN(LINE)) END_IND=J1
                IF (END_IND.GT.0) THEN
                   J2=J2+1
                   WORD=LINE(START_IND:END_IND)
                   WORDSOUT(J2)=TRIM(WORD)
                   START_IND=0
                   END_IND=0
                ENDIF
             ENDIF
             J1=J1+1
         ENDDO
      END SUBROUTINE READ_LINE

      SUBROUTINE PARSE_MIN_DATA(FNAME,NMIN,EMIN)
         USE VARS, ONLY: REAL64
         IMPLICIT NONE
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER, INTENT(IN) :: NMIN
         REAL(KIND = REAL64), INTENT(OUT) :: EMIN(NMIN)
         INTEGER :: I
         INTEGER :: MUNIT
         CHARACTER(150) :: LINE
         INTEGER , PARAMETER :: NWORDS=6
         CHARACTER(25) :: ENTRIES(NWORDS)=''

         EMIN(1:NMIN) = 0.0D0
         CALL FILE_OPEN(FNAME, MUNIT, .FALSE.)
         DO I=1,NMIN
            READ(MUNIT,'(A)') LINE
            CALL READ_LINE(LINE,NWORDS,ENTRIES)
            READ(ENTRIES(1),*) EMIN(I)
            ! the remaining entries are FRQ, HORDER, IX, IY, IZ
         END DO
         CLOSE(MUNIT) 
      END SUBROUTINE PARSE_MIN_DATA

      SUBROUTINE PARSE_TS_DATA(FNAME,NTS,ETS,CONNECTIVITY)
         USE VARS, ONLY: REAL64
         IMPLICIT NONE  
         CHARACTER(LEN=*), INTENT(IN) :: FNAME
         INTEGER, INTENT(IN) :: NTS
         REAL(KIND = REAL64), INTENT(OUT) :: ETS(NTS)
         INTEGER, INTENT(OUT) :: CONNECTIVITY(NTS,2)
         INTEGER :: I
         INTEGER :: IDX, JDX
         INTEGER :: TUNIT
         CHARACTER(150) :: LINE
         INTEGER , PARAMETER :: NWORDS=8
         CHARACTER(25) :: ENTRIES(NWORDS)=''

         ETS(1:NTS) = 0.0D0
         CALL FILE_OPEN(FNAME,TUNIT,.FALSE.)
         DO I=1,NTS
            READ(TUNIT,'(A)') LINE
            CALL READ_LINE(LINE,NWORDS,ENTRIES)
            READ(ENTRIES(1),*) ETS(I)
            READ(ENTRIES(4),*) IDX
            READ(ENTRIES(5),*) JDX
            ! the remaining entries are FRQ (2), HORDER (3), IX (6), IY (7), IZ (8)           
            IF (IDX.LT.JDX) THEN
               CONNECTIVITY(I,1) = IDX
               CONNECTIVITY(I,2) = JDX
            ELSE
               CONNECTIVITY(I,1) = JDX
               CONNECTIVITY(I,2) = IDX
            END IF 
         END DO
         CLOSE(TUNIT) 
      END SUBROUTINE PARSE_TS_DATA      

      SUBROUTINE DIJKSTRA(SRC,DISTS)
         USE VARS
         IMPLICIT NONE
         ! Variables
         INTEGER, INTENT(IN) :: SRC
         REAL(KIND=REAL64), INTENT(OUT)   :: DISTS(NNODES)                ! Shortest path lengths
         ! Local variables
         INTEGER :: VISITED(NNODES)   ! Logical variable to track which nodes we have visited
         INTEGER :: I, J, K, Q, THISNODE, IDXNODE
         LOGICAL :: FOUND
         REAL(KIND=REAL64), PARAMETER :: ELARGE = 1.0D10
         REAL(KIND=REAL64) :: NEWDIST, THISDIST
        
         ! Initialize distances and visited status
         DO I = 1, NNODES
            DISTS(I) = HUGE(0.0D0)
            VISITED(I) = 0
         END DO

         ! Initialize start node
         DISTS(SRC) = 0.0
         VISITED(SRC) = 1
         DO I = 1, NEDGES
            IF (EDGES(I,1) .EQ. SRC) THEN
               DISTS(EDGES(I,2)) = WEIGHTS(I)
            END IF
         END DO

         ! Dijkstra's algorithm
         DO J = 1, NNODES-1
            FOUND = .FALSE.
            DO K = 1, NNODES
               !need to check whether we have visited this node yet and if the distance is smaller than HUGE - then we want check the distance from this node
               IF ((VISITED(K) .EQ. 0).AND.(DISTS(K).LT.ELARGE)) THEN
                  FOUND = .TRUE.
                  THISDIST = DISTS(K)
                  THISNODE = K
                  EXIT
               END IF
            END DO
            IF (.NOT. FOUND) EXIT
            VISITED(THISNODE) = 1
            DO Q = 1, NEDGES
               IF (EDGES(Q,1) .EQ. THISNODE) THEN
                  NEWDIST = WEIGHTS(Q) + THISDIST
                  IDXNODE = EDGES(Q,2)
                  IF (DISTS(IDXNODE) .GT. NEWDIST) THEN
                     DISTS(IDXNODE) = NEWDIST
                  END IF
               END IF
            END DO
         END DO
      END SUBROUTINE DIJKSTRA

      
END PROGRAM EDIST
