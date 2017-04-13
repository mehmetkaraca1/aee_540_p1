PROGRAM SPHERE
!
! **********************************************************************
!   Ce programme determine toutes les spheres centrees a l'origine d'un
!   repere cartesien et contenant au moins un point de coordonnees
!   entieres dans ce repere, c est a dire tel que:
!                  n1**2+n2**2+n3**2=Cste=N     N entier
!   Chaque axe est discretise en Nx points indice de -Nx/2+1 a Nx/2-1.
!   En fait, on travaille sur 1/8 d espace correspondant aux trois
!   indices positifs. On compte, pour chaque valeur N=Rayon**2 contenant
!   au moins un point, le nombre de points trouve sur 1/8 de sphere et
!   on complete par symetrie. On ne garde en memoire que les spheres
!   contenant plus d un nombre donne de points.
! **********************************************************************
!
      IMPLICIT NONE
!
      INTEGER, PARAMETER :: NxMax=512, NRMAX=55000, NTRIPL_MAX=500
!
      INTEGER(KIND=4), DIMENSION(0:NRMAX,NTRIPL_MAX)         :: NN1, NN2, NN3, NGEN
      INTEGER(KIND=4)                                        :: N, M, L, N1, N2, N3, KK, KC2, NRAY, &
                                                                TEMP, TEMP0, Nx, NRG
      INTEGER(KIND=4), DIMENSION(NTRIPL_MAX)                 :: TEMP1, TEMP2, TEMP3
      INTEGER(KIND=4), DIMENSION(0:NRMAX)                    :: KSPH2, NTRIPL, NPOINTS, NPOINTS_8,  &
                                                                NTRG_8, NTRG
      REAL(KIND=8), DIMENSION(0:NRMAX)                       :: KG
      LOGICAL                                                :: TEST
!
      WRITE(6,*) 'Nx='  ; READ(5,*) Nx
      KC2=(Nx/2-1)**2
!
! -------------------------------------------- Determination de toutes -
!                                                   les spheres K=cste
!
      NRAY=0
      KSPH2(0)=0
      NTRIPL(0)=0
!
      DO N3=0,Nx/2-1
         DO N2=0,N3
            DO N1=0,N2
!
               KK=N1**2+N2**2+N3**2                ! voici un rayon possible
!
               IF ( KK <= KC2 ) THEN               ! est-il plus petit que Kc ?
!
! Test pour savoir si le rayon existe deja dans la liste
!
                  TEST=.TRUE.
                  DO N=0,NRAY
                     IF( KK == KSPH2(N)) THEN      ! ce rayon existe deja
                       TEST=.FALSE.
                       NTRIPL(N)=NTRIPL(N)+1       ! voici un nouveau triplet pour ce rayon
                       NN1(N,NTRIPL(N))=N1
                       NN2(N,NTRIPL(N))=N2         ! on stocke et on numerote les valeurs (N1,N2,N3)
                       NN3(N,NTRIPL(N))=N3         ! pour ce rayon
                     ENDIF
                  ENDDO
!
! s il n a pas encore trouve, ...
                  IF(TEST) THEN
                    NRAY=NRAY+1                 ! ca en fait un de plus
                    KSPH2(NRAY)=KK
                    NTRIPL(NRAY)=1              ! c est le 1er triplet (N1,N2,N3) pour ce rayon
                    NN1(NRAY,1)=N1              ! on stocke les valeurs
                    NN2(NRAY,1)=N2              ! de N1, N2, N3 de ce rayon, numerotees 1
                    NN3(NRAY,1)=N3
                    IF(MOD(NRAY,100) == 0) WRITE(6,*) 'NRAY=',NRAY
                  ENDIF
!
               ENDIF
!
            ENDDO
         ENDDO
      ENDDO
!
      WRITE(6,*) ' Classement par ordre croissant'
!
      DO N=1,NRAY
!
            IF(KSPH2(N) < KSPH2(N-1)) THEN   ! celui la n est pas a sa place
!
              DO M=1,N
                 IF( (KSPH2(N) > KSPH2(M-1)).AND. (KSPH2(N) < KSPH2(M)) ) THEN
! Stocke les caracteristisques du rayon M
                   TEMP=KSPH2(M)
                   TEMP0=NTRIPL(M)
                   DO L=1,NTRIPL(M)
                      TEMP1(L)=NN1(M,L)
                      TEMP2(L)=NN2(M,L)
                      TEMP3(L)=NN3(M,L)
                   ENDDO
! Ecrase le rayon M par le rayon N:
                   KSPH2(M)=KSPH2(N)
                   NTRIPL(M)=NTRIPL(N)
                   DO L=1,NTRIPL(N)
                      NN1(M,L)=NN1(N,L)
                      NN2(M,L)=NN2(N,L)
                      NN3(M,L)=NN3(N,L)
                   ENDDO
! Attribue au rayon N les caracterisitiques du rayon M
                   KSPH2(N)=TEMP
                   NTRIPL(N)=TEMP0
                   DO L=1,NTRIPL(N)
                      NN1(N,L)=TEMP1(L)
                      NN2(N,L)=TEMP2(L)
                      NN3(N,L)=TEMP3(L)
                   ENDDO
                 ENDIF
              ENDDO
!
            ENDIF
      ENDDO
!
! ------------------------- Calcul du nombre de points sur les spheres -
!
      WRITE(6,*) ' Calcul du nombre de points sur les spheres'
!
! * sur le 1er octant:
!
      DO N=0,NRAY
         NPOINTS_8(N)=NTRIPL(N)
         DO L=1,NTRIPL(N)
            N1=NN1(N,L)
            N2=NN2(N,L)
            N3=NN3(N,L)
            IF( (N1 /= N2).AND.(N2 /= N3).AND.(N1 /= N3) ) THEN
               NPOINTS_8(N)=NPOINTS_8(N)+5
               NGEN(N,L)=6
            ELSE
               IF( (N1 == N2).AND.(N2 == N3).AND.(N1 == N3) ) THEN
                  NPOINTS_8(N)=NPOINTS_8(N)
                  NGEN(N,L)=1
               ELSE
                  NPOINTS_8(N)=NPOINTS_8(N)+2
                  NGEN(N,L)=3
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
! * sur la sphere complete:
!
      DO N=0,NRAY
         NPOINTS(N)=0
         DO L=1,NTRIPL(N)
            N1=NN1(N,L)
            N2=NN2(N,L)
            N3=NN3(N,L)
            IF( (N1 /= 0).AND.(N2 /= 0).AND.(N3 /= 0)) THEN         ! pas de zeros
               NPOINTS(N)=NPOINTS(N)+NGEN(N,L)*8
            ELSE
               IF( ((N1 == 0).AND.(N2 /= 0).AND.(N3 /= 0)).OR.          &
                   ((N2 == 0).AND.(N1 /= 0).AND.(N3 /= 0)).OR.          &
                   ((N3 == 0).AND.(N1 /= 0).AND.(N2 /= 0)) ) THEN   ! 1 zero
                   NPOINTS(N)=NPOINTS(N)+NGEN(N,L)*4
               ELSE
                   IF( (N1 == 0).AND.(N2 == 0).AND.(N3 == 0)) THEN  ! 3 zeros
                       NPOINTS(N)=NPOINTS(N)+NGEN(N,L)*1
                   ELSE
                       NPOINTS(N)=NPOINTS(N)+NGEN(N,L)*2            ! 2 zeros
                   ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
!
      WRITE(6,*) 'Nombre total de spheres : ', NRAY+1
      WRITE(6,*) 'Nombre de triplets maxi : ', MAXVAL(NTRIPL)
      WRITE(6,*) '   sur la sphere numero : ', MAXLOC(NTRIPL)-1
      WRITE(6,*) '               de rayon : ', SQRT(REAL(KSPH2(MAXLOC(NTRIPL)-1), KIND=8))
      WRITE(6,*)
      READ(5,*)
!
      WRITE(6,*) ' Selection des spheres qui ont le plus de points'
!
      OPEN(UNIT=2,FILE='SPHERE.DAT',STATUS='UNKNOWN')
!
      WRITE(2,*)
!
! -------------------------------- Selection de celles qui contiennent -
!                                                    le plus de points
! * Sphere K=1:
! ----------------
      NRG=1
      KG(1)=1.D00
      NTRG_8(1)=NPOINTS_8(1)
      NTRG(1)=NPOINTS(1)
      WRITE(2,*) KG(1), NTRG_8(1)
      WRITE(2,*) NN1(1,1),NN2(1,1),NN3(1,1)
      WRITE(2,*) NN3(1,1),NN1(1,1),NN2(1,1)
      WRITE(2,*) NN2(1,1),NN3(1,1),NN1(1,1)
!
      DO N=1,NRAY
!
! * Spheres 1<K<=3:
! -----------------
         IF((KSPH2(N) > 1).AND.(KSPH2(N) <= 9)) THEN
           IF(NPOINTS(N) >= 23) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 3<K<=5:
! -----------------
         IF((KSPH2(N) > 9).AND.(KSPH2(N) <= 25)) THEN
           IF(NPOINTS(N) >= 40) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 5<K<=7:
! -----------------
         IF((KSPH2(N) > 25).AND.(KSPH2(N) <= 49)) THEN
           IF(NPOINTS(N) >= 70) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 7<K<=9:
! -----------------
         IF((KSPH2(N) > 49).AND.(KSPH2(N) <= 81)) THEN
           IF(NPOINTS(N) >= 100) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 9<K<=11:
! -----------------
         IF((KSPH2(N) > 81).AND.(KSPH2(N) <= 121)) THEN
           IF(NPOINTS(N) >= 120) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 11<K<=13:
! -----------------
         IF((KSPH2(N) > 121).AND.(KSPH2(N) <= 169)) THEN
           IF(NPOINTS(N) >= 160) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 13<K<=15:
! -----------------
         IF((KSPH2(N) > 169).AND.(KSPH2(N) <= 225)) THEN
           IF(NPOINTS(N) >= 200) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 15<K<=20:
! -----------------
         IF((KSPH2(N) > 225).AND.(KSPH2(N) <= 400)) THEN
           IF(NPOINTS(N) >= 240) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 20<K<=30:
! -----------------
         IF((KSPH2(N) > 400).AND.(KSPH2(N) <= 900)) THEN
           IF(NPOINTS(N) >= 360) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 30<K<=40:
! -----------------
         IF((KSPH2(N) > 900).AND.(KSPH2(N) <= 1600)) THEN
           IF(NPOINTS(N) >= 490) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 40<K<=50:
! -----------------
         IF((KSPH2(N) > 1600).AND.(KSPH2(N) <= 2500)) THEN
           IF(NPOINTS(N) >= 500) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 50<K<=60:
! -----------------
         IF((KSPH2(N) > 2500).AND.(KSPH2(N) <= 3600)) THEN
           IF(NPOINTS(N) >= 600) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
! * Spheres 60<K<=70:
! -----------------
         IF((KSPH2(N) > 3600).AND.(KSPH2(N) <= 4900)) THEN
           IF(NPOINTS(N) >= 700) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 70<K<=80:
!! -----------------
         IF((KSPH2(N) > 4900).AND.(KSPH2(N) <= 6400)) THEN
           IF(NPOINTS(N) >= 800) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                   WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 80<K<=90:
!! -----------------
         IF((KSPH2(N) > 6400).AND.(KSPH2(N) <= 8100)) THEN
           IF(NPOINTS(N) >= 900) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 90<K<=100:
!! -----------------
         IF((KSPH2(N) > 8100).AND.(KSPH2(N) <= 10000)) THEN
           IF(NPOINTS(N) >= 1000) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 100<K<=110:
!! -----------------
         IF((KSPH2(N) > 10000).AND.(KSPH2(N) <= 12100)) THEN
           IF(NPOINTS(N) >= 1500) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 110<K<=150:
!! ---------------------
         IF((KSPH2(N) > 12100).AND.(KSPH2(N) <= 22500)) THEN
           IF(NPOINTS(N) >= 2000) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 150<K<=200:
!! ---------------------
         IF((KSPH2(N) > 22500).AND.(KSPH2(N) <= 40000)) THEN
           IF(NPOINTS(N) >= 3000) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!!
!! * Spheres 200<K:
!! ----------------
         IF(KSPH2(N) > 40000) THEN
           IF(NPOINTS(N) >= 4000) THEN
             NRG=NRG+1
             KG(NRG)=SQRT(REAL(KSPH2(N), KIND=8))
             NTRG_8(NRG)=NPOINTS_8(N)
             NTRG(NRG)=NPOINTS(N)
             WRITE(2,*) KG(NRG), NTRG_8(NRG)
             DO L=1,NTRIPL(N)
                IF(NGEN(N,L) == 1) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                ENDIF
                IF(NGEN(N,L) == 3) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                ENDIF
                IF(NGEN(N,L) == 6) THEN
                  WRITE(2,*) NN1(N,L),NN2(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN1(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN3(N,L),NN1(N,L)
                  WRITE(2,*) NN1(N,L),NN3(N,L),NN2(N,L)
                  WRITE(2,*) NN2(N,L),NN1(N,L),NN3(N,L)
                  WRITE(2,*) NN3(N,L),NN2(N,L),NN1(N,L)
                ENDIF
             ENDDO
           ENDIF
         ENDIF
!
      ENDDO
      CLOSE(UNIT=2)
!
      OPEN(UNIT=1,FILE='SPHERE.RES',STATUS='UNKNOWN')
!
      WRITE(1,*) ' Nx                       : ', Nx
      WRITE(1,*) ' K max                    : ', SQRT(REAL(KC2, KIND=8))
      WRITE(1,*) ' Nombre de valeurs K=Cste : ', NRAY+1
      WRITE(1,*) ' Nombre de triplets maxi  : ', MAXVAL(NTRIPL)
      WRITE(1,*) '    sur la sphere numero  : ', MAXLOC(NTRIPL)-1
      WRITE(1,*) '                de rayon  : ', SQRT(REAL(KSPH2(MAXLOC(NTRIPL)-1), KIND=8))
      WRITE(1,*)
      WRITE(1,*) '************************************************'
      WRITE(1,*) '    Nombre retenu de valeurs K=Cste : ',NRG
      WRITE(1,*) '    Nombre de triplets maxi         : ',MAXVAL(NTRG_8)
      WRITE(1,*) '************************************************'
      WRITE(1,*)
      WRITE(1,*) '------------------------------------------------------'
      WRITE(1,*) ' Num  *   K    *  K**2  * Nb pts/Octant * Nb pts total'
      WRITE(1,*) '------------------------------------------------------'
!
      DO N=1,NRG
         WRITE(1,100) N,KG(N),KG(N)**2,NTRG_8(N),NTRG(N)
      ENDDO
 100  FORMAT(1X,I4,2X,F10.6,2X,F10.2,2X,I4,2X,I6)
!
      CLOSE(UNIT=1)
      !
      OPEN(UNIT=3,FILE='SPHERE.bin',STATUS='UNKNOWN')
      write(3,*) KG(1:NRG)
      close(3)
      WRITE(6,*)
      WRITE(6,*) '************************************************'
      WRITE(6,*) '    Nombre retenu de valeurs K=Cste : ',NRG
      WRITE(6,*) '    Nombre de triplets maxi         : ',MAXVAL(NTRG_8)
      WRITE(6,*) '************************************************'
      WRITE(6,*)
      WRITE(6,*) ' Ces valeurs servent a dimensionner les tableaux: '
      WRITE(6,*) '       KSPH, NTRIPL, N1SPH, N2SPH, N3SPH'
      WRITE(6,*) '  des programmes data_DNS.f90 et solveur THI '
      WRITE(6,*) 'Il faut ouvrir "a la main" le fichier SPHERE.DAT'
      WRITE(6,*) 'et recopier la valeur ',NRG,' sur la 1ere ligne'
      WRITE(6,*) 'qui est blanche. Le reste suit tout seul...'
      


      !
      STOP
      END
