!
! **********************************************************************
!
!                  PROGRAMME DE DEFINITION DES DONNEES
!                POUR LE CALCUL DE L'AUTO-AMORTISSEMENT
!                 D'UNE TURBULENCE HOMOGENE ET ISTROPE
!                      DANS UN DOMAINE CUBIQUE
!
!     Generation des Fichiers: -    THI.DAT : Parametres du calcul
!                                                et spheres K=Cste
!                              - CHINIT.DAT : Champ Spectral initial
!     -- Ivan Fedioun        
! **********************************************************************
!
!  update : Program which generates Homogeneous Isotropic Turbulent  
!                   Initial Field for Turbulence Decay Simulation.
!
!   
! Declaration des modules
!
!
      MODULE PARAM
        IMPLICIT NONE
        INTEGER(kind=4)               :: NTOTRAY, NTRMAX
        REAL   (kind=8), PARAMETER    :: PI=3.14159265358979D00, CK=1.4D00,Lx=0.014d00
        COMPLEX(kind=8), PARAMETER    :: IOTA=(0.D00, 1.D00)
      END MODULE PARAM
!
      MODULE SPHERES
        IMPLICIT NONE
        REAL   (kind=8), DIMENSION(1549)     :: KSPH
        INTEGER(kind=4), DIMENSION(1549)     :: NTRIPL
        INTEGER(kind=4), DIMENSION(1549,259) :: N1SPH, N2SPH, N3SPH
      END MODULE SPHERES
!
      MODULE ADIM
        IMPLICIT NONE
        REAL(kind=8)   :: KC,KI,REU0L
      END MODULE ADIM
!
!
      MODULE CONSTANTES
         IMPLICIT NONE
         INTEGER(kind=4), PARAMETER    :: Nx=128, Ny=Nx, Nz=Nx
         REAL   (kind=8)               :: PI=3.14159265358979D00
         COMPLEX(kind=8)               :: IOTA=(0.D00,1.D00)
      END MODULE CONSTANTES
!
!
      MODULE TEMP_FFT
         IMPLICIT NONE
         INTEGER(kind=4),                  SAVE :: N1N2N3
         REAL   (kind=8),                  SAVE :: N1N2N3M1
         INTEGER(kind=4), DIMENSION(:), POINTER :: N22, N33
         INTEGER(kind=4), DIMENSION(3),    SAVE :: Nxyz
      END MODULE TEMP_FFT
!
!
! -------------------------------------------------------------------------
!                           Programme principal
! -------------------------------------------------------------------------
!
      PROGRAM DATA_THI
!
      USE PARAM
      USE SPHERES
      USE ADIM
!
      IMPLICIT NONE
      INTEGER(kind=4)                                :: I,Nx,NTOUR,NPT,NFREQ,NRAY, &
                                                        N,NFOIS,N1,N2,N3
      INTEGER(kind=4), DIMENSION(3)                  :: Mx
      CHARACTER                                      :: REP*1
      REAL(kind=8)                                   :: DT
      COMPLEX(kind=8), DIMENSION(:,:,:), ALLOCATABLE :: UF1, UF2, UF3
      real(kind=8), DIMENSION(:,:,:), ALLOCATABLE    ::Du1D1,Du2D1,Du3D1,U1,U2,U3
!
! ---------------------------------------- Donnees du calcul numerique -
!
      WRITE(6,*) ' Nx                               : ' ; READ(5,*) Nx
      WRITE(6,*) ' Dt (fraction de L.E.turnover)    : ' ; READ(5,*) DT
      WRITE(6,*) ' Nombre de retournements calcules : ' ; READ(5,*) NTOUR
      WRITE(6,*) ' Frequence des stats./resultat    : ' ; READ(5,*) NFREQ
      WRITE(6,*) ' Reynolds (Uo,Lx)                 : ' ; READ(5,*) REU0L
!
      Mx(:)=Nx
      KC=DFLOAT(Nx/2-1)
      NPT=INT(FLOAT(NTOUR)/DT)
!
! ----------------------------------------------------- Spheres K=Cste -
!
      OPEN(UNIT=2,FILE='SPHERE.DAT')
      READ(2,*) NTOTRAY
      DO NRAY=1,NTOTRAY
         READ(2,*) KSPH(NRAY),NTRIPL(NRAY)
         DO N=1,NTRIPL(NRAY)
            READ(2,*) N1SPH(NRAY,N),N2SPH(NRAY,N),N3SPH(NRAY,N)
         ENDDO
      ENDDO
      CLOSE(UNIT=2)
      WRITE(6,*) ' Ok lecture SPHERE.DAT'
      READ(5,*)
!
! ------------------------------------------------------ Champ initial -
      OPEN(UNIT=2, FILE='COMPARE.SPE',STATUS='UNKNOWN')
      OPEN(UNIT=10,FILE='PDF.RES',    STATUS='UNKNOWN')
!
      ALLOCATE(UF1(0:Nx/2-1,-Nx/2:Nx/2-1,-Nx/2:Nx/2-1)) ;WRITE(6,*) 'U1F:',ALLOCATED(UF1)
      ALLOCATE(UF2(0:Nx/2-1,-Nx/2:Nx/2-1,-Nx/2:Nx/2-1)) ;WRITE(6,*) 'U2F:',ALLOCATED(UF2)
      ALLOCATE(UF3(0:Nx/2-1,-Nx/2:Nx/2-1,-Nx/2:Nx/2-1)) ;WRITE(6,*) 'U3F:',ALLOCATED(UF3)
      ALLOCATE(U1(0:Nx-1,0:Nx-1,0:Nx-1)) ;WRITE(6,*) 'U1:',ALLOCATED(U1)
      ALLOCATE(U2(0:Nx-1,0:Nx-1,0:Nx-1)) ;WRITE(6,*) 'U2:',ALLOCATED(U1)
      ALLOCATE(U3(0:Nx-1,0:Nx-1,0:Nx-1)) ;WRITE(6,*) 'U3:',ALLOCATED(U1)
      
!      READ(5,*)
!
      NFOIS=0
      REP='o'
!
      DO WHILE(REP == 'o')
!
         NFOIS=NFOIS+1
!
         CALL INITALEA(NFOIS,Nx,UF1,UF2,UF3)
!
         WRITE(6,*) 'Re-initialisation ? (o/n) : ' ; READ(5,'(A1)') REP
!
      ENDDO
!
      WRITE(6,*) 'Champ initial construit'
!
      CLOSE(UNIT=2)
      CLOSE(UNIT=10)
! -------------------------------------Ecriture des fichiers de sortie -
!
! * Fichier calcul:
!
      OPEN(UNIT=1,FILE='THI_DNS.DAT',STATUS='UNKNOWN')
        WRITE(1,*) NPT,DT,NFREQ
        WRITE(1,*) REU0L,Lx
        WRITE(1,*) KC,KI
        WRITE(1,*) NTOTRAY
        DO NRAY=1,NTOTRAY
           WRITE(1,*) KSPH(NRAY),NTRIPL(NRAY)
           DO N=1,NTRIPL(NRAY)
              WRITE(1,*) N1SPH(NRAY,N),N2SPH(NRAY,N),N3SPH(NRAY,N)
           ENDDO
        ENDDO
      CLOSE(UNIT=1)
!
! * Champ initial:
!
      OPEN(UNIT=1,FILE='CHINIT.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(1) UF1
          WRITE(1) UF2
          WRITE(1) UF3
      CLOSE(UNIT=1)
      CALL TFD3DR(U1,UF1,1) 
      CALL TFD3DR(U2,UF2,1) 
      CALL TFD3DR(U3,UF3,1) 
      OPEN(UNIT=2,FILE='CHINITPHY.DAT',STATUS='UNKNOWN',FORM='UNFORMATTED')
          WRITE(2) U1
          WRITE(2) U2
          WRITE(2) U3
      CLOSE(UNIT=1)

      !
      END PROGRAM DATA_THI
!
!
! **********************************************************************
!  SPG:          GENERATION D'UN CHAMP INITIAL 3D ALEATOIRE INDIVERGENT
! **********************************************************************
!
      SUBROUTINE INITALEA(NFOIS,Nx,UF1,UF2,UF3)
!
      USE PARAM
      USE SPHERES
      USE ADIM
!
      IMPLICIT NONE
      INTEGER(kind=4)                   :: Nx, I, IK, NK, IRAY, NFOIS, NN, N1, N2, N3, NKi 
      REAL(kind=8)                      :: K, AA, EPIC, EPS, EK, MODPSI,temp,KD,A,B,IntE
      REAL(kind=8), DIMENSION(400)      :: KDATA, EDATA
      REAL(kind=8), DIMENSION(NTOTRAY)  :: EINI
      COMPLEX(kind=8)                   :: PSIF(3),DIVUF
      COMPLEX(kind=8), DIMENSION(0:Nx/2-1,-Nx/2:Nx/2-1,-Nx/2:Nx/2-1) :: UF1, UF2, UF3
!
! -------------------------------------- Donnees pour le champ initial -
!
      WRITE(6,*) 'Champ initial a spectre pointu  '
!!$!
! * Construction d un spectre pointu au mode KI:
      WRITE(6,*) 'Mode excite : ' ; READ(5,*) KI
      AA=2.d0**14.d00/(105.D00*SQRT(PI)*KI**9)
      DO IRAY=1,NTOTRAY
         EINI(IRAY)=AA*KSPH(IRAY)**8*EXP(-4.D00*(KSPH(IRAY)/KI)**2)
      ENDDO
      write(124,*) KSPH
      write(123,*) EINI
      !
! ---------------------------------------- Champ aleatoire indivergent -
!
      CALL RANDOM_SEED(NFOIS)
      DO N1=0,Nx/2-1
         DO N2=-Nx/2,Nx/2-1
            DO N3=-Nx/2,Nx/2-1
!
               NN=N1*N1+N2*N2+N3*N3
               K=SQRT(FLOAT(NN))
!
! * Valeur du spectre au mode K:
!                 EK=AA*K**8*DEXP(-4.D00*(K/KI)**2)
               EK=AA*K**8*EXP(-4.D00*(K/KI)**2)
               
               !
! * Fonction de courant vectorielle spectrale aleatoire:
               IF(NN == 0) THEN
                 MODPSI=0.D00
               ELSE
                 MODPSI=DSQRT(3.D00/2.D00*EK/K**4/(6.D00*PI))
               ENDIF
               DO I=1,3
                  CALL RANDOM_NUMBER(EPS)
                  PSIF(I)=MODPSI*CDEXP(IOTA*2.D00*PI*EPS)
               ENDDO
!
! * Champ de vitesse spectral 3D (troncature isotrope):
               IF(K > KC) THEN
                 UF1(N1,N2,N3)=(0.D00,0.D00)
                 UF2(N1,N2,N3)=(0.D00,0.D00)
                 UF3(N1,N2,N3)=(0.D00,0.D00)
               ELSE
                  
                  UF1(N1,N2,N3)=IOTA*(DBLE(N2)*PSIF(3)-DBLE(N3)*PSIF(2))
                  UF2(N1,N2,N3)=IOTA*(DBLE(N3)*PSIF(1)-DBLE(N1)*PSIF(3))
                  UF3(N1,N2,N3)=IOTA*(DBLE(N1)*PSIF(2)-DBLE(N2)*PSIF(1))
               ENDIF
!
            ENDDO
         ENDDO
      ENDDO
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!$      !!!add by song
!!$      Ki=2
!!$      Kd=50
!!$      A=1
!!$      B=A*KI**(17./3.)
!!$      IntE=1./5.*A*KI**5+(-3./2.*B*KD**(-2./3.)+3./2.*B*KI**(-2./3.))
!!$      DO N1=0,Nx/2-1
!!$         DO N2=-Nx/2,Nx/2-1
!!$            DO N3=-Nx/2,Nx/2-1
!!$!
!!$               NN=N1*N1+N2*N2+N3*N3
!!$               K=SQRT(FLOAT(NN))
!!$!
!!$                  EK .........
!!$! * Fonction de courant vectorielle spectrale aleatoire:
!!$!                 ..............
!!$! * Champ de vitesse spectral 3D (troncature isotrope):
!!$!                  ...................
!!$            ENDDO
!!$         ENDDO
!!$      ENDDO
!!$      Ki=Ki/4.d00
!!$!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! --------------------------------------------- Verification div(UF)=0 -
      DIVUF=(0.D00,0.D00)
      DO N1=-Nx/2+1,Nx/2-1
         DO N2=-Nx/2+1,Nx/2-1
            DO N3=-Nx/2+1,Nx/2-1
               IF(N1 < 0) THEN
                 DIVUF=DIVUF+DCONJG(UF1(-N1,-N2,-N3))*DBLE(N1)          &
                            +DCONJG(UF2(-N1,-N2,-N3))*DBLE(N2)          &
                            +DCONJG(UF3(-N1,-N2,-N3))*DBLE(N3)
               ELSE
                 DIVUF=DIVUF+UF1(N1,N2,N3)*DBLE(N1)                     &
                            +UF2(N1,N2,N3)*DBLE(N2)                     &
                            +UF3(N1,N2,N3)*DBLE(N3)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      WRITE(6,*)
      WRITE(6,*) ' Champ de vitesse genere - div(UF): ',CDABS(DIVUF)
      WRITE(6,*)
      READ(5,*)
!
      END SUBROUTINE INITALEA
!
!
! **************************************************************************
!         INTERFACE DE CALCUL DE TRANSFORMEE DE FOURIER DISCRETE
!                        DE DONNEES REELLES 3D
!                   F(0:Nx(1)-1,0:Nx(2)-1,0:Nx(3)-1)
!                            EN COMPLEXES
!         TF(0:Nx(1)/2-1,-Nx(2)/2:Nx(2)/2-1,-Nx(3)/2:Nx(3)/2-1)
!                      PAR UNE FFT (1:64*64*64)
!
!        LES MODES   (   *   ,    *   ,-Nx(3)/2)  SONT FORCES A ZERO
!                    (   *   ,-Nx(2)/2,    *   )
!
! **************************************************************************
      !
      SUBROUTINE TFD3DR(F,TF,I_SIGN) 
        use, intrinsic ::iso_c_binding
        USE TEMP_FFT
        USE CONSTANTES
        IMPLICIT NONE 
        include 'fftw3.f03' 
        REAL   (kind=8), DIMENSION(0:Nx-1,      0:Ny-1,      0:Nz-1 ) ,INTENT(INOUT):: F
        COMPLEX(kind=8), DIMENSION(0:Nx/2-1,-Ny/2:Ny/2-1,-Nz/2:Nz/2-1),INTENT(INOUT):: TF
        type(C_PTR)                                                                 :: plan
        INTEGER(kind=4)                                                             :: I_SIGN, &
             N, N1, N2, N3 
        COMPLEX(kind=8), DIMENSION(0:Nx/2,0:Ny-1,0:Nz-1)                            :: DATA_C
        real(kind=8)                                                                :: NNNm1
        IF (I_SIGN == -1) THEN 
           NNNm1=1.d0/real(Nx*Ny*Nz,kind=8)
           ! 
           !------------------------------------ Transformee Physique --> Spectral

           !-----------plan setup
           plan = fftw_plan_dft_r2c_3d(Nz,Ny,Nx,F,DATA_C,FFTW_ESTIMATE)

           CALL fftw_execute_dft_r2c(PLAN,F,DATA_C) 


           ! 
           DO N3=0,Nz/2-1 
              DO N2=0,Ny/2-1 
                 DO N1=0,Nx/2-1 
                    TF(N1,N2,N3)=Data_c(N1,N2,N3)*NNNm1
                 ENDDO
              ENDDO
           ENDDO
           DO N3=-Nz/2,-1 
              DO N2=0,Ny/2-1 
                 DO N1=0,Nx/2-1 
                    TF(N1,N2,N3)=Data_c(N1,N2,N3+Nz)*NNNm1
                 ENDDO
              ENDDO
           ENDDO
           DO N3=0,Nz/2-1 
              DO N2=-Ny/2,-1 
                 DO N1=0,Nx/2-1 
                    TF(N1,N2,N3)=Data_c(N1,N2+Ny,N3)*NNNm1
                 ENDDO
              ENDDO
           ENDDO
           DO N3=-Nz/2,-1 
              DO N2=-Ny/2,-1 
                 DO N1=0,Nx/2-1 
                    TF(N1,N2,N3)=Data_c(N1,N2+Ny,N3+Nz)*NNNm1
                 ENDDO
              ENDDO
           ENDDO

           !* Mises a zero des modes -Nx(*)/2: 
           !* Face du 1/2 cube N2=-Nx(2)/2: 
           !TF(:,-Ny/2,:)=(0.D00,0.D00) 
           !* Face du 1/2 cube N3=-Nx(3)/2: 
           !TF(:,:,-Nz/2)=(0.D00,0.D00) 
           !
           call fftw_destroy_plan(plan)
        ELSE
           ! 
           !------------------------------------ Transformee Spectral --> Physique 
           DO N3=0,Nz/2-1 
              DO N2=0,Ny/2-1 
                 DO N1=0,Nx/2-1 
                    Data_C(N1,N2,N3)=TF(N1,N2,N3)
                 ENDDO
              ENDDO
           ENDDO
           DO N3=Nz/2,Nz-1 
              DO N2=0,Ny/2-1 
                 DO N1=0,Nx/2-1 
                    Data_C(N1,N2,N3)=TF(N1,N2,N3-Nz)
                 ENDDO
              ENDDO
           ENDDO
           DO N3=0,Nz/2-1 
              DO N2=Ny/2,Ny-1 
                 DO N1=0,Nx/2-1 
                    Data_C(N1,N2,N3)=TF(N1,N2-Ny,N3)
                 ENDDO
              ENDDO
           ENDDO
           DO N3=Nz/2,Nz-1 
              DO N2=Ny/2,Ny-1 
                 DO N1=0,Nx/2-1 
                    Data_C(N1,N2,N3)=TF(N1,N2-Ny,N3-Nz)
                 ENDDO
              ENDDO
           ENDDO


           Data_C(Nx/2,:,:)=(0.D00,0.D00)
           !write (*,*) maxval(abs(DATA_C(0:Nx/2-1,-Ny/2:Ny/2-1,-Nz/2:Nz/2-1)-TF(:,:,:)))

           plan = fftw_plan_dft_c2r_3d(Nz,Ny,Nx,DATA_C,F,FFTW_ESTIMATE)
           ! 
           CALL fftw_execute_dft_c2r(PLAN,DATA_C,F)
           call fftw_destroy_plan(plan)

        ENDIF
        ! 
      END SUBROUTINE TFD3DR

      
      SUBROUTINE TFD3DRold(F,TF,Nx,NxMax,ISIGN)
!
      INTEGER Nx(3)
      REAL*8 F(0:NxMax-1,0:NxMax-1,0:NxMax-1),DATA(300000)
      COMPLEX*16 TF(0:NxMax/2-1,-NxMax/2:NxMax/2-1,-NxMax/2:NxMax/2-1)
!
      NDIM=3
!
      IF(ISIGN.EQ.-1) THEN
! ----------------------------------- Transormee Physique --> Spectral -
        IFORM=0
        DO I3=0,Nx(3)-1
           DO I2=0,Nx(2)-1
              DO I1=0,Nx(1)-1
                 I=I1+I2*Nx(1)+I3*Nx(1)*Nx(2)
                 DATA(I+1)=F(I1,I2,I3)
              ENDDO
           ENDDO
        ENDDO
        CALL FOUR2(DATA,Nx,NDIM,ISIGN,IFORM)
        DO N1=0,Nx(1)/2-1
           DO N2=-Nx(2)/2+1,Nx(2)/2-1
              N22=N2
              IF (N2.LT.0) N22=N2+Nx(2)
              DO N3=-Nx(3)/2+1,Nx(3)/2-1
                 N33=N3
                 IF(N3.LT.0) N33=N3+Nx(3)
                 N=N1+N22*(Nx(1)/2+1)+N33*(Nx(1)/2+1)*Nx(2)
                 TF(N1,N2,N3)=DCMPLX(DATA(2*N+1),DATA(2*N+2))/DFLOAT(Nx(1)*Nx(2)*Nx(3))
              ENDDO
           ENDDO
        ENDDO
! * Mises a zero des modes -Nx(*)/2:
! * Face du 1/2 cube N2=-Nx(2)/2:
        DO N1=0,Nx(1)/2-1
           DO N3=-Nx(3)/2,Nx(3)/2-1
              TF(N1,-Nx(2)/2,N3)=(0.D00,0.D00)
           ENDDO
        ENDDO
! * Face du 1/2 cube N3=-Nx(3)/2:
        DO N1=0,Nx(1)/2-1
           DO N2=-Nx(2)/2+1,Nx(2)/2-1
              TF(N1,N2,-Nx(3)/2)=(0.D00,0.D00)
           ENDDO
        ENDDO
!
      ELSE
! ---------------------------------- Transformee Spectral --> Physique -
        IFORM=-1
        DO N1=0,Nx(1)/2
           DO N2=-Nx(2)/2,Nx(2)/2
              N22=N2
              IF(N2.LT.0) N22=N2+Nx(2)
              DO N3=-Nx(3)/2,Nx(3)/2
                 N33=N3
                 IF(N3.LT.0) N33=N3+Nx(3)
                 N=N1+N22*(Nx(1)/2+1)+N33*(Nx(1)/2+1)*Nx(2)
                 IF((N1.EQ.Nx(1)/2).OR.(N22.EQ.Nx(2)/2).OR.(N33.EQ.Nx(3)/2)) THEN
                   DATA(2*N+1)=0.D00
                   DATA(2*N+2)=0.D00
                 ELSE
                   DATA(2*N+1)=DREAL(TF(N1,N2,N3))
                   DATA(2*N+2)=DIMAG(TF(N1,N2,N3))
                 ENDIF
              ENDDO
           ENDDO
        ENDDO
        CALL FOUR2(DATA,Nx,NDIM,ISIGN,IFORM)
        DO I1=0,Nx(1)-1
           DO I2=0,Nx(2)-1
              DO I3=0,Nx(3)-1
                 I=I1+I2*Nx(1)+I3*Nx(1)*Nx(2)
                 F(I1,I2,I3)=DATA(I+1)
              ENDDO
           ENDDO
        ENDDO
!
      ENDIF
!
      RETURN
    END SUBROUTINE TFD3DRold
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!    PROGRAMME DE TRANSFORMEE DE FOURIER RAPIDE (1,...,N quelconque)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE BITRV(DATA,NPREV,N,NREM)                               !CRA08620
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*)                                                 !!CRA08630
!                                                                       !!CRA08640
      IP0=2                                                             !!CRA08650
      IP1=IP0*NPREV                                                     !!CRA08660
      IP4=IP1*N                                                         !!CRA08670
      IP5=IP4*NREM                                                      !!CRA08680
      I4REV=1                                                           !!CRA08690
      DO 60 I4=1,IP4,IP1                                                !CRA08700
      IF(I4-I4REV) 10,10,30                                             !CRA08710
   10 I1MAX=I4+IP1-IP0                                                  !CRA08720
      DO 20 I1=I4,I1MAX,IP0                                             !CRA08730
      DO 20 I5=I1,IP5,IP4                                               !CRA08740
      I5REV=I4REV+I5-I4                                                 !CRA08750
      TEMPR=DATA(I5)                                                    !CRA08760
      TEMPI=DATA(I5+1)                                                  !CRA08770
      DATA(I5)=DATA(I5REV)                                              !CRA08780
      DATA(I5+1)=DATA(I5REV+1)                                          !CRA08790
      DATA(I5REV)=TEMPR                                                 !CRA08800
   20 DATA(I5REV+1)=TEMPI                                               !CRA08810
   30 IP2=IP4/2                                                         !CRA08820
   40 IF(I4REV-IP2) 60,60,50                                            !CRA08830
   50 I4REV=I4REV-IP2                                                   !CRA08840
      IP2=IP2/2                                                         !CRA08850
      IF(IP2-IP1) 60,40,40                                              !CRA08860
   60 I4REV=I4REV+IP2                                                   !CRA08870
      RETURN                                                            !CRA08880
      END                                                               !CRA08890
      SUBROUTINE DFT2(DATA,NPREV,N,NREM,ISIGN)                          !CRA08900
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*)                                                 !CRA08910
      TWOPI=8.*DATAN(1.D0)*DFLOAT(ISIGN)                                !CRA08920
      IP0=2                                                             !CRA08930
      IP1=IP0*NPREV                                                     !CRA08940
      IP4=IP1*N                                                         !CRA08950
      IP5=IP4*NREM                                                      !CRA08960
      IP2=IP1                                                           !CRA08970
      NPART=N                                                           !CRA08980
   10 IF(NPART-2) 70,30,20                                              !CRA08990
   20 NPART=NPART/4                                                     !CRA09000
      GO TO 10                                                          !CRA09010
!                                                                      !CRA09020
   30 IF(IP2-IP4) 40,170,170                                            !CRA09030
   40 IP3=IP2*2                                                         !CRA09040
      DO 50 I1=1,IP1,IP0                                                !CRA09050
      DO 50 I5=I1,IP5,IP3                                               !CRA09060
      I3A=I5                                                            !CRA09070
      I3B=I3A+IP2                                                       !CRA09080
      TEMPR=DATA(I3B)                                                   !CRA09090
      TEMPI=DATA(I3B+1)                                                 !CRA09100
      DATA (I3B  )=DATA(I3A  )-TEMPR                                    !CRA09110
      DATA (I3B+1)=DATA(I3A+1)-TEMPI                                    !CRA09120
      DATA (I3A  )=DATA(I3A  )+TEMPR                                    !CRA09130
   50 DATA (I3A+1)=DATA(I3A+1)+TEMPI                                    !CRA09140
!                                                                      !CRA09150
   60 IP2=IP3                                                           !CRA09160
   70 IF(IP2-IP4) 80,170,170                                            !CRA09170
   80 IP3=IP2*4                                                         !CRA09180
      THETA=TWOPI/DFLOAT(IP3/IP1)                                       !CRA09190
      SINTH=DSIN(THETA/2.)                                              !CRA09200
      WSTPR=-2.*SINTH*SINTH                                             !CRA09210
      WSTPI=DSIN(THETA)                                                 !CRA09220
      WR=1.                                                             !CRA09230
      WI=0.                                                             !CRA09240
      DO 160 I2=1,IP2,IP1                                               !CRA09250
      IF(I2-1) 100,100,90                                               !CRA09260
   90 W2R=WR*WR-WI*WI                                                   !CRA09270
      W2I=2.*WR*WI                                                      !CRA09280
      W3R=W2R*WR-W2I*WI                                                 !CRA09290
      W3I=W2R*WI+W2I*WR                                                 !CRA09300
  100 I1MAX=I2+IP1-IP0                                                  !CRA09310
      DO 150 I1=I2,I1MAX,IP0                                            !CRA09320
      DO 150 I5=I1,IP5,IP3                                              !CRA09330
      I3A=I5                                                            !CRA09340
      I3B=I3A+IP2                                                       !CRA09350
      I3C=I3B+IP2                                                       !CRA09360
      I3D=I3C+IP2                                                       !CRA09370
      IF(I2-1) 120,120,110                                              !CRA09380
  110 TEMPR=DATA(I3B)                                                   !CRA09390
      DATA (I3B  )=W2R*TEMPR-W2I*DATA(I3B+1)                            !CRA09400
      DATA (I3B+1)=W2R*DATA(I3B+1)+W2I*TEMPR                            !CRA09410
      TEMPR=DATA(I3C)                                                   !CRA09420
      DATA(I3C)=WR*TEMPR-WI*DATA(I3C+1)                                 !CRA09430
      DATA(I3C+1)=WR*DATA(I3C+1)+WI*TEMPR                               !CRA09440
      TEMPR=DATA(I3D)                                                   !CRA09450
      DATA(I3D)=W3R*TEMPR-W3I*DATA(I3D+1)                               !CRA09460
      DATA(I3D+1)=W3R*DATA(I3D+1)+W3I*TEMPR                             !CRA09470
  120 T0R=DATA(I3A  )+DATA(I3B  )                                       !CRA09480
      T0I=DATA(I3A+1)+DATA(I3B+1)                                       !CRA09490
      T1R=DATA(I3A  )-DATA(I3B  )                                       !CRA09500
      T1I=DATA(I3A+1)-DATA(I3B+1)                                       !CRA09510
      T2R=DATA(I3C  )+DATA(I3D  )                                       !CRA09520
      T2I=DATA(I3C+1)+DATA(I3D+1)                                       !CRA09530
      T3R=DATA(I3C  )-DATA(I3D  )                                       !CRA09540
      T3I=DATA(I3C+1)-DATA(I3D+1)                                       !CRA09550
      DATA(I3A  )=T0R+T2R                                               !CRA09560
      DATA(I3A+1)=T0I+T2I                                               !CRA09570
      DATA(I3C  )=T0R-T2R                                               !CRA09580
      DATA(I3C+1)=T0I-T2I                                               !CRA09590
      IF(ISIGN) 130,130,140                                             !CRA09600
  130 T3R=-T3R                                                          !CRA09610
      T3I=-T3I                                                          !CRA09620
  140 DATA(I3B  )=T1R-T3I                                               !CRA09630
      DATA(I3B+1)=T1I+T3R                                               !CRA09640
      DATA(I3D  )=T1R+T3I                                               !CRA09650
  150 DATA(I3D+1)=T1I-T3R                                               !CRA09660
      TEMP=WR                                                           !CRA09670
      WR=WSTPR*WR-WSTPI*WI+WR                                           !CRA09680
      WI=WSTPR*WI+WSTPI*TEMP+WI                                         !CRA09690
      FACTR=1.5-0.5*(WR*WR+WI*WI)                                       !CRA09700
      WR=WR*FACTR                                                       !CRA09710
  160 WI=WI*FACTR                                                       !CRA09720
      GO TO 60                                                          !CRA09730
  170 RETURN                                                            !CRA09740
      END                                                               !CRA09750
      SUBROUTINE FOUR2(DATA,N,NDIM,ISIGN,IFORM)                         !CRA07580
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*),N(*)                                            !CRA07590
      IF(IABS(NDIM-1)-NDIM+1+IABS(IABS(ISIGN)-1)+IABS(IFORM)/2) 5,5,110 !CRA07600
    5 NTOT=1                                                            !CRA07610
      DO 10 IDIM=1,NDIM                                                 !CRA07620
      IF(N(IDIM)) 110,110,8                                             !CRA07630
!   8 IF(N(IDIM)-2**IFIX(ALOG(FLOAT(N(IDIM)))/ALOG(2.)+0.5)) 
    8  MDIM=N(IDIM)-2**IDINT(DLOG(DFLOAT(N(IDIM)))/DLOG(2.D0)+0.5D0)
      IF(MDIM) 130,10,130 
   10 NTOT=NTOT*N(IDIM)                                                 !CRA07650
      IF(IFORM) 70,20,20                                                !CRA07660
   20 NREM=NTOT                                                         !CRA07670
!                                                                      !CRA07680
      DO 60 IDIM=1,NDIM                                                 !CRA07690
      NREM=NREM/N(IDIM)                                                 !CRA07700
      NPREV=NTOT/(N(IDIM)*NREM)                                         !CRA07710
      NCURR=N(IDIM)                                                     !CRA07720
      IF(IDIM-1+IFORM) 30,30,40                                         !CRA07730
   30 NCURR=NCURR/2                                                     !CRA07740
   40 CALL BITRV(DATA,NPREV,NCURR,NREM)                                 !CRA07750
      CALL DFT2 (DATA,NPREV,NCURR,NREM,ISIGN)                           !CRA07760
      IF(IDIM-1+IFORM) 50,50,60                                         !CRA07770
!  50 CALL FIXRL(DATA,N(1),NREM,ISIGN,IFORM)                           !CRA07780
   50 CALL FIXRL(DATA,N,NREM,ISIGN,IFORM)
      NTOT=(NTOT/N(1))*(N(1)/2+1)                                       !CRA07790
   60 CONTINUE                                                          !CRA07800
      RETURN                                                            !CRA07810
   70 NTOT=(NTOT/N(1))*(N(1)/2+1)                                       !CRA07820
      NREM=1                                                            !CRA07830
      DO 100 JDIM=1,NDIM                                                !CRA07840
      IDIM=NDIM+1-JDIM                                                  !CRA07850
      NCURR=N(IDIM)                                                     !CRA07860
      IF(IDIM-1) 80,80,90                                               !CRA07870
   80 NCURR=NCURR/2                                                     !CRA07880
!     CALL FIXRL(DATA,N(1),NREM,ISIGN,IFORM)                           !CRA07890
       CALL FIXRL(DATA,N,NREM,ISIGN,IFORM)
      NTOT=NTOT/(N(1)/2+1)*N(1)                                         !CRA07900
   90 NPREV=NTOT/(N(IDIM)*NREM)                                         !CRA07910
      CALL BITRV(DATA,NPREV,NCURR,NREM)                                 !CRA07920
      CALL DFT2(DATA,NPREV,NCURR,NREM,ISIGN)                            !CRA07930
      NREM=NREM*N(IDIM)                                                 !CRA07940
  100 CONTINUE                                                          !CRA07950
      RETURN                                                            !CRA07960
  110 PRINT 120, NDIM,ISIGN,IFORM 
      RETURN                                                            !CRA07980
!    POUR LE CAS DE FREQUENCE DE COUPURE QUI N'EST PAS PUISSANCE DE 2  !CRA07990
  130 CALL FOURT(DATA,N,NDIM,ISIGN,IFORM)                               !CRA08000
      RETURN                                                            !CRA08010
  120 FORMAT(1X,'ERROR IN FOUR2. EITHER NDIM= ',I10,' OR ISIGN = ',I10,& !CRA08020
      ' OR IFORM = ',I10,' HAS AN ILLEGAL VALUE')                       !CRA08030
      END                                                               !CRA08040
      SUBROUTINE GOERT(DATA,NPREV,IPROD,IFACT,IREM,WORK,WMINR,WMINI,&   !CRA13410
        ROOTR,ROOTI,ROTHR,ROTHI)                                       !CRA13420
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*),WORK(*)                                         !CRA13430
      IP0=2                                                             !CRA13440
      IP1=IP0*NPREV                                                     !CRA13450
      IP2=IP1*IPROD                                                     !CRA13460
      IP3=IP2*IFACT                                                     !CRA13470
      IP5=IP3*IREM                                                      !CRA13480
      IF(WMINI)10,40,10                                                 !CRA13490
   10 WR=WMINR                                                          !CRA13500
      WI=WMINI                                                          !CRA13510
      I3MIN=1+IP2                                                       !CRA13520
      DO 30 I3=I3MIN,IP3,IP2                                            !CRA13530
      I1MAX=I3+IP1-IP0                                                  !CRA13540
      DO 20 I1=I3,I1MAX,IP0                                             !CRA13550
      DO 20 I5=I1,IP5,IP3                                               !CRA13560
      TEMPR=DATA(I5)                                                    !CRA13570
      DATA(I5)=WR*TEMPR-WI*DATA(I5+1)                                   !CRA13580
   20 DATA(I5+1)=WR*DATA(I5+1)+WI*TEMPR                                 !CRA13590
      TEMPR=WR                                                          !CRA13600
      WR=WMINR*TEMPR-WMINI*WI                                           !CRA13610
   30 WI=WMINR*WI+WMINI*TEMPR                                           !CRA13620
   40 DO 210 I1=1,IP1,IP0                                               !CRA13630
      DO 210 I5=I1,IP5,IP3                                              !CRA13640
      SUM1R=0.                                                          !CRA13650
      SUM1I=0.                                                          !CRA13660
      I3MAX=I5+IP3-IP2                                                  !CRA13670
      DO 50 I3=I5,I3MAX,IP2                                             !CRA13680
      SUM1R=SUM1R+DATA(I3  )                                            !CRA13690
   50 SUM1I=SUM1I+DATA(I3+1)                                            !CRA13700
      WORK(1)=SUM1R                                                     !CRA13710
      WORK(2)=SUM1I                                                     !CRA13720
      DO 200 IREGN=1,3                                                  !CRA13730
      IF(IREGN-2)60,70,80                                               !CRA13740
   60 IWMIN=IP0+1                                                       !CRA13750
      IWMAX=IP0*(IFACT/8)+1                                             !CRA13760
      WR=0.                                                             !CRA13770
      WI=0.                                                             !CRA13780
      GO TO 90                                                          !CRA13790
   70 IWMIN=IWMAX+IP0                                                   !CRA13800
      IWMAX=IP0*(3*(IFACT+1)/8)+1                                       !CRA13810
      WR=WR+1.                                                          !CRA13820
      GO TO 100                                                         !CRA13830
   80 IWMIN=IP0*(IFACT/2+1)+1                                           !CRA13840
      IWMAX=IP0*IFACT+2-IWMAX                                           !CRA13850
      WR=ROOTR-ROTHR                                                    !CRA13860
      WI=ROOTI-ROTHI                                                    !CRA13870
   90 IF(IWMIN-IWMAX)100,100,190                                        !CRA13880
  100 DO 180 IWORK=IWMIN,IWMAX,IP0                                      !CRA13890
      I3=I3MAX                                                          !CRA13900
      SUM1R=DATA(I3  )                                                  !CRA13910
      SUM1I=DATA(I3+1)                                                  !CRA13920
      I3=I3-IP2                                                         !CRA13930
      IF(IREGN-2)110,130,150                                            !CRA13940
  110 WR=WR+ROOTR                                                       !CRA13950
      WI=WI+ROOTI                                                       !CRA13960
      TWOWR=WR+WR                                                       !CRA13970
      SUM2R=-SUM1R                                                      !CRA13980
      SUM2I=-SUM1I                                                      !CRA13990
  120 SUM2R=SUM2R-TWOWR*SUM1R-DATA(I3  )                                !CRA14000
      SUM2I=SUM2I-TWOWR*SUM1I-DATA(I3+1)                                !CRA14010
      SUM1R=SUM1R-SUM2R                                                 !CRA14020
      SUM1I=SUM1I-SUM2I                                                 !CRA14030
      I3=I3-IP2                                                         !CRA14040
      IF(I3-I5)170,170,120                                              !CRA14050
  130 TWOWR=WR+WR                                                       !CRA14060
      SUM2R=0.                                                          !CRA14070
      SUM2I=0.                                                          !CRA14080
  140 TEMPR=SUM1R                                                       !CRA14090
      TEMPI=SUM1I                                                       !CRA14100
      SUM1R=TWOWR*SUM1R-SUM2R+DATA(I3)                                  !CRA14110
      SUM1I=TWOWR*SUM1I-SUM2I+DATA(I3+1)                                !CRA14120
      SUM2R=TEMPR                                                       !CRA14130
      SUM2I=TEMPI                                                       !CRA14140
      I3=I3-IP2                                                         !CRA14150
      IF(I3-I5)170,170,140                                              !CRA14160
  150 WR=WR-ROOTR                                                       !CRA14170
      WI=WI-ROOTI                                                       !CRA14180
      TWOWR=WR+WR                                                       !CRA14190
      SUM2R=SUM1R                                                       !CRA14200
      SUM2I=SUM1I                                                       !CRA14210
  160 SUM2R=TWOWR*SUM1R-SUM2R+DATA(I3  )                                !CRA14220
      SUM2I=TWOWR*SUM1I-SUM2I+DATA(I3+1)                                !CRA14230
      SUM1R=SUM2R-SUM1R                                                 !CRA14240
      SUM1I=SUM2I-SUM1I                                                 !CRA14250
      I3=I3-IP2                                                         !CRA14260
      IF(I3-I5) 170,170,160                                             !CRA14270
  170 TEMPR=-WI*SUM1I                                                   !CRA14280
      TEMPI=WI*SUM1R                                                    !CRA14290
      SUM1R=WR*SUM1R-SUM2R+DATA(I3  )                                   !CRA14300
      SUM1I=WR*SUM1I-SUM2I+DATA(I3+1)                                   !CRA14310
      WORK(IWORK  )=SUM1R+TEMPR                                         !CRA14320
      WORK(IWORK+1)=SUM1I+TEMPI                                         !CRA14330
      IWCNJ=IP0*(IFACT+1)-IWORK                                         !CRA14340
      WORK(IWCNJ  )=SUM1R-TEMPR                                         !CRA14350
      WORK(IWCNJ+1)=SUM1I-TEMPI                                         !CRA14360
      TEMP=WR                                                           !CRA14370
      WR=WR*ROOTR-WI*ROOTI+WR                                           !CRA14380
  180 WI=WI*ROOTR+TEMP*ROOTI+WI                                         !CRA14390
  190 WR=WR+ROOTR                                                       !CRA14400
  200 WI=WI+ROOTI                                                       !CRA14410
      IWORK=1                                                           !CRA14420
      DO 210 I3=I5,I3MAX,IP2                                            !CRA14430
      DATA(I3  )=WORK(IWORK  )                                          !CRA14440
      DATA(I3+1)=WORK(IWORK+1)                                          !CRA14450
  210 IWORK=IWORK+IP0                                                   !CRA14460
      RETURN                                                            !CRA14470
      END                                                               !CRA14480
      SUBROUTINE ASMRV(DATA,NPREV,N,NREM,IFACT,NFACT,IWORK,NWORK)       !CRA08050
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*),IFACT(*),IWORK(*),MOD(10)                       !CRA08060
      IF(NFACT-1)190,190,10                                             !CRA08070
   10 IP0=2                                                             !CRA08080
      IP1=IP0*NPREV                                                     !CRA08090
      IP2=IP1*N                                                         !CRA08100
      IP3=IP2*NREM                                                      !CRA08110
      MULT=N/IFACT(1)                                                   !CRA08120
      IPROD=MULT*IP1                                                    !CRA08130
      INVPR=IPROD                                                       !CRA08140
      DO 20 IF=2,NFACT                                                  !CRA08150
      IPROD=IPROD*IFACT(IF-1)                                           !CRA08160
      INVPR=INVPR/IFACT(IF)                                             !CRA08170
   20 MOD(IF)=IPROD-INVPR                                               !CRA08180
      I2NXT=1                                                           !CRA08190
      DO 30 IW=1,NWORK                                                  !CRA08200
      IWORK(IW)=I2NXT                                                   !CRA08210
   30 I2NXT=I2NXT+IP1                                                   !CRA08220
      I2WMX=IWORK(NWORK)                                                !CRA08230
      DO 180 I2=1,IP2,IP1                                               !CRA08240
      I2NXT=I2                                                          !CRA08250
   40 I2NXT=(I2NXT-1)*MULT+1                                            !CRA08260
      IF=NFACT                                                          !CRA08270
      GO TO 60                                                          !CRA08280
   50 I2NXT=I2NXT-MOD(IF)*((I2NXT-1)/MOD(IF))                           !CRA08290
      IF=IF-1                                                           !CRA08300
   60 IF(IF-2)70,70,50                                                  !CRA08310
   70 IQUOT=(I2NXT-1)/MOD(2)                                            !CRA08320
      IF(IQUOT-IFACT(2))90,80,80                                        !CRA08330
   80 IQUOT=IQUOT-1                                                     !CRA08340
   90 I2NXT=I2NXT-MOD(2)*IQUOT                                          !CRA08350
  100 IF(I2-I2NXT)130,180,110                                           !CRA08360
  110 IF(I2NXT-I2WMX)120,120,40                                         !CRA08370
  120 IWN=1+(I2NXT-1)/IP1                                               !CRA08380
      I2NXT=IWORK(IWN)                                                  !CRA08390
      GO TO 100                                                         !CRA08400
  130 IF(I2-I2WMX)140,140,160                                           !CRA08410
  140 IW2=1+(I2-1)/IP1                                                  !CRA08420
      ITEMP=IWORK(IW2)                                                  !CRA08430
      IWT=1+(ITEMP-1)/IP1                                               !CRA08440
      IWORK(IWT)=I2NXT                                                  !CRA08450
      IF(I2NXT-I2WMX)150,150,160                                        !CRA08460
  150 IWN=1+(I2NXT-1)/IP1                                               !CRA08470
      IWORK(IWN)=ITEMP                                                  !CRA08480
  160 I1MAX=I2+IP1-IP0                                                  !CRA08490
      DO 170 I1=I2,I1MAX,IP0                                            !CRA08500
      DO 170 I3=I1,IP3,IP2                                              !CRA08510
      I3NXT=I3+I2NXT-I2                                                 !CRA08520
      TEMPR=DATA(I3)                                                    !CRA08530
      TEMPI=DATA(I3+1)                                                  !CRA08540
      DATA(I3  )=DATA(I3NXT  )                                          !CRA08550
      DATA(I3+1)=DATA(I3NXT+1)                                          !CRA08560
      DATA(I3NXT  )=TEMPR                                               !CRA08570
  170 DATA(I3NXT+1)=TEMPI                                               !CRA08580
  180 CONTINUE                                                          !CRA08590
  190 RETURN                                                            !CRA08600
      END                                                               !CRA08610
      SUBROUTINE DFT(DATA,NPREV,N,NREM,ISIGN,IFACT,WORK)                !CRA09760
!
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*),IFACT(*),WORK(*),W(8)                           !CRA09770
      PI=DACOS(-1.D0)
!     TWOPI=6.2831853071795865*DFLOAT(ISIGN)                           !CRA09780
      TWOPI=2.D0*PI*DFLOAT(ISIGN) 
      IP0=2                                                             !CRA09790
      IP1=IP0*NPREV                                                     !CRA09800
      IP4=IP1*N                                                         !CRA09810
      IP5=IP4*NREM                                                      !CRA09820
      IF=0                                                              !CRA09830
      WSTPR=0.                                                          !CRA09840
      WSTPI=0.                                                          !CRA09850
      IP2=IP1                                                           !CRA09860
   10 IF(IP2-IP4)20,330,330                                             !CRA09870
   20 IF=IF+1                                                           !CRA09880
      IFCUR=IFACT(IF)                                                   !CRA09890
      IF(IFCUR-2)330,30,60                                              !CRA09900
   30 IF(4*IP2-IP4)40,40,70                                             !CRA09910
   40 IF(IFACT(IF+1)-2)70,50,70                                         !CRA09920
   50 IF=IF+1                                                           !CRA09930
      IFCUR=4                                                           !CRA09940
      GO TO 70                                                          !CRA09950
   60 THETA=TWOPI/DFLOAT(IFCUR)                                         !CRA09960
      ROTHR=-2.*DSIN(THETA/4.)**2                                       !CRA09970
      ROTHI=DSIN(THETA/2.)                                              !CRA09980
      ROOTR=-2.*ROTHI*ROTHI                                             !CRA09990
      ROOTI=DSIN(THETA)                                                 !CRA10000
      ROT2I=2.*ROOTI*(1.+ROOTR)                                         !CRA10010
   70 IWMAX=IP0*(IFCUR-1)                                               !CRA10020
      IP3=IP2*IFCUR                                                     !CRA10030
      W(1)=1.                                                           !CRA10040
      W(2)=0.                                                           !CRA10050
      IF(IP1-IP2)80,90,90                                               !CRA10060
   80 THETA=TWOPI/DFLOAT(IP3/IP1)                                       !CRA10070
      SINTH=DSIN(THETA/2.)                                              !CRA10080
      WSTPR=-2.*SINTH*SINTH                                             !CRA10090
      WSTPI=DSIN(THETA)                                                 !CRA10100
   90 DO 320 I2=1,IP2,IP1                                               !CRA10110
      IF(IFCUR-5)100,100,300                                            !CRA10120
  100 IF(I2-1)140,140,110                                               !CRA10130
  110 IF(IFCUR-2)330,140,120                                            !CRA10140
  120 IWMIN=IP0+1                                                       !CRA10150
      DO 130 IW=IWMIN,IWMAX,IP0                                         !CRA10160
      W(IW  )=W(IW-2)*W(1)-W(IW-1)*W(2)                                 !CRA10170
  130 W(IW+1)=W(IW-2)*W(2)+W(IW-1)*W(1)                                 !CRA10180
  140 I1MAX=I2+IP1-IP0                                                  !CRA10190
      DO 290 I1=I2,I1MAX,IP0                                            !CRA10200
      DO 290 J0=I1,IP5,IP3                                              !CRA10210
      IF(IFCUR-4)150,230,150                                            !CRA10220
  150 IF(I2-1)180,180,160                                               !CRA10230
  160 J1=J0                                                             !CRA10240
      DO 170 IW=1,IWMAX,IP0                                             !CRA10250
      J1=J1+IP2                                                         !CRA10260
      TEMPR=DATA(J1)                                                    !CRA10270
      DATA(J1)=TEMPR*W(IW)-DATA(J1+1)*W(IW+1)                           !CRA10280
  170 DATA(J1+1)=TEMPR*W(IW+1)+DATA(J1+1)*W(IW)                         !CRA10290
  180 J1=J0+IP2                                                         !CRA10300
      J2=J0+IP3-IP2                                                     !CRA10310
      IF(IFCUR-2)330,190,200                                            !CRA10320
  190 J1=J0                                                             !CRA10330
  200 TEMPR=DATA(J2  )                                                  !CRA10340
      TEMPI=DATA(J2+1)                                                  !CRA10350
      DATA(J2  )=DATA(J1  )-TEMPR                                       !CRA10360
      DATA(J2+1)=DATA(J1+1)-TEMPI                                       !CRA10370
      DATA(J1  )=DATA(J1  )+TEMPR                                       !CRA10380
      DATA(J1+1)=DATA(J1+1)+TEMPI                                       !CRA10390
      J1=J1+IP2                                                         !CRA10400
      J2=J2-IP2                                                         !CRA10410
      IF(J1-J2)200,210,210                                              !CRA10420
  210 IF(IFCUR-3)290,220,280                                            !CRA10430
  220 J1=J0+IP2                                                         !CRA10440
      J2=J1+IP2                                                         !CRA10450
      TEMPR=DATA(J0  )-0.5*DATA(J1  )                                   !CRA10460
      TEMPI=DATA(J0+1)-0.5*DATA(J1+1)                                   !CRA10470
      DATA(J0  )=DATA(J0  )+DATA(J1  )                                  !CRA10480
      DATA(J0+1)=DATA(J0+1)+DATA(J1+1)                                  !CRA10490
      DIFR=-ROOTI*DATA(J2+1)                                            !CRA10500
      DIFI=ROOTI*DATA(J2)                                               !CRA10510
      DATA(J1  )=TEMPR+DIFR                                             !CRA10520
      DATA(J1+1)=TEMPI+DIFI                                             !CRA10530
      DATA(J2  )=TEMPR-DIFR                                             !CRA10540
      DATA(J2+1)=TEMPI-DIFI                                             !CRA10550
      GO TO 290                                                         !CRA10560
  230 J1=J0+IP2                                                         !CRA10570
      J2=J1+IP2                                                         !CRA10580
      J3=J2+IP2                                                         !CRA10590
      IF(I2-1)250,250,240                                               !CRA10600
  240 TEMPR=DATA(J1)                                                    !CRA10610
      DATA(J1)=W(3)*TEMPR-W(4)*DATA(J1+1)                               !CRA10620
      DATA(J1+1)=W(3)*DATA(J1+1)+W(4)*TEMPR                             !CRA10630
      TEMPR=DATA(J2)                                                    !CRA10640
      DATA(J2)=W(1)*TEMPR-W(2)*DATA(J2+1)                               !CRA10650
      DATA(J2+1)=W(1)*DATA(J2+1)+W(2)*TEMPR                             !CRA10660
      TEMPR=DATA(J3)                                                    !CRA10670
      DATA(J3)=W(5)*TEMPR-W(6)*DATA(J3+1)                               !CRA10680
      DATA(J3+1)=W(5)*DATA(J3+1)+W(6)*TEMPR                             !CRA10690
  250 T0R=DATA(J0)+DATA(J1)                                             !CRA10700
      T0I=DATA(J0+1)+DATA(J1+1)                                         !CRA10710
      T1R=DATA(J0  )-DATA(J1  )                                         !CRA10720
      T1I=DATA(J0+1)-DATA(J1+1)                                         !CRA10730
      T2R=DATA(J2  )+DATA(J3  )                                         !CRA10740
      T2I=DATA(J2+1)+DATA(J3+1)                                         !CRA10750
      T3R=DATA(J2  )-DATA(J3  )                                         !CRA10760
      T3I=DATA(J2+1)-DATA(J3+1)                                         !CRA10770
      DATA(J0  )=T0R+T2R                                                !CRA10780
      DATA(J0+1)=T0I+T2I                                                !CRA10790
      DATA(J2  )=T0R-T2R                                                !CRA10800
      DATA(J2+1)=T0I-T2I                                                !CRA10810
      IF(ISIGN)260,260,270                                              !CRA10820
  260 T3R=-T3R                                                          !CRA10830
      T3I=-T3I                                                          !CRA10840
  270 DATA(J1  )=T1R-T3I                                                !CRA10850
      DATA(J1+1)=T1I+T3R                                                !CRA10860
      DATA(J3  )=T1R+T3I                                                !CRA10870
      DATA(J3+1)=T1I-T3R                                                !CRA10880
      GO TO 290                                                         !CRA10890
  280 J1=J0+IP2                                                         !CRA10900
      J2=J1+IP2                                                         !CRA10910
      J3=J2+IP2                                                         !CRA10920
      J4=J3+IP2                                                         !CRA10930
      SUMR=DATA(J1  )+DATA(J2  )                                        !CRA10940
      SUMI=DATA(J1+1)+DATA(J2+1)                                        !CRA10950
      DATA(J0  )=DATA(J0  )+SUMR                                        !CRA10960
      DATA(J0+1)=DATA(J0+1)+SUMI                                        !CRA10970
      SUMR=DATA(J0  )-1.25*SUMR                                         !CRA10980
      SUMI=DATA(J0+1)-1.25*SUMI                                         !CRA10990
      DIFR=(ROOTR+1.25)*(DATA(J1  )-DATA(J2  ))                         !CRA11000
      DIFI=(ROOTR+1.25)*(DATA(J1+1)-DATA(J2+1))                         !CRA11010
      T1R=SUMR+DIFR                                                     !CRA11020
      T1I=SUMI+DIFI                                                     !CRA11030
      T2R=SUMR-DIFR                                                     !CRA11040
      T2I=SUMI-DIFI                                                     !CRA11050
      SUMR=-ROT2I*DATA(J3+1)-ROOTI*DATA(J4+1)                           !CRA11060
      SUMI= ROT2I*DATA(J3  )+ROOTI*DATA(J4  )                           !CRA11070
      DIFR= ROT2I*DATA(J4+1)-ROOTI*DATA(J3+1)                           !CRA11080
      DIFI=-ROT2I*DATA(J4  )+ROOTI*DATA(J3  )                           !CRA11090
      DATA(J1  )=T1R+SUMR                                               !CRA11100
      DATA(J1+1)=T1I+SUMI                                               !CRA11110
      DATA(J2  )=T2R-DIFR                                               !CRA11120
      DATA(J2+1)=T2I-DIFI                                               !CRA11130
      DATA(J3  )=T2R+DIFR                                               !CRA11140
      DATA(J3+1)=T2I+DIFI                                               !CRA11150
      DATA(J4  )=T1R-SUMR                                               !CRA11160
      DATA(J4+1)=T1I-SUMI                                               !CRA11170
  290 CONTINUE                                                          !CRA11180
      GO TO 310                                                         !CRA11190
  300 CALL GOERT(DATA(I2),NPREV,IP2/IP1,IFCUR,IP5/IP3,WORK,W(1),W(2),&  !CRA11200
        ROOTR,ROOTI,ROTHR,ROTHI)                                       !CRA11210
  310 TEMP=W(1)                                                         !CRA11220
      W(1)=TEMP*WSTPR-W(2)*WSTPI+TEMP                                   !CRA11230
  320 W(2)=W(2)*WSTPR+TEMP*WSTPI+W(2)                                   !CRA11240
      IP2=IP3                                                           !CRA11250
      GO TO 10                                                          !CRA11260
  330 RETURN                                                            !CRA11270
      END                                                               !CRA11280
      SUBROUTINE FACTR(N,IFACT,NFACT,ISYM,IFSYM,NFSYM,ICENT,IFCNT,NFCNT,&!CRA11290
      IFMAX)                                                           !CRA11300
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION IFACT(*),IFSYM(*),IFCNT(*)                              !CRA11310
      IF=0                                                              !CRA11320
      IFS=0                                                             !CRA11330
      IFC=0                                                             !CRA11340
      ISYM=1                                                            !CRA11350
      ICENT=1                                                           !CRA11360
      NREM=N                                                            !CRA11370
      IDIV=2                                                            !CRA11380
      GO TO 30                                                          !CRA11390
   10 IDIV=1                                                            !CRA11400
   20 IDIV=IDIV+2                                                       !CRA11410
   30 IQUOT=NREM/IDIV                                                   !CRA11420
      IF(NREM-IQUOT*IDIV)70,40,70                                       !CRA11430
   40 NREM=IQUOT                                                        !CRA11440
      IQUOT=NREM/IDIV                                                   !CRA11450
      IF(NREM-IQUOT*IDIV)60,50,60                                       !CRA11460
   50 NREM=IQUOT                                                        !CRA11470
      IF=IF+1                                                           !CRA11480
      IFACT(IF)=IDIV                                                    !CRA11490
      GO TO 30                                                          !CRA11500
   60 IFC=IFC+1                                                         !CRA11510
      IFCNT(IFC)=IDIV                                                   !CRA11520
      ICENT=ICENT*IDIV                                                  !CRA11530
   70 IF(IQUOT-IDIV)90,80,80                                            !CRA11540
   80 IF(IDIV-2)10,10,20                                                !CRA11550
   90 IFMAX=IDIV                                                        !CRA11560
      IF(NREM-1)110,110,100                                             !CRA11570
  100 IFC=IFC+1                                                         !CRA11580
      IFCNT(IFC)=NREM                                                   !CRA11590
      ICENT=ICENT*NREM                                                  !CRA11600
      IFMAX=NREM                                                        !CRA11610
  110 NFSYM=IF                                                          !CRA11620
      KFS=NFSYM                                                         !CRA11630
      NFCNT=IFC                                                         !CRA11640
      IF(NFCNT)140,140,120                                              !CRA11650
  120 DO 130 IFC=1,NFCNT                                                !CRA11660
      IF=IF+1                                                           !CRA11670
  130 IFACT(IF)=IFCNT(IFC)                                              !CRA11680
      KFS=KFS+1                                                         !CRA11690
      IFSYM(NFSYM+1)=ICENT                                              !CRA11700
  140 IF(NFSYM)170,170,150                                              !CRA11710
  150 DO 160 IFS=1,NFSYM                                                !CRA11720
      IFSYM(IFS)=IFACT(IFS)                                             !CRA11730
      ISYM=ISYM*IFACT(IFS)                                              !CRA11740
      IF=IF+1                                                           !CRA11750
      JFS=NFSYM+1-IFS                                                   !CRA11760
      IFACT(IF)=IFACT(JFS)                                              !CRA11770
      KFS=KFS+1                                                         !CRA11780
  160 IFSYM(KFS)=IFACT(JFS)                                             !CRA11790
  170 NFACT=IF                                                          !CRA11800
      NFSYM=KFS                                                         !CRA11810
      RETURN                                                            !CRA11820
      END                                                               !CRA11830
      SUBROUTINE FIXRL(DATA,N,NREM,ISIGN,IFORM)                         !CRA11840
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*)                                                 !CRA11850
      PI=DACOS(-1.D0)
      TWOPI=2.D0*PI*DFLOAT(ISIGN)  
!     TWOPI=6.28318 53071 79586 5 *FLOAT(ISIGN)                        !CRA11860
      IP0=2                                                             !CRA11870
      IP1=IP0*(N/2)                                                     !CRA11880
      IP2=IP1*NREM                                                      !CRA11890
      IF(IFORM) 10,70,70                                                !CRA11900
   10 J1=IP1+1                                                          !CRA11910
      DATA(2)=DATA(J1)                                                  !CRA11920
      IF(NREM-1) 70,70,20                                               !CRA11930
   20 J1=J1+IP0                                                         !CRA11940
      I2MIN=IP1+1                                                       !CRA11950
      DO 60 I2=I2MIN,IP2,IP1                                            !CRA11960
      DATA(I2)=DATA(J1)                                                 !CRA11970
      J1=J1+IP0                                                         !CRA11980
      IF(N-2) 50,50,30                                                  !CRA11990
   30 I1MIN=I2+IP0                                                      !CRA12000
      I1MAX=I2+IP1-IP0                                                  !CRA12010
      DO 40 I1=I1MIN,I1MAX,IP0                                          !CRA12020
      DATA(I1  )=DATA(J1  )                                             !CRA12030
      DATA(I1+1)=DATA(J1+1)                                             !CRA12040
   40 J1=J1+IP0                                                         !CRA12050
   50 DATA(I2+1)=DATA(J1)                                               !CRA12060
   60 J1=J1+IP0                                                         !CRA12070
   70 DO 80 I2=1,IP2,IP1                                                !CRA12080
      TEMPR=DATA(I2)                                                    !CRA12090
      DATA(I2)=DATA(I2)+DATA(I2+1)                                      !CRA12100
   80 DATA(I2+1)=TEMPR-DATA(I2+1)                                       !CRA12110
      IF(N-2) 200,200,90                                                !CRA12120
   90 THETA=TWOPI/FLOAT(N)                                              !CRA12130
      SINTH=DSIN(THETA/2.)                                              !CRA12140
      ZSTPR=-2.*SINTH*SINTH                                             !CRA12150
      ZSTPI=DSIN(THETA)                                                 !CRA12160
      ZR=(1.-ZSTPI)/2.                                                  !CRA12170
      ZI=(1.+ZSTPR)/2.                                                  !CRA12180
      IF(IFORM) 100,110,110                                             !CRA12190
  100 ZR=1.-ZR                                                          !CRA12200
      ZI=-ZI                                                            !CRA12210
  110 I1MIN=IP0+1                                                       !CRA12220
      I1MAX=IP0*(N/4)+1                                                 !CRA12230
      DO 190 I1=I1MIN,I1MAX,IP0                                         !CRA12240
      DO 180 I2=I1,IP2,IP1                                              !CRA12250
      I2CNJ=I2+IP0*(N/2)-2*(I1-1)                                       !CRA12260
      IF(I2-I2CNJ) 150,120,120                                          !CRA12270
  120 IF(ISIGN*(2*IFORM+1)) 130,140,140                                !CRA12280
  130 DATA(I2+1)=-DATA(I2+1)                                            !CRA12290
  140 IF(IFORM) 170,180,180                                             !CRA12300
  150 DIFR=DATA(I2)-DATA(I2CNJ)                                         !CRA12310
      DIFI=DATA(I2+1)+DATA(I2CNJ+1)                                     !CRA12320
      TEMPR=DIFR*ZR-DIFI*ZI                                             !CRA12330
      TEMPI=DIFR*ZI+DIFI*ZR                                             !CRA12340
      DATA(I2  )=DATA(I2  )-TEMPR                                       !CRA12350
      DATA(I2+1)=DATA(I2+1)-TEMPI                                       !CRA12360
      DATA(I2CNJ  )=DATA(I2CNJ  )+TEMPR                                 !CRA12370
      DATA(I2CNJ+1)=DATA(I2CNJ+1)-TEMPI                                 !CRA12380
      IF(IFORM) 160,180,180                                             !CRA12390
  160 DATA(I2CNJ  )=DATA(I2CNJ  )+DATA(I2CNJ  )                         !CRA12400
      DATA(I2CNJ+1)=DATA(I2CNJ+1)+DATA(I2CNJ+1)                         !CRA12410
  170 DATA(I2  )=DATA(I2  )+DATA(I2  )                                  !CRA12420
      DATA(I2+1)=DATA(I2+1)+DATA(I2+1)                                  !CRA12430
  180 CONTINUE                                                          !CRA12440
      TEMP=ZR-0.5                                                       !CRA12450
      ZR=ZSTPR*TEMP-ZSTPI*ZI+ZR                                         !CRA12460
      ZI=ZSTPR*ZI+ZSTPI*TEMP+ZI                                         !CRA12470
      FACTR=1.5-2.*((ZR-0.5)*(ZR-0.5)+ZI*ZI)                            !CRA12480
      ZR=0.5+(ZR-0.5)*FACTR                                             !CRA12490
  190 ZI=ZI*FACTR                                                       !CRA12500
  200 IF(IFORM) 270,210,210                                             !CRA12510
  210 I2=IP2+1                                                          !CRA12520
      I1=I2                                                             !CRA12530
      J1=IP0*(N/2+1)*NREM+1                                             !CRA12540
      GO TO 250                                                         !CRA12550
  220 DATA(J1)=DATA(I1)                                                 !CRA12560
      DATA(J1+1)=DATA(I1+1)                                             !CRA12570
      I1=I1-IP0                                                         !CRA12580
      J1=J1-IP0                                                         !CRA12590
  230 IF(I2-I1) 220,240,240                                             !CRA12600
  240 DATA(J1)=DATA(I1)                                                 !CRA12610
      DATA(J1+1)=0.                                                     !CRA12620
  250 I2=I2-IP1                                                         !CRA12630
      J1=J1-IP0                                                         !CRA12640
      DATA(J1)=DATA(I2+1)                                               !CRA12650
      DATA(J1+1)=0.                                                     !CRA12660
      I1=I1-IP0                                                         !CRA12670
      J1=J1-IP0                                                         !CRA12680
      IF(I2-1) 260,260,230                                              !CRA12690
  260 DATA(2)=0.                                                        !CRA12700
  270 RETURN                                                            !CRA12710
      END                                                               !CRA12720
      SUBROUTINE FOURT(DATA,N,NDIM,ISIGN,IFORM)                         !CRA12730
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*),N(*)                                            !CRA12740
      DIMENSION IFACT(32),IFSYM(32),IFCNT(10),WORK(2,50000),IWORK(1)    !CRA12750
      EQUIVALENCE (WORK(1,1),IWORK(1))                                  !CRA12760
      NWORK=50000                                                       !CRA12770
      IF(IFORM)10,10,40                                                 !CRA12780
   10 IF(N(1)-2*(N(1)/2))20,40,20                                       !CRA12790
   20 PRINT 30,IFORM,N(1)                                               !CRA12800
   30 FORMAT('ya un pb: IFORM =' ,I2,' N(1) = ',I5,' NO TRANSFORM DONE')       !CRA12820
      RETURN                                                            !CRA12830
   40 NTOT=1                                                            !CRA12840
      DO 70 IDIM=1,NDIM                                                 !CRA12850
      IF(N(IDIM))50,50,70                                               !CRA12860
   50 PRINT 60,IDIM,N(IDIM)                                             !CRA12870
   60 FORMAT(20H0ERROR IN FOURT.  N(,I2,4H) = ,I5,41H IS ZERO OR NEGATIVE.  NO TRANSFORM DONE.)                                           !CRA12890
      RETURN                                                            !CRA12900
   70 NTOT=NTOT*N(IDIM)                                                 !CRA12910
!...                                                                    !CRA12920
!    SQRTN=SQRT(FLOAT(NTOT))                                           !CRA12930
!...                                                                    !CRA12940
      NREM=NTOT                                                         !CRA12950
      IMAX=(IFORM+1)*NTOT                                               !CRA12960
      IF(IFORM)80,90,90                                                 !CRA12970
   80 NREM=1                                                            !CRA12980
      NTOT=(NTOT/N(1))*(N(1)/2+1)                                       !CRA12990
      IMAX=2*NTOT                                                       !CRA13000
!...                                                                    !CRA13010
! 90 DO 100 I=1,IMAX                                                   !CRA13020
!100 DATA(I)=DATA(I)/SQRTN                                             !CRA13030
   90 CONTINUE                                                          !CRA13040
!...                                                                    !CRA13050
      DO 260 JDIM=1,NDIM                                                !CRA13060
      IF(IFORM)110,120,120                                              !CRA13070
  110 IDIM=NDIM+1-JDIM                                                  !CRA13080
      GO TO 130                                                         !CRA13090
  120 IDIM=JDIM                                                         !CRA13100
      NREM=NREM/N(IDIM)                                                 !CRA13110
  130 NCURR=N(IDIM)                                                     !CRA13120
      IF(IDIM-1)140,140,170                                             !CRA13130
  140 IF(IFORM)150,160,170                                              !CRA13140
  150 CALL FIXRL(DATA,N(1),NREM,ISIGN,IFORM)                            !CRA13150
      NTOT=(NTOT/(N(1)/2+1))*N(1)                                       !CRA13160
  160 NCURR=NCURR/2                                                     !CRA13170
  170 IF(NCURR-1)220,220,180                                            !CRA13180
  180 CALL FACTR(NCURR,IFACT,NFACT,ISYM,IFSYM,NFSYM,ICENT,IFCNT,NFCNT,IFMAX)                                                            !CRA13200
      IF(IFMAX-NWORK)210,210,190                                        !CRA13210
  190 PRINT 200,NWORK,IFMAX,IDIM,N(IDIM)                                !CRA13220
  200 FORMAT(I3,I5,I1,I5)                                       !CRA13250
      RETURN                                                            !CRA13260
  210 NPREV=NTOT/(N(IDIM)*NREM)                                         !CRA13270
      CALL SYMRV(DATA,NPREV,NCURR,NREM,IFSYM,NFSYM)                     !CRA13280
      CALL ASMRV(DATA,NPREV*ISYM,ICENT,ISYM*NREM,IFCNT,NFCNT,IWORK,2*NWORK)                                                         !CRA13300
      CALL DFT(DATA,NPREV,NCURR,NREM,ISIGN,IFACT,WORK)                  !CRA13310
  220 IF(IFORM)230,240,260                                              !CRA13320
  230 NREM=NREM*N(IDIM)                                                 !CRA13330
      GO TO 260                                                         !CRA13340
  240 IF(IDIM-1)250,250,260                                             !CRA13350
  250 CALL FIXRL(DATA,N(1),NREM,ISIGN,IFORM)                            !CRA13360
      NTOT=NTOT/N(1)*(N(1)/2+1)                                         !CRA13370
  260 CONTINUE                                                          !CRA13380
      RETURN                                                            !CRA13390
      END                                                               !CRA13400
      SUBROUTINE SYMRV(DATA,NPREV,N,NREM,IFACT,NFACT)                   !CRA14490
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      DIMENSION DATA(*),IFACT(*)                                        !CRA14500
      IF(NFACT-1)80,80,10                                               !CRA14510
   10 IP0=2                                                             !CRA14520
      IP1=IP0*NPREV                                                     !CRA14530
      IP4=IP1*N                                                         !CRA14540
      IP5=IP4*NREM                                                      !CRA14550
      I4REV=1                                                           !CRA14560
      DO 70 I4=1,IP4,IP1                                                !CRA14570
      IF(I4-I4REV)20,40,40                                              !CRA14580
   20 I1MAX=I4+IP1-IP0                                                  !CRA14590
      DO 30 I1=I4,I1MAX,IP0                                             !CRA14600
      DO 30 I5=I1,IP5,IP4                                               !CRA14610
      I5REV=I4REV+I5-I4                                                 !CRA14620
      TEMPR=DATA(I5)                                                    !CRA14630
      TEMPI=DATA(I5+1)                                                  !CRA14640
      DATA(I5)=DATA(I5REV)                                              !CRA14650
      DATA(I5+1)=DATA(I5REV+1)                                          !CRA14660
      DATA(I5REV)=TEMPR                                                 !CRA14670
   30 DATA(I5REV+1)=TEMPI                                               !CRA14680
   40 IP3=IP4                                                           !CRA14690
      DO 60 IF=1,NFACT                                                  !CRA14700
      IP2=IP3/IFACT(IF)                                                 !CRA14710
      I4REV=I4REV+IP2                                                   !CRA14720
      IF(I4REV-IP3)70,70,50                                             !CRA14730
   50 I4REV=I4REV-IP3                                                   !CRA14740
   60 IP3=IP2                                                           !CRA14750
   70 CONTINUE                                                          !CRA14760
   80 RETURN                                                            !CRA14770
      END                                                               !CRA14780
