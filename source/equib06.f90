! For input atomic data file formats see Readme file attached at the end
! **********************************************************************
!          Program EQUIB  (FORTRAN 77)
!
!    Programming history:
!
!    1981 May 3    IDH    Version 1
!    1981 May 5    IDH    Minibug fixed!
!    1981 May 7    IDH    Now takes collision rates or strengths
!    1981 Aug 3    SA     Interpolates collision strengths
!    1981 Aug 7    SA     Input method changed
!    1984 Nov 19   RESC   SA files entombed in scratch disk. Logical
!                         filenames given to SA's data files.
!    1995 Aug      DPR    Changed input file format. Increased matrices.
!    1996 Feb      XWL    Tidy up. SUBROUTINES SPLMAT, HGEN, CFY and CFD
!                         modified such that matrix sizes (i.e. maximum
!                         of Te and maximum no of levels) can now be cha
!                         by modifying the parameters NDIM1, NDIM2 and N
!                         in the Main program. EASY!
!                         Now takes collision rates as well.
!                         All variables are declared explicitly
!                         Generate two extra files (ionpop.lis and ionra
!                         of plain stream format for plotting
!    1996 June     CJP    Changed input data format for cases IBIG=1,2.
!                         Fixed readin bug for IBIG=2 case.
!                         Now reads reformatted upsilons (easier to see
!                         and the 0 0 0 data end is excluded for these c
!                         The A values have a different format for IBIG=
!    2006           BE    Converted to F90
!    2015           RW    Misc updates and improvements
!    2019           RW    allocatable arrays using actual number of
!                         temperatures and levels instead of hard coding
!
! ***** N.B!!  NO TRAPS FOR BAD DATA!!  TAKE CARE!! ****C
!

      USE mod_subroutines
      IMPLICIT NONE
      INTEGER, PARAMETER :: DP = KIND(1.D0)

      INTEGER :: GX, I, I1, I2, J, K, L, KK, LL, JT, JJD,       &
     & NLINES, NLEV, NTEMP, IBIG, IRATS, NTRA, ITEMP,            &
     & IN, NLEV1, KP1, INT, IND, IOPT, IT, IM1, JM1, IP1,               &
     & IAPR, IBPR, ICPR, IKT, IA, IB, IC, IA1, IA2, IB1, IB2, IC1, IC2

      INTEGER,DIMENSION(2) :: ID,JD
      INTEGER,DIMENSION(:),ALLOCATABLE :: G
      INTEGER,DIMENSION(:,:),ALLOCATABLE :: ITRANA,ITRANB,ITRANC

      REAL(KIND=DP) :: TEMPI, TINC, DENSI, DINC, DENS, DLOGD, TEMP, TLOGT,       &
     & TEMP2, DD, DELTEK, EXPE, VALUE, SUMN, TTT, TTP, AHB, EJI, WAV,   &
     & RLINT, FINT, SUMA, SUMB, SUMC, QX, AX, EX, FRAT, DEE
      REAL(KIND=DP),DIMENSION(:),ALLOCATABLE :: N,N2,WAVA,WAVB,WAVC,QQ,E,T,ROOTT,Y,Y2,YKEEP,D,GH
      REAL(KIND=DP),DIMENSION(:,:),ALLOCATABLE :: TDRAT,TNIJ,FINTIJ,CS,QEFF,A,X,X2,XKEEP,HMH
      REAL(KIND=DP),DIMENSION(:,:,:),ALLOCATABLE :: Q,QOM

      CHARACTER(LEN=20) :: ION
      CHARACTER(LEN=20),DIMENSION(:),ALLOCATABLE :: LABEL
      CHARACTER(LEN=1),DIMENSION(78) :: LTEXT

      WRITE(6,1000)                                    !Write out ions available
      ION = '                    '
      WRITE(6,1001)
      READ(5,1002) ION                                  !Interrogate for input

OPEN(UNIT=1,STATUS='OLD',file='/usr/share/equib06/'//trim(ion)//'.dat')
!     + file='atomic_data/'//trim(ion)//'.dat')

      READ(1,*) NLINES                                !Read in no. comment lines
      DO I = 1, NLINES
        READ(1,1003) LTEXT                                             !Comments
!       WRITE(6,1003) LTEXT
      ENDDO
      READ (1,*) NLEV, NTEMP              !Read no. of levels (max=NDIM2) NLEV,

!! allocate arrays

      ALLOCATE(G(NLEV))
      ALLOCATE(N(NLEV))
      ALLOCATE(N2(NLEV))
      ALLOCATE(WAVA(NLEV))
      ALLOCATE(WAVB(NLEV))
      ALLOCATE(WAVC(NLEV))
      ALLOCATE(E(NLEV))
      ALLOCATE(LABEL(NLEV))
      ALLOCATE(D(NTEMP))
      ALLOCATE(QQ(NTEMP))
      ALLOCATE(Y(NTEMP))
      ALLOCATE(Y2(NTEMP))
      ALLOCATE(YKEEP(NTEMP))
      ALLOCATE(T(NTEMP))
      ALLOCATE(ROOTT(NTEMP))

      ALLOCATE(HMH(NTEMP,NTEMP))
      ALLOCATE(ITRANA(2,NLEV))
      ALLOCATE(ITRANB(2,NLEV))
      ALLOCATE(ITRANC(2,NLEV))
      ALLOCATE(TNIJ(NLEV,NLEV))
      ALLOCATE(FINTIJ(NLEV,NLEV))
      ALLOCATE(CS(NLEV,NLEV))
      ALLOCATE(QEFF(NLEV,NLEV))
      ALLOCATE(A(NLEV,NLEV))
      ALLOCATE(X(NLEV,NLEV))
      ALLOCATE(X2(NLEV,NLEV))
      ALLOCATE(XKEEP(NLEV,NLEV))

      ALLOCATE(Q(NTEMP,NLEV,NLEV))
      ALLOCATE(QOM(NTEMP,NLEV,NLEV))

      ALLOCATE(GH(NTEMP*3))

! initialisations

      G=0
      ITRANA=0
      ITRANB=0
      ITRANC=0

      IAPR=0
      IBPR=0
      ICPR=0
      JJD=1
      K=0

      DO I = 1, NLEV                      !no. of Te (max=NDIM1) NTEMP and the
         READ (1,1002) LABEL(I)           !input format (cf Readme)
      ENDDO
!     be
      ibig=0
!     WRITE(6,1004) (I,LABEL(I),I=1,NLEV)   !Tell user what levels are available
      DO I = 1, NTEMP                       !Read in Te's where coll. strengths
           READ (1,*) T(I)                  !are tabulated
           T(I) = LOG10 (T(I))
           ROOTT(I) = SQRT(T(I))
      ENDDO

      READ(1,*) IRATS        !If IRATS=0, what tabulated are collision strengths
!                            !Else Coll. rates = tabulated values * 10 ** IRATS
      IF(IBIG.EQ.0) THEN
   10   READ (1,*) ID(2), JD(2), QX
        IF (QX.EQ.0.0) GOTO 20
        IF (ID(2).EQ.0) THEN
          ID(2) = ID(1)
          K = K + 1
        ELSE
          ID(1) = ID(2)
          K = 1
        ENDIF
        IF (JD(2).EQ.0) THEN
          JD(2) = JD(1)
        ELSE
          JD(1) = JD(2)
        ENDIF
        I = ID(2)
        J = JD(2)
        QOM(K,I,J) = QX
        GO TO 10
      ENDIF
   20 IF(IBIG.EQ.1.OR.IBIG.EQ.2) THEN
        READ(1,*) NTRA
        DO IN = 1, NTRA
          READ(1,*) I,J,(QOM(ITEMP,I,J),ITEMP=1,NTEMP)
        ENDDO
      ENDIF
      NLEV1 = NLEV - 1                            !Read transition probabilities
      IF (IBIG.EQ.1) THEN
       READ(1,7000) ((I,J,A(J,I),L=K+1,NLEV),K=1,NLEV1)
      ELSE
      DO K = 1, NLEV1
        KP1 = K + 1
          DO L = KP1, NLEV
            READ (1,*) I, J, AX
            A(J,I) = AX
          ENDDO
      ENDDO
      ENDIF
      DO J = 1, NLEV             !Read statistical weights, energy levels (cm-1)
        READ (1,*) I, GX, EX
        G(I) = GX
        E(I) = EX
      ENDDO
      CLOSE (UNIT=1)

      WRITE(6,1010)             !Get levels for ratio
      READ(5,*) ((ITRANA(LL,KK),LL=1,2),KK=1,NLEV)          !150 large enough
      WRITE(6,1011)
      READ(5,*) ((ITRANB(LL,KK),LL=1,2),KK=1,NLEV)
      WRITE(6,1012)
      READ(5,*) ((ITRANC(LL,KK),LL=1,2),KK=1,NLEV)
      WRITE(6,1005)                            !Read in Te and Ne where the line
      READ(5,*) TEMPI,TINC,INT                 !ratio is to be calculated
      INT=INT+1
      WRITE(6,1006)
      READ(5,*) DENSI,DINC,IND
      IND=IND+1
!      OPEN(UNIT=3,file='equib_pop.dat')
!      OPEN(UNIT=4,file='equib_rat.dat')
!      OPEN(UNIT=7,file='equib_pop.lis')
!      OPEN(UNIT=8,file='equib_rat.lis')
     OPEN(UNIT=3,file=trim(ion)//'pop.dat',STATUS='UNKNOWN')
     OPEN(UNIT=4,file=trim(ion)//'rat.dat',STATUS='UNKNOWN')
     OPEN(UNIT=7,file=trim(ion)//'pop.lis',STATUS='UNKNOWN')
     OPEN(UNIT=8,file=trim(ion)//'rat.lis',STATUS='UNKNOWN')

      ALLOCATE(TDRAT(2,IND))

      DO JT = 1, INT                                           !Start of Te loop
        TEMP=TEMPI+(JT-1)*TINC
!       IF(TEMPI.LT.30.0) THEN
!         TEMP=10.0**TEMP
!       ENDIF
        DO JJD = 1, IND                                        !Start of Ne loop
          DENS=DENSI+(JJD-1)*DINC
          IF(DENSI.LT.30.0) THEN
            DENS=10.0**DENS
          ENDIF
          IF (TEMP.LE.0.0.OR.DENS.LE.0.0) THEN
            WRITE (6,6100)
            STOP
          ENDIF
          DLOGD = LOG10 (DENS)
          TLOGT = LOG10 (TEMP)
          TEMP2= SQRT (TEMP)

! intialisations
          X=0
          CS=0
          QEFF=0
          TNIJ=0
          Y=0

          IOPT=0
          IF (NTEMP.EQ.1) THEN
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 1 Te only - assuming const'
          ELSEIF (NTEMP.EQ.2) THEN
            WRITE (6,*)
            WRITE (6,*)                                                 &
     &      'Coll. strengths available for 2 Te only - linear interp'
          ELSE
            CALL SPLMAT(T, NTEMP, IOPT, HMH,GH,Y)
            CALL CFD(TLOGT,T,NTEMP,HMH,D)
          ENDIF
          DO I = 2, NLEV
            DO J = I, NLEV
              DELTEK = (E(I-1)-E(J))*1.43884630           !Negative!
              EXPE = EXP(DELTEK/TEMP)
              DO IT = 1, NTEMP
                IF (IRATS.EQ.0.D+00) THEN
                  QQ(IT) = QOM(IT,I-1,J)
                ELSE
                  QQ(IT) = QOM(IT,I-1,J) / EXPE       !Take out the exp. depend.
                ENDIF                                 !before interpolation
              ENDDO
              IF (NTEMP.EQ.1) THEN
                 DD = QQ(1)
              ELSEIF (NTEMP.EQ.2) THEN
                 DD = QQ(1) +                                           &
     &            (QQ(2) - QQ(1))/(T(2) - T(1)) * (TLOGT - T(1))
              ELSE
                CALL CFY(TLOGT, DD, T, QQ, NTEMP, D)
              ENDIF
              IF (IRATS.EQ.0.D+00) THEN
                CS(I-1,J) = DD
              ELSE
                CS(I-1,J) = DD * EXPE
              ENDIF
              IF (IRATS .EQ. 0.D+00) THEN
                QEFF(I-1,J) = 8.63E-06*CS(I-1,J) * EXPE / (G(I-1)*TEMP2)
                QEFF(J,I-1) = 8.63E-06 * CS(I-1,J) / (G(J)*TEMP2)
              ELSE
                QEFF(I-1,J) = CS(I-1,J) * 10. ** IRATS
                QEFF(J,I-1) = G(I-1) * QEFF(I-1,J) / (EXPE * G(J))   !Be careful
              ENDIF                                                  !G integer!
            ENDDO
          ENDDO
          DO I = 2, NLEV
            DO J = 1, NLEV
              IF (J.NE.I) THEN
                X(I,J) = X(I,J) + DENS * QEFF(J,I)
                X(I,I) = X(I,I) - DENS * QEFF(I,J)
                IF (J.GT.I) THEN
                  X(I,J) = X(I,J) + A(J,I)
                ELSE
                  X(I,I) = X(I,I) - A(I,J)
                ENDIF
              ENDIF
            ENDDO
          ENDDO
          DO I = 2, NLEV
            IM1 = I - 1
            VALUE = 0.0 - X(I,1)
            Y(IM1) = VALUE
            Y2(IM1) = VALUE
            YKEEP(IM1) = VALUE
            DO J = 2, NLEV
              JM1 = J - 1
              VALUE = X(I,J)
              X(IM1,JM1) = VALUE
              X2(IM1,JM1) = VALUE
              XKEEP(IM1,JM1) = VALUE
            ENDDO
          ENDDO
          CALL LUSLV(X,Y,NLEV1)         !Solve matrices for populations
          DO I = NLEV, 2, -1
            N(I) = Y(I-1)
          ENDDO
          SUMN = 1.0
          DO I = 2, NLEV
            SUMN = SUMN + N(I)
          ENDDO
          DO I = 2, NLEV
            N(I) = N(I) / SUMN
          ENDDO
          N(1) = 1.0 / SUMN
          WRITE (3,3000) ION,TEMP,TLOGT,DENS,DLOGD                !Output data
          DO I = 1, NLEV
            WRITE (3,3100) I, LABEL(I), N(I)
          ENDDO
          WRITE (3,3200)
          TTT=TEMP*1.0E-4
          TTP=TTT**(-0.870)
          AHB=3.036E-14*TTP            !Eff. recombination coef. of Hb
          DO I = 1, NLEV1
            IP1 = I + 1
            DO J = IP1, NLEV
               IF (A(J,I).NE.0.0) THEN
                 EJI = E(J) - E(I)
                 WAV = 1.E8 / EJI
                 RLINT = A(J,I) * EJI
                 RLINT = RLINT *N(J)
                 TNIJ(I,J)=RLINT
                 FINT=N(J)*A(J,I)*4861.33/(DENS*AHB*WAV)
                 FINTIJ(I,J)=FINT
                 WRITE (3,3300) I,J,WAV,RLINT,FINT
               ENDIF
            ENDDO
          ENDDO
          SUMA=0.0  !Search ITRANA, ITRANB & ITRANC for transitions & sum up
          SUMB=0.0
          SUMC=0.0
          IAPR=0
          IBPR=0
          ICPR=0
          DO IKT = 1, NLEV
            IA1=ITRANA(1,IKT)
            IA2=ITRANA(2,IKT)
            IF(IA1.NE.0.AND.IA2.NE.0) THEN
              SUMA=SUMA+TNIJ(IA1,IA2)
              IAPR=IAPR+1
            ENDIF
            IB1=ITRANB(1,IKT)
            IB2=ITRANB(2,IKT)
            IF(IB1.NE.0.AND.IB2.NE.0) THEN
             IBPR=IBPR+1
             SUMB=SUMB+TNIJ(IB1,IB2)
            ENDIf

            IC1=ITRANC(1,IKT)
            IC2=ITRANC(2,IKT)
            IF(IC1.NE.0.AND.IC2.NE.0) THEN
             ICPR=ICPR+1
             SUMC=SUMC+FINTIJ(IC1,IC2)
            ENDIf
          ENDDO
          FRAT=SUMA/SUMB
          if (sumc.ne.0.d0) SUMC = 1./SUMC
          TDRAT(1,JJD)=DENS
          TDRAT(2,JJD)=FRAT
!          write(6,*),jd,suma,sumb,sumc,dens,frat
          WRITE(7,1017) TEMP, DENS, SUMC
          WRITE(8,1017) TEMP, DENS, FRAT
        ENDDO                                          !End of the Ne loop
        DO IA = 1, IAPR
          I1=ITRANA(1,IA)
          I2=ITRANA(2,IA)
          DEE=E(I2)-E(I1)
          WAVA(IA)=1.E8/DEE
        ENDDO
        DO IB = 1, IBPR
          I1=ITRANB(1,IB)
          I2=ITRANB(2,IB)
          DEE=E(I2)-E(I1)
          WAVB(IB)=1.E8/DEE
        ENDDO
        DO IC = 1, ICPR
          I1=ITRANC(1,IC)
          I2=ITRANC(2,IC)
          DEE=E(I2)-E(I1)
          WAVC(IC)=1.E8/DEE
        ENDDO
        IF(JT.EQ.1) THEN
          WRITE(4,1018) ION,(WAVA(IA),IA=1,IAPR),(WAVB(IB),IB=1,IBPR)
          WRITE(4,*) (TDRAT(1,KK),KK=1,IND) !1013
        ENDIF
        WRITE(4,1014) TEMP,(TDRAT(2,KK),KK=1,IND)
      ENDDO                                              !End of the Te loop
      WRITE(6,1015) (WAVC(IC),IC=1,ICPR)
      CLOSE(UNIT=3)
      CLOSE(UNIT=4)
      CLOSE(UNIT=7)
      CLOSE(UNIT=8)
      STOP
 1000 FORMAT(/,60('*'),/,5X,' Welcome to the EQUIB program',            &
     & /5X,' Data is available for the following ions ( and others) : ',&
     & /,9X,' AlIII, AlV, AlVI, AlVII, AlVIII, AlIX',                   &
     & /,9X,' ArII, ArIII, ArIV, ArV, ArVI, ArX, ArXI, ArXIV',          &
     & /,9X,' CI, CII, CIII, CIV, CaIV, CaVII, CaVIII, CaXII, CaXVI',   &
     & /,9X,' ClIII, ClIX, ClX, CoVI, CoXI, CrIII, CrVIII, CrXVI',      &
     & /,9X,' FeII, FeV, FeVI, FeVII, FeX, FeXVIII, FeXXII',            &
     & /,9X,' FII, FIV, KIII, KVII, KXI',                               &
     & /,9X,' MgII, MgIV, MgV, MgVI, MgVII, MgVIII',                    &
     & /,9X,' MnIV, MnIX, MnXVII, NaIII, NaIV, NaVI, NaVII',            &
     & /,9X,' NeII, NeIII, NeIV, NeV, NeVI, Ni2, Ni7, Ni12',            &
     & /,9X,' NI, NII, NIII, NIV, NV, OI, OII, OIII, OIV, OV',          &
     & /,9X,' PVII, PVIII, PX, PXI, ScV, ScXIII',                       &
     & /,9X,' Si2, Si3, Si4, Si6, Si7, Si8, Si9, Si10',                 &
     & /,9X,' SII, SIII, SIV, SVIII, SIX, SXI, SXII',                   &
     & /,9X,' Ti6, Ti14, V2, V7, V15',                                  &
     & /1x,60('*'),/)
 1001 FORMAT(1X,'Enter name of ion : ',$)
 1002 FORMAT(A20)
 1003 FORMAT(78A1)
! 1004 FORMAT(/1X,'Data tabulated for the following levels ;',/,         &
!     &(2X,I2,1X,A20))
 1005 FORMAT(1X,'Enter initial Temperature, Increment, No of incre',    &
     & 'ments',/1X,' [linear only] : ',$)
 1006 FORMAT(1X,'Enter initial Density, Increment, No of incre',        &
     & 'ments',/1X,' [ log (< 30) or linear ] : ',$)
 1010 FORMAT(1X,'Enter data for line ratio A/B ',/,                     &
     & 1X,'Transitions for line A : ',$)
 1011 FORMAT(1X,'Transitions for line B : ',$)
 1012 FORMAT(1X,'Transitions for which A value is required : ',$)
! 1013 FORMAT(//1X,' TEMPorDENS | ',(1X,1PE9.3),                         &
!     & /1X,' ---------+-',('----------'))
 1014 FORMAT(1H ,1PE9.3,' | ',(1X,1PE9.3))
 1015 FORMAT(/1X,'The A value, N(ion)/N(H+) = A * I(sum)/I(Hb), is',    &
     & /1X,'for the total intensity of the following lines',            &
     & /1X,(F9.1,' ')/)
 1017 FORMAT(3(1X,1PE10.3))
 1018 FORMAT(/1X,A20,' ( ',(F9.1,','),' ) / ( ',                        &
     & (F9.1,','),' )'/)
 3000 FORMAT (//1H ,10X,A20/                                            &
     & 10X,' T =',F9.0,', LOG T =',F7.3/                                &
     & 10X,' D =',1PD9.2,', LOG D =',0PF7.3/                            &
     & ' POPULATIONS:')
 3100 FORMAT (1H ,I4,3X,A20,1PD12.4)
 3200 FORMAT (1H ,'TRANSITION',4X,'LAMDA(A)',5X,'INTENSITY',            &
     & 4X,' I(line)*N(H+)/ I(Hbeta)*N(ion)')
 3300 FORMAT (1H ,I4,I3,5X,F12.2,2(1PD13.3))
 6100 FORMAT (' PROCESSING COMPLETED'/                                  &
     & ' GOODBYE!!'///)
 7000 FORMAT (4(2I4,2X,1PE10.3))
      END
