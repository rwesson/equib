module mod_subroutines
contains
!
!---- PROC LUSLV
      SUBROUTINE LUSLV(A,B,N)                        !Solving linear equations
      IMPLICIT NONE
      INTEGER :: N
      REAL,DIMENSION(:,:) :: A
      REAL,DIMENSION(:) :: B
      CALL LURED(A,N)
      CALL RESLV(A,B,N)
      RETURN
      END
!
!---- PROC LURED
      SUBROUTINE LURED(A,N)
      IMPLICIT NONE
      INTEGER :: N, NM1, I, J, K, IP1
      REAL :: FACT
      REAL,DIMENSION(:,:) :: A
      IF(N.EQ.1) RETURN
      NM1=N-1
      DO I=1,NM1
        IP1=I+1
        DO K=IP1,N
          FACT=A(K,I)/A(I,I)
          DO J=IP1,N
            A(K,J)=A(K,J)-A(I,J)*FACT
          ENDDO
        ENDDO
      ENDDO
      RETURN
      END
!
!---- PROC RESLV
      SUBROUTINE RESLV(A,B,N)                               !Resolve A with B
      IMPLICIT NONE
      INTEGER N, NM1, I, J, K, L, IP1
      REAL,DIMENSION(:) :: B
      REAL,DIMENSION(:,:) :: A
      IF(N.EQ.1) GOTO 1
      NM1=N-1
      DO I=1,NM1
        IP1=I+1
        DO J=IP1,N
          B(J)=B(J)-B(I)*A(J,I)/A(I,I)
        ENDDO
      ENDDO
      B(N)=B(N)/A(N,N)
      DO I=1,NM1
        K=N-I
        L=K+1
        DO J=L,N
          B(K)=B(K)-B(J)*A(K,J)
        ENDDO
        B(K)=B(K)/A(K,K)
      ENDDO
      RETURN
    1 B(N)=B(N)/A(N,N)
      RETURN
      END
!
!---- PROC SPLMAT
      SUBROUTINE SPLMAT(XX,NPT,IOPT,HMH,GH,Y)
      IMPLICIT NONE
      INTEGER :: NPT, IOPT, NPM, NELEM
      REAL,DIMENSION(:) :: XX,GH,Y
      REAL,DIMENSION(:,:) :: HMH
      NPM=NPT-2
      CALL GHGEN(GH,XX,NPT,IOPT)
      NELEM=3*NPM-2
      CALL ELU(GH,NPM)
      CALL HGEN(XX,GH,Y,NPT,IOPT,HMH)
      RETURN
      END
!
!---- PROC DERIV
!     Calculate the first derivative of the lagrangian interpolator
!     of a function F, tabulated at the N points XY(I), I=1 to N.
!     The derivative is given as the coefficients of F(I), I=1 to N,
!     in the array D(I), I=1 to N.
      SUBROUTINE DERIV(XY,D,X,N)
      IMPLICIT NONE
      INTEGER :: N, I, J, K
      REAL :: X, P1, P2, S
      REAL,DIMENSION(:) :: XY,D
      DO I=1,N
        P1=1.
        S=0.
        DO J=1,N
          IF(J.NE.I) THEN
            P1=P1*(XY(I)-XY(J))
            P2=1.
            DO K=1,N
              IF(K.NE.I.AND.K.NE.J) P2=P2*(X-XY(K))
            ENDDO
            S=S+P2
          ENDIF
        ENDDO
        D(I)=S/P1
      ENDDO
      RETURN
      END
!
!---- PROC HGEN
!     Cubic spline interpolation
!     The equation for the second derivatives at internal points
!     is of the form G*YPP=B, where G has been evaluated and LU
!     decomposed.
!     this routine writes B=HMH*Y and then solves YPP=G**(-1)*HMH*Y,
!     =HMH*Y.
!     Three options are provided for boundary conditions-
!     IOPT = 0  YPP=0 at end points
!     IOPT = 1  YP=0  at end points
!     IOPT = 2  YP at end points from lagarnge interpolant of a set of
!     internal points.
      SUBROUTINE HGEN(XX,GH,Y,NPT,IOPT,HMH)
      IMPLICIT NONE
      INTEGER :: NPT, IOPT, NDIM3, NIP, I, J, K, NPM, INDX
      REAL,DIMENSION(:,:) :: HMH
      REAL,DIMENSION(:) :: XX,GH,Y
      REAL,DIMENSION(5) :: XY,D
      REAL,DIMENSION(2,5) :: C
      REAL :: A0, AN1, H1, H2

      IF(IOPT.EQ.2) THEN           !Case of derivative boundary condition, with
        NDIM3=5                    !derivatives from NIP-point Lagrange at
        NIP=3                      !internal points
        DO J=1,2
          DO I=1,NIP
            K=(NPT-NIP)*(J-1)
            XY(I)=XX(K+I)
          ENDDO
          K=1+(NPT-1)*(J-1)
          CALL DERIV(XY,D,XX(K),NIP)
          DO I=1,NIP
            C(J,I)=D(I)
          ENDDO
        ENDDO
      ENDIF
      A0=XX(2)-XX(1)                         !Set up matrix equation G*YPP=HMH*Y
      AN1=XX(NPT)-XX(NPT-1)
      NPM=NPT-2
      DO I=1,NPM
        H1=6./(XX(I+1)-XX(I))
        H2=6./(XX(I+2)-XX(I+1))
        DO J=1,NPT
          HMH(I,J)=0.
          IF(J.EQ.I) HMH(I,J)=H1
          IF(J.EQ.I+2) HMH(I,J)=H2
          IF(J.EQ.I+1) HMH(I,J)=-H1-H2
        ENDDO
      ENDDO
      IF(IOPT.EQ.1.OR.IOPT.EQ.2) THEN            !Correct matrix for case of
        HMH(1,1)=HMH(1,1)+3/A0                   !derivative boundary conditions
        HMH(1,2)=HMH(1,2)-3/A0
        HMH(NPM,NPT-1)=HMH(NPM,NPT-1)-3/AN1
        HMH(NPM,NPT)=HMH(NPM,NPT)+3/AN1
      ENDIF
      IF(IOPT.EQ.2) THEN
        DO J=1,NIP
          HMH(1,J)=HMH(1,J)+3*C(1,J)
          K=NPT+J-NIP
          HMH(NPM,K)=HMH(NPM,K)-3*C(2,J)
        ENDDO
      ENDIF
      DO I=1,NPM
      ENDDO
      DO I=1,NPT                 !Solve matrix equation with results in the form
        Y(1)=HMH(1,I)            !YPP=HMH*Y. matrix g has been LU decomposed
        INDX=0
        DO J=2,NPM
          INDX=INDX+3
          Y(J)=HMH(J,I)-GH(INDX)*Y(J-1)
        ENDDO
        INDX=INDX+1
        Y(NPM)=Y(NPM)/GH(INDX)
        DO J=2,NPM
          K=NPM-J+1
          INDX=INDX-3
          Y(K)=(Y(K)-GH(INDX+1)*Y(K+1))/GH(INDX)
        ENDDO
        DO J=1,NPM
        HMH(J+1,I)=Y(J)
        ENDDO
        HMH(1,I)=0.                 !Insert values for second derivatives at end
        HMH(NPT,I)=0.               !points: first and last rows of the matrix
      ENDDO
      IF(IOPT.GT.0) THEN                 !Case of derivative boundary conditions
        DO J=1,NPT
          HMH(1,J)=-0.5*HMH(2,J)
          HMH(NPT,J)=-0.5*HMH(NPT-1,J)
        ENDDO
        HMH(1,1)=HMH(1,1)-3/(A0*A0)
        HMH(1,2)=HMH(1,2)+3/(A0*A0)
        HMH(NPT,NPT-1)=HMH(NPT,NPT-1)+3/(AN1*AN1)
        HMH(NPT,NPT)=HMH(NPT,NPT)-3/(AN1*AN1)
      ENDIF
      IF(IOPT.EQ.2) THEN
        DO J=1,NIP
          HMH(1,J)=HMH(1,J)-3*C(1,J)/A0
          K=NPT+J-NIP
          HMH(NPT,K)=HMH(NPT,K)+3*C(2,J)/AN1
        ENDDO
      ENDIF
      DO I=1,NPT
      ENDDO
      RETURN
      END
!
!---- PROC GHGEN
      SUBROUTINE GHGEN(GH,XX,NPT,IOPT)
      IMPLICIT NONE
      INTEGER :: NPT, IOPT, INDX, NPTM, I, J, IP, JP, IK
      REAL,DIMENSION(:) :: XX,GH
      INDX=0
      NPTM=NPT-1
      DO I=2,NPTM
        IP=I-1
        DO J=1,3
          JP=IP+J-2
          IF(JP.GE.1.AND.JP.LE.NPTM-1) THEN
            INDX=INDX+1
            IF(J.EQ.2) THEN
              GH(INDX)=2*(XX(I+1)-XX(I-1))
            ELSE
              IK=I+(J-1)/2
              GH(INDX)=XX(IK)-XX(IK-1)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      IF(IOPT.GE.1) THEN
        GH(1)=GH(1)-(XX(2)-XX(1))/2.
        GH(INDX)=GH(INDX)-(XX(NPT)-XX(NPT-1))/2.
      ENDIF
      RETURN
      END
!
!---- PROC ELU
      SUBROUTINE ELU(GH,N)
      IMPLICIT NONE
      INTEGER :: N, INDX, I, J, JP
      REAL,DIMENSION(:) :: GH
      INDX=0
      DO I=1,N
        DO J=1,3
          JP=I+J-2
          IF(JP.GE.1.AND.JP.LE.N) THEN
            INDX=INDX+1
            IF(I.GT.1) THEN
              IF(J.EQ.1) THEN
                GH(INDX)=GH(INDX)/GH(INDX-2)
              ENDIF
              IF(J.EQ.2) THEN
                GH(INDX)=GH(INDX)-GH(INDX-1)*GH(INDX-2)
              ENDIF
            ENDIF
          ENDIF
        ENDDO
      ENDDO
      RETURN
      END
!
!---- PROC CFY
      SUBROUTINE CFY(X,Y,XX,YY,NPT,D)
      IMPLICIT NONE
      INTEGER :: NPT, J
      REAL :: X, Y, TT
      REAL,DIMENSION(:) :: XX,YY,D

      IF(X.LT.XX(1)) THEN
        Y=YY(1)
        WRITE(6,400) XX(1)
      ENDIF
      IF(X.GT.XX(NPT)) THEN
        Y=YY(NPT)
        WRITE(6,401) XX(NPT)
      ENDIF
      TT=0.
      DO J=1,NPT
        TT=TT+D(J)*YY(J)
      ENDDO
      Y=TT
      RETURN
  400 FORMAT(2X,'TEMPERATURE CHOSEN IS BELOW THE AVAILABLE RANGE,',/,   &
     &2X,'THE FIRST TEMP-COLLISION STRENGTH IS BEING USED,',/,          &
     &2X,'FOR T=',F8.4)
  401 FORMAT(2X,'TEMPERATURE CHOSEN IS ABOVE THE AVAILABLE RANGE,',/,   &
     &2X,'THE HIGHEST TEMP-COLLISION STRENGTH IS BEING USED,',/,        &
     &2X,'FOR T=',F8.4)
      END
!
!---- PROC CFD
      SUBROUTINE CFD(X,XX,NPT, HMH, D)
      IMPLICIT NONE
      INTEGER :: NPT, NPTM, I, J
      REAL :: X, X1, X2, A1, A2, HI
      REAL,DIMENSION(:) :: XX, D
      REAL,DIMENSION(:,:) :: HMH
      IF(X.LT.XX(1)) THEN
        WRITE(6,400) XX(1)
        RETURN
      ENDIF
      IF(X.GT.XX(NPT)) THEN
        WRITE(6,401) XX(NPT)
        RETURN
      ENDIF
      NPTM=NPT-1
      DO I=1,NPTM
        IF(X.LT.XX(I+1)) THEN
          X1=XX(I+1)-X
          X2=X-XX(I)
          HI=XX(I+1)-XX(I)
          A1=X1*(X1*X1/(6*HI)-HI/6)
          A2=X2*(X2*X2/(6*HI)-HI/6)
          DO J=1,NPT
            D(J)=(A1*HMH(I,J)+A2*HMH(I+1,J))
          ENDDO
          D(I)=D(I)+X1/HI
          D(I+1)=D(I+1)+X2/HI
          RETURN
        ENDIF
      ENDDO
  400 FORMAT(2X,'TEMPERATURE CHOSEN IS BELOW THE AVAILABLE RANGE,',/,   &
     &2X,'THE FIRST TEMP-COLLISION STRENGTH IS BEING USED,',/,          &
     &2X,'FOR T=',F8.4)
  401 FORMAT(2X,'TEMPERATURE CHOSEN IS ABOVE THE AVAILABLE RANGE,',/,   &
     &2X,'THE HIGHEST TEMP-COLLISION STRENGTH IS BEING USED,',/,        &
     &2X,'FOR T=',F8.4)
      END
end module mod_subroutines