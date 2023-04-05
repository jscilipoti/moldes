!==========================================================================
      SUBROUTINE GENESTINT (IFAM,MDV,NUD,J)
!C
      PARAMETER (NSVA=20,NMAX=60000)
      IMPLICIT real*8 (A-H,O-Z)
      INTEGER NUD(NSVA,NMAX)
!c
!c		NUD es la matriz donde se va almacenando la composición 
!c		por subgrupos de cada nueva estructura intermedia generada
!C
!C     GENERATION OF INTERMEDIATE STRUCTURES
!C
      DO 101 J=1,NMAX
         DO 101 I=1,MDV
            NUD(I,J)=0
101   CONTINUE
      J=0
!C
!C     GENERATION OF ONE GROUP STRUCTURES
!C
      IF ((IFAM.EQ.0).OR.(IFAM.EQ.2).OR.(IFAM.EQ.3).OR.(IFAM.EQ.4)) THEN
         DO 130 I1=1,MDV
            J=J+1
            NUD(I1,J)=1
130      CONTINUE
      END IF
!C
!C     GENERATION OF TWO GROUP STRUCTURES
!C
      IF ((IFAM.EQ.0).OR.(IFAM.EQ.3).OR.(IFAM.EQ.4)) THEN
         DO 140 I1=1,MDV
            DO 140 I2=I1,MDV
               J=J+1
               NUD(I1,J)=1
               NUD(I2,J)=NUD(I2,J)+1
140      CONTINUE
      END IF
!C
!C        GENERATION OF THREE GROUP STRUCTURES
!C
      IF ((IFAM.EQ.3).OR.(IFAM.EQ.4)) THEN
         DO 150 I1=1,MDV
            DO 150 I2=I1,MDV
               DO 150 I3=I2,MDV
                  J=J+1
                  NUD(I1,J)=1
                  NUD(I2,J)=NUD(I2,J)+1
                  NUD(I3,J)=NUD(I3,J)+1
150      CONTINUE
      END IF
      IF ((IFAM.EQ.3).OR.(IFAM.EQ.4).OR.(IFAM.EQ.5)) THEN
!C   
!C        GENERATION OF FOUR GROUP STRUCTURES
!C
         DO 160 I1=1,MDV
            DO 160 I2=I1,MDV
               DO 160 I3=I2,MDV
                  DO 160 I4=I3,MDV
                     J=J+1
                     NUD(I1,J)=1
                     NUD(I2,J)=NUD(I2,J)+1
                     NUD(I3,J)=NUD(I3,J)+1
                     NUD(I4,J)=NUD(I4,J)+1
160      CONTINUE
!C
!C        GENERATION OF FIVE GROUP STRUCTURES
!C
         INIP = J + 1
         DO 170 I1=1,MDV
            DO 170 I2=I1,MDV
               DO 170 I3=I2,MDV
                  DO 170 I4=I3,MDV
                     DO 170 I5=I4,MDV
                        J=J+1
                        NUD(I1,J)=1
                        NUD(I2,J)=NUD(I2,J)+1
                        NUD(I3,J)=NUD(I3,J)+1
                        NUD(I4,J)=NUD(I4,J)+1
                        NUD(I5,J)=NUD(I5,J)+1
170      CONTINUE
         INIF = J
      END IF
      IF ((IFAM.EQ.1).OR.(IFAM.EQ.3).OR.(IFAM.EQ.4).OR.(IFAM.EQ.5)) THEN
!C
!C        GENERATION OF SIX GROUP STRUCTURES
!C
         DO 180 I1=1,MDV
            DO 180 I2=I1,MDV
               DO 180 I3=I2,MDV
                  DO 180 I4=I3,MDV
                     DO 180 I5=I4,MDV
                        DO 180 I6=I5,MDV
                           J=J+1
                           NUD(I1,J)=1
                           NUD(I2,J)=NUD(I2,J)+1
                           NUD(I3,J)=NUD(I3,J)+1
                           NUD(I4,J)=NUD(I4,J)+1
                           NUD(I5,J)=NUD(I5,J)+1
                           NUD(I6,J)=NUD(I6,J)+1
180      CONTINUE
      END IF
!c
      RETURN
      END
      
!==========================================================================
      subroutine alif_fact (ngdv,ngsdv,mdv1,mdv2,nud,jisl,ims_test)
!---------------------------------------------------
!     Esta subrutina elimina todas las estructuras 
!     alifáticas intermedias no factibles contenidas
!     en el vector nud. 
!     Variables de entrada:
!       -ngdv: vector que contiene los grupos interme-
!       dios seleccionados.
!       -ngsdv: vector con la información de los 
!       tipos de enlace de cada grupo guardado en el
!       vector ngdv.
!       -mdv1: cantidad de grupos aromáticos en el 
!       vector ngdv.
!       -mdv2: cantidad de grupos alifáticos en el 
!       vector ngdv.
!     Variables de entrada/salida:
!       -nud:vector que contiene todas las estructu-
!       ras alifáticas intermedias.
!       -jisl: cantidad de estructuras en el vector 
!       nud. 
!---------------------------------------------------
! Definición de variables      
      implicit none
      integer, parameter::NSVA=20, NMAX=60000, na1=150
      integer, intent (in):: ngdv(0:na1),ngsdv(NSVA,5),mdv1,mdv2
      integer, intent (out):: ims_test(2,NMAX)
      integer, intent (inout):: jisl
      integer,intent (inout):: nud(NSVA,NMAX)
      integer i,j,jtest,ktest,ngrup
!     
      i=0
      do while (i<jisl)
        i=i+1
        ngrup=0
        jtest=0
        ktest=0
        do j=1,mdv2
            if(nud(j,i)==0)cycle
            jtest = nud(j,i)*ngsdv(j+mdv1,2)+jtest
			ktest = nud(j,i)*ngsdv(j+mdv1,3)+ktest
			ngrup=ngrup+nud(j,i)
        enddo
        ims_test(1,i)=jtest
        ims_test(2,i)=ktest
        if(ngrup==1)cycle
        if(ktest>2)then
            if((2*ktest)>(jtest+2)) call eliminar_estructura(nud,jisl,i)
        else
            if(ktest>jtest) call eliminar_estructura(nud,jisl,i)
        endif
      enddo
!      
      endsubroutine
      
!==========================================================================
     SUBROUTINE AROMI (IPAREQ,MOP,ISF,NGSDV1,SOLDES,natot)
!C-----------------------------------------------------------------------
!C     Esta subrutina realiza el test de factibilidad y controla los casos
!C     especiales de estructuras aromaticas intermedias
!C-----------------------------------------------------------------------
      PARAMETER (NMOD=3,NCOM=3,NSCM=10,NSVA=20)
      implicit real*8 (A-H,O-Z)
      LOGICAL SOLDES,INV
      DIMENSION NGSDV1(NSVA,5)
      INTEGER HTEST,ITEST,natot
      COMMON/GRUPAR/IACH(NMOD),IACCH2(NMOD),IACCL(NMOD),IACNH2(NMOD),IACNO2(NMOD)
      COMMON/GRUPES/ICH(NMOD),IACOH(NMOD),IACCH(NMOD),IAC(NMOD),&
                   IOH(NMOD),ICOOH(NMOD),IACCOO(NMOD)
      COMMON/US/MS(NCOM,NSCM,2)
      COMMON/EVAL/INV
!C
!C
      ITEST=0.0
      HTEST=0.0
      DO 36 I=1,ISF
         ITEST=MS(3,I,2)*NGSDV1(I,4)+ITEST !cantidad de enlaces I
         HTEST=MS(3,I,2)*NGSDV1(I,5)+HTEST !cantidad de enlaces H
  36  CONTINUE
      IF (MOP.EQ.1) THEN !extracción líquido-líquido
         !IF ((ITEST.LE.4).OR.(HTEST.GT.2)) THEN
         IF (ITEST.LT.3) THEN
            SOLDES = .TRUE.
         END IF
      ELSE
         IF ((ITEST.LT.3).OR.(HTEST.GT.3)) THEN
            SOLDES = .TRUE.
         END IF
      END IF
      IF (.NOT.(SOLDES)) THEN !Restricciones operativas
         NACOH = 0	!número del grupo ACOH
         NAC = 0	!número del grupo AC
         NACCH2 = 0	!número del grupo ACCH2
         NACCH = 0	!número del grupo ACCH
         DO 37 I=1,ISF
            IF (MS(3,I,1).EQ.IACOH(IPAREQ)) THEN
               NACOH = NACOH + MS(3,I,2)
            ELSE IF (MS(3,I,1).EQ.IAC(IPAREQ)) THEN
               NAC = NAC + MS(3,I,2)
            ELSE IF (MS(3,I,1).EQ.IACCH2(IPAREQ)) THEN
               NACCH2 = NACCH2 + MS(3,I,2)
            ELSE IF (MS(3,I,1).EQ.IACCH(IPAREQ)) THEN
               NACCH = NACCH + MS(3,I,2)
            ELSE IF (MS(3,I,1).EQ.IACCOO(IPAREQ)) THEN
               NACCOO = NACCOO + MS(3,I,2)
            END IF
  37     CONTINUE
        NATOT = NAC + NACCH2 + NACCOO + 2*NACCH
        IF (.NOT.INV) THEN
            if((MOP.EQ.1).AND.(NACOH.GT.1))then
                SOLDES = .True.
            elseif((MOP.EQ.2).AND.(NACOH.GT.2))then
                SOLDES = .TRUE.
            elseif((NAC.GT.2))then
                SOLDES = .TRUE.
            elseif((NACCH2.GT.2))then
                SOLDES = .TRUE.
            elseif(((NACCH2.GT.0).OR.(NAC.GT.0)).AND.(NACCH.GT.0))then
                SOLDES = .TRUE.
            elseif((NACCH.GT.1))then
                SOLDES = .TRUE.
            elseif((NACCH2.GT.1).AND.(NAC.EQ.1))then
                SOLDES = .TRUE.
            elseif((NACCH2.EQ.1).AND.(NAC.GT.1))then
                SOLDES = .TRUE.
            endif
        ENDIF
      END IF
      RETURN
      endsubroutine

!==========================================================================
    SUBROUTINE COMB_TERM (MSV,NUS,KST,JST,IST)
!C
      PARAMETER (NMAX=60000,NSVA=20)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION NUS(NSVA,NMAX)
	K = 0
!En el siguiente bucle se generan todas las combinaciones pares entre los MSV
!grupos terminales seleccionados
	DO 1510 I1=1,MSV !MSV=cant de subgrupos terminales seleccionados
		DO 1510 I2=I1,MSV
			K = K + 1
			NUS(I1,K) = 1
			NUS(I2,K)= NUS(I2,K) + 1
1510	CONTINUE
	KST = K !cantidad de combinación de pares posibles
!En el siguiente bucle se generan todas las combinaciones ternarias entre los MSV
!grupos terminales seleccionados
      DO 1520 I1=1,MSV
          DO 1520 I2=I1,MSV
              DO 1520 I3=I2,MSV
              K = K + 1
              NUS(I1,K)=1
              NUS(I2,K)=NUS(I2,K)+1
              NUS(I3,K)=NUS(I3,K)+1
1520  CONTINUE
	JST = K !cantidad de combinaciones ternarias posibles
!En el siguiente bucle se generan todas las combinaciones cuaternarias entre los MSV
!grupos terminales seleccionados
      DO 1530 I1=1,MSV
          DO 1530 I2=I1,MSV
              DO 1530 I3=I2,MSV
                  DO 1530 I4=I3,MSV
                  K = K + 1
                  NUS(I1,K)=1
                  NUS(I2,K)=NUS(I2,K)+1
                  NUS(I3,K)=NUS(I3,K)+1
                  NUS(I4,K)=NUS(I4,K)+1
1530  CONTINUE
      IST = K
      RETURN
      END

!==========================================================================
      SUBROUTINE COMB_INTER (JISL,NUR,KSR,JSR,ISR,HSR)

      PARAMETER (NMAX=60000,NA=150)
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION NUR(NA,NMAX)
      INTEGER HSR

	K = 0
      DO 130 I1=1,JISL
          K=K+1
          NUR(I1,K)=1
130   CONTINUE
      KSR = K !Cant de combinaciones unitarias posibles
	DO 1510 I1=1,JISL !JISL=cant de estructuras intermedias generadas
		DO 1510 I2=I1,JISL
			K = K + 1
			NUR(I1,K) = 1
			NUR(I2,K)= NUR(I2,K) + 1
1510	CONTINUE
	JSR = K !cantidad de combinación de pares posibles
      DO 1520 I1=1,JISL
          DO 1520 I2=I1,JISL
              DO 1520 I3=I2,JISL
              K = K + 1
              NUR(I1,K)=1
              NUR(I2,K)=NUR(I2,K)+1
              NUR(I3,K)=NUR(I3,K)+1
1520  CONTINUE
	ISR = K !cantidad de combinaciones ternarias posibles
      DO 1530 I1=1,JISL
          DO 1530 I2=I1,JISL
              DO 1530 I3=I2,JISL
                  DO 1530 I4=I3,JISL
                  K = K + 1
                  NUR(I1,K)=1
                  NUR(I2,K)=NUR(I2,K)+1
                  NUR(I3,K)=NUR(I3,K)+1
                  NUR(I4,K)=NUR(I4,K)+1
1530  CONTINUE
	HSR = K !cantidad de combinaciones cuaternarias posibles
      RETURN
      END

!========================================================================
      SUBROUTINE AROMN (IPAREQ,ISF,NAC,NACCH,NACCH2,NATOT)
!C-----------------------------------------------------------------------
!C     Esta subrutina determina el numero de grupos terminales que deben 
!C     agregarse a una estructura aromatica intermedia.
!C-----------------------------------------------------------------------

      PARAMETER (NMOD=3,NCOM=3,NSCM=10)
      implicit real*8 (A-H,O-Z)
      COMMON/GRUPAR/ IACH(NMOD),IACCH2(NMOD),IACCL(NMOD),IACNH2(NMOD),IACNO2(NMOD)
      COMMON/GRUPES/ ICH(NMOD),IACOH(NMOD),IACCH(NMOD),IAC(NMOD),&
                    IOH(NMOD),ICOOH(NMOD),IACCOO(NMOD)
      COMMON/IUS/IMS(NSCM,2)
      NAC = 0
      NACC = 0
      NACCH2 = 0
      NACCH = 0
      NACCOO = 0
      DO 47 I=1,ISF
         IF (IMS(I,1).EQ.IAC(IPAREQ)) THEN
            NAC = NAC + IMS(I,2)
!            NACC = NACC + IMS(I,2)
         ELSE IF (IMS(I,1).EQ.IACCH2(IPAREQ)) THEN
            NACCH2 = NACCH2 + IMS(I,2)
!            NACC = NACC + IMS(I,2)
         ELSE IF (IMS(I,1).EQ.IACCH(IPAREQ)) THEN
            NACCH = NACCH + IMS(I,2)
!            NACC = 
         ELSE IF (IMS(I,1).EQ.IACCOO(IPAREQ)) THEN
            NACCOO = NACCOO + IMS(I,2)
!            NACC = NACC + IMS(I,2)
         END IF
!         IF ((NAC.EQ.4).OR.(NACC.EQ.4).OR.(NACCH2.EQ.4).OR.(NACCH.EQ.1)
!     +   .OR.(NACCOO.EQ.2))THEN
!            GO TO 48
!         END IF
  47  CONTINUE
  48  NATOT = NAC + NACCH2 + NACCOO + 2*NACCH
      RETURN
      END

!===========================================================================
	SUBROUTINE Ingresar_EstInter (MDV,MS,NINTER,ISF,NGDV,NUM,iCH2)
!---------------------------------------------------------------------------
!	Esta subrutina ingresa NUM grupos NTERM a la IMS de modo tal que
!	su estructura de grupos en MS(3,_,1) y MS(3,_,2) quede ordenada	
!	en forma ascendente según el número de grupo y pueda luego ser
!	comparada con la base de datos PROP.PDB para encontrar isómeros.
!	Para ello primero se define el LUGAR dentro del vector donde 
!	debe ingresarse el grupo terminal. Luego "se corren" un lugar 
!	todos los grupos que se encontraban desde esa posición hasta la
!	 última (ISF) y finalmente se ingresa el grupo terminal.
!---------------------------------------------------------------------------
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (NCOM=3,NSCM=10,NA=150)
	DIMENSION MS(NSCM,2),NINTER(20),NGDV(0:NA),GrupAgre(20)
	logical bool
	LUGAR=ISF
	NINTER(:)=NINTER(:)*NUM
!Las siguientes tres variables se utilizarán para averiguar si el grupo CH2 se agregó solo o 
!acompañado por otro(s) grupo(s).
	k=0
	iCH2=0
	nCH2=0
	DO I=1,20
	  IF (NINTER(I).EQ.0) cycle
!Averigua si el grupo ya fue cargado previamente
        bool=.false.
        j=1
        do while (MS(j,1).ne.0)
            if(MS(j,1).ne.NGDV(I+MDV)) then
                j=j+1
                cycle
            endif
            MS(j,2)=MS(j,2)+NINTER(I)
            bool=.true.
            k=k+(1*NUM)
            if(NGDV(I+MDV).eq.2)nCH2=nCH2+(1*NUM)
            exit
        enddo
        if(bool)cycle
	  LUGAR=LUGAR+1
	  MS(LUGAR,1)=NGDV(I+MDV)
	  MS(LUGAR,2)=NINTER(I)
	  k=k+(1*NUM)
        if(NGDV(I+MDV).eq.2)nCH2=nCH2+(1*NUM)
      ENDDO
      CALL ORDENAR (MS)
      if(nCH2.eq.k)iCH2=1
 	RETURN 
	endsubroutine
	
!==========================================================================
   SUBROUTINE INGRESAR_TERMINAL (MS,ISF,NTERM,NUM)
!--------------------------------------------------------------------------
!	Esta subrutina ingresa NUM grupos NTERM a la IMS de modo tal que
!	su estructura de grupos en MS(3,_,1) y MS(3,_,2) quede ordenada	
!	en forma ascendente según el número de grupo y pueda luego ser
!	comparada con la base de datos PROP.PDB para encontrar isómeros.
!	Para ello primero se define el LUGAR dentro del vector donde 
!	debe ingresarse el grupo terminal. Luego "se corren" un lugar 
!	todos los grupos que se encontraban desde esa posición hasta la
!	última (ISF) y finalmente se ingresa el grupo terminal.
!--------------------------------------------------------------------------
      PARAMETER (NCOM=3,NSCM=10)
	DIMENSION MS(NCOM,NSCM,2)
	LUGAR=ISF+1
      IF (MS(3,1,1).eq.0) THEN
		LUGAR=1
	ELSE
		DO 11 N=1,ISF
			IF (MS(3,N,1).GT.NTERM) THEN
				LUGAR=N
				GOTO 12
			END IF
  11		CONTINUE
	GOTO 14 ! cuando LUGAR=ISF+1
  12		DO 13 I=ISF,LUGAR,-1
			MS(3,I+1,1)=MS(3,I,1)
			MS(3,I+1,2)=MS(3,I,2)
  13		CONTINUE
      END IF
  14  MS(3,LUGAR,1)=NTERM
	MS(3,LUGAR,2)=NUM
 	RETURN 
	endsubroutine

!==========================================================================
    SUBROUTINE AROM1 (IPAREQ,IST,NACCH,NACCH2,SOLDES)
!C----------------------------------------------------------------------
!C     Esta subrutina controla la existencia de casos especiales para
!C     compuestos aromaticos.
!C----------------------------------------------------------------------

      PARAMETER (NMOD=3,NCOM=3,NSCM=10)
      implicit real*8 (A-H,O-Z)
      LOGICAL SOLDES
      COMMON/GRUPES/ ICH(NMOD),IACOH(NMOD),IACCH(NMOD),IAC(NMOD),&
                    IOH(NMOD),ICOOH(NMOD),IACCOO(NMOD)
      COMMON/GRUPAL/ICH3(NMOD),ICH2(NMOD)
      COMMON/US/MS(NCOM,NSCM,2)

      N=0
	DO 7 I=1,IST
	  IF (MS(3,I,1).EQ.ICH3(IPAREQ).OR.MS(3,I,1).EQ.IOH(IPAREQ)&
         .OR.MS(3,I,1).EQ.ICOOH(IPAREQ)) THEN
	  N = N+1
      END IF
   7	CONTINUE   
	IF (N.GT.2*NACCH+NACCH2) THEN
          SOLDES = .TRUE.
      END IF
      RETURN
   endsubroutine

!==========================================================================   
   SUBROUTINE AROM2 (K,NACCH,NSV,NSV2,NGSSV1,SOLDES)
!C----------------------------------------------------------------------
!C     Esta subrutina controla la existencia de casos especiales para
!C     compuestos aromaticos.
!C----------------------------------------------------------------------

      PARAMETER (NMAX=60000,NSVA=20)
      implicit real*8 (A-H,O-Z)
      DIMENSION NSV2(3,NMAX),NGSSV1(NSVA,5)
      LOGICAL SOLDES
      IF (NACCH.EQ.1) THEN
        JTEST = 2
        MTEST = 0
        KTEST = 0
        MSST=0
        KSST=0
	  DO 52 I=1,NSV
            MSST=NSV2(I,K)*NGSSV1(I,1)+MSST
            KSST=NSV2(I,K)*NGSSV1(I,3)+KSST
  52    CONTINUE
        MS1=MSST+MTEST
        KS1=KSST+KTEST
        MJT=MS1+JTEST/2
        IF (KS1.GT.MJT) THEN
            SOLDES = .TRUE.
        END IF
      END IF
      RETURN
   endsubroutine

!==========================================================================   
	SUBROUTINE eliminar_repetidos(nact,TERM,isfi,IST,rep)
    PARAMETER (NSCM=10,NMAX=60000)
	LOGICAL rep
	integer TERM(NMAX,NSCM,2)
	rep = .FALSE.
	do i=1,nact-1
		do j=isfi+1,IST	! j=No.de grupo terminal en la term. nact
			do k=isfi+1,IST	! k=No.de grupo terminal en la term. i
				if(TERM(nact,j,1).eq.TERM(i,k,1)) then
     					if(TERM(nact,j,2).eq.TERM(i,k,2)) go to 225
					go to 227
				end if
 			end do
			go to 227
 225		end do
		rep = .TRUE.
		go to 229
 227	end do
 229	RETURN
      endsubroutine

!==========================================================================         
      SUBROUTINE AROMATICOS_ESPECIALES (ISUBFAM,IPAREQ,MS,IST)
!--------------------------------------------------------------------------
!C	Cuando ifam=4 esta subrutina agrega la parte aromática a la 
!c	IMS alifática ya generada por GENESTINT.
!--------------------------------------------------------------------------
      PARAMETER (NMOD=3,NCOM=3,NSCM=10)
      implicit real*8 (A-H,O-Z)
      DIMENSION MS(NCOM,NSCM,2),IACH(NMOD),IACCH2(NMOD),IACCL(NMOD)
      DIMENSION IACNH2(NMOD),IACNO2(NMOD)
      COMMON/GRUPAR/IACH,IACCH2,IACCL,IACNH2,IACNO2
!C
      IF (ISUBFAM.EQ.1)THEN
!c            MS(3,ISF+2,1)=IACH(IPAREQ)
!c            MS(3,ISF+2,2)=5
!c            MS(3,ISF+3,1)=IACCH2(IPAREQ)
!c            MS(3,ISF+3,2)=1
            CALL INGRESAR_TERMINAL (MS,IST,IACH(IPAREQ),5)
            CALL INGRESAR_TERMINAL (MS,IST+1,IACCH2(IPAREQ),1)
            IST=IST+2
      ELSE IF(ISUBFAM.EQ.2)THEN
!C            MS(3,ISF+2,1)=IACH(IPAREQ)
!C            MS(3,ISF+2,2)=4
!C            MS(3,ISF+3,1)=IACCH2(IPAREQ)
!C            MS(3,ISF+3,2)=1
!C            MS(3,ISF+4,1)=IACCL(IPAREQ)
!C            MS(3,ISF+4,2)=1
            CALL INGRESAR_TERMINAL (MS,IST,IACH(IPAREQ),4)
            CALL INGRESAR_TERMINAL (MS,IST+1,IACCH2(IPAREQ),1)
            CALL INGRESAR_TERMINAL (MS,IST+2,IACCL(IPAREQ),1)
            IST=IST+3
      ELSE IF(ISUBFAM.EQ.3)THEN
!C            MS(3,ISF+2,1)=IACH(IPAREQ)
!C            MS(3,ISF+2,2)=4
!C            MS(3,ISF+3,1)=IACCH2(IPAREQ)
!C            MS(3,ISF+3,2)=1
!C            MS(3,ISF+4,1)=IACNH2(IPAREQ)
!C            MS(3,ISF+4,2)=1
            CALL INGRESAR_TERMINAL (MS,IST,IACH(IPAREQ),4)
            CALL INGRESAR_TERMINAL (MS,IST+1,IACCH2(IPAREQ),1)
            CALL INGRESAR_TERMINAL (MS,IST+2,IACNH2(IPAREQ),1)
            IST=IST+3
      ELSE IF(ISUBFAM.EQ.4)THEN
!C            MS(3,ISF+2,1)=IACH(IPAREQ)
!C            MS(3,ISF+2,2)=4
!C            MS(3,ISF+3,1)=IACCH2(IPAREQ)
!C            MS(3,ISF+3,2)=1
!C            MS(3,ISF+4,1)=IACNO2(IPAREQ)
!C            MS(3,ISF+4,2)=1
            CALL INGRESAR_TERMINAL (MS,IST,IACH(IPAREQ),4)
            CALL INGRESAR_TERMINAL (MS,IST+1,IACCH2(IPAREQ),1)
            CALL INGRESAR_TERMINAL (MS,IST+2,IACNO2(IPAREQ),1)
            IST=IST+3
      END IF
      RETURN
  endsubroutine
  
!=======================================================================  
  SUBROUTINE ALIF2 (IPAREQ,SOLDES)
!C----------------------------------------------------------------------
!C     Esta subrutina controla la existencia de casos especiales para
!C     compuestos alifaticos.
!C----------------------------------------------------------------------
!C
      PARAMETER (NMOD=3,NCOM=3,NSCM=10)
      implicit real*8 (A-H,O-Z)
      LOGICAL SOLDES
      COMMON/GRUPAL/ICH3(NMOD),ICH2(NMOD)
      COMMON/US/MS(NCOM,NSCM,2)
      IF (MS(3,1,1).EQ.ICH2(IPAREQ).AND.MS(3,1,2).EQ.2.AND.&
         MS(3,2,1).EQ.ICH3(IPAREQ).AND.MS(3,2,2).EQ.2) THEN
         SOLDES = .TRUE.
      END IF
      RETURN
      endsubroutine

!UTILIZADAS POR SUBROUTINES DE ESTE ARCHIVO
!==============================================================
      subroutine eliminar_estructura(nud,jisl,i)
!Deficinión de variables
      implicit none
      integer, parameter:: NSVA=20,NMAX=60000
      integer,intent(inout)::jisl,i
      integer,intent(inout)::nud(NSVA,NMAX)
      integer j
!
      do j=i,jisl-1
        nud(:,j)=nud(:,j+1)
      enddo
      nud(:,jisl)=0
      jisl=jisl-1
      i=i-1
!
      endsubroutine







!      SUBROUTINE CICLII (IPAREQ,MOP,ISF,NGSDV1,SOLDES)
!C-----------------------------------------------------------------------
!C     Esta subrutina controla los casos especiales de estructuras ciclicas
!C     intermedias (finales, ya que no hay ramificaciones)
!C-----------------------------------------------------------------------
!      PARAMETER (NMOD=3,NCOM=3,NSCM=10,NSVA=20)
!      implicit real*8 (A-H,O-Z)
!      LOGICAL SOLDES
!      DIMENSION NGSDV1(NSVA,5)
!      COMMON/GRUPES/ICH(NMOD),IACOH(NMOD),IACCH(NMOD),IAC(NMOD),
!     +              IOH(NMOD),ICOOH(NMOD),IACCOO(NMOD)
!      COMMON/GRUPAL/ICH3(NMOD),ICH2(NMOD)
!      COMMON/US/MS(NCOM,NSCM,2)
!C
!C
!      KTEST=0
!      MTEST=0
!      JTEST=0
!      DO 38 I=1,ISF
!         MTEST = MS(3,I,2)*NGSDV1(I,1)+MTEST
!         JTEST = MS(3,I,2)*NGSDV1(I,2)+JTEST
!         KTEST = MS(3,I,2)*NGSDV1(I,3)+KTEST
!  38  CONTINUE
!      MJI = MTEST + JTEST/2
!      IF (KTEST.GT.MJI) THEN
!         SOLDES = .TRUE.
!      END IF
!      IF (.NOT.(SOLDES)) THEN
!         NCH2 = 0
!         NCH = 0
!         DO 37 I=1,ISF
!            IF (MS(3,I,1).EQ.ICH2(IPAREQ)) THEN
!               NCH2 = NCH2 + MS(3,I,2)
!            ELSE IF (MS(3,I,1).EQ.ICH(IPAREQ)) THEN
!               NCH = NCH + MS(3,I,2)
!            END IF
!            IF (((MOP.EQ.1).AND.(NCH2.EQ.3).AND.(NCH.EQ.1)).OR.
!     +          (NCH.GT.1)) 
!     +      THEN
!               SOLDES = .TRUE.
!               GO TO 39
!            END IF
!  37     CONTINUE
!      END IF
!  39  RETURN
!      END


