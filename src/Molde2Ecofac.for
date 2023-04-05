      SUBROUTINE ESCRIBIR_RESULTADOS(MOP,JIST,IS,NSOL,MST,ISOL,
     *                               RDER,TITULO,IDEV,OPT)
!-----------------------------------------------------------------------
!
!-----VARIABLES DE ENTRADA
!     MOP:  Operación de separación (0-Predicción de propiedades, 
!           1-Extracción líq-líq,2-Dest.Extract.)
!     JIST: Cantidad de compuestos a imprimir
!     IS:   Número de compuestos por página
!     NSOL: Máximo número de compuestos a ser escritos
!     MST:  Arreglo que contiene todos los compuestos
!     NGSL: Arreglo con la información de la cantidad de grupos de cada
!           compuesto
!     ISOL: Arreglo con la información de las posiciones de cada 
!           compuesto según su Coef. de Dist.
!     DIFBP:Diferencia de puntos de ebullición entre el compuesto generado
!           el componente a ser recuperado
!     VREL: RELATIVE VOLATILITY
!     Hvap:  VAP. HEAT
!     RDER: DENS. RAT.
!     TITULO:Caracter cadena que guarda texto informativo sobre el tipo 
!           de estructura a escribir
!     IDEV: Variable entera que informa si se llama a la subrutina para 
!           escribir en pantalla o en archivo.
!-----VARIABLES DE SALIDA
!     OPT:  Variable caracter donde se almacena la elección del usuario
!           cada vez que se imprime una página
!      
! ----------------------------------------------------------------------
      use PropertiesData
      PARAMETER (NMOD=3,NMSV=20,NGMAX=84,NSCM=10,NMG=150,
     *            NA=150)
      implicit real*8 (A-H,O-Z)
      DIMENSION MST(NMAX,NSCM,2),ISOL(NMAX,2)
      real*8 RDER(NMAX)
	integer identsolv(NMAX,8),nisomsolv(NMAX)
	character*35 nombresolv(NMAX,8),formulasolv(NMAX,8) 
      real*8 TB(NMAX,8),DENST(NMAX,8),VREL(NMAX),
     *VRELSOLV(NMAX,8)
      CHARACTER*70 TITULO
      CHARACTER*8 FS(NA)
      CHARACTER*1 SALTO,OPT
      real*8 VLIQ(NMAX,8),TCDB(NMAX,8),PCDB(NMAX,8),
     *VCDB(NMAX,8),ZCDB(NMAX,8),TFDB(NMAX,8),DENSDB(NMAX,8)
c	CHARACTER*1 SIGUE
      COMMON/PUNSUB/NPUNT(NA),NGRUP(NA),NCANT
      COMMON/EXTDIS/PSAT1,PSAT2,TAZEO,X1AZEO,ERROR
      COMMON/NOM/FS
	COMMON/SOLVISOM/NISOMSOLV,IDENTSOLV,NOMBRESOLV,FORMULASOLV,TB,
     *      DENST,VRELSOLV
	COMMON/gru/GRUPO1,NGRMAX
	COMMON/PAREQ/IPAREQ
      COMMON/DBPROP/TCDB,PCDB,VCDB,ZCDB,TFDB,VLIQ
C
!      SALTO = CHAR(12)
!      IF (JIST.EQ.IS) THEN  !define la cantidad máxima de páginas que 
!        NPAG = 1            !se imprimirán
!      ELSE
!        NPAG = JIST/IS + 1
!      END IF
!      L = 0
!      DO 705 K1=1,NPAG !en cada iteración se escribe una página
!        IF (IDEV.NE.6) THEN
!            WRITE (IDEV,150) SALTO 
!        END IF
!        IF (MOP.EQ.1) THEN ! ESCRIBE ENCABEZADO
!            IF (IDEV.EQ.6) THEN 
!                WRITE (IDEV,111)
!            ELSE
!                WRITE (IDEV,110)
!            ENDIF
!        ELSEIF(MOP.EQ.2)THEN
!            IF (TAZEO.NE.1) THEN
!                WRITE (IDEV,120)
!            ELSE
!                WRITE (IDEV,170)
!            END IF
!        ELSEIF(MOP.EQ.0)THEN
!            WRITE (IDEV,112)
!        END IF
!        IF ((MOP.EQ.1.AND.IDEV.EQ.6).OR.MOP.EQ.0) THEN ! **********....
!            WRITE(IDEV,693)
!        ELSE
!            WRITE(IDEV,692)
!        ENDIF
!        IF (TITULO.NE.'    ') THEN !Escribe tipo de estructuras a presentar
!            WRITE(IDEV,100) TITULO
!        END IF
!        IF (K1*IS.GT.JIST) THEN !indica la cantidad de compuestos que 
!            NSULT = JIST        !se imprimirán por página (sólo la 
!        ELSE                    !última página puede tener menos de IS 
!            NSULT = K1*IS       !compuestos)
!        END IF
!        DO 10 K=(K1-1)*IS+1,NSULT !en cada iteración se escribe un compuesto
!            L=L+1
!            !Escribe la posición y la fórmula química de la molécula
!            WRITE(IDEV,130)ISOL(L,2),(FS(NPUNT(MST(ISOL(K,1),I,1))),
!     *            MST(ISOL(K,1),I,2),I=1,NGSL(ISOL(K,1)))
!            IF (MOP.EQ.1) THEN ! Aquí se escriben las propiedades predichas
!                IF (IDEV.EQ.6) THEN
!                    WRITE(IDEV,141)DistCoef(ISOL(K,1)),Sel(ISOL(K,1)),
!     *              SolPow(ISOL(K,1)),SLost(ISOL(K,1)),
!     *              BP(ISOL(K,1)),DIFBP(ISOL(K,1)),MW(ISOL(K,1))
!!                    WRITE(IDEV,141)Sel(ISOL(K,1)),SolPow(ISOL(K,1))
!!     *             ,SLost(ISOL(K,1)),MW(ISOL(K,1)),BP(ISOL(K,1))
!!     *             ,DIFBP(ISOL(K,1)),VREL(ISOL(K,1))
!                ELSE
!                    WRITE(IDEV,140)Sel(ISOL(K,1)),SolPow(ISOL(K,1))
!     *             ,SLost(ISOL(K,1)),MW(ISOL(K,1)),BP(ISOL(K,1))
!     *             ,DIFBP(ISOL(K,1)),RDER(ISOL(K,1)),Hvap(ISOL(K,1)),
!     *		        DistCoef(ISOL(K,1)),Visc(ISOL(K,1)),VREL(ISOL(K,1))
!	          ENDIF
!	      ELSEIF(MOP.EQ.2)THEN
!                IF (TAZEO.NE.1) THEN
!	              WRITE(IDEV,180) Sel(ISOL(K,1)),SolPow(ISOL(K,1)),
!     +                MW(ISOL(K,1)),BP(ISOL(K,1)),Hvap(ISOL(K,1)),
!     +		        SLost(ISOL(K,1)),DistCoef(ISOL(K,1)),Visc(ISOL(K,1))
!                else 
!	              WRITE(IDEV,190) Sel(ISOL(K,1)),SolPow(ISOL(K,1)),
!     +                MW(ISOL(K,1)),BP(ISOL(K,1)),Hvap(ISOL(K,1)),
!     +		        DistCoef(ISOL(K,1)),Visc(ISOL(K,1))
!	          END IF
!	      ELSE
!	          WRITE(IDEV,141)MW(ISOL(K,1)),BP(ISOL(K,1)),
!     *              TC(ISOL(K,1)),PC(ISOL(K,1)),VC(ISOL(K,1)),
!     *             Visc(ISOL(K,1)),Dens(ISOL(K,1)),Kow(ISOL(K,1))
!            END IF
!			IF (NISOMSOLV(ISOL(K,1)).GT.0) THEN ! propiedades prop.pdb
!				DO 703 N=1,NISOMSOLV(ISOL(K,1))
!				    IF(MOP.EQ.0)THEN
!				        WRITE (IDEV,306)
!				        WRITE (IDEV,309)IDENTSOLV(ISOL(K,1),N),
!     +			        NOMBRESOLV(ISOL(K,1),N),FORMULASOLV(ISOL(K,1),N)
!                        WRITE(IDEV,310)
!                        WRITE(IDEV,311)TB(ISOL(K,1),N),TCDB(ISOL(K,1),N)
!     *                        ,PCDB(ISOL(K,1),N),VCDB(ISOL(K,1),N)*1000,
!     *                        (0.001*MW(ISOL(K,1))/(VLIQ(ISOL(K,1),n)))
! 310  FORMAT(7X,'B.P.[K]',5X,'  TC[K]',5X,'PC[atm]',1X,'VC[cm3/mol]',
!     *2X,'DENS[g/mL]')
! 311  FORMAT(2X,6F12.3)
!				    ELSE
!			            WRITE (IDEV,307)			!!!evaluacion de solventes	    
!				        WRITE (IDEV,308)IDENTSOLV(ISOL(K,1),N),
!     +			        NOMBRESOLV(ISOL(K,1),N),FORMULASOLV(ISOL(K,1),N),
!     +			        TB(ISOL(K,1),N)!,VRELSOLV(ISOL(K,1),N)
!                    ENDIF
! 703				CONTINUE
!			END IF
! 704	      CONTINUE
!            IF (K.EQ.NSOL) THEN
!		        IF (IDEV.EQ.6) THEN
! 220		            WRITE (IDEV,160)
!                    READ (5,150) OPT
!                    I1 = INDEX ('1234',OPT)
!                    IF (I1.EQ.0) THEN
!                        GO TO 220
!                    ELSE IF ((I1.EQ.2).OR.(I1.EQ.3)) THEN
!                        GO TO 1000
!                    ELSE IF (I1.EQ.4) THEN
!                        WRITE (6,5020)
!		                IF (MOP.EQ.1) THEN
!                            WRITE (IDEV,5030)
!                        ELSE
!                            IF (TAZEO.NE.1) THEN
!                                WRITE (IDEV,5040)
!                            ELSE
!                                WRITE (IDEV,5050)
!                            END IF
!                        END IF
!                        CALL PAUSA 
!                        OPT = '2'
!                        GO TO 1000
!                    END IF
!                END IF
!                GO TO 706
!            END IF
!  10    CONTINUE
!        IF (IDEV.EQ.6) THEN
!210         WRITE (IDEV,160) !escribe opciones
!            READ (5,150) OPT
!            I1 = INDEX ('12345',OPT)
!            IF (I1.EQ.0) THEN
!                GO TO 210
!            ELSE IF ((I1.EQ.2).OR.(I1.EQ.3).OR.I1.EQ.5) THEN
!                GO TO 1000
!            ELSE IF (I1.EQ.4) THEN
!                WRITE (6,5020)
!                IF (MOP.EQ.1) THEN
!                    WRITE (IDEV,5030)
!                ELSE
!                    IF (TAZEO.NE.1) THEN
!                       WRITE (IDEV,5040)
!                    ELSE
!                       WRITE (IDEV,5050)
!                    END IF
!                END IF
!                CALL PAUSA 
!                OPT = '2'
!                GO TO 1000
!            END IF
!        END IF
!705   CONTINUE
!706   CONTINUE    
!      call ci_ban1   
C
C------FORMATS
 306	FORMAT (2X,'PROP.PDB contains the following isomers for this ',
     +		'solvent molecular strucure:',/,
     +		2X,'R.N.',x,'Name',23x,'Formula')
 307	FORMAT (2X,'PROP.PDB contains the following isomers for this ',
     +		'solvent molecular strucure:',/,
     +		2X,'R.N.',x,'Name',23x,'Formula',20x,'B.P.(K)',3x)
 308	FORMAT (2X,I4,X,A27,A27,X,F5.1,4X,F8.2)	
 309	FORMAT (2X,I4,X,A27,A27,X,F5.1)
5020  format 
     *(' The  list  of  solvents  is shown according to the following',
     * ' nomenclature:',//,
     * '                UNIFAC solvent structure.     ',/,
     * '                           |',/,
     * '               ---------------------------',/,
     * '      13**     (CH2) 1(CH2COO) 1   (CH3) 2 ',/,
     * '      ----',/,
     * '       |')
5030  FORMAT (
     * '    Solvent ranking based on distribution coefficient.',///,
     * ' SELECT.   : solvent selectivity [ wt ].',/,
     * ' SOL.POW.  : solvent power [ wt % ].',/,
     * ' SOL.LOS.  : solvent loss [ wt % ].',/,
     * ' M.W.      : solvent molecular weight.',/,
     * ' B.P.      : solvent boiling point [ K ].',/,
     * ' B.P.DIFF. : solvent and solute boiling points difference',
     * ' [ K ].',/,
     * ' DENS.RAT. : solvent density ratio.',/,
     * ' VAP.HEAT  : solvent latent heat of vaporization [cal/molgr].'
     + ,/,
     * ' DIST.COEF.: solvent distribution coefficient [ wt ].',/,
     * ' VISC.SOLV.: viscosidad del solvente [mPa*seg].' )
5040  FORMAT (
     * '    Solvent ranking based on PI index.',///,
     * ' REL.VOL.  : relative volatility [ wt ].',/,
     * ' SOL.POW.  : solvent power [ wt % ].',/,
     * ' M.W.      : solvent molecular weight.',/,
     * ' B.P.      : solvent boiling point [ K ].',/,
     * ' HEAT VAP. : solvent latent heat of vaporization [cal/molgr].'
     + ,/,
     * ' MINX3.    : solvent concentration to break the azeotrope [mol',
     + '%].',/,
     * ' P.I.      : performance index.',/,
     * ' VISC.SOLV.: viscosidad del solvente [mPa*seg].' )
5050  FORMAT (
     * '    Solvent ranking based on PI index.',///,
     * ' REL.VOL.  : relative volatility [ wt ].',/,
     * ' SOL.POW.  : solvent power [ wt % ].',/,
     * ' M.W.      : solvent molecular weight.',/,
     * ' B.P.      : solvent boiling point [ K ].',/,
     * ' HEAT VAP. : solvent latent heat of vaporization [cal/molgr].'
     + ,/,
     * ' P.I.      : performance index.',/,
     * ' VISC.SOLV.: viscosidad del solvente [mPa*seg].' )
c
160    FORMAT (' ',/,3X,'CONTINUE: 1',5X,'BEGIN AGAIN: 2',5X,
     *         'EXIT: 3',5X,'HELP: 4',5X,'SEARCH: 5',5X,'> ',$)
150    FORMAT(A)
110    FORMAT(/,' SELECT.',2X,'SOL.POW.',2X,'SOL.LOS.',4X
     *         ,'M.W.',4X,'B.P.',3X,'B.P.DIFF.',2X,'DENS.RAT.',2X,
     *	       'VAP.HEAT',2X,'DIS.COEF.',2X,'VISC.SOLV.',2X,'REL.VOL.',
     *	       /,2X,'[WT]',4X,'[WT %]',
     +	       4X,'[WT %]',14X,'[K]',6X,'[K]',16X,
     *	       '[CAL/GMOL]',3X,'[WT]',
     *	       5X,'[mPa*seg]')
111    FORMAT(/,1X,' DIS.COEF.',1X,'   SELECT.',1X,'  SOL.POW.',1X
     *         ,' SOL.LOST.',1X,'      B.P.',1X,' B.P.DIFF.',1X,
     *	       '      M.W.',/,'    [ WT ]',1X,'    [ WT ]',1X,
     *	       '  [ WT % ]',1X,'  [ WT % ]',1X,'       [K]',1X,
     *         '       [K]')
112    FORMAT(/,1X,'  M.W.',1X,'     B.P.',1X,'       TC ',1X
     *         ,'      PC ',1X,'      VC ',1X,'    VISC',1X,
     *	       '     DENS',1X,'   LogPow',/,'         ',1X,'    [K]',1X,
     *	       '      [K]',1X,'     [atm]',1X,'[cm3/mol]',1X,
     *         ' [mPa*s]',1x,'   [g/mL]',' (298.15K)')
130    FORMAT(1X,/,1X,I4,'*',10(1A8,1I2))
692    FORMAT(110('*'))
693    FORMAT(80('*'))
100    FORMAT (A70)
140    FORMAT (F6.2,F10.2,E13.3,F8.1,F8.1,F9.1
     *	     ,F10.2,F13.2,F9.3,F11.4,F11.2)
141    FORMAT (F7.1,7(F10.2))
120	 FORMAT(' REL.VOL.',2X,'SOL.POW.',2X,' M.W. ',2X,'B.P. ',2X,
     +	      'VAP.HEAT',4X,' MINX3  ',2X,'  P.I. (REL.VOL./M.W.)',
     *	 'VISC. SOLV.',/,' ',
     *	 10X,'[ WT % ]',10X,'[K]',3X,'[CAL/GMOL]',3X,'[MOL %]',9X,
     *	 '(MINX3 M.W.))',4x,'[mPa*seg]')
170	 FORMAT(' REL.VOL.',2X,'SOL.POW.',2X,' M.W. ',2X,'B.P. ',2X,
     +	 'VAP.HEAT',11X,'  P.I. ',3x,'VISC.SOLV.',/,' ',10X,
     *	 '[ WT % ]',10X,'[K]',3X,'[CAL/GMOL]',5X,'(REL.VOL./M.W.)',
     *	 '[mPa*seg]')
180	 FORMAT(' ',F8.2,2X,F8.2,2X,F6.1,2X,F5.1,2X,F8.2,2X,
     *	      F8.2,2X,F8.4,X,F8.4)
190	 FORMAT(' ',F8.2,2X,F8.2,2X,F6.1,2X,F5.1,2X,F8.2,2X,
     *          F8.4,x,F8.4)
1000   RETURN
       END



