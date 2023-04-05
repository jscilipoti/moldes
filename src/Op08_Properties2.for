      SUBROUTINE PROPERTIES2
      USE GRPS
      use PropertiesData
      use input_data, only:seleccion_grupos
      use SubGrupos
      use CONSTANTES
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (INS1=10,NSCM=10,NG=70)!,NA=150,NSVA=20,
!     *,NCOM=3,NMODEL=3)
!      real*8 LB
      CHARACTER*8 FS(NMG)
      CHARACTER*21 family(8)
      logical::SALIR
!      character*17 tabla(nmodel)
!      character*63 cartel(8),concar(2),carte2(6)
!      character*64 icartl(3),icart2(3)
!      CHARACTER*36 cartel3(7)
!      character*35 nombresolv(NMAX,8),NOMBRE(8),formulasolv(NMAX,8),
!     *formula(8),ANAME,FORM2
!	CHARACTER*20 FORM1
!      CHARACTER TITULO*70,OPT*1
!      real*8 TB(NMAX,8),VLIQ(NMAX,8),TCDB(NMAX,8),PCDB(NMAX,8),
!     *VCDB(NMAX,8),ZCDB(NMAX,8),TFDB(NMAX,8)
      DIMENSION MAINGR(NMG),MST(NMAX,NSCM,2),KGRUP(NSEL),ngr1(ins1),
     *ngr2(ins1),ngr3(ins1),ngr4(ins1),MGR(NSEL)
!-----COMMONS
      COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NUM
      COMMON/PUNGRU/NPINT(NG),NINTT(NG),NUMINT
!     *,
!     *NGSSV(NSVA,5),NGSDV(NSVA,5),dato(8),dato2(5),idato(3),NPUNT(NA),
!     *NGRUP(NA),NGSL(NMAX),PROV2(NMAX),dato3(14)
!      DIMENSION IGR0(NMAX,2),IGR1(NMAX,2),IGR2(NMAX,2),IGR3(NMAX,2) 
!      INTEGER PROV(NMAX,2) !POSICIONES EN LA ESCRITURA
!      integer comp(20),nisomsolv(NMAX),ident(8),UNIF(20),
!     *identsolv(NMAX,8)
!      LOGICAL SALIR
      LOGICAL PROPER
!      COMMON/PUNSUB/NPUNT,NGRUP,NCANT
      COMMON/PAREQ/IPAREQ
!      common /cant/ dato,dato2,dato3,idato,family,tabla,cartel,icartl,
!     *       concar,carte2,icart2,cartel3
!      COMMON/SOLVISOM/NISOMSOLV,IDENTSOLV,NOMBRESOLV,FORMULASOLV,TB,
!     *DENS,VRELSOLV
!      COMMON/DBPROP/TCDB,PCDB,VCDB,ZCDB,TFDB,VLIQ
      COMMON/PROPIEDADES/PROPER
!      COMMON/NOM/FS
!      COMMON/US/MS(NCOM,NSCM,2),nms(NCOM)
      IDEV=6
      IDEVR=5
      PROPER=.TRUE.
!      DO I=1,NMAX
!        PROV(I,1)=I
!        PROV(I,2)=I
!      ENDDO

	!write(6,2900)
c---Seleccion de Parametros Unifac
! 60   IPAREQ=2     
 !60   CALL TABLA_PARAMETROS(IPAREQ)	
!---Give the family of solvents to be generated:
1200  write (idev,590) 
      read (idevr,*,err=1200) ifam
      if ((ifam.lt.1).or.(ifam.gt.4)) then
		go to 1200
      end if
!---Selección cotas
      icod=0
      WRITE (*,802)
 802  FORMAT (/,'Choose lower and upper limits for the following 
     *properties:',/)
 !	write (idev,830)
!	do 710 i=1,7 !Opciones del 1 al 8
!	  LB=0.
!	  UB=0.
!		write (idev,800) cartel3(i) !,dato3(i*2-1),dato3(i*2)
!		WRITE (IDEV,803)'Lower limit: '
!  		READ(IDEVR,*)LB
!  		WRITE (IDEV,803)'Upper limit: '
!  		READ(IDEVR,*)UB
!  		IF(LB.NE.0.)DATO3(I*2-1)=LB
!	  IF(UB.NE.0.)DATO3(I*2)=UB
! 710	continue
! 850	write (idev,840)
!	read (idevr,*,err=850) icod
!3090	if (icod.ne.0) then
!		write (idev,800) cartel3(icod)
!		WRITE (IDEV,803)'Lower limit: '
!		READ(IDEVR,801)LB
!		WRITE (IDEV,803)'Upper limit: '
!		READ(IDEVR,801)UB
!		IF(LB.NE.0.)DATO3(I*2-1)=LB
!	  IF(UB.NE.0.)DATO3(I*2)=UB
!		go to 850
!	end if
!---Cotas
!      PMMALB=DATO3(1)
!      PMMA=DATO3(2)
!      BPLB=DATO3(3)
!      BPUB=DATO3(4)
!      TcLB=DATO3(5)
!      TcUB=DATO3(6)
!      PcLB=DATO3(7)
!      PcUB=DATO3(8)
!      VcLB=DATO3(9)
!      VcUB=DATO3(10)
!      VisLB=DATO3(11)
!      VisUB=DATO3(12)
!      DensLB=DATO3(13)
!      DensUB=DATO3(14)

      PMMALB=0
      PMMA=9000
      BPLB=0
      BPUB=9000
      TcLB=0
      TcUB=9000
      PcLB=0
      PcUB=9000
      VcLB=0
      VcUB=9000
      VisLB=0
      VisUB=9000
      DensLB=0
      DensUB=9000
      CALL SELECCION_GRUPOS (IMPRIM,KGRUP,IFIN2,NGR1,
     *                          NGR2,NGR3,NGR4,MGR,family)
      I=1
      DO WHILE (NGDV(I).NE.0)           !Carga de grupos intermedios ingresados por el usuario
        CALL CR_PUNTF(NGDV(I),NPUNT,NGRUP,NPINT,NINTT)
        I=I+1
      ENDDO
      I=1
      DO WHILE (NGSV(I).NE.0)           !Carga de grupos terminales ingresados por el usuario
        CALL CR_PUNTF(NGSV(I),NPUNT,NGRUP,NPINT,NINTT)
        I=I+1
      ENDDO
      I=1
      DO WHILE (NGSV1(I).NE.0)
        CALL CR_PUNTF(NGSV1(I),NPUNT,NGRUP,NPINT,NINTT)
        I=I+1
      ENDDO
!      ms(1,1,1) = 1
!	ms(1,1,2) = 1
!	ms(1,2,1) = 2
!	ms(1,2,2) = 7
!	ms(1,3,1) = 15
!	ms(1,3,2) = 1
!	if(ipareq.eq.1) ms(1,3,1) = 14
!	ms(2,1,1) = 17
!	ms(2,1,2) = 1
!      CALL CR_PUNT (MS,3,1,NGDV,MDV,NGSV,MSV,NGSV1,MSV1,
!     *                  IPAREQ,IFAM,IALAR,NALAR)
      CALL Store_Pr (IPAREQ)                     !Lee parámetros grupales
    !  CALL STRUCTURE_GENERATOR (JIST,MST,MS,3,NGSSV,
      !*                  TEMP1,NGSDV,SALIR,.TRUE.)
!-----Checkea existencia de estructuras en prop.pdb---------------------
!     	write(6,*) ' Checking solvent structures with data base'
!
!
!	DO 711 J=1,JIST
!	  do 9 i=1,19,2
!		    k =	(i+1)/2
!			comp(i)=MST(J,k,1)
!			comp(i+1)=MST(J,k,2)
!   9	  continue
!		call checkdatos (comp,nisom,ident,nombre,formula)
!		NISOMSOLV(j)=NISOM
!		IF (NISOM.GT.0) THEN
!			DO 222 N=1,NISOM
!				IDENTSOLV(J,N)=IDENT(N)
!				NOMBRESOLV(J,N)=NOMBRE(N)
!				FORMULASOLV(J,N)=FORMULA(N)
!				imp=7
!				READ(IMP,333,REC=IDENT(N)) NI,ANAME,FORM1,FORM2,
!     *               (UNIF(K),K=1,20),TCDB(J,N),PCDB(J,N),VCDB(J,N),
!     *                ZCDB(J,N),TFDB(J,N),TB(J,N),VLIQ(j,n)
!                PCDB(J,N)=PCDB(J,N)/100000
! 333 				FORMAT(I4,A35,A20,A35,20I3,7D13.6)
! 222			CONTINUE
!		END IF
! 711	CONTINUE
!      PRINT*,CHAR(7)
!	!-----------------------------------------------------------------
!!-----Ordena por grupos funcionales
!      CALL OR_POR_GRUPOS_FUNCIONALES (IPAREQ,MST,NGSL,JIST,
!     *                                         PROV,
!     *                                         INGR0,INGR1,INGR2,INGR3,
!     *                                         IGR0,IGR1,IGR2,IGR3)
!!----------------------------------
!      !COMPUESTOS SIN GRUPOS FUNCIONALES
!      TITULO = '**** SOLVENTS WITHOUT FUNCTIONAL GROUPS *****'
!      CALL ESCRIBIR_RESULTADOS (0,INGR0,5,100,MST,NGSL,IGR0,
!     *					             PROV2,PROV2,TITULO,6,OPT)
!      !COMPUESTOS CON UN GRUPO FUNCIONAL
!      TITULO = '**** SOLVENTS WITH ONE FUNCTIONAL GROUP '
!      CALL ESCRIBIR_RESULTADOS (0,INGR1,5,100,MST,NGSL,IGR1,
!     *					             PROV2,PROV2,TITULO,6,OPT)   
!      !COMPUESTOS CON DOS GRUPOS FUNCIONALES
!      TITULO = '**** SOLVENTS WITH TWO FUNCTIONAL GROUPS '
!      CALL ESCRIBIR_RESULTADOS (0,INGR2,5,100,MST,NGSL,IGR2,
!     *					             PROV2,PROV2,TITULO,6,OPT)  
!      !COMPUESTOS CON MÁS DE DOS GRUPOS FUNCIONALES
!      TITULO = '**** SOLVENTS WITH MORE THAN TWO '//
!     *                     'FUNCTIONAL GROUPS ****'
!      CALL ESCRIBIR_RESULTADOS (0,INGR3,5,100,MST,NGSL,IGR3,
!     *					             PROV2,PROV2,TITULO,6,OPT)
!-----FORMATS      
 20   format (' Choose the separation operation:',//,
     *        29x,'Liquid-Liquid Extraction: 1',/,
     *        29x,'Extractive  Distillation: 2')
 30   format (1x,/,60x,'> ',$)  
 590  format (/,' Give the family of solvents to be generated:',///,
     *        10x,'-Aromatic solvents                  : 1',//
     *        10x,'-Single substance groups            : 2',//
     *        10x,'-No-aromatics with up to 12 groups in',/
     *        10x,' the final structure                : 3',//
     *        10x,'-Cyclic solvents                    : 4',///
     *	      51x,' > ',$)    
 800  format (/,' ',A36)
 801  format(E)
 803  FORMAT(A13,$)
 830  format (' ','Code',3x,'Parameter',30x,'min',7x,'max'/)
 840  format (' ',/,' If you want to change a parameter value give code'
     *        ,' number.','  If not 0.',/,71x,'> ',$)
 860  format (' ',a63,7x)
      END 