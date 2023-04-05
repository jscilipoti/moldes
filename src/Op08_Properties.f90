SUBROUTINE PROPERTIES1
      USE GRPS
      use PropertiesData
      use blockdatas
      use StructuresDesign
      use input_data, only:seleccion_grupos,checkdatos
     ! use PROP_PDB, only:checkdatos
      use design
      use SubGrupos
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (NA=150,INS1=10,NSVA=20,NSCM=10)
      real*8 LB
      CHARACTER*8 FS(NMG)
      character*35 nombresolv(NMAX,8),NOMBRE(8),formulasolv(NMAX,8),formula(8),ANAME,FORM2
	  CHARACTER*20 FORM1
      CHARACTER TITULO*70,OPT*1
      real*8 TB(NMAX,8),VLIQ(NMAX,8),TCDB(NMAX,8),PCDB(NMAX,8),VCDB(NMAX,8),ZCDB(NMAX,8),TFDB(NMAX,8)
      DIMENSION MAINGR(NMG),KGRUP(NSEL),ngr1(ins1),ngr2(ins1),ngr3(ins1),ngr4(ins1),MGR(NSEL),MST(NMAX,NSCM,2),&
      NGSSV(NSVA,5),NGSDV(NSVA,5),NPUNT(NA),NGRUP(NA),PROV2(NMAX)
      DIMENSION IGR0(NMAX,2),IGR1(NMAX,2),IGR2(NMAX,2),IGR3(NMAX,2) 
      INTEGER PROV(NMAX,2) !POSICIONES EN LA ESCRITURA
      integer comp(20),nisomsolv(NMAX),ident(8),UNIF(20),identsolv(NMAX,8)
      LOGICAL SALIR
      LOGICAL PROPER
      COMMON/PUNSUB/NPUNT,NGRUP,NCANT
      COMMON/PAREQ/IPAREQ
      COMMON/SOLVISOM/NISOMSOLV,IDENTSOLV,NOMBRESOLV,FORMULASOLV,TB,DENST,VRELSOLV
      COMMON/DBPROP/TCDB,PCDB,VCDB,ZCDB,TFDB,VLIQ
      COMMON/PROPIEDADES/PROPER
      COMMON/NOM/FS
      COMMON/US/MS(NCOM,NSCM,2),nms(NCOM)
      IDEV=6
      IDEVR=5
      PROPER=.TRUE.
      DO I=1,NMAX
        PROV(I,1)=I
        PROV(I,2)=I
      ENDDO
!---Apertura del banco de datos BADAUN
      call ab_ban
	!write(6,2900)
!---Seleccion de Parametros Unifac
 60   IPAREQ=2     
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
 802  FORMAT (/,'Choose lower and upper limits for the following properties:',/)
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
      CALL SELECCION_GRUPOS (IMPRIM,KGRUP,IFIN2,NGR1,NGR2,NGR3,NGR4,MGR,family)
      ms(1,1,1) = 1
	ms(1,1,2) = 1
	ms(1,2,1) = 2
	ms(1,2,2) = 7
	ms(1,3,1) = 15
	ms(1,3,2) = 1
	if(ipareq.eq.1) ms(1,3,1) = 14
	ms(2,1,1) = 17
	ms(2,1,2) = 1
      CALL CR_PUNT (MS,3,1,NGDV,MDV,NGSV,MSV,NGSV1,MSV1,IPAREQ,IFAM,IALAR,NALAR)
      CALL Store_Pr (IPAREQ)
      CALL STRUCTURE_GENERATOR (JIST,MS,3,NGSSV,NGSDV,SALIR,.TRUE.)
!-----Checkea existencia de estructuras en prop.pdb---------------------
     	write(6,*) ' Checking solvent structures with data base'


	DO 711 J=1,JIST
		call checkdatos (MST(j,:,:),nisom,ident,nombre,formula)
		NISOMSOLV(j)=NISOM
		IF (NISOM.GT.0) THEN
			DO 222 N=1,NISOM
				IDENTSOLV(J,N)=IDENT(N)
				NOMBRESOLV(J,N)=NOMBRE(N)
				FORMULASOLV(J,N)=FORMULA(N)
				
				READ(7,333,REC=IDENT(N)) NI,ANAME,FORM1,FORM2,&
                    (UNIF(K),K=1,20),TCDB(J,N),PCDB(J,N),VCDB(J,N),&
                     ZCDB(J,N),TFDB(J,N),TB(J,N),VLIQ(j,n)
                PCDB(J,N)=PCDB(J,N)/100000
 333 				FORMAT(I4,A35,A20,A35,20I3,7D13.6)
 222			CONTINUE
		END IF
 711	CONTINUE
      PRINT*,CHAR(7)
!----------------------------------
      !COMPUESTOS SIN GRUPOS FUNCIONALES
      TITULO = '**** SOLVENTS WITHOUT FUNCTIONAL GROUPS *****'
      CALL ESCRIBIR_RESULTADOS (0,INGR0,5,100,MST,IGR0,PROV2,TITULO,6,OPT)
      !COMPUESTOS CON UN GRUPO FUNCIONAL
      TITULO = '**** SOLVENTS WITH ONE FUNCTIONAL GROUP '
      CALL ESCRIBIR_RESULTADOS (0,INGR1,5,100,MST,IGR1,PROV2,TITULO,6,OPT)   
      !COMPUESTOS CON DOS GRUPOS FUNCIONALES
      TITULO = '**** SOLVENTS WITH TWO FUNCTIONAL GROUPS '
      CALL ESCRIBIR_RESULTADOS (0,INGR2,5,100,MST,IGR2,PROV2,TITULO,6,OPT)  
      !COMPUESTOS CON MÁS DE DOS GRUPOS FUNCIONALES
      TITULO = '**** SOLVENTS WITH MORE THAN TWO '//'FUNCTIONAL GROUPS ****'
      CALL ESCRIBIR_RESULTADOS (0,INGR3,5,100,MST,IGR3,PROV2,TITULO,6,OPT)
!-----FORMATS      
 20   format (' Choose the separation operation:',//,&
             29x,'Liquid-Liquid Extraction: 1',/,&
             29x,'Extractive  Distillation: 2')
 30   format (1x,/,60x,'> ',$)  
 590  format (/,' Give the family of solvents to be generated:',///,&
             10x,'-Aromatic solvents                  : 1',//&
             10x,'-Single substance groups            : 2',//&
             10x,'-No-aromatics with up to 12 groups in',/&
             10x,' the final structure                : 3',//&
             10x,'-Cyclic solvents                    : 4',///&
     	      51x,' > ',$)    
 800  format (/,' ',A36)
 801  format(E)
 803  FORMAT(A13,$)
 830  format (' ','Code',3x,'Parameter',30x,'min',7x,'max'/)
 840  format (' ',/,' If you want to change a parameter value give code'&
             ,' number.','  If not 0.',/,71x,'> ',$)
 860  format (' ',a63,7x)
      END


SUBROUTINE PROPERTIES
      use input_data,only:Leer_Comp
      use SubGrupos
      use CONSTANTES
   !   USE PROPI
      IMPLICIT real*8 (A-H,O-Z)
!INTERFACES
interface
    subroutine ab_ban1(mod)
    use input_data, only:Model_S
    integer, intent(out),optional::mod
    endsubroutine ab_ban1
endinterface
      PARAMETER (NMAX=60000,NG=10,NFILE=100)
      CHARACTER IFILE*35,FICHSEL(30)*35,FS(NMG)*8,nomcom*37
      integer::idr(nsel),nyr(nsel),comp(DiffStructGroups,2)
      DIMENSION isv(nsel),isd(nsel),idv(nsel),icyc(nsel),&
      iar(nsel),igrp(nsel),i3v(nsel),i4v(nsel)
      INTEGER UNIF(NMAX,20),RESP
      INTEGER GRUPOS(NMG),MAINGR(NMG),MOP
      real*8,DIMENSION(:),ALLOCATABLE::PMOL,TC ,PC,VC,BPOINT,VISC,DENSID,HVAP
      INTEGER,DIMENSION(:),ALLOCATABLE::CANTGRUP
      INTEGER,DIMENSION(:,:,:),ALLOCATABLE:: MSS
      INTEGER::isvt(NMG),isdt(NMG),&
      idvt(NMG),icyct(NMG),iart(NMG),igrpt(NMG),i3vt(NMG),i4vt(NMG),nsvt,&
      nsdt,ndvt,n3vt,n4vt,ncyct,nart,ngrpt
      COMMON /CONT/ PMG(NMG),DELTC(NMG),DELV(NMG),DELPC(NMG),&
     		    DELVI(NMG),DELN(NMG),DELTCA(NMG)
      COMMON/PUNSUB/NPUNT(NMG),NGRUP(NMG),NCANT
      LOGICAL LOG
      
      CALL AB_BAN1(model)
      IPROP=1
      IF(IPROP.EQ.1)THEN
       !INGRESO DE DATOS
  110   WRITE(6,600)
        READ(5,500,ERR=110) nopt
        IF(nopt.ne.1.and.nopt.ne.2)THEN
            GOTO 110
        ELSEIF(nopt.eq.2)THEN !Ingreso datos desde archivo
 115        WRITE(6,610)
            NFICH=0
            FICHSEL(:)=''
            DO
                NFICH=NFICH+1
 130            CONTINUE
                IFILE=''
                WRITE(6,620)NFICH
                READ (5,510,ERR=130)IFILE
                IF (IFILE=='') EXIT
                !IF (IFICHSEL.LT.1.OR.IFICHSEL.GT.NF) GOTO 130
                CALL BUSCARCADENA(FICHSEL,NFICH,IFILE,LOG)
                IF (LOG) THEN
                    WRITE (*,*) "This file has been chosen"
                    GOTO 130
                ENDIF
                FICHSEL(NFICH)=IFILE
            ENDDO
            NFICH=NFICH-1
           !Muestra ficheros seleccionados
            WRITE(6,630)
            DO I=1,NFICH
                WRITE(6,510) FICHSEL(I)
            ENDDO      
            WRITE(6,640)
            READ (5,500,ERR=112) RESP
            IF (RESP==1) GOTO 115
 112        CONTINUE
!           Averigua cantidad de componentes
            DO I=1,NFICH
                IFILE=FICHSEL(I)
                OPEN(unit=NFILE+I,file=IFILE,status='old',&
                access='direct', form='formatted',RECL=62,ERR=126)
                GOTO 127
               !En caso de error
 126            CONTINUE
                WRITE(6,650) FICHSEL(I)
                CALL PAUSA
                GOTO 1000
!
 127            CONTINUE
                DO J=1,NMAX
                   READ (NFILE+I,502,rec=j,err=123) IN
                ENDDO
 123            NCANTCOM=NCANTCOM+J-1
            ENDDO
            ALLOCATE(MSS(NCANTCOM,NG,2))
            MSS(:,:,:)=0
           !Lectura de datos
            K=1
            DO I=1, NFICH
                DO J=1,NMAX
                   READ (NFILE+I,503,rec=j,err=113)(MSS(K,L,1),L=1,10)
                   READ (NFILE+I,504,rec=j,err=113)(MSS(K,L,2),L=1,10)
                   K=K+1
                ENDDO
 113            CLOSE (UNIT=NFILE+I)
            ENDDO
        ELSEIF(nopt.eq.1)THEN
           !Apertura del banco de datos
            call ab_ban1(model)						 
!
           !Seleccion de Parametros Unifac
            !CALL TABLA_PARAMETROS(IPAREQ)
            ipareq=2
            !Lectura de los datos segun la seleccion de parametros
!            write (idev,240) !*** Loading data ***
            ngrmax = Size_LSubGroups()
            DO i=1,ngrmax
                GRUPOS(i) = i
                MAINGR(i) = mainsg (i,ipareq)
            ENDDO
            call car_car (ipareq,isvt,isdt,idvt,i3vt,i4vt,icyct,iart,igrpt,&
                                    nsvt,nsdt,ndvt,n3vt,n4vt,ncyct,nart,ngrpt)
!
            NCANTCOM=1
            NOMCOM="COMPOUND"
            write (idev,660)
            CALL leer_comp (ipareq,nomcom,comp)
!                           (ipareq,fs,grupo1,ngrmax,nomcom,idr,nyr,msol)
            ALLOCATE(MSS(NCANTCOM,NG,2))
            MSS(:,:,:)=0
      		DO i=1,10
			    MSS(1,I,1)=idr(I)
			    MSS(1,I,2)=nyr(I) 
            ENDDO
        ENDIF
      ELSE
        call pausa
      ENDIF
      ALLOCATE(CANTGRUP(NCANTCOM))
      CANTGRUP(:)=0
      !Acomodación de grupos
      CALL CR_PUNT3(NCANTCOM,MSS,NCANT,NPUNT,NGRUP,CANTGRUP)
      !Propiedades de grupos
      DO 10 I=1,NCANT
        CALL Leer_Pr(NGRUP(I),2,FS(I),MM,MJ,MK,MI,MH,PMG(I),&
       RPAR,QPAR,DELTC(I),DELPC(I),DELV(I),delvi(I),deln(I),DELTCA(I))
  10  CONTINUE
      !Imprime compuestos seleccionados
      WRITE(6,603) 
      DO I=1, NCANTCOM
        WRITE(6,601) (FS(NPUNT(MSS(I,J,1))),MSS(I,J,2),J=1,CANTGRUP(I))
      ENDDO
      WRITE(6,618) NCANTCOM
      CALL PAUSA
!-----CALCULO DE PROPIEDADES DE COMPUESTO PURO
!     Inicialización de vectores donde se almacenarán las propiedades
      N=NCANTCOM
      ALLOCATE (PMOL(N),TC(N),PC(N),VC(N),BPOINT(N),VISC(N),DENSID(N),HVAP(N))
      PMOL(:)=0.
      TC(:)=0.
      PC(:)=0.
      VC(:)=0.
      BPOINT(:)=0.
      VISC(:)=0.
      DENSID(:)=0.
      HVAP(:)=0.      
      DO I=1,NCANTCOM
        !CALL PROPI2 (MSS(I,:,:),1,ICOMP,PMOL(I),TC(I))
      ENDDO
      CALL PAUSA
!-----FORMATS
 500  FORMAT(I1)
 501 	FORMAT (20I3)
 503  FORMAT (10(I3,3X))
 504  FORMAT (10(3X,I3))
 502  FORMAT(I3)
 510  FORMAT(a)
 600  FORMAT(////' ','Component read from:',&
            //,11x,' keyboard:           :  1',&
            /,11x,' input file          :  2',&
     	     //,60x,'> ',$)
 601  FORMAT (/,"- ",1X,10(A9,I2)) 
 603  FORMAT (///, "These compounds have been entered: ",/) 
 610  FORMAT(//,"Enter the name(s) of the file(s) to open (<ret> to finish)",/)
 618  FORMAT (//,1X,"*",I3," COMPOUNDS *")
 620  FORMAT(/,"File N°",I3,":   "$)
 630  FORMAT (//,"These files have been selected: ",/)
 640  FORMAT (//,"If it is OK: <ret>. If not: 1 ",10X,">",$)
 650  FORMAT (//," *ERROR* al abrir archivo",A35," - To exit: 1" ,/,70X,'<RET>',$)
 660  format (' Give the group composition of the following ',&
             'components: id1,ny1,id2,ny2,etc.',/,' Where: ',/,&
             10x,'id1,id2: subgroup identification number',/,&
             10x,'ny1,ny2: number of subgroups id1,id2,etc',/)
1000  CONTINUE
      END
!

!
!

!
!
!

