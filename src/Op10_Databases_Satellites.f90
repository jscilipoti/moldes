SUBROUTINE VER_COM(NR)
    use SubGrupos
    IMPLICIT real*8 (A-H,O-Z)
	INTEGER NR,NI,UNIF(20),K1,NO,OB,PARINT,BD,grupo1(150),maingr(150)
	INTEGER IPAREQN(3),IPAREQLN(4)
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1,IPAREQT(3)
	CHARACTER*4  TIPO,CAR*3
	character*8 fsi(150)
	character fs*8,MG*8,caract*3
	logical dat, IPAREQL
	integer::opt
	COMMON/OPC/NO,OB,PARINT,BD
	COMMON/GRUPOSRAM/nlgr,nl,MG,fs,item,caract,mm,&
     						  mj,mk,mi,mh,rpar,qpar,&
     						  pmg,deltc,delpc,delv,delvi,deln	
	COMMON/NOP1/K1
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO
      data IPAREQT /'liquid-liquid       ',&
                   'liquid-vapor        ',&
                   'infinite dilution   '/   

      IF (OB.EQ.1) THEN
  30	  IF(NR.EQ.0) THEN
	      dat=.true.
  	      WRITE(6,*) ' '
	      WRITE(6,*) 'NUMERO IDENTIFICATORIO DEL COMPONENTE',' QUE SE DESEA VER = '
	      READ(5,*)  NR
	  else
	      dat=.false.
	  END IF

	  CALL LEC_COM(NR)

	  IF (NI.EQ.0) THEN
	      WRITE(6,*)'   '
  	      WRITE(6,120) 'NO EXISTE NINGUN COMPONENTE CON NUMERO '
            WRITE(6,122) 'IDENTIFICATORIO = ',NR,' EN "PROP.PDB"'
 120	      FORMAT(10X,A39)
 122	      FORMAT(10X,A18,I4,A14)
	      nr=1
	      GOTO 30
	  END IF

	  CALL TITULO
	  CALL PAGINA1

	  CALL PAUSA

	  CALL TITULO
	  CALL PAGINA2

	  CALL PAUSA

	  CALL TITULO
	  CALL PAGINA3	

	  CALL PAUSA

	  if(dat) then
      	    WRITE(6,*) ' '
	      WRITE(6,*) 'INDICAR :'
	      WRITE(6,*) '    1- VER OTRO COMPONENTE'
	      WRITE(6,*) '    2- ACCEDER A OTRA OPCION'
	      WRITE(6,*) '    3- TERMINAR'
	      WRITE(6,*)'    '
	      READ(5,*)  K1
	      IF (K1.EQ.1) then
	          nr=1
	          GOTO 30
	      end if
	  end if
	ELSE IF (OB.EQ.2) THEN
        CALL TABLA_PARAMETROS(IPAREQ)
        CALL MOSTRAR_TABLA (IPAREQ,OPT)
	ELSE IF (OB.EQ.3) THEN
 303	  CONTINUE
        IPAREQL=.FALSE.
        IPAREQ=4
        ngrmax = Size_LSubGroups()
        do 21 i=1,ngrmax
            grupo1(i) = i
            maingr(i) = mainsg (i,ipareq)
            call carac (i,ipareq,fsi(i),car)
  21    continue
        write(6,*) 'Choose group:                    '
        write(6,*)' '
        call Groups_Present ()
        WRITE(6,6)
        read(5,*) nlgr
!c
!c-----Lectura en gruposram de la informaci�n del grupo seleccionado
!c
 124		irec2=nlgr+(ipareq-1)*150
		read (14,54,rec=irec2)nl,MG,fs,item,item,caract,mm,&
     						  mj,mk,mi,mh,rpar,qpar,&
     						  pmg,deltc,delpc,delv,delvi,deln
        IPAREQN(:) = 0
        IPAREQLN(1:3)= 0
        IPAREQLN(4)=IREC2
        DO I=1,3
            CALL BUSCAR (I,IPAREQL,FS,NUM)
            IF (IPAREQL) THEN
                IPAREQN(I)=I
                IPAREQLN(I)=NUM
            ENDIF
        ENDDO
        WRITE(6,16)
        DO I=1,3
            IF (IPAREQN(I).NE.0) THEN
                WRITE (6,17) I,IPAREQT(IPAREQN(I))
            ENDIF
        ENDDO  
        CALL PAUSA
!c-----Presentaci�n por pantalla de la informaci�n del grupo y
!c-----selecci�n del par�metro a modificar

	  WRITE(6,11)

	  CALL NUM_PRO
	  WRITE(6,*) ' '
	  WRITE(6,*) 'Choose :'
	  WRITE(6,*) '    1- See another component'
	  WRITE(6,*) '    2- Access another choose'
	  WRITE(6,*) '    3- Exit'
	  WRITE(6,6)
	  READ(5,*)  K1
	  IF (K1.EQ.1) GOTO 303
	ENDIF
!C-----FORMATOS
  6   FORMAT (/,52x,'> ',$)
  11  FORMAT (/,' The current properties are as follows :',/)
  16  FORMAT (/,' The select goup exists in the',' following tables of parameters:',/)
  17  FORMAT (I2,' - ',A20)     
  54	format (i4,2a8,2i4,a3,5i2,8d15.8)
	RETURN
 	END

SUBROUTINE MOD_COM

      use CONSTANTES
      use SubGrupos
      implicit real*8 (A-H,O-Z)
      parameter (NG=70)
	INTEGER NI,UNIF(20),K3,NP,K10,NO,OB,PARINT,BD,grupo1(150),maingr(150)
	INTEGER IPAREQN(3),IPAREQLN(4),OPT,GRUPOS(NMG)
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
      character*37 nomco1
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1,IPAREQT(3)
	CHARACTER  TIPO*4,car*3
	character*8 fsi(150)
	dimension idrf(nsel),nyrf(nsel),AINTV(NG)
	COMMON/OPC/NO,OB,PARINT,BD
	COMMON/NOP2/K10
	character fs*8,MG*8,caract*3,RESP*1,RESP1
	CHARACTER*8 MGV(NMG)
	LOGICAL IPAREQL
	COMMON/GRUPOSRAM/nlgr,nl,MG,fs,item,caract,mm,&
     						  mj,mk,mi,mh,rpar,qpar,&
     						  pmg,deltc,delpc,delv,delvi,deln
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO
      data IPAREQT /'liquid-liquid       ','liquid-vapor        ','infinite dilution   '/     

      idev=6
      idevr=5
      IF (OB.EQ.1) THEN
 300    WRITE(6,10)
	  READ(5,*) NR

	  CALL LEC_COM(NR)

	  IF(NI.NE.0) GOTO 400
        WRITE (6,15)NR
 !200	  FORMAT(10X,A39)
 !222	  FORMAT(10X,A18,I4,A14)
	  GOTO 300

	  NI=NR

 400    WRITE(6,11)

 	  CALL NUM_PRO

 500	  WRITE(6,12)
	  READ(5,*) K3
	  DO 101 I=1,K3
	      WRITE(6,*)' '
	      WRITE(6,*)'PROPIEDAD A CORREGIR O AGREGAR ( 1-50 ) ='
	      READ(5,*) NP

	      CALL IN_PRO(NP)

  101	  CONTINUE
	  WRITE(6,13)

	  CALL NUM_PRO

        WRITE(6,14)
	  READ(5,*) K10
	  IF(K10.EQ.1) GOTO 500

	  WRITE(7,100,REC=NI) NI,ANAME,FORM1,FORM2,(UNIF(K),K=1,20),&
                           TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,TIPO,&
                           (ANT(KI),KI=1,5),(CPG(KJ),KJ=1,5),&
                           (CPL(KL),KL=1,5),CIANT,CSANT,CICPG,CSCPG,&
                           CICPL,CSCPL,H0G,DHVAP
 100 	  FORMAT(I4,A35,A20,A35,20I3,12D13.6,A4,23D13.6)
	  IF(K10.EQ.2) GOTO 300
	  GOTO 600
!C     
      ELSE IF (OB.EQ.2) THEN 
 270    CALL TABLA_PARAMETROS(IPAREQ)
 280    WRITE(6,*)'Would you like to see parameters? (y/n)'
        READ (5,*) RESP
        IF (RESP.NE.'y'.AND.RESP.NE.'n') GOTO 280
        IF (RESP.EQ.'y') THEN
            CALL MOSTRAR_TABLA (IPAREQ,OPT)
        ENDIF
        IF (OPT.EQ.2) GOTO 600
        ngrmax = max_sub_int (ipareq)
        DO I=1,NGRMAX
            irec =I+(ipareq-1)*70
			grupos(i) = i
			READ(BD,24,rec=irec)MGV(I)
        ENDDO
        CALL Groups_Present ()
        WRITE (6,23)
        WRITE (6,20) 
        READ  (5,*) M
        WRITE (6,21) 
        READ  (5,*) N
        CALL MODIFICAR_INTRCN(M,N,IPAREQ)
        WRITE (6, 25)
  281   READ (5,18) RESP1
        IF (RESP1.NE.'y'.AND.RESP1.NE.'n') GOTO 281
        IF (RESP1.EQ.'y') GOTO 270
      ELSE IF (OB.EQ.3) THEN 
!C
  303	  CONTINUE
        IPAREQL=.FALSE.
        IPAREQ=4
        ngrmax = Size_LSubGroups()
        do 31 i=1,ngrmax
            grupo1(i) = i
            maingr(i) = mainsg (i,ipareq)
            call carac (i,ipareq,fsi(i),car)
  31    continue
        write(6,*) 'Choose group to modify:                    '
        write(6,*)' '
        call Groups_Present ()
        WRITE(6,6)
        read(5,*) nlgr
!c
!c-----Lectura en gruposram de la informaci�n del grupo seleccionado
!c
 124		irec2=nlgr+(ipareq-1)*150
		read (14,54,rec=irec2)nl,MG,fs,item,item,caract,mm,&
     						  mj,mk,mi,mh,rpar,qpar,&
     						  pmg,deltc,delpc,delv,delvi,deln
        IPAREQN(:) = 0
        IPAREQLN(1:3)= 0
        IPAREQLN(4)=IREC2
        DO I=1,3
            CALL BUSCAR (I,IPAREQL,FS,NUM)
            IF (IPAREQL) THEN
                IPAREQN(I)=I
                IPAREQLN(I)=NUM
            ENDIF
        ENDDO
        WRITE(6,16)
        DO I=1,3
            IF (IPAREQN(I).NE.0) THEN
                WRITE (6,17) I,IPAREQT(IPAREQN(I))
            ENDIF
        ENDDO  
        CALL PAUSA  
       
      ! 54	format (i4,2a8,2i4,a3,5i2,8d15.8)
!c
!c-----Presentaci�n por pantalla de la informaci�n del grupo y
!c-----selecci�n del par�metro a modificar
!c
 125	  WRITE(6,11)
!C
	  CALL NUM_PRO
!C 
 503	  WRITE(6,12)
	  READ(5,*) K3 
		DO 33 I=1,K3
 504	      WRITE(6,*)' '
	      WRITE(6,*)'Property to modify or add ( 1-14 ) ='
	      READ(5,*) NP
	      IF (NP.LT.1.OR.NP.GT.14) GOTO 504
!C
	      CALL IN_PRO(NP)
!C
  33	  CONTINUE
  	  WRITE(6,13)
!C
	  CALL NUM_PRO
!C
        WRITE(6,19)
  505   READ(5,18)RESP
        IF (RESP.NE.'y'.AND.RESP.NE.'n') GOTO 505
        IF (RESP.EQ.'n') GOTO 125
	  WRITE(6,14)
	  READ(5,*) K10
	  IF(K10.EQ.1) GOTO 503
	  
!C
        DO I=1,4
            IF (IPAREQLN(I).NE.0) THEN
  	          WRITE (14,9,rec=IPAREQLN(I))nl,MG,fs,item,item,&
			              caract,mm,mj,mk,mi,mh,rpar,qpar,&
						  pmg,deltc,delpc,delv,delvi,deln
            ENDIF
        ENDDO
!C
	  IF(K10.EQ.2) GOTO 303
	  GOTO 600
!C
  	ENDIF
!C-----FORMATOS
  10  FORMAT (/,' NUMERO IDENTIFICATORIO DEL COMPONENTE ','QUE DESEA MODIFICARSE =')
  11  FORMAT (/,' The current properties are as follows :',/)
  12  FORMAT (/,' Total number of data to modify and/or add =')
  13  FORMAT (/,' The corrected data are as follows:',/)
  14  FORMAT (/,' Choose :',/,&
               '    1- Modify more',/,&
               '    2- Modify another component',/,&
               '    3- Access another choose',/,&
               '    4- Exit',/) 
  15  FORMAT (/,'NO EXISTE NINGUN COMPONENTE CON NUMERO ','IDENTIFICATORIO ',I4,' EN "PROP.PDB"')
  16  FORMAT (/,' The select goup exists in the', ' following tables of parameters:',/)
  17  FORMAT (I2,' - ',A20)
  18  FORMAT (a1)
  19  FORMAT (/,'Is that correct? (y/n):',//,52x,'> ',$)
  20  FORMAT (/,' m = ',$)
  21  FORMAT (/,' n = ',$)
  22  FORMAT (a8,70d12.5)
  23  FORMAT (//,'SELECT PARAMETER a(mn)')
  24  FORMAT (a8)
  25  FORMAT (/,2X,'Would you like to modify some other parameter? (y/n)',//,52x,'> ',$)
  6   FORMAT (/,52x,'> ',$)
  54	format (i4,2a8,2i4,a3,5i2,8d15.8)
   9	format (i4,2a8,2i4,a3,5i2,8d15.8)
 600	RETURN
	END 
 
 	SUBROUTINE TITULO

	INTEGER NI,UNIF(20),I,IU,IU1,U2,U21,NL1,NL2
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	CHARACTER*1  NAME(35),FORM(20)
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO

	DO 10 I=1,35
	   NAME(I)=ANAME(I:I)
  10	CONTINUE
	DO 20 I=1,20
	   FORM(I)=FORM1(I:I)
  20	CONTINUE

	NL1=0
	DO 30 I=35,1,-1
	   IF (NAME(I).EQ.' ') THEN
	       NL1=NL1+1
	   ELSE 
	      GOTO 35
	   END IF
  30	CONTINUE
  35	NL2=0
	DO 40 I=20,1,-1
	   IF (FORM(I).EQ.' ') THEN
	       NL2=NL2+1
	   END IF
  40	CONTINUE

	IU1=35-NL1
	IU=IU1+1
	U21=20-NL2
	U2=U21+1

	WRITE(6,1000)
	WRITE(6,1100)
	WRITE(6,1200) 'Numero identificatorio : ',NI,(FORM(I),I=U2,20),&
                      (FORM(I),I=1,U21)
	WRITE(6,1300) (NAME(I),I=IU,35),(NAME(I),I=1,IU1)
	WRITE(6,1400) 'Formula estructural : ',FORM2
	WRITE(6,1000)
	call pausa

1000	FORMAT(80('#'))
1100	FORMAT(' ','#',77X,'#')
1200	FORMAT(' ','#',3X,A25,I4,23X,20A1,2X,'#')
1300	FORMAT(' ','#',40X,35A1,2X,'#')
1400	FORMAT(' ','#',3X,A22,A35,17X,'#')

	RETURN
	END

	SUBROUTINE PAGINA1

	INTEGER NI,UNIF(20)
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO

	WRITE(6,1000)
	WRITE(6,1100) 'PROPIEDAD','UNIDADES','VALOR'
	WRITE(6,1000)
	WRITE(6,1200)
	WRITE(6,1300) 'Peso molecular','kg/kmol',PM
	WRITE(6,1400) 'Temperatura critica','K',TC
	WRITE(6,1500) 'Presion critica','Pa',PC
	WRITE(6,1600) 'Volumen critico','m**3/kmol',VC
	WRITE(6,1700) 'Factor de compresibilidad critico',ZC
	WRITE(6,1200)
	WRITE(6,1800) 'Punto de fusion','K',TF
	WRITE(6,1900) 'Punto normal de ebullicion','K',TEB
	WRITE(6,2000) 'Volumen molar liquido','m**3/kmol',VLIQ
	WRITE(6,1000)
	WRITE(6,1200)

1000	FORMAT(' ','|',37X,'|',16X,'|',22X,'|')
1100	FORMAT(' ','|',14X,A9,14X,'|',4X,A8,4X,'|',9X,A5,8X,'|')
1200	FORMAT(80('-'))
1300	FORMAT(' ','|',2X,A14,21X,'|',4X,A7,5X,'|',5X,G13.6,4X,'|')
1400	FORMAT(' ','|',2X,A19,16X,'|',8X,A1,7X,'|',5X,G13.6,4X,'|')
1500	FORMAT(' ','|',2X,A15,20X,'|',7X,A2,7X,'|',5X,G13.6,4X,'|')
1600	FORMAT(' ','|',2X,A15,20X,'|',4X,A9,3X,'|',5X,G13.6,4X,'|')
1700	FORMAT(' ','|',2X,A33,2X,'|',16X,'|',5X,G13.6,4X,'|')
1800	FORMAT(' ','|',2X,A15,20X,'|',8X,A1,7X,'|',5X,G13.6,4X,'|')
1900	FORMAT(' ','|',2X,A26,9X,'|',8X,A1,7X,'|',5X,G13.6,4X,'|')
2000	FORMAT(' ','|',2X,A21,14X,'|',4X,A9,3X,'|',5X,G13.6,4X,'|')

	RETURN
	END


	SUBROUTINE PAGINA2

!	EXTERNAL NOM_SUB
	INTEGER NI,UNIF(20),I,ND,ND1,N
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	CHARACTER*8  FORMUNI(19)
	CHARACTER*8  NOM_SUB
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO

	ND=0
	DO 10 I=1,19,2
	   IF (UNIF(I).NE.0) THEN
	      ND=ND+1
	   END IF
  10	CONTINUE
	
	ND1=2*ND-1
!	DO 20 I=1,ND1,2
!	   FORMUNI(I)=NOM_SUB(UNIF(I),2)
!  20	CONTINUE

	WRITE(6,950)
	WRITE(6,900) 'PROPIEDAD','UNIDADES','VALOR'
	WRITE(6,950)
	WRITE(6,1400)
	WRITE(6,1000) 'Factor acentrico',W
	WRITE(6,1100) 'Radio de giro','m',RD
	WRITE(6,1200) 'Parametro de solubilidad','(J/m**3)**.5',PS
	WRITE(6,1300) 'Momento dipolar','C*m',DMU
	WRITE(6,1350) 'Tipo de componente',TIPO
	WRITE(6,1400) 
	WRITE(6,1600) 'CONFIGURACION UNIFAC'
	WRITE(6,1700)
	IF(ND.EQ.0) WRITE(6,1750)' No puede representarse '
	IF(ND.EQ.1) WRITE(6,1800) FORMUNI(1),UNIF(2)
	IF(ND.EQ.2) WRITE(6,1900) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4)
	IF(ND.EQ.3) WRITE(6,2000) FORMUNI(1),UNIF(2),FORMUNI(3),&
     	                          UNIF(4),FORMUNI(5),UNIF(6)
	IF(ND.EQ.4) WRITE(6,2100) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),&
                                 FORMUNI(5),UNIF(6),FORMUNI(7),UNIF(8)
	IF(ND.EQ.5) WRITE(6,2200) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),&
     	       FORMUNI(5),UNIF(6),FORMUNI(7),UNIF(8),FORMUNI(9),UNIF(10)
	IF(ND.EQ.6) THEN
	WRITE(6,2200) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),FORMUNI(5),&
                     UNIF(6),FORMUNI(7),UNIF(8),FORMUNI(9),UNIF(10)
	WRITE(6,1800) FORMUNI(11),UNIF(12)
	END IF
	IF(ND.EQ.7) THEN
	WRITE(6,2200) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),FORMUNI(5),&
                     UNIF(6),FORMUNI(7),UNIF(8),FORMUNI(9),UNIF(10)
	WRITE(6,1900) FORMUNI(11),UNIF(12),FORMUNI(13),UNIF(14)
	END IF
	IF(ND.EQ.8) THEN
	WRITE(6,2200) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),FORMUNI(5),&
                     UNIF(6),FORMUNI(7),UNIF(8),FORMUNI(9),UNIF(10)
	WRITE(6,2000) FORMUNI(11),UNIF(12),FORMUNI(13),UNIF(14),&
     	              FORMUNI(15),UNIF(16)
	END IF
	IF(ND.EQ.9) THEN
	WRITE(6,2200) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),FORMUNI(5),&
                     UNIF(6),FORMUNI(7),UNIF(8),FORMUNI(9),UNIF(10)
	WRITE(6,2100) FORMUNI(11),UNIF(12),FORMUNI(13),UNIF(14),&
     	              FORMUNI(15),UNIF(16),FORMUNI(17),UNIF(18)
	END IF
	IF(ND.EQ.10) THEN
	WRITE(6,2200) FORMUNI(1),UNIF(2),FORMUNI(3),UNIF(4),FORMUNI(5),&
                     UNIF(6),FORMUNI(7),UNIF(8),FORMUNI(9),UNIF(10)
	WRITE(6,2200) FORMUNI(11),UNIF(12),FORMUNI(13),UNIF(14),&
        FORMUNI(15),UNIF(16),FORMUNI(17),UNIF(18),FORMUNI(19),UNIF(20)
	END IF
	WRITE(6,1400)

 950	FORMAT(' ','|',37X,'|',16X,'|',22X,'|')
 900	FORMAT(' ','|',14X,A9,14X,'|',4X,A8,4X,'|',9X,A5,8X,'|')
1000	FORMAT(' ','|',2X,A16,19X,'|',16X,'|',5X,G13.6,4X,'|')
1100	FORMAT(' ','|',2X,A13,22X,'|',8X,A1,7X,'|',5X,G13.6,4X,'|')
1200	FORMAT(' ','|',2X,A24,11X,'|',2X,A11,3X,'|',5X,G13.6,4X,'|')
1300	FORMAT(' ','|',2X,A15,20X,'|',7X,A3,6X,'|',5X,G13.6,4X,'|')
1350	FORMAT(' ','|',2X,A18,17X,'|',16X,'|',9X,A4,9X,'|')
1400	FORMAT(80('-'))
1600	FORMAT(' ','|',29X,A20,28X,'|')
1700	FORMAT(' ','|',77X,'|')
1750    FORMAT(' ','|',27X,A24,26X,'|')
1800	FORMAT(' ','|',6X,A8,I3,60X,'|')
1900	FORMAT(' ','|',6X,2(A8,I3,3X),43X,'|')
2000	FORMAT(' ','|',6X,3(A8,I3,3X),29X,'|')
2100	FORMAT(' ','|',6X,4(A8,I3,3X),15X,'|')
2200	FORMAT(' ','|',6X,5(A8,I3,3X),X,'|')

	RETURN
	END


	SUBROUTINE PAGINA3

	INTEGER NI,UNIF(20),I
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO

	WRITE(6,1000) 'Presion de vapor','T min (',CIANT,')','T max (',CSANT,')'
	WRITE(6,1100)
	WRITE(6,1200) (ANT(I),I=1,5)
	WRITE(6,1100)
	WRITE(6,1300) 'Capacidad calorifica liquido','T min (',CICPL,')','T max (',CSCPL,')'
	WRITE(6,1100)
	WRITE(6,1200) (CPL(I),I=1,5)
	WRITE(6,1100)
 	WRITE(6,1400) 'Capacidad calorifica gas','T min (',CICPG,')','T max (',CSCPG,')'
	WRITE(6,1100)
	WRITE(6,1200) (CPG(I),I=1,5)
	WRITE(6,1100)
	WRITE(6,1500) 'Entalpia referencia gas',H0G
	WRITE(6,1600) 'Calor de vaporizacion',DHVAP
	WRITE(6,1150)
	WRITE(6,1100)

1000	FORMAT(' ','|',2X,A16,18X,A7,F7.2,A1,6X,A7,F7.2,A1,5X,'|')
1100	FORMAT(80('-'))
1150	FORMAT(' ','|',77X,'|')
1200	FORMAT(' ','|',4(G13.6,1X,'|',1X),G13.6,'|')
1300	FORMAT(' ','|',2X,A28,6X,A7,F7.2,A1,6X,A7,F7.2,A1,5X,'|')
1400	FORMAT(' ','|',2X,A24,10X,A7,F7.2,A1,6X,A7,F7.2,A1,5X,'|')
1500	FORMAT(' ','|',2X,A23,15X,G13.6,24X,'|')
1600	FORMAT(' ','|',2X,A21,17X,G13.6,24X,'|')

	RETURN
	END

SUBROUTINE NUM_PRO

      IMPLICIT real*8 (A-H,O-Z)
	INTEGER NI,UNIF(20),I,I1,I2,I3,I4,NO,OB,PARINT,BD
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	COMMON/OPC/NO,OB,PARINT,BD
	CHARACTER fs*8,MG*8,caract*3
	COMMON/GRUPOSRAM/nlgr,nl,MG,fs,item,caract,mm,&
     						  mj,mk,mi,mh,rpar,qpar,&
     						  pmg,deltc,delpc,delv,delvi,deln
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO

      IF (OB.EQ.1) THEN
	WRITE(6,*)'  1 - NUMERO IDENTIFICATORIO  = ',NI
	WRITE(6,200)'   2 - NOMBRE COMPUESTO        = ',ANAME
	WRITE(6,210)'   3 - FORMULA QUIMICA         = ',FORM1
	WRITE(6,200)'   4 - FORMULA ESTRUCTURAL     = ',FORM2
 200	FORMAT(A33,A35)
 210	FORMAT(A33,A20)
	WRITE(6,*)' '
	WRITE(6,*)'        CONFIGURACION UNIFAC'
	WRITE(6,*)' '
	DO 30 I=1,10
	   I1=2*I-1
	   I2=2*I
	   I3=I+4
	   WRITE(6,220)'  ',I3,' - GRUPO ',I,'   = ',UNIF(I1),UNIF(I2)
 220	   FORMAT(A2,I2,A9,I2,A5,5X,I3,5X,I3)
  30    CONTINUE

	WRITE(6,*)'  '
	WRITE(6,*)' 15 - TEMPERATURA CRITICA     = ',TC

	CALL PAUSA

	WRITE(6,*)' 16 - PRESION CRITICA         = ',PC
	WRITE(6,*)' 17 - VOLUMEN CRITICO         = ',VC
	WRITE(6,*)' 18 - FACTOR COMPRES. CRITICO = ',ZC
	WRITE(6,*)' 19 - PUNTO DE FUSION         = ',TF
	WRITE(6,*)' 20 - PUNTO NORMAL EBULLICION = ',TEB
	WRITE(6,*)' 21 - VOL. MOLAR LIQUIDO      = ',VLIQ
	WRITE(6,*)' 22 - FACTOR ACENTRICO        = ',W
	WRITE(6,*)' 23 - PARAMETRO SOLUBILIDAD   = ',PS
	WRITE(6,*)' 24 - PESO MOLECULAR          = ',PM
	WRITE(6,*)' 25 - RADIO DE GIRO           = ',RD
	WRITE(6,*)' 26 - MOMENTO DIPOLAR         = ',DMU
	WRITE(6,*)' 27 - TIPO DE COMPONENTE      = ',TIPO 	  	
	WRITE(6,*)'  '
	WRITE(6,*)'     PRESION DE VAPOR'
	DO 40 I=1,5
	   I4=I+27
	   WRITE(6,250)'     ',I4,' - COEFIC. ',I,'    = ',ANT(I)
  40 	CONTINUE
 250    FORMAT(A5,I2,A10,I1,A6,D13.6)
	WRITE(6,*)'    33 - TEMPERATURA MIN     = ',CIANT
	WRITE(6,*)'    34 - TEMPERATURA MAX     = ',CSANT

	CALL PAUSA

	WRITE(6,*)'     CAPACIDAD CALORIFICA GAS'
	WRITE(6,*)'  '
	DO 50 I=1,5
	   I4=I+34
	   WRITE(6,250)'     ',I4,' - COEFIC. ',I,'    = ',CPG(I)   
  50	CONTINUE 
	WRITE(6,*)'    40 - TEMPERATURA MIN     = ',CICPG
	WRITE(6,*)'    41 - TEMPERATURA MAX     = ',CSCPG
	WRITE(6,*)'  '
	WRITE(6,*)'     CAPACIDAD CALORIFICA LIQUIDO'
	WRITE(6,*)'   '
	DO 60 I=1,5
	   I4=I+41
	   WRITE(6,250)'     ',I4,' - COEFIC. ',I,'    = ',CPL(I)
  60	CONTINUE
	WRITE(6,*)'    47 - TEMPERATURA MIN     = ',CICPL
	WRITE(6,*)'    48 - TEMPERATURA MAX     = ',CSCPL
	WRITE(6,*)'  '
	WRITE(6,*)' 49 - ENTALPIA REFERENCIA GAS   = ',H0G
	WRITE(6,*)' 50 - CALOR DE VAPORIZACION     = ',DHVAP

	CALL PAUSA
      ELSE IF (OB.EQ.3) THEN
  	  write (6,7)nlgr,fs,caract,mm,mj,mk,mi,mh
        call pausa
        write (6,8)rpar,qpar,pmg,deltc,delpc,delv,delvi,deln
        call pausa
	ENDIF
!C
!C-----FORMATOS  
  7	format (2x,'Group ',i3,':  ',a8,/,&
     		'  1- Caract:',54x,a3,//,&
     		'Combination properties',/,&
             '  2- Number of m attachments:',39x,i2,/,&
             '  3- Number of j attachments:',39x,i2,/,&
             '  4- Number of k attachments:',39x,i2,/,&
             '  5- Number of i attachments:',39x,i2,/,&
             '  6- Number of h attachments:',39x,i2,//)
 
 8   	format ('Group contributions to properties estimates',/,&
     		'  7- Rpar:',44x,d15.8,/,&
     		'  8- Qpar:',44x,d15.8,/,&
     		'  9- Pmg :',44x,d15.8,/,&
     		' 10- Deltc:',43x,d15.8,/,&
     		' 11- Delpc:',43x,d15.8,/,&
     		' 12- Delv:',44x,d15.8,/,&
     		' 13- Delvi:',43x,d15.8,/,&
     		' 14- Deln:',44x,d15.8,//,70x)

	RETURN
ENDsubroutine



	SUBROUTINE IN_PRO(NP)

	INTEGER NI,UNIF(20),NP,I,I1,I2,I3,NO,OB,PARINT,BD
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	COMMON/OPC/NO,OB,PARINT,BD
	character fs*8,MG*8,caract*3
	COMMON/GRUPOSRAM/nlgr,nl,MG,fs,item,caract,mm,&
     						  mj,mk,mi,mh,rpar,qpar,&
     						  pmg,deltc,delpc,delv,delvi,deln
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO
!C
	IF (OB.EQ.1) THEN
	  IF(NP.EQ.1) THEN
	      WRITE(6,*)'NUMERO IDENTIFICATORIO ='
	      READ(5,*) NI
	  END IF
      	IF(NP.EQ.2) THEN
	    WRITE(6,*) 'NOMBRE COMPUESTO= '
	    READ(5,200)   ANAME
	  END IF
	  IF(NP.EQ.3) THEN
	      WRITE(6,*) 'FORMULA QUIMICA= '
	      READ(5,200)   FORM1
	  END IF
	  IF(NP.EQ.4) THEN
	      WRITE(6,*) 'FORMULA ESTRUCTURAL= '
	      READ(5,200)   FORM2
	  END IF
 200   	FORMAT(A)
	  DO 32 I=1,10
	   I1=I+4
	   I2=I*2-1
	   I3=2*I
	   IF(NP.EQ.I1) THEN
	     WRITE(6,*)'GRUPO ',I,' (ident.,cantidad) ='
	     READ(5,*) UNIF(I2),UNIF(I3)
	   END IF
  32	  CONTINUE
	  IF(NP.EQ.15) THEN
	      WRITE(6,*) 'TEMPERATURA CRITICA='
	      READ(5,*)   TC
	  END IF
	  IF(NP.EQ.16) THEN
	      WRITE(6,*) 'PRESION CRITICA='
	      READ(5,*)   PC
	  END IF
	  IF(NP.EQ.17) THEN
	      WRITE(6,*) 'VOLUMEN CRITICO='
	      READ(5,*)   VC
	  END IF
	  IF(NP.EQ.18) THEN
	      WRITE(6,*) 'FACTOR DE COMPRESIBILIDAD='
	      READ(5,*)   ZC
	  END IF
	  IF(NP.EQ.19) THEN
	      WRITE(6,*) 'PUNTO DE FUSION='
	      READ(5,*)   TF
	  END IF
	  IF(NP.EQ.20) THEN
	      WRITE(6,*) 'PUNTO NORMAL DE EBULLICION='
	      READ(5,*)   TEB
	  END IF
	  IF(NP.EQ.21) THEN
	      WRITE(6,*) 'VOLUMEN MOLAR LIQUIDO='
	      READ(5,*)   VLIQ
	  END IF
	  IF(NP.EQ.22) THEN
	      WRITE(6,*) 'FACTOR ACENTRICO='
	      READ(5,*)   W
	  END IF
	  IF(NP.EQ.23) THEN
	      WRITE(6,*) 'PARAMETRO DE SOLUBILIDAD='
	      READ(5,*)   PS
	  END IF
	  IF(NP.EQ.24) THEN
	      WRITE(6,*)'PESO MOLECULAR ='
	      READ(5,*) PM 
	  END IF
	  IF(NP.EQ.25) THEN
	      WRITE(6,*) 'RADIO DE GIRO ='
	      READ(5,*) RD
	  END IF
	  IF(NP.EQ.26) THEN
	      WRITE(6,*)'MOMENTO DIPOLAR ='
	      READ(5,*) DMU
	  END IF
	  IF(NP.EQ.27) THEN
	      WRITE(6,*)'TIPO DE COMPONENTE ='
	      READ(5,222) TIPO
 222	  FORMAT(A4)
	  END IF
	  DO 36 I=1,5
	   I1=I+27
	   I2=I+34
	   I3=I+41
	   IF(NP.EQ.I1) THEN
	     WRITE(6,500)'COEF. ANT ',I,' ='
 500	     FORMAT(' ',A10,3X,I1,A2)
	     READ(5,*) ANT(I)
	   END IF
	   IF(NP.EQ.I2) THEN
	     WRITE(6,510)'COEF. CP GAS ',I,' ='
 510	     FORMAT(' ',A13,2X,I1,A2) 
	     READ(5,*) CPG(I)
	   END IF
	   IF(NP.EQ.I3) THEN
	     WRITE(6,520)'COEF. CP LIQUIDO ',I,' ='
 520	     FORMAT(' ',A17,2X,I1,A2)
	     READ(5,*) CPL(I)
	   END IF
  36	  CONTINUE
	  IF(NP.EQ.33) THEN
	      WRITE(6,*)'TEMP. MIN. ANT. ='
	      READ(5,*) CIANT
	  END IF
	  IF(NP.EQ.34) THEN
	      WRITE(6,*)'TEMP. MAX. ANT ='
	      READ(5,*) CSANT
	  END IF
	  IF(NP.EQ.40) THEN
	      WRITE(6,*)'TEMP. MIN. CP GAS ='
	      READ(5,*) CICPG
	  END IF
	  IF(NP.EQ.41) THEN
	      WRITE(6,*) 'TEMP. MAX. CP GAS ='
	      READ(5,*) CSCPG
	  END IF
	  IF(NP.EQ.47) THEN
	      WRITE(6,*)'TEMP. MIN. CP LIQ. ='
	      READ(5,*) CICPL
	  END IF
	  IF(NP.EQ.48) THEN
	      WRITE(6,*)'TEMP. MAX. CP LIQ. ='
	      READ(5,*) CSCPL
	  END IF
	  IF(NP.EQ.49) THEN
	      WRITE(6,*) 'ENTALPIA REFERENCIA GAS ='
	      READ(5,*) H0G
	  END IF
	  IF(NP.EQ.50) THEN
	      WRITE(6,*)'CALOR DE VAPORIZACION ='
	      READ(5,*) DHVAP
	  END IF
      ELSE IF (OB.EQ.3) THEN
        if (np.eq.1) then 
 114			write (6,8)fs
			write (6,3)
			read (5,*) caract
			if (caract.ne.'AR1'.and.caract.ne.'GR1'.and.&
     		caract.ne.'DV1'.and.caract.ne.'DV2'.and.caract&
     		.ne.'SV2'.and.caract.ne.'SV1'.and.caract&
     		.ne.'0'.and.caract.ne.'3V1'.and.caract.ne.'4V1') then	
				goto 114
			end if
		else if (NP.eq.2) then 
			write (6,*) '  Give the number of m attachments for group ',fs
			write (6,3)
			read (5,*) mm
		else if (NP.eq.3) then 
			write (6,*) '  Give the number of j attachments for group ',fs
			write (6,3)
			read (5,*) mj
		else if (NP.eq.4) then 
			write (6,*) '  Give the number of k attachments for group ',fs
			write (6,3)
			read (5,*) mk
		else if (NP.eq.5) then 
			write (6,*) '  Give the number of i attachments for group ',fs
			write (6,3)
			read (5,*) mi
		else if (NP.eq.6) then 
			write (6,*) '  Give the number of h attachments for group ',fs
			write (6,3)
			read (5,*) mh
		else if (NP.eq.7) then 
			write (6,*) '  Give the new value of Rpar for group ',fs
			write (6,3)
			read (5,*) rpar
		else if (NP.eq.8) then 
			write (6,*) '  Give the new value of Qpar for group ',fs
			write (6,3)
			read (5,*) qpar
		else if (NP.eq.9) then 
			write (6,*) '  Give the new value of Pmg for group ',fs
			write (6,3)
			read (5,*) pmg
		else if (NP.eq.10) then 
			write (6,*) '  Give the new value of Deltc for group ',fs
			write (6,3)
			read (5,*) deltc
		else if (NP.eq.11) then 
			write (6,*) '  Give the new value of Delpc for group ',fs
			write (6,3)
			read (5,*) delpc
		else if (NP.eq.12) then 
			write (6,*) '  Give the new value of Delv for group ',fs
			write (6,3)
			read (5,*) delv
		else if (NP.eq.13) then 
			write (6,*) '  Give the new value of Delvi for group ',fs
			write (6,3)
			read (5,*) delvi
		else 
			write (6,*) '  Give the new value of Deln for group ',fs
			write (6,3)
			read (5,*) deln
		end if
      ENDIF
!C
!C-----FORMATOS
   3  format (1x,/,60x,'> ',$)
   8	format (' Give the new caract to use the group ',a8,&
     		/,'  only in the synthesis of:',&
     		/,'  Aromatic solvents:',42x,'AR1',&
     		/,'  Single substance groups:',36x,'GR1',&
     		/,'  Aliphatic and Mixed aromatic-aliphatic solvents: '&
     		 ,11x,'DV1',/,'  Aliphatic, Mixed aromatic-aliphatic and',&
     		  ' Cyclic solvents:',4x,'DV2',&
     		/,' or as a terminal group:',&
     		/,'  when there is not any dual valence group in its ',&
     		  'maingroup:',2X,'SV2',/,'  when there are any dual valence ',&
     		  'group in its maingroup:',5x,'SV1',/,&
     		  ' If the group will not be taken into account:',19x,'0')

	RETURN
	END         
!C
      SUBROUTINE BUSCAR (I,IPAREQL,FF,NUM)
      CHARACTER*8 FS,FF,MG
      LOGICAL IPAREQL
!C
      IPAREQL=.FALSE.
      num=(I-1)*150
      fs= ''
      MG= ''
      DO WHILE (MG.NE.'fin')
	  num=num+1
	  read (14,10,rec=num)MG,fs
	  IF (FS.EQ.FF) THEN
	      IPAREQL=.TRUE.
	      EXIT
	  ENDIF
	enddo
!c
!c
  10  format (4x,2a8) 
      return
      end
!C
      SUBROUTINE MOSTRAR_TABLA (ipareq,OPT)
!C
      IMPLICIT real*8 (A-H,O-Z)
      integer,parameter :: ALTOP=9,ANCHOP=5
      CHARACTER MG*8
	  character*8,dimension(ANCHOP)::MGN
      INTEGER NGRMAX,ANCHO,ALTO,PAGN,PAGL,RESN,RESL
      INTEGER OPT,OPT1,MGNAN(ANCHOP),MGNAL(ALTOP),NO,OB,PARINT,BD
      real*8 FILA(ANCHOP)
      COMMON/OPC/NO,OB,PARINT,BD
!C
 11   ngrmax = max_sub_int (ipareq)
      ANCHO=ANCHOP
      ALTO=ALTOP
      RESN = MOD (NGRMAX,ANCHO)
      RESL = MOD (NGRMAX,ALTO)
      PAGN = NGRMAX/ANCHO
      PAGL = NGRMAX/ALTO
      DO L=1,PAGL+1
        ANCHO=ANCHOP
        IF(L.EQ.PAGL+1) ALTO=RESL
        DO K=1,PAGN+1
            IF(K.EQ.PAGN+1) ANCHO=RESN
            DO H=1,ANCHO
                F=(K-1)*ANCHOP+H
                MGNAN(H)=F
                irec =F+(ipareq-1)*70
                READ (BD,201,rec=irec)MGN(H)
            ENDDO
            WRITE (6,104) MGNAN(1:ANCHO)
            WRITE (6,101) MGN(1:ANCHO)
            DO J=1,ALTO
                DO I=1,ANCHO
                    F=(K-1)*ANCHOP+I
                    N=(L-1)*ALTOP+J
                    irec =N+(ipareq-1)*70
                    read (BD,200,rec=irec)MG,(AINT,M=1,F)
                    FILA(I)= AINT
                ENDDO
                WRITE (6,100)N, MG, FILA(1:ANCHO)
            ENDDO
            IF (NO.EQ.1) THEN
 13             WRITE (6,105)
                READ (5,202,ERR=13) OPT
                IF (OPT.GT.4.OR.OPT.LT.2) GOTO 13
            ELSE
 10             WRITE (6,102)
                READ (5,202,ERR=10) OPT
                IF (OPT.GT.4.OR.OPT.LT.1) GOTO 10
            ENDIF
            IF (OPT.EQ.3) THEN
                GOTO 11
            ELSE IF (OPT.LE.2) THEN
                GOTO 1000
            ENDIF
        ENDDO
      ENDDO
      IF (NO.EQ.1) THEN
 14     WRITE (6,106) 
        READ (5,202,ERR=14) OPT
        IF (OPT.GT.3.OR.OPT.LT.2) GOTO 14 
      ELSE    
 12     WRITE (6,103) 
        READ (5,202,ERR=12) OPT
        IF (OPT.GT.3.OR.OPT.LT.1) GOTO 12
      ENDIF
      IF (OPT.EQ.3) GOTO 11
!C-----FORMATOS
 100  FORMAT (/,I2,1X,a8,2x,5d13.5)
 101  FORMAT (/,13X,5(5X,A8))
 102  FORMAT (' ',//,10X,'SELECT: 1',5X,'EXIT: 2',5X,'BEGIN AGAIN: 3',5X,'CONTINUE: 4',5X,'> ',$)
 103  FORMAT (' ',//,10X,'SELECT: 1',5X,'EXIT: 2',5X,'BEGIN AGAIN: 3',5X,'> ',$)
 104  FORMAT (///,13X,5(10X,I3))
 105  FORMAT (' ',//,24X,'EXIT: 2',5X,'BEGIN AGAIN: 3',5X,'CONTINUE: 4',5X,'> ',$)
 106  FORMAT (' ',//,24X,'EXIT: 2',5X,'BEGIN AGAIN: 3',5X,'> ',$)
 200  FORMAT (a8,70d12.5)
 201  FORMAT (a8)
 202  FORMAT (I1)
1000  RETURN
      END
!C
!C
      SUBROUTINE MODIFICAR_INTRCN (M,N,IPAREQ)
      PARAMETER (NG=70)
!C      
      IMPLICIT real*8 (A-H,O-Z)
      DIMENSION AINTV(NG)
      CHARACTER*8 MGM,MGN
      CHARACTER*1 RESP, RESP1
      INTEGER NO, OB,PARINT,BD
	COMMON/OPC/NO,OB,PARINT,BD
!C
      irec = M+(ipareq-1)*70
      irecN= N+(ipareq-1)*70
      READ (BD,10,rec=irec)MGM,(AINTV(F),F=1,NG)
      READ (BD,11,rec=irecN)MGN
 101  WRITE (6,50)MGN,MGM,AINTV(N)
      WRITE (6,54)
 102  READ (5,12) RESP1
      IF (RESP1.NE.'y'.AND.RESP1.NE.'n') GOTO 102
      IF (RESP1.EQ.'n') GOTO 1000
      WRITE (6,51)
      READ (5,*) AINTV(N)
      WRITE (6,52)MGN,MGM,AINTV(N)
      WRITE (6,53)
 100  READ (5,12) RESP
      IF (RESP.NE.'y'.AND.RESP.NE.'n') GOTO 100
      IF (RESP.EQ.'n') GOTO 101
      WRITE (BD,10,rec=irec)MGM,(AINTV(F),F=1,NG)
!C-----FORMATOS
  10  FORMAT (a8,70d12.5)
  11  FORMAT (a8)
  12  FORMAT (a1)
  50  FORMAT (//,2X,'Current Value =',//15x,a8,//,1x,a8,2x,d12.5)
  51  FORMAT (//,2X,'Enter new value = ',$) 
  52  FORMAT (///,2X,'Input value =',//15x,a8,//,1x,a8,2x,d12.5) 
  53  FORMAT (///,2X,'Is That Correct? (y/n)',//,52x,'> ',$)
  54  FORMAT (/,2X,'Would you like to modify this? (y/n)',//,52x,'> ',$)
1000  RETURN
      END
      
SUBROUTINE LEC_COM(NR)

	INTEGER UNIF(20),K,K1,K2,K3
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),&
     	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,&
                  W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,&
                  CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO

	IF (NR.GT.0) THEN
	  READ(7,100,REC=NR,ERR=20) NI,ANAME,FORM1,FORM2,&
                                   (UNIF(K),K=1,20),TC,PC,VC,ZC,TF,TEB,&
                                   VLIQ,W,PS,PM,RD,DMU,TIPO,&
                                   (ANT(K1),K1=1,5),(CPG(K2),K2=1,5),&
                                   (CPL(K3),K3=1,5),CIANT,CSANT,CICPG,&
                                   CSCPG,CICPL,CSCPL,H0G,DHVAP
 100 	  FORMAT(I4,A35,A20,A35,20I3,12D13.6,A4,23D13.6)
	  write(6,*)'ni=',ni
	END IF
	goto 500

  20	WRITE(6,*) '*ERROR* EN LA LECTURA DE "PROP.PDB"'

	call pausa

 500	RETURN 
	END