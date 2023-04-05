	SUBROUTINE IN_COM
C
      IMPLICIT real*8 (A-H,O-Z)
      INTEGER NO,OB,PARINT,BD,NGRMAX
      character resp*1,MGC*8
	COMMON/OPC/NO,OB,PARINT,BD
C
      IF (OB.EQ.1) THEN
        CALL IN_PROP
      ELSE IF (OB.EQ.2) THEN
        CALL IN_INTRCN
      ELSE IF (OB.EQ.3) THEN
        resp='n'
        CALL IN_GRUPOS (resp,ipareq)
      ENDIF
      
      RETURN
	END 
C
      SUBROUTINE IN_PROP
C      
      INTEGER NR,NI,UNIF(20),K1,K2,IMP,K3,I,NP
	real*8  TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,ANT(5),CPG(5),
     *	        CPL(5),CIANT,CSANT,CICPG,CSCPG,CICPL,CSCPL,H0G,DHVAP
 	CHARACTER*35 ANAME,FORM2
	CHARACTER*20 FORM1
	CHARACTER*4  TIPO
	COMMON/NOP3/K2
	COMMON/PROP/UNIF,ANT,CPG,CPL,TC,PC,VC,ZC,TF,TEB,VLIQ,
     *             W,PS,PM,RD,DMU,CIANT,CSANT,CICPG,CSCPG,CICPL,
     *             CSCPL,H0G,DHVAP,NI,ANAME,FORM2,FORM1,TIPO
C     
 300	WRITE(6,*)'  '
	WRITE(6,*)' NUMERO IDENTIFICATORIO DEL COMPONENTE ',
     *            ' QUE DESEA INGRESARSE = '
	READ(5,*) NR
C 
	CALL LEC_COM(NR)
C
	IF(NI.NE.0) THEN
	  WRITE(6,*)'  '
	  WRITE(6,260)' YA EXISTE UN COMPONENTE CON NUMERO '
          WRITE(6,270)' IDENTIFICATORIO = ',NR,' EN "PROP.PDB"'  
 260      FORMAT(10X,A35)
 270	  FORMAT(10X,A19,I4,A15)
	  GOTO 300
	END IF

	NI=NR
C
	WRITE(6,*)'  '
	WRITE(6,*)' INGRESAR LAS SIGUIENTES PROPIEDADES :'
	WRITE(6,*)'  '
	WRITE(6,*) 'NOMBRE COMPUESTO = '
	READ(5,200)   ANAME
	WRITE(6,*) 'FORMULA QUIMICA = '
	READ(5,210)   FORM1
	WRITE(6,*) 'FORMULA ESTRUCTURAL = '
	READ(5,200)   FORM2
 200  FORMAT(A35)
 210	FORMAT(A20)
	WRITE(6,*) '                                               '
	WRITE(6,*) 'Para la CONFIGURACION UNIFAC pueden ingresarse ',
     *	           '10 grupos como maximo. Los mismos seran pedidos',
     *	           ' de a uno a continuacion.'
	WRITE(6,*) '                       '
	WRITE(6,*) 'NUMERO TOTAL DE GRUPOS DISTINTOS ='
	READ(5,*)   NG
	WRITE(6,*) '                             '
	DO 12 I=1,NG
	   WRITE(6,*) 'GRUPO ',I,
     *   ' = (identificacion, cantidad de este grupo)'
	   I1=2*I-1
	   I2=2*I
	   READ(5,*)  UNIF(I1),UNIF(I2)
  12	CONTINUE
	NG2=NG+1
	DO 14 I=NG2,10
	   I1=2*I-1
	   I2=2*I
	   UNIF(I1)=0
	   UNIF(I2)=0
  14  CONTINUE
	WRITE(6,*) 'TEMPERATURA CRITICA ='
	READ(5,*)   TC
	WRITE(6,*) 'PRESION CRITICA ='
	READ(5,*)   PC
	WRITE(6,*) 'VOLUMEN CRITICO ='
	READ(5,*)   VC
	WRITE(6,*) 'FACTOR DE COMPRESIBILIDAD ='
	READ(5,*)   ZC
	WRITE(6,*) 'PUNTO DE FUSION ='
	READ(5,*)   TF
	WRITE(6,*) 'PUNTO NORMAL DE EBULLICION ='
	READ(5,*)   TEB
	WRITE(6,*) 'VOLUMEN MOLAR LIQUIDO ='
	READ(5,*)   VLIQ
	WRITE(6,*) 'FACTOR ACENTRICO ='
	READ(5,*)   W
	WRITE(6,*) 'PARAMETRO DE SOLUBILIDAD ='
	READ(5,*)   PS
	WRITE(6,*) 'PESO MOLECULAR ='
	READ(5,*)   PM
	WRITE(6,*) 'RADIO DE GIRO ='
	READ(5,*)   RD
	WRITE(6,*) 'MOMENTO DIPOLAR ='
	READ(5,*)   DMU
	WRITE(6,*) 'TIPO DE COMPONENTE ='
	READ(5,222)   TIPO
 222	FORMAT(A4)
	WRITE(6,*) '  '
	WRITE(6,*) 'Para la PRESION DE VAPOR deben ingresarse las ',
     *             'cinco constantes de ANTOINE (se pediran de a una) ',
     *	           'y las cotas inferior y superior'
	WRITE(6,*)' '
	DO 75 I=1,5
	   WRITE(6,*)'   COEFF. ',I,' ='
	   READ(5,*)  ANT(I)
  75	CONTINUE
	WRITE(6,*) 'TEMP. MIN.='
	READ(5,*) CIANT
	WRITE(6,*) 'TEMP. MAX.=' 
	READ(5,*) CSANT
	WRITE(6,*)'  '
	WRITE(6,*)'Ingresar los coeficientes para el calculo de ',
     *            'la CAPACIDAD CALORIFICA DEL GAS de a uno y las ',
     *            'temperaturas min. y max.'
	WRITE(6,*)'  '
	DO 76 I=1,5
	   WRITE(6,*)'   COEF. ',I,' ='
	   READ(5,*) CPG(I)
  76	CONTINUE
	WRITE(6,*)'TEMP. MIN.='
	READ(5,*) CICPG
	WRITE(6,*)'TEMP. MAX.='
	READ(5,*) CSCPG
	WRITE(6,*)'  '
	WRITE(6,*)'Ingreasar los coeficientes para el calculo de ',
     *            'la CAPACIDAD CALORIFICA DEL LIQUIDO de a uno y ',
     *            'las temperaturas min. y max.'
	WRITE(6,*)'  '
	DO 78 I=1,5
	   WRITE(6,*)'   COEF. ',I,' ='
	   READ(5,*)  CPL(I)
  78	CONTINUE
	WRITE(6,*)'TEMP. MIN.='
	READ(5,*)  CICPL
	WRITE(6,*)'TEMP. MAX.='
	READ(5,*) CSCPL
	WRITE(6,*)' '
	WRITE(6,*)'ENTALPIA DE REFERNCIA DEL GAS ='
	READ(5,*) H0G
	WRITE(6,*)'CALOR DE VAPORIZACION ='
	READ(5,*) DHVAP
C
 600	Continue
C
 	WRITE(6,*)'  '
	WRITE(6,*)' LOS DATOS INGRESADOS SON LOS SIGUIENTES :'
	WRITE(6,*)'  '
C
	CALL NUM_PRO
C
	WRITE(6,*)'  '
	WRITE(6,*)' SON TODOS LOS DATOS CORRECTOS ? (1=SI,0=NO)'
	READ(5,*)K1
	IF(K1.EQ.1) GOTO 700	
	WRITE(6,*)' '
	WRITE(6,*)' NUMERO TOTAL DE DATOS A CORREGIR ='
	READ(5,*) K3
	DO 20 I=1,K3
	   WRITE(6,*)' ' 
	   WRITE(6,*)' PROPIEDAD A CORREGIR ( 1-50 ) ='
	   READ(5,*) NP
C
	   CALL IN_PRO(NP)
C  
20	CONTINUE
	GOTO 600
 700	WRITE(IMP,100,REC=NI) NI,ANAME,FORM1,FORM2,(UNIF(K),K=1,20),
     *                      TC,PC,VC,ZC,TF,TEB,VLIQ,W,PS,PM,RD,DMU,TIPO,
     *                      (ANT(KI),KI=1,5),(CPG(KJ),KJ=1,5),
     *                      (CPL(KL),KL=1,5),CIANT,CSANT,CICPG,CSCPG,
     *                      CICPL,CSCPL,H0G,DHVAP
 100 	 FORMAT(I4,A35,A20,A35,20I3,12D13.6,A4,23D13.6)
c	CLOSE(UNIT=IMP)
C 
	WRITE(6,*)' '
	WRITE(6,*)' INDICAR :'
	WRITE(6,*)'     1- INGRESAR OTRO COMPONENTE'
	WRITE(6,*)'     2- ACCEDER A OTRA OPCION'
	WRITE(6,*)'     3- TERMINAR'
	WRITE(6,*)'   '
	READ(5,*) K2
	IF(K2.EQ.1) GOTO 300
c	CLOSE(UNIT=IMP)
	GOTO 500
C 
 10	WRITE(6,*)' *ERROR* EN LA APERTURA DE "PROP.PDB"'
 500	RETURN
	END 
C
	SUBROUTINE IN_INTRCN
C
      PARAMETER (NG=70)
      IMPLICIT real*8 (A-H,O-Z)
      CHARACTER*8 MGC,MG,MGV(NG),FS
      character*1 resp
      INTEGER GRUPOS(NG)
      DIMENSION PINTER(NG)
      COMMON/GRUPOSRAM/nlgr,nl,MG,fs,item,caract,mm,
     *						  mj,mk,mi,mh,rpar,qpar,
     *						  pmg,deltc,delpc,delv,delvi,deln
C
      CALL TABLA_PARAMETROS(IPAREQ)
      write (6,110)
	read (5,201) mgc
	call ajustar_nombre (mgc,mg)
C-----GENERA NUEVA LINEA
      CALL GENERA_LINEA_INTRCN (IPAREQ,IREC)
C-----AGREGA EL NUEVO GRUPO
      read (13,200,rec=irec)    MGC,(PINTER(F),F=1,NG)
      WRITE (13,200,REC=IREC)   MG,(PINTER(F),F=1,NG)
C-----AGREGA PARÁMETROS DE INTERACCIÓN
      CALL CARGA_PARAMETROS_INTRCN (IPAREQ)
      WRITE(6,120)
      CALL PAUSA
      FS=MG
      RESP='y'
      CALL IN_GRUPOS (resp,ipareq)
      
C-----FORMAT
 110  FORMAT (/,'Give name of main group: ',$)
 120  FORMAT (//,'*** Now need to modify GRUPOSRAM.MDS ***')
 200  FORMAT (a8,70d12.5)
 201  FORMAT (a8)
      RETURN
      END	
C
      SUBROUTINE IN_GRUPOS (resp,ipareq)
C
      use SubGrupos
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (NMG=150,NG=70)
      CHARACTER*8 MG,FS,mgv(NMG),FSC,MGC
      CHARACTER*3 CARACT,caractc
      CHARACTER*1 resp1,resp,resp2
      INTEGER GRUPOS(NMG)
      LOGICAL lgroup
      DIMENSION PINTER(NG)
	COMMON/GRUPOSRAM/nlgr,nl,MG,fs,item,caract,mm,
     *						  mj,mk,mi,mh,rpar,qpar,
     *						  pmg,deltc,delpc,delv,delvi,deln
  
      idev = 6
      idevr= 5  
      IF (resp.eq.'y') THEN
		call GENERA_LINEA_GRUPOSRAM (MG,ipareq,IREC3,resp)
		resp1=resp
		GOTO 13
      ENDIF
      CALL TABLA_PARAMETROS(IPAREQ)
      
  10  write (6,*) 'Is a main group (y/n)?'
	read (5,107,ERR=10) resp1
 	if (resp1.eq.'y') then
		write (6,*) 'Give name of main group'
		read (5,200) mgc
		call ajustar_nombre (mgc,mg)
		write (6,*) 'Give name of new sub-group'
		read (5,200) fsc
		call ajustar_nombre (fsc,fs)	
		call GENERA_LINEA_GRUPOSRAM (MG,ipareq,IREC3,resp1)	
      else if(resp1.eq.'n') then
!----lectura de grupos principales
        k=0
		do while (mg.ne.'fin')
			k=k+1
			irec2=k+(ipareq-1)*70
			read (13,200,rec=irec2) mg
			mgv(k)=mg
			grupos(k)=k 
		enddo 
C----Selecciona el grupos principal al que pertenece el subgrupo
		write (6,102)
		call Groups_Present ()
  12	  write (6,3)
		read (5,*,err=12)i1
		if (i1.gt.k.or.i1.lt.0)goto 12
        i=0
        do
            i=i+1
            irec3=i+(ipareq-1)*150
            read (14,203,rec=irec3)ic,mg
            if(i1.eq.ic) exit
        enddo
C----Nombre del nuevo Subgrupo
		write (6,*) 'Give name of new sub-group'
		read (5,200) fsc
        call ajustar_nombre(fsc,fs)
!-------Genera una nueva linea donde se escribira el nuevo grupo
        call GENERA_LINEA_GRUPOSRAM (MG,ipareq,IREC3,resp1)
      else
        goto 10
      endif
C----Pide al usuario los valores de las propiedades del nuevo grupo
 13	write (6,*) 'Give the following values:'
      write (6,*) ' '
 17	write (6,8) MG
	write (6,3)
	read (5,103) caract
	if (caract.ne.'AR1'.and.caract.ne.'GR1'.and.
     *caract.ne.'DV1'.and.caract.ne.'DV2'.and.caract
     *.ne.'SV2'.and.caract.ne.'SV1'.and.caract
     *.ne.'0'.and.caract.ne.'3V1'.and.caract.ne.'4V1') goto 17
	write (6,*) 'number of m attachments'
	write (6,3)
	read (5,104) mm
	write (6,*) 'number of j attachments'
	write (6,3)
	read (5,104) mj
	write (6,*) 'number of k attachments'
	write (6,3)
	read (5,104) mk
	write (6,*) 'number of i attachments'
	write (6,3)
	read (5,104) mi 
	write (6,*) 'number of h attachments'
	write (6,3)
	read (5,104) mh
	write (6,*) 'Rpar '
	write (6,3)
	read (5,*) rpar
	write (6,*) 'Qpar '
	write (6,3)
	read (5,*) qpar
	write (6,*) 'Pmg '
	write (6,3)
	read (5,*) pmg
	write (6,*) 'Deltc '
	write (6,3)
	read (5,*) deltc
	write (6,*) 'Delpc '
	write (6,3)
	read (5,*) delpc
	write (6,*) 'Delv '
	write (6,3)
	read (5,*) delv
	write (6,*) 'Delvi '
	write (6,3)
	read (5,*) delvi
	write (6,*) 'Deln '
	write (6,3)
	read (5,*) deln
c-----Escribe el nuevo grupo y sus valores en la tabla de parámetros correspondiente
      if (resp1.eq.'y') then
        read (14,105,rec=irec3-1)i1
        i1=i1+1
     	  write (14,201,rec=irec3)i1,mg,fs,item,item,caract,mm,
     *				          mj,mk,mi,mh,rpar,qpar,
     *						  pmg,deltc,delpc,delv,delvi,deln 
        IF (RESP.EQ.'y') GOTO 18
        !Modifica INTRCN
        CALL GENERA_LINEA_INTRCN (IPAREQ,IREC)
        !AGREGA EL NUEVO GRUPO A INTRCN
        read (13,205,rec=irec)    MGC,(PINTER(F),F=1,NG)
        WRITE (13,205,REC=IREC)   MG,(PINTER(F),F=1,NG)
        !AGREGA PARÁMETROS DE INTERACCIÓN  
  16    WRITE (6,106)
        READ  (5,107) RESP2
        IF (RESP2.NE.'y'.AND.RESP2.NE.'n') GOTO 16
        IF (RESP2.EQ.'y') CALL CARGA_PARAMETROS_INTRCN (IPAREQ)
      else
	  write (14,201,rec=irec3+1)i1,mg,fs,item,item,caract,mm,
     *				          mj,mk,mi,mh,rpar,qpar,
     *						  pmg,deltc,delpc,delv,delvi,deln
      endif
!-----Checkea si el grupo existe en las otras bases y si no existe lo agrega
  18	lgroup=.false.
	i=0
	do
		i=i+1
		irec5=i+(4-1)*150
		read (14,202,rec=irec5)mgc,fsc
		if (mgc.eq.'fin') exit 
		if (resp1.eq.'y')then
		    if (mgc.eq.mg) then 
		        lgroup=.true. 
		        exit
		    endif
		else
		    if (fsc.eq.fs) then 
		        lgroup=.true. 
		        exit
		    endif
		endif
	enddo
	if (lgroup)goto 15
!---------Genera una nueva linea donde se escribira el nuevo grupo en la última tabla
      CALL GENERA_LINEA_GRUPOSRAM (MG,4,IREC4,resp1)
      if(resp1.eq.'y')then
        write (14,201,rec=irec4)0,mg,fs,item,item,caract,mm,
     *			              mj,mk,mi,mh,rpar,qpar,
     *						  pmg,deltc,delpc,delv,delvi,deln
      else
	  write (14,201,rec=irec4+1)0,mg,fs,item,item,caract,mm,
     *			              mj,mk,mi,mh,rpar,qpar,
     *						  pmg,deltc,delpc,delv,delvi,deln
      endif
  15	WRITE(6,*) '***The new group was added successfully!***'

C-----FORMATOS
   3  FORMAT (1x,/,60x,'> ',$)
   8	format (' Give the caract to use the group ',a8,
     +		/,'  only in the synthesis of:',
     +		/,'  Aromatic solvents:',42x,'AR1',
     +		/,'  Single substance groups:',36x,'GR1',
     +		/,'  Aliphatic and Mixed aromatic-aliphatic solvents: '
     +		 ,11x,'DV1',/,'  Aliphatic, Mixed aromatic-aliphatic and',
     +		  ' Cyclic solvents:',4x,'DV2',
     +		/,' or as a terminal group:',
     +		/,'  when there is not any dual valence group in its ',
     +		  'maingroup:',2X,'SV2',/,'  when there are any dual valence ',
     +		  'group in its maingroup:',5x,'SV1',/,
     +		  ' If the group will not be taken into account:',19x,'0')
 101	FORMAT (/,' Option not available in the actual version ')
 102  FORMAT(/,' Elija el grupo principal al que pertenece', 
     *          ' el subgrupo:')
 103  FORMAT (a3)
 104  format (i2)
 105  FORMAT (i4)
 106  FORMAT (/,'The new group was added successfully to GRUPOSRAM.MDS 
     +and INTRCN.MDS',//,'Do you want to enter some interaction 
     +parameter for the new group? (y/n)',//,52x,'>',$)
 107  FORMAT (a1)
 200  FORMAT (a8)
 201	FORMAT (i4,2a8,2i4,a3,5i2,8d15.8)
 202  FORMAT (4x,2a8)
 203  FORMAT (i4,a8)
 204  FORMAT (4x,a8)
 205  FORMAT (a8,70d12.5)
      RETURN
      END
C
      SUBROUTINE GENERA_LINEA_INTRCN (IPAREQ,IREC)
      IMPLICIT real*8 (A-H,O-Z)
      PARAMETER (NG=70)
      CHARACTER*8 MGC
      DIMENSION PINTER(NG)
      
      ngrmax = max_sub_int (ipareq)+1
      irec=NGRMAX+(ipareq-1)*70
      read (13,10,rec=irec)    MGC,(PINTER(F),F=1,NG)
      WRITE (13,10,rec=irec+1) MGC,(PINTER(F),F=1,NG)
!-----FORMATOS
  10  FORMAT (a8,70d12.5)
      RETURN
      END
C
      SUBROUTINE GENERA_LINEA_GRUPOSRAM (MG,IPAREQ,IREC4,resp)
      use SubGrupos
      IMPLICIT real*8 (A-H,O-Z)
      character*8 mgc,fsc,mg
      character*3 caractc
      character*1 resp
      lin = Size_LSubGroups()+1
      if (resp.eq.'y') then
        irec4=lin+(ipareq-1)*150
		read (14,201,rec=irec4)   ni,mgc,fsc,item,item,caractc,mmc,
     *				              mjc,mkc,mic,mhc,rparc,qparc,
     *							  pmgc,deltcc,delpcc,delvc,delvic,delnc
        write (14,201,rec=irec4+1)ni,mgc,fsc,item,item,caractc,mmc,
     *				              mjc,mkc,mic,mhc,rparc,qparc,
     *							  pmgc,deltcc,delpcc,delvc,delvic,delnc
      else
	  do
		    irec4=lin+(ipareq-1)*150
		    read (14,201,rec=irec4)   ni,mgc,fsc,item,item,caractc,mmc,
     *				              mjc,mkc,mic,mhc,rparc,qparc,
     *							  pmgc,deltcc,delpcc,delvc,delvic,delnc
		    if(mg.eq.mgc)exit
		    write (14,201,rec=irec4+1)ni,mgc,fsc,item,item,caractc,mmc,
     *				              mjc,mkc,mic,mhc,rparc,qparc,
     *							  pmgc,deltcc,delpcc,delvc,delvic,delnc
     		    lin=lin-1		
	  enddo
	endif
201	FORMAT (i4,2a8,2i4,a3,5i2,8d15.8)
	return
	end
C
      SUBROUTINE CARGA_PARAMETROS_INTRCN (IPAREQ)
      USE CONSTANTES
      use SubGrupos
      IMPLICIT real*8 (A-H,O-Z)
      CHARACTER*8 MGV(NMG)
      INTEGER GRUPOS(NMG)
      
  11  WRITE(6,25)
      READ (5,*) NP
      ngrmax = max_sub_int (ipareq)
      DO I=1,NGRMAX
        irec =I+(ipareq-1)*70
		grupos(i) = i
		READ(13,24,rec=irec)MGV(I)
      ENDDO   
      write (6,26) 
      DO I=1,NP
        CALL Groups_Present ()
  12    WRITE (6,23)
        WRITE (6,21) 
        READ  (5,*) N
        IF (N.GE.NGRMAX) THEN
            WRITE (6,22) ngrmax
            GOTO 12
        ENDIF
        CALL MODIFICAR_INTRCN(NGRMAX,N,IPAREQ)
      ENDDO  
!-----FORMATOS
  21  FORMAT (/,' n = ',$)
  22  FORMAT (/,'The value of N should be less than = ',I2)
  23  FORMAT (//,'Select parameter a(mn)')
  24  FORMAT (a8)    
  25  FORMAT (//,'How many parameters would you like 
     +add to this group? ',$)
  26  FORMAT (////)
      RETURN      
      END
C
      SUBROUTINE AJUSTAR_NOMBRE (NC,NB)
C
      CHARACTER*8 NC,NB
      
      	nc = adjustr (nc)
		nb = ')'
		nb = adjustr (nb)
		i=8
		do while (nc(i:i).ne.' ')
		    j=i-1
		    nb(j:j)=nc(i:i)
		    i=i-1
		enddo
		j=i-1
		nb(j:j)='('
	return
	end
