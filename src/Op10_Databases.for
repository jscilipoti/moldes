	subroutine DATOS
C 
	INTEGER NO, OB,PARINT,BD
	LOGICAL ASOC
	COMMON/OPC/NO,OB,PARINT,BD
	COMMON/NOP1/K1
	COMMON/NOP2/K10
	COMMON/NOP3/K2
	COMMON/AS/ASOC
C
  33	Continue
	CALL OPCION
	CALL AB_BAN
	CALL PRESEN
	CALL PAUSA
C
	IF(NO.EQ.1) THEN
		NR=0
	  CALL VER_COM(NR)
	  IF(K1.EQ.2) GOTO 33
	END IF
	IF(NO.EQ.2) THEN
	  CALL MOD_COM
	  IF(K10.EQ.3) GOTO 33
	END IF
	IF(NO.EQ.3) THEN
	  CALL IN_COM
	  IF(K2.EQ.2) GOTO 33
	END IF
C
	CALL CI_BAN1

	GOTO 500
 	WRITE(6,*)'  '
	WRITE(6,*)'OK'
C

	call pausa
 500	return
	END
C
C
C
	SUBROUTINE PRESEN
C
	INTEGER NO, OB,PARINT,BD
	COMMON/OPC/NO,OB,PARINT,BD
C
!	WRITE(6,2000)
!	WRITE(6,2100)
!	WRITE(6,2200)'PLANTA PILOTO DE INGENIERIA QUIMICA'
!2200	FORMAT(' ','#',2X,A35,40X,'#')
!	WRITE(6,2300)'PROP.PDB','FECHA : 12/91'
!2300	FORMAT(' ','#',2X,A8,37x,A13,17X,'#')
!	WRITE(6,2400)'VERSION : 1'
!2400	FORMAT(' ','#',47X,A11,19X,'#')
!	WRITE(6,2500)'PROGRAMACION : A. Mengarelli'
!	WRITE(6,2500)'               E. Pretel    '
!2500	FORMAT(' ','#',47X,A28,2X,'#')
!	WRITE(6,2100)
      IF (OB==1) THEN
	  WRITE(6,2000)
	  WRITE(6,2100)
	  WRITE(6,2100)
	  WRITE(6,2600)'   PROP.PDB is a physical properties       '
	  WRITE(6,2600)'      of pure components database          '
      	WRITE(6,2100)
	  WRITE(6,2100)
	  WRITE(6,2700)'REFERENCIA : Physical and Thermodynamic Properties '
     *               ,'of Pure Chemicals'
        WRITE(6,2800)'Data Compilation'
	  WRITE(6,2800)'   T.E.Daubert  '
	  WRITE(6,2800)'   R.P.Danner   '  
      	WRITE(6,2100)
	  WRITE(6,2000)
	ELSE IF (OB==2) THEN
		WRITE(6,2000)
	  WRITE(6,2100)
	  WRITE(6,2100)
	  IF(PARINT.EQ.1)THEN
	      WRITE(6,2600)'      INTRCN.MDS is a UNIFAC group         '
	      WRITE(6,2600)'     interaction parameters database       '
	  ELSEIF(PARINT.EQ.2)THEN
	      WRITE(6,2600)'     INTRCNAS.MDS is a A-UNIFAC group      '
	      WRITE(6,2600)'     interaction parameters database       '
	  ELSEIF(PARINT.EQ.3)THEN
	      WRITE(6,2600)'     PARENEAS.MDS is a A-UNIFAC energy     '
	      WRITE(6,2600)'     of association parameters database    '
	  ELSE
	      WRITE(6,2600)'     PARVOLAS.MDS is a A-UNIFAC volume     '
	      WRITE(6,2600)'     of association parameters database    '
	  ENDIF	  
      	WRITE(6,2100)
      	WRITE(6,2100)
	  WRITE(6,2000)
	ELSE IF (OB==3) THEN
	  WRITE(6,2000)
	  WRITE(6,2100)
	  WRITE(6,2100)
	  WRITE(6,2600)'    GRUPOSRAM.MDS is a properties and      '
	  WRITE(6,2600)'characterizations of UNIFAC groups database'
      	WRITE(6,2100)
      	WRITE(6,2100)
	  WRITE(6,2000)
	ENDIF
C-----FORMATOS
2000	FORMAT(80('#'))
2100	FORMAT('#',78X,'#')
2600	FORMAT('#',17X,A43,18X,'#')
2700	FORMAT(' ','#',4X,A51,A17,5X,'#')
2800	FORMAT(' ','#',30X,A16,31X,'#')
C
	RETURN
	END
C
C
C
	SUBROUTINE OPCION
C
	INTEGER NO, OB, PARINT,BD
	COMMON/OPC/NO,OB,PARINT,BD
	COMMON/AS/ASOC
	LOGICAL ASOC
c
C.....Selección del tipo de base de datos
      WRITE (6,1000)
  10  WRITE (6,1010)
	READ (5,2000,ERR=10) OB
	IF (OB<1.OR.OB>3) GOTO 10
C
C.....Selección del tipo de base de datos de interacción
c.....si OB=2
      IF (OB.EQ.2)THEN
        WRITE (6,1001)
  30    WRITE (6,1010)
  	  READ (5,2000,ERR=30) PARINT
	  IF (PARINT<1.OR.PARINT>4) GOTO 30
	  IF(PARINT.EQ.1)THEN
	      ASOC=.FALSE.
	  ELSE
	      ASOC=.TRUE.
	  ENDIF
	  IF(PARINT.EQ.1.OR.PARINT.EQ.2)BD=13
	  IF(PARINT.EQ.3)BD=15
	  IF(PARINT.EQ.4)BD=16
	ENDIF
C.....Selección de la acción a realizar sobre la base de datos
      WRITE(6,1020)
      IF (OB.EQ.1) THEN
        WRITE(6,1100)
	ELSE IF (OB.EQ.2) THEN
	  WRITE(6,1200)
	ELSE IF (OB.EQ.3) THEN
	  WRITE(6,1300)
	ENDIF
C
  20  WRITE(6,1010)
	READ(5,2000,ERR=20) NO
	IF (NO<1.OR.NO>3) GOTO 20
C
C-----FORMATOS
1000  FORMAT(///,' Select the database:',///,
     +      6X,'1 - Prop.pdb',//,
     +      6X,'2 - Databases interaction parameters',//,
     +      6X,'3 - Gruposram.mds',//)
1001  FORMAT(///,' Select the database interaction parameters:',///,
     +      6X,'1 - Intrcn.mds',//,
     +      6X,'2 - Intrcnas.mds',//,
     +      6X,'3 - Pareneas.mds',//,
     +      6X,'4 - Parvolas.mds',//)
1010  FORMAT(/,60X,'> ',$)
1020  FORMAT (/,' Choose how to work on the database:',//)
1100  FORMAT(6X,'1 - View the properties for a component',//,
     +       6X,'2 - Modify and/or add properties for a component',//,
     +       6X,'3 - Enter a new component to the database',/)
1200  FORMAT(6X,'1 - Check parameters',//,
     +       6X,'2 - Modify and/or add parameters',//,
     +       6X,'3 - Enter a new group',/)
1300  FORMAT(6X,'1 - View the properties for a subgroup',//,
     +       6X,'2 - Modify and/or add properties for a subgroup',//,
     +       6X,'3 - Enter a new subgroup to the database',/)    
2000  FORMAT(I1)
	RETURN
	END 
