program ECOFAS
!
    use CONSTANTES
    use input_data,only:CHECK_DATABASE
    use Evaluation, only:Initialization
    use Input
    
    implicit none
    
!...Integers
	integer I,nfunc
!...real*8
!...Characters
    character*1 ret
	character*6 nomb
!...Logical
    logical::NOACE,oprop
!...Commons
    common/EVAL/INV	
    logical INV

!...Sentencias
    write (6,"(2X,'****************************************************************'&
                    ,//,27X,'MolDeS 2.1',//,13x,&
                ' Version prepared by J. A. Scilipoti',//,13X,&
                ' December, 2021.',//,13X,&
                ' IPQA (UNC-CONICET), 5016 Cordoba, Argentina.',//,2X,&
                '****************************************************************',////)")
    call PAUSA 
! segunda prueba git
    call initialization()
    call open_proppdb(oprop)
    if(oprop)stop
    nfunc=1
    do while (.true.)    
        write (6,600)
        read (5,'(i2)') nfunc
	    if (nfunc.eq.1) then !Run a molecular design problem
            call Moldes(.True.)
        elseif (nfunc.eq.2) then !Evaluate specific solvents performances
		    call Moldes(.False.)
	    elseif (nfunc.eq.3) then !n-octanol/water partition coefficients
            call log_Pow()
        elseif (nfunc.eq.4) then !distribution coefficients in other systems
3200	    write (6,'(/," Give output file name (maximum 6 characters)")')
            read (5,'(a)',err=3200) nomb
            open (unit=1,file=nomb//'.pgr',form='formatted')
	        write (1,*) ' ECOFAC 1.0 '
	        write (1,*) ' Case: '
            write (1,'(6x,a6)') nomb
	        call distcoef
	        write(6,*) ' Find results in the output file .pgr'
	        close (unit=1)
        elseif (nfunc.eq.5) then !solubilities in water
	        write(6,*)' Please use the option 6'
            call PAUSA 
	    elseif (nfunc.eq.6) then ! other solubilities
	        call solubilities ()
        elseif (nfunc.eq.7) then ! activity coefficients
            call activity_coefficients
	    elseif (nfunc.eq.8) then ! Check databases
            CALL CHECK_DATABASE
	    elseif (nfunc.eq.9) then !View, modify or add components in databases
	        !open (unit=45,file='Getting_Started.pdf')
	        call datos
        elseif(NFUNC.EQ.10)then !Property prediction
            call properties1
        elseif(nfunc==11)then ! Design of solvents for chemical reactions
            call solvrea()
        elseif(nfunc==0)then
            exit
	    endif
    enddo

    close(unit=7)
!-----FORMATOS
 30   format (1x,/,60x,'> ',$)
 981  FORMAT (I1)
 600  format (////' ','Select the function to perform with ECOFAC :',/&
             //,11x,' Run a molecular design problem             : 1',&
             //,11x,' Evaluate specific solvents performances    : 2',&
             //,11x,' Estimate',&
              /,11x,'     n-octanol/water partition coefficients : 3',&
              /,11x,' distribution coefficients in other systems : 4',&
              /,11x,'	                 solubilities in water : 5',  &
              /,11x,'	                    other solubilities : 6',  &
              /,11x,'                      activity coefficients : 7',&
             //,11x,' Check databases                            : 8',&
             //,11x,' View, modify or add components in databases: 9',&
             //,11x,' Property prediction                        :10',&
             //,11x,' Design of solvents for chemical reactions  :11',&             
             //,11x,' Exit                                       : 0',&
             //,60x,'> ',$)

endprogram ECOFAS





!
!FORMATO PARA ESCRIBIR SUBRUTINAS
!
!    subroutine name
!!-------------------------------------------------
!!   Descripci�n
!!   - Variables de entrada
!!       vars:
!!   - Variables de salida
!!       vars:
!!-------------------------------------------------
!    parameter()
!    implicit none
!!...Integers
!
!!...real*8
!
!!...Characters
!
!!...Logical
!
!!...Commons
!
!!...Sentencias
!
!    endsubroutine
