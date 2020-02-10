! ***************************************************************
! ***************************************************************
! * fichier include pour les sorties de modele en netcdf        *
! *								*
! * auteur         :    jfleroux                                *
! * org            :    IFREMER                                 *
! * date creation  :    11/01/02                                *
! * derniere modif :    06/06/02 par jfleroux                   *
! *-------------------------------------------------------------*
! * Modifications :     					*
! *								*
! ***************************************************************
! ***************************************************************

        implicit none
	
! MODIF JFLR 06/06/02
! dimensions du tableau temporaire de travail

	integer		ionc_dimibid
	parameter 	(ionc_dimibid=1200)
	integer		ionc_dimjbid
	parameter 	(ionc_dimjbid=600)
	integer		ionc_dimkbid
	parameter 	(ionc_dimkbid=40)
! FIN MODIF JFLR 06/06/02

	real*4		ionc_valmanque
	parameter 	(ionc_valmanque=1.e+37)
	integer 	ionc_longtabfich
	parameter 	(ionc_longtabfich=1000)
	integer 	ionc_longnomfich
	parameter 	(ionc_longnomfich=200)
	character*(*)	ionc_originet
	parameter 	(ionc_originet='01-JAN-0000 00:00:00')
	character*(*)	ionc_longnamet
	parameter 	(ionc_longnamet='Secondes depuis la date origine')
	character*(*)	ionc_unitst
	parameter 	(ionc_unitst='seconds')

	integer 	ionc_limi,ionc_lima,ionc_pasi
	integer 	ionc_ljmi,ionc_ljma,ionc_pasj
	integer 	ionc_lkmi,ionc_lkma,ionc_pask
	integer 	ionc_nfich
	character*(ionc_longnomfich) ionc_nomfich(ionc_longtabfich)
	character*(30)  ionc_rout
	integer		ionc_idfich(ionc_longtabfich)
! MODIF JFLR 06/06/02
	real*4		ionc_bid(ionc_dimibid * ionc_dimjbid * ionc_dimkbid)
	common /ionc_var/ ionc_bid
! FIN MODIF JFLR 06/06/02
	common /ionc_var/ ionc_nomfich,ionc_idfich,ionc_nfich,ionc_rout
	common /ionc_mod/ ionc_limi,ionc_lima,ionc_ljmi,ionc_ljma
        common /ionc_mod/ ionc_lkmi,ionc_lkma,ionc_pasi
        common /ionc_mod/ ionc_pasj,ionc_pask

!
! codes d'erreur:
!
      integer ionc_errcreation
      integer ionc_errfich
      integer ionc_err_i
      integer ionc_err_j
      integer ionc_err_k
      integer ionc_err_t

      parameter (ionc_errcreation = -1000)
      parameter (ionc_errfich     = -1001)
      parameter (ionc_err_i       = -1002)
      parameter (ionc_err_j       = -1003)
      parameter (ionc_err_k       = -1004)
      parameter (ionc_err_t       = -1005)
