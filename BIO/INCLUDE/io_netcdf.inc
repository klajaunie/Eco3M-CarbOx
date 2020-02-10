c ***************************************************************
c ***************************************************************
c * fichier include pour les sorties de modele en netcdf        *
c *                                                             *
c * auteur         :    jfleroux                                *
c * org            :    IFREMER                                 *
c * date creation  :    21/02/01                                *
c * derniere modif :    09/12/03 par jfleroux                   *
c *-------------------------------------------------------------*
c * Modifications :                                             *
c *     09/01/02  : var SIG -> Z, ajout routine ionc_ouvre      *
c *                 ajout ionc_ecrit_lon, ionc_ecrit_lat        *
c *                 ajout ionc_lit_lon, ionc_lit_lat            *
c *     11/01/02  :                                             *
c *       att units de lat degree_north -> degrees_north	*
c *       att units de lon degree_east -> degrees_east	        *
c *       ajout fonction ionc_flimax ~ ionc_limax               *
c *       ajout fonction ionc_fljmax ~ ionc_ljmax               *
c *       ajout fonction ionc_flkmax ~ ionc_lkmax               *
c *       ajout fonction ionc_fltmax ~ ionc_ltmax               *
c *       ajout routine ionc_creer8                             *
c *       ajout routine ionc_ecritr8                            *
c *       ajout routine ionc_litr8                              *
c *	14/01/02  : correction bug ionc_liti4                   *
c *	24/01/02  : correction bug ionc_litr4_xy                *
c *       correction bug ionc_litr4_xyt                         *
c *       ajout routine ionc_creer8_t                           *
c *       ajout routine ionc_ecritr8_t                          *
c *       ajout routine ionc_litr8_t                            *
c *       ajout routine ionc_creer4_xt                          *
c *       ajout routine ionc_ecritr4_xt                         *
c *       ajout routine ionc_litr4_xt                           *
c *	15/02/02  : correction ionc_litr4_zxy                   *
c *	18/02/02  : correction ionc_creer4_zxyt                 *
c *       correction ionc_litr4_zxyt                            *
c *	28/02/02  : ajout ionc_litr4_xt                         *
c *       ajout fonction ionc_fnlimax  	                        *
c *       ajout fonction ionc_fnljmax  	                        *
c *       ajout fonction ionc_fnlkmax  	                        *
c *       ajout fonction ionc_fnltmax 	                        *
c *	05/03/02  : ajout ionc_litr4_zxt                        *
c *	15/03/02  : ajout ionc_litr4_subxy                      *
c *	  ajout ionc_litr4_subxyt                               *
c *	   04/02  : ajout ionc_creer_sta                        *
c *	            ajout ionc_creer4_sta                       *
c *	            ajout ionc_ecritr4_sta                      *
c *	            ajout ionc_creer4_stat                      *
c *	            ajout ionc_ecritr4_stat                     *
c *	            ajout ionc_litr4_stat                       *
c *	            ajout ionc_creer4_zstat                     *
c *	            ajout ionc_ecritr4_zstat                    *
c *	            ajout ionc_litr4_zstat                      *
c *	16/05/02  : ajout ionc_fnlsta                           *
c *	04/07/02  : correction   ionc_litr4_subxy               *
c *	  correction ionc_litr4_subxyt                          *
c *	27/08/02  : ajout ionc_creer4_z                         *
c *	            ajout ionc_ecritr4_z                        *
c *	            ajout ionc_putr4_z                          *
c *	            ajout ionc_litr4_z                          *
c *	            ajout ionc_ecritr4_sig                      *
c *	            ajout ionc_putr4_sig                        *
c *	24/10/02  : ajout ionc_littorigine                      *
c *	25/10/02  : ajout ionc_litdimi                          *
c *	            ajout ionc_litdimj                          *
c *	26/10/02  : ajout ionc_showvar                          *
c *	13/03/03  : variables lat,lon en double precision	*
c *	25/03/03  : lecture et ecriture des lat, lon avec       *
c *                 des variables en real*4              	*
c *		    ajout ionc_litr4_lon, ionc_litr4_lat        *
c *		    ajout ionc_ecritr4_lon, ionc_ecritr4_lat	*
c *	18/04/03  : creation avec variables lat,lon real        *
c *		    ajout ionc_creerr4,ionc_creerr4_sta         *
c *	25/08/03  :                                             *
c *                 ajout routine ionc_creer8_xy                *
c *                 ajout routine ionc_ecritr8_xy               *
c *                 ajout routine ionc_putr8_xy                 *
c *                 ajout routine ionc_litr8_xy                 *
c *	26/08/03  :                                             *
c *                 ajout routine ionc_creer8_zxyt              *
c *                 ajout routine ionc_ecritr8_zxyt             *
c *                 ajout routine ionc_putr8_zxyt               *
c *                 ajout routine ionc_litr8_zxyt               *
c *                 ajout routine ionc_creer8_zxy               *
c *                 ajout routine ionc_ecritr8_zxy              *
c *                 ajout routine ionc_putr8_zxy                *
c *                 ajout routine ionc_litr8_zxy                *
c *                 ajout routine ionc_creer8_xyt               *
c *                 ajout routine ionc_ecritr8_xyt              *
c *                 ajout routine ionc_putr8_xyt                *
c *                 ajout routine ionc_litr8_xyt                *
c *	09/12/03  : ajout ionc_litr4_subzxyt                    *
c *	10/12/03  : ajout ionc_litdimk                          *
c *                 ajout ionc_litdiml                          *
c *	12/12/03  : ajout ionc_litr4_sublxyt                    *
c *	16/12/03  : ajout ionc_litr4_sublzxyt                   *
c *                 ajout ionc_ecritr4_sublxyt                  *
c *                 ajout ionc_ecritr4_sublzxyt                 *
c ***************************************************************
c ***************************************************************
c *     valeur manquante par defaut : 1.e+37                    *
c *     nombre de fichiers          : 10                        *
c *     taille de nom du fichier    : 200                       *
c *     origine des temps par defaut: '01-JAN-0000 00:00:00'    *
c *     nom long du temps par defaut:                           *
c *                          'Secondes depuis la date origine'  *
c *     unite de temps par defaut   :	'seconds'               *
c ***************************************************************

c ***************************************************************
c * subroutine ionc_init             				*
c * Role : initialisation du module ionc 			*
c *								*
c ***************************************************************
c      ionc_init()

c ***************************************************************
c * SUBROUTINE ionc_creer             				*
c *								*
c * Role : creer le fichier Netcdf qui contiendra les donnees	*
c *	   Initialiser les donnees ne dependant pas du time :	*
c *			- lon,lat				*
c *			- SIG(k)				*
c *	   initialise les variables correspondants à la loupe	*
c *	   initialise les unites de temps et l'origine          *
c *        aux valeurs par defaut                               *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier : chemin complet du fichier a 	*
c *				creer				*
c *		- dg   :   	   	   	   	   	*
c *		- g0   :   	   	   	   	   	*
c *		- dfi  :   	   	   	   	   	*
c *		- fi0  :   	   	   	   	   	*
c *		- limi : limite min de l'indice i de la loupe   *
c *		- lima : limite max de l'indice i de la loupe   *
c *		- pasi : pas de l'indice i de la loupe          *     
c *		- ljmi : limite min de l'indice j de la loupe   *
c *		- ljma : limite max de l'indice j de la loupe	*
c *		- pasj : pas de l'indice j de la loupe          *
c *		- lkmi : limite min de l'indice k de la loupe   *
c *		- lkma : limite max de l'indice k de la loupe   *
c *		- pask : pas de l'indice k de la loupe          *
c *		- kmax : borne max pour k                       *
c *		- dim_time : nombre de pas de time que le 	*
c *			     fichier contiendra si dim_time <=0 *
c *			     la dimension time est unlimited    *
c *			     si dim_time = -1 pas de dim temps  *
c *		- SIG : niveaux					*
c ***************************************************************
c      ionc_creer(nom_fichier,dg,g0,dfi,fi0,
c                 limi,lima,pasi,ljmi,ljma,pasj,
c                 lkmi,lkma,pask,kmax,dim_time,sig)
c ******** parametres de la subroutine *******
c      character*(*)   nom_fichier
c      double	       dg,g0,dfi,fi0
c      integer         limi,lima,pasi,ljmi,ljma,pasj,lkmi,lkma,pask,kmax
c      integer         dim_time
c      real	       sig(0:kmax+1)

c ***************************************************************
c * SUBROUTINE ionc_creerr4           				*
c *								*
c * Role : creer le fichier Netcdf qui contiendra les donnees	*
c *	   Initialiser les donnees ne dependant pas du time :	*
c *			- lon,lat				*
c *			- SIG(k)				*
c *	   initialise les variables correspondants à la loupe	*
c *	   initialise les unites de temps et l'origine          *
c *        aux valeurs par defaut                               *
c *  les variables lat et lon sont ecrites en real		*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier : chemin complet du fichier a 	*
c *				creer				*
c *		- dg   :   	   	   	   	   	*
c *		- g0   :   	   	   	   	   	*
c *		- dfi  :   	   	   	   	   	*
c *		- fi0  :   	   	   	   	   	*
c *		- limi : limite min de l'indice i de la loupe   *
c *		- lima : limite max de l'indice i de la loupe   *
c *		- pasi : pas de l'indice i de la loupe          *     
c *		- ljmi : limite min de l'indice j de la loupe   *
c *		- ljma : limite max de l'indice j de la loupe	*
c *		- pasj : pas de l'indice j de la loupe          *
c *		- lkmi : limite min de l'indice k de la loupe   *
c *		- lkma : limite max de l'indice k de la loupe   *
c *		- pask : pas de l'indice k de la loupe          *
c *		- kmax : borne max pour k                       *
c *		- dim_time : nombre de pas de time que le 	*
c *			     fichier contiendra si dim_time <=0 *
c *			     la dimension time est unlimited    *
c *			     si dim_time = -1 pas de dim temps  *
c *		- SIG : niveaux					*
c ***************************************************************
c      ionc_creerr4(nom_fichier,dg,g0,dfi,fi0,
c                 limi,lima,pasi,ljmi,ljma,pasj,
c                 lkmi,lkma,pask,kmax,dim_time,sig)
c ******** parametres de la subroutine *******
c      character*(*)   nom_fichier
c      real	       dg,g0,dfi,fi0
c      integer         limi,lima,pasi,ljmi,ljma,pasj,lkmi,lkma,pask,kmax
c      integer         dim_time
c      real	       sig(0:kmax+1)

c ***************************************************************
c * SUBROUTINE ionc_creer_sta             			*
c *								*
c * Role : creer le fichier Netcdf qui contiendra les donnees   *
c *        non maillees (stationnelles)                         *
c *	   Initialiser les donnees ne dependant pas du time :	*
c *			- sta				*
c *			- SIG(k)				*
c *	   initialise les unites de temps et l'origine          *
c *        aux valeurs par defaut                               *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier : chemin complet du fichier a 	*
c *				creer				*
c *		- lkmi : limite min de l'indice k de la loupe   *
c *		- lkma : limite max de l'indice k de la loupe   *
c *		- pask : pas de l'indice k de la loupe          *
c *		- kmax : borne max pour k                       *
c *			0 si 2D					*
c *		- dim_sta : nombre de stations                  *
c *		- dim_time : nombre de pas de time que le 	*
c *			     fichier contiendra si dim_time =0  *
c *			     la dimension time est unlimited    *
c *			     si dim_time = -1 pas de dim temps  *
c *		- SIG : niveaux					*
c ***************************************************************
c       ionc_creer_sta(nom_fichier,lkmi,lkma,pask,kmax,
c                      dim_sta,dim_time,sig,lat,lon)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)	nom_fichier
c	integer 	lkmi,lkma,pask
c	integer 	dim_time
c	integer 	dim_sta
c	integer 	kmax
c        real 		sig(0:kmax+1)
c	double          lat(dim_sta)
c	double          lon(dim_sta)

c ***************************************************************
c * SUBROUTINE ionc_creerr4_sta             			*
c *								*
c * Role : creer le fichier Netcdf qui contiendra les donnees   *
c *        non maillees (stationnelles)                         *
c *	   Initialiser les donnees ne dependant pas du time :	*
c *			- sta				*
c *			- SIG(k)				*
c *	   initialise les unites de temps et l'origine          *
c *        aux valeurs par defaut                               *
c *  les variables lat et lon sont ecrites en real		*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier : chemin complet du fichier a 	*
c *				creer				*
c *		- lkmi : limite min de l'indice k de la loupe   *
c *		- lkma : limite max de l'indice k de la loupe   *
c *		- pask : pas de l'indice k de la loupe          *
c *		- kmax : borne max pour k                       *
c *			0 si 2D					*
c *		- dim_sta : nombre de stations                  *
c *		- dim_time : nombre de pas de time que le 	*
c *			     fichier contiendra si dim_time =0  *
c *			     la dimension time est unlimited    *
c *			     si dim_time = -1 pas de dim temps  *
c *		- SIG : niveaux					*
c ***************************************************************
c       ionc_creerr4_sta(nom_fichier,lkmi,lkma,pask,kmax,
c                      dim_sta,dim_time,sig,lat,lon)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)	nom_fichier
c	integer 	lkmi,lkma,pask
c	integer 	dim_time
c	integer 	dim_sta
c	integer 	kmax
c       real 		sig(0:kmax+1)
c	real            lat(dim_sta)
c	real            lon(dim_sta)

c ***************************************************************
c * subroutine ionc_titre             			        *
c *								*
c * Role : attribut titre d'un fichier netcdf                   *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- titre       :  titre du fichier	        *
c ***************************************************************
c      ionc_titre(nom_fichier,titre)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  titre

c ***************************************************************
c * subroutine ionc_history            			        *
c *								*
c * Role : attribut history d'un fichier netcdf                 *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- history    : historique de creation du fichier*
c ***************************************************************
c      ionc_history(nom_fichier,history)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  history

c ***************************************************************
c * subroutine ionc_tunite            			        *
c *								*
c * Role : modifie l'attribut units de la variable time         *
c *        d'un fichier netcdf                                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- tunite      :  unite de temps                 *
c ***************************************************************
c      ionc_tunite(nom_fichier,tunite)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  tunite

c ***************************************************************
c * subroutine ionc_torigine            			*
c *								*
c * Role : modifie l'attribut timr_origin de la variable time   *
c *        d'un fichier netcdf                                  *
c *	   entrer l'orine sous la forme JJ-MMM-AAAA hh:mm:ss	*
c *	   avec MMM trois premieres lettres du mois en anglais	*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- torigine    :  origine du temps               *
c ***************************************************************
c      ionc_torigine(nom_fichier,torigine)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  torigine

c ***************************************************************
c * subroutine ionc_vatt_range             			*
c *								*
c * Role : ajoute l'attribut valid_range a une variable         * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- valid_min   :  valeur valide min              *
c * 		- valid_max   :  valeur valide max              *
c ***************************************************************
c      ionc_vatt_range(nom_fichier,var_name,valid_min,valid_max)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*) var_name
c      real*4 valid_min, valid_max

c ***************************************************************
c * subroutine ionc_vatt_fformat             			*
c *								*
c * Role : ajoute l'attribut Fortran_format a une variable      * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- fformat     :  format fortran                 *
c ***************************************************************
c      ionc_vatt_fformat(nom_fichier,var_name,fformat)
c ******** PARAMETRES DE LA SUBROUTINE *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      character*(*)  fformat

c ***************************************************************
c * subroutine ionc_vatt_fill             			*
c *								*
c * Role : ajoute l'attribut _FillValue a une variable          * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- fillvalue   :  valeur de remplissage          *
c ***************************************************************
c      ionc_vatt_fill(nom_fichier,var_name,fillvalue)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real           fillvalue

c ***************************************************************
c * subroutine ionc_vatt_missing             			*
c *								*
c * Role : ajoute l'attribut missing_value a une variable       * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- mvalue      :  valeur manquante               *
c ***************************************************************
c      ionc_vatt_missing(nom_fichier,var_name,mvalue)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real           mvalue

c ***************************************************************
c * subroutine ionc_ecrit_lon             			*
c *								*
c * Role : ecrit la variable lon                		*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c ***************************************************************
c      ionc_ecrit_lon(nom_fichier,var,imin,imax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        imin,imax
c       double         var(imin:imax) 

c ***************************************************************
c * subroutine ionc_ecrit_lat             			*
c *								*
c * Role : ecrit la variable lat                		*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var         :  variable	                *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j	        *
c ***************************************************************
c     ionc_ecrit_lat(nom_fichier,var,jmin,jmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        jmin,jmax
c       double         var(jmin:jmax) 

c ***************************************************************
c * subroutine ionc_ecritr4_lon             			*
c *								*
c * Role : ecrit la variable lon a partir d'une variable real   *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c ***************************************************************
c      ionc_ecritr4_lon(nom_fichier,var,imin,imax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        imin,imax
c       real         var(imin:imax) 

c ***************************************************************
c * subroutine ionc_ecritr4_lat             			*
c *								*
c * Role : ecrit la variable lat a partir d'une variable real   *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var         :  variable	                *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j	        *
c ***************************************************************
c     ionc_ecritr4_lat(nom_fichier,var,jmin,jmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        jmin,jmax
c       real         var(jmin:jmax) 

c ***************************************************************
c * subroutine ionc_lit_lon             			*
c *								*
c * Role : lit la variable lon                		        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_lit_lon(nom_fichier,var,imin,imax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        imin,imax
c       double         var(imin:imax) 

c ***************************************************************
c * subroutine ionc_lit_lat             			*
c *								*
c * Role : lit la variable lat                		        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j	        *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_lit_lat(nom_fichier,var,jmin,jmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        jmin,jmax
c       double         var(jmin:jmax) 

c ***************************************************************
c * subroutine ionc_litr4_lon             			*
c *								*
c * Role : lit la variable lon dans une variable real           *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_lon(nom_fichier,var,imin,imax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        imin,imax
c       real           var(imin:imax) 

c ***************************************************************
c * subroutine ionc_litr4_lat             			*
c *								*
c * Role : lit la variable lat dans une variable real           *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j	        *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_lat(nom_fichier,var,jmin,jmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c       integer        jmin,jmax
c       real           var(jmin:jmax) 

c ***************************************************************
c * subroutine ionc_creei4             				*
c *								*
c * Role : creer une variable entiere 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la varible	        *
c ***************************************************************
c      ionc_creei4(nom_fichier,var_name,units,long_name)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name

c ***************************************************************
c * subroutine ionc_ecriti4             			*
c *								*
c * Role : ecrit une variable entiere 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable        	        *
c ***************************************************************
c      ionc_ecriti4(nom_fichier,var_name,var)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        var

c ***************************************************************
c * subroutine ionc_liti4             			        *
c *								*
c * Role : lit une variable entiere dans le fichier 	        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 	sortie:	- var         :  variable        	        *
c ***************************************************************
c      ionc_liti4(nom_fichier,var_name,var)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        var

c ***************************************************************
c * subroutine ionc_creer4             				*
c *								*
c * Role : creer une variable reelle 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c ***************************************************************
c      ionc_creer4(nom_fichier,var_name,units,long_name)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*) var_name,units,long_name

c ***************************************************************
c * subroutine ionc_creer4_sta            			*
c *								*
c * Role : creer une variable reelle 1D = nb stations           * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c    ionc_creer4_sta(nom_fichier,var_name,units,
c                      long_name,scale_factor,add_offset)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*) nom_fichier
c	character*(*) var_name,units,long_name
c        real*4 scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4             			*
c *								*
c * Role : ecrit une variable reelle 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c ***************************************************************
c      ionc_ecritr4(nom_fichier,var_name,var)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real           var

c ***************************************************************
c * subroutine ionc_ecritr4_sta             			*
c *								*
c * Role : ecrit une variable reelle 1D = nb stations           *
c *        a partir d'une variable reelle 1D = nb stations      * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c     ionc_ecritr4_sta(nom_fichier,var_name,var,
c                          imin,imax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax
c        real           var(imin:imax)

c ***************************************************************
c * subroutine ionc_litr4             			*
c *								*
c * Role : ecrit une variable reelle 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 	sortie:	- var         :  variable        	        *
c ***************************************************************
c      ionc_litr4(nom_fichier,var_name,var)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real           var

c ***************************************************************
c * subroutine ionc_creer8             				*
c *								*
c * Role : creer une variable double 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c ***************************************************************
c      ionc_creer8(nom_fichier,var_name,units,long_name)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*) var_name,units,long_name

c ***************************************************************
c * subroutine ionc_ecritr8             			*
c *								*
c * Role : ecrit une variable double 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c ***************************************************************
c      ionc_ecritr8(nom_fichier,var_name,var)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real*8           var

c ***************************************************************
c * subroutine ionc_litr8             			        *
c *								*
c * Role : ecrit une variable double 			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 	sortie:	- var         :  variable        	        *
c ***************************************************************
c      ionc_litr8(nom_fichier,var_name,var)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real*8           var

c ***************************************************************
c * subroutine ionc_ecrit_t             			*
c *								*
c * Role : ecrit le temps        			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c * 		- temps       :  valeur du temps	        *
c ***************************************************************
c      ionc_ecrit_t(nom_fichier,nrec,temps)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      real*8         temps
c      integer        nrec

c ***************************************************************
c * subroutine ionc_lit_t               			*
c *								*
c * Role : lit le temps          			        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- nrec        :  numero d'enregistrement	*
c * 	sortie:	- temps       :  variable temps     	        *
c ***************************************************************
c      ionc_lit_t(nom_fichier,nrec,temps)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      real*8         temps
c      integer        nrec

c ***************************************************************
c * fonction ionc_lit_dimt            	 			*
c *								*
c * Role : lit la dimension du temps        			*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * Resultat   :						*
c *	sortie:	- temps	      :  variable        	        *
c ***************************************************************
c      integer function ionc_lit_dimt(nom_fichier)
c ******** parametres de la fonction *******
c *	character*(*)  nom_fichier

c ***************************************************************
c * subroutine ionc_creer4_t             			*
c *								*
c * Role : creer une variable reelle 1D = temps                 *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer4_t(nom_fichier,var_name,units,long_name,
c                    scale_factor,add_offset)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name
c      real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_t             			*
c *								*
c * Role : ecrit une variable reelle 1D = temps 	        *
c *	   l'indice du temps est fixe par nrec                  *
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      ionc_ecritr4_t(nom_fichier,var_name,var,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real           var
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_t             			*
c *								*
c * Role : ecrit une variable reelle 1D = temps		        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- nrec        :  numero d'enregistrement        *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_t(nom_fichier,var_name,var,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real           var
c      integer        nrec

c ***************************************************************
c * subroutine ionc_creer8_t             			*
c *								*
c * Role : creer une variable double (real*8) 1D = temps        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer8_t(nom_fichier,var_name,units,long_name,
c                    scale_factor,add_offset)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name
c      real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr8_t             			*
c *								*
c * Role : ecrit une variable double (real*8) 1D = temps  	*
c *	   l'indice du temps est fixe par nrec                  *
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      ionc_ecritr8_t(nom_fichier,var_name,var,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real*8           var
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr8_t             			*
c *								*
c * Role : lit une variable double (real*8) 1D = temps		*
c *	   l'indice du temps est fixe par nrec                  *
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- nrec        :  numero d'enregistrement        *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr8_t(nom_fichier,var_name,var,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      real*8           var
c      integer        nrec

c ***************************************************************
c * subroutine ionc_creer4_xt             			*
c *								*
c * Role : creer une variable reelle 3D = lon, temps            * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      subroutine ionc_creer4_xt(nom_fichier,var_name,units,long_name,
c                            scale_factor,add_offset)
c
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*) var_name,units,long_name
c       real*4 scale_factor,add_offset


c ***************************************************************
c * subroutine ionc_ecritr4_xt             			*
c *								*
c * Role : ecrit une variable reelle 3D = lon,  temps           *
c *        a partir d'une variable reelle 2D = lon              * 
c *        et de l'indice du temps                              * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      subroutine ionc_ecritr4_xt(nom_fichier,var_name,var,
c                                 imin,imax,nrec)
c
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c       integer        imin,imax,jmin,jmax
c       real           var(imin:imax,jmin:jmax)
c       integer        nrec


c ***************************************************************
c * subroutine ionc_creer4_xy             			*
c *								*
c * Role : creer une variable reelle 2D = lon, lat 	        * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer4_xy(nom_fichier,var_name,units,long_name,
c                     scale_factor,add_offset)
c ******** parametres de la subroutine *******
c	character*(*)  nom_fichier
c	character*(*)  var_name,units,long_name
c       real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_xy             			*
c *								*
c * Role : ecrit une variable reelle 2D = lon,lat 	        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c ***************************************************************
c      ionc_ecritr4_xy(nom_fichier,var_name,var,
c                      imin,imax,jmin,jmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_litr4_xy             			*
c *								*
c * Role : lit une variable reelle 2D = lon,lat 	        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_xy(nom_fichier,var_name,var,
c                    imin,imax,jmin,jmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_creer8_xy             			*
c *								*
c * Role : creer une variable double 2D = lon, lat 	        * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer8_xy(nom_fichier,var_name,units,long_name,
c                     scale_factor,add_offset)
c ******** parametres de la subroutine *******
c	character*(*)  nom_fichier
c	character*(*)  var_name,units,long_name
c       real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr8_xy             			*
c *								*
c * Role : ecrit une variable double 2D = lon,lat 	        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c ***************************************************************
c      ionc_ecritr8_xy(nom_fichier,var_name,var,
c                      imin,imax,jmin,jmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real*8         var(imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_litr8_xy             			*
c *								*
c * Role : lit une variable double 2D = lon,lat 	        *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr8_xy(nom_fichier,var_name,var,
c                    imin,imax,jmin,jmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real*8         var(imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_litr4_subxy             			*
c *								*
c * Role : lit une sous matrice(imin,imax,jmin,jmax)            *
c *         d'une variable reelle 2D lon, lat                   *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_subxy(nom_fichier,var_name,var,
c                    imin,imax,jmin,jmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_creer4_xyt             			*
c *								*
c * Role : creer une variable reelle 3D = lon, lat, temps       * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer4_xyt(nom_fichier,var_name,units,long_name,
c                      scale_factor,add_offset)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name
c      real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_xyt             			*
c *								*
c * Role : ecrit une variable reelle 3D = lon, lat, temps       *
c *        a partir d'une variable reelle 2D = lon, lat         * 
c *        et de l'indice du temps                              * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      ionc_ecritr4_xyt(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_ecritr4_sublxyt             		*
c *								*
c * Role : ecrit une variable reelle 3D = lon, lat, temps       *
c *        a partir d'une variable reelle 2D = lon, lat         * 
c *        (ilmin:ilmin+(imax-imin),jlmim:jlmin+(jmax-jmin))    *
c *        et de l'indice du temps                              * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c * 		- ilmin       :  borne min var d'entree	        *
c * 		- jlmin       :  borne min var d'entree	        *
c ***************************************************************
c ionc_ecritr4_sublxyt(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,nrec,ilmin,jlmin)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,jmin,jmax,ilmin,jlmin	
c        real           var(ilmin:ilmin+(imax-imin),
c                           jlmin:jlmin+(jmax-jmin))
c        integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_xyt             			*
c *								*
c * Role : lit une variable reelle 3D = lon, lat, temps         *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enrtegistrement       *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_xyt(nom_fichier,var_name,var,
c                     imin,imax,jmin,jmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_creer8_xyt             			*
c *								*
c * Role : creer une variable double 3D = lon, lat, temps       * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer8_xyt(nom_fichier,var_name,units,long_name,
c                      scale_factor,add_offset)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name
c      real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr8_xyt             			*
c *								*
c * Role : ecrit une variable double 3D = lon, lat, temps       *
c *        a partir d'une variable double 2D = lon, lat         * 
c *        et de l'indice du temps                              * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      ionc_ecritr8_xyt(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real*8         var(imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr8_xyt             			*
c *								*
c * Role : lit une variable double 3D = lon, lat, temps         *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enrtegistrement       *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr8_xyt(nom_fichier,var_name,var,
c                     imin,imax,jmin,jmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real*8         var(imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_subxyt             			*
c *								*
c * Role : lit une sous matrice d'une variable reelle 3D        *
c *                                     lon, lat, temps         *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enrtegistrement       *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_subxyt(nom_fichier,var_name,var,
c                     imin,imax,jmin,jmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_sublxyt             			*
c *								*
c * Role : lit une sous matrice d'une variable reelle 3D        *
c * lon, lat, temps dans une variable(ilmin:ilmin+(imax-imin),*
c *				      jlmim:jlmin+(jmax-jmin))*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- nrec        :  numero d'enrtegistrement       *
c * 		- ilmin       :  borne imin variable de lecture *
c * 		- jlmin       :  borne jmin variable de lecture *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c   ionc_litr4_sublxyt(nom_fichier,var_name,var,
c                                imin,imax,jmin,jmax,nrec,
c                                ilmin,ilmax,jlmin,jlmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,jmin,jmax
c        real           var(ilmin:ilmin+(imax-imin),
c                           jlmin:jlmin+(jmax-jmin))
c        integer        nrec

c ***************************************************************
c * subroutine ionc_creer4_zxy             			*
c *								*
c * Role : creer une variable reelle 3D = lon, lat, z           * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer4_zxy(nom_fichier,var_name,units,long_name,
c                      scale_factor,add_offset)
c ******** parametres de la subroutine *******
c	character*(*) nom_fichier
c	character*(*) var_name,units,long_name
c       real*4        scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_zxy             			*
c *								*
c * Role : ecrit une variable reelle 3D = lon, lat, z           *
c *        a partir d'une variable reelle 3D = z , lon, lat     * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c ***************************************************************
c      ionc_ecritr4_zxy(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,kmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(kmax,imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_litr4_zxy             			*
c *								*
c * Role : lit une variable reelle 3D = lon, lat, z             *
c *        dans une variable reelle 3D = z , lon, lat           * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_zxy(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,kmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(kmax,imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_creer8_zxy             			*
c *								*
c * Role : creer une variable double 3D = lon, lat, z           * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer8_zxy(nom_fichier,var_name,units,long_name,
c                      scale_factor,add_offset)
c ******** parametres de la subroutine *******
c	character*(*) nom_fichier
c	character*(*) var_name,units,long_name
c       real*4        scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_zxy             			*
c *								*
c * Role : ecrit une variable double 3D = lon, lat, z           *
c *        a partir d'une variable double 3D = z , lon, lat     * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c ***************************************************************
c      ionc_ecritr8_zxy(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,kmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real*8         var(kmax,imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_litr8_zxy             			*
c *								*
c * Role : lit une variable double 3D = lon, lat, z             *
c *        dans une variable double 3D = z , lon, lat           * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr8_zxy(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,kmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real*8         var(kmax,imin:imax,jmin:jmax)

c ***************************************************************
c * subroutine ionc_ecritr4_xyz             			*
c *								*
c * Role : ecrit une variable reelle 3D = lon, lat, z           *
c *        a partir d'une variable reelle 3D = lon, lat,z       * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c ***************************************************************
c      ionc_ecritr4_xyz(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,kmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax,kmax)

c ***************************************************************
c * subroutine ionc_litr4_xyz             			*
c *								*
c * Role : lit une variable reelle 3D = lon, lat, z             *
c *        dans une variable reelle 3D = lon, lat,z             * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_xyz(nom_fichier,var_name,var,
c                       imin,imax,jmin,jmax,kmax)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax
c      real           var(imin:imax,jmin:jmax,kmax)

c ***************************************************************
c * subroutine ionc_creer4_zxyt             			*
c *								*
c * Role : creer une variable reelle 4D = lon, lat, z, temps    * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer4_zxyt(nom_fichier,var_name,units,long_name,
c                       scale_factor,add_offset)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name
c      real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_zxyt             			*
c *								*
c * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
c *        a partir d'une variable reelle 3D = z , lon, lat     * 
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      ionc_ecritr4_zxyt(nom_fichier,var_name,var,
c                        imin,imax,jmin,jmax,kmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax,kmax
c      real           var(kmax,imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_ecritr4_sublzxyt             		*
c *								*
c * Role : ecrit une variable reelle 4D = lon, lat, z, temps    *
c *        a partir d'une variable reelle 3D = z , lon, lat     * 
c *        (klmin:klmin+(kmax-kmin),                            *
c *         ilmin:ilmin+(imax-imin),                            *
c *         jlmim:jlmin+(jmax-jmin))                            *
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmin        :  borne min dim k	        *
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c * 		- ilmin       :  borne min var d'entree	        *
c * 		- jlmin       :  borne min var d'entree	        *
c * 		- klmin       :  borne min var d'entree	        *
c ***************************************************************
c     ionc_ecritr4_sublzxyt(nom_fichier,var_name,var,
c                                imin,imax,jmin,jmax,kmin,kmax,nrec,
c                                ilmin,jlmin,klmin)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,jmin,jmax,kmin,kmax
c        integer        ilmin,jlmin,klmin
c        real           var(klmin:klmin+(kmax-kmin),
c                          ilmin:ilmin+(imax-imin),
c                          jlmin:jlmin+(jmax-jmin))
c        integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_zxyt             			*
c *								*
c * Role : lit une variable reelle 4D = lon, lat, z, temps      *
c *        dans une variable reelle 3D = z , lon, lat           * 
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_zxyt(nom_fichier,var_name,var,
c                      imin,imax,jmin,jmax,kmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax,kmax
c      real           var(kmax,imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_subzxyt             			*
c *								*
c * Role : lit une sous matrice (k,i,j) d'une variable reelle 4D*
c *                                     lon, lat, z, temps      *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmin        :  borne min dim k	        *
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enrtegistrement       *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_subzxyt(nom_fichier,var_name,var,
c                           imin,imax,jmin,jmax,kmin,kmax,nrec)
c
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,jmin,jmax, kmin, kmax
c        real           var(kmin:kmax,imin:imax,jmin:jmax)
c        integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_sublzxyt             			*
c *								*
c * Role : lit une sous matrice d'une variable reelle 3D        *
c * lon, lat, temps dans une variable(klmin:klmin+(kmax-kmin),  *
c *                                   ilmin:ilmin+(imax-imin),  *
c *				      jlmim:jlmin+(jmax-jmin))  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmin        :  borne min dim k	        *
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enrtegistrement       *
c * 		- ilmin       :  borne imin variable de lecture *
c * 		- jlmin       :  borne jmin variable de lecture *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c       ionc_litr4_sublzxyt(nom_fichier,var_name,var,
c                                imin,imax,jmin,jmax,kmin,kmax,nrec,
c                                ilmin,jlmin,klmin)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,jmin,jmax,kmin,kmax
c        integer        ilmin,jlmin,klmin
c        real           var(klmin:klmin+(kmax-kmin),
c     	                   ilmin:ilmin+(imax-imin),
c                          jlmin:jlmin+(jmax-jmin))
c        integer        nrec

c ***************************************************************
c * subroutine ionc_creer8_zxyt             			*
c *								*
c * Role : creer une variable double 4D = lon, lat, z, temps    * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer8_zxyt(nom_fichier,var_name,units,long_name,
c                       scale_factor,add_offset)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name,units,long_name
c      real*4         scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr8_zxyt             			*
c *								*
c * Role : ecrit une variable double 4D = lon, lat, z, temps    *
c *        a partir d'une variable double 3D = z , lon, lat     * 
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c      ionc_ecritr8_zxyt(nom_fichier,var_name,var,
c                        imin,imax,jmin,jmax,kmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax,kmax
c      real*8         var(kmax,imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr8_zxyt             			*
c *								*
c * Role : lit une variable double 4D = lon, lat, z, temps      *
c *        dans une variable double 3D = z , lon, lat           * 
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- jmin        :  borne min dim j	        *
c * 		- jmax        :  borne max dim j		*
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr8_zxyt(nom_fichier,var_name,var,
c                      imin,imax,jmin,jmax,kmax,nrec)
c ******** parametres de la subroutine *******
c      character*(*)  nom_fichier
c      character*(*)  var_name
c      integer        imin,imax,jmin,jmax,kmax
c      real*8         var(kmax,imin:imax,jmin:jmax)
c      integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_zxt             			*
c *								*
c * Role : lit une variable reelle 3D = lon, z, temps           * 
c *        dans une variable reelle 2D = z , lon                * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enrtegistrement       *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      subroutine ionc_litr4_zxt(nom_fichier,var_name,var,
c                                imin,imax,kmax,nrec)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*)  var_name
c       integer        imin,imax,kmax
c       real           var(kmax,imin:imax)
c       integer        nrec

c ***************************************************************
c * subroutine ionc_creer4_stat             			*
c *								*
c * Role : creer une variable reelle 2D = lon, temps            * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c     ionc_creer4_stat(nom_fichier,var_name,units,long_name,
c                         scale_factor,add_offset)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*) var_name,units,long_name
c        real*4 scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_stat             			*
c *								*
c * Role : ecrit une variable reelle 2D = lon,  temps           *
c *        a partir d'une variable reelle 2D = lon              * 
c *        et de l'indice du temps                              * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c    ionc_ecritr4_stat(nom_fichier,var_name,var,
c                               imin,imax,nrec)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax
c        real           var(imin:imax)
c        integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_stat             			*
c *								*
c * Role : lit une variable reelle 2D = lon, temps              * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- nrec        :  numero d'enrtegistrement       *
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_stat(nom_fichier,var_name,var,
c                              imin,imax,nrec)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax
c        real           var(imin:imax)
c        integer        nrec

c ***************************************************************
c * subroutine ionc_creer4_zstat             			*
c *								*
c * Role : creer une variable reelle 3D = nb station, z, temps  * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c      ionc_creer4_zstat(nom_fichier,var_name,units,long_name,
c                       scale_factor,add_offset)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*) var_name,units,long_name
c        real*4 scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_zstat             			*
c *								*
c * Role : ecrit une variable reelle 3D = nb stations, z, temps *
c *        a partir d'une variable reelle 2D = z , nb stations  * 
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c * 		                 0 si t est illimité		*
c ***************************************************************
c     ionc_ecritr4_zstat(nom_fichier,var_name,var,
c                           imin,imax,kmax,nrec)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,kmax
c        real           var(kmax,imin:imax)
c        integer        nrec

c ***************************************************************
c * subroutine ionc_litr4_zstat             			*
c *								*
c * Role : lit une variable reelle 3D = nb stations, z, temps   *
c *        dans une variable reelle 2D = z , nb stations        * 
c *	   l'indice du temps est fixe par nrec                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- imin        :  borne min dim i	        *
c * 		- imax        :  borne max dim i	        *
c * 		- kmax        :  borne max dim k		*
c * 		- nrec        :  numero d'enregistrement	*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c      ionc_litr4_zstat(nom_fichier,var_name,var,
c                           imin,imax,kmax,nrec)
c ******** PARAMETRES DE LA SUBROUTINE *******
c
c	character*(*)  nom_fichier
c	character*(*)  var_name
c        integer        imin,imax,kmax
c        real           var(kmax,imin:imax)
c	integer        nrec

c ***************************************************************
c * subroutine ionc_loupe                			*
c *								*
c * Role : lit les parametres de loupe dans le fichier netcdf   *
c *								*
c * Parametres :						*
c * 	entree:	- nc_id :  id du fichier		        *
c *	sortie:	- limi : limite min de l'indice i de la loupe   *
c *		- lima : limite max de l'indice i de la loupe   *
c *		- pasi : pas de l'indice i de la loupe          *     
c *		- ljmi : limite min de l'indice j de la loupe   *
c *		- ljma : limite max de l'indice j de la loupe	*
c *		- pasj : pas de l'indice j de la loupe          *
c *		- lkmi : limite min de l'indice k de la loupe   *
c *		- lkma : limite max de l'indice k de la loupe   *
c *		- pask : pas de l'indice k de la loupe          *
c ***************************************************************
c      ionc_litloupe(nc_id)
c                    limi,lima,pasi,ljmi,ljma,pasj,
c                    lkmi,lkma,pask)
c ******** parametres de la subroutine *******
c      integer nc_id
c      integer limi,lima,pasi
c      integer ljmi,ljma,pasj
c      integer lkmi,lkma,pask

c ***************************************************************
c * fonction ionc_flimax            			        *
c *								*
c * Role : lit l'indice max des longitudes                      *
c *								*
c * Parametres :						*
c * 	entree:	- nc_id :  id du fichier		        *
c * Resultat :							*
c * 	indice longitude max					*
c ***************************************************************
c      integer function ionc_flimax(nc_id)
c ******** parametres de la fonction *******
c      integer nc_id
       integer ionc_flimax
       external ionc_flimax

c ***************************************************************
c * fonction ionc_fljmax            			        *
c *								*
c * Role : lit l'indice max des latitudes                       *
c *								*
c * Parametres :						*
c * 	entree:	- nc_id :  id du fichier		        *
c * Resultat :							*
c * 	indice latitude max					*
c ***************************************************************
c      integer function ionc_fljmax(nc_id)
c ******** parametres de la fonction *******
c      integer nc_id
       integer ionc_fljmax
       external ionc_fljmax

c ***************************************************************
c * fonction ionc_flkmax            			        *
c *								*
c * Role : lit l'indice max des niveaux                         *
c *								*
c * Parametres :						*
c * 	entree:	- nc_id :  id du fichier		        *
c * Resultat :							*
c * 	indice niveaux max					*
c ***************************************************************
c      integer function ionc_flkmax(nc_id)
c ******** parametres de la fonction *******
c      integer nc_id
       integer ionc_flkmax
       external ionc_flkmax

c ***************************************************************
c * fonction ionc_fltmax            			        *
c *								*
c * Role : lit l'indice max des temps                           *
c *								*
c * Parametres :						*
c * 	entree:	- nc_id :  id du fichier		        *
c * Resultat :							*
c * 	indice temps max					*
c ***************************************************************
c      integer function ionc_fltmax(nc_id)
c ******** parametres de la fonction *******
c      integer nc_id
        integer ionc_fltmax
        external ionc_fltmax

c ***************************************************************
c * function ionc_fnlimax            			        *
c *								*
c * Role : lit l'indice max des longitudes                      *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * Resultat :							*
c * 	indice longitude max					*
c ***************************************************************
c      integer function ionc_fnlimax(nom_fichier)
c ******** parametres de la fonction *******
c        character*(*) nom_fichier
       integer ionc_fnlimax
       external ionc_fnlimax

c ***************************************************************
c * fonction ionc_fnljmax            			        *
c *								*
c * Role : lit l'indice max des latitudes                       *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * Resultat :							*
c * 	indice latitude max					*
c ***************************************************************
c      integer function ionc_fnljmax(nom_fichier)
c ******** parametres de la fonction *******
c        character*(*) nom_fichier
       integer ionc_fnljmax
       external ionc_fnljmax

c ***************************************************************
c * fonction ionc_fnlkmax            			        *
c *								*
c * Role : lit l'indice max des niveaux                         *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * Resultat :							*
c * 	indice niveau max					*
c ***************************************************************
c      integer function ionc_flkmax(nom_fichier)
c ******** parametres de la fonction *******
c        character*(*) nom_fichier
       integer ionc_fnlkmax
       external ionc_fnlkmax

c ***************************************************************
c * function ionc_fnltmax            			        *
c *								*
c * Role : lit l'indice max des temps                           *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * Resultat :							*
c * 	indice temps max					*
c ***************************************************************
c      integer function ionc_fnltmax(nom_fichier)
c ******** parametres de la fonction *******
c        character*(*) nom_fichier
        integer ionc_fnltmax
        external ionc_fnltmax

c ***************************************************************
c * function ionc_fnlsta                                        *
c *                                                             *
c * Role : lit l'indice max des stations                        *
c *                                                             *
c * Parametres :                                                *
c *     entree: - nom_fichier :  nom du fichier                 *
c * Resultat :                                                  *
c *     indice nombre de stations max                           *
c ***************************************************************
c      integer function ionc_fnlsta(nom_fichier)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        character*(*) nom_fichier
        integer ionc_fnlsta
        external ionc_fnlsta


c ***************************************************************
c * subroutine ionc_ouvrel             				*
c *								*
c * Role : ouverture d'un fichier netcdf en lecture		*
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c ***************************************************************
c      ionc_ouvrel(nom_fichier)
c ******** parametres de la subroutine *******
c      character*(*) nom_fichier

c ***************************************************************
c * subroutine ionc_ouvre             				*
c *								*
c * Role : ouverture d'un fichier netcdf en lecture/ecriture	*
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c ***************************************************************
c      ionc_ouvre(nom_fichier)
c ******** parametres de la subroutine *******
c      character*(*) nom_fichier

c ***************************************************************
c * subroutine ionc_fermer             				*
c *								*
c * Role : fermeture d'un fichier netcdf en lecture		*
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c ***************************************************************
c      ionc_fermer(nom_fichier)

c ***************************************************************
c * subroutine ionc_fermertout             			*
c *								*
c * Role : fermeture de tout les fichiers netcdf         	*
c *								*
c * Parametres :						*
c ***************************************************************
c      ionc_fermertout()

c ***************************************************************
c * subroutine ionc_creer4_z             			*
c *								*
c * Role : creer une variable reelle 1D = z                     * 
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- units       :  unite			        *
c * 		- long_name   :  nom long de la variable        *
c * 		- scale_factor:  coeff de multiplication        *
c * 		- add_offset  :  offset			        *
c ***************************************************************
c ionc_creer4_z(nom_fichier,var_name,units,long_name,
c               scale_factor,add_offset)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*) var_name,units,long_name
c       real*4 scale_factor,add_offset

c ***************************************************************
c * subroutine ionc_ecritr4_z             			*
c *								*
c * Role : ecrit une variable reelle 1D =  z(kmax)              *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- kmax        :  borne max dim k		*
c ***************************************************************
c    ionc_ecritr4_z(nom_fichier,var_name,var,kmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*)  var_name
c       integer        kmax
c       real           var(kmax)

c ***************************************************************
c * subroutine ionc_putr4_z	          			*
c *								*
c * Role :   ecrit une variable reelle 1D =  z(kmax)            * 
c *								*
c * Parametres :						*
c ***************************************************************
c   ionc_putr4_z(nc_id,var_id,var,kmax,
c                dimk,scale_factor,add_offset)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        integer        nc_id,var_id
c        integer        kmax
c        integer        dimk
c        real           var(kmax)
c        real*4         add_offset,scale_factor

c ***************************************************************
c * subroutine ionc_litr4_z            	        		*
c *								*
c * Role : lit une variable reelle 1D =  z 	                *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- kmax        :  borne max dim k		*
c *     sortie: - var         :  variable	                *
c ***************************************************************
c    ionc_litr4_z(nom_fichier,var_name,var,kmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*)  var_name
c       integer        kmax
c       real           var(kmax)

c ***************************************************************
c * subroutine ionc_ecritr4_sig             			*
c *								*
c * Role : ecrit une variable reelle 1D =  z(0:kmax+1)          *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 		- var_name    :  nom de la variable	        *
c * 		- var         :  variable	                *
c * 		- kmax        :  borne max dim k		*
c ***************************************************************
c   ionc_ecritr4_sig(nom_fichier,var_name,var,kmax)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*)  var_name
c       integer        kmax
c       real           var(0:kmax+1)

c ***************************************************************
c * subroutine ionc_putr4_sig	          			*
c *								*
c * Role :   ecrit une variable reelle 1D =  z(0:kmax+1)        * 
c *								*
c * Parametres :						*
c ***************************************************************
c    ionc_putr4_sig(nc_id,var_id,var,kmax,
c                   dimk,scale_factor,add_offset)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        integer        nc_id,var_id
c        integer        kmax
c        integer        dimk
c        real           var(0:kmax+1)
c        real*4         add_offset,scale_factor

c ***************************************************************
c * subroutine ionc_littorigine            			*
c *								*
c * Role : lecture de l'attribut time_origin de la varible time *
c *        d'un fichier netcdf                                  *
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c * 	sortie:	- torigine    :  origine du temps               *
c ***************************************************************
c       ionc_littorigine(nom_fichier,torigine)
c ******** PARAMETRES DE LA SUBROUTINE *******
c	character*(*)  nom_fichier
c	character*(*)  torigine

c ***************************************************************
c * subroutine ionc_litdimi            			        *
c *								*
c * Role : lit la taille de la grille en i                      *
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c * 	sortie:	- dimi  :  taille de la grille en i	        *
c ***************************************************************
c      ionc_litdimi(nom_fichier,dimi)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        character*(*) nom_fichier
c        integer dimi

c ***************************************************************
c * subroutine ionc_litdimj            			        *
c *								*
c * Role : lit la taille de la grille en j                      *
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c * 	sortie:	- dimj  :  taille de la grille en j	        *
c ***************************************************************
c      ionc_litdimj(nom_fichier,dimj)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        character*(*) nom_fichier
c        integer dimj

c ***************************************************************
c * subroutine ionc_litdimk            			        *
c *								*
c * Role : lit la taille de la grille en k                      *
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c * 	sortie:	- dimk  :  taille de la grille en k	        *
c ***************************************************************
c      ionc_litdimk(nom_fichier,dimk)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        character*(*) nom_fichier
c        integer dimk

c ***************************************************************
c * subroutine ionc_litdiml            			        *
c *								*
c * Role : lit la taille de la grille en l (temps)              *
c *								*
c * Parametres :						*
c * 	entree:	- nomfich   :  nom du fichier    	        *
c * 	sortie:	- dimt  :  taille de la grille en l (temps)     *
c ***************************************************************
c      ionc_litdiml(nom_fichier,diml)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        character*(*) nom_fichier
c        integer diml

c ***************************************************************
c * subroutine ionc_showvar            			        *
c *								*
c * Role : affichage des variables du fichier                   *
c *	sans les variables lat, lon, z et time			*
c *								*
c * Parametres :						*
c * 	entree:	- nom_fichier :  nom du fichier		        *
c ***************************************************************
c    ionc_showvar(nom_fichier)
c ******** PARAMETRES DE LA SUBROUTINE *******
c        character*(*) nom_fichier
