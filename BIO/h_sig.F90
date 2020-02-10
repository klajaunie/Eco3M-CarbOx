      Subroutine h_sig(i_arg,j_arg,k_arg,ZSUP,ZINF)
!***********************************************************************
!-- calcul des bornes superieures et inférieures de la couche sigma k_arg
!
!********************************************************************
!
!     objet    :  initialisation du modèle biogeochimique
!
!
!     appelant         : f_ppb_vincent.inc (inclus dans mod_process )
!                      :
!     appeles 1        : eco3m, sub_init, sub_init_mod,Sub_var_ini,sub_phys_ini
!     appeles 2        : sortie_bio, sortie_temp,sortie_bio_autre
!     appeles 3        : geom_bio, sub_affiche
!     variables in     :  mode_bio ( 0= INIT BIO ; 1 = CALCUL BIO )
!***********************************************************************

 use PHY_BIO
 use parametres
 use comvars2d, only: h0      
 use comvars3d, only: sig

implicit none    


 integer i_arg,j_arg,k_arg
 real    ZINF,ZSUP,H_eau


      H_eau=h0(i_arg,j_arg)+xe3bio(i_arg,j_arg)
!ajout 13/04
if(H_eau.gt.0.0) then        
      IF (k_arg.EQ.KMAX) THEN

              !ZINF=H_eau*(sig(k_arg)-1.)+xe3bio(i_arg,j_arg)
              ZINF=H_eau*(sig(k_arg-1)-1.)+xe3bio(i_arg,j_arg)
              !ZSUP=H_eau+xe3bio(i_arg,j_arg)
              ZSUP=xe3bio(i_arg,j_arg)

      ELSEIF(k_arg.EQ.1) THEN
              ZSUP=H_eau*(sig(k_arg)-1.)+xe3bio(i_arg,j_arg)
              ZINF=-1*H_eau+xe3bio(i_arg,j_arg)
      ELSE
              ZINF=H_eau*(sig(k_arg)-1.)+xe3bio(i_arg,j_arg)
              ZSUP=H_eau*(sig(k_arg+1)-1.)+xe3bio(i_arg,j_arg)

      END IF
endif


      return
      end subroutine h_sig

