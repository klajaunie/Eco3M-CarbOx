!***************************************************************************
!***************************************************************************
!Copyright or © or Copr. CNRS/IRD/Université de la Méditerranée
!contributor(s) : Melika BAKLOUTI & Vincent FAURE (10/10/2006)
!
!m.baklouti@univmed.fr; vincent.faure@univmed.fr
!
!This software (Eco3M) is a computer program whose purpose is to perform 
!biogeochemical or coupled physical-biogeochemical modelling.
!
!This software is governed by the CeCILL license under French law and
!abiding by the rules of distribution of free software. You can  use, 
!modify and/ or redistribute the software under the terms of the CeCILL
!license as circulated by CEA, CNRS and INRIA at the following URL
!"http://www.cecill.info". 
!
!As a counterpart to the access to the source code and  rights to copy,
!modify and redistribute granted by the license, users are provided only
!with a limited warranty  and the software''s author,  the holder of the
!economic rights,  and the successive licensors  have only  limited
!liability. 
!
!In this respect, the user''s attention is drawn to the risks associated
!with loading,  using,  modifying and/or developing or reproducing the
!software by the user in light of its specific status of free software,
!that may mean  that it is complicated to manipulate,  and  that  also
!therefore means  that it is reserved for developers  and  experienced
!professionals having in-depth computer knowledge. Users are therefore
!encouraged to load and test the software''s suitability as regards their
!requirements in conditions enabling the security of their systems and/or 
!data to be ensured and,  more generally, to use and operate it in the 
!same conditions as regards security. 
!
!The fact that you are presently reading this means that you have had
!knowledge of the CeCILL license and that you accept its terms.
!***************************************************************************
!***************************************************************************
!
!-------------------------------------------------------------------------------
! 
                         SUBROUTINE Sub_resbio_ascii

! Routine permettant d'ecrire au format ASCII l'evolution temporelle des 
! concentrations biogeochimiques calculees par le modele 
!
! derniere modification : 27/06/07
!
!-------------------------------------------------------------------------------
Use DEF_TYPE
Use VAR_USER
Use VAR_GLOBAL
Use COUPLEUR_PHY_BIO
Use MOD_FCHAIN

Implicit None

!-- Arguments     
real(8):: t_arg

! Variables locales
integer :: i,j,k,ivar,ifich,id,idorg_old
integer :: nfich1,nfich2
 Character(10)  :: fich1,fich2

id=0
idorg_old=-100

!-- Au premier appel au module bio --> Ouverture du fichier de sortie
if (nbcallbio == 0) then
do 
   id=id+1
   if (id > nbvar) exit
   if (VAR(id)%idorg /= idorg_old) then
     idorg_old=VAR(id)%idorg
!
     fich1 = VAR(id)%comp
     nfich1 = len_trim(adjustl(fich1))
!
     fich2 = VAR(id)%scomp
     nfich2 = len_trim(adjustl(fich2))
!

  !--Definition du nom du fichier de sortie et entete
     resu_bio = trim(rep_sortie) // trim(suffbio) //'_'// fich1(1:nfich1)//'_'//fich2(1:nfich2)//'.OUT'
    
     ifich = 213 + VAR(id)%idorg
     OPEN(ifich,FILE=resu_bio)
     write(ifich,'(/,a22,$)') '%  time(h) m  n  l ||'
   endif
       write(ifich,'($,a5,$a1,$a3,$a2)') VAR(id)%scomp, '/', VAR(id)%elmt, '||'
  enddo
OPEN(213,FILE=trim(rep_sortie) // trim(suffbio) //'_BP.OUT')
 write(213,'(/,a25,$)') '%  time(h) m  n  l || BP'
! marion rajout fichier sortie des limitation par les sel nut
OPEN(212,FILE=trim(rep_sortie) // trim(suffbio) //'_LIM_PICO.OUT')
 write(212,'(/,a180,$)') '%  time(h) m  n  l ||  phyto_lim_N phyto_lim_P phyto_lim_fin lim_T_graz lim_T_ppb lim_I'
OPEN(211,FILE=trim(rep_sortie) // trim(suffbio) //'_LIM_DIA.OUT')
 write(211,'(/,a180,$)') '%  time(h) m  n  l ||  phyto_lim_N phyto_lim_P phyto_lim_fin lim_T_graz lim_T_ppb lim_I'
endif 


!-- Pour les autres appels:

idorg_old=-100

#ifdef COUPL
do i=nx_min+1,nx_max-1,nint((nx_max-nx_min)/2.)-2
 do j=ny_min+1,ny_max-1,nint((ny_max-ny_min)/2.)-2
  do k=1,nz_max-1,nz_max/2
#endif

#ifdef NCOUPL
do i=nx_min,nx_max
 do j=ny_min,ny_max
  do k=1,nz_max
#endif
      id=0
      do
       id=id+1
        if (id > nbvar) exit
          ifich = 213 + VAR(id)%idorg
          if (VAR(id)%idorg /= idorg_old) then
           WRITE(ifich,'(/,F8.2,3I4,$)') tps / 3600., i,j,k
           idorg_old=VAR(id)%idorg
          endif
            WRITE(ifich,'($1X,E15.4)') VAR(id)%conc(i,j,k)
        enddo
           WRITE(240,'(/,F8.2,3I4,$)') tps / 3600., i,j,k
           write(240,'($1X,E15.4)') PPB(1)%val(i,j,k),PPB(2)%val(i,j,k),hQno(1)%val(i,j,k),hQno(2)%val(i,j,k),pr(1)%val(i,j,k),pr(2)%val(i,j,k),br(1)%val(i,j,k),PPB(1)%val(i,j,k)+PPB(2)%val(i,j,k)+2*hQno(1)%val(i,j,k)+2*hQno(2)%val(i,j,k)-pr(1)%val(i,j,k)-pr(2)%val(i,j,k)-br(1)%val(i,j,k)
           WRITE(213,'(/,F8.2,3I4,$)') tps / 3600., i,j,k
            WRITE(213,'($1X,E15.4)') bp(1)%val(i,j,k)
    ! marion ajout ecriture fichier limitation sel nut   
           WRITE(212,'(/,F8.2,3I4,$)') tps / 3600., i,j,k
            WRITE(212,'($1X,E15.4)') limn(1)%val(i,j,k),limp(1)%val(i,j,k),lim_nut(1)%val(i,j,k),limT_graz(1)%val(i,j,k),limT_ppb(1)%val(i,j,k),lim_lum(1)%val(i,j,k)
           WRITE(211,'(/,F8.2,3I4,$)') tps / 3600., i,j,k
            WRITE(211,'($1X,E15.4)') limn(2)%val(i,j,k),limp(2)%val(i,j,k),lim_nut(2)%val(i,j,k),limT_graz(2)%val(i,j,k),limT_ppb(2)%val(i,j,k),lim_lum(2)%val(i,j,k)
      enddo
    enddo
enddo


#ifdef MODTEST
   write(*,*) 'fin de la subroutine  sub_resbio_ascii'
#endif

End Subroutine sub_resbio_ascii
