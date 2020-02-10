 !-- Ce fichier est genere automatiquement lors de l execution du programme
 !-- a partir des informations du fichier config.ini   --
 !-- Il contient les instructions necessaires a la sauvegarde --
 !-- des differentes composantes d'un flux entre deux variables --
  
  
 dimfonc = size(FLUX_VAL(1,22)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_ppb_geider(1,22, &
      FLUX_PAR(1)%valpar(      1 ), & 
      FLUX_PAR(1)%valpar(      2 ), & 
      FLUX_PAR(1)%valpar(      3 ), & 
      FLUX_PAR(1)%valpar(      4 ), & 
      FLUX_PAR(1)%valpar(      5 ), & 
      FLUX_PAR(1)%valpar(      6 ), & 
      FLUX_PAR(1)%valpar(      7 ), & 
      FLUX_PAR(1)%valpar(      8 ), & 
      FLUX_PAR(1)%valpar(      9 ), & 
      FLUX_PAR(1)%valpar(     10 ), & 
      FLUX_PAR(1)%valpar(     11 ) )
  
 ss_flux(           2 ,:,:,:) =   &
(+1) * &
 f_resp(1,22, &
      FLUX_PAR(2)%valpar(      1 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')1
 write(         203 ,'(I4,$)')22
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(1,22)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(5,22)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_ppb_geider(5,22, &
      FLUX_PAR(3)%valpar(      1 ), & 
      FLUX_PAR(3)%valpar(      2 ), & 
      FLUX_PAR(3)%valpar(      3 ), & 
      FLUX_PAR(3)%valpar(      4 ), & 
      FLUX_PAR(3)%valpar(      5 ), & 
      FLUX_PAR(3)%valpar(      6 ), & 
      FLUX_PAR(3)%valpar(      7 ), & 
      FLUX_PAR(3)%valpar(      8 ), & 
      FLUX_PAR(3)%valpar(      9 ), & 
      FLUX_PAR(3)%valpar(     10 ), & 
      FLUX_PAR(3)%valpar(     11 ) )
  
 ss_flux(           2 ,:,:,:) =   &
(+1) * &
 f_resp(5,22, &
      FLUX_PAR(4)%valpar(      1 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')5
 write(         203 ,'(I4,$)')22
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(5,22)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(2,18)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_droop(2,18, &
      FLUX_PAR(5)%valpar(      1 ), & 
      FLUX_PAR(5)%valpar(      2 ), & 
      FLUX_PAR(5)%valpar(      3 ), & 
      FLUX_PAR(5)%valpar(      4 ), & 
      FLUX_PAR(5)%valpar(      5 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')2
 write(         203 ,'(I4,$)')18
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(2,18)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(2,19)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_droop(2,19, &
      FLUX_PAR(6)%valpar(      1 ), & 
      FLUX_PAR(6)%valpar(      2 ), & 
      FLUX_PAR(6)%valpar(      3 ), & 
      FLUX_PAR(6)%valpar(      4 ), & 
      FLUX_PAR(6)%valpar(      5 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')2
 write(         203 ,'(I4,$)')19
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(2,19)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(3,20)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_droop(3,20, &
      FLUX_PAR(7)%valpar(      1 ), & 
      FLUX_PAR(7)%valpar(      2 ), & 
      FLUX_PAR(7)%valpar(      3 ), & 
      FLUX_PAR(7)%valpar(      4 ), & 
      FLUX_PAR(7)%valpar(      5 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')3
 write(         203 ,'(I4,$)')20
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(3,20)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(6,18)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_droop(6,18, &
      FLUX_PAR(8)%valpar(      1 ), & 
      FLUX_PAR(8)%valpar(      2 ), & 
      FLUX_PAR(8)%valpar(      3 ), & 
      FLUX_PAR(8)%valpar(      4 ), & 
      FLUX_PAR(8)%valpar(      5 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')6
 write(         203 ,'(I4,$)')18
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(6,18)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(6,19)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_droop(6,19, &
      FLUX_PAR(9)%valpar(      1 ), & 
      FLUX_PAR(9)%valpar(      2 ), & 
      FLUX_PAR(9)%valpar(      3 ), & 
      FLUX_PAR(9)%valpar(      4 ), & 
      FLUX_PAR(9)%valpar(      5 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')6
 write(         203 ,'(I4,$)')19
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(6,19)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(7,20)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_droop(7,20, &
      FLUX_PAR(10)%valpar(      1 ), & 
      FLUX_PAR(10)%valpar(      2 ), & 
      FLUX_PAR(10)%valpar(      3 ), & 
      FLUX_PAR(10)%valpar(      4 ), & 
      FLUX_PAR(10)%valpar(      5 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')7
 write(         203 ,'(I4,$)')20
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(7,20)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(1,15)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out(1,15) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')1
 write(         203 ,'(I4,$)')15
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(1,15)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(5,15)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out(5,15) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')5
 write(         203 ,'(I4,$)')15
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(5,15)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(2,16)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out(2,16) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')2
 write(         203 ,'(I4,$)')16
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(2,16)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(6,16)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out(6,16) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')6
 write(         203 ,'(I4,$)')16
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(6,16)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(3,17)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out(3,17) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')3
 write(         203 ,'(I4,$)')17
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(3,17)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(7,17)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out(7,17) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')7
 write(         203 ,'(I4,$)')17
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(7,17)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(9,22)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(+1) * &
 f_resp_bac(9,22, &
      FLUX_PAR(17)%valpar(      1 ), & 
      FLUX_PAR(17)%valpar(      2 ), & 
      FLUX_PAR(17)%valpar(      3 ), & 
      FLUX_PAR(17)%valpar(      4 ), & 
      FLUX_PAR(17)%valpar(      5 ), & 
      FLUX_PAR(17)%valpar(      6 ), & 
      FLUX_PAR(17)%valpar(      7 ), & 
      FLUX_PAR(17)%valpar(      8 ), & 
      FLUX_PAR(17)%valpar(      9 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')9
 write(         203 ,'(I4,$)')22
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(9,22)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(9,12)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(9,12, &
      FLUX_PAR(18)%valpar(      1 ), & 
      FLUX_PAR(18)%valpar(      2 ), & 
      FLUX_PAR(18)%valpar(      3 ), & 
      FLUX_PAR(18)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')9
 write(         203 ,'(I4,$)')12
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(9,12)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(9,15)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(9,15, &
      FLUX_PAR(19)%valpar(      1 ), & 
      FLUX_PAR(19)%valpar(      2 ), & 
      FLUX_PAR(19)%valpar(      3 ), & 
      FLUX_PAR(19)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')9
 write(         203 ,'(I4,$)')15
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(9,15)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(10,13)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(10,13, &
      FLUX_PAR(20)%valpar(      1 ), & 
      FLUX_PAR(20)%valpar(      2 ), & 
      FLUX_PAR(20)%valpar(      3 ), & 
      FLUX_PAR(20)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')10
 write(         203 ,'(I4,$)')13
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(10,13)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(10,16)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(10,16, &
      FLUX_PAR(21)%valpar(      1 ), & 
      FLUX_PAR(21)%valpar(      2 ), & 
      FLUX_PAR(21)%valpar(      3 ), & 
      FLUX_PAR(21)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')10
 write(         203 ,'(I4,$)')16
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(10,16)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(10,18)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(10,18, &
      FLUX_PAR(22)%valpar(      1 ), & 
      FLUX_PAR(22)%valpar(      2 ), & 
      FLUX_PAR(22)%valpar(      3 ), & 
      FLUX_PAR(22)%valpar(      4 ) )
  
 ss_flux(           2 ,:,:,:) =   &
(+1) * &
 f_remin(10,18, &
      FLUX_PAR(23)%valpar(      1 ), & 
      FLUX_PAR(23)%valpar(      2 ), & 
      FLUX_PAR(23)%valpar(      3 ), & 
      FLUX_PAR(23)%valpar(      4 ), & 
      FLUX_PAR(23)%valpar(      5 ), & 
      FLUX_PAR(23)%valpar(      6 ), & 
      FLUX_PAR(23)%valpar(      7 ), & 
      FLUX_PAR(23)%valpar(      8 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')10
 write(         203 ,'(I4,$)')18
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(10,18)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(11,14)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(11,14, &
      FLUX_PAR(24)%valpar(      1 ), & 
      FLUX_PAR(24)%valpar(      2 ), & 
      FLUX_PAR(24)%valpar(      3 ), & 
      FLUX_PAR(24)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')11
 write(         203 ,'(I4,$)')14
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(11,14)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(11,17)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(11,17, &
      FLUX_PAR(25)%valpar(      1 ), & 
      FLUX_PAR(25)%valpar(      2 ), & 
      FLUX_PAR(25)%valpar(      3 ), & 
      FLUX_PAR(25)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')11
 write(         203 ,'(I4,$)')17
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(11,17)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(11,20)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_upt_bact_droop(11,20, &
      FLUX_PAR(26)%valpar(      1 ), & 
      FLUX_PAR(26)%valpar(      2 ), & 
      FLUX_PAR(26)%valpar(      3 ), & 
      FLUX_PAR(26)%valpar(      4 ) )
  
 ss_flux(           2 ,:,:,:) =   &
(+1) * &
 f_remin(11,20, &
      FLUX_PAR(27)%valpar(      1 ), & 
      FLUX_PAR(27)%valpar(      2 ), & 
      FLUX_PAR(27)%valpar(      3 ), & 
      FLUX_PAR(27)%valpar(      4 ), & 
      FLUX_PAR(27)%valpar(      5 ), & 
      FLUX_PAR(27)%valpar(      6 ), & 
      FLUX_PAR(27)%valpar(      7 ), & 
      FLUX_PAR(27)%valpar(      8 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')11
 write(         203 ,'(I4,$)')20
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(11,20)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(18,19)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(+1) * &
 f_nit(18,19, &
      FLUX_PAR(28)%valpar(      1 ), & 
      FLUX_PAR(28)%valpar(      2 ), & 
      FLUX_PAR(28)%valpar(      3 ), & 
      FLUX_PAR(28)%valpar(      4 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')18
 write(         203 ,'(I4,$)')19
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(18,19)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(21,22)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(+1) * &
 f_out_mod(21,22, &
      FLUX_PAR(29)%valpar(      1 ), & 
      FLUX_PAR(29)%valpar(      2 ), & 
      FLUX_PAR(29)%valpar(      3 ), & 
      FLUX_PAR(29)%valpar(      4 ), & 
      FLUX_PAR(29)%valpar(      5 ), & 
      FLUX_PAR(29)%valpar(      6 ), & 
      FLUX_PAR(29)%valpar(      7 ), & 
      FLUX_PAR(29)%valpar(      8 ), & 
      FLUX_PAR(29)%valpar(      9 ), & 
      FLUX_PAR(29)%valpar(     10 ), & 
      FLUX_PAR(29)%valpar(     11 ), & 
      FLUX_PAR(29)%valpar(     12 ), & 
      FLUX_PAR(29)%valpar(     13 ), & 
      FLUX_PAR(29)%valpar(     14 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')21
 write(         203 ,'(I4,$)')22
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(21,22)%val(1,1,1)
 
 deallocate (ss_flux)
 
 dimfonc = size(FLUX_VAL(22,26)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
  
  
 ss_flux(           1 ,:,:,:) =   &
(-1) * &
 f_disscarbo(22,26, &
      FLUX_PAR(30)%valpar(      1 ) )
  
 ss_flux(           2 ,:,:,:) =   &
(+1) * &
 f_precipcarbo(22,26, &
      FLUX_PAR(31)%valpar(      1 ), & 
      FLUX_PAR(31)%valpar(      2 ) )
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')22
 write(         203 ,'(I4,$)')26
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   FLUX_VAL(22,26)%val(1,1,1)
 
 deallocate (ss_flux)
 
 
 !--- 2/ Expressions des flux (non nuls) n'impliquant qu'une 
 !       seule variable d'état
 
 dimfonc = size(SELF_VAL(1)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_phy(1,0, &
      SELF_PAR(1)%valpar(      1 ), & 
      SELF_PAR(1)%valpar(      2 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')1
 write(         203 ,'(I4,$)')1
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(1)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(5)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_phy(5,0, &
      SELF_PAR(2)%valpar(      1 ), & 
      SELF_PAR(2)%valpar(      2 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')5
 write(         203 ,'(I4,$)')5
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(5)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(2)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_phy(2,0, &
      SELF_PAR(3)%valpar(      1 ), & 
      SELF_PAR(3)%valpar(      2 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')2
 write(         203 ,'(I4,$)')2
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(2)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(6)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_phy(6,0, &
      SELF_PAR(4)%valpar(      1 ), & 
      SELF_PAR(4)%valpar(      2 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')6
 write(         203 ,'(I4,$)')6
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(6)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(3)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_phy(3,0, &
      SELF_PAR(5)%valpar(      1 ), & 
      SELF_PAR(5)%valpar(      2 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')3
 write(         203 ,'(I4,$)')3
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(3)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(7)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_phy(7,0, &
      SELF_PAR(6)%valpar(      1 ), & 
      SELF_PAR(6)%valpar(      2 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')7
 write(         203 ,'(I4,$)')7
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(7)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(9)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_bact(9,0, &
      SELF_PAR(7)%valpar(      1 ), & 
      SELF_PAR(7)%valpar(      2 ), & 
      SELF_PAR(7)%valpar(      3 ), & 
      SELF_PAR(7)%valpar(      4 ), & 
      SELF_PAR(7)%valpar(      5 ), & 
      SELF_PAR(7)%valpar(      6 ), & 
      SELF_PAR(7)%valpar(      7 ), & 
      SELF_PAR(7)%valpar(      8 ), & 
      SELF_PAR(7)%valpar(      9 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')9
 write(         203 ,'(I4,$)')9
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(9)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(10)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_bact(10,0, &
      SELF_PAR(8)%valpar(      1 ), & 
      SELF_PAR(8)%valpar(      2 ), & 
      SELF_PAR(8)%valpar(      3 ), & 
      SELF_PAR(8)%valpar(      4 ), & 
      SELF_PAR(8)%valpar(      5 ), & 
      SELF_PAR(8)%valpar(      6 ), & 
      SELF_PAR(8)%valpar(      7 ), & 
      SELF_PAR(8)%valpar(      8 ), & 
      SELF_PAR(8)%valpar(      9 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')10
 write(         203 ,'(I4,$)')10
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(10)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(11)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-1) * &
   f_gra_bact(11,0, &
      SELF_PAR(9)%valpar(      1 ), & 
      SELF_PAR(9)%valpar(      2 ), & 
      SELF_PAR(9)%valpar(      3 ), & 
      SELF_PAR(9)%valpar(      4 ), & 
      SELF_PAR(9)%valpar(      5 ), & 
      SELF_PAR(9)%valpar(      6 ), & 
      SELF_PAR(9)%valpar(      7 ), & 
      SELF_PAR(9)%valpar(      8 ), & 
      SELF_PAR(9)%valpar(      9 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')11
 write(         203 ,'(I4,$)')11
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(11)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(12)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_det(12,0, &
      SELF_PAR(10)%valpar(      1 ), & 
      SELF_PAR(10)%valpar(      2 ), & 
      SELF_PAR(10)%valpar(      3 ), & 
      SELF_PAR(10)%valpar(      4 ), & 
      SELF_PAR(10)%valpar(      5 ), & 
      SELF_PAR(10)%valpar(      6 ), & 
      SELF_PAR(10)%valpar(      7 ), & 
      SELF_PAR(10)%valpar(      8 ), & 
      SELF_PAR(10)%valpar(      9 ), & 
      SELF_PAR(10)%valpar(     10 ), & 
      SELF_PAR(10)%valpar(     11 ), & 
      SELF_PAR(10)%valpar(     12 ) )  
   
  
 ss_flux(           2 ,:,:,:) =  &
(-1) * &
   f_gra_det(12,0, &
      SELF_PAR(11)%valpar(      1 ), & 
      SELF_PAR(11)%valpar(      2 ) )  
   
  
 ss_flux(           3 ,:,:,:) =  &
(+1) * &
   f_p_zoo(12,0, &
      SELF_PAR(12)%valpar(      1 ), & 
      SELF_PAR(12)%valpar(      2 ), & 
      SELF_PAR(12)%valpar(      3 ), & 
      SELF_PAR(12)%valpar(      4 ), & 
      SELF_PAR(12)%valpar(      5 ), & 
      SELF_PAR(12)%valpar(      6 ), & 
      SELF_PAR(12)%valpar(      7 ), & 
      SELF_PAR(12)%valpar(      8 ), & 
      SELF_PAR(12)%valpar(      9 ), & 
      SELF_PAR(12)%valpar(     10 ), & 
      SELF_PAR(12)%valpar(     11 ), & 
      SELF_PAR(12)%valpar(     12 ), & 
      SELF_PAR(12)%valpar(     13 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')12
 write(         203 ,'(I4,$)')12
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           3 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(12)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(13)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_det(13,0, &
      SELF_PAR(13)%valpar(      1 ), & 
      SELF_PAR(13)%valpar(      2 ), & 
      SELF_PAR(13)%valpar(      3 ), & 
      SELF_PAR(13)%valpar(      4 ), & 
      SELF_PAR(13)%valpar(      5 ), & 
      SELF_PAR(13)%valpar(      6 ), & 
      SELF_PAR(13)%valpar(      7 ), & 
      SELF_PAR(13)%valpar(      8 ), & 
      SELF_PAR(13)%valpar(      9 ), & 
      SELF_PAR(13)%valpar(     10 ), & 
      SELF_PAR(13)%valpar(     11 ), & 
      SELF_PAR(13)%valpar(     12 ) )  
   
  
 ss_flux(           2 ,:,:,:) =  &
(-1) * &
   f_gra_det(13,0, &
      SELF_PAR(14)%valpar(      1 ), & 
      SELF_PAR(14)%valpar(      2 ) )  
   
  
 ss_flux(           3 ,:,:,:) =  &
(+1) * &
   f_p_zoo(13,0, &
      SELF_PAR(15)%valpar(      1 ), & 
      SELF_PAR(15)%valpar(      2 ), & 
      SELF_PAR(15)%valpar(      3 ), & 
      SELF_PAR(15)%valpar(      4 ), & 
      SELF_PAR(15)%valpar(      5 ), & 
      SELF_PAR(15)%valpar(      6 ), & 
      SELF_PAR(15)%valpar(      7 ), & 
      SELF_PAR(15)%valpar(      8 ), & 
      SELF_PAR(15)%valpar(      9 ), & 
      SELF_PAR(15)%valpar(     10 ), & 
      SELF_PAR(15)%valpar(     11 ), & 
      SELF_PAR(15)%valpar(     12 ), & 
      SELF_PAR(15)%valpar(     13 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')13
 write(         203 ,'(I4,$)')13
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           3 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(13)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(14)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_det(14,0, &
      SELF_PAR(16)%valpar(      1 ), & 
      SELF_PAR(16)%valpar(      2 ), & 
      SELF_PAR(16)%valpar(      3 ), & 
      SELF_PAR(16)%valpar(      4 ), & 
      SELF_PAR(16)%valpar(      5 ), & 
      SELF_PAR(16)%valpar(      6 ), & 
      SELF_PAR(16)%valpar(      7 ), & 
      SELF_PAR(16)%valpar(      8 ), & 
      SELF_PAR(16)%valpar(      9 ), & 
      SELF_PAR(16)%valpar(     10 ), & 
      SELF_PAR(16)%valpar(     11 ), & 
      SELF_PAR(16)%valpar(     12 ) )  
   
  
 ss_flux(           2 ,:,:,:) =  &
(-1) * &
   f_gra_det(14,0, &
      SELF_PAR(17)%valpar(      1 ), & 
      SELF_PAR(17)%valpar(      2 ) )  
   
  
 ss_flux(           3 ,:,:,:) =  &
(+1) * &
   f_p_zoo(14,0, &
      SELF_PAR(18)%valpar(      1 ), & 
      SELF_PAR(18)%valpar(      2 ), & 
      SELF_PAR(18)%valpar(      3 ), & 
      SELF_PAR(18)%valpar(      4 ), & 
      SELF_PAR(18)%valpar(      5 ), & 
      SELF_PAR(18)%valpar(      6 ), & 
      SELF_PAR(18)%valpar(      7 ), & 
      SELF_PAR(18)%valpar(      8 ), & 
      SELF_PAR(18)%valpar(      9 ), & 
      SELF_PAR(18)%valpar(     10 ), & 
      SELF_PAR(18)%valpar(     11 ), & 
      SELF_PAR(18)%valpar(     12 ), & 
      SELF_PAR(18)%valpar(     13 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')14
 write(         203 ,'(I4,$)')14
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           3 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(14)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(15)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_mod(15,0, &
      SELF_PAR(19)%valpar(      1 ), & 
      SELF_PAR(19)%valpar(      2 ), & 
      SELF_PAR(19)%valpar(      3 ), & 
      SELF_PAR(19)%valpar(      4 ), & 
      SELF_PAR(19)%valpar(      5 ), & 
      SELF_PAR(19)%valpar(      6 ), & 
      SELF_PAR(19)%valpar(      7 ), & 
      SELF_PAR(19)%valpar(      8 ), & 
      SELF_PAR(19)%valpar(      9 ), & 
      SELF_PAR(19)%valpar(     10 ), & 
      SELF_PAR(19)%valpar(     11 ), & 
      SELF_PAR(19)%valpar(     12 ), & 
      SELF_PAR(19)%valpar(     13 ), & 
      SELF_PAR(19)%valpar(     14 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')15
 write(         203 ,'(I4,$)')15
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(15)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(16)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_mod(16,0, &
      SELF_PAR(20)%valpar(      1 ), & 
      SELF_PAR(20)%valpar(      2 ), & 
      SELF_PAR(20)%valpar(      3 ), & 
      SELF_PAR(20)%valpar(      4 ), & 
      SELF_PAR(20)%valpar(      5 ), & 
      SELF_PAR(20)%valpar(      6 ), & 
      SELF_PAR(20)%valpar(      7 ), & 
      SELF_PAR(20)%valpar(      8 ), & 
      SELF_PAR(20)%valpar(      9 ), & 
      SELF_PAR(20)%valpar(     10 ), & 
      SELF_PAR(20)%valpar(     11 ), & 
      SELF_PAR(20)%valpar(     12 ), & 
      SELF_PAR(20)%valpar(     13 ), & 
      SELF_PAR(20)%valpar(     14 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')16
 write(         203 ,'(I4,$)')16
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(16)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(17)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_mod(17,0, &
      SELF_PAR(21)%valpar(      1 ), & 
      SELF_PAR(21)%valpar(      2 ), & 
      SELF_PAR(21)%valpar(      3 ), & 
      SELF_PAR(21)%valpar(      4 ), & 
      SELF_PAR(21)%valpar(      5 ), & 
      SELF_PAR(21)%valpar(      6 ), & 
      SELF_PAR(21)%valpar(      7 ), & 
      SELF_PAR(21)%valpar(      8 ), & 
      SELF_PAR(21)%valpar(      9 ), & 
      SELF_PAR(21)%valpar(     10 ), & 
      SELF_PAR(21)%valpar(     11 ), & 
      SELF_PAR(21)%valpar(     12 ), & 
      SELF_PAR(21)%valpar(     13 ), & 
      SELF_PAR(21)%valpar(     14 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')17
 write(         203 ,'(I4,$)')17
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(17)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(18)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_mod(18,0, &
      SELF_PAR(22)%valpar(      1 ), & 
      SELF_PAR(22)%valpar(      2 ), & 
      SELF_PAR(22)%valpar(      3 ), & 
      SELF_PAR(22)%valpar(      4 ), & 
      SELF_PAR(22)%valpar(      5 ), & 
      SELF_PAR(22)%valpar(      6 ), & 
      SELF_PAR(22)%valpar(      7 ), & 
      SELF_PAR(22)%valpar(      8 ), & 
      SELF_PAR(22)%valpar(      9 ), & 
      SELF_PAR(22)%valpar(     10 ), & 
      SELF_PAR(22)%valpar(     11 ), & 
      SELF_PAR(22)%valpar(     12 ), & 
      SELF_PAR(22)%valpar(     13 ), & 
      SELF_PAR(22)%valpar(     14 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')18
 write(         203 ,'(I4,$)')18
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(18)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(20)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_out_mod(20,0, &
      SELF_PAR(23)%valpar(      1 ), & 
      SELF_PAR(23)%valpar(      2 ), & 
      SELF_PAR(23)%valpar(      3 ), & 
      SELF_PAR(23)%valpar(      4 ), & 
      SELF_PAR(23)%valpar(      5 ), & 
      SELF_PAR(23)%valpar(      6 ), & 
      SELF_PAR(23)%valpar(      7 ), & 
      SELF_PAR(23)%valpar(      8 ), & 
      SELF_PAR(23)%valpar(      9 ), & 
      SELF_PAR(23)%valpar(     10 ), & 
      SELF_PAR(23)%valpar(     11 ), & 
      SELF_PAR(23)%valpar(     12 ), & 
      SELF_PAR(23)%valpar(     13 ), & 
      SELF_PAR(23)%valpar(     14 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')20
 write(         203 ,'(I4,$)')20
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(20)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(21)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(-2.00   ) * &
   f_nit(21,0, &
      SELF_PAR(24)%valpar(      1 ), & 
      SELF_PAR(24)%valpar(      2 ), & 
      SELF_PAR(24)%valpar(      3 ), & 
      SELF_PAR(24)%valpar(      4 ) )  
   
  
 ss_flux(           2 ,:,:,:) =  &
(+1) * &
  f_prodO2(21,0) 
   
  
 ss_flux(           3 ,:,:,:) =  &
(+1) * &
   f_aera(21,0, &
      SELF_PAR(26)%valpar(      1 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')21
 write(         203 ,'(I4,$)')21
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           3 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(21)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(22)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+1) * &
   f_aera(22,0, &
      SELF_PAR(27)%valpar(      1 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')22
 write(         203 ,'(I4,$)')22
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 , '(A10,$)' )      ' NaN '    
 write(         203 ,'(E12.3E3)')   SELF_VAL(22)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
 dimfonc = size(SELF_VAL(25)%idproc)
 Allocate (ss_flux(dimfonc,nx_min:nx_max,ny_min:ny_max,nz_max))
  
 ss_flux(           1 ,:,:,:) =  &
(+2.     ) * &
   f_disscarbo(25,0, &
      SELF_PAR(28)%valpar(      1 ) )  
   
  
 ss_flux(           2 ,:,:,:) =  &
(-2.     ) * &
   f_precipcarbo(25,0, &
      SELF_PAR(29)%valpar(      1 ), & 
      SELF_PAR(29)%valpar(      2 ) )  
   
  
 ss_flux(           3 ,:,:,:) =  &
(-2.     ) * &
   f_nit(25,0, &
      SELF_PAR(30)%valpar(      1 ), & 
      SELF_PAR(30)%valpar(      2 ), & 
      SELF_PAR(30)%valpar(      3 ), & 
      SELF_PAR(30)%valpar(      4 ) )  
   
  
 ss_flux(           4 ,:,:,:) =  &
(+1) * &
   f_upt_droop(25,0, &
      SELF_PAR(31)%valpar(      1 ), & 
      SELF_PAR(31)%valpar(      2 ), & 
      SELF_PAR(31)%valpar(      3 ), & 
      SELF_PAR(31)%valpar(      4 ), & 
      SELF_PAR(31)%valpar(      5 ) )  
   
  
 ss_flux(           5 ,:,:,:) =  &
(+1) * &
   f_remin(25,0, &
      SELF_PAR(32)%valpar(      1 ), & 
      SELF_PAR(32)%valpar(      2 ), & 
      SELF_PAR(32)%valpar(      3 ), & 
      SELF_PAR(32)%valpar(      4 ), & 
      SELF_PAR(32)%valpar(      5 ), & 
      SELF_PAR(32)%valpar(      6 ), & 
      SELF_PAR(32)%valpar(      7 ), & 
      SELF_PAR(32)%valpar(      8 ) ) 
  
  
 write(         203 ,'(E10.2,$)') tps/3600. 
 write(         203 ,'(I4,$)')25
 write(         203 ,'(I4,$)')25
 write(         203 ,"(E10.2,$)") ss_flux(           1 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           2 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           3 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           4 ,1,1,1)
 write(         203 ,"(E10.2,$)") ss_flux(           5 ,1,1,1)
 write(         203 ,'(E12.3E3)')   SELF_VAL(25)%val(1,1,1)
 
 ss_flux=0.d0
 deallocate(ss_flux)
 
