Implicit none 
Integer:: err,nb,i,j
Character(25) :: fichpar,chain1,chain2
Character(120) :: chain
Character(50) :: chaine
real(8) :: num
Character(25),Allocatable::chaine1(:),chaine2(:)
logical::present

Write(*,*) '-------------------------------------------------------------------'
Write(*,*) 'A LIRE ---- A LIRE --- A LIRE --- A LIRE --- A LIRE --- A LIRE ----'
Write(*,*) 'Ce programme permet de convertir les noms des parametres  utilises'
Write(*,*) 'dans config.ini_txt en nombres suivant les associations nom/valeur'
Write(*,*) 'donnees dans le fichier de parametres (dont le nom doit etre '
Write(*,*) 'introduit au clavier). '
Write(*,*) 'Pour que la conversion réussisse, il est IMPORTANT :'
Write(*,*) '1/ q''aucun chiffre n''apparaisse dans les NOMS des parametres (par'
Write(*,*) 'exemple  AMMO peut etre utilise a la place de NH4)'
Write(*,*) '2/ que le fichier config.ini texte soit renomme config.ini_txt'
Write(*,*) 'avant de lancer le programme'
Write(*,*) 'Le resultat de la conversion est donne dans le fichier'
Write(*,*) 'config.ini_conv'
Write(*,*) 'A LIRE ---- A LIRE --- A LIRE --- A LIRE --- A LIRE --- A LIRE ----'
Write(*,*) '-------------------------------------------------------------------'

Read(*,*)

!-- fichier de commandes bash
open (10,file='txt2nb.sh')

!-- fichier de parametres
write(*,*)'nom du fichier contenant les valeurs des parametres ?'
read(*,*) fichpar
open (20,file=fichpar)

!-- 1ere lecture pour longueur du fichier de parametres:
nb=0
do
! read(20,*,IOSTAT=err) chain,num
 read(20,*,IOSTAT=err) chain1,chain2
 if (err /=0) exit
 nb=nb+1
enddo

Inquire(file='config.ini_txt',Exist=Present)
if (.not. present) then
 write(*,*) 'le fichier config.ini_txt n existe pas, ce programme va être arreté '
 stop
endif

write(10,*) 'cp config.ini_txt fc'

Allocate(chaine1(nb),chaine2(nb))
rewind(20)

!-- Classement des parametres par ordre de chaines numeriques
!   decroissantes:

do i=1,nb
  read(20,*,IOSTAT=err) chaine1(i),chaine2(i)
          chain1 = trim(adjustl(chaine1(i)))
          chain2 = trim(adjustl(chaine2(i)))
  if (i > 1) then
          do j= i-1,1,-1
          if (len_trim(chain1) > len_trim(chaine1(j))) then
                chaine2(j+1)=chaine2(j)
                chaine1(j+1)=chaine1(j)
                chaine2(j)=chain2
                chaine1(j)=chain1
          endif
          enddo
  endif          
enddo
do i=1,nb
  chain1 = trim(adjustl(chaine1(i)))
  chain2 = trim(adjustl(chaine2(i)))
  write(*,*) chain1,chain2
  chain='perl -pi"*.bak" -e "s/>'//trim(chain1)//'>/>'//trim(chain2)//'>/g unless /^#/"'//' fc'
  write(10,*) chain
  chain='perl -pi"*.bak" -e "s/'//trim(chain1)//'>/'//trim(chain2)//'>/g unless /^#/"'//' fc'
  write(10,*) chain
  chain='perl -pi"*.bak" -e "s/>'//trim(chain1)//'/>'//trim(chain2)//'/g unless /^#/"'//' fc'
  write(10,*) chain
  chain='perl -pi"*.bak" -e "s/'//trim(chain1)//'/'//trim(chain2)//'/g unless /^#/"'//' fc'
  write(10,*) chain
!  write(10,*)'perl -pi"*.bak" -e ''s/',chain1,'/',chain2,'/'' config.ini_symb'
enddo

write(10,*) 'cp fc config.ini_conv'

call SYSTEM('chmod 777 txt2nb.sh')
call SYSTEM('./txt2nb.sh')
end
