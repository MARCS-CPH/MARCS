program dmarcs

integer :: idriftok, io, i,n
character :: line*27
real*8  :: delta_t, f_opac, delta_opac, d_t, dt_old
character(len=10) :: file_id
character(len=50) :: file_name, opac_file1, opac_file2 
character(len=50) :: dust_file, chem1_file, chem2_file, chem3_file, sw_log
character(len=50) :: nuclea_file, thermo_file

!set the starting cloud fraction in f_opac. Generally set to 0.0.
f_opac = 0.0
!delta_opac is fractional increment to be added to cloud opacity in every iteration.
delta_opac = 0.05

i = 0

open(unit=987, file='opac_frac.in', status='replace')
write(987,*) f_opac
close(987)
write(file_id, '(i0)') i

do while (f_opac < 1.01)
  write(file_id, '(i0)') i
  print*, "************************** AT " , f_opac*100 , "% OF CLOUD OPACITY ****************"
  if (i>0) then
  call system('cp marcs2drift.dat ./INPUT')
  print *, 'Calling Static Weather.'
  call system('./static_weather16 model_MARCS.in > sw.out')
  call system('rm ./out_default/restart.dat')
  dust_file = './sw_models/dust_out_' // trim(adjustl(file_id)) // '.dat'
  chem1_file = './dust_chem/chem1_out_' // trim(adjustl(file_id)) // '.dat'
  chem2_file = './dust_chem/chem2_out_' // trim(adjustl(file_id)) // '.dat'
  chem3_file = './dust_chem/chem3_out_' // trim(adjustl(file_id)) // '.dat'
  nuclea_file = './dust_nuclea/nuclea_out_' // trim(adjustl(file_id)) // '.dat'
  thermo_file = './dust_thermo/thermo_out_' // trim(adjustl(file_id)) // '.dat'
  call system ('cp ./out_default/out3_dust.dat ' // dust_file)
  call system ('cp ./out_default/out3_chem1.dat ' // chem1_file)
  call system ('cp ./out_default/out3_chem2.dat ' // chem2_file)
  call system ('cp ./out_default/out3_chem3.dat ' // chem3_file)
  call system ('cp ./out_default/out3_nuclea.dat ' // nuclea_file)
  call system ('cp ./out_default/out3_thermo.dat ' // thermo_file)
  call system ('cp ./out_default/drift2marcs.dat ./sw_models/drift2marcs_' // trim(adjustl(file_id)) // '.dat')

  print *, 'Checking Static Weather output.'
  open(unit=2,file='sw.out',readonly)
  idriftok = 0
  do
    read(2,'(a27)',iostat=io) line(1:27)
    if(io < 0) exit
    if(index(line,' regular end of integration')) then
      idriftok = 1
      print*, "Static Weather converged."
      exit
    end if
  end do
  close(2)
  sw_log = './sw_logs/sw_' // trim(adjustl(file_id)) // '.out'
  call system ('cp sw.out ' // sw_log)
  call system ('rm sw.out')
  if(idriftok /= 1) then
      print *, 'Static Weather failed.'
      exit
  end if
  end if
  print *, 'Calling MARCS.'
  if (i.eq.0) then
    call system('cp marcsdrift_nodust.input  mxms7.input')
    call system('cp marcs_init.arciv arcivaaa.dat')
    call system('cp marcs_init.dat old.dat')
  else 
    call system('cp marcsdrift_dust.input  mxms7.input')
    call system('cp arcivaab.dat arcivaaa.dat')
    call system('cp mxmodel.dat old.dat')
  end if
  call system('./marcs')
  call system('rm fort*')
  
  call system('cp mxmodel.dat new.dat')
  file_name = './marcs_models/marcs_model_' // trim(adjustl(file_id)) // '.dat'
  call system('cp mxmodel.dat ' // file_name)
  file_name = './marcs_models/marcs_model_' // trim(adjustl(file_id)) // '.arciv'
  call system('cp marcs2drift.dat ./marcs_models/marcs2drift_' // trim(adjustl(file_id)) //'.dat')
  call system('cp arcivaab.dat ' // file_name)
  call system('cp gas_opac.dat  ./opacities/gas_opac_'  // trim(adjustl(file_id)) //'.dat')
  if (i>0) then
  
  call system('cp dust_opac.dat ./opacities/dust_opac_' // trim(adjustl(file_id)) //'.dat')
  print *, 'Checking model difference.'
  call system('python3 MD_convergence.py build')
  open(unit=1,file='check_conv.txt',readonly)
    read(1,'(f30.1)') d_t
  close(1)
  if (i==1) dt_old=d_t
  print*, "Model difference is ", d_t
  call system('rm check_conv.txt')
  ! if ((d_t> 1.50*dt_old).or. d_t>1.0e-2) then
  !     delta_opac = delta_opac/2.
  !     print*, "delta opac is now ", delta_opac
  ! end if
  dt_old = d_t
  end if
  
  !if (i>200) call abort
  
  call system('rm opac_frac.in')
  call system('rm  ./out_default/drift2marcs.dat')
  ! if (d_t > delta_t) then
  !  if  (d_t > 1.0e-2) delta_opac = delta_opac / 2.
  ! end if
  ! delta_t = d_t
  ! print*, "delta opac is ", delta_opac*100.0 , " %"
  f_opac = f_opac + delta_opac
  open(unit=987, file='opac_frac.in', status='replace')
  write(987,*) f_opac
  close(987)
  i = i+1
  
end do

m = 0
 f_opac = 1.0
 delta_t = 1.0
 call system('rm opac_frac.in')
 open(unit=987, file='opac_frac.in', status='replace')
 write(987,*) f_opac
 close(987)
 print*, "At 100% cloud, check for convergence once again"
 do while (delta_t > (1.0e-2) .and. m < 10)
   write(file_id, '(i0)') i
   call system('cp marcs2drift.dat ./INPUT')
   print *, 'Calling Static Weather.'
   call system('./static_weather16 model_MARCS.in > sw.out')
   call system('rm ./out_default/restart.dat')
   dust_file = './sw_models/dust_out_' // trim(adjustl(file_id)) // '.dat'
   chem1_file = './dust_chem/chem1_out_' // trim(adjustl(file_id)) // '.dat'
   chem2_file = './dust_chem/chem2_out_' // trim(adjustl(file_id)) // '.dat'
   chem3_file = './dust_chem/chem3_out_' // trim(adjustl(file_id)) // '.dat'
   call system ('cp ./out_default/out3_dust.dat ' // dust_file)
   call system ('cp ./out_default/out3_chem1.dat ' // chem1_file)
   call system ('cp ./out_default/out3_chem2.dat ' // chem2_file)
   call system ('cp ./out_default/out3_chem3.dat ' // chem3_file)
   call system('rm ./out_default/out3_dust.dat ')
   call system ('cp ./out_default/drift2marcs.dat ./sw_models/drift2marcs_' // trim(adjustl(file_id)) // '.dat')
   print *, 'Checking Static Weather output.'
   open(unit=2,file='sw.out',readonly)
   idriftok = 0
   do
     read(2,'(a27)',iostat=io) line(1:27)
     if(io < 0) exit
     if(index(line,' regular end of integration')) then
       idriftok = 1
       print*, "Static Weather converged."
       exit
     end if
   end do
   close(2)
   if(idriftok /= 1) then
       print *, 'Static Weather failed.'
       exit
   end if
   print *, 'Calling MARCS.'
  
  if (i==0) then
   call system('cp marcsdrift_nodust.input  mxms7.input')
   call system('cp marcs_init.arciv arcivaaa.dat')
   call system('cp marcs_init.dat old.dat')

   else
     call system('cp arcivaab.dat arcivaaa.dat')
     call system('cp mxmodel.dat old.dat')
  end if

   call system('./marcs')
   call system('rm fort*')

  
   call system('cp mxmodel.dat new.dat')
   file_name = './marcs_models/marcs_model_' // trim(adjustl(file_id)) // '.dat'
   call system('cp mxmodel.dat ' // file_name)
   file_name = './marcs_models/marcs_model_' // trim(adjustl(file_id)) // '.arciv'
   call system('cp marcs2drift.dat ./marcs_models/marcs2drift_' // trim(adjustl(file_id)) //'.dat')
   call system('cp arcivaab.dat ' // file_name)
   call system('cp gas_opac.dat  ./opacities/gas_opac_'  // trim(adjustl(file_id)) //'.dat')
   call system('cp dust_opac.dat ./opacities/dust_opac_' // trim(adjustl(file_id)) //'.dat')

   print *, 'Checking MARCS-SW convergence.'
   call system('python3 MD_convergence.py build')
   open(unit=1,file='check_conv.txt',readonly)
     read(1,'(f30.1)') delta_t
   close(1)
   call system('rm check_conv.txt')
   if ((delta_t <= 1.0e-2)) then 
     print*, "MARCS-SW converged, delta T of ", delta_t
     call system('cp arcivaab.dat marcs_init.arciv')
     call system('cp mxmodel.dat marcs_init.dat')
     exit 
   else 
     print*, "MARCS-SW did not converge yet..."
     print*, "Maximum T difference was: ", delta_t
   end if
  
   i = i+1
   m = m+1
end do 
 call system('cp arcivaab.dat drift_testp.arciv')
 call system('cp mxmodel.dat drift_testp.dat')

end
