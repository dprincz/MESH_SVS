!-------------------------------------- LICENCE BEGIN ------------------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------------------

subroutine inisoili_svs2(ni, trnch)
   use sfcbus_mod
   use svs_configs
   use sfc_options
   implicit none
!!!#include <arch_specific.hf>

   integer ni, trnch

   !@Author  Maria Abrahamowicz, Stephane Belair , Vincent Fortin (20xx)
   !@Object  Compute soil properties for given soil texture. Compute these properties on 
   !         native levels of database providing soil texture, and then map properties 
   !         unto SVS levels.
   !@Arguments
   !             - Input -
   ! NI          longueur d'une tranche horizontale

   integer :: i, k, kk
   REAL b, usb, fb, crit1_wfcint, crit2_wfcint, soillayer_depth, beta_soc
   REAL :: rho_silt = 2692., & ! Density of silt (kg/m3)
           rho_sand = 2656., & ! Density of sand (kg/m3)
           rho_clay = 2761.    ! Density of clay (kg/m3)
   
   ! "geo" variables are on the levels of the geophysical soil texture datbase
   REAL, dimension(ni,nl_stp) :: wsat_geo, wwilt_geo, wfc_geo, b_geo, psisat_geo, &
           ksat_geo, wfcint_geo, fb_geo ,quartz_geo,rhosoil_geo,conddry_geo,condsld_geo, &
           soilhcapz_dry_geo, fvom_geo, rbom_geo
   REAL :: fmom_geo,   rho_sms_geo, rho_bms_geo, rho_b_geo, rho_bom_geo

   ! 100% soil organic content variables 
   REAL, dimension(ni,nl_svs) :: rbom_soc
   REAL :: wsat_soc, wwilt_soc, wfc_soc, b_soc, psisat_soc, &
           ksat_soc, wfcint_soc, fb_soc ,conddry_soc,condsld_soc, hcap_soc, fs_vom
   real, pointer, dimension(:) :: zcgsat, zgrkef, zdraindens, zslop

   ! variables on the levels of SVS
   real, pointer, dimension(:,:) :: zbcoef, zclay, zfbcof, zksat, zpsisat, zsand, zwfc, zwfcint, zwsat, zwwilt, zfvom 
   real, pointer, dimension(:,:) :: zconddry, zcondsld, zquartz, zrhosoil, zsoilhcapz_dry, zgravel, zbulksoil, zoc, zsilt


  
#define MKPTR1D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni) => busptr(vd%NAME2%i)%ptr(:,trnch)
#define MKPTR2D(NAME1,NAME2) nullify(NAME1); if (vd%NAME2%i > 0 .and. associated(busptr(vd%NAME2%i)%ptr)) NAME1(1:ni,1:vd%NAME2%mul*vd%NAME2%niveaux) => busptr(vd%NAME2%i)%ptr(:,trnch)

   MKPTR1D(zcgsat, cgsat)
   MKPTR1D(zdraindens, draindens)
   MKPTR1D(zgrkef, grkef)
   MKPTR1D(zslop, slop)


   MKPTR2D(zbcoef, bcoef)
   MKPTR2D(zbulksoil , bulksoil)
   MKPTR2D(zclay, clay)
   MKPTR2D(zconddry, conddry)
   MKPTR2D(zcondsld, condsld)  
   MKPTR2D(zfbcof, fbcof)
   MKPTR2D(zfvom , fvom) ! volumetric fraction of soil organic content
   MKPTR2D(zgravel , gravel)
   MKPTR2D(zksat, ksat)
   MKPTR2D(zoc , oc)
   MKPTR2D(zpsisat , psisat)
   MKPTR2D(zquartz, quartz)  
   MKPTR2D(zrhosoil, rhosoil)  
   MKPTR2D(zsand, sand)
   MKPTR2D(zsilt, silt)
   MKPTR2D(zwfc, wfc)
   MKPTR2D(zwfcint, wfcint)
   MKPTR2D(zwsat, wsat)
   MKPTR2D(zwwilt , wwilt)
   MKPTR2D(zsoilhcapz_dry , soilhcapz_dry)


   !call subroutine to compute layer thicknesses
   call layer_thickness()

   ! calculate soil parameters on native GEO layers, and then map them unto model layers. 
   ! calculate weights to be used in phybusinit.... because here... we are
   ! re-doing the calculation for each row of domain...
   ! but the weights are the same !

   ! Compute averaged values for sand and clay for each layer when SOILGRID is used

   if ( soiltext == "SOILGRIDS" ) then   
      do i=1,ni
         do k=1,(nl_stp-1)
            zsand(i,k) = (zsand(i,k)+zsand(i,k+1))/2.
            zclay(i,k) = (zclay(i,k)+zclay(i,k+1))/2.
         enddo
      enddo
   endif



   !     Computer soil properties for GEO layers
   do i=1,ni
      do k=1,nl_stp

         wsat_geo  (i,k)  =  -0.00126   * zsand(i,k) + 0.489
         wwilt_geo (i,k)  =  37.1342e-3 * sqrt(max(1.,zclay(i,k)))
         wfc_geo   (i,k)  =  89.0467e-3 * max(1.,zclay(i,k))**0.3496
         psisat_geo(i,k)  =  0.01 * ( 10.0**(-0.0131 * zsand(i,k) + 1.88) )
         ksat_geo  (i,k)  =  ( 10.0**(0.0153 * zsand(i,k) - 0.884) ) * 7.0556E-6

         b                 =  0.137 * zclay(i,k)  + 3.501
         b_geo     (i,k)  =  b
         usb               =  1./b
         fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
         fb_geo(i,k)      =  fb

         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation

         crit1_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*(wsat_geo(i,k)/wwilt_geo(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*psisat_geo(i,k)*fb**b

         if (abs(zslop(i)).gt.crit1_wfcint) then
            wfcint_geo(i,k) = wwilt_geo(i,k)        
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            wfcint_geo(i,k) = wsat_geo(i,k) 
         elseif (zslop(i).ne.0.0) then
            wfcint_geo(i,k) = wsat_geo(i,k) * fb * &
                 ( psisat_geo(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            wfcint_geo(i,k) = wfc_geo(i,k)
         endif
         
         ! Compute soil thermal properties for heat diffusion

!       Quartz content (ref : NL95 & PL98)): 
        quartz_geo(i,k)  = 0.038 + 0.0095*zsand(i,k)

!       Soil dry density (PL98):
        rhosoil_geo(i,k) = (1.0-wsat_geo(i,k))*2700.

!       Soil solid conductivity:
        if (quartz_geo(i,k).gt.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (2.0**(1.0-quartz_geo(i,k))) 
        endif
        if (quartz_geo(i,k).le.0.20) then
           condsld_geo(i,k) = (7.7**quartz_geo(i,k)) *  &
                           (3.0**(1.0-quartz_geo(i,k)))
        endif

!       Soil dry conductivity:
        conddry_geo(i,k) = (0.135*rhosoil_geo(i,k) + 64.7) / &
                        (2700. - 0.947*rhosoil_geo(i,k))   

!       Soil heat capacity
        if ((zsand(i,k)+zclay(i,k)) .gt. epsilon_svs) then
            ! From Lawrence and Slater 2008, 
            soilhcapz_dry_geo(i,k) = (2.128*zsand(i,k) + 2.385*zclay(i,k))/(zsand(i,k)+zclay(i,k))*1.E6
        else    
            ! Previous formulation used in SVS - constant value
            soilhcapz_dry_geo(i,k) = 2700. * 733.
        endif


      enddo
   enddo

   !     Computer fraction of soil organic matter for GEO layers
   if (read_oc) then

      ! Initialize rbom_soc calculated on SVS levels
      rbom_soc = 0.

      do i=1,ni
         do k=1,nl_stp
            ! Theory and variables explained in Decharme (2025) 
            ! zoc(i,k)*1E-2 is SOC mass fraction [kg/kg]
            ! zclay(i,k)*1E-2 is clay mass fraction [kg/kg], same for sand, and sild
            ! zbulksoil(i,k) is bulk density of the fine earth fraction [kg/m3]
            ! zgravel(i,k)*1E-2 is volumetric fraction of coarse fragments [m3/m3]


            if (zbulksoil(i,k) .gt. epsilon_svs .and. zoc(i,k) .gt. epsilon_svs) then 

               fmom_geo = 1.848 * (zoc(i,k)*1E-2)**(0.967) 
               fmom_geo = max(min(fmom_geo, 1.),0.)

               rho_sms_geo = 1./(zclay(i,k)*1E-2 / rho_clay + zsilt(i,k)*1E-2 / rho_silt + zsand(i,k)*1E-2 / rho_sand)
               rho_bms_geo = (1. - wsat_geo(i,k)) * rho_sms_geo 
               rho_b_geo   = zbulksoil(i,k) * (1.- zgravel(i,k)*1E-2)

               if (1. / rho_b_geo - (1.-fmom_geo)/rho_bms_geo .gt. epsilon_svs) then
                  rho_bom_geo = fmom_geo / (1. / rho_b_geo - (1.-fmom_geo)/rho_bms_geo)

                  fvom_geo(i,k) =  fmom_geo * rho_b_geo / rho_bom_geo
                  fvom_geo(i,k) = max(min(fvom_geo(i,k), 1.),0.)

                  rbom_geo(i,k) = rho_bom_geo *1E-3
                  rbom_geo(i,k) = max(min(rbom_geo(i,k), 1.),0.)

               else
                  fvom_geo(i,k) = 0.
                  rbom_geo(i,k) = 0.
               endif

            else
               rbom_geo(i,k) = 0.
               fvom_geo(i,k) = 0.
            endif

         enddo
      enddo
   else
      ! Force variables for soil organic to 0
      rbom_geo = 0.
      rbom_soc = 0.
      fvom_geo = 0.
   endif

   ! "Map" GEO soil properties unto model soil layers
   Do i = 1 , ni
      Do k = 1, nl_svs
         do kk = 1 , nl_stp
            
            zwsat  (i,k)  = zwsat  (i,k) + wsat_geo  (i,kk)  * weights( k , kk)
            zwwilt (i,k)  = zwwilt (i,k) + wwilt_geo (i,kk)  * weights( k , kk)
            
            zwfc   (i,k)  = zwfc   (i,k) + wfc_geo   (i,kk)  * weights( k , kk)
            zbcoef (i,k)  = zbcoef (i,k) + b_geo     (i,kk)  * weights( k , kk)
            zfbcof (i,k)  = zfbcof (i,k) + fb_geo    (i,kk)  * weights( k , kk)
            zpsisat(i,k)  = zpsisat(i,k) + psisat_geo(i,kk)  * weights( k , kk)
            zksat  (i,k)  = zksat  (i,k) + ksat_geo  (i,kk)  * weights( k , kk)
            zwfcint(i,k)  = zwfcint(i,k) + wfcint_geo(i,kk)  * weights( k , kk)
            zconddry  (i,k)  = zconddry  (i,k) + conddry_geo  (i,kk)  * weights( k , kk)
            zcondsld  (i,k)  = zcondsld  (i,k) + condsld_geo  (i,kk)  * weights( k , kk)
            zquartz   (i,k)  = zquartz   (i,k) + quartz_geo   (i,kk)  * weights( k , kk)
            zrhosoil  (i,k)  = zrhosoil  (i,k) + rhosoil_geo  (i,kk)  * weights( k , kk)
            zsoilhcapz_dry  (i,k)  = zsoilhcapz_dry  (i,k) + soilhcapz_dry_geo  (i,kk)  * weights( k , kk)

            ! Soil organic matter variables
            zfvom  (i,k)  = zfvom  (i,k) + fvom_geo  (i,kk)  * weights( k , kk)
            rbom_soc  (i,k)  = rbom_soc  (i,k) + rbom_geo  (i,kk)  * weights( k , kk)

         enddo
      enddo
      ! compute thermal coeff. 
      ! for 1st model layer only --- here simply use 1st GEO soil texture !!! Do not map !
      zcgsat (i)  = ( -1.557e-2 * zsand(i,1) &
           -  1.441e-2 * zclay(i,1) + 4.7021 ) * 1.E-6 
      
      ! Compute effective parameter for watdrain
      zgrkef(i)   = 2.* zdraindens(i) * zslop(i)

   enddo

   !     Computer 100% soil organic matter properties (Decharme 2025)
   if (read_oc) then
      
      do i=1,ni
         do k=1,nl_svs

            if (k .EQ. 1) then
               soillayer_depth = delz(k) * 0.5
            else
               soillayer_depth = soillayer_depth + (delz(k) + delz(k-1))*0.5
            endif


            wsat_soc    =  0.95 - 0.437 * rbom_soc(i,k)
            b_soc       =  2.933 + 0.442 * rbom_soc(i,k)**(0.463) + exp(1.321 * rbom_soc(i,k))
            psisat_soc  = abs((101.633 * rbom_soc(i,k)**4 - 46.913 * rbom_soc(i,k)**5 - 61.625 * rbom_soc(i,k)**(2.635)) * 0.0168**(rbom_soc(i,k)))
            ksat_soc    = 10**(-7.955 - 1.89 * LOG10( min(soillayer_depth,3.) + 0.068) - 2.96 * LOG10( min(0.25, rbom_soc(i,k)) + 0.045)) 

            wfc_soc     = 3.1486 * 0.12**rbom_soc(i,k) * rbom_soc(i,k)**0.70
            wwilt_soc   = 0.9355 * 0.20**rbom_soc(i,k) * rbom_soc(i,k)**0.71

            hcap_soc    = rbom_soc(i,k)*1E3 / (1. - wsat_soc) * 1972.
            condsld_soc = 0.25 
            conddry_soc = 0.05   

            fs_vom =  zfvom(i,k) * (1. - wsat_soc) / (1. - zwsat  (i,k))


            ! Calculating soil properties accounting for soil organic matter

           ! Arithmetic mean
            zwsat  (i,k)  = (1.-zfvom(i,k)) * zwsat  (i,k) + zfvom(i,k) * wsat_soc 
            zwwilt (i,k)  = (1.-zfvom(i,k)) * zwwilt  (i,k) + zfvom(i,k) * wwilt_soc 
            zwfc   (i,k)  = (1.-zfvom(i,k)) * zwfc  (i,k) + zfvom(i,k) * wfc_soc 
            zpsisat(i,k)  = (1.-zfvom(i,k)) * zpsisat  (i,k) + zfvom(i,k) * psisat_soc  
            zsoilhcapz_dry (i,k) = (1.-fs_vom) * zsoilhcapz_dry  (i,k) +fs_vom * hcap_soc 

            ! Geometric mean
            zcondsld(i,k) = zcondsld(i,k)**(1.-zfvom(i,k)) * condsld_soc**zfvom(i,k)
            zksat  (i,k)  = zksat(i,k)**(1.-zfvom(i,k)) * ksat_soc**zfvom(i,k)
            zconddry  (i,k)  = ( (zconddry(i,k) * conddry_soc)**0.5  / (conddry_soc**0.5 *(1.-zfvom(i,k)) + zconddry(i,k)**0.5 * zfvom(i,k)) )**2

            b                 = (1.-zfvom(i,k)) * zbcoef  (i,k) + zfvom(i,k) * b_soc 
            zbcoef     (i,k)  =  b
            usb               =  1./b
            fb                =  b**usb/(b-1.) * ((3.*b+2.)**(1.-usb)-(2.*b+2.)**(1.-usb))
            zfbcof(i,k)      =  fb

         enddo
      enddo

   endif

   Do i = 1 , ni
      Do k = 1, nl_svs

         ! Compute water content at field capacity along sloping aquifer based on Soulis et al. 2012
         ! Ensure that wc at fc stays between wilting point and saturation

         crit1_wfcint   = 2.*zdraindens(i)*zpsisat(i,k)*(zwsat(i,k)/zwwilt(i,k)*fb)**b
         crit2_wfcint   = 2.*zdraindens(i)*zpsisat(i,k)*fb**b

         if (abs(zslop(i)).gt.crit1_wfcint) then
            zwfcint(i,k) = zwwilt(i,k)        
         elseif (abs(zslop(i)).lt.crit2_wfcint) then
            zwfcint(i,k) = zwsat(i,k) 
         elseif (zslop(i).ne.0.0) then
            zwfcint(i,k) = zwsat(i,k) * fb * &
                 ( zpsisat(i,k)/ABS(zslop(i)) *2. * zdraindens(i) )**usb
         else
            zwfcint(i,k) = zwfc(i,k)
         endif

      enddo

   enddo

   return
 end subroutine inisoili_svs2
