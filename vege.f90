MODULE VEGE
 
USE PRECISION_PARAMETERS
USE GLOBAL_CONSTANTS
USE MESH_POINTERS
USE TRAN
USE PART
USE MEMORY_FUNCTIONS, ONLY:CHKMEMERR
USE TYPES, ONLY: PARTICLE_TYPE, PARTICLE_CLASS_TYPE, PARTICLE_CLASS! WALL_TYPE,SURFACE_TYPE 
IMPLICIT NONE
PRIVATE
PUBLIC INITIALIZE_RAISED_VEG, DEALLOCATE_VEG_ARRAYS, RAISED_VEG_MASS_ENERGY_TRANSFER, LEVEL_SET_FIRESPREAD, &
       GET_REV_vege, BNDRY_VEG_MASS_ENERGY_TRANSFER
TYPE (PARTICLE_TYPE), POINTER :: LP=>NULL()
TYPE (PARTICLE_CLASS_TYPE), POINTER :: PC=>NULL()
!TYPE (WALL_TYPE), POINTER :: WC
!TYPE (SURFACE_TYPE), POINTER :: SF 
CHARACTER(255), PARAMETER :: vegeid='$Id: vege.f90 9718 2011-12-30 17:49:06Z drjfloyd $'
CHARACTER(255), PARAMETER :: vegerev='$Revision: 9718 $'
CHARACTER(255), PARAMETER :: vegedate='$Date: 2011-12-30 09:49:06 -0800 (Fri, 30 Dec 2011) $'
LOGICAL, ALLOCATABLE, DIMENSION(:,:,:) :: VEG_PRESENT_FLAG,CELL_TAKEN_FLAG
INTEGER :: IZERO,NLP_VEG_FUEL,NCONE_TREE,NXB,NYB
REAL(EB) :: RCELL,R_TREE,XCELL,XI,YJ,YCELL,ZCELL,ZK,DYN_SR_MAX
 
CONTAINS
 

SUBROUTINE INITIALIZE_RAISED_VEG(NM)

USE MEMORY_FUNCTIONS, ONLY: RE_ALLOCATE_PARTICLES
USE COMP_FUNCTIONS, ONLY: GET_FILE_NUMBER
USE TRAN, ONLY : GET_IJK
REAL(EB) CROWN_LENGTH,CROWN_VOLUME,TANGENT,CROWN_WIDTH
REAL(EB) DX_RING,DZ_RING,INNER_RADIUS,OUTER_RADIUS,R_CTR_CYL,  &
         RING_BOTTOM,RING_TOP,SLANT_WIDTH,TFLAG
REAL(EB) V_CELL,XLOC,YLOC,ZLOC,X_EXTENT,Y_EXTENT,Z_EXTENT
INTEGER NCT,NLP_TREE,NLP_RECT_VEG,N_TREE,NXB,NYB,NZB,IPC,IJK_CT(6)
INTEGER N_CFCR_TREE,N_FRUSTUM_TREE,N_RECT_TREE,N_CUST_TREE,N_RING_TREE,N_IGN
INTEGER I,II,I_OUTER_RING,JJ,KK,K_BOTTOM_RING,LU_CBD
INTEGER, INTENT(IN) :: NM
CHARACTER(30) :: CBDFILE,NMFILE,NTFILE
LOGICAL :: EX

!IF (.NOT. TREE) RETURN !Exit if there are no trees anywhere
IF (.NOT. TREE_MESH(NM)) RETURN !Exit routine if no raised veg in mesh
IF (EVACUATION_ONLY(NM)) RETURN  ! Don't waste time if an evac mesh
CALL POINT_TO_MESH(NM)

ALLOCATE(VEG_PRESENT_FLAG(0:IBP1,0:JBP1,0:KBP1))
CALL ChkMemErr('VEGE','VEG_PRESENT_FLAG',IZERO)
ALLOCATE(CELL_TAKEN_FLAG(0:IBP1,0:JBP1,0:KBP1))
CALL ChkMemErr('VEGE','CELL_TAKEN_FLAG',IZERO)

!Diagnostic files
!IF (NM == NMESHES) THEN
!OPEN(9999,FILE='total_PARTICLE_mass.out',STATUS='REPLACE')
! OPEN(9998,FILE='diagnostics.out',STATUS='REPLACE')
!ENDIF

TREE_MESH(NM)          = .FALSE. 
CONE_TREE_PRESENT      = .FALSE.
FRUSTUM_TREE_PRESENT   = .FALSE.
CYLINDER_TREE_PRESENT  = .FALSE.
RING_TREE_PRESENT      = .FALSE.
RECTANGLE_TREE_PRESENT = .FALSE.
CUSTOM_TREE_PRESENT    = .FALSE.

TREE_LOOP: DO NCT=1,N_TREES

   VEG_PRESENT_FLAG = .FALSE. ; CELL_TAKEN_FLAG = .FALSE.
   IPC = TREE_PARTICLE_CLASS(NCT)
   PC=>PARTICLE_CLASS(IPC)
   PC%KILL_RADIUS = 0.5_EB/PC%VEG_SV !radius bound below which fuel elements are removed
! 
! Build a conical volume of solid (vegetation) fuel
!
   IF_CONE_VEGETATION: IF(VEG_FUEL_GEOM(NCT) == 'CONE') THEN
!
   CONE_TREE_PRESENT = .TRUE.
   N_CFCR_TREE = TREE_CFCR_INDEX(NCT)
   CROWN_WIDTH  = CROWN_W(N_CFCR_TREE)
   CROWN_LENGTH = TREE_H(N_CFCR_TREE) - CROWN_B_H(N_CFCR_TREE)
   IF(CROWN_LENGTH <= 0.0_EB) THEN
     PRINT*,'ERROR CONE TREE: Crown base height >= tree height for (maybe) tree number ', NCT
     PRINT*,'CONE TREE_HEIGHT = ',TREE_H(N_CFCR_TREE)
     PRINT*,'CONE CROWN_BASE_HEIGHT = ',CROWN_B_H(N_CFCR_TREE)
     STOP
   ENDIF
   TANGENT = 0.5_EB*CROWN_W(N_CFCR_TREE)/CROWN_LENGTH
   CROWN_VOLUME = PI*CROWN_WIDTH**2*CROWN_LENGTH/12._EB
 
   NLP_TREE = 0

   DO NZB=1,KBAR
     IF (Z(NZB)>=Z_TREE(N_CFCR_TREE)+CROWN_B_H(N_CFCR_TREE) .AND. & 
         Z(NZB)<=Z_TREE(N_CFCR_TREE)+TREE_H(N_CFCR_TREE)) THEN
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
!      R_TREE = TANGENT*(TREE_H(N_CFCR_TREE)+Z_TREE(N_CFCR_TREE)-Z(NZB)+0.5_EB*DZ(NZB))
      R_TREE = TANGENT*(TREE_H(N_CFCR_TREE)+Z_TREE(N_CFCR_TREE)-Z(NZB))
      DO NXB = 1,IBAR
       DO NYB = 1,JBAR
        RCELL = SQRT((X(NXB)-X_TREE(N_CFCR_TREE))**2 + (Y(NYB)-Y_TREE(N_CFCR_TREE))**2)
        IF (RCELL <= R_TREE) THEN
         NLP  = NLP + 1
         NLP_TREE = NLP_TREE + 1
         IF (NLP>NLPDIM) THEN
          CALL RE_ALLOCATE_PARTICLES(1,NM,0,1000)
          PARTICLE=>MESHES(NM)%PARTICLE
         ENDIF
         LP=>PARTICLE(NLP)
         LP%VEG_VOLFRACTION = 1._EB
         LP%TAG = PARTICLE_TAG
         LP%X = REAL(NXB,EB)
         LP%Y = REAL(NYB,EB)
         LP%Z = REAL(NZB,EB)
         LP%CLASS = IPC
         LP%PWT   = 1._EB  ! This is not used, but it is necessary to assign a non-zero weight factor to each particle
         VEG_PRESENT_FLAG(NXB,NYB,NZB) = .TRUE.
        ENDIF
       ENDDO   
      ENDDO 
     ENDIF
   ENDDO
   NLP_VEG_FUEL = NLP_TREE
!
   ENDIF IF_CONE_VEGETATION
! 
! Build a frustum volume of solid (vegetation) fuel
!
   IF_FRUSTUM_VEGETATION: IF(VEG_FUEL_GEOM(NCT) == 'FRUSTUM') THEN
!
   FRUSTUM_TREE_PRESENT = .TRUE.
   N_CFCR_TREE          = TREE_CFCR_INDEX(NCT)
   N_FRUSTUM_TREE       = TREE_FRUSTUM_INDEX(NCT)
   CROWN_LENGTH         = TREE_H(N_CFCR_TREE) - CROWN_B_H(N_CFCR_TREE)
   IF(CROWN_LENGTH <= 0.0_EB) THEN
     PRINT*,'ERROR FRUSTUM TREE: Crown base height >= tree height for (maybe) tree number ', NCT
     PRINT*,'FRUSTUM TREE_HEIGHT = ',TREE_H(N_CFCR_TREE)
     PRINT*,'FRUSTUM CROWN_BASE_HEIGHT = ',CROWN_B_H(N_CFCR_TREE)
     STOP
   ENDIF
   R_CTR_CYL    = 0.5*MIN(CROWN_W_TOP(N_FRUSTUM_TREE),CROWN_W_BOTTOM(N_FRUSTUM_TREE))
   SLANT_WIDTH  = 0.5*ABS(CROWN_W_TOP(N_FRUSTUM_TREE) - CROWN_W_BOTTOM(N_FRUSTUM_TREE))

   TANGENT = SLANT_WIDTH/CROWN_LENGTH
   CROWN_VOLUME = PI*CROWN_LENGTH*(CROWN_W_BOTTOM(N_FRUSTUM_TREE)**2 + & 
                  CROWN_W_TOP(N_FRUSTUM_TREE)*CROWN_W_TOP(N_FRUSTUM_TREE) + CROWN_W_TOP(N_FRUSTUM_TREE)**2)/3._EB
 
   NLP_TREE = 0

   DO NZB=1,KBAR
     IF (Z(NZB)>=Z_TREE(N_CFCR_TREE)+CROWN_B_H(N_CFCR_TREE) .AND. & 
                 Z(NZB)<=Z_TREE(N_CFCR_TREE)+TREE_H(N_CFCR_TREE)) THEN
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      IF(CROWN_W_TOP(N_FRUSTUM_TREE) <= CROWN_W_BOTTOM(N_FRUSTUM_TREE)) & 
                 R_TREE = R_CTR_CYL + TANGENT*(TREE_H(N_CFCR_TREE)+Z_TREE(N_CFCR_TREE)-Z(NZB))
      IF(CROWN_W_TOP(N_FRUSTUM_TREE) >  CROWN_W_BOTTOM(N_FRUSTUM_TREE)) &
                 R_TREE = R_CTR_CYL + TANGENT*(Z(NZB)-Z_TREE(N_CFCR_TREE)-CROWN_B_H(N_CFCR_TREE))
      DO NXB = 1,IBAR
       DO NYB = 1,JBAR
        RCELL = SQRT((X(NXB)-X_TREE(N_CFCR_TREE))**2 + (Y(NYB)-Y_TREE(N_CFCR_TREE))**2)
        IF (RCELL <= R_TREE) THEN
         NLP  = NLP + 1
         NLP_TREE = NLP_TREE + 1
         IF (NLP>NLPDIM) THEN
          CALL RE_ALLOCATE_PARTICLES(1,NM,0,1000)
          PARTICLE=>MESHES(NM)%PARTICLE
         ENDIF
         LP=>PARTICLE(NLP)
         LP%VEG_VOLFRACTION = 1._EB
         LP%TAG = PARTICLE_TAG
         LP%X = REAL(NXB,EB)
         LP%Y = REAL(NYB,EB)
         LP%Z = REAL(NZB,EB)
         LP%CLASS = IPC
         LP%PWT   = 1._EB  ! This is not used, but it is necessary to assign a non-zero weight factor to each particle
         VEG_PRESENT_FLAG(NXB,NYB,NZB) = .TRUE.
        ENDIF
       ENDDO   
      ENDDO 
     ENDIF
   ENDDO
   NLP_VEG_FUEL = NLP_TREE
!
   ENDIF IF_FRUSTUM_VEGETATION
!
! Build a cylindrical volume of vegetative fuel
!
   IF_CYLINDRICAL_VEGETATION: IF (VEG_FUEL_GEOM(NCT) == 'CYLINDER') THEN
!
   CYLINDER_TREE_PRESENT = .TRUE.
   N_CFCR_TREE           = TREE_CFCR_INDEX(NCT)
   CROWN_WIDTH           = CROWN_W(N_CFCR_TREE)
   R_TREE                = 0.5*CROWN_WIDTH
   CROWN_LENGTH          = TREE_H(N_CFCR_TREE) - CROWN_B_H(N_CFCR_TREE)
   IF(CROWN_LENGTH <= 0.0_EB) THEN
     PRINT*,'ERROR CYLINDER TREE: Crown base height >= tree height for (maybe) tree number ', NCT
     PRINT*,'CYLINDER TREE_HEIGHT = ',TREE_H(N_CFCR_TREE)
     PRINT*,'CYLINDER CROWN_BASE_HEIGHT = ',CROWN_B_H(N_CFCR_TREE)
     STOP
   ENDIF
   CROWN_VOLUME = 0.25*PI*CROWN_WIDTH**2*CROWN_LENGTH
   NLP_TREE = 0

   DO NZB=1,KBAR
     IF (Z(NZB)>=Z_TREE(N_CFCR_TREE)+CROWN_B_H(N_CFCR_TREE) .AND. Z(NZB)<=Z_TREE(N_CFCR_TREE)+TREE_H(N_CFCR_TREE)) THEN
      PARTICLE_TAG = PARTICLE_TAG + NMESHES
      DO NXB = 1,IBAR
       DO NYB = 1,JBAR
        RCELL = SQRT((X(NXB)-X_TREE(N_CFCR_TREE))**2 + (Y(NYB)-Y_TREE(N_CFCR_TREE))**2)
        IF (RCELL <= R_TREE) THEN
         NLP  = NLP + 1
         NLP_TREE = NLP_TREE + 1
         IF (NLP>NLPDIM) THEN
          CALL RE_ALLOCATE_PARTICLES(1,NM,0,1000)
          PARTICLE=>MESHES(NM)%PARTICLE
         ENDIF
         LP=>PARTICLE(NLP)
         LP%VEG_VOLFRACTION = 1._EB
         LP%TAG = PARTICLE_TAG
         LP%X = REAL(NXB,EB)
         LP%Y = REAL(NYB,EB)
         LP%Z = REAL(NZB,EB)
         LP%CLASS = IPC
         LP%PWT   = 1._EB  ! This is not used, but it is necessary to assign a non-zero weight factor to each particle
         VEG_PRESENT_FLAG(NXB,NYB,NZB) = .TRUE.
        ENDIF
       ENDDO   
      ENDDO 
     ENDIF
   ENDDO
   NLP_VEG_FUEL = NLP_TREE
!
   ENDIF IF_CYLINDRICAL_VEGETATION
!
! Build a rectangular volume containing vegetation
!
   IF_RECTANGULAR_VEGETATION:IF (VEG_FUEL_GEOM(NCT) == 'RECTANGLE')THEN
       RECTANGLE_TREE_PRESENT = .TRUE.
       N_RECT_TREE            = TREE_RECT_INDEX(NCT)
       NLP_RECT_VEG           = 0
       DO NZB=0,KBAR-1
        ZLOC = Z(NZB) + 0.5_EB*DZ(NZB)
        IF (ZLOC>=ZS_RECT_VEG(N_RECT_TREE) .AND. ZLOC<=ZF_RECT_VEG(N_RECT_TREE)) THEN
         DO NXB = 0,IBAR-1
          XLOC = X(NXB) + 0.5_EB*DX(NXB)
          IF (XLOC >= XS_RECT_VEG(N_RECT_TREE) .AND. XLOC <= XF_RECT_VEG(N_RECT_TREE)) THEN
           DO NYB = 0,JBAR-1
            YLOC = Y(NYB) + 0.5_EB*DY(NYB)
            IF (YLOC >= YS_RECT_VEG(N_RECT_TREE) .AND. YLOC <= YF_RECT_VEG(N_RECT_TREE)) THEN
             NLP  = NLP + 1
             NLP_RECT_VEG = NLP_RECT_VEG + 1
             IF (NLP>NLPDIM) THEN
              CALL RE_ALLOCATE_PARTICLES(1,NM,0,1000)
              PARTICLE=>MESHES(NM)%PARTICLE
             ENDIF
             LP=>PARTICLE(NLP)
             LP%TAG = PARTICLE_TAG
             LP%X = REAL(NXB,EB)
             LP%Y = REAL(NYB,EB)
             LP%Z = REAL(NZB,EB)
             LP%CLASS = IPC
             LP%PWT   = 1._EB  ! This is not used, but it is necessary to assign a non-zero weight factor to each particle
             VEG_PRESENT_FLAG(NXB,NYB,NZB) = .TRUE.
             X_EXTENT = XF_RECT_VEG(N_RECT_TREE) - XS_RECT_VEG(N_RECT_TREE)
             Y_EXTENT = YF_RECT_VEG(N_RECT_TREE) - YS_RECT_VEG(N_RECT_TREE)
             Z_EXTENT = ZF_RECT_VEG(N_RECT_TREE) - ZS_RECT_VEG(N_RECT_TREE)
             IF(X_EXTENT <= 0.0_EB .OR. Y_EXTENT <= 0.0_EB .OR. Z_EXTENT <= 0.0_EB) THEN
               PRINT*,'ERROR RECTANGULAR TREE: for (maybe) tree number ', NCT
               PRINT*,'ZERO OR NEGATIVE TREE WIDTH IN ONE OR MORE DIRECTIONS'
               PRINT*,'X LENGTH = ',X_EXTENT
               PRINT*,'Y LENGTH = ',Y_EXTENT
               PRINT*,'Z LENGTH = ',Z_EXTENT
               STOP
             ENDIF
             LP%VEG_VOLFRACTION = 1._EB
!            IF (X_EXTENT < DX(NXB)) LP%VEG_VOLFRACTION = LP%VEG_VOLFRACTION*X_EXTENT/DX(NXB)
!            IF (Y_EXTENT < DY(NYB)) LP%VEG_VOLFRACTION = LP%VEG_VOLFRACTION*Y_EXTENT/DY(NYB)
             IF (Z_EXTENT < DZ(NZB)) LP%VEG_VOLFRACTION = LP%VEG_VOLFRACTION*Z_EXTENT/DZ(NZB)
!            print*,'veg_volfraction',z_extent,dz(nzb),LP%veg_volfraction
!            print*,'veg_volfraction',xs_rect_veg(nct),xf_rect_veg(nct),ys_rect_veg(nct),yf_rect_veg(nct), &
!                                     zs_rect_veg(nct),zf_rect_veg(nct),z_extent,dz(nzb),LP%VEG_VOLFRACTION
            ENDIF
           ENDDO   
          ENDIF
         ENDDO 
        ENDIF
       ENDDO
       NLP_VEG_FUEL = NLP_RECT_VEG
   ENDIF IF_RECTANGULAR_VEGETATION  
   
!
! Build a custom volume containing vegetation
!
   IF_CUSTOM_VEGETATION:IF (VEG_FUEL_GEOM(NCT) == 'CUSTOM')THEN
       WRITE(NMFILE,*) NM
       WRITE(NTFILE,*) NCT
       CBDFILE = 'CBD_'//TRIM(ADJUSTL(NTFILE))//'_'//TRIM(ADJUSTL(NMFILE))//'.txt'
       INQUIRE(FILE=CBDFILE,EXIST=EX)
       IF (EX) THEN 
           LU_CBD=GET_FILE_NUMBER()
           CUSTOM_TREE_PRESENT = .TRUE.
           N_CUST_TREE            = TREE_CUSTOM_INDEX(NCT)
           NLP_TREE               = 0
           OPEN(LU_CBD,FILE=CBDFILE,FORM='FORMATTED')
           READ(LU_CBD,*) IJK_CT
           !diagnostic check 
           !print*,IJK_CT(1),IJK_CT(2),IJK_CT(3),IJK_CT(4),IJK_CT(5),IJK_CT(6)
           DO NZB=IJK_CT(5),IJK_CT(6)
              DO NXB=IJK_CT(1),IJK_CT(2)
                 DO NYB=IJK_CT(3),IJK_CT(4)
                  READ(LU_CBD,*) TFLAG
                  IF (TFLAG.NE.0) THEN
                   NLP  = NLP + 1
                   NLP_TREE = NLP_TREE + 1
                   IF (NLP>NLPDIM) THEN
                    CALL RE_ALLOCATE_PARTICLES(1,NM,0,1000)
                    PARTICLE=>MESHES(NM)%PARTICLE
                   ENDIF
                   LP=>PARTICLE(NLP)
                   LP%TAG = PARTICLE_TAG
                   !Decide whether to show particle in smokeview
                   IF (MOD(NLP,PC%SAMPLING)==0) THEN
                      LP%SHOW=.TRUE.
                   ELSE
                      LP%SHOW=.FALSE.
                   ENDIF
                   LP%X = REAL(NXB,EB)
                   LP%Y = REAL(NYB,EB)
                   LP%Z = REAL(NZB,EB)
                   LP%CLASS = IPC
                   LP%PWT   = 1._EB  ! This is not used, but it is necessary to assign a non-zero weight factor to each particle
                   VEG_PRESENT_FLAG(NXB,NYB,NZB) = .TRUE.
                   LP%VEG_VOLFRACTION = 1._EB
                   LP%VEG_FUEL_MASS = TFLAG
                  ENDIF
               ENDDO 
             ENDDO
           ENDDO
           CLOSE(LU_CBD)
           NLP_VEG_FUEL = NLP_TREE
      ENDIF
   ENDIF IF_CUSTOM_VEGETATION 
   
!
! Build a ring volume of vegetation fuel
!
   IF_RING_VEGETATION_BUILD: IF (VEG_FUEL_GEOM(NCT) == 'RING') THEN
       RING_TREE_PRESENT = .TRUE.
       N_CFCR_TREE = TREE_CFCR_INDEX(NCT)
       N_RING_TREE = TREE_RING_INDEX(NCT)
       K_BOTTOM_RING = 0
       DZ_RING       = 0.0_EB
       OUTER_RADIUS  = 0.5_EB*CROWN_W(N_CFCR_TREE)
       RING_BOTTOM   = Z_TREE(N_CFCR_TREE) + CROWN_B_H(N_CFCR_TREE)
       RING_TOP      = Z_TREE(N_CFCR_TREE) + TREE_H(N_CFCR_TREE)
!  print*,'--------- NM = ',nm
!  print*,outer_radius
       DO II=1,IBAR-1
        IF(X(II) <= OUTER_RADIUS .AND. X(II+1) > OUTER_RADIUS) I_OUTER_RING = II
       ENDDO
!  print*,i_outer_ring,nct
!  print*,dx(i_outer_ring),ring_thickness_veg(nct)
!  DX_RING = MAX(DX(I_OUTER_RING),RING_THICKNESS_VEG(N_RING_TREE))
       DX_RING = DX(1)
       INNER_RADIUS = OUTER_RADIUS - DX_RING
       DO KK=1,KBAR-1
        IF(Z(KK) <= RING_BOTTOM .AND. Z(KK+1) > RING_BOTTOM) K_BOTTOM_RING = KK
       ENDDO
       IF (K_BOTTOM_RING > 0) DZ_RING  = MAX(DZ(K_BOTTOM_RING),RING_TOP-RING_BOTTOM)
       RING_TOP = RING_BOTTOM + DZ_RING
       NLP_TREE = 0
!
       DO NZB=1,KBAR
         IF (Z(NZB)>=RING_BOTTOM .AND. Z(NZB)<=RING_TOP) THEN
          PARTICLE_TAG = PARTICLE_TAG + NMESHES
          DO NXB = 1,IBAR
           DO NYB = 1,JBAR
            RCELL = SQRT((X(NXB)-X_TREE(N_CFCR_TREE))**2 + (Y(NYB)-Y_TREE(N_CFCR_TREE))**2)
            IF (RCELL <= OUTER_RADIUS .AND. RCELL >= INNER_RADIUS) THEN
             NLP  = NLP + 1
             NLP_TREE = NLP_TREE + 1
             IF (NLP>NLPDIM) THEN
              CALL RE_ALLOCATE_PARTICLES(1,NM,0,1000)
              PARTICLE=>MESHES(NM)%PARTICLE
             ENDIF
             LP=>PARTICLE(NLP)
             LP%VEG_VOLFRACTION = 1._EB
             LP%TAG = PARTICLE_TAG
             LP%X = REAL(NXB,EB)
             LP%Y = REAL(NYB,EB)
             LP%Z = REAL(NZB,EB)
             LP%CLASS = IPC
             LP%PWT   = 1._EB  ! This is not used, but it is necessary to assign a non-zero weight factor to each particle
             VEG_PRESENT_FLAG(NXB,NYB,NZB) = .TRUE.
            ENDIF
           ENDDO   
          ENDDO 
         ENDIF
       ENDDO
       NLP_VEG_FUEL = NLP_TREE
   ENDIF IF_RING_VEGETATION_BUILD
!
! For the current vegetation type (particle class) assign one fuel 
! element (PARTICLE) to each grid cell and initialize PARTICLE properties
! (this is precautionary needs more test to determine its necessity)
!
   REP_VEG_ELEMS: DO I=NLP-NLP_VEG_FUEL+1,NLP
    LP=>PARTICLE(I)
    LP%IGNITOR = .FALSE.
    DO NZB=0,KBAR
     DO NXB=0,IBAR
      GRID_LOOP: DO NYB=0,JBAR
      IF (.NOT. VEG_PRESENT_FLAG(NXB,NYB,NZB)) CYCLE GRID_LOOP
       IF (REAL(NXB,EB)==LP%X .AND. REAL(NYB,EB)==LP%Y .AND. REAL(NZB,EB)==LP%Z) THEN 
        IF(CELL_TAKEN_FLAG(NXB,NYB,NZB)) THEN
         LP%R = 0.0001_EB*PC%KILL_RADIUS
         CYCLE REP_VEG_ELEMS
        ENDIF
        CELL_TAKEN_FLAG(NXB,NYB,NZB) = .TRUE.
        LP%X = X(NXB) - 0.5_EB*DX(NXB)
        LP%Y = Y(NYB) - 0.5_EB*DY(NYB)
        LP%Z = Z(NZB) - 0.5_EB*DZ(NZB)
        IF (VEG_FUEL_GEOM(NCT) == 'RECTANGLE'.OR.VEG_FUEL_GEOM(NCT) == 'CUSTOM')THEN
         LP%X = X(NXB) + 0.5_EB*DX(NXB)
         LP%Y = Y(NYB) + 0.5_EB*DY(NYB)
         LP%Z = Z(NZB) + 0.5_EB*DZ(NZB)
        ENDIF
        TREE_MESH(NM) = .TRUE.
        LP%SHOW = .TRUE.
        LP%T   = 0.
        LP%U = 0.
        LP%V = 0.
        LP%W = 0.
        LP%R =  2./PC%VEG_SV !cylinder, Porterie
        LP%IOR = 0
        IF (VEG_FUEL_GEOM(NCT).NE.'CUSTOM') LP%VEG_FUEL_MASS=PC%VEG_BULK_DENSITY
        LP%VEG_MOIST_MASS = PC%VEG_MOISTURE*LP%VEG_FUEL_MASS
!       LP%VEG_CHAR_MASS  = PC%VEG_BULK_DENSITY*PC%VEG_CHAR_FRACTION
        LP%VEG_CHAR_MASS  = 0.0_EB
        LP%VEG_ASH_MASS   = 0.0_EB
        LP%VEG_PACKING_RATIO = LP%VEG_FUEL_MASS/PC%VEG_DENSITY 
        LP%VEG_SV            = PC%VEG_SV 
        LP%VEG_KAPPA = 0.25*PC%VEG_SV*LP%VEG_FUEL_MASS/PC%VEG_DENSITY
        LP%TMP = PC%VEG_INITIAL_TEMPERATURE
        LP%VEG_IGNITED = .FALSE.
        IF(IGN_ELEMS(NCT)) THEN
          IGNITOR_PRESENT = .TRUE.
          LP%TMP = TMPA
          LP%IGNITOR = .TRUE.
          N_IGN = TREE_IGN_INDEX(NCT)
          LP%VEG_IGN_TON      = TON_IGN_ELEMS(N_IGN)
          LP%VEG_IGN_TOFF     = TOFF_IGN_ELEMS(N_IGN)
          LP%VEG_IGN_TRAMPON  = T_RAMPON_IGN_ELEMS(N_IGN)
          LP%VEG_IGN_TRAMPOFF = T_RAMPOFF_IGN_ELEMS(N_IGN)
        ENDIF
        LP%VEG_EMISS = 4._EB*SIGMA*LP%VEG_KAPPA*LP%TMP**4
        LP%VEG_DIVQR = 0.0_EB
        LP%VEG_N_TREE_OUTPUT = 0
!       TREE_MESH_OUT(NM) = .FALSE.
        IF (N_TREE_OUT(NCT) /= 0) THEN
!print*,'in vege1'
         LP%VEG_N_TREE_OUTPUT = N_TREE_OUT(NCT)
         LP%IOR = 0 !airborne static PARTICLE
!        TREE_MESH_OUT(NM) = .TRUE.
        ENDIF
        CYCLE REP_VEG_ELEMS
       ENDIF
      ENDDO GRID_LOOP
     ENDDO
    ENDDO
   ENDDO REP_VEG_ELEMS
!
!print*,'in vege 2: NM,NCT,N_TREE_OUT(NCT)', NM,NCT,N_TREE_OUT(NCT)
!print*,'in vege 2: NLP,NM,TREE_MESH_OUT(NM),NCT',NLP,NM,TREE_MESH_OUT(NM),NCT
ENDDO TREE_LOOP

CALL REMOVE_PARTICLES(0._EB,NM)

!Fill veg output arrays with initial values
IF (N_TREES_OUT > 0) THEN 
  CALL POINT_TO_MESH(NM)
  TREE_OUTPUT_DATA(:,:,NM) = 0._EB
  PARTICLE_LOOP: DO I=1,NLP
   LP=>PARTICLE(I)
   N_TREE = LP%VEG_N_TREE_OUTPUT
   IF (N_TREE /= 0) THEN
     CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
     V_CELL = DX(II)*DY(JJ)*DZ(KK)
     TREE_OUTPUT_DATA(N_TREE,1,NM) = TREE_OUTPUT_DATA(N_TREE,1,NM) + LP%TMP - 273._EB !C
     TREE_OUTPUT_DATA(N_TREE,2,NM) = TREE_OUTPUT_DATA(N_TREE,2,NM) + TMPA   - 273._EB !C
     TREE_OUTPUT_DATA(N_TREE,3,NM) = TREE_OUTPUT_DATA(N_TREE,3,NM) + LP%VEG_FUEL_MASS*V_CELL !kg
     TREE_OUTPUT_DATA(N_TREE,4,NM) = TREE_OUTPUT_DATA(N_TREE,4,NM) + LP%VEG_MOIST_MASS*V_CELL !kg
     TREE_OUTPUT_DATA(N_TREE,5,NM) = TREE_OUTPUT_DATA(N_TREE,5,NM) + LP%VEG_CHAR_MASS*V_CELL !kg
     TREE_OUTPUT_DATA(N_TREE,6,NM) = TREE_OUTPUT_DATA(N_TREE,6,NM) + LP%VEG_ASH_MASS*V_CELL !kg
     TREE_OUTPUT_DATA(N_TREE,7,NM) = TREE_OUTPUT_DATA(N_TREE,7,NM) + LP%VEG_DIVQR*V_CELL*0.001_EB !kW
     TREE_OUTPUT_DATA(N_TREE,8,NM) = TREE_OUTPUT_DATA(N_TREE,8,NM) + LP%VEG_DIVQR*V_CELL*0.001_EB !kW
     TREE_OUTPUT_DATA(N_TREE,10,NM) = 0.0_EB !kg
     TREE_OUTPUT_DATA(N_TREE,11,NM) = 0.0_EB !kW
   ENDIF
  ENDDO PARTICLE_LOOP
ENDIF

!Deallocate arrays 
DEALLOCATE(VEG_PRESENT_FLAG)
DEALLOCATE(CELL_TAKEN_FLAG)
END SUBROUTINE INITIALIZE_RAISED_VEG

SUBROUTINE DEALLOCATE_VEG_ARRAYS
!Deallocate arrays used to initialize vegetation particles

IF (CONE_TREE_PRESENT .OR. FRUSTUM_TREE_PRESENT .OR. CYLINDER_TREE_PRESENT .OR. RING_TREE_PRESENT) THEN
 DEALLOCATE(TREE_CFCR_INDEX) 
 DEALLOCATE(X_TREE)
 DEALLOCATE(Y_TREE)
 DEALLOCATE(Z_TREE)
 DEALLOCATE(TREE_H)
 DEALLOCATE(CROWN_B_H)
 DEALLOCATE(CROWN_W)
ENDIF
IF (FRUSTUM_TREE_PRESENT) THEN
 DEALLOCATE(TREE_FRUSTUM_INDEX)
 DEALLOCATE(CROWN_W_TOP)
 DEALLOCATE(CROWN_W_BOTTOM)
ENDIF

IF (RECTANGLE_TREE_PRESENT) THEN
 DEALLOCATE(TREE_RECT_INDEX)
 DEALLOCATE(XS_RECT_VEG)
 DEALLOCATE(XF_RECT_VEG)
 DEALLOCATE(YS_RECT_VEG)
 DEALLOCATE(YF_RECT_VEG)
 DEALLOCATE(ZS_RECT_VEG)
 DEALLOCATE(ZF_RECT_VEG)
ENDIF

IF (CUSTOM_TREE_PRESENT) THEN
 DEALLOCATE(TREE_CUSTOM_INDEX)
ENDIF

IF (IGNITOR_PRESENT) THEN
 DEALLOCATE(IGN_ELEMS)
 DEALLOCATE(TREE_IGN_INDEX)
 DEALLOCATE(TON_IGN_ELEMS)
 DEALLOCATE(TOFF_IGN_ELEMS)
 DEALLOCATE(T_RAMPOFF_IGN_ELEMS)
 DEALLOCATE(T_RAMPON_IGN_ELEMS)
ENDIF
END SUBROUTINE DEALLOCATE_VEG_ARRAYS



SUBROUTINE RAISED_VEG_MASS_ENERGY_TRANSFER(T,NM)
    
! Mass and energy transfer between gas and raised vegetation fuel elements 

USE PHYSICAL_FUNCTIONS, ONLY : GET_MASS_FRACTION,GET_SPECIFIC_HEAT
USE MATH_FUNCTIONS, ONLY : AFILL2
USE TRAN, ONLY: GET_IJK
!arrays for debugging
REAL(EB), POINTER, DIMENSION(:,:,:) :: HOLD1,HOLD2,HOLD3,HOLD4
REAL(EB), POINTER, DIMENSION(:,:,:) :: UU,VV,WW !,RHOP

REAL(EB) :: RE_D,RCP_GAS,CP_GAS
REAL(EB) :: RDT,T,V_CELL,V_VEG
REAL(EB) :: CP_ASH,CP_H2O,CP_CHAR,H_VAP_H2O,TMP_H2O_BOIL
REAL(EB) :: K_AIR,MASS_GAS,MU_AIR,RHO_GAS,RRHO_GAS_NEW,TMP_GAS,UBAR,VBAR,WBAR,UREL,VREL,WREL
REAL(EB) :: CHAR_FCTR,CHAR_FCTR2,CP_VEG,DTMP_VEG,MPV_MOIST,MPV_MOIST_MIN,DMPV_VEG,MPV_VEG,MPV_VEG_MIN, &
            SV_VEG,TMP_VEG,TMP_VEG_NEW
REAL(EB) :: TMP_IGNITOR
REAL(EB) :: MPV_ADDED,MPV_MOIST_LOSS,MPV_VOLIT,MPV_MOIST_LOSS_MAX,MPV_VOLIT_MAX
REAL(EB) :: QCON_VEG,QNET_VEG,QRAD_VEG,QREL,TMP_GMV,Q_FOR_DRYING,Q_VOLIT,Q_FOR_VOLIT, &
            Q_UPTO_VOLIT
REAL(EB) :: H_SENS_VEG_VOLIT,Q_ENTHALPY,Q_VEG_MOIST,Q_VEG_VOLIT,Q_VEG_CHAR
REAL(EB) :: MW_AVERAGE,MW_VEG_MOIST_TERM,MW_VEG_VOLIT_TERM
REAL(EB) :: XI,YJ,ZK
REAL(EB) :: A_H2O_VEG,E_H2O_VEG,A_PYR_VEG,E_PYR_VEG,H_PYR_VEG,R_H_PYR_VEG
REAL(EB) :: A_CHAR_VEG,E_CHAR_VEG,BETA_CHAR_VEG,NU_CHAR_VEG,NU_ASH_VEG,NU_O2_CHAR_VEG, &
            MPV_ASH,MPV_ASH_MAX,MPV_CHAR,MPV_CHAR_LOSS,MPV_CHAR_MIN,MPV_CHAR_CO2,Y_O2, &
            H_CHAR_VEG ,ORIG_PACKING_RATIO,CP_VEG_FUEL_AND_CHAR_MASS,CP_MASS_VEG_SOLID,     &
            TMP_CHAR_MAX
REAL(EB) :: ZZ_GET(0:N_TRACKED_SPECIES)
INTEGER :: I,II,JJ,KK,IIX,JJY,KKZ,IPC,N_TREE,I_FUEL
INTEGER, INTENT(IN) :: NM
LOGICAL :: VEG_DEGRADATION_LINEAR,VEG_DEGRADATION_ARRHENIUS
INTEGER :: IDT,NDT_CYCLES
REAL(EB) :: FCTR_DT_CYCLES,FCTR_RDT_CYCLES,Q_VEG_CHAR_TOTAL,MPV_CHAR_CO2_TOTAL,MPV_CHAR_LOSS_TOTAL, &
            MPV_MOIST_LOSS_TOTAL,MPV_VOLIT_TOTAL
REAL(EB) :: VEG_CRITICAL_MASSFLUX,VEG_CRITICAL_MASSSOURCE

!place holder
REAL(EB) :: RCP_TEMPORARY

!Debug
REAL(EB)TOTAL_BULKDENS_MOIST,TOTAL_BULKDENS_DRY_FUEL,TOTAL_MASS_DRY_FUEL,TOTAL_MASS_MOIST


!IF (.NOT. TREE) RETURN !Exit if no raised veg anywhere
IF (.NOT. TREE_MESH(NM)) RETURN !Exit if raised veg is not present in mesh
CALL POINT_TO_MESH(NM)

!IF (PREDICTOR) THEN
    UU => U
    VV => V
    WW => W
!   RHOP => RHO
!ELSE
!   UU => US
!   VV => VS
!   WW => WS
!   RHOP => RHOS
!ENDIF

! Initializations

RDT    = 1._EB/DT
!RCP_TEMPORARY = 1._EB/CP_GAMMA
RCP_TEMPORARY = 1._EB/1010._EB

!Critical mass flux (kg/(s m^2)
VEG_CRITICAL_MASSFLUX = 0.0025_EB !kg/s/m^2 for qradinc=50 kW/m^2, M=4% measured by McAllister Fire Safety J., 61:200-206 2013
!VEG_CRITICAL_MASSFLUX = 0.0035_EB !kg/s/m^2 largest measured by McAllister Fire Safety J., 61:200-206 2013
!VEG_CRITICAL_MASSFLUX = 999999._EB !kg/s/m^2 for testing

!Constants for Arrhenius pyrolyis and Arrhenius char oxidation models
!are from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005)
CP_H2O       = 4190._EB !J/kg/K specific heat of water
TMP_H2O_BOIL = 373.15_EB
!H_VAP_H2O    = 2259._EB*1000._EB !J/kg/K heat of vaporization of water
TMP_CHAR_MAX = 1300._EB !K

!Kinetic constants used by multiple investigators from Porterie or Morvan papers
!A_H2O_VEG      = 600000._EB !1/s sqrt(K)
!E_H2O_VEG      = 5800._EB !K
!A_PYR_VEG      = 36300._EB !1/s
!E_PYR_VEG      = 7250._EB !K
!E_CHAR_VEG     = 9000._EB !K
!BETA_CHAR_VEG  = 0.2_EB
!!NU_CHAR_VEG    = 0.3_EB
!!NU_ASH_VEG     = 0.1_EB
!NU_O2_CHAR_VEG = 1.65_EB

CP_ASH         = 800._EB !J/kg/K

!Kinetic constants used by Morvan and Porterie mostly obtained from Grishin
!H_PYR_VEG      = 418000._EB !J/kg 
!A_CHAR_VEG     = 430._EB !m/s 
!H_CHAR_VEG     = -12.0E+6_EB ! J/kg

!Kinetic constants used by Yolanda and Paul
!H_PYR_VEG      = 418000._EB !J/kg 
!A_CHAR_VEG     = 215._EB !m/s Yolanda, adjusted from Morvan, Porterie values based on HRR exp
!H_CHAR_VEG     = -32.74E+6_EB !J/kg via Susott

!Kinetic constants used by Shankar
!H_PYR_VEG      = 418._EB !J/kg Shankar
!A_CHAR_VEG     = 430._EB !m/s Porterie, Morvan
!H_CHAR_VEG     = -32.74E+6_EB !J/kg Shankar via Susott

!Kinetic constants used by me for ROS vs Slope excelsior experiments
!H_PYR_VEG      = 711000._EB !J/kg excelsior Catchpole et al. (via Susott)
!A_CHAR_VEG     = 430._EB !m/s Porterie, Morvan
!H_CHAR_VEG     = -32.74E+6_EB !J/kg via Susott

!R_H_PYR_VEG    = 1._EB/H_PYR_VEG

!D_AIR  = 2.6E-5_EB  ! Water Vapor - Air binary diffusion (m2/s at 25 C, Incropera & DeWitt, Table A.8) 
!SC_AIR = 0.6_EB     ! NU_AIR/D_AIR (Incropera & DeWitt, Chap 7, External Flow)
!PR_AIR = 0.7_EB     

! Working arrays
IF(N_TREES_OUT > 0) TREE_OUTPUT_DATA(:,:,NM) = 0._EB !for output of veg data
!DMPVDT_FM_VEG  = 0.0_EB

!Clear arrays and scalars
HOLD1 => WORK4 ; WORK4 = 0._EB
HOLD2 => WORK5 ; WORK5 = 0._EB
HOLD3 => WORK6 ; WORK6 = 0._EB
HOLD4 => WORK7 ; WORK7 = 0._EB
TOTAL_BULKDENS_MOIST    = 0.0_EB
TOTAL_BULKDENS_DRY_FUEL = 0.0_EB
TOTAL_MASS_MOIST    = 0.0_EB
TOTAL_MASS_DRY_FUEL = 0.0_EB
V_VEG               = 0.0_EB

!print*,'vege h-m transfer: NM, NLP',nm,nlp

PARTICLE_LOOP: DO I=1,NLP

 LP => PARTICLE(I)
 IPC = LP%CLASS
 PC=>PARTICLE_CLASS(IPC)
 IF (.NOT. PC%TREE) CYCLE PARTICLE_LOOP !Ensure grid cell has vegetation
 IF (PC%MASSLESS) CYCLE PARTICLE_LOOP   !Skip PARTICLE if massless

 THERMAL_CALC: IF (.NOT. PC%VEG_STEM) THEN   !compute heat transfer, etc if thermally thin

!Quantities for sub-cycling the thermal degradation time stepping
 NDT_CYCLES  = PC%VEG_NDT_SUBCYCLES !number of thermal degradation time stepping loops within one gas phase DT
 FCTR_DT_CYCLES   = 1._EB/REAL(NDT_CYCLES,EB)
 FCTR_RDT_CYCLES  = REAL(NDT_CYCLES,EB)

! Intialize quantities
 LP%VEG_MLR     = 0.0_EB
 Q_VEG_CHAR     = 0.0_EB
 Q_VEG_MOIST    = 0.0_EB
 Q_VEG_VOLIT    = 0.0_EB
 Q_UPTO_VOLIT   = 0.0_EB
 Q_VOLIT        = 0.0_EB
 MPV_MOIST_LOSS = 0.0_EB
 MPV_CHAR_LOSS  = 0.0_EB
 MPV_CHAR_CO2   = 0.0_EB
 MPV_VOLIT      = 0.0_EB
 MPV_ADDED      = 0.0_EB
 MW_VEG_MOIST_TERM = 0.0_EB
 MW_VEG_VOLIT_TERM = 0.0_EB
 CP_VEG_FUEL_AND_CHAR_MASS = 0.0_EB
 CP_MASS_VEG_SOLID         = 0.0_EB
 VEG_DEGRADATION_LINEAR    = .FALSE.
 VEG_DEGRADATION_ARRHENIUS = .FALSE.
 MPV_CHAR_CO2_TOTAL   = 0.0_EB !needed for time subcycling
 MPV_CHAR_LOSS_TOTAL  = 0.0_EB !needed for time subcycling
 MPV_MOIST_LOSS_TOTAL = 0.0_EB !needed for time subcycling
 MPV_VOLIT_TOTAL  = 0.0_EB !needed for time subcycling
 Q_VEG_CHAR_TOTAL = 0.0_EB !needed for time subcyling

! Vegetation variables
 NU_CHAR_VEG        = PC%VEG_CHAR_FRACTION
 NU_ASH_VEG         = PC%VEG_ASH_FRACTION
 CHAR_FCTR          = 1._EB - PC%VEG_CHAR_FRACTION !factor used to determine volatile mass
 CHAR_FCTR2         = 1._EB/CHAR_FCTR !factor used to determine char mass
 SV_VEG             = LP%VEG_SV !surface-to-volume ration 1/m
 TMP_VEG            = LP%TMP
 MPV_VEG            = LP%VEG_FUEL_MASS !bulk density of dry veg kg/m^3
 MPV_CHAR           = LP%VEG_CHAR_MASS !bulk density of char
 MPV_ASH            = LP%VEG_ASH_MASS  !bulk density of ash 
 MPV_MOIST          = LP%VEG_MOIST_MASS !bulk density of moisture in veg
 MPV_VEG_MIN        = PC%VEG_FUEL_MPV_MIN
 MPV_CHAR_MIN       = PC%VEG_CHAR_MPV_MIN
 MPV_MOIST_MIN      = PC%VEG_MOIST_MPV_MIN
 MPV_ASH_MAX        = PC%VEG_ASH_MPV_MAX   !maxium ash bulk density
 MPV_MOIST_LOSS_MAX = PC%VEG_DEHYDRATION_RATE_MAX*DT*FCTR_DT_CYCLES
 MPV_VOLIT_MAX      = PC%VEG_BURNING_RATE_MAX*DT*FCTR_DT_CYCLES
 ORIG_PACKING_RATIO = PC%VEG_BULK_DENSITY/PC%VEG_DENSITY 
 H_VAP_H2O          = PC%VEG_H_H2O !J/kg/K heat of vaporization of water
 A_H2O_VEG          = PC%VEG_A_H2O !1/s sqrt(K)
 E_H2O_VEG          = PC%VEG_E_H2O !K
 H_PYR_VEG          = PC%VEG_H_PYR !J/kg 
 A_PYR_VEG          = PC%VEG_A_PYR !1/s
 E_PYR_VEG          = PC%VEG_E_PYR !K
 H_CHAR_VEG         = PC%VEG_H_CHAR ! J/kg
 A_CHAR_VEG         = PC%VEG_A_CHAR !m/s 
 E_CHAR_VEG         = PC%VEG_E_CHAR !K
 BETA_CHAR_VEG      = PC%VEG_BETA_CHAR
 NU_O2_CHAR_VEG     = PC%VEG_NU_O2_CHAR

! Thermal degradation approach parameters
 IF(PC%VEG_DEGRADATION == 'LINEAR') VEG_DEGRADATION_LINEAR = .TRUE.
 IF(PC%VEG_DEGRADATION == 'ARRHENIUS') VEG_DEGRADATION_ARRHENIUS = .TRUE.

 R_H_PYR_VEG    = 1._EB/H_PYR_VEG

!Bound on volumetric mass flux
 VEG_CRITICAL_MASSSOURCE = VEG_CRITICAL_MASSFLUX*SV_VEG*LP%VEG_PACKING_RATIO

! Determine grid cell quantities of the vegetation fuel element
 CALL GET_IJK(LP%X,LP%Y,LP%Z,NM,XI,YJ,ZK,II,JJ,KK)
 IIX = FLOOR(XI+0.5_EB)
 JJY = FLOOR(YJ+0.5_EB)
 KKZ = FLOOR(ZK+0.5_EB)
 V_CELL = DX(II)*DY(JJ)*DZ(KK)

! Gas velocities in vegetation grid cell
 UBAR = AFILL2(UU,II-1,JJY,KKZ,XI-II+1,YJ-JJY+.5_EB,ZK-KKZ+.5_EB)
 VBAR = AFILL2(VV,IIX,JJ-1,KKZ,XI-IIX+.5_EB,YJ-JJ+1,ZK-KKZ+.5_EB)
 WBAR = AFILL2(WW,IIX,JJY,KK-1,XI-IIX+.5_EB,YJ-JJY+.5_EB,ZK-KK+1)
 UREL = LP%U - UBAR
 VREL = LP%V - VBAR
 WREL = LP%W - WBAR
 QREL = MAX(1.E-6_EB,SQRT(UREL*UREL + VREL*VREL + WREL*WREL))

! Gas thermophysical quantities
 TMP_GAS  = TMP(II,JJ,KK)
 RHO_GAS  = RHO(II,JJ,KK)
 MASS_GAS = RHO_GAS*V_CELL
 MU_AIR   = MU_Z(MIN(5000,NINT(TMP_GAS)),0)*SPECIES_MIXTURE(0)%MW
 K_AIR    = CPOPR*MU_AIR !W/m.K

TIME_SUBCYCLING_LOOP: DO IDT=1,NDT_CYCLES

! Veg thermophysical properties
 TMP_GMV  = TMP_GAS - TMP_VEG
 CP_VEG   = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !J/kg/K Ritchie IAFSS 1997:177-188
 CP_CHAR  = 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494

! Divergence of convective and radiative heat fluxes
!print*,'---- NM=',NM
!print*,rho_gas,qrel,sv_veg,mu_air
 RE_D     = RHO_GAS*QREL*2._EB/(SV_VEG*MU_AIR)
 !IF (TMP_VEG >= TMP_GAS )QCON_VEG = SV_VEG*(0.5_EB*K_AIR*0.683_EB*RE_D**0.466_EB)*0.5_EB*TMP_GMV !W/m^2 from Porterie
 !IF (TMP_VEG <  TMP_GAS ) QCON_VEG = TMP_GMV*1.42_EB*(ABS(TMP_GMV)/DZ(KK))**0.25_EB !Holman
!QCON_VEG = SV_VEG*(0.5_EB*K_AIR*0.683_EB*RE_D**0.466_EB)*0.5_EB*TMP_GMV !W/m^2 from Porterie
 QCON_VEG = TMP_GMV*1.42_EB*(ABS(TMP_GMV)/DZ(KK))**0.25_EB !Holman
 QCON_VEG = SV_VEG*LP%VEG_PACKING_RATIO*QCON_VEG !W/m^3
 LP%VEG_DIVQC = QCON_VEG
 QRAD_VEG = LP%VEG_DIVQR

! Divergence of net heat flux
 QNET_VEG = QCON_VEG + QRAD_VEG !W/m^3

! Update temperature of vegetation
!CP_VEG_FUEL_AND_CHAR_MASS = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR
!DTMP_VEG    = DT*QNET_VEG/(CP_VEG_FUEL_AND_CHAR_MASS + CP_H2O*MPV_MOIST)
 CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
 DTMP_VEG    = FCTR_DT_CYCLES*DT*QNET_VEG/(CP_MASS_VEG_SOLID + CP_H2O*MPV_MOIST)
!print*,'vege:tmpveg,qnet_veg,cp_mass_veg_solid',tmp_veg,qnet_veg,cp_mass_veg_solid
 TMP_VEG_NEW = TMP_VEG + DTMP_VEG

! Set temperature of inert ignitor elements
 IF(LP%IGNITOR) THEN
  TMP_IGNITOR = PC%VEG_INITIAL_TEMPERATURE
  TMP_VEG_NEW = TMP_GAS
  IF(T>=LP%VEG_IGN_TON .AND. T<=LP%VEG_IGN_TON+LP%VEG_IGN_TRAMPON) THEN
    TMP_VEG_NEW = &
      TMPA + (TMP_IGNITOR-TMPA)*(T-LP%VEG_IGN_TON)/LP%VEG_IGN_TRAMPON
  ENDIF  
  IF(T>LP%VEG_IGN_TON+LP%VEG_IGN_TRAMPON) TMP_VEG_NEW = TMP_IGNITOR
  IF(T>=LP%VEG_IGN_TOFF .AND. T<=LP%VEG_IGN_TOFF+LP%VEG_IGN_TRAMPOFF)THEN 
    TMP_VEG_NEW = &
      TMP_IGNITOR - (TMP_IGNITOR-TMP_GAS)*(T-LP%VEG_IGN_TOFF)/LP%VEG_IGN_TRAMPOFF
  ENDIF
  IF(T > LP%VEG_IGN_TOFF+LP%VEG_IGN_TRAMPOFF) THEN
   LP%R = 0.0001_EB*PC%KILL_RADIUS !remove ignitor element
   TMP_VEG_NEW = TMP_GAS
  ENDIF
 ENDIF

!      ************** Fuel Element Linear Pyrolysis Degradation model *************************
! Drying occurs if qnet > 0 with Tveg held at 100 c
! Pyrolysis occurs if qnet > 0 according to Morvan & Dupuy empirical formula. Linear
! temperature dependence with qnet factor. 
! Char oxidation occurs if qnet > 0 (user must request char ox) after pyrolysis is completed.
!
 IF_VEG_DEGRADATION_LINEAR: IF(VEG_DEGRADATION_LINEAR) THEN
 IF_NET_HEAT_INFLUX: IF (QNET_VEG > 0.0_EB .AND. .NOT. LP%IGNITOR) THEN !dehydrate or pyrolyze 

! Drying of fuel element vegetation 
   IF_DEHYDRATION: IF (MPV_MOIST > MPV_MOIST_MIN .AND. TMP_VEG_NEW > TMP_H2O_BOIL) THEN
     Q_FOR_DRYING   = (TMP_VEG_NEW - TMP_H2O_BOIL)/DTMP_VEG * QNET_VEG
     MPV_MOIST_LOSS = MIN(DT*Q_FOR_DRYING/H_VAP_H2O,MPV_MOIST-MPV_MOIST_MIN)
     MPV_MOIST_LOSS = LP%VEG_VOLFRACTION*MPV_MOIST_LOSS !accounts for veg not filling grid cell in z
     MPV_MOIST_LOSS = MIN(MPV_MOIST_LOSS,MPV_MOIST_LOSS_MAX) !use specified max
     TMP_VEG_NEW       = TMP_H2O_BOIL
     LP%VEG_MOIST_MASS = MPV_MOIST - MPV_MOIST_LOSS !kg/m^3
     IF (LP%VEG_MOIST_MASS <= MPV_MOIST_MIN) LP%VEG_MOIST_MASS = 0.0_EB
     Q_VEG_MOIST       = MPV_MOIST_LOSS*CP_H2O*(TMP_VEG_NEW - TMPA)
     MW_VEG_MOIST_TERM = MPV_MOIST_LOSS/MW_H2O
!    IF (I == 1) print*,MPV_MOIST,MPV_MOIST_LOSS
   ENDIF IF_DEHYDRATION

! Volitalization of fuel element vegetation
  IF_VOLITALIZATION: IF(MPV_MOIST <= MPV_MOIST_MIN) THEN

    IF_MD_VOLIT: IF(MPV_VEG > MPV_VEG_MIN .AND. TMP_VEG_NEW >= 400._EB) THEN !Morvan & Dupuy volitalization
     Q_UPTO_VOLIT = CP_MASS_VEG_SOLID*MAX((400._EB-TMP_VEG),0._EB)
     Q_FOR_VOLIT  = DT*QNET_VEG - Q_UPTO_VOLIT
     Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(MIN(500._EB,TMP_VEG)-400._EB)

!    MPV_VOLIT    = Q_VOLIT*R_H_PYR_VEG
     MPV_VOLIT    = CHAR_FCTR*Q_VOLIT*R_H_PYR_VEG
     MPV_VOLIT    = MAX(MPV_VOLIT,0._EB)
     MPV_VOLIT    = LP%VEG_VOLFRACTION*MPV_VOLIT !accounts for veg not filling grid cell in z
     MPV_VOLIT    = MIN(MPV_VOLIT,MPV_VOLIT_MAX) !user specified max

     DMPV_VEG     = CHAR_FCTR2*MPV_VOLIT
     DMPV_VEG     = MIN(DMPV_VEG,(MPV_VEG - MPV_VEG_MIN))
     MPV_VEG      = MPV_VEG - DMPV_VEG

     MPV_VOLIT    = CHAR_FCTR*DMPV_VEG
!    MPV_CHAR     = MPV_CHAR + NU_CHAR_VEG*MPV_VOLIT !kg/m^3
!    MPV_CHAR     = MPV_CHAR + PC%VEG_CHAR_FRACTION*MPV_VOLIT !kg/m^3
     MPV_CHAR     = MPV_CHAR + PC%VEG_CHAR_FRACTION*DMPV_VEG !kg/m^3
     Q_VOLIT      = MPV_VOLIT*H_PYR_VEG
     CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR 
     TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/CP_MASS_VEG_SOLID
     TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB) !set to high pyrol temp if too hot

!Handle veg. fuel elements if element mass <= prescribed minimum
     IF (MPV_VEG <= MPV_VEG_MIN) THEN
       MPV_VEG = MPV_VEG_MIN
       IF(PC%VEG_REMOVE_CHARRED .AND. .NOT. PC%VEG_CHAR_OXIDATION) & 
            LP%R = 0.0001_EB*PC%KILL_RADIUS !fuel element will be removed
     ENDIF
!Enthalpy of fuel element volatiles using Cp,volatiles(T) from Ritchie
     H_SENS_VEG_VOLIT = 0.0445_EB*(TMP_VEG**1.5_EB - TMP_GAS**1.5_EB) - 0.136_EB*(TMP_VEG - TMP_GAS)
     H_SENS_VEG_VOLIT = H_SENS_VEG_VOLIT*1000._EB !J/kg
     Q_VEG_VOLIT      = CHAR_FCTR*MPV_VOLIT*H_SENS_VEG_VOLIT !J/m^3
     MW_VEG_VOLIT_TERM= MPV_VOLIT/SPECIES(FUEL_INDEX)%MW
    ENDIF IF_MD_VOLIT
    LP%VEG_FUEL_MASS = MPV_VEG
    LP%VEG_CHAR_MASS = MPV_CHAR !kg/m^3
  ENDIF IF_VOLITALIZATION

 ENDIF IF_NET_HEAT_INFLUX

!Char oxidation of fuel element with the Linear pyrolysis model from Morvan and Dupuy, Comb.
!Flame, 138:199-210 (2004)
!(note that this can be handled only approximately with the conserved
!scalar based gas-phase combustion model - the oxygen is consumed by
!the char oxidation reaction is not accounted for since it would be inconsistent with the state
!relation for oxygen that is based on the conserved scalar approach used for gas phase
!combustion)
   IF_CHAR_OXIDATION_LIN: IF (PC%VEG_CHAR_OXIDATION .AND. MPV_MOIST <= MPV_MOIST_MIN) THEN

      ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
      CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
      MPV_CHAR_LOSS = DT*RHO_GAS*Y_O2*A_CHAR_VEG/NU_O2_CHAR_VEG*SV_VEG*LP%VEG_PACKING_RATIO*  &
                       EXP(-E_CHAR_VEG/TMP_VEG)*(1+BETA_CHAR_VEG*SQRT(2._EB*RE_D))
      MPV_CHAR      = MAX(MPV_CHAR - MPV_CHAR_LOSS,0.0_EB)
      MPV_CHAR_LOSS = LP%VEG_CHAR_MASS - MPV_CHAR
      MPV_CHAR_CO2  = (1._EB + NU_O2_CHAR_VEG - NU_ASH_VEG)*MPV_CHAR_LOSS
      LP%VEG_CHAR_MASS  = MPV_CHAR !kg/m^3
      CP_MASS_VEG_SOLID = MPV_VEG*CP_VEG + MPV_CHAR*CP_CHAR + MPV_ASH*CP_ASH

! Reduce fuel element size based on char consumption
!     IF (MPV_VEG <= MPV_VEG_MIN) THEN !charring reduce veg elem size
       LP%VEG_PACKING_RATIO = LP%VEG_PACKING_RATIO - MPV_CHAR_LOSS/(PC%VEG_DENSITY*PC%VEG_CHAR_FRACTION)
       LP%VEG_SV     = PC%VEG_SV*(ORIG_PACKING_RATIO/LP%VEG_PACKING_RATIO)**0.333_EB 
       LP%VEG_KAPPA  = 0.25_EB*LP%VEG_SV*LP%VEG_PACKING_RATIO
!     ENDIF

!remove fuel element if char ox is complete
      IF (MPV_CHAR <= MPV_CHAR_MIN .AND. MPV_VEG <= MPV_VEG_MIN) THEN 
        CP_MASS_VEG_SOLID = MPV_CHAR_MIN*CP_CHAR + MPV_ASH*CP_ASH
        LP%VEG_CHAR_MASS = 0.0_EB
        IF(PC%VEG_REMOVE_CHARRED) LP%R = 0.0001_EB*PC%KILL_RADIUS
      ENDIF

      Q_VEG_CHAR   = MPV_CHAR_LOSS*H_CHAR_VEG 
      TMP_VEG_NEW  = TMP_VEG_NEW - PC%VEG_CHAR_ENTHALPY_FRACTION*Q_VEG_CHAR/CP_MASS_VEG_SOLID
      TMP_VEG_NEW  = MIN(TMP_CHAR_MAX,TMP_VEG_NEW)
!     print*,'vege: q_veg_char,temp_veg_new,',q_veg_char,tmp_veg_new
!          print*,'------------------'
!    ENDIF IF_CHAR_OXIDATION_LIN_2

   ENDIF IF_CHAR_OXIDATION_LIN
  
! ENDIF IF_VOLITALIZATION

!ENDIF IF_NET_HEAT_INFLUX
 ENDIF IF_VEG_DEGRADATION_LINEAR

!      ************** Fuel Element Arrehnius Degradation model *************************
! Drying and pyrolysis of fuel element occur according to Arrehnius expressions obtained 
! from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005
! Predicting wildland fire behavior and emissions using a fine-scale physical
! model
!
 IF_VEG_DEGRADATION_ARRHENIUS: IF(VEG_DEGRADATION_ARRHENIUS) THEN

  IF_NOT_IGNITOR1: IF (.NOT. LP%IGNITOR) THEN !dehydrate or pyrolyze 

! Drying of fuel element vegetation 
   IF_DEHYDRATION_2: IF (MPV_MOIST > MPV_MOIST_MIN) THEN
     MPV_MOIST_LOSS = MIN(FCTR_DT_CYCLES*DT*MPV_MOIST*A_H2O_VEG*EXP(-E_H2O_VEG/TMP_VEG)/SQRT(TMP_VEG), &
                          MPV_MOIST-MPV_MOIST_MIN)
     MPV_MOIST_LOSS = MIN(MPV_MOIST_LOSS,MPV_MOIST_LOSS_MAX) !use specified max
     MPV_MOIST      = MPV_MOIST - MPV_MOIST_LOSS
     LP%VEG_MOIST_MASS = MPV_MOIST !kg/m^3
     IF (MPV_MOIST <= MPV_MOIST_MIN) LP%VEG_MOIST_MASS = 0.0_EB
     MW_VEG_MOIST_TERM = MPV_MOIST_LOSS/MW_H2O
     Q_VEG_MOIST  = MPV_MOIST_LOSS*CP_H2O*(TMP_VEG - TMPA)
!    IF (I == 1) print*,MPV_MOIST,MPV_MOIST_LOSS
   ENDIF IF_DEHYDRATION_2

! Volitalization of fuel element vegetation
  IF_VOLITALIZATION_2: IF(MPV_VEG > MPV_VEG_MIN) THEN
     MPV_VOLIT    = FCTR_DT_CYCLES*DT*CHAR_FCTR*MPV_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG)
     MPV_VOLIT    = MIN(MPV_VOLIT,MPV_VOLIT_MAX) !user specified max

     DMPV_VEG     = CHAR_FCTR2*MPV_VOLIT
     DMPV_VEG     = MIN(DMPV_VEG,(MPV_VEG - MPV_VEG_MIN))
     MPV_VEG      = MPV_VEG - DMPV_VEG

     MPV_VOLIT    = CHAR_FCTR*DMPV_VEG 
     MPV_CHAR     = MPV_CHAR + PC%VEG_CHAR_FRACTION*DMPV_VEG !kg/m^3
!    MPV_CHAR     = MPV_CHAR + PC%VEG_CHAR_FRACTION*MPV_VOLIT !kg/m^3
     CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH

!Yolanda's
!    MPV_VOLIT    = FCTR_DT_CYCLES*DT*MPV_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG)
!!   MPV_VOLIT    = MAX(MPV_VOLIT,0._EB)
!    MPV_VOLIT    = MIN(MPV_VOLIT,CHAR_FCTR*MPV_VOLIT_MAX) !user specified max
!    MPV_VOLIT    = MIN(MPV_VOLIT,(MPV_VEG-MPV_VEG_MIN))
!    MPV_VEG      = MPV_VEG - MPV_VOLIT
!    MPV_CHAR     = MPV_CHAR + PC%VEG_CHAR_FRACTION*MPV_VOLIT !kg/m^3
!    CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
!    MPV_VOLIT    = CHAR_FCTR*MPV_VOLIT

!Handle veg. fuel elements if original element mass <= prescribed minimum
     IF (MPV_VEG <= MPV_VEG_MIN) THEN
!      MPV_VEG = MPV_VEG_MIN
       MPV_VEG = 0.0_EB
       CP_MASS_VEG_SOLID = CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
       IF(PC%VEG_REMOVE_CHARRED .AND. .NOT. PC%VEG_CHAR_OXIDATION) LP%R = 0.0001_EB*PC%KILL_RADIUS !remove part
     ENDIF
!Enthalpy of fuel element volatiles using Cp,volatiles(T) from Ritchie
     H_SENS_VEG_VOLIT = 0.0445_EB*(TMP_VEG**1.5_EB - TMP_GAS**1.5_EB) - 0.136_EB*(TMP_VEG - TMP_GAS)
     H_SENS_VEG_VOLIT = H_SENS_VEG_VOLIT*1000._EB !J/kg
     Q_VEG_VOLIT      = MPV_VOLIT*H_SENS_VEG_VOLIT !J
     MW_VEG_VOLIT_TERM= MPV_VOLIT/SPECIES(FUEL_INDEX)%MW
  ENDIF IF_VOLITALIZATION_2
  LP%VEG_FUEL_MASS = MPV_VEG
  LP%VEG_CHAR_MASS = MPV_CHAR


!Char oxidation or fuel element within the Arrhenius pyrolysis model
!(note that this can be handled only approximately with the conserved
!scalar based gas-phase combustion model - no gas phase oxygen is consumed by
!the char oxidation reaction since it would be inconsistent with the state
!relation for oxygen based on the conserved scalar approach for gas phase
!combustion)
! MPV_CHAR_LOSS = 0.0_EB
  IF_CHAR_OXIDATION: IF (PC%VEG_CHAR_OXIDATION) THEN
     ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
     CALL GET_MASS_FRACTION(ZZ_GET,O2_INDEX,Y_O2)
     MPV_CHAR_LOSS = FCTR_DT_CYCLES*DT*RHO_GAS*Y_O2*A_CHAR_VEG/NU_O2_CHAR_VEG*SV_VEG*LP%VEG_PACKING_RATIO*  &
                      EXP(-E_CHAR_VEG/TMP_VEG)*(1+BETA_CHAR_VEG*SQRT(2._EB*RE_D))
     MPV_CHAR_LOSS = MIN(MPV_CHAR,MPV_CHAR_LOSS)
     MPV_CHAR      = MPV_CHAR - MPV_CHAR_LOSS
     MPV_ASH       = MPV_ASH + NU_ASH_VEG*MPV_CHAR_LOSS
     MPV_CHAR_CO2  = (1._EB + NU_O2_CHAR_VEG - NU_ASH_VEG)*MPV_CHAR_LOSS
     CP_MASS_VEG_SOLID = CP_VEG*MPV_VEG + CP_CHAR*MPV_CHAR + CP_ASH*MPV_ASH
     LP%VEG_CHAR_MASS = MPV_CHAR !kg/m^3
     LP%VEG_ASH_MASS  = MPV_ASH

! Reduce veg element size based on char consumption
!    LP%VEG_PACKING_RATIO = LP%VEG_PACKING_RATIO - MPV_CHAR_LOSS/(PC%VEG_DENSITY*PC%VEG_CHAR_FRACTION)
!    LP%VEG_SV     = PC%VEG_SV*(ORIG_PACKING_RATIO/LP%VEG_PACKING_RATIO)**0.333_EB 
!    LP%VEG_KAPPA  = 0.25_EB*LP%VEG_SV*LP%VEG_PACKING_RATIO

     IF (MPV_CHAR <= MPV_CHAR_MIN .AND. MPV_VEG <= MPV_VEG_MIN) THEN 
!    IF (MPV_ASH >= MPV_ASH_MAX .AND. MPV_VEG <= MPV_VEG_MIN) THEN 
!      CP_MASS_VEG_SOLID = CP_CHAR*MPV_CHAR_MIN
       CP_MASS_VEG_SOLID = CP_ASH*MPV_ASH
       LP%VEG_CHAR_MASS = 0.0_EB
       IF(PC%VEG_REMOVE_CHARRED) LP%R = 0.0001_EB*PC%KILL_RADIUS !fuel element will be removed
     ENDIF
!  ENDIF IF_CHAR_OXIDATION_2
  ENDIF IF_CHAR_OXIDATION

  Q_VEG_CHAR        = MPV_CHAR_LOSS*H_CHAR_VEG 
  Q_VEG_CHAR_TOTAL  = Q_VEG_CHAR_TOTAL + Q_VEG_CHAR
  TMP_VEG_NEW  = TMP_VEG_NEW - (MPV_MOIST_LOSS*H_VAP_H2O + MPV_VOLIT*H_PYR_VEG + & 
                                PC%VEG_CHAR_ENTHALPY_FRACTION*Q_VEG_CHAR) / &
                               (LP%VEG_MOIST_MASS*CP_H2O + CP_MASS_VEG_SOLID)
  TMP_VEG_NEW  = MIN(TMP_CHAR_MAX,TMP_VEG_NEW)
  IF (MPV_VEG <= MPV_VEG_MIN) MPV_VOLIT = 0.0_EB
 ENDIF IF_NOT_IGNITOR1
 ENDIF IF_VEG_DEGRADATION_ARRHENIUS

 LP%TMP = TMP_VEG_NEW
 LP%VEG_EMISS = 4.*SIGMA*LP%VEG_KAPPA*LP%TMP**4 !used in RTE solver

 MPV_CHAR_LOSS_TOTAL  = MPV_CHAR_LOSS_TOTAL  + MPV_CHAR_LOSS !needed for subcycling
 MPV_MOIST_LOSS_TOTAL = MPV_MOIST_LOSS_TOTAL + MPV_MOIST_LOSS !needed for subcycling
 MPV_VOLIT_TOTAL      = MPV_VOLIT_TOTAL      + MPV_VOLIT !needed for subcycling
 MPV_CHAR_CO2_TOTAL   = MPV_CHAR_CO2_TOTAL   + MPV_CHAR_CO2
 MPV_ADDED = MPV_ADDED + MPV_MOIST_LOSS + MPV_VOLIT + MPV_CHAR_CO2

! Check if critical mass flux condition is met
!IF (MPV_ADDED*RDT < VEG_CRITICAL_MASSSOURCE .AND. .NOT. LP%VEG_IGNITED) THEN
! MPV_ADDED      = 0.0_EB
! MW_AVERAGE     = 0.0_EB
! MPV_MOIST_LOSS = 0.0_EB
! MPV_VOLIT      = 0.0_EB
! Q_VEG_MOIST    = 0.0_EB
! Q_VEG_VOLIT    = 0.0_EB
!ELSE
! LP%VEG_IGNITED = .TRUE.
!ENDIF

! Add affects of fuel element thermal degradation of vegetation to velocity divergence

 ZZ_GET(1:N_TRACKED_SPECIES) = ZZ(II,JJ,KK,1:N_TRACKED_SPECIES)
 CALL GET_SPECIFIC_HEAT(ZZ_GET,CP_GAS,TMP_GAS)
 RCP_GAS    = 1._EB/CP_GAS
 !MW_TERM    = MW_VEG_MOIST_TERM + MW_VEG_VOLIT_TERM
 MW_AVERAGE = R0/RSUM(II,JJ,KK)/RHO_GAS*(MW_VEG_MOIST_TERM + MW_VEG_VOLIT_TERM)
 Q_ENTHALPY = Q_VEG_MOIST + Q_VEG_VOLIT - (1.0_EB - PC%VEG_CHAR_ENTHALPY_FRACTION)*Q_VEG_CHAR
 !D_LAGRANGIAN(II,JJ,KK) = D_LAGRANGIAN(II,JJ,KK)  +           & 
 !                         (-FCTR_RDT_CYCLES*QCON_VEG*RCP_GAS + Q_ENTHALPY*RCP_GAS)/(RHO_GAS*TMP_GAS) + &
 !                         RDT*MW_AVERAGE 
 D_LAGRANGIAN(II,JJ,KK) = D_LAGRANGIAN(II,JJ,KK)  +           & 
                          (Q_ENTHALPY*RCP_GAS)/(RHO_GAS*TMP_GAS) + &
                          RDT*MW_AVERAGE 


 TMP_VEG   = TMP_VEG_NEW
 IF (MPV_MOIST <= MPV_MOIST_MIN) THEN !for time sub cycling
  MPV_MOIST = 0.0_EB
  MW_VEG_MOIST_TERM = 0.0_EB
  Q_VEG_MOIST = 0.0_EB
 ENDIF
 IF (MPV_VEG <= MPV_VEG_MIN) THEN
  MPV_VOLIT = 0.0_EB
  MPV_VEG   = 0.0_EB
  MW_VEG_VOLIT_TERM = 0.0_EB
  Q_VEG_VOLIT = 0.0_EB
 ENDIF

ENDDO TIME_SUBCYCLING_LOOP

 D_LAGRANGIAN(II,JJ,KK) = D_LAGRANGIAN(II,JJ,KK) + (-QCON_VEG*RCP_GAS)/(RHO_GAS*TMP_GAS)

 IF_NOT_IGNITOR2: IF (.NOT. LP%IGNITOR) THEN !add fuel,H2O,CO2 to mixture factions

! Add water vapor, fuel vapor, and CO2 mass to total density
! MPV_ADDED     = MPV_MOIST_LOSS + MPV_VOLIT + MPV_CHAR_CO2
  LP%VEG_MLR    = MPV_ADDED*RDT !kg/m^3/s used in FVX,FVY,FVZ along with drag in part.f90
  RHO(II,JJ,KK) = RHO_GAS + MPV_ADDED
  RRHO_GAS_NEW  = 1._EB/RHO(II,JJ,KK)
! print*,'NM =',NM
! print*,'** ',rho(ii,jj,kk)

! Add gas species created by degradation of vegetation Yi_new = (Yi_old*rho_old + change in rho_i)/rho_new
! Add water vapor mass from drying to water vapor mass fraction
  IF (I_WATER > 0) THEN 
!  ZZ(II,JJ,KK,I_WATER) = (MPV_MOIST_LOSS + ZZ(II,JJ,KK,I_WATER)*RHO_GAS)/(MPV_MOIST_LOSS + RHO_GAS)
!  ZZ(II,JJ,KK,I_WATER) = ZZ(II,JJ,KK,I_WATER) +  MPV_MOIST_LOSS*RRHO_GAS_NEW
   ZZ(II,JJ,KK,I_WATER) = ZZ(II,JJ,KK,I_WATER) + (MPV_MOIST_LOSS_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_WATER))*RRHO_GAS_NEW
!  ZZ(II,JJ,KK,I_WATER) = MIN(1._EB,ZZ(II,JJ,KK,I_WATER))
!  DMPVDT_FM_VEG(II,JJ,KK,I_WATER) = DMPVDT_FM_VEG(II,JJ,KK,I_WATER) + RDT*MPV_MOIST_LOSS
  ENDIF

! Add fuel vapor mass from pyrolysis to fuel mass fraction
  I_FUEL = REACTION(1)%FUEL_SMIX_INDEX
  IF (I_FUEL /= 0) THEN 
!  ZZ(II,JJ,KK,I_FUEL) = (MPV_VOLIT + ZZ(II,JJ,KK,I_FUEL)*RHO_GAS)/(MPV_VOLIT + RHO_GAS)
!  ZZ(II,JJ,KK,I_FUEL) = ZZ(II,JJ,KK,I_FUEL) + MPV_VOLIT*RRHO_GAS_NEW
   ZZ(II,JJ,KK,I_FUEL) = ZZ(II,JJ,KK,I_FUEL) + (MPV_VOLIT_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_FUEL))*RRHO_GAS_NEW
!  ZZ(II,JJ,KK,I_FUEL) = MIN(1._EB,ZZ(II,JJ,KK,I_FUEL))
!  DMPVDT_FM_VEG(II,JJ,KK,I_FUEL) = DMPVDT_FM_VEG(II,JJ,KK,I_FUEL) + RDT*MPV_VOLIT
  ENDIF

! Add CO2 vapor mass from char oxidation mass to CO2 mass fraction
  IF (I_CO2 /= 0 .AND. PC%VEG_CHAR_OXIDATION) THEN 
   ZZ(II,JJ,KK,I_CO2) = ZZ(II,JJ,KK,I_CO2) + (MPV_CHAR_CO2_TOTAL - MPV_ADDED*ZZ(II,JJ,KK,I_CO2))*RRHO_GAS_NEW
  ENDIF

 ENDIF IF_NOT_IGNITOR2

! WRITE(9998,'(A)')'T,TMP_VEG,QCON_VEG,QRAD_VEG'
!IF (II==0.5*IBAR .AND. JJ==0.5*JBAR .AND. KK==0.333*KBAR) THEN
!IF (II==12 .AND. JJ==12 .AND. KK==4) THEN 
!IF (II==20 .AND. JJ==20 .AND. KK==25) THEN !M=14% and 49% element burnout
!IF (II==27 .AND. JJ==20 .AND. KK==7) THEN !M=49% not full element burnout
! WRITE(9998,'(9(ES12.4))')T,TMP_GAS,TMP_VEG,QCON_VEG,QRAD_VEG,LP%VEG_MOIST_MASS,LP%VEG_FUEL_MASS, &
!                          MPV_MOIST_LOSS_MAX*RDT,MPV_VOLIT_MAX*RDT
!ENDIF

! V_VEG               = V_VEG + V_CELL
! TOTAL_MASS_MOIST    = TOTAL_MASS_MOIST + LP%VEG_MOIST_MASS*V_CELL
! TOTAL_MASS_DRY_FUEL = TOTAL_MASS_DRY_FUEL + LP%VEG_FUEL_MASS*V_CELL

 ENDIF THERMAL_CALC  ! end of thermally thin heat transfer, etc. calculations

 N_TREE = LP%VEG_N_TREE_OUTPUT
 IF (N_TREE /= 0) THEN
  TREE_OUTPUT_DATA(N_TREE,1,NM) = TREE_OUTPUT_DATA(N_TREE,1,NM) + LP%TMP - 273._EB !C
  TREE_OUTPUT_DATA(N_TREE,2,NM) = TREE_OUTPUT_DATA(N_TREE,2,NM) + TMP_GAS - 273._EB !C
  TREE_OUTPUT_DATA(N_TREE,3,NM) = TREE_OUTPUT_DATA(N_TREE,3,NM) + LP%VEG_FUEL_MASS*V_CELL !kg
  TREE_OUTPUT_DATA(N_TREE,4,NM) = TREE_OUTPUT_DATA(N_TREE,4,NM) + LP%VEG_MOIST_MASS*V_CELL !kg
  TREE_OUTPUT_DATA(N_TREE,5,NM) = TREE_OUTPUT_DATA(N_TREE,5,NM) + LP%VEG_CHAR_MASS*V_CELL !kg
  TREE_OUTPUT_DATA(N_TREE,6,NM) = TREE_OUTPUT_DATA(N_TREE,6,NM) + LP%VEG_ASH_MASS*V_CELL !kg
  TREE_OUTPUT_DATA(N_TREE,7,NM) = TREE_OUTPUT_DATA(N_TREE,7,NM) + LP%VEG_DIVQC*V_CELL*0.001_EB !kW
  TREE_OUTPUT_DATA(N_TREE,8,NM) = TREE_OUTPUT_DATA(N_TREE,8,NM) + LP%VEG_DIVQR*V_CELL*0.001_EB !kW
  TREE_OUTPUT_DATA(N_TREE,9,NM) = TREE_OUTPUT_DATA(N_TREE,9,NM) + 1._EB !number of particles
  TREE_OUTPUT_DATA(N_TREE,10,NM) = TREE_OUTPUT_DATA(N_TREE,10,NM) + MPV_CHAR_LOSS_TOTAL*V_CELL !kg 
  TREE_OUTPUT_DATA(N_TREE,11,NM) = TREE_OUTPUT_DATA(N_TREE,11,NM) - Q_VEG_CHAR_TOTAL*V_CELL*RDT*0.001_EB !kW

! TREE_OUTPUT_DATA(N_TREE,4,NM) = TREE_OUTPUT_DATA(N_TREE,4,NM) + LP%VEG_PACKING_RATIO
! TREE_OUTPUT_DATA(N_TREE,5,NM) = TREE_OUTPUT_DATA(N_TREE,5,NM) + LP%VEG_SV

 ENDIF

ENDDO PARTICLE_LOOP

!print*,'--------------------------------'
!print*,'VEGE: NM, TREE_OUTPUT_DATA(1,1,NM),(1,2,NM)',nm, tree_output_data(1,1,nm),tree_output_data(1,2,nm)

! Write out total bulk
!TOTAL_BULKDENS_MOIST = TOTAL_MASS_MOIST/V_VEG
!TOTAL_BULKDENS_DRY_FUEL = TOTAL_MASS_DRY_FUEL/V_VEG
!WRITE(9999,'(5(ES12.4))')T,TOTAL_BULKDENS_DRY_FUEL,TOTAL_BULKDENS_MOIST,TOTAL_MASS_DRY_FUEL,TOTAL_MASS_MOIST

!VEG_TOTAL_DRY_MASS(NM)   = TOTAL_MASS_DRY_FUEL
!VEG_TOTAL_MOIST_MASS(NM) = TOTAL_MASS_MOIST

! Remove vegetation that has completely burned (i.e., LP%R set equal to zero)
CALL REMOVE_PARTICLES(T,NM)
 
END SUBROUTINE RAISED_VEG_MASS_ENERGY_TRANSFER

! ***********************************************************************************************

SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER(T,NM)
!
! Issues:
! 1. Are SF%VEG_FUEL_FLUX_L and SF%VEG_MOIST_FLUX_L needed in linear degradation model?
REAL(EB) :: DT_BC,RDT_BC,T
INTEGER, INTENT(IN) :: NM
INTEGER  ::  IW
INTEGER  ::  I,IIG,JJG,KKG
REAL(EB) :: CP_MOIST_AND_VEG,DZVEG_L,ETAVEG_H,H_CONV_L, &
            KAPPA_VEG,K_AIR,MU_AIR,QRADM_INC,QRADP_INC,RHO_GAS, &
            TMP_BOIL,TMP_G,DTMP_L,RE_H,RE_VEG_PART,U2,V2
!REAL(EB) :: H_CONV_FDS_WALL,DTMP_FDS_WALL,QCONF_FDS_WALL,LAMBDA_AIR,TMPG_A
INTEGER  IIVEG_L,IVEG_L,J,LBURN,NVEG_L,I_FUEL
!REAL(EB), ALLOCATABLE, DIMENSION(:) :: VEG_DIV_QRNET_EMISS,VEG_DIV_QRNET_INC,
!         VEG_QRNET_EMISS,VEG_QRNET_INC,VEG_QRM_EMISS,VEG_QRP_EMISS, VEG_QRM_INC,VEG_QRP_INC
REAL(EB) :: VEG_DIV_QRNET_EMISS(50),VEG_DIV_QRNET_INC(50),VEG_QRNET_EMISS(0:50),VEG_QRNET_INC(0:50), &
            VEG_QRM_EMISS(0:50),VEG_QRP_EMISS(0:50), VEG_QRM_INC(0:50),VEG_QRP_INC(0:50)
REAL(EB) :: A_H2O_VEG,E_H2O_VEG,A_PYR_VEG,E_PYR_VEG,H_PYR_VEG,RH_PYR_VEG
REAL(EB) :: A_CHAR_VEG,E_CHAR_VEG,BETA_CHAR_VEG,NU_CHAR_VEG,NU_ASH_VEG,NU_O2_CHAR_VEG,H_CHAR_VEG 
REAL(EB) :: CP_CHAR,CP_H2O,CP_VEG,CP_TOTAL,DTMP_VEG,H_VAP_H2O,TMP_VEG,TMP_VEG_NEW
REAL(EB) :: CHAR_FCTR,CHAR_FCTR2,MPA_MOIST,MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX,MPA_MOIST_MIN,DMPA_VEG, &
            MPA_CHAR,MPA_VEG,MPA_VEG_MIN,MPA_VOLIT,MPA_VOLIT_LOSS_MAX
REAL(EB) :: DETA_VEG,ETA_H,ETAFM_VEG,ETAFP_VEG
REAL(EB) :: QCONF_L,Q_FOR_DRYING,Q_VEG_MOIST,Q_VEG_VOLIT,QNET_VEG,Q_FOR_VOLIT,Q_VOLIT,Q_UPTO_VOLIT
LOGICAL  :: VEG_DEGRADATION_ARRHENIUS,VEG_DEGRADATION_LINEAR
LOGICAL  :: H_VERT_CYLINDER_LAMINAR,H_CYLINDER_RE
logical  :: fuel_elem_degrad,fds4_degrad

INTEGER  :: IC,II,IOR,JJ,KK,IW_CELL

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

TYPE (WALL_TYPE),    POINTER :: WC1 =>NULL() !to handle qrad on slopes
TYPE (SURFACE_TYPE), POINTER :: SF1 =>NULL() !to handle qrad on slopes

CALL POINT_TO_MESH(NM)

TMP_BOIL  = 373._EB
CP_H2O    = 4190._EB !J/kg/K specific heat of water
H_VAP_H2O = 2259._EB*1000._EB !J/kg/K heat of vaporization of water
H_PYR_VEG = 416000._EB !J/kg Morvan
!H_PYR_VEG = 2640000._EB !J/kg Drysdale,Doug Fir
RH_PYR_VEG = 1._EB/H_PYR_VEG
DT_BC     = T - VEG_CLOCK_BC
!DT_BC     = MESHES(NM)%DT
RDT_BC    = 1.0_EB/DT_BC

IF (N_REACTIONS>0) I_FUEL = REACTION(1)%FUEL_SMIX_INDEX

! Thermal degradation approach parameters
! VEG_DEGRADATION_LINEAR    = .TRUE.
! VEG_DEGRADATION_ARRHENIUS = .FALSE.
  fuel_elem_degrad = .true.
  fds4_degrad      = .false.
!
! Loop through vegetation wall cells and burn
!
VEG_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE VEG_WALL_CELL_LOOP

    SF  => SURFACE(WC%SURF_INDEX)
!
  IF (.NOT. SF%VEGETATION) CYCLE VEG_WALL_CELL_LOOP

  VEG_DEGRADATION_LINEAR    = SF%VEG_LINEAR_DEGRAD
  VEG_DEGRADATION_ARRHENIUS = SF%VEG_ARRHENIUS_DEGRAD
! IF(SF%VEG_NO_BURN) THEN
!  print*,'No Burn Veg'
!  print*,'VEG_DEGRADATION_LINEAR    = ',SF%VEG_LINEAR_DEGRAD
!  print*,'VEG_DEGRADATION_ARRHENIUS = ',SF%VEG_ARRHENIUS_DEGRAD
! ENDIF
! IF(.NOT. SF%VEG_NO_BURN) THEN
!  print*,'BURN VEG'
!  print*,'VEG_DEGRADATION_LINEAR    = ',SF%VEG_LINEAR_DEGRAD
!  print*,'VEG_DEGRADATION_ARRHENIUS = ',SF%VEG_ARRHENIUS_DEGRAD
! STOP
! ENDIF


  IIG = WC%IIG
  JJG = WC%JJG
  KKG = WC%KKG
  TMP_G = TMP(IIG,JJG,KKG)
  CHAR_FCTR  = 1._EB - SF%VEG_CHARFRAC
  CHAR_FCTR2 = 1._EB/CHAR_FCTR
  IF(SF%VEG_NO_BURN) WC%VEG_HEIGHT = SF%VEG_HEIGHT
  VEG_DRAG(IIG,JJG) = SF%VEG_DRAG_INI*(SF%VEG_CHARFRAC + CHAR_FCTR*WC%VEG_HEIGHT/SF%VEG_HEIGHT)

  IF(SF%VEG_NO_BURN) CYCLE VEG_WALL_CELL_LOOP

! Initialize quantities
  Q_VEG_MOIST     = 0.0_EB
  Q_VEG_VOLIT     = 0.0_EB
  Q_UPTO_VOLIT    = 0.0_EB
  Q_VOLIT         = 0.0_EB
  MPA_MOIST_LOSS  = 0.0_EB
  MPA_VOLIT       = 0.0_EB
  SF%VEG_DIVQNET_L = 0.0_EB
  SF%VEG_MOIST_FLUX_L = 0.0_EB
  SF%VEG_FUEL_FLUX_L  = 0.0_EB
  WC%MASSFLUX(I_FUEL) = 0.0_EB 
  WC%QCONF           = 0.0_EB
  IF (I_WATER > 0) WC%MASSFLUX(I_WATER) = 0.0_EB

! Vegetation variables and minimum bounds
  NVEG_L = SF%NVEG_L
  LBURN  = 0
  MPA_VEG_MIN   = 0.001_EB*SF%VEG_LOAD/REAL(NVEG_L,EB) !kg/m^2
  MPA_MOIST_MIN = SF%VEG_MOISTURE*MPA_VEG_MIN !ks/m^2
  IF (SF%VEG_MOISTURE == 0.0_EB) MPA_MOIST_MIN = MPA_VEG_MIN
  DZVEG_L   = SF%VEG_HEIGHT/REAL(NVEG_L,EB)
  KAPPA_VEG = SF%VEG_KAPPA
  DETA_VEG  = DZVEG_L*KAPPA_VEG

! Find top of vegetation which burns downward from the top
  DO IVEG_L = 1,NVEG_L 
    IF(WC%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN) LBURN = IVEG_L
  ENDDO
! LBURN = 0
  WC%VEG_HEIGHT = REAL(NVEG_L-LBURN,EB)*DZVEG_L
! LBURN = 0 !keep charred veg
  !FIRELINE_MLR_MAX = w*R*(1-ChiChar)
  MPA_VOLIT_LOSS_MAX = SF%FIRELINE_MLR_MAX*DT_BC*DZVEG_L 
  MPA_MOIST_LOSS_MAX = MPA_VOLIT_LOSS_MAX
! MPA_MOIST_LOSS_MAX = 9999999._EB
! MPA_VOLIT_LOSS_MAX = 9999999._EB

! Factors for computing divergence of incident and self emission radiant fluxes
! in vegetation fuel bed. These need to be recomputed as the height of the
! vegetation surface layer decreases with burning

! Factors for computing decay of +/- incident fluxes
  SF%VEG_FINCM_RADFCT_L(:) =  0.0_EB
  SF%VEG_FINCP_RADFCT_L(:) =  0.0_EB
! ETA_H = KAPPA_VEG*WC%VEG_HEIGHT
  ETA_H = KAPPA_VEG*REAL(NVEG_L-LBURN,EB)*DZVEG_L
  DO IVEG_L = 0,SF%NVEG_L - LBURN
    ETAFM_VEG = IVEG_L*DETA_VEG
    ETAFP_VEG = ETA_H - ETAFM_VEG
    SF%VEG_FINCM_RADFCT_L(IVEG_L) = EXP(-ETAFM_VEG)
    SF%VEG_FINCP_RADFCT_L(IVEG_L) = EXP(-ETAFP_VEG)
  ENDDO

!  Integrand for computing +/- self emission fluxes
  SF%VEG_SEMISSP_RADFCT_L(:,:) = 0.0_EB
  SF%VEG_SEMISSM_RADFCT_L(:,:) = 0.0_EB
! q+
  DO IIVEG_L = 0,SF%NVEG_L-LBURN !veg grid coordinate
    DO IVEG_L = IIVEG_L,SF%NVEG_L-1-LBURN !integrand index
!    ETAG_VEG = IIVEG_L*DETA_VEG
!    ETAI_VEG =  IVEG_L*DETA_VEG
!    SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-(ETAI_VEG-ETAG_VEG))
     ETAFM_VEG = (IVEG_L-IIVEG_L)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSP_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
  ENDDO
! q-
  DO IIVEG_L = 0,SF%NVEG_L-LBURN
    DO IVEG_L = 1,IIVEG_L
!    ETAG_VEG = IIVEG_L*DETA_VEG
!    ETAI_VEG =  IVEG_L*DETA_VEG
!    SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-(ETAG_VEG-ETAI_VEG))
     ETAFM_VEG = (IIVEG_L-IVEG_L)*DETA_VEG
     ETAFP_VEG = ETAFM_VEG + DETA_VEG
     SF%VEG_SEMISSM_RADFCT_L(IVEG_L,IIVEG_L) = EXP(-ETAFM_VEG) - EXP(-ETAFP_VEG)
    ENDDO
  ENDDO
!
! -----------------------------------------------
! compute CONVECTIVE HEAT FLUX on vegetation
! -----------------------------------------------
   H_VERT_CYLINDER_LAMINAR = .TRUE.
   H_CYLINDER_RE           = .FALSE.
! cylinder heat transfer coefficient, hc, from Albini CST, assumes
! lambda ~ rho*cp*T^1.5/p where cp (of air) is assumed to be 
! independent of temperature. Flux is from Morvan and Dupuy assuming
! constant physical properties and integrating vertically over fuel
! bed to get a factor of h multiplying their qc'''
! DTMP*BETA*sigma*h*hc*(T-Ts)
! hc = 0.350*(sigma/4)*lambda in Albini CST 1985 assumes quiescent air
! hc = 0.683*(sigma/4)*lambda*Re^0.466 ; Re=|u|r/nu, r=2/sigma
!      used by Porterie, cylinders in air flow
! lambda = lambda0*(rho/rho0)(T/T0)^a; a=1.5 below
! TMPG_A     = (TMP_G*0.0033_EB)**1.5
! LAMBDA_AIR = 0.026_EB*RHO(IIG,JJG,KKG)*0.861_EB*TMPG_A
!Albini assumes quiescent air
! H_CONV_L = 0.35*LAMBDA_AIR*SF%VEG_SVRATIO*0.25
!Holman "Heat Transfer",5th Edition, McGraw-Hill, 1981 p.285 
!assumes vertical cylinder laminar air flow
! H_CONV_L = 1.42*(DTMP/VEG_HEIGHT_S(SURF_INDEX))**0.25 !W/m^2/C
!Porterie allow for air flow
! U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
! V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
! RE_VEG_PART = SQRT(U2 + V2)*2./SF%VEG_SVRATIO/TMPG_A/15.11E-6
! H_CONV_L = 0.5*LAMBDA_AIR*0.683*RE_VEG_PART**0.466*0.5*SF%VEG_SVRATIO
!
! DTMP_FDS_WALL   = TMP_G - WALL(IW)%TMP_F
! H_CONV_FDS_WALL = 1.42_EB*(ABS(DTMP_FDS_WALL)/DZVEG_L)**0.25
! QCONF_FDS_WALL  = H_CONV_FDS_WALL*DTMP_FDS_WALL
! QCONF(IW)       = QCONF_FDS_WALL !W/m^2
! print*,'dtmp_fds_wall,qconf',dtmp_fds_wall,qconf(iw)
! print*,'tmp_g,tmp_f(iw)',tmp_g,tmp_f(iw)
! SF%VEG_DIVQNET_L(1) = SF%VEG_PACKING*SF%VEG_SVRATIO*QCONF_L*DZVEG_L !W/m^2
!
! Quantities following fuel element model approach
  IF (H_CYLINDER_RE) THEN
   RHO_GAS  = RHO(IIG,JJG,KKG)
   MU_AIR   = MU_Z(MIN(5000,NINT(TMP_G)),0)*SPECIES_MIXTURE(0)%MW
   K_AIR    = CPOPR*MU_AIR !W/m.K
   U2 = 0.25*(U(IIG,JJG,KKG)+U(IIG-1,JJG,KKG))**2
   V2 = 0.25*(V(IIG,JJG,KKG)+V(IIG,JJG-1,KKG))**2
   RE_VEG_PART = 4._EB*RHO_GAS*SQRT(U2 + V2 + W(IIG,JJG,1)**2)/SF%VEG_SVRATIO/MU_AIR
   RE_H = RE_VEG_PART**0.466_EB
  ENDIF

! Divergence of convective and radiative heat fluxes
!print*,'---- NM=',NM
!print*,rho_gas,qrel,sv_veg,mu_air

  DO I=1,NVEG_L-LBURN
    DTMP_L = TMP_G - WC%VEG_TMP_L(I+LBURN)

!Convective heat correlation for laminar flow (Holman see ref above) 
    IF (H_VERT_CYLINDER_LAMINAR) H_CONV_L = 1.42_EB*(ABS(DTMP_L)/DZVEG_L)**0.25

!Convective heat correlation that accounts for air flow (used by Porteried)
    IF(H_CYLINDER_RE) H_CONV_L = 0.25*0.683*SF%VEG_SVRATIO*K_AIR*RE_H !W/m^2 from Porterie(DeWitt)
!
    QCONF_L  = H_CONV_L*DTMP_L
    SF%VEG_DIVQNET_L(I) = SF%VEG_PACKING*SF%VEG_SVRATIO*QCONF_L*DZVEG_L !W/m^2
  ENDDO
!
  WALL(IW)%QCONF = SUM(SF%VEG_DIVQNET_L) !*RDN(IW)*WC%VEG_HEIGHT
! qconf(iw) = 0.0_EB
!
! -----------------------------------------------
! Compute +/- radiation fluxes and their divergence due to self emission within vegetation
! -----------------------------------------------
  LAYER_RAD_FLUXES: IF (LBURN < NVEG_L) THEN
    VEG_QRP_EMISS   = 0.0_EB ; VEG_QRM_EMISS = 0.0_EB 
    VEG_QRNET_EMISS = 0.0_EB ; VEG_DIV_QRNET_EMISS = 0.0_EB
! qe+
    DO J=0,NVEG_L-LBURN !veg grid coordinate loop
      DO I=J,NVEG_L-LBURN !integrand loop 
         VEG_QRP_EMISS(J) =  VEG_QRP_EMISS(J) + SF%VEG_SEMISSP_RADFCT_L(I,J)*WC%VEG_TMP_L(I+LBURN)**4
      ENDDO
    ENDDO
! qe-
    DO J=0,NVEG_L-LBURN  !veg grid coordinate
      DO I=0,J           !integrand for q-
         VEG_QRM_EMISS(J) = VEG_QRM_EMISS(J) + SF%VEG_SEMISSM_RADFCT_L(I,J)*WC%VEG_TMP_L(I+LBURN)**4
      ENDDO
    ENDDO
    VEG_QRP_EMISS =  VEG_QRP_EMISS*SIGMA
    VEG_QRM_EMISS =  VEG_QRM_EMISS*SIGMA
!
    DO I=0,NVEG_L-LBURN
      VEG_QRNET_EMISS(I) = VEG_QRP_EMISS(I)-VEG_QRM_EMISS(I)
    ENDDO
!    DO I=1,NVEG_L-LBURN
!      VEG_QRNET_EMISS(I)  = VEG_QRNET_EMISS(I) - VEG_QRM_EMISS(I)
!    ENDDO
!
    DO I=1,NVEG_L-LBURN
      VEG_DIV_QRNET_EMISS(I) = VEG_QRNET_EMISS(I-1) - VEG_QRNET_EMISS(I)
    ENDDO
!
! Compute +/- radiation fluxes and their divergence due to incident fluxes on boundaries
    QRADM_INC = WALL(IW)%QRADIN/WALL(IW)%E_WALL !sigma*Ta^4 + flame
!   QRADM_INC = QRADIN(IW)/E_WALL(IW) + SIGMA*TMP_F(IW)**4 ! as done in FDS4
!   print*,'vege: QRADIN(IW)',qradin(iw)

! Adjust incident radiant flux to account for sloped terrain
! assumes user put VEGETATION_NO_BURN=.TRUE. for vertical faces
! sets qrad on cell downspread of vertical face = qrad on cell face upspread of vertical face
!   QRADM_INC = QRADM_INC*1.0038_EB !adjustment for horizontal faces assuming 5 degree slope
!   II = WC%II
!   JJ = WC%JJ
!   KK = WC%KK
!   IC = CELL_INDEX(II-1,JJ,KK)
!   IOR = 1
!   IW_CELL = WALL_INDEX(IC,IOR) 
!   WC1 => WALL(IW_CELL)
!   SF1 => SURFACE(WC1%SURF_INDEX)
!print*,'vege: i,j,k,iw,sf',ii,jj,kk,iw,sf1%veg_no_burn
!   IF(SF1%VEG_NO_BURN) THEN
!print*,'vege: in vertical face qrad determination'
!!   QRADM_INC_SLOPE_VERTFACE = QRADM_INC_SLOPE_VERTFACE + WALL(IW_CELL)%RADIN/WALL(IW_CELL)%E_WALL
!!   QRADM_INC_SLOPE_VERTFACE = QRADM_INC_SLOPE_VERTFACE*0.0872_EB !assumes 5 degree slope
!!   QRADM_INC = QRADM_INC + QRADM_INC_SLOPE_VERTFACE !adjustment for adjacent vertical faces

!   IOR = -3
!   IW_CELL = WALL_INDEX(IC,IOR)
!adjustment for horizontal faces downspread of vertical face
!set flux = to max of flux up or downspread 
!print*,'vege: i,j,k,iw,qr',ii,jj,kk,wall(iw_cell)%qradin,wall(iw)%qradin
!   WALL(IW)%QRADIN = MAX(WALL(IW_CELL)%QRADIN,WALL(IW)%QRADIN) 
!   QRADM_INC = 1.0038_EB*WALL(IW)%QRADIN/WALL(IW_CELL)%E_WALL !assumes 5 degree slope!!!
!print*,'vege: qradm_inc,wallqrad',qradm_inc,wall(iw)%qradin
!   ENDIF

    ETAVEG_H  = (NVEG_L - LBURN)*DETA_VEG
    !this QRADP_INC ensures zero net radiant fluxes at bottom of vegetation
    IF(SF%VEG_GROUND_ZERO_RAD) QRADP_INC = QRADM_INC*SF%VEG_FINCM_RADFCT_L(NVEG_L-LBURN) + VEG_QRM_EMISS(NVEG_L-LBURN)
    !this QRADP_INC assumes the ground stays at user specified temperature
    IF(.NOT. SF%VEG_GROUND_ZERO_RAD) QRADP_INC = SIGMA*SF%VEG_GROUND_TEMP**4
!   QRADP_INC = SIGMA*WC%VEG_TMP_L(NVEG_L)**4*EXP(-ETAVEG_H) + VEG_QRM_EMISS(NVEG_L-LBURN) !fds4
    VEG_QRM_INC   = 0.0_EB ; VEG_QRP_INC = 0.0_EB 
    VEG_QRNET_INC = 0.0_EB ; VEG_DIV_QRNET_INC = 0.0_EB
    DO I=0,NVEG_L-LBURN
      VEG_QRM_INC(I)   = QRADM_INC*SF%VEG_FINCM_RADFCT_L(I)
      VEG_QRP_INC(I)   = QRADP_INC*SF%VEG_FINCP_RADFCT_L(I)
      VEG_QRNET_INC(I) = VEG_QRP_INC(I)-VEG_QRM_INC(I)
    ENDDO
    DO I=1,NVEG_L-LBURN
      VEG_DIV_QRNET_INC(I) = VEG_QRNET_INC(I-1) - VEG_QRNET_INC(I)
    ENDDO
  ENDIF LAYER_RAD_FLUXES
!
! Add divergence of net radiation flux to divergence of convection flux
  DO I=1,NVEG_L-LBURN
    SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - (VEG_DIV_QRNET_INC(I) + VEG_DIV_QRNET_EMISS(I))
!   SF%VEG_DIVQNET_L(I)= SF%VEG_DIVQNET_L(I) - VEG_DIV_QRNET_INC(I)
  ENDDO
!
!
!      ************** Boundary Fuel Non-Arrehnius (Linear in temp) Degradation model *************************
! Drying occurs if qnet > 0 with Tveg held at 100 c
! Pyrolysis occurs according to Morvan & Dupuy empirical formula. Linear
! temperature dependence with qnet factor
!

  IF_VEG_DEGRADATION_LINEAR: IF (VEG_DEGRADATION_LINEAR) THEN

  if_fds4_degrad: if(fds4_degrad) then
!
! compute mass flux of H20 vapor or fuel gas
!
!     VEG_MOIST_FLUX_L   = 0.0
!     VEG_FUEL_FLUX_L    = 0.0
!     WC%MASSFLUX(IFUEL) = 0.0
      DO IVEG_L = LBURN+1,NVEG_L

       IF(SF%VEG_DIVQNET_L(IVEG_L-LBURN) > 0._EB) THEN 

!                                 -- boiling 
        IF(WC%VEG_TMP_L(IVEG_L)>=TMP_BOIL .AND. WC%VEG_MOISTMASS_L(IVEG_L)>MPA_MOIST_MIN) THEN
          SF%VEG_MOIST_FLUX_L(IVEG_L) = SF%VEG_DIVQNET_L(IVEG_L-LBURN)/H_VAP_H2O
          WC%VEG_MOISTMASS_L(IVEG_L) = WC%VEG_MOISTMASS_L(IVEG_L) - DT_BC*SF%VEG_MOIST_FLUX_L(IVEG_L)
!         print*,'&& layer, water mass',
!     .  iveg_l,veg_water_mass_per_area_l(iw,iveg_l)
        ENDIF
!                                 -- pyrolysis multiple layers
        IF (WC%VEG_MOISTMASS_L(IVEG_L)<=MPA_MOIST_MIN .AND. WC%VEG_FUELMASS_L(IVEG_L)>MPA_VEG_MIN) THEN

        IF(WC%VEG_TMP_L(IVEG_L)>= 400._EB .AND. WC%VEG_TMP_L(IVEG_L)<= 500._EB)  &
         SF%VEG_FUEL_FLUX_L(IVEG_L) = SF%VEG_DIVQNET_L(IVEG_L-LBURN)*0.0000025_EB*(WC%VEG_TMP_L(IVEG_L)-400._EB)*0.01_EB 
!
        WC%VEG_FUELMASS_L(IVEG_L) = WC%VEG_FUELMASS_L(IVEG_L) - DT_BC*SF%VEG_FUEL_FLUX_L(IVEG_L)
        WC%MASSFLUX(I_FUEL)= SF%VEG_FUEL_FLUX_L(IVEG_L) + WC%MASSFLUX(I_FUEL)
        ENDIF !pyrolysis models
        WC%VEG_FUELMASS_L(IVEG_L) = MAX(WC%VEG_FUELMASS_L(IVEG_L),MPA_VEG_MIN)
       ENDIF !qnetflux_l > 0
       ENDDO !boil off and pyrolysis
!
! Compute temperature of vegetation
!
      VEG_TEMP_LOOP: DO IVEG_L = LBURN+1,NVEG_L
!
      IF (WC%VEG_MOISTMASS_L(IVEG_L) < MPA_MOIST_MIN) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
!
      CP_VEG = (0.01_EB + 0.0037_EB*WC%VEG_TMP_L(IVEG_L))*1000._EB !W/kg/K
      CP_MOIST_AND_VEG = CP_H2O*WC%VEG_MOISTMASS_L(IVEG_L) + CP_VEG*WC%VEG_FUELMASS_L(IVEG_L)
      WC%VEG_TMP_L(IVEG_L) = WC%VEG_TMP_L(IVEG_L) + DT_BC*( SF%VEG_DIVQNET_L(IVEG_L-LBURN)  &
                             - SF%VEG_MOIST_FLUX_L(IVEG_L)*H_VAP_H2O &
                             - SF%VEG_FUEL_FLUX_L(IVEG_L)*H_PYR_VEG )/CP_MOIST_AND_VEG

!  Set veg. temp to boiling temp if appropriate
      IF (WC%VEG_MOISTMASS_L(IVEG_L)>=MPA_MOIST_MIN .AND. WC%VEG_TMP_L(IVEG_L)>=TMP_BOIL)  &
           WC%VEG_TMP_L(IVEG_L) = TMP_BOIL

!  Set veg. temp to pyroysis temp if appropriate
!     IF (SF%VEG_FUEL_FLUX_L(IVEG_L)>0._EB .AND. WC%VEG_TMP_L(IVEG_L)>500._EB) WC%VEG_TMP_L(IVEG_L) = 500.
      IF (WC%VEG_FUELMASS_L(IVEG_L)>MPA_VEG_MIN .AND. WC%VEG_TMP_L(IVEG_L)>500._EB) WC%VEG_TMP_L(IVEG_L) = 500._EB

      ENDDO VEG_TEMP_LOOP

    WC%VEG_TMP_L(LBURN) = WC%VEG_TMP_L(LBURN+1)

  endif if_fds4_degrad


  IF_FUEL_ELEM_DEGRAD: IF (FUEL_ELEM_DEGRAD) THEN

    LAYER_LOOP1: DO IVEG_L = LBURN+1,NVEG_L
!
! Compute temperature of vegetation
!
      MPA_CHAR    = WC%VEG_CHARMASS_L(IVEG_L)
      MPA_VEG     = WC%VEG_FUELMASS_L(IVEG_L)
      MPA_MOIST   = WC%VEG_MOISTMASS_L(IVEG_L)
      TMP_VEG     = WC%VEG_TMP_L(IVEG_L)
      QNET_VEG    = SF%VEG_DIVQNET_L(IVEG_L-LBURN)
      CP_VEG      = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !J/kg/K
      CP_CHAR     = 420._EB + 2.09_EB*TMP_VEG + 6.85E-4_EB*TMP_VEG**2 !J/kg/K Park etal. C&F 2010 147:481-494
      CP_TOTAL    = CP_H2O*MPA_MOIST +  CP_VEG*MPA_VEG + CP_CHAR*MPA_CHAR
      DTMP_VEG    = DT_BC*QNET_VEG/CP_TOTAL
!     DTMP_VEG    = DT_BC*QNET_VEG/CP_MOIST_AND_VEG
      TMP_VEG_NEW = TMP_VEG + DTMP_VEG 
!     IF(TMP_VEG_NEW >= 800._EB .AND. MPA_VEG <= MPA_VEG_MIN) TMP_VEG_NEW = 800._EB

      IF_DIVQ_L_GE_0: IF(QNET_VEG > 0._EB) THEN 

! -- drying of veg layer
      IF(MPA_MOIST > MPA_MOIST_MIN .AND. TMP_VEG_NEW >= TMP_BOIL) THEN
        Q_FOR_DRYING   = (TMP_VEG_NEW - TMP_BOIL)/DTMP_VEG * QNET_VEG
        MPA_MOIST_LOSS = MIN(DT_BC*Q_FOR_DRYING/H_VAP_H2O,MPA_MOIST_LOSS_MAX)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST-MPA_MOIST_MIN)
        TMP_VEG_NEW    = TMP_BOIL
        WC%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST - MPA_MOIST_LOSS !kg/m^2
        IF( WC%VEG_MOISTMASS_L(IVEG_L) <= MPA_MOIST_MIN ) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
        IF (I_WATER > 0) WC%MASSFLUX(I_WATER) = WC%MASSFLUX(I_WATER) + RDT_BC*MPA_MOIST_LOSS
!       WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW
      ENDIF

! -- pyrolysis multiple layers
      IF_VOLITIZATION: IF (MPA_MOIST <= MPA_MOIST_MIN) THEN

        IF(TMP_VEG_NEW >= 400._EB .AND. MPA_VEG > MPA_VEG_MIN) THEN
          Q_UPTO_VOLIT = (CP_VEG*MPA_VEG + CP_CHAR*MPA_CHAR)*MAX((400._EB-TMP_VEG),0.0_EB)
          Q_FOR_VOLIT  = DT_BC*QNET_VEG - Q_UPTO_VOLIT
!         Q_FOR_VOLIT  = DT_BC*(TMP_VEG_NEW - 400._EB)/DTMP_VEG * QNET_VEG
          Q_VOLIT      = Q_FOR_VOLIT*0.01_EB*(TMP_VEG-400._EB)
!         Q_VOLIT      = Q_FOR_VOLIT

          MPA_VOLIT    = CHAR_FCTR*Q_VOLIT*RH_PYR_VEG
          MPA_VOLIT    = MAX(MPA_VOLIT,0._EB)
          MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX) !user specified max

          DMPA_VEG     = CHAR_FCTR2*MPA_VOLIT
          DMPA_VEG     = MIN(DMPA_VEG,(MPA_VEG-MPA_VEG_MIN))
          MPA_VEG      = MPA_VEG - DMPA_VEG
          MPA_CHAR     = MPA_CHAR + SF%VEG_CHARFRAC*DMPA_VEG !**

          MPA_VOLIT    = CHAR_FCTR*DMPA_VEG
          Q_VOLIT      = MPA_VOLIT*H_PYR_VEG 

          TMP_VEG_NEW  = TMP_VEG + (Q_FOR_VOLIT-Q_VOLIT)/(MPA_VEG*CP_VEG + MPA_CHAR*CP_CHAR)
          TMP_VEG_NEW  = MIN(TMP_VEG_NEW,500._EB)
!         TMP_VEG_NEW  = MIN(TMP_VEG_NEW,900._EB)
          WC%VEG_CHARMASS_L(IVEG_L) = MPA_CHAR
          WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
          IF( WC%VEG_FUELMASS_L(IVEG_L) <= MPA_VEG_MIN ) WC%VEG_FUELMASS_L(IVEG_L) = 0.0_EB !**
          WC%MASSFLUX(I_FUEL)= WC%MASSFLUX(I_FUEL) + RDT_BC*MPA_VOLIT
!         WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW
        ENDIF        

      ENDIF IF_VOLITIZATION

      ENDIF IF_DIVQ_L_GE_0
      
!     IF(MPA_VEG <= MPA_VEG_MIN) TMP_VEG_NEW = TMP_G
      WC%VEG_TMP_L(IVEG_L) = TMP_VEG_NEW

    ENDDO LAYER_LOOP1

!   WC%VEG_TMP_L(LBURN) = WC%VEG_TMP_L(LBURN+1)
    WC%VEG_TMP_L(LBURN) = TMP_G

  ENDIF IF_FUEL_ELEM_DEGRAD

  ENDIF  IF_VEG_DEGRADATION_LINEAR

!      ************** Boundary Fuel Arrehnius Degradation model *************************
! Drying and pyrolysis occur according to Arrehnius expressions obtained 
! from the literature (Porterie et al., Num. Heat Transfer, 47:571-591, 2005
! Predicting wildland fire behavior and emissions using a fine-scale physical
! model

  IF_VEG_DEGRADATION_ARRHENIUS: IF(VEG_DEGRADATION_ARRHENIUS) THEN
    A_H2O_VEG      = 600000._EB !1/s sqrt(K)
    E_H2O_VEG      = 5800._EB !K
    A_PYR_VEG      = 36300._EB !1/s
    E_PYR_VEG      = 7250._EB !K
    H_PYR_VEG      = 418000._EB !J/kg
    A_CHAR_VEG     = 430._EB !m/s
    E_CHAR_VEG     = 9000._EB !K
    BETA_CHAR_VEG  = 0.2_EB
!   NU_CHAR_VEG    = 0.3_EB
    NU_CHAR_VEG    = SF%VEG_CHARFRAC
    NU_ASH_VEG     = 0.1_EB
    NU_O2_CHAR_VEG = 1.65_EB
    H_CHAR_VEG     = -12.0E+6_EB !J/kg

    LAYER_LOOP2: DO IVEG_L = LBURN+1,NVEG_L

      MPA_MOIST = WC%VEG_MOISTMASS_L(IVEG_L)
      MPA_VEG   = WC%VEG_FUELMASS_L(IVEG_L)
      TMP_VEG   = WC%VEG_TMP_L(IVEG_L)

      TEMP_THRESEHOLD: IF (WC%VEG_TMP_L(IVEG_L) > 323._EB) THEN
              !arbitrary thresehold to prevent low-temp hrr reaction
              !added for drainage runs

! Drying of vegetation (Arrhenius)
      IF_DEHYDRATION_2: IF (MPA_MOIST > MPA_MOIST_MIN) THEN
        MPA_MOIST_LOSS = MIN(DT_BC*MPA_MOIST*A_H2O_VEG*EXP(-E_H2O_VEG/TMP_VEG)/SQRT(TMP_VEG), &
                         MPA_MOIST-MPA_MOIST_MIN)
        MPA_MOIST_LOSS = MIN(MPA_MOIST_LOSS,MPA_MOIST_LOSS_MAX) !user specified max
        MPA_MOIST      = MPA_MOIST - MPA_MOIST_LOSS
        WC%VEG_MOISTMASS_L(IVEG_L) = MPA_MOIST !kg/m^2
        IF (MPA_MOIST <= MPA_MOIST_MIN) WC%VEG_MOISTMASS_L(IVEG_L) = 0.0_EB
      ENDIF IF_DEHYDRATION_2

! Volitalization of vegetation(Arrhenius)
      IF_VOLITALIZATION_2: IF(MPA_VEG > MPA_VEG_MIN) THEN
        MPA_VOLIT    = MAX(CHAR_FCTR*DT_BC*MPA_VEG*A_PYR_VEG*EXP(-E_PYR_VEG/TMP_VEG),0._EB)
        MPA_VOLIT    = MIN(MPA_VOLIT,MPA_VOLIT_LOSS_MAX) !user specified max
        MPA_VOLIT    = MIN(MPA_VOLIT,(MPA_VEG-MPA_VEG_MIN))
        MPA_VEG      = MPA_VEG - MPA_VOLIT
        WC%VEG_FUELMASS_L(IVEG_L) = MPA_VEG
      ENDIF IF_VOLITALIZATION_2

      WC%MASSFLUX(I_FUEL)= WC%MASSFLUX(I_FUEL) + MPA_VOLIT*RDT_BC
      IF (I_WATER > 0) WC%MASSFLUX(I_WATER) = WC%MASSFLUX(I_WATER) + MPA_MOIST*RDT_BC

      ENDIF TEMP_THRESEHOLD

! Vegetation temperature (Arrhenius)
      CP_VEG = (0.01_EB + 0.0037_EB*TMP_VEG)*1000._EB !W/kg/K
      CP_MOIST_AND_VEG = CP_H2O*WC%VEG_MOISTMASS_L(IVEG_L) +  CP_VEG*WC%VEG_FUELMASS_L(IVEG_L)

      WC%VEG_TMP_L(IVEG_L) = WC%VEG_TMP_L(IVEG_L) + (DT_BC*SF%VEG_DIVQNET_L(IVEG_L-LBURN) - &
                             (MPA_MOIST_LOSS*H_VAP_H2O + MPA_VOLIT*H_PYR_VEG) )/CP_MOIST_AND_VEG
      WC%VEG_TMP_L(IVEG_L) = MAX( WC%VEG_TMP_L(IVEG_L), TMPA)

    ENDDO LAYER_LOOP2

  ENDIF IF_VEG_DEGRADATION_ARRHENIUS
  
  WC%VEG_TMP_L(LBURN) = MAX(TMP_G,TMPA)
  WC%MASSFLUX_ACTUAL(I_FUEL) = WC%MASSFLUX(I_FUEL)
  IF (I_WATER > 0) WC%MASSFLUX_ACTUAL(I_WATER) = WC%MASSFLUX(I_WATER)
 
! Temperature boundary condtions 
! Mass boundary conditions are determine in subroutine SPECIES_BC in wall.f90 for case SPECIFIED_MASS_FLUX
! TMP_F(IW) = WC%VEG_TMP_L(NVEG_L)
! IF (LBURN < NVEG_L)  TMP_F(IW) = WC%VEG_TMP_L(1+LBURN)
  IF (LBURN < NVEG_L) THEN
    WALL(IW)%TMP_F = WC%VEG_TMP_L(1+LBURN)
!   TMP_F(IW) = ((VEG_QRP_INC(0)+VEG_QRP_EMISS(0))/SIGMA)**.25 !as done in FDS4
  ELSE
    WALL(IW)%TMP_F = MAX(TMP_G,TMPA) !Tveg=Tgas if veg is completely burned
!   TMP_F(IW) = TMPA  !Tveg=Tambient if veg is completely burned
  ENDIF
! TMP_F(IW) = MAX(TMP_F(IW),TMPA)

ENDDO VEG_WALL_CELL_LOOP

VEG_CLOCK_BC = T

END SUBROUTINE BNDRY_VEG_MASS_ENERGY_TRANSFER
!
! ***********************************************************************************************
SUBROUTINE LEVEL_SET_FIRESPREAD(NM)
!
! Level set based modeling of fire spread across terrain. Currently, no computation of the wind field 
! is needed. Instead, U0 and V0 which are specified on the MISC line of the input file, are used for 
! the wind field direction. Does use the extent of the vegetation as defined in the fds input file. 
! Level ground spread rates for the head, flank, and back fires are user defined.
! Spread rate dependence on slope is according to McArthur's rules. 
!
! Issues:
! 1) Need to make level set computation mesh dependent so the the LS slice file
!    is created only where fire is expected
! 2) Need to use multiprocessors
!
!
INTEGER, INTENT(IN) :: NM
INTEGER :: J_FLANK,I,IM1,IM2,IIG,IP1,IP2,IW,J,JJ,JJG,JM1,JP1,LU_SLCF_LS,N_FINAL,N_STEPS_OUT
INTEGER :: IDUM,JDUM,ITER
LOGICAL :: IGNITION = .FALSE.
REAL(EB) :: LX,SR_MAX,T_FINAL,U_AMBIENT,UMAG,UMAX_LS, V_AMBIENT,VMAX_LS,YMIN_LS,YMAX_LS
REAL(EB) :: G_EAST,G_WEST,G_SOUTH,G_NORTH,DT_COEF
REAL(EB) :: HEAD_WIDTH_FCTR,IGNITION_WIDTH_Y,ROS_HEAD1,TIME_LS_LAST
REAL(EB) :: DT_OUTPUT,CPUTIME,LS_T_BEG,LS_T_END,SUMTIME,SUM_T_SLCF,PHI_CHECK
REAL(EB), ALLOCATABLE, DIMENSION(:,:) :: PHI_TEMP

!---Vars for Farsite emulation (phi_w calculation)---
REAL(EB) :: B,C,E,BETA_OP,UMF_TMP,SIGMA,PHX,PHY,MAG_PHI,DZT_DUM, PIO2
REAL(EB) :: PHI_W_X,PHI_W_Y,PHI_S_X,PHI_S_Y,MAG_PHI_S,UMF_X,UMF_Y
!--------------------------------

REAL(EB), ALLOCATABLE, DIMENSION(:) :: X_LS,Y_LS
!
REAL(FB) :: TIME_LS_OUT
!REAL(FB) :: MAG_SR
REAL(FB), ALLOCATABLE, DIMENSION(:,:) :: PHI_OUT, TOA
CHARACTER(30) :: SMOKEVIEW_LABEL,SMOKEVIEW_BAR_LABEL,UNITS

REAL(EB), POINTER, DIMENSION(:,:) :: ZT => NULL()

TYPE (WALL_TYPE),    POINTER :: WC =>NULL()
TYPE (SURFACE_TYPE), POINTER :: SF =>NULL()

CALL CPU_TIME(CPUTIME)
LS_T_BEG = CPUTIME

CALL POINT_TO_MESH(NM)

ZT => LS_Z_TERRAIN

!print*,'level set: z(*)',z
!print*,'level set: ls_z_terrain(1,1)',ls_z_terrain(:,100)
!
!-Initialize variables
!
!-- Domain specification (meters) from input file (assumes there's only one mesh)
!
 LX = XF - XS ; NX_LS = IBAR ; DX_LS = LX/REAL(NX_LS,EB)
 LX = YF - YS ; NY_LS = JBAR ; DY_LS = LX/REAL(NY_LS,EB)
 T_FINAL = T_END
 U_AMBIENT = U0 ; V_AMBIENT = V0
 YMIN_LS = YS
 YMAX_LS = YF
 PIO2 = PI/2
!
! Define spread rate across domain (including no burn areas)
!
ALLOCATE(HEAD_WIDTH(NX_LS,NY_LS))  ; CALL ChkMemErr('VEGE:LEVEL SET','HEAD_WIDTH',IZERO)
ALLOCATE(ROS_HEAD(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_HEAD',IZERO)
ALLOCATE(ROS_FLANK(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_FLANK',IZERO)
ALLOCATE(ROS_BACKU(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','ROS_BACKU',IZERO)
ALLOCATE(FLANKFIRE_LIFETIME(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','FLANKFIRE_LIFETIME',IZERO)
ALLOCATE(WIND_EXP(NX_LS,NY_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','WIND_EXP',IZERO)

!----------wind factor and wind at mean flame height for Farsite emulation -----------
ALLOCATE(PHI_WS(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)
ALLOCATE(PHI_S(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_S',IZERO)
ALLOCATE(PHI_W(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_W',IZERO)
ALLOCATE(UMF(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','UMF',IZERO)
ALLOCATE(THETA_ELPS(NX_LS,NY_LS))    ; CALL ChkMemErr('VEGE:LEVEL SET','THETA_ELPS',IZERO)
!-----------------------------------------------------------------------------------


FLANKFIRE_LIFETIME = 0.0_EB !handles finite length (lifetime) flankfires. For
!                            quenching flanks with lifetimes > TIME_FLANKFIRE_QUENCH
TIME_FLANKFIRE_QUENCH = 20.0_EB !flankfire lifetime in seconds

! Assign spread rates (i.e., vegetation types) to locations on terrain

SUMTIME = 0._EB
SUM_T_SLCF = 0._EB
DT_LS = 0.1_EB
ITER = 0
DT_COEF = 0.5

ROS_HEAD  = 0.0_EB
ROS_HEAD1 = 0.0_EB !needed when dependence on head width is computed
ROS_FLANK = 0.0_EB
ROS_BACKU = 0.0_EB
WIND_EXP  = 1.0_EB
PHI_WS    = 0.0_EB
LSET_ELLIPSE = .FALSE.
THETA_ELPS   = 0.0_EB
LSET_TAN2    = .FALSE.
HEAD_WIDTH   = 1.0_EB

!print*,'surface ros',surface%veg_lset_ros_head
!print*,'surface wind_exp',surface%veg_lset_wind_exp
!
!C_F = 0.2_EB
!
! -- Flux limiter
!LIMITER_LS = 1 !MINMOD
!LIMITER_LS = 2 !SUPERBEE
!LIMITER_LS = 3 !First order upwinding
!
!
LIMITER_LS = FLUX_LIMITER
IF (LIMITER_LS > 3) LIMITER_LS = 1

! -- Output file
DT_OUTPUT  = 0.5_EB
IF (DT_SLCF > 0._EB) DT_OUTPUT = DT_SLCF

TIME_LS    = 0._EB
LU_SLCF_LS = 9999
SMOKEVIEW_LABEL = 'phifield'
SMOKEVIEW_BAR_LABEL = 'phifield'
UNITS  = 'C'
OPEN(LU_SLCF_LS,FILE=TRIM(CHID)//'_lsfs.sf',FORM='UNFORMATTED',STATUS='REPLACE')
WRITE(LU_SLCF_LS) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS) SMOKEVIEW_LABEL(1:30)
WRITE(LU_SLCF_LS) UNITS(1:30)
WRITE(LU_SLCF_LS)1,NX_LS,1,NY_LS,1,1
!
!
! =============== end of case specifications ========================
!
!-- Allocate arrays
ALLOCATE(U_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','U_LS',IZERO) ; U_LS = 0._EB
ALLOCATE(V_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','V_LS',IZERO) ; V_LS = 0._EB
ALLOCATE(PHI_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_LS',IZERO)
ALLOCATE(PHI0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI0_LS',IZERO)
ALLOCATE(PHI1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI1_LS',IZERO)
ALLOCATE(PHI_TEMP(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','PHI_TEMP',IZERO)
ALLOCATE(PHI_OUT(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','PHI_OUT',IZERO) ; PHI_OUT = 0.0
ALLOCATE(SR_X_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_X_LS',IZERO)
ALLOCATE(SR_Y_LS(NX_LS,NY_LS)) ; CALL ChkMemErr('VEGE:LEVEL SET','SR_Y_LS',IZERO)
ALLOCATE(FLUX0_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX0_LS',IZERO)
ALLOCATE(FLUX1_LS(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','FLUX1_LS',IZERO)
ALLOCATE(DZTDX(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTX',IZERO)
ALLOCATE(DZTDY(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','DZDTY',IZERO)
ALLOCATE(MAG_ZT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_ZT',IZERO)
ALLOCATE(MAG_SR_OUT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','MAG_SR_OUT',IZERO)
ALLOCATE(ASPECT(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','ASPECT',IZERO)
!
!-- Computational grid
ALLOCATE(X_LS(NX_LS))   ; CALL ChkMemErr('VEGE:LEVEL SET','X_LS',IZERO)
ALLOCATE(Y_LS(NY_LS+1)) ; CALL ChkMemErr('VEGE:LEVEL SET','Y_LS',IZERO)

!'Time of arrival' grid 02dec11
ALLOCATE(TOA(NX_LS,NY_LS)); CALL ChkMemErr('VEGE:LEVEL SET','LSTOA',IZERO)
TOA = -1.0_EB

DO I = 0,NX_LS-1
!X_LS(I+1) = -0.5_EB*LX + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
 X_LS(I+1) = XS + 0.5_EB*DX_LS + DX_LS*REAL(I,EB)
ENDDO
!
DO J = 0,NY_LS
!Y_LS(J+1) = YMIN_LS + DY_LS*REAL(J,EB)
 Y_LS(J+1) = YS + DY_LS*REAL(J,EB)
ENDDO

!---- Compute components of slope gradient and magnitude of gradient

GRADIENT_ILOOP: DO I = 1,NX_LS
 IM1=I-1 ; IM2=I-2
 IP1=I+1 ; IP2=I+2

 IF (I==1) IM1 = I
 IF (I==NX_LS) IP1 = I
 
 DO J = 1,NY_LS
   JM1=J-1
   JP1=J+1
    
   IF (J==1) JM1 = J
   IF (J==NX_LS) JP1 = J
   
! Changed to GIS-type slope calculation on 5-Feb_2013 -asb
   !G_EAST  = ZT(IP1,JP1) + 2._EB * ZT(IP1,J) + ZT(IP1,JM1) 
   !G_WEST  = ZT(IM1,JP1) + 2._EB * ZT(IM1,J) + ZT(IM1,JM1) 
   !G_NORTH = ZT(IM1,JP1) + 2._EB * ZT(I,JP1) + ZT(IP1,JP1) 
   !G_SOUTH = ZT(IM1,JM1) + 2._EB * ZT(I,JM1) + ZT(IP1,JM1) 
   !
   !DZTDX(I,J) = (G_EAST-G_WEST) / (8._EB*DX_LS)
   !DZTDY(I,J) = (G_NORTH-G_SOUTH) / (8._EB*DY_LS)
   !
   
   G_EAST  = 0.5*( ZT(I,J) + ZT(IP1,J) )
   G_WEST  = 0.5*( ZT(I,J) + ZT(IM1,J) )
   G_NORTH = 0.5*( ZT(I,J) + ZT(I,JP1) )
   G_SOUTH = 0.5*( ZT(I,J) + ZT(I,JM1) )

   DZTDX(I,J) = (G_EAST-G_WEST) / DX_LS
   DZTDY(I,J) = (G_NORTH-G_SOUTH) / DY_LS


   MAG_ZT(I,J) = SQRT(DZTDX(I,J)**2 + DZTDY(I,J)**2)
   
   ASPECT(I,J) = PIO2 - ATAN2(-DZTDY(I,J),-DZTDX(I,J)) 
   IF (ASPECT(I,J) < 0.0_EB) ASPECT(I,J) = 2*PI + ASPECT(I,J)
    
        


 ENDDO

ENDDO GRADIENT_ILOOP

!
!-- Build wind field ------------------------------------------------------
!U_LS = U(1:NX_LS,1:NY_LS,1) ; use this if computed wind is to be used
!V_LS = V(1:NX_LS,1:NY_LS,1)

U_LS = U0
V_LS = V0

!U_LS = U_LS + U_AMBIENT
!V_LS = V_LS + V_AMBIENT
!
! Compute time step
UMAX_LS  = MAXVAL(ABS(U_LS))
VMAX_LS  = MAXVAL(ABS(V_LS))
UMAG     = SQRT(UMAX_LS**2 + VMAX_LS**2)



ROS_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE ROS_WALL_CELL_LOOP

  SF  => SURFACE(WC%SURF_INDEX)

  IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE ROS_WALL_CELL_LOOP



  IIG = WC%IIG
  JJG = WC%JJG
!print*,'IIG,JJG',iig,jjg
!print*,'ROS_HEAD',SF%VEG_LSET_ROS_HEAD
!print*,'ROS_HEAD,ROS_FLANK,ROS_BACK',SF%VEG_LSET_ROS_HEAD,SF%VEG_LSET_ROS_FLANK,SF%VEG_LSET_ROS_BACK

!print*,'IIG,JJG',iig,jjg
!print*,'ROS_HEAD',SF%VEG_LSET_ROS_HEAD

  HEAD_WIDTH(IIG,JJG)= DX_LS
  ROS_HEAD(IIG,JJG)  = SF%VEG_LSET_ROS_HEAD
  ROS_FLANK(IIG,JJG) = SF%VEG_LSET_ROS_FLANK
  ROS_BACKU(IIG,JJG) = SF%VEG_LSET_ROS_BACK
  WIND_EXP(IIG,JJG)  = SF%VEG_LSET_WIND_EXP
  
  !If any surfaces uses tan^2 function for slope, tan^2 will be used throughout simulation
  IF (SF%VEG_LSET_TAN2) LSET_TAN2=.TRUE.

  
  IF (SF%VEG_LSET_ELLIPSE) THEN    
    
      ROS_HEAD(IIG,JJG) = SF%VEG_LSET_ELLIPSE_HEAD 
      !If any surfaces set to ellipse, then elliptical model used for all calculations
      IF (.NOT. LSET_ELLIPSE) LSET_ELLIPSE=.TRUE.
    
      SF%VEG_LSET_HT = MAX(0.001,SF%VEG_LSET_HT)
      
    !Wind at midflame height (UMF). Factor 60 converts to m/min used in Farsite. 
    !UMF_TMP = 1._EB / LOG((20.0_EB + 1.18 * SF%VEG_LSET_HT) /(0.43 * SF%VEG_LSET_HT))
    
    ! From Andrews 2012, USDA FS Gen Tech Rep. RMRS-GTR-266 (with added SI conversion)
    UMF_TMP = 1.83_EB / LOG((20.0_EB + 1.18 * SF%VEG_LSET_HT) /(0.43 * SF%VEG_LSET_HT))

    UMF_X = UMF_TMP * U_LS(IIG,JJG) * 60.0_EB
    UMF_Y = UMF_TMP * V_LS(IIG,JJG) * 60.0_EB
  
     SIGMA = SF%VEG_LSET_SIGMA ! SAV ratio, converted to 1/cm in read.f90 
     !Variables used in Phi_W formulas below
     B = 0.15988 * (SIGMA**0.54)
     C = 7.47 * EXP(-0.8711 * (SIGMA**0.55))
     E = 0.715 * EXP(-0.01094 * SIGMA)
     BETA_OP = 0.20395 * (SIGMA**(-0.8189))! Optimum packing ratio
     
    
     
     ! Components of wind factor - affects spread rate
     PHI_W_X = C * ((3.281 * ABS(UMF_X))**B) * (SF%VEG_LSET_BETA / BETA_OP)**(-E)
     PHI_W_X = SIGN(PHI_W_X,UMF_X)
     
     PHI_W_Y = C * ((3.281 * ABS(UMF_Y))**B) * (SF%VEG_LSET_BETA / BETA_OP)**(-E)
     PHI_W_Y = SIGN(PHI_W_Y,UMF_Y)

     PHI_W(IIG,JJG) =  SQRT(PHI_W_X**2 + PHI_W_Y**2)
     
     !Limit effect to slope lte 60 degrees
     !Phi_s_x,y are slope factors
     DZT_DUM = MIN(1.73,ABS(DZTDX(IIG,JJG))) ! 1.73 ~ tan 60 deg
     !DZT_DUM = SIGN(DZT_DUM,DZTDX)
     PHI_S_X = 5.275 * ((SF%VEG_LSET_BETA)**(-0.3)) * DZT_DUM**2
     PHI_S_X = SIGN(PHI_S_X,DZTDX(IIG,JJG))
     
     DZT_DUM = MIN(1.73,ABS(DZTDY(IIG,JJG))) ! 1.73 ~ tan 60 deg
     !DZT_DUM = SIGN(DZT_DUM,DZTDY)
     PHI_S_Y = 5.275 * ((SF%VEG_LSET_BETA)**(-0.3)) * DZT_DUM**2
     PHI_S_Y = SIGN(PHI_S_Y,DZTDY(IIG,JJG))
     
     MAG_PHI_S = SQRT(PHI_S_X**2 + PHI_S_Y**2)
     
     PHI_S(IIG,JJG) = MAG_PHI_S  !5.275 * MAG_ZT(I,J)**2 * (SF%VEG_LSET_BETA)**-0.3
     
     
     IF (MAG_PHI_S > 0.0_EB) THEN      
        
        PHX = PHI_W_X + PHI_S_X
        PHY = PHI_W_Y + PHI_S_Y
        MAG_PHI = SQRT(PHX**2 + PHY**2)
        
        !Total phi (phi_w + phi_s) for use in spread rate section
        PHI_WS(IIG,JJG) = MAG_PHI
        
        !Theta_elps is angle of direction (0 to 2pi) of highest spread rate
        !0<=theta_elps<=2pi as measured clockwise from Y-axis 
        THETA_ELPS(IIG,JJG) = ATAN2(PHY,PHX)
        
        !"Effective midflame windspeed" used in length-to-breadth ratio calculation (spread rate routine)
        ! is the wind + slope effect obtained by solving Phi_w eqs. above for UMF
        UMF(IIG,JJG) = (((MAG_PHI * (SF%VEG_LSET_BETA / BETA_OP)**E)/C)**(1/B))/3.281
        
       !IF (MAG_ZT(I,J)>0._EB) THEN
        ! Adjust phi_s according to aspect and direction of highest ros (theta_elps)   
       ! PHI_S(I,J) = PHI_S(I,J) * COS(ASPECT(I,J) - PI - THETA_ELPS(I,J)) 
       ! PHI_S(I,J) = SIGN(PHI_S(I,J),ASPECT(I,J)-THETA_ELPS(I,J))
       !ENDIF
       
     
     ELSE
     
        PHI_WS(IIG,JJG) = SQRT (PHI_W_X**2 + PHI_W_Y**2)
        !IF (PHY == 0._EB) PHY = 1.E-6_EB
        !0<= Theta_elps <=2pi as measured clockwise from Y-axis 
        THETA_ELPS(IIG,JJG) = ATAN2(PHI_W_Y,PHI_W_X)    
        UMF(IIG,JJG) = SQRT(UMF_X**2 + UMF_Y**2)
        
     ENDIF
       
     !The following two lines convert ATAN2 output to compass system (0 to 2 pi CW from +Y-axis)
     THETA_ELPS(IIG,JJG) = PIO2 - THETA_ELPS(IIG,JJG)
     IF (THETA_ELPS(IIG,JJG) < 0.0_EB) THETA_ELPS(IIG,JJG) = 2.0*PI + THETA_ELPS(IIG,JJG)
  
  ENDIF
  
ENDDO ROS_WALL_CELL_LOOP

!print*,'before assign ROS'
print*,'ROS_HEAD max',MAXVAL(ROS_HEAD)
ROS_HEAD1 = MAXVAL(ROS_HEAD)
print*,'ROS_HEAD1',ROS_HEAD1



!SR_X_MAX = ( MAXVAL(ROS_HEAD_U0_INFW) + MAXVAL(ROS_HEAD_U_INFW) )*UMAX_LS
!SR_Y_MAX = ( MAXVAL(ROS_HEAD_U0_INFW) + MAXVAL(ROS_HEAD_U_INFW) )*VMAX_LS
!SR_MAX   = MAX(SR_Y_MAX,SR_X_MAX)
!SR_MAX   = MAX(SR_MAX,ROS_HEADS)
!SR_X_MAX =  MAXVAL(ROS_HEAD)
SR_MAX   = MAXVAL(ROS_HEAD)
DYN_SR_MAX = 0._EB

!print*,'SR_MAX',sr_max
SR_MAX   = MAX(SR_MAX,MAXVAL(ROS_FLANK))

!Flank rate not available when ellipse/farsite model is used
IF (LSET_ELLIPSE) THEN
    SR_MAX   = MAXVAL(ROS_HEAD) * (1._EB + MAXVAL(PHI_S) + MAXVAL(PHI_W)) ! + MAXVAL(PHI_S))
ENDIF
print*,'Phi_S max',MAXVAL(PHI_S)
print*,'Phi_W max',MAXVAL(PHI_W)
print*,'UMF max',MAXVAL(UMF)
print*,'Mag_zt max',MAXVAL(MAG_ZT)
print*,'SR_MAX',sr_max
IF (.NOT. LSET_ELLIPSE) SR_MAX   = 2._EB*SR_MAX !rough accounting for upslope spread aligned with wind
!DT_LS    = MIN(DX_LS,DY_LS)/SR_MAX
DT_LS = 0.25_EB*MIN(DX_LS,DY_LS)/SR_MAX
!DT_LS = 0.1603_EB !to make AU F19 ignition sequence work
N_STEPS_OUT = DT_OUTPUT/DT_LS
N_STEPS_OUT = MAX(N_STEPS_OUT,1)
N_FINAL     = MAX(T_FINAL/DT_LS,1._EB)
print*,'vege: t_final,dt_ls',t_final,dt_ls
print*,'flux limiter= ',LIMITER_LS
!
!
!-- Initialize level set field. Fireline is at PHI_LS = 0 -------------------
!
!
!--- Finite width ignition line
PHI_MIN_LS = -1._EB
PHI_MAX_LS = 1._EB
PHI_LS = PHI_MIN_LS
!
!
!-- Time step solution using second order Runge-Kutta -----------------------
!

print*,'vege: n_final',n_final

DO WHILE (TIME_LS < T_FINAL)

    ITER = ITER + 1
!
!-- Find flank-to-flank distance at base of fire assume symmetry about ymid and
!   define spread rate based on AU head fire width dependence
!IF (.NOT. LSET_ELLIPSE) THEN

!  ROS_WALL_CELL_LOOP2: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
!    WC  => WALL(IW)
!    IF (WC%BOUNDARY_TYPE==NULL_BOUNDARY) CYCLE ROS_WALL_CELL_LOOP2
!    SF  => SURFACE(WC%SURF_INDEX)
!    IF (.NOT. SF%VEG_LSET_SPREAD) CYCLE ROS_WALL_CELL_LOOP2
!    IF (.NOT. SF%VEG_LSET_HEADWIDTH_DEPENDENCE) CYCLE ROS_WALL_CELL_LOOP2 
!    I = WC%IIG
!    J = WC%JJG
!Case C064     
!    IF(TIME_LS > 0._EB .AND. TIME_LS < 27._EB)  HEAD_WIDTH(I,J) = 2._EB*0.9_EB*TIME_LS !Ignition procdure 0.9 m/s rate
!    IF(TIME_LS >= 27._EB) HEAD_WIDTH(I,J) = HEAD_WIDTH(I,J) + 2._EB*ROS_FLANK(I,J)*(TIME_LS-TIME_LS_LAST)
!Case F19
!During Field Ignition procedure 1.54 m/s rate
!    IF(TIME_LS > 0._EB .AND. TIME_LS < 57._EB)  HEAD_WIDTH(I,J) = 2._EB*1.54_EB*TIME_LS 
!After field ignition procedure
!    IF(TIME_LS >= 57._EB .AND. TIME_LS < 100._EB) &
!                              HEAD_WIDTH(I,J) = HEAD_WIDTH(I,J) + 2._EB*ROS_FLANK(I,J)*(TIME_LS-TIME_LS_LAST)
!    IF(TIME_LS >= 100._EB) HEAD_WIDTH(I,J) = 100000._EB
!    HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH(I,J))
!    IF(ROS_HEAD(I,J) > 0.0_EB) ROS_HEAD(I,J)=ROS_HEAD1*HEAD_WIDTH_FCTR

!   ENDDO ROS_WALL_CELL_LOOP2
!ENDIF
     TIME_LS_LAST = TIME_LS


!    IF(TIME_LS > 0._EB .AND. TIME_LS < 24._EB)  HEAD_WIDTH = 2._EB*1._EB*TIME_LS !Ignition procdure 1 m/s rate
!    IF(TIME_LS >= 24._EB) HEAD_WIDTH = HEAD_WIDTH + 2._EB*ROS_FLANK1*(TIME_LS - TIME_LS_LAST)
!    TIME_LS_LAST = TIME_LS
!    HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH)
!     DO J = 1,NY_LS
!      DO I = 1,NX_LS
!       IF(ROS_HEAD(I,J) > 0.0_EB) ROS_HEAD(I,J)=1.48*HEAD_WIDTH_FCTR
!      ENDDO
!     ENDDO

   IF (SF%VEG_LSET_HEADWIDTH_DEPENDENCE) THEN 
     IGNITION_WIDTH_Y = 8
     J_FLANK = 0
     DO JJ = NY_LS/2,NY_LS
   !  IF(PHI_LS(26,JJ) <= 0.0_EB .AND. J_FLANK==0) J_FLANK = JJ
      IF(PHI_LS(19,JJ) > 0.0_EB) J_FLANK = J_FLANK + 1
     ENDDO
!    HEAD_WIDTHV = 2.0_EB*J_FLANK*DY_LS
!    IF (HEAD_WIDTHV < IGNITION_WIDTH_Y) HEAD_WIDTHV = IGNITION_WIDTH_Y
!    HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTHV)
     DO J = 1,NY_LS
      DO I = 1,NX_LS
       HEAD_WIDTH(I,J) = HEAD_WIDTH(I,J) + 2._EB*ROS_FLANK(I,J)*(TIME_LS-TIME_LS_LAST)
       HEAD_WIDTH_FCTR = EXP(-(0.859_EB + 2.036_EB*UMAG)/HEAD_WIDTH(I,J))
       IF(ROS_HEAD(I,J) > 0.0_EB) ROS_HEAD(I,J)=ROS_HEAD1*HEAD_WIDTH_FCTR
      ENDDO
     ENDDO
   ENDIF

IF (((TIME_LS<=10.0) .OR. (SUMTIME > 100._EB)) .OR. (SUMTIME > 10._EB**FLOOR(LOG10(TIME_LS)))) THEN
 SUMTIME = 0._EB
 print*,'vege:LS:-------------------------------------'
 print*,'vege:LS:time_ls',time_ls
 print*,'vege:ls:dt',dt_ls
 print*,'vege:LS:HW,ros_h',head_width(nx_ls/2,ny_ls/2),ros_head(nx_ls/2,ny_ls/2)
 print*,'vege:LS:ros_f',ros_flank(nx_ls/2,ny_ls/2)
 print*,'vege:LS:max_ros',dyn_sr_max
ENDIF

! -- Ignite landscape at user specified location(s) and time(s)
! ** change to go into the wall cell loop only if time corresponds to an
! ** ignition time. Need to create a separate array with sorted ignition times

IGNITOR_WALL_CELL_LOOP: DO IW=1,N_EXTERNAL_WALL_CELLS+N_INTERNAL_WALL_CELLS
  WC  => WALL(IW)
  SF  => SURFACE(WC%SURF_INDEX)


  IIG = WC%IIG
  JJG = WC%JJG

  IF (SF%VEG_LSET_IGNITE_T >= TIME_LS .AND. SF%VEG_LSET_IGNITE_T <= TIME_LS + DT_LS) PHI_LS(IIG,JJG) = PHI_MAX_LS 

ENDDO IGNITOR_WALL_CELL_LOOP

IF (ANY(PHI_LS==PHI_MAX_LS)) IGNITION = .TRUE.

  
  !DT_LS = 0.3_eb ! Override
!
PHI_TEMP = PHI_LS 
!--- RK Stage 1
 RK2_PREDICTOR_LS = .TRUE.
 CALL LEVEL_SET_SPREAD_RATE
 CALL LEVEL_SET_ADVECT_FLUX
 PHI1_LS = PHI_LS - DT_LS*FLUX0_LS
 !PHI1_LS = MAX(PHI1_LS,PHI_MIN_LS)
 !PHI1_LS = MIN(PHI1_LS,PHI_MAX_LS)

!--- RK Stage2
 RK2_PREDICTOR_LS = .FALSE.
 MAG_SR_OUT       = 0.0
 CALL LEVEL_SET_SPREAD_RATE 
 CALL LEVEL_SET_ADVECT_FLUX
 PHI_LS = PHI_LS - 0.5_EB*DT_LS*(FLUX0_LS + FLUX1_LS)
 !PHI_LS = MAX(PHI_LS,PHI_MIN_LS)
 !PHI_LS = MIN(PHI_LS,PHI_MAX_LS)
 

 ! **********************************************************
 ! **************** Flexible time step **********************
 
  PHI_CHECK = MAXVAL(ABS(PHI_LS - PHI_TEMP))
 
 !write(*,*)"Phi check ",phi_check
 !write(*,*)"dt_coef",dt_coef

 IF (IGNITION) THEN
     ! If any phi values change by more than 0.5, or all change less
     ! than 0.1 (both values are arbitrary), then adjust time step
     
     IF (PHI_CHECK > 0.5) THEN
         ! Cut time step in half and cycle the do-while loop
         DT_COEF = 0.5 * DT_COEF 
         DT_LS = DT_COEF * MIN(DX_LS,DY_LS)/DYN_SR_MAX
         DT_LS = MIN(DT_LS,100._EB)
         PHI_LS = PHI_TEMP ! Exchange for previous phi and cycle
         WRITE(*,*)"Halving time step."
         CYCLE 
     ENDIF
 
     ! Increas time step by 1/4 if changes are small
     IF (PHI_CHECK < 0.1) DT_COEF = DT_COEF * 1.25
     
     DYN_SR_MAX = MAX(DYN_SR_MAX,0.01) ! dyn_sr_max must be g.t. zero
     DT_LS = DT_COEF * MIN(DX_LS,DY_LS)/DYN_SR_MAX
     DT_LS = MIN(DT_LS,100._EB)
     
 ENDIF
 
 ! **********************************************************



 !
 !Construct Time of Arrival (TOA) grid matching the horizontal domain grid
 !Phi_ls goes from -1 to 1
 DO IDUM=1,NX_LS
     DO JDUM=1,NY_LS
         IF (PHI_LS(IDUM,JDUM)>=0._EB .AND. TOA(IDUM,JDUM)<=-1.0_EB) TOA(IDUM,JDUM)=TIME_LS
     ENDDO
 ENDDO
 

 TIME_LS = TIME_LS + DT_LS
 SUMTIME = SUMTIME + DT_LS
 SUM_T_SLCF = SUM_T_SLCF + DT_LS

 !--- Output slice file for smokeview
 !IF (MOD(N_TIME,N_STEPS_OUT) < 0.0001_EB) THEN
 IF (SUM_T_SLCF >= DT_OUTPUT) THEN    
  SUM_T_SLCF = 0._EB
  PHI_OUT = PHI_LS
  TIME_LS_OUT = TIME_LS
  WRITE(LU_SLCF_LS) TIME_LS_OUT
!negative for consistency with wall thickness output from wfds
  WRITE(LU_SLCF_LS) ((-PHI_OUT(I,J),I=1,NX_LS),J=1,NY_LS) 

! output magnitude of spread rate
!MAG_SR_OUT = 0.0
!DO I = 1,NX_LS
! DO J = 1,NY_LS
!!   MAG_SR = SQRT(SR_X_LS(I,J)**2 + SR_Y_LS(I,J)**2)
!!   IF(MAG_SR  > PHI_OUT(I,J)) PHI_OUT(I,J) = MAG_SR
!!   PHI_OUT(I,J) = MAG_SR
!  IF(PHI_LS(I,J) <= 0.2_EB .AND. PHI_LS(I,J) >= -0.2) MAG_SR_OUT(I,J) =-FLUX1_LS(I,J)
!  IF(PHI_OUT(I,J) < MAG_SR_OUT(I,J)) PHI_OUT(I,J) = MAG_SR_OUT(I,J)
! ENDDO
!ENDDO
!! WRITE(LU_SLCF_LS) ((PHI_OUT(I,J),I=1,NX_LS),J=1,NY_LS) 
!! PHI_OUT = SR_X_LS + SR_Y_LS
! WRITE(LU_SLCF_LS) ((MAG_SR_OUT(I,J),I=1,NX_LS),J=1,NY_LS) 

 ENDIF
!
ENDDO !While loop

CALL CPU_TIME(CPUTIME)
LS_T_END = CPUTIME
print*,'CPU Time: ',LS_T_END - LS_T_BEG

! ******  Write array to file **************

OPEN(9998,FILE='time_of_arrival.toa',STATUS='REPLACE')
WRITE(9998,'(I5)') NX_LS,NY_LS
WRITE(9998,'(F7.2)') XS,XF,YS,YF
!Write across row (TOA(1,1), TOA(1,2), ...) to match Farsite output
WRITE(9998,'(F7.2)') ((TOA(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
CLOSE(9998)

!OPEN(9998,FILE='Phi_S.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((PHI_S(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
!!
!OPEN(9998,FILE='Phi_W.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((PHI_W(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
!
!OPEN(9998,FILE='alt.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F10.5)') ((ZT(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
!
!OPEN(9998,FILE='DZTDX.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((DZTDX(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
!
!OPEN(9998,FILE='DZTDY.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((DZTDY(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
!
!OPEN(9998,FILE='Theta_Ellipse.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((Theta_Elps(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
!
!OPEN(9998,FILE='UMF.txt',STATUS='REPLACE')
!WRITE(9998,'(I5)') NX_LS,NY_LS
!WRITE(9998,'(F7.2)') XS,XF,YS,YF
!WRITE(9998,'(F7.2)') ((UMF(IDUM,JDUM),JDUM=1,NY_LS),IDUM=1,NX_LS)
!CLOSE(9998)
! ************************************************************

!
CLOSE(LU_SLCF_LS)
!
END SUBROUTINE LEVEL_SET_FIRESPREAD


SUBROUTINE LEVEL_SET_SPREAD_RATE
!
! Compute components of spread rate vector
!
INTEGER :: I,J,IM1,IP1,JM1,JP1
REAL(EB) :: COS_THETA_WIND,COS_THETA_SLOPE,COS_THETA_WIND_H,COS_THETA_WIND_B, &
            COS_THETA_SLOPE_H,COS_THETA_SLOPE_B,DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH, &
            GRAD_SLOPE_DOT_NORMAL_FIRELINE,MAG_F,MAG_SR,MAG_U,WIND_DOT_NORMAL_FIRELINE,NEXP_WIND
REAL(EB) :: RAD_TO_DEGREE,DEGREES_SLOPE,SLOPE_FACTOR
!Variables for Farsite emulation
REAL(EB) :: COS_THETA,SIN_THETA,XSF,YSF,UMF_DUM,PIO2
REAL(EB) :: A_ELPS,A_ELPS2,B_ELPS2,B_ELPS,C_ELPS,DENOM,ROS_TMP,LB,LBD,HB
REAL(EB), DIMENSION(:)   :: NORMAL_FIRELINE(2)
 
RAD_TO_DEGREE = 90._EB/ASIN(1._EB)
PIO2 = PI/2
!NEXP_WIND = 2


IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS
SR_X_LS = 0.0_EB ; SR_Y_LS = 0.0_EB
DYN_SR_MAX = 0.0_EB

FLUX_ILOOP: DO I = 1,NX_LS
!
 !IM1=I-1; IF (IM1<1) IM1=IM1+NX_LS
 !IM2=I-2; IF (IM2<1) IM2=IM2+NX_LS
 !
 !IP1=I+1; IF (IP1>NX_LS) IP1=IP1-NX_LS
 !IP2=I+2; IF (IP2>NX_LS) IP2=IP2-NX_LS
  
  IM1=I-1 
  IP1=I+1 
  IF (I==1) IM1 = I
  IF (I==NX_LS) IP1 = I
  
  DO J = 1,NY_LS
    
    JM1=J-1
    JP1=J+1
    IF (J==1) JM1 = J
    IF (J==NX_LS) JP1 = J
  
      !IF (J==1) THEN
      !  JP1 = J+1     
      !  F_EAST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
      !  F_WEST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
      !  F_NORTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
      !  F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI_MIN_LS ) 
      !ELSEIF (J==NY_LS) THEN
      !  JM1 = J-1
      !  F_EAST =  0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
      !  F_WEST =  0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
      !  F_NORTH = 0.5*( PHI0_LS(I,J) + PHI_MIN_LS )
      !  F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )
      !ELSE
      !  JM1=J-1
      !  JP1=J+1
      !  F_EAST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
      !  F_WEST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
      !  F_NORTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
      !  F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )
      !
      !ENDIF

    F_EAST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IP1,J) )
    F_WEST  = 0.5*( PHI0_LS(I,J) + PHI0_LS(IM1,J) )
    F_NORTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JP1) )
    F_SOUTH = 0.5*( PHI0_LS(I,J) + PHI0_LS(I,JM1) )
       

   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   
   !------------------------------------------------------------------------
   !    The two loops below check for rare cases of symmetrically merging fire lines
   !     where the central difference may be zero but there are forward or 
   !     backward differences
   !------------------------------------------------------------------------
   IF ((DPHIDX == 0._EB) .AND. (F_EAST > 0._EB)) THEN
        IF ((F_EAST == F_WEST) .AND. (PHI0_LS(I,J) < 0._EB)) THEN
            DPHIDX =  (PHI0_LS(IP1,J) - PHI0_LS(I,J))/DX_LS
        ENDIF
   ENDIF
   
   IF ((DPHIDY == 0._EB) .AND. (F_NORTH > 0._EB)) THEN
        IF ((F_NORTH == F_SOUTH) .AND. (PHI0_LS(I,J) < 0._EB)) THEN
            DPHIDY = (PHI0_LS(I,JP1) - PHI0_LS(I,J) )/DY_LS
        ENDIF
   ENDIF
   
   MAG_F = SQRT(DPHIDX**2 + DPHIDY**2)
   IF (MAG_F > 0._EB) THEN   !components of unit vector normal to PHI contours
        NORMAL_FIRELINE(1) = -DPHIDX/MAG_F
        NORMAL_FIRELINE(2) = -DPHIDY/MAG_F
        ! Eulerian normal approximation from Rehm and Mcdermott 2009
        ! For elliptical front calculation
        XSF = (DPHIDY / MAG_F ) !* DX_LS 
        YSF = (-DPHIDX / MAG_F) !* DY_LS    
   ELSE
        NORMAL_FIRELINE = 0._EB
        XSF=0._EB
        YSF=0._EB
   ENDIF

   GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*(DPHIDY/MAG_F) + DZTDY(I,J)*(-DPHIDY/MAG_F) ! ???
   COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB
   
   IF (MAG_ZT(I,J) > 0.0_EB) THEN
       COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)
       !XSF = XSF * COS_THETA_SLOPE
       !YSF = YSF * COS_THETA_SLOPE
   ENDIF
   
   DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE
    
   
   IF (LSET_ELLIPSE) THEN
       
       
       ! Effective wind direction (theta) is clockwise from y-axis (Richards 1990)
       COS_THETA = COS(THETA_ELPS(I,J)) !V_LS(I,J) / MAG_U
       SIN_THETA = SIN(THETA_ELPS(I,J)) !U_LS(I,J) / MAG_U

       !ROS_TMP = ROS_HEAD(I,J) * (1.0_EB + PHI_S(I,J) + PHI_W(I,J))
       ROS_TMP = ROS_HEAD(I,J) * (1.0_EB + PHI_WS(I,J))
       
       !Mag of wind speed at midflame ht must be in units of m/s here   
       UMF_DUM = UMF(I,J)/60.0_EB
       
       !Length to breadth ratio of ellipse based on effective umf
       LB = 0.936 * EXP(0.2566 * UMF_DUM) + 0.461 * EXP(-0.1548 * UMF_DUM) - 0.397 
       
       !Constraint LB max = 8 from Finney 2004
       LB = MAX(1.0_EB,MIN(LB,8.0_EB))
       
       LBD = SQRT(LB**2 - 1.0_EB)
       
       !Head to back ratio based on LB
       HB = (LB + LBD) / (LB - LBD)
       
       !*** A_ELPS and B_ELPS are *opposite* in notation from Farsite and Richards
       A_ELPS =  0.5_EB * (ROS_TMP + ROS_TMP/HB)
       A_ELPS2 = A_ELPS**2
       B_ELPS =  A_ELPS / LB !0.5_EB * (ROS_TMP + ROS_TMP/HB) / LB
       B_ELPS2=  B_ELPS**2
       C_ELPS =  A_ELPS - (ROS_TMP/HB)
       
  
       ! Denominator used in spread rate equation from Richards 1990 (also in Farsite)
       DENOM = B_ELPS2 * (YSF * COS_THETA + XSF * SIN_THETA)**2 + &
                A_ELPS2 * (XSF * COS_THETA - YSF * SIN_THETA)**2
           
       ! Finney's formulation
       !DENOM = B_ELPS2 * (XS * SIN_THETA - YS * COS_THETA)**2 - &
       !A_ELPS2 * (XS * COS_THETA + YS * SIN_THETA)**2
             
       IF (DENOM > 0_EB) THEN                 
        DENOM = 1._EB / SQRT(DENOM)        
       ELSE
        DENOM = 0._EB
       ENDIF
        
  
       IF (PHI_W(I,J) > 0._EB .OR. MAG_ZT(I,J) > 0._EB) THEN
   
        SR_X_LS(I,J) = DENOM * (B_ELPS2 * COS_THETA * (XSF * SIN_THETA + YSF * COS_THETA) -&
                        A_ELPS2 * SIN_THETA * (XSF * COS_THETA - YSF * SIN_THETA)) + C_ELPS * SIN_THETA
  
       
        SR_Y_LS(I,J) = DENOM * (-B_ELPS2 * SIN_THETA * (XSF * SIN_THETA + YSF * COS_THETA) -&
                        A_ELPS2 * COS_THETA * (XSF * COS_THETA - YSF * SIN_THETA)) + C_ELPS * COS_THETA
        
    
       ELSE
   
            !For no-wind, no-slope case
            SR_X_LS(I,J) = ROS_HEAD(I,J) * NORMAL_FIRELINE(1)
            SR_Y_LS(I,J) = ROS_HEAD(I,J) * NORMAL_FIRELINE(2)
        
       ENDIF  
       
       ! Spread rates are calculated for flat ground. Lines below are equivalent to projecting
       ! delta_x and delta_y from the slope to the plane surface, e.g., say dx/dz>0 and dy/dz=0. On the slope, 
       ! ROS_slp=dx_slp/dt and dx_surf=dx_slp*cos(theta) --> R_surf = dx_surf/dt = dx_slp/dt * cos(theta).
       
       IF (ABS(DZTDX(I,J)) > 0._EB) SR_X_LS(I,J) = SR_X_LS(I,J) * ABS(COS(ATAN(DZTDX(I,J))))
       IF (ABS(DZTDY(I,J)) > 0._EB) SR_Y_LS(I,J) = SR_Y_LS(I,J) * ABS(COS(ATAN(DZTDY(I,J))))
       
       MAG_SR = SQRT(SR_X_LS(I,J)**2 + SR_Y_LS(I,J)**2)   
   
   
   ELSE !McArthur Spread Model
             
        
           WIND_DOT_NORMAL_FIRELINE = U_LS(I,J)*NORMAL_FIRELINE(1) + V_LS(I,J)*NORMAL_FIRELINE(2)
           MAG_U  = SQRT(U_LS(I,J)**2 + V_LS(I,J)**2)

           COS_THETA_WIND = 0.0_EB ; COS_THETA_WIND_H = 0.0_EB ; COS_THETA_WIND_B = 0.0_EB
           IF(MAG_U > 0.0_EB) COS_THETA_WIND = WIND_DOT_NORMAL_FIRELINE/MAG_U

           GRAD_SLOPE_DOT_NORMAL_FIRELINE = DZTDX(I,J)*NORMAL_FIRELINE(1) + DZTDY(I,J)*NORMAL_FIRELINE(2) 
           COS_THETA_SLOPE = 0.0_EB ; COS_THETA_SLOPE_H = 0.0_EB ; COS_THETA_SLOPE_B = 0.0_EB
   
           IF (MAG_ZT(I,J) > 0.0_EB) THEN
               COS_THETA_SLOPE = GRAD_SLOPE_DOT_NORMAL_FIRELINE/MAG_ZT(I,J)
           ENDIF
   
           DEGREES_SLOPE = ATAN(MAG_ZT(I,J))*RAD_TO_DEGREE
    
           SLOPE_FACTOR  = MAG_ZT(I,J)**2
           IF (SLOPE_FACTOR > 3._EB) SLOPE_FACTOR = 3._EB
        
        
           ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
            IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
             IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
             IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)

              MAG_SR    = 0.0_EB
              ROS_HEADS = 0.0_EB
              ROS_BACKS = 0.0_EB

             NEXP_WIND = WIND_EXP(I,J)
        
  
        ! Spread with the wind and upslope
           IF(COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS = 0.33_EB*ROS_HEAD(I,J)
                IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =         ROS_HEAD(I,J)
                IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  3._EB*ROS_HEAD(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                    ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
           ENDIF
        !  IF(ABS(COS_THETA_WIND) < 0.5_EB .AND. MAG_F > 0._EB) MAG_SR = 0.0_EB
        !  IF(ABS(COS_THETA_WIND) < 0.5_EB .AND. MAG_F > 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
        !  IF(FLANKFIRE_LIFETIME(I,J) > TIME_FLANKFIRE_QUENCH) MAG_SR = 0.0_EB

        ! Spread with the wind and downslope
           IF(COS_THETA_WIND >= 0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_HEADS =  0.33_EB*ROS_HEAD(I,J)
            IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_HEADS =  0.50_EB*ROS_HEAD(I,J)
            IF(DEGREES_SLOPE >= 20._EB)                              ROS_HEADS =  0.75_EB*ROS_HEAD(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB + COS_THETA_WIND*COS_THETA_SLOPE) + &
                     (ROS_HEAD(I,J) - ROS_FLANK(I,J))*COS_THETA_WIND**NEXP_WIND + &
                     (ROS_HEADS     - ROS_FLANK(I,J))*COS_THETA_SLOPE  !magnitude of spread rate
        !   if(cos_theta_wind == 0._EB) FLANKFIRE_LIFETIME(I,J) = FLANKFIRE_LIFETIME(I,J) + DT_LS
        !   if(flankfire_lifetime(i,j) > time_flankfire_quench) mag_sr = 0.0_EB
           ENDIF

        ! Spread against the wind and upslope
           IF(COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE >= 0._EB) THEN
            IF (.NOT. LSET_TAN2) THEN
                IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = -0.33_EB*ROS_BACKU(I,J)
                IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS =         -ROS_BACKU(I,J)
                IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = -3.0_EB*ROS_BACKU(I,J)
            ELSEIF (DEGREES_SLOPE > 0._EB) THEN
                ROS_HEADS = ROS_HEAD(I,J) * SLOPE_FACTOR !Dependence on TAN(slope)^2
            ENDIF
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
           ENDIF

        ! Spread against the wind and downslope
           IF(COS_THETA_WIND <  0._EB .AND. COS_THETA_SLOPE < 0._EB) THEN
            IF(DEGREES_SLOPE >= 5._EB .AND. DEGREES_SLOPE < 10._EB)  ROS_BACKS = 0.33_EB*ROS_BACKU(I,J)
            IF(DEGREES_SLOPE >= 10._EB .AND. DEGREES_SLOPE < 20._EB) ROS_BACKS = 0.50_EB*ROS_BACKU(I,J)
            IF(DEGREES_SLOPE >= 20._EB)                              ROS_BACKS = 0.75_EB*ROS_BACKU(I,J)
            MAG_SR = ROS_FLANK(I,J)*(1._EB - ABS(COS_THETA_WIND)**NEXP_WIND*COS_THETA_SLOPE) + &
                     (ROS_FLANK(I,J) - ROS_BACKU(I,J))*(-ABS(COS_THETA_WIND)**NEXP_WIND) + &
                     (ROS_FLANK(I,J) - ROS_BACKS)*COS_THETA_SLOPE  !magnitude of spread rate
           ENDIF



        !  MAG_SR = ROS_FLANK(I,J) + ROS_HEAD(I,J)*COS_THETA_WIND**1.5 !magnitude of spread rate
        !  MAG_SR = ROS_FLANK(I,J) + ROS_HEAD(I,J)*MAG_U*COS_THETA_WIND**1.5 !magnitude of spread rate
           SR_X_LS(I,J) = MAG_SR*NORMAL_FIRELINE(1) !spread rate components
           SR_Y_LS(I,J) = MAG_SR*NORMAL_FIRELINE(2) 
        !  MAG_SR_OUT(I,J) = MAG_SR
  
   ENDIF !Ellipse or McArthur Spread 
   
   DYN_SR_MAX = MAX(DYN_SR_MAX,MAG_SR) 

  ENDDO

ENDDO FLUX_ILOOP




END SUBROUTINE LEVEL_SET_SPREAD_RATE 

!--------------------------------------------------------------------
!
SUBROUTINE LEVEL_SET_ADVECT_FLUX
!
! Use the spread rate [SR_X_LS,SR_Y_LS] to compute the limited scalar gradient
! and take dot product with spread rate vector to get advective flux

INTEGER :: I,IM1,IM2,IP1,IP2,J,JM1,JM2,JP1,JP2
REAL(EB), DIMENSION(:) :: Z(4)
REAL(EB), DIMENSION(:,:) :: FLUX_LS(NX_LS,NY_LS)
REAL(EB) :: DPHIDX,DPHIDY,F_EAST,F_WEST,F_NORTH,F_SOUTH
REAL(EB) :: PHIMAG

IF (RK2_PREDICTOR_LS) PHI0_LS = PHI_LS
IF (.NOT. RK2_PREDICTOR_LS) PHI0_LS = PHI1_LS

ILOOP: DO I = 1,NX_LS

 IM1=I-1; IF (IM1<1) IM1=IM1+NX_LS
 IM2=I-2; IF (IM2<1) IM2=IM2+NX_LS

 IP1=I+1; IF (IP1>NX_LS) IP1=IP1-NX_LS
 IP2=I+2; IF (IP2>NX_LS) IP2=IP2-NX_LS

 DO J = 1,NY_LS
   
   JM1=J-1
   JM2=J-2
   JP1=J+1
   JP2=J+2

!-- east face
   Z(1) = PHI0_LS(IM1,J)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(IP1,J)
   Z(4) = PHI0_LS(IP2,J)
   F_EAST = SCALAR_FACE_VALUE(SR_X_LS(I,J),Z,LIMITER_LS)
   
!-- west face
   Z(1) = PHI0_LS(IM2,J)
   Z(2) = PHI0_LS(IM1,J)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(IP1,J)
   F_WEST = SCALAR_FACE_VALUe(SR_X_LS(I,J),Z,LIMITER_LS)
   
 IF (J<2 .OR. J>(NY_LS-2)) THEN  
   
      IF (J==1) THEN
    !    north face
        Z(1) = PHI_MAX_LS
        Z(2) = PHI0_LS(I,J)
        Z(3) = PHI0_LS(I,JP1)
        Z(4) = PHI0_LS(I,JP2)
        F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

    !    south face
        Z(1) = PHI_MAX_LS
        Z(2) = PHI_MAX_LS
        Z(3) = PHI0_LS(I,J)
        Z(4) = PHI0_LS(I,JP1)
        F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

       ELSEIF (j==2) THEN
    !    north face
        Z(1) = PHI0_LS(I,JM1)
        Z(2) = PHI0_LS(I,J)
        Z(3) = PHI0_LS(I,JP1)
        Z(4) = PHI0_LS(I,JP2)
        F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

    !    south face
        Z(1) = PHI_MAX_LS
        Z(2) = PHI0_LS(I,JM1)
        Z(3) = PHI0_LS(I,J)
        Z(4) = PHI0_LS(I,JP1)
        F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)
   
   
        ELSEIF (J == NY_LS-1) THEN
    !    north face
        Z(1) = PHI0_LS(I,JM1)
        Z(2) = PHI0_LS(I,J)
        Z(3) = PHI0_LS(I,JP1)
        Z(4) = PHI_MIN_LS
        F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

    !    south face
        Z(1) = PHI0_LS(I,JM2)
        Z(2) = PHI0_LS(I,JM1)
        Z(3) = PHI0_LS(I,J)
        Z(4) = PHI0_LS(I,JP1)
        F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

       ELSE ! must be J == NY_LS
    !    north face
        Z(1) = PHI0_LS(I,JM1)
        Z(2) = PHI0_LS(I,J)
        Z(3) = PHI_MIN_LS
        Z(4) = PHI_MIN_LS
        F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

    !    south face
        Z(1) = PHI0_LS(I,JM2)
        Z(2) = PHI0_LS(I,JM1)
        Z(3) = PHI0_LS(I,J)
        Z(4) = PHI_MIN_LS
        F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)
    
       ENDIF !IF (J==1) 
   
   ELSE
!    north face
   Z(1) = PHI0_LS(I,JM1)
   Z(2) = PHI0_LS(I,J)
   Z(3) = PHI0_LS(I,JP1)
   Z(4) = PHI0_LS(I,JP2)
   F_NORTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)

!    south face
   Z(1) = PHI0_LS(I,JM2)
   Z(2) = PHI0_LS(I,JM1)
   Z(3) = PHI0_LS(I,J)
   Z(4) = PHI0_LS(I,JP1)
   F_SOUTH = SCALAR_FACE_VALUE(SR_Y_LS(I,J),Z,LIMITER_LS)
   
   ENDIF !IF (J<2 .OR. J>(NY_LS-2) 
        
   DPHIDX = (F_EAST-F_WEST)/DX_LS
   DPHIDY = (F_NORTH-F_SOUTH)/DY_LS
   FLUX_LS(I,J) = SR_X_LS(I,J)*DPHIDX + SR_Y_LS(I,J)*DPHIDY
   
   PHIMAG          = SQRT(DPHIDX**2 + DPHIDY**2)
   MAG_SR_OUT(I,J) = 0.0_EB
   IF(PHIMAG > 0.0_EB) MAG_SR_OUT(I,J) = FLUX_LS(I,J)/PHIMAG
        
!  fx = (f_east-f_west)/dx
!  fy = (f_north-f_south)/dy
!       phi(i,j) = phi0(i,j) - dt*[Fx(i,j) Fy(i,j)]*[fx fy]
 ENDDO


 FLUX_LS(:,1) = FLUX_LS(:,2)



ENDDO ILOOP

IF (RK2_PREDICTOR_LS) FLUX0_LS = FLUX_LS
IF (.NOT. RK2_PREDICTOR_LS) FLUX1_LS = FLUX_LS


END SUBROUTINE LEVEL_SET_ADVECT_FLUX 
!
! ----------------------------------------------------
REAL(EB) FUNCTION SCALAR_FACE_VALUE(SR_XY,Z,LIMITER)
!
! From Randy 7-11-08
! This function computes the scalar value on a face.
! The scalar is denoted Z, and the velocity is denoted U.
! The gradient (computed elsewhere) is a central difference across 
! the face subject to a flux limiter.  The flux limiter choices are:
! 
! limiter = 1 implements the MINMOD limiter
! limiter = 2 implements the SUPERBEE limiter of Roe
! limiter = 3 implements first-order upwinding (monotone)
!
!
!                    location of face
!                            
!                            f
!    |     o     |     o     |     o     |     o     |
!                     SRXY        SRXY
!                 (if f_east)  (if f_west)
!         Z(1)        Z(2)        Z(3)        Z(4)
!
INTEGER :: LIMITER
REAL(EB) :: SR_XY
REAL(EB), INTENT(IN), DIMENSION(4) :: Z
REAL(EB) :: B,DZLOC,DZUP,R,ZUP,ZDWN

IF (SR_XY > 0._EB) THEN
!     the flow is left to right
 DZLOC = Z(3)-Z(2)
 DZUP  = Z(2)-Z(1)

 IF (ABS(DZLOC) > 0._EB) THEN
  R = DZUP/DZLOC
 ELSE
  R = 0._EB
 ENDIF
 ZUP  = Z(2)
 ZDWN = Z(3)
ELSE
!     the flow is right to left
 DZLOC = Z(3)-Z(2)
 DZUP  = Z(4)-Z(3)

 IF (ABS(DZLOC) > 0._EB) THEN
  R = DZUP/DZLOC
 ELSE
  r = 0._EB
 ENDIF
  ZUP  = Z(3)
  ZDWN = Z(2)
ENDIF

! flux limiter
IF (LIMITER==1) THEN
!     MINMOD
    B = MAX(0._EB,MIN(1._EB,R))
ELSEIF (limiter==2) THEN
!     SUPERBEE
    B = MAX(0._EB,MIN(2._EB*R,1._EB),MIN(R,2._EB))
ELSEIF (limiter==3) THEN
!     first-order upwinding
    B = 0._EB
ENDIF

SCALAR_FACE_VALUE = ZUP + 0.5*B*( ZDWN - ZUP )

END FUNCTION SCALAR_FACE_VALUE

SUBROUTINE GET_REV_vege(MODULE_REV,MODULE_DATE)
INTEGER,INTENT(INOUT) :: MODULE_REV
CHARACTER(255),INTENT(INOUT) :: MODULE_DATE

WRITE(MODULE_DATE,'(A)') vegerev(INDEX(vegerev,':')+1:LEN_TRIM(vegerev)-2)
READ (MODULE_DATE,'(I5)') MODULE_REV
WRITE(MODULE_DATE,'(A)') vegedate

END SUBROUTINE GET_REV_vege


END MODULE VEGE
