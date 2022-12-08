      FUNCTION USIA( Model, nodenumber, dumy) RESULT(xvel)
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,xvel
            INTEGER :: nodenumber
            REAL(kind=dp) :: x,y
            REAL(kind=dp) :: n, u_s

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            A = 4.9e-25  ! in s^-1Pa^-3
            rho = 910
            g = 9.8
            alpha = 1.0e-3
            H = 2000.
            n = 3.0

            yis = 365.24 * 24. * 60. * 60.

            u_s = 0.6666666666666666666666

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            xvel = sign(u_s * (1.0_dp - ((H - y) / H) ** (n + 1)), x)
            RETURN 
        END

        ! Deal with units and combining basal heat terms
      FUNCTION FricAndGeo(  Model, Node, BetaAndF)RESULT(totalheat)
      ! A little helper in dealing with the different heat terms at the bed
      ! we want both frictional heating and geothermal
      ! it is easiest to combine these in fortran rather than dealing
      ! with separate terms in Elmer/Ice
      USE Types
      USE DefUtils

      IMPLICIT NONE
      
      !------------------------------------------------------------------------------
      TYPE(Model_t) :: Model
      INTEGER :: Node
      REAL(KIND=dp) :: BetaAndF(2), frictionHeat
      !----------------------------------------------------------------------------
      
      INTEGER :: DIM, i, j
      REAL(KIND=dp), POINTER :: FlowValues(:),NormalValues(:)
      REAL(KIND=dp) :: normal(3), velo(3), un, ut
      REAL(KIND=dp) :: geo, totalheat
      INTEGER, POINTER :: FlowPerm(:), NormalPerm(:)
      LOGICAL :: FirstTime=.TRUE.,UnFoundFatal=.TRUE.
      TYPE(Variable_t), POINTER :: FlowSol,NormalVar
      CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
      
      SAVE DIM, FunctionName
      
      IF (FirstTime) THEN
         WRITE(FunctionName, '(A)') 'getFrictionHeat'
         DIM = CoordinateSystemDimension()
         FirstTime = .FALSE.
      END IF
      
      ! Get the variable velocity
      !---------------------------
      FlowSol => VariableGet( Model % Variables, 'AIFlow',UnFoundFatal=UnFoundFatal)
      FlowPerm    => FlowSol % Perm
      FlowValues  => FlowSol % Values
      
      ! Get the variable for normal vector
      NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
      NormalPerm => NormalVar % Perm
      NormalValues => NormalVar % Values
      
      DO i=1, DIM
         normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
         velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
      END DO
      
      !Tangential velocity
      un = SUM(velo(1:DIM)*(normal(1:DIM))) 
      ut = SQRT(SUM( (velo(1:DIM))**2.0_dp ) - (un**2.0_dp) )

      ! First mW/m^2, then elmer's conversion
      geo = abs(BetaAndF(2)*1.0e-3*1.0e-6*365.25*24.*60.*60.)

      ! the force is (ut*beta)^2
      frictionHeat = min(1.0, abs(ut*ut*BetaAndF(1)*BetaAndF(1)))
      totalheat = geo+frictionHeat
    END FUNCTION FricAndGeo

    FUNCTION Fric(  Model, Node, BetaAndF)RESULT(totalheat)
      ! A little helper in dealing with the different heat terms at the bed
      ! we want both frictional heating and geothermal
      ! it is easiest to combine these in fortran rather than dealing
      ! with separate terms in Elmer/Ice
      USE Types
      USE DefUtils

      IMPLICIT NONE
      
      !------------------------------------------------------------------------------
      TYPE(Model_t) :: Model
      INTEGER :: Node
      REAL(KIND=dp) :: BetaAndF, frictionHeat
      !----------------------------------------------------------------------------
      
      INTEGER :: DIM, i, j
      REAL(KIND=dp), POINTER :: FlowValues(:),NormalValues(:)
      REAL(KIND=dp) :: normal(3), velo(3), un, ut
      REAL(KIND=dp) :: geo, totalheat
      INTEGER, POINTER :: FlowPerm(:), NormalPerm(:)
      LOGICAL :: FirstTime=.TRUE.,UnFoundFatal=.TRUE.
      TYPE(Variable_t), POINTER :: FlowSol,NormalVar
      CHARACTER(LEN=MAX_NAME_LEN) :: FunctionName
      
      SAVE DIM, FunctionName
      
      IF (FirstTime) THEN
         WRITE(FunctionName, '(A)') 'getFrictionHeat'
         DIM = CoordinateSystemDimension()
         FirstTime = .FALSE.
      END IF
      
      ! Get the variable velocity
      !---------------------------
      FlowSol => VariableGet( Model % Variables, 'AIFlow',UnFoundFatal=UnFoundFatal)
      FlowPerm    => FlowSol % Perm
      FlowValues  => FlowSol % Values
      
      ! Get the variable for normal vector
      NormalVar =>  VariableGet(Model % Variables,'Normal Vector',UnFoundFatal=UnFoundFatal)
      NormalPerm => NormalVar % Perm
      NormalValues => NormalVar % Values
      
      DO i=1, DIM
         normal(i) = NormalValues(DIM*(NormalPerm(Node)-1) + i)      
         velo(i) = FlowValues( (DIM+1)*(FlowPerm(Node)-1) + i )
      END DO
      
      !Tangential velocity
      un = SUM(velo(1:DIM)*(normal(1:DIM))) 
      ut = SQRT(SUM( (velo(1:DIM))**2.0_dp ) - (un**2.0_dp) )

      ! the force is (ut*beta)^2
      frictionHeat = ut*ut*BetaAndF*BetaAndF

      ! We will return this in mW
      totalheat = frictionHeat / (1.0e-3*1.0e-6*365.25*24.*60.*60.)
    END FUNCTION Fric

      FUNCTION AgeFromCFRaymond(Model, nodenumber, inputs) RESULT(age)
          ! Calculate an age by integrating the veloicty profile
          ! Defined by VCFRaymond
          ! I feel like it is unclear how to pick between this and the divide
          ! profile, but I think this is more self-consistent
          USE Types
          IMPLICIT NONE
          TYPE(Model_t) :: Model
          REAL(kind=dp) :: yvel, age
          REAL(kind=dp) :: acc, surf, bed
          REAL(kind=dp) :: inputs
          INTEGER :: nodenumber, i
          REAL(kind=dp) :: y, frac
          REAL(kind=dp) :: n=3.0_dp, dy=1.0_dp, agelim=2.0e5

          acc = 0.3
          surf = 2000.0_dp
          bed = 0.0_dp
          y=Model % Nodes % y (nodenumber)
          age = 0

          ! simple do loop, age good to nearest meter should be fine
          DO i=1,NINT((surf - y) / dy)
                frac = FLOAT(i) / (surf - bed)
                yvel = abs(acc * (1.0_dp - frac * ((n + 2.0_dp) / (n + 1.0_dp) - 1.0_dp / (1.0_dp + n) * frac ** (n + 1))))
                age = age + dy / yvel
          END DO
          age = MIN(age, agelim)

          RETURN 
      END

