      FUNCTION X_sps( Model, nodenumber, dumy) RESULT(betav)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,betav
            INTEGER :: nodenumber, nentries, i, ind
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                WRITE(beta_fn,"(A,I0.3,A)") '../edc_data/x_sps.txt'

                OPEN(10, file=beta_fn)
                READ(10, *) header, nentries

                ALLOCATE(dist(nentries), betaarr(nentries))

                DO i=1,nentries
                    READ(10, *) dist(i), betaarr(i)
                END DO
                CLOSE(10)

                dx = (dist(nentries) - dist(1)) / (nentries - 1) 
            END IF

            x=Model % Nodes % x (nodenumber)
            ind = floor((x-dist(1)) / dx) + 1
            ind = max(ind, 1)
            ind = min(ind, nentries - 1)
            betav = betaarr(ind) + (x - dist(ind)) * &
                (betaarr(ind+1) - betaarr(ind)) / (dist(ind+1) - dist(ind))

            RETURN 
        END

      FUNCTION Y_sps( Model, nodenumber, dumy) RESULT(betav)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,betav
            INTEGER :: nodenumber, fln, nentries, i, ind
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                WRITE(beta_fn,"(A,I0.3,A)") '../edc_data/y_sps.txt'

                OPEN(10, file=beta_fn)
                READ(10, *) header, nentries

                ALLOCATE(dist(nentries), betaarr(nentries))

                DO i=1,nentries
                    READ(10, *) dist(i), betaarr(i)
                END DO
                CLOSE(10)

                dx = (dist(nentries) - dist(1)) / (nentries - 1) 
            END IF

            x=Model % Nodes % x (nodenumber)
            ind = floor((x-dist(1)) / dx) + 1
            ind = max(ind, 1)
            ind = min(ind, nentries - 1)
            betav = betaarr(ind) + (x - dist(ind)) * &
                (betaarr(ind+1) - betaarr(ind)) / (dist(ind+1) - dist(ind))

            RETURN 
        END


        ! Helper function to deal with interpolating unvenly spaced data
        ! Used for the unven spacing of the ice-core records
        FUNCTION UnevenInterp(x0, x, y, nx) RESULT(y0)
            USE Types
            INTEGER :: nx, i
            REAL(KIND=dp) :: x(nx), y(nx)
            REAL(KIND=dp) :: x0, y0
            i = 1
            IF (x0.LE.x(1)) THEN
                y0 = y(1)
            ELSE IF (x0.GE.x(nx)) THEN
                y0 = y(nx)
            ELSE
                DO WHILE (x(i).LT.x0)
                   i = i + 1
                END DO
                y0 = (y(i-1) * (x(i) - x0) + y(i) * (x0 - x(i - 1))) / (x(i) - x(i - 1))
            END IF
            RETURN
         END


      ! Rescale the surface temperature by a history
      FUNCTION TDSurfTemp( Model, nodenumber, dumy) RESULT(f)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(Variable_t), POINTER :: TimeVar
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy,f
            INTEGER :: nodenumber, nentries, i, ind, nentries1
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n, time_offset
            REAL(kind=dp) :: time, alwaystime
            REAL(kind=dp) :: UnevenInterp, divtempnow
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr, times, dividetemp
            LOGICAL :: gotit, FIRSTTIME=.TRUE., fixedtscale=.FALSE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME, nentries1, dividetemp, times, fixedtscale, alwaystime, time_offset

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                Material => GetMaterial()
                alwaystime=GetInteger(Material, "FixedTime", gotit)
                IF (.not.gotit) THEN
                    fixedtscale=.FALSE.
                    alwaystime=0.
                ELSE
                    fixedtscale=.TRUE.
                END IF
                time_offset=GetInteger(Material, "TimeOffset", gotit)
                IF (.not.gotit) THEN
                    time_offset=0.
                ELSE
                    WRITE(Message,'(A, f8.0)') 'Time offset is', time_offset
                    CALL INFO('TDSurfTemp',Message,Level=3)
                END IF

                WRITE(beta_fn,"(A,I0.3, A)")'../model_inputs/surftemp_c.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentries
                ALLOCATE(dist(nentries), betaarr(nentries))
                DO i=1,nentries
                    READ(10, *) dist(i), betaarr(i)
                END DO
                CLOSE(10)
                dx = (dist(nentries) - dist(1)) / (nentries - 1) 

                beta_fn = '../model_inputs/GRIP_temperature.tab'
                OPEN(10, file=beta_fn)
                READ(10, *) nentries1
                ALLOCATE(times(nentries1), dividetemp(nentries1))
                DO i=1,nentries1
                    READ(10, *) times(i), dividetemp(i)
                END DO
                CLOSE(10)
            END IF


            x=Model % Nodes % x (nodenumber)
            ind = floor((x-dist(1)) / dx) + 1
            ind = max(ind, 1)
            ind = min(ind, nentries - 1)
            f = betaarr(ind) + (x - dist(ind)) * &
                (betaarr(ind+1) - betaarr(ind)) / (dist(ind+1) - dist(ind))
            f = f + 273.15

            IF (fixedtscale) THEN
                time = alwaystime
            ELSE
                TimeVar => VariableGet(Model % Mesh % Variables, 'Time')
                Time = TimeVar % Values(1)
            END IF
            time = time + time_offset
            divtempnow = UnevenInterp(time, times, dividetemp, nentries1)
            f = f + divtempnow - betaarr(1)
            RETURN 
        END


      ! Rescale the surface temperature by a history
      FUNCTION SpatConstTemp( Model, nodenumber, dumy) RESULT(f)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(Variable_t), POINTER :: TimeVar
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy,f
            INTEGER :: nodenumber, fln, nentries, i, ind, nentries1
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n, time_offset
            REAL(kind=dp) :: time, alwaystime
            REAL(kind=dp) :: UnevenInterp, divtempnow, dumb, dumber, dumbest
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr, times, dividetemp
            LOGICAL :: gotit, FIRSTTIME=.TRUE., fixedtscale=.FALSE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME, nentries1, dividetemp, times, fixedtscale, alwaystime, time_offset

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                Material => GetMaterial()
                alwaystime=GetInteger(Material, "FixedTime", gotit)
                IF (.not.gotit) THEN
                    fixedtscale=.FALSE.
                    alwaystime=0.
                ELSE
                    fixedtscale=.TRUE.
                END IF

                time_offset=GetInteger(Material, "TimeOffset", gotit)
                IF (.not.gotit) THEN
                    time_offset=0.
                ELSE
                    WRITE(Message,'(A, f8.0)') 'Time offset is', time_offset
                    CALL INFO('TDSurfTemp',Message,Level=3)
                END IF

                beta_fn = '../edc_data/jouzel_2007_EDC_dD_temp.tab'
                OPEN(10, file=beta_fn)
                READ(10, *) nentries1
                ALLOCATE(times(nentries1), dividetemp(nentries1))
                DO i=1,nentries1
                    READ(10, *) dumb, times(i), dumber, dividetemp(i), dumbest
                    times(i) = times(i) * 1000.0
                END DO
                CLOSE(10)
            END IF


            
            ! Take the average temperature above 100 m in the borehole as mean
            f = -54.0 + 273.15

            IF (fixedtscale) THEN
                time = alwaystime
            ELSE
                TimeVar => VariableGet(Model % Mesh % Variables, 'Time')
                Time = TimeVar % Values(1)
            END IF
            time = time_offset - time
            divtempnow = UnevenInterp(time, times, dividetemp, nentries1)
            f = f + divtempnow
            RETURN 
        end


      ! Rescale the surface temperature by a history
      FUNCTION SpatConstAcc( Model, nodenumber, dumy) RESULT(f)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(Variable_t), POINTER :: TimeVar
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy,f
            INTEGER :: nodenumber, fln, nentries, i, ind, nentries1
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n, time_offset
            REAL(kind=dp) :: time, alwaystime
            REAL(kind=dp) :: UnevenInterp, accnow, dumb, dumber, dumbest
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr, times, accarr
            LOGICAL :: gotit, FIRSTTIME=.TRUE., fixedtscale=.FALSE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME, nentries1, accarr, times, fixedtscale, alwaystime, time_offset

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                Material => GetMaterial()
                alwaystime=GetInteger(Material, "FixedTime", gotit)
                IF (.not.gotit) THEN
                    fixedtscale=.FALSE.
                    alwaystime=0.
                ELSE
                    fixedtscale=.TRUE.
                END IF

                time_offset=GetInteger(Material, "TimeOffset", gotit)
                IF (.not.gotit) THEN
                    time_offset=0.
                ELSE
                    WRITE(Message,'(A, f8.0)') 'Time offset is', time_offset
                    CALL INFO('TDSurfTemp',Message,Level=3)
                END IF

                beta_fn = '../edc_data/AICC2012_acc.csv'
                OPEN(10, file=beta_fn)
                READ(10, *) nentries1
                ALLOCATE(times(nentries1), accarr(nentries1))
                DO i=1,nentries1
                    READ(10, *) times(i), accarr(i)
                END DO
                CLOSE(10)
            END IF

            IF (fixedtscale) THEN
                time = alwaystime
            ELSE
                TimeVar => VariableGet(Model % Mesh % Variables, 'Time')
                Time = TimeVar % Values(1)
            END IF
            time = time_offset - time
            f = UnevenInterp(time, times, accarr, nentries1)
            RETURN 
        end

      FUNCTION SpatialAcc( Model, nodenumber, dumy) RESULT(betav)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,betav
            INTEGER :: nodenumber, fln, nentries, i, ind
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                WRITE(beta_fn,"(A,I0.3,A)") '../edc_data/acc.txt'

                OPEN(10, file=beta_fn)
                READ(10, *) header, nentries

                ALLOCATE(dist(nentries), betaarr(nentries))

                DO i=1,nentries
                    READ(10, *) dist(i), betaarr(i)
                END DO
                CLOSE(10)

                dx = (dist(nentries) - dist(1)) / (nentries - 1) 
            END IF

            x=Model % Nodes % x (nodenumber)
            ind = floor((x-dist(1)) / dx) + 1
            ind = max(ind, 1)
            ind = min(ind, nentries - 1)
            betav = betaarr(ind) + (x - dist(ind)) * &
                (betaarr(ind+1) - betaarr(ind)) / (dist(ind+1) - dist(ind))

            RETURN 
        END

      FUNCTION DualVarAcc( Model, NodeNumber, dumy) RESULT(ACC)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: NodeNumber
            REAL(kind=dp) :: dumy, Acc, SpatConstAcc, SpatialAcc
            ACC = SpatConstAcc( Model, NodeNumber, dumy ) * &
                    SpatialAcc( Model, NodeNumber, dumy )
      END

      FUNCTION SpatialTemp( Model, nodenumber, dumy) RESULT(betav)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            REAL(kind=dp) :: dumy,betav
            INTEGER :: nodenumber, fln, nentries, i, ind
            REAL(kind=dp) :: x, dx
            REAL(kind=dp) :: n
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: dist, betaarr
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            REAL(kind=dp) :: A, rho, g, H, alpha, yis

            SAVE dx, dist, nentries, betaarr, FIRSTTIME

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                WRITE(beta_fn,"(A,I0.3,A)") '../edc_data/ts.txt'

                OPEN(10, file=beta_fn)
                READ(10, *) header, nentries

                ALLOCATE(dist(nentries), betaarr(nentries))

                DO i=1,nentries
                    READ(10, *) dist(i), betaarr(i)
                END DO
                CLOSE(10)

                dx = (dist(nentries) - dist(1)) / (nentries - 1) 
            END IF

            x=Model % Nodes % x (nodenumber)
            ind = floor((x-dist(1)) / dx) + 1
            ind = max(ind, 1)
            ind = min(ind, nentries - 1)
            betav = betaarr(ind) + (x - dist(ind)) * &
                (betaarr(ind+1) - betaarr(ind)) / (dist(ind+1) - dist(ind))

            RETURN 
        END

      FUNCTION DualVarTemp( Model, NodeNumber, dumy) RESULT(Temp)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: NodeNumber
            REAL(kind=dp) :: dumy, Temp, SpatConstTemp, SpatialTemp
            temp = SpatConstTemp( Model, NodeNumber, dumy ) + &
                     SpatialTemp( Model, NodeNumber, dumy )
      END

      FUNCTION BuizertTemp( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth
            CHARACTER(len=1) :: header

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)

            SAVE depth, nentries, betaarr, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/buizert_2021_boreholeT.csv'
                OPEN(10, file=beta_fn)
                READ(10, *) nentries
                ALLOCATE(depth(nentries), betaarr(nentries))
                DO i=1,nentries
                    READ(10, *) depth(i), betaarr(i)
                    depth(i) = -depth(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries) + 273.15
            RETURN 
        END


       FUNCTION BuizertTempC( Model, nodenumber, dumy) RESULT(T)
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            INTEGER :: nodenumber
            REAL(kind=dp) :: dumy, T, BuizertTemp
            T = BuizertTemp( Model, nodenumber, dumy) - 273.15
            RETURN
       END


      FUNCTION AiccAge( Model, nodenumber, dumy) RESULT(age)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, age
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, nentries, betaarr, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/AICC2012_ages.csv'
                OPEN(10, file=beta_fn)
                READ(10, *) nentries
                ALLOCATE(depth(nentries), betaarr(nentries))
                DO i=1,nentries
                    READ(10, *) depth(i), betaarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            age = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION ZsTop( Model, nodenumber, dumy) RESULT(zs)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, age
            INTEGER :: nodenumber, nentriessurf, i, ind
            REAL(kind=dp) :: x, zs
            REAL(kind=dp) :: n

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn
            CHARACTER(len=1) :: header

            SAVE FIRSTTIME, xsurf, surfarr, nentriessurf

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)
            END IF

            x=Model % Nodes % x (nodenumber)

            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            RETURN 
        END


      FUNCTION nlm1( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth
            CHARACTER(len=1) :: header

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)

                nentries = 3
                depth = (/ 0.0_dp, 2500.0_dp, 3300.0_dp /)
                betaarr = (/ 0.0_dp, -0.27039194_dp, -0.27039194_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)


            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            IF (realdepth.lt.250.0_dp)  T = 0.0_dp
            RETURN 
        END

      FUNCTION nlm3( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)

                nentries = 3
                depth = (/ 0.0_dp, 2500.0_dp, 3300.0_dp /)
                betaarr = (/ 0.0_dp, -0.2207741_dp, -0.2207741_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END


      FUNCTION nlm5( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)

                nentries = 3
                depth = (/ 0.0_dp, 2000.0_dp, 3300.0_dp /)
                betaarr = (/ 0.0_dp, -0.27039194_dp, -0.27039194_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm2weak( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth
            CHARACTER(len=1) :: header

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)

                nentries = 3
                depth = (/ 0.0_dp, 2000.0_dp, 3300.0_dp /)
                betaarr = (/ 0.0_dp, -0.027_dp, -0.027_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm4weak( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)

                nentries = 3
                depth = (/ 0.0_dp, 2000.0_dp, 3300.0_dp /)
                betaarr = (/ 0.0_dp, -0.022_dp, -0.022_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END


      FUNCTION nlm6weak( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                beta_fn = '../edc_data/surf.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriessurf
                ALLOCATE(xsurf(nentriessurf), surfarr(nentriessurf))
                DO i=1,nentriessurf
                    READ(10, *) xsurf(i), surfarr(i)
                END DO
                CLOSE(10)

                beta_fn = '../edc_data/bed.txt'
                OPEN(10, file=beta_fn)
                READ(10, *) header, nentriesbed
                ALLOCATE(xbed(nentriesbed), bedarr(nentriesbed))
                DO i=1,nentriesbed
                    READ(10, *) xbed(i), bedarr(i)
                END DO
                CLOSE(10)

                nentries = 3
                depth = (/ 0.0_dp, 2000.0_dp, 3300.0_dp /)
                betaarr = (/ 0.0_dp, -0.027_dp, -0.027_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = UnevenInterp(x, xbed, bedarr, nentriesbed)
            zs = UnevenInterp(x, xsurf, surfarr, nentriessurf)

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm2sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, -0.0_dp, -0.269882_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm4sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, 0.0_dp, -0.220358_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm6sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                nentries = 3
                depth = (/ 100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, -0.0_dp, -0.269882_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm7sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.
                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, 0.0_dp, 0.186445_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm9sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, 0.0_dp, 0.140939_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END

      FUNCTION nlm11sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, 0.0_dp, 0.133706_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END


      FUNCTION nlm13sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, 0.0_dp, 0.140939_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END


      FUNCTION nlm15sf( Model, nodenumber, dumy) RESULT(T)
            ! Temp with depth
            USE DefUtils
            USE types

            IMPLICIT NONE
            TYPE(Model_t) :: Model
            TYPE(ValueList_t), POINTER :: Material
            REAL(kind=dp) :: dumy, T
            INTEGER :: nodenumber, fln, nentries, nentriessurf, nentriesBed, i, ind
            REAL(kind=dp) :: x, y, dy, zb, zs, PercFromTop, dxsurf, dxbed
            REAL(kind=dp) :: n, width, realdepth

            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: depth, betaarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xsurf, surfarr
            REAL(kind=dp), DIMENSION(:), ALLOCATABLE :: xbed, bedarr
            REAL(kind=dp) :: UnevenInterp
            LOGICAL :: gotit, FIRSTTIME=.TRUE.
            CHARACTER(len=MAX_NAME_LEN) :: beta_fn

            REAL(kind=dp) :: A, rho, g, H, alpha, yis
            TYPE(Variable_t), POINTER :: ZsTopSol, ZsBotSol
            REAL(KIND=dp), POINTER :: ZsTop(:), ZsBot(:)
            INTEGER, POINTER :: ZsTopPerm(:), ZsBotPerm(:)
            CHARACTER(len=1) :: header

            SAVE depth, betaarr, nentries, FIRSTTIME, xbed, xsurf,&
                 bedarr, surfarr, nentriessurf, nentriesbed

            IF (FIRSTTIME) THEN
                FIRSTTIME = .FALSE.

                nentries = 3
                depth = (/ -100.0_dp, 0.0_dp, 2000.0_dp /)
                betaarr = (/ 0.0_dp, 0.0_dp, 0.186445_dp /)
            END IF

            x=Model % Nodes % x (nodenumber)
            y=Model % Nodes % y (nodenumber)

            zb = 0.0_dp
            zs = 2000.0_dp

            realdepth = (zs - y) / ( zs - zb) * depth(nentries)

            T = UnevenInterp(realdepth, depth, betaarr, nentries)
            RETURN 
        END
