MODULE nswclib

  USE nswcconstant
  use precision
  use globals
  
  IMPLICIT NONE

  !
  !  These routines are part of the Naval Surface Warfare Center (NSWC) Mathematics Library
  !  as provided by This Fortran 90 version has been compiled by Alan Miller
  !
  !      http://www.ozemail.com.au/~milleraj/
  !  
  !  I extracted the bessel function routines which is all that I needed for
  !  the code to work.
  !
  


CONTAINS



  SUBROUTINE bssli(mo, a, in, w)

    !     ******************************************************************
    !     FORTRAN SUBROUTINE FOR MODIFIED BESSEL FUNCTION OF INTEGRAL ORDER
    !     ******************************************************************
    !     MO = MODE OF OPERATION
    !     A  = ARGUMENT (COMPLEX NUMBER)
    !     IN = ORDER (INTEGER)
    !     W  = FUNCTION OF FIRST KIND (COMPLEX NUMBER)
    !     -------------------

    COMPLEX(cdp), INTENT(IN)  :: a
    INTEGER(idp), INTENT(IN)  :: mo, in

    COMPLEX(cdp), INTENT(OUT) :: w
    REAL(fdp), PARAMETER :: cd(30) = (/ 0.00000000000000,  &
         -1.64899505142212E-2, -7.18621880068536E-2, -1.67086878124866E-1,  &
         -3.02582250219469E-1, -4.80613945245927E-1, -7.07075239357898E-1,  &
         -9.92995790539516E-1, -1.35583925612592, -1.82105907899132,  &
         -2.42482175310879, -3.21956655708750, -4.28658077248384,  &
         -5.77022816798128, -8.01371260952526, 0.00000000000000,  &
         -5.57742429879505E-3, -4.99112944172476E-2, -1.37440911652397E-1,  &
         -2.67233784710566E-1, -4.40380166808682E-1, -6.61813614872541E-1,  &
         -9.41861077665017E-1, -1.29754130468326, -1.75407696719816,  &
         -2.34755299882276, -3.13041332689196, -4.18397120563729,  &
         -5.65251799214994, -7.87863959810677 /), ce(30)  &
         = (/ 0.00000000000000, -4.80942336387447E-3,  &
         -1.31366200347759E-2, -1.94843834008458E-2, -2.19948900032003E-2,  &
         -2.09396625676519E-2, -1.74600268458650E-2, -1.27937813362085E-2,  &
         -8.05234421796592E-3, -4.15817375002760E-3, -1.64317738747922E-3,  &
         -4.49175585314709E-4, -7.28594765574007E-5, -5.38265230658285E-6,  &
         -9.93779048036289E-8, 0.00000000000000, 7.53805779200591E-2,  &
         7.12293537403464E-2, 6.33116224228200E-2, 5.28240264523301E-2,  &
         4.13305359441492E-2, 3.01350573947510E-2, 2.01043439592720E-2,  &
         1.18552223068074E-2, 5.86055510956010E-3, 2.25465148267325E-3,  &
         6.08173041536336E-4, 9.84215550625747E-5, 7.32139093038089E-6,  &
         1.37279667384666E-7 /)
    REAL(cdp)    :: an, aq(2), az(2), fi(2), pm, pn, qm, qn, qz(2), rz(2), sz(2), &
         zr(2), ts(2), tm(2), rm(4), sm(4), sn, ss, qf(2), zm, zs
    INTEGER(idp) :: i, m, n


    az(1) = REAL(a)
    az(2) = AIMAG(a)
    zs = az(1) * az(1) + az(2) * az(2)
    zm = SQRT(zs)
    pn = ABS(in)
    sn = +1.0
    IF (az(1) < 0) THEN
       qz(1) = -az(1)
       qz(2) = -az(2)
       IF (in == in/2*2) GO TO 10
       sn = -1.0
    ELSE
       qz(1) = az(1)
       qz(2) = az(2)
    END IF
10  IF (zm > 17.5+0.5*pn*pn) THEN
       qn = pn
    ELSE
       qn = 0.5 * zm - 0.5 * ABS(qz(1)) + 0.5 *ABS(0.5*zm-ABS(qz(1)))
       IF (pn > qn) THEN
          qn = +AINT(0.0625*zs)
          IF (pn <= qn) GO TO 130
          qn = pn
          GO TO 130
       END IF
       IF (zm > 17.5) THEN
          qn = +AINT(SQRT(2.0*(zm-17.5)))
       ELSE
          IF (zs >= 1.0) THEN
             IF (-ABS(az(1))+0.096*az(2)*az(2) >= 0) GO TO 20
          END IF
          qn = AINT(0.0625*zs)
          IF (pn <= qn) GO TO 130
          qn = pn
          GO TO 130
20        qn = 0.0
       END IF
    END IF
    sz(1) = qz(1)
    sz(2) = qz(2)
    qm = sn * 0.398942280401433
    zr(1) = SQRT(sz(1)+zm)
    zr(2) = sz(2) / zr(1)
    zr(1) = 0.707106781186548 * zr(1)
    zr(2) = 0.707106781186548 * zr(2)
    qf(1) = +qm * zr(1) / zm
    qf(2) = -qm * zr(2) / zm
    IF (zm > 17.5) THEN
       rz(1) = +0.5 * qz(1) / zs
       rz(2) = -0.5 * qz(2) / zs
       an = qn * qn - 0.25
       sm(1) = 0.0
       sm(2) = 0.0
       sm(3) = 0.0
       sm(4) = 0.0
       tm(1) = 1.0
       tm(2) = 0.0
       pm = 0.0
       GO TO 40
30     an = an - 2.0 * pm
       pm = pm + 1.0
       ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
       ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
       tm(1) = an * ts(1) / pm
       tm(2) = an * ts(2) / pm
40     sm(1) = sm(1) + tm(1)
       sm(2) = sm(2) + tm(2)
       an = an - 2.0 * pm
       pm = pm + 1.0
       ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
       ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
       tm(1) = an * ts(1) / pm
       tm(2) = an * ts(2) / pm
       IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
          IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 50
       END IF
       sm(3) = sm(3) + tm(1)
       sm(4) = sm(4) + tm(2)
       IF (pm < 35.0) GO TO 30
50     ts(1) = sm(1) + sm(3)
       ts(2) = sm(2) + sm(4)
       sm(1) = sm(1) - sm(3)
       sm(2) = sm(2) - sm(4)
       sm(3) = ts(1)
       sm(4) = ts(2)
    ELSE
       sm(1) = 1.0
       sm(2) = 0.0
       sm(3) = 1.0
       sm(4) = 0.0
       m = 15.0 * qn + 2.0
       n = 15.0 * qn + 15.0
       DO i = m, n
          ts(1) = -qz(1) - cd(i)
          ts(2) = -qz(2)
          ss = ts(1) * ts(1) + ts(2) * ts(2)
          tm(1) = +ce(i) * ts(1) / ss
          tm(2) = -ce(i) * ts(2) / ss
          sm(1) = sm(1) + tm(1)
          sm(2) = sm(2) + tm(2)
          ts(1) = qz(1) - cd(i)
          ts(2) = qz(2)
          ss = ts(1) * ts(1) + ts(2) * ts(2)
          tm(1) = +ce(i) * ts(1) / ss
          tm(2) = -ce(i) * ts(2) / ss
          sm(3) = sm(3) + tm(1)
          sm(4) = sm(4) + tm(2)
       END DO
    END IF
    rm(1) = sm(1)
    rm(2) = sm(2)
    IF (qz(1) < 17.5) THEN
       aq(1) = -2.0 * qz(1)
       IF (qz(2) < 0) THEN
          aq(2) = -2.0 * qz(2) - 3.14159265358979 * (qn+0.5)
       ELSE
          aq(2) = -2.0 * qz(2) + 3.14159265358979 * (qn+0.5)
       END IF
       qm = EXP(aq(1))
       ts(1) = qm * COS(aq(2))
       ts(2) = qm * SIN(aq(2))
       rm(1) = rm(1) + ts(1) * sm(3) - ts(2) * sm(4)
       rm(2) = rm(2) + ts(1) * sm(4) + ts(2) * sm(3)
    END IF
    IF (qn /= pn) THEN
       rm(3) = rm(1)
       rm(4) = rm(2)
       qn = qn + 1.0
       IF (zm > 17.5) THEN
          an = qn * qn - 0.25
          sm(1) = 0.0
          sm(2) = 0.0
          sm(3) = 0.0
          sm(4) = 0.0
          tm(1) = 1.0
          tm(2) = 0.0
          pm = 0.0
          GO TO 80
70        an = an - 2.0 * pm
          pm = pm + 1.0
          ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
          ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
          tm(1) = an * ts(1) / pm
          tm(2) = an * ts(2) / pm
80        sm(1) = sm(1) + tm(1)
          sm(2) = sm(2) + tm(2)
          an = an - 2.0 * pm
          pm = pm + 1.0
          ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
          ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
          tm(1) = an * ts(1) / pm
          tm(2) = an * ts(2) / pm
          IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
             IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 90
          END IF
          sm(3) = sm(3) + tm(1)
          sm(4) = sm(4) + tm(2)
          IF (pm < 35.0) GO TO 70
90        ts(1) = sm(1) + sm(3)
          ts(2) = sm(2) + sm(4)
          sm(1) = sm(1) - sm(3)
          sm(2) = sm(2) - sm(4)
          sm(3) = ts(1)
          sm(4) = ts(2)
       ELSE
          sm(1) = 1.0
          sm(2) = 0.0
          sm(3) = 1.0
          sm(4) = 0.0
          m = 15.0 * qn + 2.0
          n = 15.0 * qn + 15.0
          DO i = m, n
             ts(1) = -qz(1) - cd(i)
             ts(2) = -qz(2)
             ss = ts(1) * ts(1) + ts(2) * ts(2)
             tm(1) = +ce(i) * ts(1) / ss
             tm(2) = -ce(i) * ts(2) / ss
             sm(1) = sm(1) + tm(1)
             sm(2) = sm(2) + tm(2)
             ts(1) = +qz(1) - cd(i)
             ts(2) = +qz(2)
             ss = ts(1) * ts(1) + ts(2) * ts(2)
             tm(1) = +ce(i) * ts(1) / ss
             tm(2) = -ce(i) * ts(2) / ss
             sm(3) = sm(3) + tm(1)
             sm(4) = sm(4) + tm(2)
          END DO
       END IF
       rm(1) = sm(1)
       rm(2) = sm(2)
       IF (qz(1) >= 17.5) GO TO 120
       aq(1) = -2.0 * qz(1)
       IF (qz(2) < 0) THEN
          aq(2) = -2.0 * qz(2) - 3.14159265358979 * (qn+0.5)
       ELSE
          aq(2) = -2.0 * qz(2) + 3.14159265358979 * (qn+0.5)
       END IF
       qm = EXP(aq(1))
       ts(1) = qm * COS(aq(2))
       ts(2) = qm * SIN(aq(2))
       rm(1) = rm(1) + ts(1) * sm(3) - ts(2) * sm(4)
       rm(2) = rm(2) + ts(1) * sm(4) + ts(2) * sm(3)
       GO TO 120
110    tm(1) = -2.0 * qn * qz(1) / zs
       tm(2) = +2.0 * qn * qz(2) / zs
       ts(1) = tm(1) * rm(1) - tm(2) * rm(2) + rm(3)
       ts(2) = tm(1) * rm(2) + tm(2) * rm(1) + rm(4)
       rm(3) = rm(1)
       rm(4) = rm(2)
       rm(1) = ts(1)
       rm(2) = ts(2)
       qn = qn + 1.0
120    IF (qn < pn) GO TO 110
    END IF
    IF (mo == 0) THEN
       qm = EXP(qz(1))
       tm(1) = qm * COS(qz(2))
       tm(2) = qm * SIN(qz(2))
       ts(1) = tm(1) * rm(1) - tm(2) * rm(2)
       ts(2) = tm(1) * rm(2) + tm(2) * rm(1)
       rm(1) = ts(1)
       rm(2) = ts(2)
    END IF
    fi(1) = qf(1) * rm(1) - qf(2) * rm(2)
    fi(2) = qf(1) * rm(2) + qf(2) * rm(1)
    w = CMPLX(fi(1),fi(2))
    RETURN
130 sz(1) = 0.25 * (qz(1)*qz(1)-qz(2)*qz(2))
    sz(2) = 0.5 * qz(1) * qz(2)
    an = qn
    sm(1) = 0.0
    sm(2) = 0.0
    sm(3) = 0.0
    sm(4) = 0.0
    tm(1) = 1.0
    tm(2) = 0.0
    pm = 0.0
140 an = an + 1.0
    ts(1) = tm(1) / an
    ts(2) = tm(2) / an
    sm(3) = sm(3) + ts(1)
    sm(4) = sm(4) + ts(2)
    tm(1) = ts(1) * sz(1) - ts(2) * sz(2)
    tm(2) = ts(1) * sz(2) + ts(2) * sz(1)
    pm = pm + 1.0
    tm(1) = tm(1) / pm
    tm(2) = tm(2) / pm
    IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
       IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 150
    END IF
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
    GO TO 140
150 sm(1) = sm(1) + 1.0
    an = qn + 1.0
    sm(3) = an * sm(3)
    sm(4) = an * sm(4)
    GO TO 170
160 an = qn * (qn+1.0)
    tm(1) = sz(1) / an
    tm(2) = sz(2) / an
    ts(1) = +tm(1) * sm(3) - tm(2) * sm(4)
    ts(2) = +tm(1) * sm(4) + tm(2) * sm(3)
    sm(3) = sm(1)
    sm(4) = sm(2)
    sm(1) = sm(1) + ts(1)
    sm(2) = sm(2) + ts(2)
    qn = qn - 1.0
170 IF (qn > pn) GO TO 160
    qf(1) = sn
    qf(2) = 0.0
    qn = 0.0
    GO TO 190
180 qn = qn + 1.0
    tm(1) = qf(1) * qz(1) - qf(2) * qz(2)
    tm(2) = qf(1) * qz(2) + qf(2) * qz(1)
    qf(1) = 0.5 * tm(1) / qn
    qf(2) = 0.5 * tm(2) / qn
190 IF (qn < pn) GO TO 180
    IF (mo /= 0) THEN
       qm = EXP(-qz(1))
       tm(1) = qm * COS(-qz(2))
       tm(2) = qm * SIN(-qz(2))
       ts(1) = tm(1) * qf(1) - tm(2) * qf(2)
       ts(2) = tm(1) * qf(2) + tm(2) * qf(1)
       qf(1) = ts(1)
       qf(2) = ts(2)
    END IF
    fi(1) = qf(1) * sm(1) - qf(2) * sm(2)
    fi(2) = qf(1) * sm(2) + qf(2) * sm(1)

    w = CMPLX(fi(1),fi(2))

  END SUBROUTINE bssli



  SUBROUTINE bsslj(a, in, w)

    !     ******************************************************************
    !     FORTRAN SUBROUTINE FOR ORDINARY BESSEL FUNCTION OF INTEGRAL ORDER
    !     ******************************************************************
    !     A  = ARGUMENT (COMPLEX NUMBER)
    !     IN = ORDER (INTEGER)
    !     W  = FUNCTION OF FIRST KIND (COMPLEX NUMBER)
    !     -------------------

    COMPLEX(cdp), INTENT(IN)  :: a
    INTEGER(idp), INTENT(IN)  :: in

    COMPLEX(cdp), INTENT(OUT) :: w
    REAL(fdp), PARAMETER :: cd(30) = (/ 0.00000000000000,  &
         -1.64899505142212E-2, -7.18621880068536E-2, -1.67086878124866E-1,  &
         -3.02582250219469E-1, -4.80613945245927E-1, -7.07075239357898E-1,  &
         -9.92995790539516E-1, -1.35583925612592, -1.82105907899132,  &
         -2.42482175310879, -3.21956655708750, -4.28658077248384,  &
         -5.77022816798128, -8.01371260952526, 0.00000000000000,  &
         -5.57742429879505E-3, -4.99112944172476E-2, -1.37440911652397E-1,  &
         -2.67233784710566E-1, -4.40380166808682E-1, -6.61813614872541E-1,  &
         -9.41861077665017E-1, -1.29754130468326, -1.75407696719816,  &
         -2.34755299882276, -3.13041332689196, -4.18397120563729,  &
         -5.65251799214994, -7.87863959810677 /), ce(30)  &
         = (/ 0.00000000000000, -4.80942336387447E-3,  &
         -1.31366200347759E-2, -1.94843834008458E-2, -2.19948900032003E-2,  &
         -2.09396625676519E-2, -1.74600268458650E-2, -1.27937813362085E-2,  &
         -8.05234421796592E-3, -4.15817375002760E-3, -1.64317738747922E-3,  &
         -4.49175585314709E-4, -7.28594765574007E-5, -5.38265230658285E-6,  &
         -9.93779048036289E-8, 0.00000000000000, 7.53805779200591E-2,  &
         7.12293537403464E-2, 6.33116224228200E-2, 5.28240264523301E-2,  &
         4.13305359441492E-2, 3.01350573947510E-2, 2.01043439592720E-2,  &
         1.18552223068074E-2, 5.86055510956010E-3, 2.25465148267325E-3,  &
         6.08173041536336E-4, 9.84215550625747E-5, 7.32139093038089E-6,  &
         1.37279667384666E-7 /)
    REAL(fdp)    :: an, aq(2), az(2), fj(2), pm, pn, qf(2), qm, qn, qz(2), rz(2),  &
         sz(2), rm(4), sm(4), sn, ss, ts(2), tm(2), zm, zr(2), zs
    INTEGER(idp) :: i, m, n


    az(1) = REAL(a)
    az(2) = AIMAG(a)
    zs = az(1) * az(1) + az(2) * az(2)
    zm = SQRT(zs)
    pn = ABS(in)
    sn = 1.0
    IF (in < 0) THEN
       IF (in /= in/2*2) THEN
          sn = -1.0
       END IF
    END IF
    IF (az(1) < 0) THEN
       qz(1) = -az(1)
       qz(2) = -az(2)
       IF (in == in/2*2) GO TO 10
       sn = -sn
    ELSE
       qz(1) = +az(1)
       qz(2) = +az(2)
    END IF
10  IF (zm > 17.5+0.5*pn*pn) THEN
       qn = pn
    ELSE
       qn = 0.5 * zm - 0.5 * ABS(qz(2)) + 0.5 *ABS(0.5*zm-ABS(qz(2)))
       IF (pn > qn) THEN
          qn = +AINT(0.0625*zs)
          IF (pn <= qn) GO TO 130
          qn = pn
          GO TO 130
       END IF
       IF (zm > 17.5) THEN
          qn = +AINT(SQRT(2.0*(zm-17.5)))
       ELSE
          IF (zs >= 1.0) THEN
             IF (-ABS(az(2))+0.096*az(1)*az(1) >= 0) GO TO 20
          END IF
          qn = +AINT(0.0625*zs)
          IF (pn <= qn) GO TO 130
          qn = pn
          GO TO 130
20        qn = 0.0
       END IF
    END IF
    sz(1) = qz(1)
    sz(2) = qz(2)
    qm = sn * 0.797884560802865
    zr(1) = SQRT(sz(1)+zm)
    zr(2) = sz(2) / zr(1)
    zr(1) = 0.707106781186548 * zr(1)
    zr(2) = 0.707106781186548 * zr(2)
    qf(1) = +qm * zr(1) / zm
    qf(2) = -qm * zr(2) / zm
    IF (zm > 17.5) THEN
       rz(1) = +0.5 * qz(1) / zs
       rz(2) = -0.5 * qz(2) / zs
       an = qn * qn - 0.25
       sm(1) = 0.0
       sm(2) = 0.0
       sm(3) = 0.0
       sm(4) = 0.0
       tm(1) = 1.0
       tm(2) = 0.0
       pm = 0.0
       GO TO 40
30     an = an - 2.0 * pm
       pm = pm + 1.0
       ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
       ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
       tm(1) = -an * ts(1) / pm
       tm(2) = -an * ts(2) / pm
40     sm(1) = sm(1) + tm(1)
       sm(2) = sm(2) + tm(2)
       an = an - 2.0 * pm
       pm = pm + 1.0
       ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
       ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
       tm(1) = +an * ts(1) / pm
       tm(2) = +an * ts(2) / pm
       IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
          IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 60
       END IF
       sm(3) = sm(3) + tm(1)
       sm(4) = sm(4) + tm(2)
       IF (pm < 35.0) GO TO 30
    ELSE
       sm(1) = 1.0
       sm(2) = 0.0
       sm(3) = 1.0
       sm(4) = 0.0
       m = 15.0 * qn + 2.0
       n = 15.0 * qn + 15.0
       DO i = m, n
          ts(1) = +qz(2) - cd(i)
          ts(2) = -qz(1)
          ss = ts(1) * ts(1) + ts(2) * ts(2)
          tm(1) = +ce(i) * ts(1) / ss
          tm(2) = -ce(i) * ts(2) / ss
          sm(1) = sm(1) + tm(1)
          sm(2) = sm(2) + tm(2)
          ts(1) = -qz(2) - cd(i)
          ts(2) = +qz(1)
          ss = ts(1) * ts(1) + ts(2) * ts(2)
          tm(1) = +ce(i) * ts(1) / ss
          tm(2) = -ce(i) * ts(2) / ss
          sm(3) = sm(3) + tm(1)
          sm(4) = sm(4) + tm(2)
       END DO
       ts(1) = +0.5 * (sm(2)-sm(4))
       ts(2) = -0.5 * (sm(1)-sm(3))
       sm(1) = +0.5 * (sm(1)+sm(3))
       sm(2) = +0.5 * (sm(2)+sm(4))
       sm(3) = ts(1)
       sm(4) = ts(2)
    END IF
60  aq(1) = qz(1) - 1.57079632679490 * (qn+0.5)
    aq(2) = qz(2)
    ts(1) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
    ts(2) = -SIN(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
    tm(1) = sm(1) * ts(1) - sm(2) * ts(2)
    tm(2) = sm(1) * ts(2) + sm(2) * ts(1)
    ts(1) = +SIN(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
    ts(2) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
    rm(1) = tm(1) - sm(3) * ts(1) + sm(4) * ts(2)
    rm(2) = tm(2) - sm(3) * ts(2) - sm(4) * ts(1)
    IF (qn /= pn) THEN
       rm(3) = rm(1)
       rm(4) = rm(2)
       qn = qn + 1.0
       IF (zm > 17.5) THEN
          an = qn * qn - 0.25
          sm(1) = 0.0
          sm(2) = 0.0
          sm(3) = 0.0
          sm(4) = 0.0
          tm(1) = 1.0
          tm(2) = 0.0
          pm = 0.0
          GO TO 80
70        an = an - 2.0 * pm
          pm = pm + 1.0
          ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
          ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
          tm(1) = -an * ts(1) / pm
          tm(2) = -an * ts(2) / pm
80        sm(1) = sm(1) + tm(1)
          sm(2) = sm(2) + tm(2)
          an = an - 2.0 * pm
          pm = pm + 1.0
          ts(1) = tm(1) * rz(1) - tm(2) * rz(2)
          ts(2) = tm(1) * rz(2) + tm(2) * rz(1)
          tm(1) = +an * ts(1) / pm
          tm(2) = +an * ts(2) / pm
          IF (ABS(sm(3))+ABS(tm(1)) == ABS(sm(3))) THEN
             IF (ABS(sm(4))+ABS(tm(2)) == ABS(sm(4))) GO TO 100
          END IF
          sm(3) = sm(3) + tm(1)
          sm(4) = sm(4) + tm(2)
          IF (pm < 35.0) GO TO 70
       ELSE
          sm(1) = 1.0
          sm(2) = 0.0
          sm(3) = 1.0
          sm(4) = 0.0
          m = 15.0 * qn + 2.0
          n = 15.0 * qn + 15.0
          DO i = m, n
             ts(1) = +qz(2) - cd(i)
             ts(2) = -qz(1)
             ss = ts(1) * ts(1) + ts(2) * ts(2)
             tm(1) = +ce(i) * ts(1) / ss
             tm(2) = -ce(i) * ts(2) / ss
             sm(1) = sm(1) + tm(1)
             sm(2) = sm(2) + tm(2)
             ts(1) = -qz(2) - cd(i)
             ts(2) = +qz(1)
             ss = ts(1) * ts(1) + ts(2) * ts(2)
             tm(1) = +ce(i) * ts(1) / ss
             tm(2) = -ce(i) * ts(2) / ss
             sm(3) = sm(3) + tm(1)
             sm(4) = sm(4) + tm(2)
          END DO
          ts(1) = +0.5 * (sm(2)-sm(4))
          ts(2) = -0.5 * (sm(1)-sm(3))
          sm(1) = +0.5 * (sm(1)+sm(3))
          sm(2) = +0.5 * (sm(2)+sm(4))
          sm(3) = ts(1)
          sm(4) = ts(2)
       END IF
100    aq(1) = qz(1) - 1.57079632679490 * (qn+0.5)
       aq(2) = qz(2)
       ts(1) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
       ts(2) = -SIN(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
       tm(1) = sm(1) * ts(1) - sm(2) * ts(2)
       tm(2) = sm(1) * ts(2) + sm(2) * ts(1)
       ts(1) = +SIN(aq(1)) * 0.5 * (EXP(+aq(2))+EXP(-aq(2)))
       ts(2) = +COS(aq(1)) * 0.5 * (EXP(+aq(2))-EXP(-aq(2)))
       rm(1) = tm(1) - sm(3) * ts(1) + sm(4) * ts(2)
       rm(2) = tm(2) - sm(3) * ts(2) - sm(4) * ts(1)
       GO TO 120
110    tm(1) = +2.0 * qn * qz(1) / zs
       tm(2) = -2.0 * qn * qz(2) / zs
       ts(1) = tm(1) * rm(1) - tm(2) * rm(2) - rm(3)
       ts(2) = tm(1) * rm(2) + tm(2) * rm(1) - rm(4)
       rm(3) = rm(1)
       rm(4) = rm(2)
       rm(1) = ts(1)
       rm(2) = ts(2)
       qn = qn + 1.0
120    IF (qn < pn) GO TO 110
    END IF
    fj(1) = qf(1) * rm(1) - qf(2) * rm(2)
    fj(2) = qf(1) * rm(2) + qf(2) * rm(1)
    w = CMPLX(fj(1), fj(2))
    RETURN
130 sz(1) = 0.25 * (qz(1)*qz(1) - qz(2)*qz(2))
    sz(2) = 0.5 * qz(1) * qz(2)
    an = qn
    sm(1) = 0.0
    sm(2) = 0.0
    sm(3) = 0.0
    sm(4) = 0.0
    tm(1) = 1.0
    tm(2) = 0.0
    pm = 0.0
140 an = an + 1.0
    ts(1) = +tm(1) / an
    ts(2) = +tm(2) / an
    sm(3) = sm(3) + ts(1)
    sm(4) = sm(4) + ts(2)
    tm(1) = -ts(1) * sz(1) + ts(2) * sz(2)
    tm(2) = -ts(1) * sz(2) - ts(2) * sz(1)
    pm = pm + 1.0
    tm(1) = tm(1) / pm
    tm(2) = tm(2) / pm
    IF (ABS(sm(1))+ABS(tm(1)) == ABS(sm(1))) THEN
       IF (ABS(sm(2))+ABS(tm(2)) == ABS(sm(2))) GO TO 150
    END IF
    sm(1) = sm(1) + tm(1)
    sm(2) = sm(2) + tm(2)
    GO TO 140
150 sm(1) = sm(1) + 1.0
    an = qn + 1.0
    sm(3) = an * sm(3)
    sm(4) = an * sm(4)
    GO TO 170
160 an = qn * (qn+1.0)
    tm(1) = sz(1) / an
    tm(2) = sz(2) / an
    ts(1) = -tm(1) * sm(3) + tm(2) * sm(4)
    ts(2) = -tm(1) * sm(4) - tm(2) * sm(3)
    sm(3) = sm(1)
    sm(4) = sm(2)
    sm(1) = sm(1) + ts(1)
    sm(2) = sm(2) + ts(2)
    qn = qn - 1.0
170 IF (qn > pn) GO TO 160
    qf(1) = sn
    qf(2) = 0.0
    qn = 0.0
    GO TO 190
180 qn = qn + 1.0
    tm(1) = qf(1) * qz(1) - qf(2) * qz(2)
    tm(2) = qf(1) * qz(2) + qf(2) * qz(1)
    qf(1) = 0.5 * tm(1) / qn
    qf(2) = 0.5 * tm(2) / qn
190 IF (qn < pn) GO TO 180
    fj(1) = qf(1) * sm(1) - qf(2) * sm(2)
    fj(2) = qf(1) * sm(2) + qf(2) * sm(1)

    w = CMPLX(fj(1),fj(2))

  END SUBROUTINE bsslj



END MODULE nswclib
