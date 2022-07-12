C**** Subroutines in misc.f (copied from MISCSUBS.FOR (VAX), D.J.Kerridge)
C     
C     List last revised 25 January 1999 : E.Clarke
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
C     SUBROUTINE DMS(JD,JM,JS,DEG,IC)
C     SUBROUTINE DMSP(JD1,JM1,JS1,JD2,JM2,JS2,JD3,JM3,JS3,DEG3)
C     SUBROUTINE DMSM(JD1,JM1,JS1,JD2,JM2,JS2,JD3,JM3,JS3,DEG3)
C     SUBROUTINE DAMON(ID,IM,IY,ND,MODE)
C     SUBROUTINE ETRAP(ECODE,ID,IM,IY,ND,MODE)
C     SUBROUTINE DAYSIN(MON,IYEAR,NDAYS)
C     SUBROUTINE DAYSINMON(IY,NDM)
C     FUNCTION   IDAYSINYR(IY)
C     SUBROUTINE DOW(ID,M,IY,DAY)
C     SUBROUTINE MNVAR(ARR,NPTS,FLAG,RMN,VAR,NMISS)
C     SUBROUTINE MNMX(ARR,NPTS,FLAG,RMIN,RMAX)
C     SUBROUTINE IDIVREM(NUM,IDIV,IQUO,IREM)
C     SUBROUTINE GETLEN(STRING,NCHAR)
C     SUBROUTINE CHCASE(STRING,CASE)
C     SUBROUTINE RMEAN(ARR,AFIL,NARR,RMISS,FILT,NFILT)
C     SUBROUTINE BROT(IDAY,MON,IYR,IROT,IDNO)
C     FUNCTION ROUND(X1,X2)
C     FUNCTION IROUND(I1,I2)
C     SUBROUTINE FDMON(IY,IFDM)
C     SUBROUTINE AUTOC(X,NX,LAG,ACF,XVR)
C     SUBROUTINE LASTMON(ID,IM,IY,ND,NDYR)
C     SUBROUTINE THISMON(ID,IM,IY,ND,NDYR)
C     SUBROUTINE PURGER(FILE)
C     SUBROUTINE DATERD(ID,IM,IY,ND,NDYR)
c     subroutine tdaterd(id,im,iy,nd,ndyr)
c     subroutine stat(x,n,flag,xmin,xmax,xrms,xsd,nmiss)
C
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine DMS to convert from degrees (JD), minutes (JM) and
C     seconds (JS) to decimal degrees (DEG) when IC = 1, and from decimal
C     degrees to degrees minutes and seconds when IC = 2.  
C
C**** No 'allowed range' checks are made on the input data, so the input
C     can be in minutes or seconds for conversion to decimal degrees.
C     eg. JD=JS=0 and JM=420 will lead to DEG=7.0
C
C     The output format for negative angles is:
C     a. -12 34 56 when the degrees are non-zero
C     b.  00-34 56 when degrees = 0
C     c   00 00-56 when degrees = 0 and minutes = 0
C
C**** The input parameters are unchanged on exit.
C
C**** Last revised 5th August 1987  DJK
C
      SUBROUTINE DMS(JD,JM,JS,DEG,IC)
      INTEGER IC,ISIG,JD,JM,JS
      REAL DEG,DG,DR,EM,SC,SIG
C
      IF(IC .EQ. 2) GOTO 1
      IF(IC .NE. 1) THEN
         WRITE(6,*) 'IC should be 1 or 2, it has been passed as',IC
         STOP
      END IF
C
C**** IC = 1, so convert from dms to decimal degrees
      SIG = 1.0
      IF(JD .LT. 0 .OR. JM .LT. 0 .OR. JS .LT. 0) SIG = -1.0
      DG  = IABS(JD)
      EM  = IABS(JM)
      SC  = IABS(JS)
      DEG = (DG + EM/60. + SC/3600.)*SIG
      RETURN
C
C**** IC = 2, so convert from decimal degrees to dms
    1 SIG = SIGN(1.1,DEG)
      DR  = ABS(DEG)
      JD  = INT(DR)
      DG  = JD
      JM  = INT(60.*(DR-DG))
      EM  = JM
      JS  = INT((DR - DG - EM/60.)*3600. + .5)
C
      IF(JS .EQ. 60) THEN
         JS = 0
         JM = JM + 1
      END IF
C
      IF(JM .EQ. 60) THEN
         JM = 0
         JD = JD+1
      END IF
C
C**** Assign the correct sign to the result
C
      ISIG = INT(SIG)
      IF(JD .NE. 0) THEN
         JD = JD*ISIG
      ELSE
             IF(JM .NE. 0) THEN
                JM = JM*ISIG
             ELSE
                JS = JS*ISIG
             END IF
      END IF
C
      RETURN
      END
C
C**** Subroutine DMSP to add two angles supplied in degrees, minutes and
C     seconds.  The result, (angle3 = angle1+angle2), is returned in both
C     degrees minutes and seconds and decimal degrees.
C
C**** This subroutine calls subroutine DMS
C
C**** Last revised 22nd June 1987 : DJK
C
      SUBROUTINE DMSP(JD1,JM1,JS1,JD2,JM2,JS2,JD3,JM3,JS3,DEG3)
      INTEGER JD1,JM1,JS1,JD2,JM2,JS2,JD3,JM3,JS3
      REAL DEG1,DEG2,DEG3
C
      CALL DMS(JD1,JM1,JS1,DEG1,1)
      CALL DMS(JD2,JM2,JS2,DEG2,1)
      DEG3 = DEG1 + DEG2
      CALL DMS(JD3,JM3,JS3,DEG3,2)
      RETURN
      END
C
C**** Subroutine DMSM to subtract two angles supplied in degrees, minutes
C     and seconds.  The result, (angle3 = angle1-angle2), is returned in
C     both degrees minutes and seconds and decimal degrees.
C
C     This subroutine calls subroutine DMS
C
C**** Last revised 22nd June 1987 : DJK
C
      SUBROUTINE DMSM(JD1,JM1,JS1,JD2,JM2,JS2,JD3,JM3,JS3,DEG3)
      INTEGER JD1,JM1,JS1,JD2,JM2,JS2,JD3,JM3,JS3
      REAL DEG1,DEG2,DEG3
C
      CALL DMS(JD1,JM1,JS1,DEG1,1)
      CALL DMS(JD2,JM2,JS2,DEG2,1)
      DEG3 = DEG1 - DEG2
      CALL DMS(JD3,JM3,JS3,DEG3,2)
      RETURN
      END
C
C**** Subroutine DAMON to convert day (ID), month (IM), year (IY) to day
C     number when MODE=1 and to compute day and month from day number and
C     year when MODE=2.  The full year should be specified, eg. 1988 not
C     just 88.  In the Gregorian calendar, introduced in 1582 but not
C     adopted in English-speaking countries until 1752, years which are
C     are divisible by 100 are leap years only if they are also divisible
C     by 400.
C
C**** Last revised 11th February 1988 : DJK (adapted from PRR subroutine)
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      SUBROUTINE DAMON(ID,IM,IY,ND,MODE)
      INTEGER MD,ID,IM,IY,ND,MODE
      DIMENSION MD(13)
      DATA MD /0,31,59,90,120,151,181,212,243,273,304,334,365/
C
      IF(MODE .NE. 1 .AND. MODE .NE. 2) CALL ETRAP(1,ID,IM,IY,ND,MODE)
C
      IF(MODE .EQ. 1) THEN
         IF(ID .GT. 31 .OR. ID .LT. 1)  CALL ETRAP(2,ID,IM,IY,ND,MODE)
         IF(IM .GT. 12 .OR. IM .LT. 1)  CALL ETRAP(3,ID,IM,IY,ND,MODE)
      	 IF(IY .LT. 1582) CALL ETRAP(4,ID,IM,IY,ND,MODE)
         IF(ID .GT. MD(IM+1)-MD(IM)) THEN
            IF(IM .NE. 2 .OR. (IM .EQ. 2 .AND. ID .NE. 29))
     &         CALL ETRAP(5,ID,IM,IY,ND,MODE)
            IF((IM .EQ. 2 .AND. ID .EQ. 29) .AND. (MOD(IY,100) .EQ. 0
     &          .AND. MOD(IY,400) .NE. 0)) 
     &          CALL ETRAP(5,ID,IM,IY,ND,MODE)
         END IF
      ELSE
         IF(ND .GT. 366 .OR. ND .LT. 1) CALL ETRAP(6,ID,IM,IY,ND,MODE)
      END IF
C
      IF(MODE .EQ. 1) THEN
         ND  = MD(IM) + ID
         IF(MOD(IY,4) .EQ. 0 .AND. IM .GT. 2) THEN
            IF(MOD(IY,100) .EQ. 0 .AND. MOD(IY,400) .NE. 0) THEN
            ELSE
              ND = ND+1
            END IF
         END IF            
      ELSE
         L = 0
         IF(MOD(IY,4) .EQ. 0 .AND. ND .GT. 59) THEN
            IF(MOD(IY,100) .EQ. 0 .AND. MOD(IY,400) .NE. 0) THEN
            ELSE
              L = 1
            END IF
         END IF
         NND = ND-L
         IF(NND .EQ. 366) CALL ETRAP(7,ID,IM,IY,ND,MODE)
         DO 1 IL=12,1,-1
         IM  = IL
         ID  = NND-MD(IM)
         IF(ID .GE. 1) GOTO 2
    1    CONTINUE
    2    IF(IM .EQ. 2) ID = ID+L
      END IF
C
      RETURN
      END
C
C**** Subroutine ETRAP for DAMON
C
      SUBROUTINE ETRAP(ECODE,ID,IM,IY,ND,MODE)
      INTEGER ECODE,MODE,ID,IM,IY,ND
C
      GOTO (10,20,30,40,50,60,70) ECODE
C
   10 WRITE(6,*) ' MODE must be either 1 or 2'
      WRITE(6,*) ' MODE =',MODE,'  has been passed to DAMON'
      STOP
   20 WRITE(6,*) 'The day of the month has been passed as',ID
      STOP
   30 WRITE(6,*) 'The month has been passed as',IM
      STOP
   40 WRITE(6,*) 'The Gregorian calendar began in 1582'
      WRITE(6,*) 'The year has been passed as',IY
      STOP
   50 WRITE(6,*) 'The day number and month are inconsistent',ID,IM
      STOP
   60 WRITE(6,*) 'The day number has been passed as',ND
      STOP
   70 WRITE(6,*) 'There is no day 366 in year',IY
      STOP
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine DAYSIN to return the number of days in month MON in
C     year IYEAR
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      SUBROUTINE DAYSIN(MON,IYEAR,NDAYS)
      INTEGER N
      DIMENSION MDAYS(12)
      DATA MDAYS /31,-99,31,30,31,30,31,31,30,31,30,31/
C
      N=365
      CALL DAMON(31,12,IYEAR,N,1)
      IF(N .EQ. 365) THEN
         MDAYS(2) = 28
      ELSE IF (N .EQ. 366) THEN
         MDAYS(2) = 29
      ELSE
         WRITE(6,'(1X,A)') 'Error in subroutine DAYSIN'
         STOP
      END IF
      NDAYS = MDAYS(MON)
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine DAYSINMON to return the number of days in each month in 
C     year IY in the array NDM.  Function IDAYSINYR is called.
C
C**** 26th September 1990
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      SUBROUTINE DAYSINMON(IY,NDM)
      IMPLICIT NONE
      INTEGER NDM,IY,ND,IDAYSINYR
      DIMENSION NDM(12)
C
      NDM(1)  =  31
      NDM(2)  = -99
      NDM(3)  =  31
      NDM(4)  =  30
      NDM(5)  =  31
      NDM(6)  =  30
      NDM(7)  =  31
      NDM(8)  =  31
      NDM(9)  =  30
      NDM(10) =  31
      NDM(11) =  30
      NDM(12) =  31
C
      ND = IDAYSINYR(IY)
      IF(ND .EQ. 365) THEN
         NDM(2) = 28
      ELSE IF(ND .EQ. 366) THEN
         NDM(2) = 29
      ELSE
         WRITE(6,'(/1X,A/)') '>>> Error in subroutine DAYSINMON <<<'
         STOP
      END IF
C
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine FDMON to return the day number for the first day in each
C     month in year IY in the array NDM.  Subroutine DAYSINMON is called.
C
C**** 26th September 1990
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      SUBROUTINE FDMON(IY,IFDM)
      IMPLICIT NONE
      INTEGER IFDM,NDM,IY,IL
      DIMENSION NDM(12),IFDM(12)
C
      CALL DAYSINMON(IY,NDM)
      IFDM(1) = 1
      DO 10 IL = 2,12
         IFDM(IL) = IFDM(IL-1) + NDM(IL-1)
   10 CONTINUE
C
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine DOW to compute the day of the week given the date in
C     day (ID), month (M) and year (IY) format.  The JULDAY function
C     from 'Numerical Recipes' is used.  The day is written into the
C     character variable DAY.  
C     (Note: the JULDAY function can be used independently of the 
C      subroutine to obtain the Julian day number.)
C
C**** Last revised 26th September 1988 : DJK
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      SUBROUTINE DOW(ID,M,IY,DAY)
      INTEGER ID,IY,M,JDAY
      CHARACTER*9 DAYS(7),DAY
      DATA DAYS /'Sunday','Monday','Tuesday','Wednesday','Thursday',
     &           'Friday','Saturday'/
C
      JDAY = JULDAY(ID,M,IY)
      IREM = MOD(JDAY+1,7) + 1
      WRITE(DAY,'(A)') DAYS(IREM)
      RETURN
      END
C
      FUNCTION JULDAY(ID,MM,IYYY)
      PARAMETER (IGREG=15+31*(10+12*1582))
      IF(IYYY.EQ.0) THEN
         WRITE(6,'(1X,A)') 'There is no Year Zero'
         STOP
      END IF
C
      IF (IYYY.LT.0) IYYY=IYYY+1
      IF (MM.GT.2) THEN
        JY=IYYY
        JM=MM+1
      ELSE
        JY=IYYY-1
        JM=MM+13
      ENDIF
      JULDAY=INT(365.25*JY) + INT(30.6001*JM) + ID + 1720995
      IF (ID+31*(MM+12*IYYY) .GE. IGREG) THEN
        JA=INT(0.01*JY)
        JULDAY=JULDAY + 2 - JA + INT(0.25*JA)
      ENDIF
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Function IDAYSINYR to return the number of days in the year IY.
C     Subroutine DAMON is used.
C
C**** 26th September 1990
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      FUNCTION IDAYSINYR(IY)
      IMPLICIT NONE
      INTEGER IDAYSINYR,IY,IDAY,MON,MODE,ND
      DATA IDAY,MON,MODE /31,12,1/
C
      CALL DAMON(IDAY,MON,IY,ND,MODE)
      IDAYSINYR = ND
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine MNVAR to calculate the mean (RMN) and variance (VAR)
C     of the NPTS elements passed to ARR(NPTS).
C     FLAG is the value used to flag missing data.
C
C*CC**** Last revised 12-06-1989 : DJK
C
      SUBROUTINE MNVAR(ARR,NPTS,FLAG,RMN,VAR,NMISS)
      IMPLICIT NONE
      INTEGER IL,NPTS,NST,IOUT,NMISS
      REAL ARR,DEL,RMN,VAR,FLAG,AMN,TINY
      DIMENSION ARR(NPTS)
      PARAMETER (IOUT = 6, TINY = 0.01)
C      
      IF(NPTS .EQ. 1) THEN
         RMN =  ARR(1)
         VAR =  FLAG
      ELSE
         RMN   = 0.0
         NMISS =   0
         NST   =   0
   10    NST   = NST + 1
         IF(NST .EQ. NPTS+1) THEN
            RMN = FLAG
            VAR = FLAG
            RETURN
         END IF
         DEL   = ABS(ARR(NST) - FLAG)
         IF(DEL .LT. TINY) THEN
            NMISS = NMISS + 1
            GOTO 10
         ELSE
            AMN   = ARR(NST)
         END IF
C
         DO 20 IL = NST+1,NPTS
            DEL = ABS(ARR(IL) - FLAG)
            IF(DEL .GT. TINY) THEN
               RMN = RMN + ARR(IL)-AMN
            ELSE
               NMISS = NMISS + 1
            END IF
   20    CONTINUE
         RMN = AMN + RMN/REAL(NPTS-NMISS)
C
         VAR = 0.0
         DO 30 IL=1,NPTS
            DEL = ABS(ARR(IL) - FLAG)
            IF(DEL .GT. 0.01) VAR = VAR + (RMN-ARR(IL))*(RMN-ARR(IL))
   30    CONTINUE
         IF(NPTS-NMISS-1 .GT. 0) THEN
            VAR = VAR/REAL(NPTS-NMISS-1)
         ELSE
            VAR = FLAG
         END IF
C         IF(NMISS .GT. NPTS/10) THEN
C             VAR = FLAG
C             RMN = FLAG
C          ENDIF
      END IF
C
      RETURN
      END
C
C**** Subroutine MNMX to find the smallest (RMIN) and largest (RMAX)
C     elements in the NPTS passed into array ARR(NPTS).  
C     FLAG is the value used to flag missing data.
C
C**** Last revised 27-09-1990 : DJK
C
      SUBROUTINE MNMX(ARR,NPTS,FLAG,RMIN,RMAX)
      IMPLICIT NONE
      INTEGER IL,NPTS,NDAT
      REAL ARR,DEL,RMIN,RMAX,FLAG,TINY
      DIMENSION ARR(NPTS)
      DATA TINY /0.1/
C
      RMIN = FLAG
      RMAX = FLAG
      NDAT = 0
      DO 10 IL=1,NPTS
         DEL  = ABS(ARR(IL) - FLAG)
         IF(DEL .GT. TINY) THEN
            NDAT = NDAT + 1
            IF(NDAT .EQ. 1) THEN
               RMIN = ARR(IL)
               RMAX = ARR(IL)
            ELSE
               RMIN  = MIN(RMIN,ARR(IL))
               RMAX  = MAX(RMAX,ARR(IL))
            END IF
         END IF
   10 CONTINUE
C
      RETURN
      END
C
C**** Subroutine IDIVREM to perform integer division and return a 
C     remainder.  IDIV goes into NUM IQUO times with remainder IREM.
C
C**** 29th January 1988 : DJK
C
      SUBROUTINE IDIVREM(NUM,IDIV,IQUO,IREM)
      INTEGER NUM,IDIV,IQUO,IREM
C
      IQUO = NUM/IDIV
      IREM = MOD(NUM,IDIV)
C
      RETURN
      END
C
C**** Subroutine GETLEN to find the number of 'significant' characters,
C     NCHAR, in STRING.  NCHAR is found by counting back from the end
C     of STRING until a non-space character is found.
C
C**** 23rd August1988 : DJK
C
      SUBROUTINE GETLEN(STRING,NCHAR)
      INTEGER NCHAR
      CHARACTER*(*) STRING
      NCHAR = LEN(STRING) + 1
   10 NCHAR = NCHAR - 1
      IF(STRING(NCHAR:NCHAR) .NE. ' ') RETURN
      GOTO 10
      END
C
C**** Subroutine CHCASE to change the case of a character variable
C     STRING according to the specifier CASE.  Only alphabetic 
C     characters are affected.  A string with a mixture of upper and
C     lower case on entry will leave with all characters as specified
C     by CASE.
C
C**** Last revised 11th August 1988 : DJK
C
      SUBROUTINE CHCASE(STRING,CASE)
      CHARACTER*(*) STRING
      CHARACTER*5 CASE
C
      IF(CASE .NE. 'UPPER' .AND. CASE .NE. 'lower') THEN
         WRITE(6,'(/1X,4A/)') 'Case conversion specifier passed to',
     &   ' CHCASE as ''',CASE,''''
         WRITE(6,'(1X,A/)') 'It must be either ''UPPER'' or ''lower'''
         STOP
      END IF
C
      NCHAR = LEN(STRING)
      IF(CASE .EQ. 'UPPER') THEN
         DO 10 IL = 1,NCHAR
            IASC  = ICHAR(STRING(IL:IL))
            IF(IASC .GE. 97 .AND. IASC .LE. 122) THEN
               IASC = IASC - 32
               STRING(IL:IL) = CHAR(IASC)
            END IF
   10       CONTINUE
      ELSE
         DO 20 IL = 1,NCHAR
            IASC  = ICHAR(STRING(IL:IL))
            IF(IASC .GE. 65 .AND. IASC .LE. 90) THEN
               IASC = IASC + 32
               STRING(IL:IL) = CHAR(IASC)
            END IF
   20       CONTINUE
      END IF
C
      RETURN
      END
C
C**** Subroutine RMEAN to filter (running mean) the data in array ARR(NARR).
C     The filter coefficients are supplied in array FILT(NFILT) and the filtered
C     data are written into array AFIL(NARR).  RMISS is the value used to show
C     missing data and this value is put into AFIL in the element at the 
C     centre of the running mean.
C
C**** Last revised 12th September 1988 : DJK
C
      SUBROUTINE RMEAN(ARR,AFIL,NARR,RMISS,FILT,NFILT)
      INTEGER IL,ISW,IWT,JL,NMISS,NFILT,NST,NFN,NARR
      REAL ARR,AFIL,FILT,RM,RMISS
      DIMENSION ARR(NARR),AFIL(NARR),FILT(NFILT)
C
C**** Check for a symmetric filter
      IF(MOD(NFILT,2) .EQ. 0) THEN
         WRITE(6,'(/1X,2A)') 'The number of filter coefficients ',
     &                       'must be odd'
         WRITE(6,'(1X,A,I4/)') 'The value supplied is',NFILT
         STOP
      END IF
C
      NMISS = NFILT/2
      NST   = NMISS + 1
      NFN   = NARR  - NMISS
      DO 20 IL = 1,NARR
         RM   = 0.0
         ISW  =   0
         IWT  =   0
         IF(IL.GE.NST .AND. IL.LE.NFN) THEN
            DO 10 JL = IL-NMISS,IL+NMISS
               IWT = IWT + 1
               IF(ABS(ARR(JL)-RMISS) .LT. 0.1) ISW = 1
               RM  = RM + ARR(JL)*FILT(IWT)
   10          CONTINUE
         AFIL(IL)  = RM
C**** ISW=1 indicates that data is missing in ARR within the range of the
C     filter. The filtered value is written as RMISS.
         IF(ISW .EQ. 1) AFIL(IL) = RMISS
         ELSE
            AFIL(IL) = RMISS
         END IF
   20    CONTINUE
C
      RETURN
      END
C
C**** Subroutine BROT to calculate the Bartels rotation number and the
C     day number in the rotation for the date supplied.
C
C     Day 1 of rotation 1 is 8th February 1832.
C
C**** Last revised 13th February 1990 : DJK
C**** Corrected 6th Jan 1992 (it had done multiples of 27 incorrectly)
C
C**** Y2K initial check : 25th January 1999, E.Clarke
C**** Y2K final check   : 
C
      SUBROUTINE BROT(IDAY,MON,IYR,IROT,IDNO)
C
      INTEGER NDAY0,NDAY1,NDAY2

      NDAY0=0
      NDAY1=0
      NDAY2=0
      
      CALL DAMON(7,2,1832,NDAY0,1)
      IYR1 = IYR
      CALL DAMON(IDAY,MON,IYR,NDAY2,1)
      IF(IYR1 .EQ. 1832) THEN
         NDAY1 = 0
         NDAY2 = NDAY2 - NDAY0
      ELSE
         CALL DAMON(31,12,1832,NDAY1,1)
         NDAY1 = NDAY1 - NDAY0
      END IF
C
      NDAY4 = 0
      IF(IYR .GT. 1832) THEN
         DO 10 IL = 1833,IYR-1
            CALL DAMON(31,12,IL,NDAY3,1)
            NDAY4 = NDAY4 + NDAY3
   10    CONTINUE
      END IF
C
      NDAYT = NDAY1 + NDAY4 + NDAY2
      CALL IDIVREM(NDAYT,27,IROT,IREM)
      IF(IREM .EQ. 0) THEN
         IDNO = 27
      ELSE
         IROT = IROT + 1
         IDNO = IREM
      END IF
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Function ROUND to return the closest multiple of X2 to X1,
C     (eg. X1=1234.5, X2=10, then ROUND=1230.0)
C
C**** Last revised 20th April 1990
C
      FUNCTION ROUND(X1,X2)
      REAL ROUND,X1,X2,X3,S
      S  = SIGN(1.,X1)
      X1 = ABS(X1)
      X3 = MOD(X1,X2)
      ROUND = X1-X3
      IF(X3 .GT. X2/2) ROUND = ROUND+X2
      X1    = X1*S
      ROUND = ROUND*S
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Function IROUND to return the closest multiple of I2 to I1,
C     (eg. I1=1234, I2=10, then IROUND=1230)
C
C**** Last revised 20th April 1990
C
      FUNCTION IROUND(I1,I2)
      INTEGER IROUND,I1,I2,I3,IS
      IS = SIGN(1,I1)
      I1 = ABS(I1)
      I3 = MOD(I1,I2)
      IROUND = I1-I3
      IF(I3 .GT. I2/2) IROUND = IROUND+I2
      I1     = I1*IS
      IROUND = IROUND*IS
      RETURN
      END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
C
C**** Subroutine AUTOC to compute the autocorrelation coefficients for 
C     the NX data in array X up from lag 1 to LAG.  The coefficients are 
C     returned in array ACF
C
C**** 21 Nov 1991 : DJK
C
      SUBROUTINE AUTOC(X,NX,LAG,ACF,XVR)
      IMPLICIT NONE
      INTEGER IL,KL,NX,LAG,NT
      DOUBLE PRECISION XMN,XVR,AMN,ACOV
      REAL X,ACF
      DIMENSION X(NX),ACF(LAG)
C
C**** Calculate the mean and variance of the data
      XMN = 0.D0
      XVR = 0.D0
      AMN = DBLE(X(1))
      DO 10 IL = 2,NX
         XMN   = XMN + DBLE(X(IL))-AMN
   10    CONTINUE
         XMN   = AMN + XMN/DBLE(NX)
C
      DO 20 IL = 1,NX
         XVR   = XVR + (DBLE(X(IL))-XMN)**2
   20 CONTINUE
C
C**** Now compute the autocorrelation coefficients
      DO 40 KL  = 1,LAG
        ACOV    = 0.D0
        NT      = NX - KL
           DO 30 IL = 1,NT
              ACOV  = ACOV + (DBLE(X(IL))-XMN)*(DBLE(X(IL+KL))-XMN)
   30      CONTINUE
        ACF(KL) = REAL(ACOV/XVR)
   40   CONTINUE
C
        RETURN
        END
C>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

C*******************************************************************************
C Routine:   LEAP2     Author: AWP Thomson         Date: 20-07-94 
C
C Function:  Gets number of days in specified year 
C
C Notes:     FIX FOR Y2K - WON'T WORK IN 2090-2100 !! WILL WORK FOR 4 DIGIT YEAR
C
C*******************************************************************************
      SUBROUTINE LEAP2(YEAR,IDAYSINYR)
      IMPLICIT   NONE
      INTEGER    YEAR,IDAYSINYR
      LOGICAL    LEAPYEAR
      IF(YEAR.LT.100)  YEAR=YEAR+1900
CY2K
      IF(YEAR.LE.1990) YEAR=YEAR+100
CY2K
      LEAPYEAR=.FALSE.
      IF(MOD(YEAR,4).NE.0) THEN
           LEAPYEAR =.FALSE.
      ELSE
           IF(MOD(YEAR,100).NE.0) THEN
            LEAPYEAR=.TRUE.
           ELSE 
            IF(MOD(YEAR,400).NE.0) THEN
              LEAPYEAR=.FALSE.
            ELSE
              LEAPYEAR=.TRUE.
            END IF
           END IF
      END IF
      IDAYSINYR=365
      IF(LEAPYEAR) IDAYSINYR=366
      RETURN
      END
C*******************************************************************************
C Routine:   DAYNUMTOMON     Author: AWP Thomson         Date: 20-07-94
C
C Function:  Gets month and day from year and daynumber
C
C Notes:     Checks supplied daynumber fits expected number of days for year.
C Notes:     FIX FOR Y2K - MAY NOT WORK IN 2099-2100 !!!
C
C*******************************************************************************
      SUBROUTINE DAYNUMTOMON(YERIN,DAYIN,YER,MON,DAY)
      IMPLICIT   NONE
      INTEGER    YER,DAYNUM,YERNUM,MON,DAY,IDAYSINYR,I,DAYINMON1,
     1           YERIN,DAYIN,DAYINMON
      DIMENSION  DAYINMON1(12),DAYINMON(12)
      DATA       DAYINMON1/31,59,90,120,151,181,
     1                    212,243,273,304,334,365/

      YERNUM=YERIN
      DAYNUM=DAYIN
      YER=0
      MON=0
      DAY=0
      CALL LEAP2(YERNUM,IDAYSINYR)
      IF(DAYNUM.LE.0) THEN
       YERNUM=YERNUM-1
CY2K
       IF(YERNUM.LT.0) YERNUM=YERNUM+100
CY2K
       CALL LEAP2(YERNUM,IDAYSINYR)
       DAYNUM=DAYNUM+IDAYSINYR
      END IF
      IF(DAYNUM.GT.IDAYSINYR) THEN
       YERNUM=YERNUM+1
CY2K
       IF((YERNUM.GE.100).AND.(YERNUM.LT.1900))YERNUM=YERNUM-100
CY2K
       DAYNUM=DAYNUM-IDAYSINYR
      END IF
      YER=YERNUM
      CALL LEAP2(YERNUM,IDAYSINYR)
      DAYINMON(1)=DAYINMON1(1)
      IF(IDAYSINYR.GT.365) THEN
       DO 10 I=2,12
        DAYINMON(I)=DAYINMON1(I)+1
10     CONTINUE
      ELSE
       DO 11 I=2,12
        DAYINMON(I)=DAYINMON1(I)
11     CONTINUE
      END IF
      DO 20 I=1,12
       IF(DAYNUM.LE.DAYINMON(I)) THEN
        MON=I
        IF(I.GT.1) THEN
         DAY=DAYNUM-DAYINMON(I-1)
        ELSE
         DAY=DAYNUM
        END IF
        GOTO 30
       END IF
20    CONTINUE
30    CONTINUE
      RETURN
      END
C*******************************************************************************
C Routine:   MONTODAYNUM     Author: AWP Thomson         Date: 20-07-94
C
C Function:  Gets year and daynumber from day, month and year
C
C Notes:    
C
C*******************************************************************************
      SUBROUTINE MONTODAYNUM(YER,MON,DAY,YEROUT,DAYOUT)
      IMPLICIT   NONE
      INTEGER    YER,DAYNUM,YERNUM,MONNUM,
     1           MON,DAY,IDAYSINYR,I,DAYINMON1,
     1           YEROUT,DAYOUT,DAYINMON
      DIMENSION  DAYINMON1(12),DAYINMON(12)
      DATA       DAYINMON1/31,59,90,120,151,181,
     1                    212,243,273,304,334,365/

      YERNUM=YER
      MONNUM=MON
      DAYNUM=DAY
      YEROUT=0
      DAYOUT=0
      CALL LEAP2(YERNUM,IDAYSINYR)
      DAYINMON(1)=DAYINMON1(1)
      IF(IDAYSINYR.GT.365) THEN
       DO 10 I=2,12
        DAYINMON(I)=DAYINMON1(I)+1
10     CONTINUE
      ELSE
       DO 11 I=2,12
        DAYINMON(I)=DAYINMON1(I)
11     CONTINUE
      END IF
      YEROUT=YERNUM
      IF(MONNUM.GT.1) THEN
       DAYOUT=DAYINMON(MONNUM-1)+DAYNUM
      ELSE
       DAYOUT=DAYNUM
      ENDIF
      RETURN
      END
