C
C**********************************************************************
C Automated fault diagnosis of chemical process plants using model-based reasoning
C Download:
C http://dspace.mit.edu/bitstream/handle/1721.1/14194/22192007.pdf?sequence=1
C Author: Finch, Francis Eric
C Citable URI: http://hdl.handle.net/1721.1/14194
C Other Contributors: Massachusetts Institute of Technology.
C Dept. of Chemical Engineering.
C Advisor: Mark A. Kramer.
C Department: Massachusetts Institute of Technology.
C Dept. of Chemical Engineering.
C Publisher: Massachusetts Institute of Technology
C Date Issued: 1989
C Description:
C Thesis (Ph. D.)--Massachusetts Institute of Technology,
C Dept. of Chemical Engineering, 1989.
C Science hard copy bound in 2 v.Includes bibliographical
C references C (leaves 300-307).
C URI: http://hdl.handle.net/1721.1/14194
C Keywords: Chemical Engineering.
C**********************************************************************
C
C
C
C Jacketed CSTR Simulation Program:
C
C23456789        ** JACKETED CSTR DYNAMIC SIMULATION PROGRAM **
C
C          F. ERIC FINCH
C          DEPARTMENT OF CHEMICAL ENGINEERING
C            MASSACHUSETTS INSTITUTE OF TECHNOLOGY
C          CAMBRIDGE, MA 02139
C
C           REVISED 5/89
C
C          ** VARIABLE DEFINITIONS **
C
C ALPHA1                - PRIMARY RATE CONSTANT PRE-EXPONENTIAL FACTOR (1/MIN)
C BETA1                 - PRIMARY ACTIVATION ENERGY (KJ/KMOL)
C ALPHA2                - SECONDARY RATE CONSTANT PRE-EXPONENTIAL FACTOR (1/MIN)
C BETA2                 - SECONDARY ACTIVATION ENERGY (KJ/KMOL)
C B(3,3)                - ARRAY OF CONTROLLER CONSTANTS FOR PID CONTROL
C                         (IF B(1,2)=0 & B(1,3)=0 THEN P CONTROL,
C                         IF B(1,3)=0 THEN PI CONTROL)
C CA0                   - FEED CONCENTRATION (A) (KMOL/M3)
C CA                    - REACTOR CONCENTRATION (A) (KMOL/M3)
C CB                    - REACTOR CONCENTRATION (B) (KMOL/M3)
C CC                    - REACTOR CONCENTRATION (C) (KMOL/M3)
C CNT(3)                - CONTROLLER OUTPUT VECTOR
C CP,CPOLD              - HEAT CAPACITY OF REACTOR CONTENTS (KJ/KG C)
C CP1                   - FEED HEAT CAPACITY (KJ/KG C)
C CP2                   - COOLING WATER HEAT CAPACITY (KJ/KG C)
C CWPD                  - COOLING WATER PRESSURE DROP CONSTRAINT RESIDUAL (M3/MIN)
C DATAFILE              - FILE NAME FOR OUTPUT
C DAY                   - DAY STAMP
C DELAY(3)              - DELAY UNTIL FAULT INITIATION (MIN)
C DERROR(3)             - VECTOR OF DIFFERENCES BETWEEN ERROR AT
C                         SUCCESSIVE TIME STEPS
C DEXT(3)               - FAULT EXTENT DIFFERENCE AT CURRENT TIME STEP (MAX 3 FAULTS)
C DEXT0(3)              - DIFFERENCE BETWEEN FINAL EXTENT AND INITIAL (MAX 3 FAULTS)
C                         EXTENT OF VARIABLE AFFECTED BY FAULT
C DSEED                 - SEED FOR GAUSSIAN RANDOM NUMBER GENERATOR
C DT                    - TIME INCREMENT OF MAIN LOOP (MIN)
C ERROR(3)              - VECTOR OF CONTROLLER ERRORS (MAX 3 FAULTS)
C EXTENT0(3)            - ULTIMATE FAULT EXTENT (CONTEXT DEPENDENT!!) (MAX 3 FAULTS)
C EPD                   - EFFLUENT PRESSURE DROP CONSTRAINT RESIDUAL (M3/MIN)
C F(3)                  - FAULT CODES (CAN HAVE UP TO 3 SIMULTANEOUS FAULTS)
C FF                    - ADJUSTABLE PARAMETER
C FINTEG                - INTEGRAL OF FLOW RESIDUAL (M3)
C FLOW(8)               - VECTOR OF PROCESS FLOWRATES (M3/MIN)
C FMODE(2)              - CHARACTER ARRAY OF SENSOR FAILURE MODES
C FTYPE(19)             - CHARACTER ARRAY OF FAULT TYPES
C HOUR                  - HOUR STAMP
C ISEED                 - SEED FOR RANDOM NUMBER GENERATOR
C LEVEL                 - CSTR LEVEL (M)
C M                     - NUMBER OF FAULTS TO BE SIMULATED
C MASSBAL               - RESIDUAL OF INVENTORY CONSTRAINT (M3)
C MEAS(18,2)            - ARRAY OF PROCESS MEASUREMENTS; COL 2 IS
C                          MOST RECENT, COL1 IS EWMA HISTORY
C MINUTE                - MINUTE STAMP
C MOLBAL                - MOL BALANCE CONSTRAINT RESIDUAL (KMOL)
C MOLIN                 - CUMULATIVE MOL INFLUX (KMOL)
C MOLOUT                - CUMULATIVE MOL OUTFLUX (KMOL)
C MON                   - MONTH STAMP
C N                     - TOTAL NUMBER OF SENSORS AND CONSTRAINTS
C NORMVAL(18)           - ARRAY OF NORMAL STEADY STATE PROCESS MEASUREMENTS
C PB,PBCOMP             - PRESSURE AT TANK OUTLET (KG/M2)
C PCW                   - COOLING WATER SUPPLY PRESSURE (KG/M2)
C PP,PP0                - PUMP DIFFERENTIAL PRESSURE (KG/M2)
C QEXT                  - EXTERNAL HEAT SOURCE(SINK) TO CSTR (KJ/MIN)   (external heat loss due to fault)
C QJAC                  - HEAT TRANSFERED IN HEAT EXCHANGER (KJ/MIN)
C QREAC1                - PRIMARY HEAT OF REACTION (KJ/KMOL)
C QREAC2                - SECONDARY HEAT OF REACTION (KJ/KMOL)
C R(10)                 - VECTOR OF FLOW RESISTANCES (MIN KG^0.5/M4)
C RAND(3)               - RANDOM NUMBER VECTOR
C RCOMP(10)             - COMPUTED VALUES FOR FLOW RESISTANCE
C                         BASED ON SENSOR MEASUREMENT DATA
C REG1(3),REG2(3)       - STORAGE REGISTERS FOR PAST CONTROLLER ERRORS
C RHO,RHOLD             - DENSITY OF REACTOR CONTENTS (KG/M3)
C RHO1                  - DENSITY OF FEED (KG/M3)
C RHO2                  - DENSITY OF COOLING WATER (KG/M3)
C R0(10)                - NOMINAL VALUES FOR FLOW RESISTANCE (MIN KG^0.5/M4)
C RRATE1                - RATE OF PRIMARY REACTION IN CSTR (KMOL/M3 MIN)
C RRATE2                - RATE OF SECONDARY REACTION IN CSTR (KMOL/M3 MIN)
C RC                    - OVERALL RESISTANCE FOR CW STREAM (MIN KG^0.5/M4)
C RE                    - OVERALL RESISTANCE FOR EFFLUENT STREAM (MIN KG^0.5/M4)
C S0,S1,S2,S3,
C S4,S5,S6              - PROGRAM CONTROL PARAMETERS
C SDEV(18)              - STANDARD DEVIATION OF RANDOM SENSOR NOISE
C SEC                   - SECOND STAMP
C SENSORS(18)           - CHARACTER ARRAY OF SENSOR NAMES
C SP(3)                 - VECTOR OF CONTROLLER SETPOINTS
C T(4)                  - VECTOR OF PROCESS TEMPERATURES (C)
C TC(3)                 - EXPONENTIAL TIME CONSTANT FOR FAULT EXTENT
C TAREA                 - FLOOR AREA OF CSTR (M2)
C TEMP                  - A TEMPORARY STORAGE REGISTER
C TH                    - TIME HORIZON OF SIMULATION (MIN)
C THETA                 - EWMA (EXP. FILTER) PARMETER
C TIME                  - SIMULATION TIME (MIN)
C TUNE                  - CONTROLLER TUNING LOGICAL CONTROL
C TUNEI,TUNEJ           - CONTROLLER TUNING ARRAY ELEMENTS
C UA                    - HEAT EXCHANGER CONSTANT (KJ/MIN C)
C UNITS(18)             - CHARACTER ARRAY OF SENSOR UNI1S
C V(2)                  - VECTOR OF CONTROL VALVE POSITIONS [0,100]
C VOL,VOLD              - CSTR LIQUID VOLUME (M3)
C YEAR                  - YEAR STAMP
C ZCOUNT                - COUNTER FOR PRINTING OUTPUT
C ZLIM                  - PRINT INTERVAL (MIN)
C iter					- Counter for iteration
C
C           **BEGIN PROGRAM**
C
C DECLARE AND DIMENSION ALL VARIABLES
C
C REAL VARIABLES (ALL DOUBLE PRECISION)
C
      IMPLICIT NONE

      integer prot,dataout,nsmp
      parameter(prot=17,dataout=18)
      character *40 labelstr,csvsep*1
      real *8 sumConcABC,mol,dmolin,dmolout,measold



      REAL *8  ALPHA1,ALPHA2,BETA1,BETA2,B(3,3),CA0,CA,CB,CC,CCNOM,
     1         CNT(3),CP,CPOLD,CP1,CP2,CWPD,DELAY(3),DERROR(3),
     1         DEXT(3),DEXT0(3),DSEED,DT,EPD,
     1         ERROR(3),EXTENT0(3),FF,FINTEG,FLOW(8),JEP,
     1         LEVEL,MASSBAL,MEAS(18,2),MOLBAL,MOLIN,MOLOUT,
     1         NORMVAL(18),PB,PBCOMP,   aux,DFRAC, ! Additional func
     1         PCW,PP,PP0,QEXT,QJAC,QREAC1,QREAC2,R(10),RAND(3),
     1         RCOMP(10),REG1(3),REG2(3),REP,RHO,RHOLD,RHO1,RHO2,
     1         R0(10),RRATE1,RRATE2,RC,RE,RV(18),SP(3),SDEV(18),
     1         T(4),TAREA,TC(3),TEMP,TH,THETA,TIME,UA,V(2),VOL,VOLD,ZLIM
C
C INTEGER VARIABLES
C
      INTEGER I,J,K
      INTEGER *2 DAY,F(3),HOUR,MINUTE,M,MON,RO(18,2),SAMP(3,23),S0,
     1            S1(3),S2(3),S3,S4,S5,S6,S7,S8,SEC,SS,TUNE,TUNEI,
     1            TUNEJ
      INTEGER *4 N,YEAR,ZCOUNT,iter
      INTEGER *4 ISEED
C
C CHARACTER VARIABLES
C
      CHARACTER *8   UNITS(18)
      CHARACTER *12  FMODE(2)
      CHARACTER *40  SENSORS(18),SENSORACRONYMS(18),DATAFILE
      CHARACTER *40  FTYPE(23)
C
5     CONTINUE
C
C
C INITIALIZE CHARACTER DATA
C
C
C INITIALIZE FAULT NAMES
C
      DATA FTYPE /'NO FAULT                                ',
     1            'BLOCKAGE AT TANK OUTLET                 ', ! Fault 2
     1            'BLOCKAGE IN JACKET                      ', ! Fault 3
     1            'JACKET LEAK TO ENVIRONMENT              ', ! Fault 4
     1            'JACKET LEAK TO TANK                     ', ! Fault 5
     1            'LEAK FROM PUMP                          ', ! Fault 6
     1            'LOSS OF PUMP PRESSURE                   ', ! Fault 7
     1            'JACKET EXCHANGE SURFACE FOULING         ', ! Fault 8
     1            'EXTERNAL HEAT SOURCE (SINK)             ', ! Fault 9
     1            'PRIMARY REACTION ACTIVATION ENERGY      ', ! Fault 10
     1            'SECONDARY REACTION ACTIVATION ENERGY    ', ! Fault 11
     1            'ABNORMAL FEED FLOWRATE                  ', ! Fault 12, nominal value is 0.25, defined in initialization of FLOW(...)
     1            'ABNORMAL FEED TEMPERATURE               ', ! Fault 13
     1            'ABNORMAL FEED CONCENTRATION             ', ! Fault 14
     1            'ABNORMAL COOLING WATER TEMPERATURE      ', ! Fault 15
     1            'ABNORMAL COOLING WATER PRESSURE         ', ! Fault 16
     1            'ABNORMAL JACKET EFFLUENT PRESSURE       ', ! Fault 17
     1            'ABNORMAL REACTOR EFFLUENT PRESSURE      ', ! Fault 18
     1            'ABNORMAL LEVEL CONTROLLER SETPOINT      ', ! Fault 19
     1            'ABNORMAL TEMPERATURE CONTROLLER SETPOINT', ! Fault 20
     1            'CONTROL VALVE (CV-1) STUCK              ', ! Fault 21
     1            'CONTROL VALVE (CV-2) STUCK              ', ! Fault 22
     1            'SENSOR FAULT                            '/ ! Fault 23
C
C INITIALIZE SENSOR NAMES AND ACRONYMS
C
      DATA SENSORS /'FEED_CONCENTRATION_SENSOR       ',
     1              'FEED_FLOWRATE_SENSOR            ',
     1              'FEED_TEMPERATURE_SENSOR         ',
     1              'REACTOR_LEVEL_SENSOR            ',
     1              'CONCENTRATION_A_SENSOR          ',
     1              'CONCENTRATION_B_SENSOR          ',
     1              'REACTOR_TEMPERATURE_SENSOR      ',
     1              'COOLING_WATER_FLOWRATE_SENSOR   ',
     1              'PRODUCT_FLOWRATE_SENSOR         ',
     1              'COOLING_WATER_TEMPERATURE_SENSOR',
     1              'COOLING_WATER_PRESSURE_SENSOR   ',
     1              'LEVEL_CONTROLLER_OUTPUT_SIGNAL  ',
     1              'CW_FLOW_CONTROLLER_OUTPUT_SIGNAL',
     1              'COOLING_WATER_SETPOINT          ',
     1              'INVENTORY CONSTRAINT            ',
     1              'CW_PRESSURE_DROP_CONSTRAINT     ',
     1              'EFFL_PRESSURE_DROP_CONSTRAINT   ',
     1              'MOL_BALANCE_CONSTRAINT          '/
C
      DATA SENSORACRONYMS /'c_{A0}','Q_{1}','T_{1}','L','c_{A}',
     1                     'c_{B}','T_{2}','Q_{5}','Q_{4}','T_{3}',
     1                     '\mathrm{PCW}','\mathrm{CNT}_{1}',
     1                     '\mathrm{CNT}_{3}=SP_{Q_{5}}',
     1                     '\mathrm{CNT}_{2}','\mathrm{INV\ CON}',
     1                     '\mathrm{CW\ DP\ CON}',
     1                     '\mathrm{EFFL\ DP\ CON}',
     1                     '\mathrm{MOLBAL\ CON}'/

C
C
C
C INITIALIZE SENSOR UNITS
C
      DATA UNITS /'KMOL/M3 ',
     1            'M3/MIN  ',
     1            'C       ',
     1            'M       ',
     1            'KMOL/M3 ',
     1            'KMOL/M3 ',
     1            'C       ',
     1            'M3/MIN  ',
     1            'M3/MIN  ',
     1            'C       ',
     1            'KG/M2   ',
     1            '% OPEN  ',
     1            '% OPEN  ',
     1            'M3/MIN  ',
     1            'M3      ',
     1            'M3/MIN  ',
     1            'M3/MIN  ',
     1            'KMOL    '/
C
C INITIALIZE SENSOR FAULT MODE NAMES
C
      DATA FMODE /'FIXED BIAS ','FIXED VALUE '/
C
C INITIALIZE NUMERICAL DATA
C
      DATA CNT /74.7D0,0.9D0,59.3D0/      ! Controller outputs of three controllers
      DATA FLOW /0.25D0,0.25D0,0.00D0,0.25D0,0.9D0,0.00D0,0.00D0,0.9D0/
C Ph.D. Thesis Finch, p. 322 FLOW=Q
C	Q1=Feed flowrate (constant flow into the tank, nominal value 0.25), associated fault nr: 12 (abnormal feed flowrate)
C	Q2=flowrate between tank exit and division into leak flow Q3 and effluent flow Q4
C	Q3=leak flowrate, modeled by resistance R2 (high resistance=small leak)
C	Q4=Effluent flowrate, nominal value 0.25 (incoming into tank=outgoing product flowrate from system, without leaking)
C
C	Cooling system:
C	Q5=Cooling flowrate (constant flow into the tank, nominal value 0.9), calculated
C	Q6=Flowrate from cooling circuit through jacket into reactor, associated fault nr: 5 (jacket leak to tank), diminish resistance R7
C	Q7=Flowrate from cooling circuit into exterior area, associated fault nr: 4 (jacket leak to environment), diminish resistance R8
C	Q8=Effluent cooling circuit flowrate, nominal value 0.9 (incoming into jacket=outgoing from jacket, without leaking into tank or environment)
C
C
C
C
      DATA (SAMP(1,I),I=1,19) /1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,
     1  16,17,18,18/
      DATA (SAMP(2,I),I=1,19) /2,4,8,9,11,15,16,17,0,0,0,0,0,0,
     1                          0,0,0,0,8/
      DATA (SAMP(3,I),I=1,19) /3,7,10,12,13,14,0,0,0,0,0,0,0,0,0,0,
     1  0,0,6/
C
      DATA S1 /0, 0, 0/
      DATA S2 /0, 0, 0/
      DATA F /0, 0, 0/
C
      DATA SP  /2.00D0,80.0D0,0.9D0/  ! Set points of the controllers: REACTOR LEVEL, REACTOR TEMP, COOLING WATER FLOW RATE
      DATA SDEV /0.15D0,0.002D0,0.15D0,0.01D0,0.02D0,0.14D0,0.15D0,
     1           0.003D0,0.002D0, 0.15D0,400.0D0,0.01D0,0.01D0,
     1           0.0001D0,0.005D0,0.005D0,0.01D0,0.5D0/
      DATA ERROR /0.00D0,0.00D0,0.00D0/
      DATA DERROR /0.00D0,0.00D0,0.00D0/
C
      DATA (EXTENT0(I),I=1,3) /0.0D0,0.0D0,0.0D0/
      DATA (DELAY(I),I=1,3) /0.0D0,0.0D0,0.0D0/
      DATA (DEXT(I),I=1,3)  /0.0D0,0.0D0,0.0D0/
      DATA (DEXT0(I),I=1,3) /0.0D0,0.0D0,0.0D0/
      DATA (TC(I),I=1,3)  /0.0D0,0.0D0,0.0D0/
C
      DATA R      /100.0D0,1000000.0D0,0.00D0,500.0D0,72.0D0,0.00D0,
     1             1000000.0D0,1000000.0D0,0.00D0,65.0D0/
      DATA R0     /100.0D0,1000000.0D0,0.00D0,500.0D0,72.0D0,0.00D0,
     1             1000000.0D0,1000000.0D0,0.00,65.0D0/
      DATA REG1 /0.00D0,0.00D0,0.00D0/
      DATA REG2 /0.00D0,0.00D0,0.00D0/
      DATA T     /30.0D0,80.0D0,20.0D0,40.0D0/
      DATA V     /74.7D0,59.3D0/
C
      DATA (B(I,1),I=1,3) /35.0D0,-0.040D0,-25.0D0/
      DATA (B(I,2),I=1,3) /5.00D0,-0.020D0,-75.0D0/
      DATA (B(I,3),I=1,3) /0.00D0,0.00D0,0.00D0/
C
C <<<<<<<<<<<< CORRECTION: Normal CNT2 41.7 is wrong, since in first simlation loop changes to 40.7
      DATA (MEAS(I,1),I=1,18)  /20.0D0,0.25D0,30.0D0,2.00D0,2.85D0,
     1                          17.11D0,80.0D0,0.9D0,0.25D0,20.0D0,
     1                          56250.0D0,25.3D0,40.7D0,0.9D0,0.0D0,
     1                          0.0D0,0.0D0,0.0D0/
C
      DATA (MEAS(I,2),I=1,18)  /20.0D0,0.25D0,30.0D0,2.00D0,2.85D0,
     1                          17.11D0,80.0D0,0.9D0,0.25D0,20.0D0,
     1                          56250.0D0,25.3D0,40.7D0,0.9D0,0.0D0,
     1                          0.0D0,0.0D0,0.0D0/
C
      DATA (NORMVAL(I),I=1,18)  /20.0D0,0.25D0,30.0D0,2.00D0,2.85D0,
     1                           17.11D0,80.0D0,0.9D0,0.25D0,20.0D0,
     1                           56250.0D0,25.3D0,40.7D0,0.9D0,0.0D0,
     1                           0.0D0,0.0D0,0.0D0/
C
      data LEVEL /2.0D0/  ! CSTR level at t=0
      data MASSBAL /0.0D0/
      data EPD /0.0D0/
      data CWPD /0.0D0/
      data MOLBAL /0.0D0/

      open(unit=prot,file='./output/log.txt',status='replace')
      write(prot,*) 'Opening protocol file ''./output/log.txt''...'
      write(*,*) 'Opening protocol file ''log.txt''...'
      open(unit=dataout,file='./output/X.csv',status='replace')
      write(*,*) 'Opening output data file ''./output/X.csv''...'
      csvsep = char(9)      ! char(9)=TAB
      csvsep = ';'



C     TOOLDIAG OUTPUT FORMAT: First line = number of features
C     Number of features
C      write(dataout,'(A2)') '18'

C     Python Pandas (or CSV)
      write(dataout, '(18(A30,A1),A10)')
     1              (SENSORACRONYMS(I),csvsep,I=1,18),'CLASS'
C
C INITIALIZE SCALAR VARIABLES
C
      ALPHA1      =     2500.0D0
      BETA1       =     25000.0D0
      ALPHA2      =     3000.0D0
      BETA2       =     45000.0D0
      CA0         =     20.0D0
      CA          =     2.850D0  ! CA+CB+CC = 19.9866 --- Shouldn't it be 20.0?
      CB          =     17.114D0
      CCNOM       =     0.0226D0
      CC          =     CCNOM
      CP          =     4.20D0
      CP1         =     4.20D0
      CP2         =     4.20D0
      DT          =     0.02D0
      FF          =     0.10D0
      FINTEG      =     0.0D0
      JEP         =     0.0D0
      MOLIN       =     0.0D0
      MOLOUT      =     0.0D0
      N           =     18
      PB          =     2000.0D0
      PCW         =     56250.0D0
      PP          =     48000.0D0
      PP0         =     48000.0D0
      QREAC1      =     30000.0D0
      QREAC2      =     -10000.0D0
      QEXT        =     0.0D0
      REP         =     0.0D0
      RHO         =     1000.0D0
      RHO1        =     1000.0D0
      RHO2        =     1000.0D0
      S0          =     0
      S3          =     0
      S4          =     0
      S7          =     0
      S8          =     0
      TAREA       =     1.5D0
      TIME        =     0.0D0
      UA          =     1901.0D0
      VOL         =     3.0D0
C
C PLOT QUERY
C
      WRITE (*,800)           ! *** JACKETED CSTR DYNAMIC SIMULATION ***
      WRITE (*,735)           ! PLOT RESULTS OF PREVIOUS RUN? [0=NO , 1=YES]
      READ (*,*) S6
C      IF (S6 .EQ. 1) GOTO 550
      IF (S6 .EQ. 1) THEN
          CALL SIMPLOT(DATAFILE,N,NORMVAL,MEAS,SDEV,SENSORS,UNITS)
          STOP
      ENDIF
C
C OPEN OUTPUT FILE FOR RAW DATA DUMP (MIDAS FORMAT)
C
      WRITE (*,700)           ! ENTER OUTPUT FILE NAME [DEFAULT=DUMP.RAW]
      READ (*,710) DATAFILE
C      IF (DATAFILE .EQ. '') THEN
      IF (DATAFILE(1:1) .EQ. ' ') THEN
        DATAFILE='DUMP.RAW'
      ENDIF
C
      OPEN (UNIT=14,FILE='./output/'//DATAFILE)
C
C ENTER DATE STAMP DATA
C
C  WRITE (*,720)              ! ENTER DATE STAMP [MM/DD/YYYY]
C  READ (*,730) MON,DAY,YEAR
      MON=07
      DAY=3
      YEAR=2024
C
C CONTROLLER TUNING (OPTIONAL)
C
600   CONTINUE
      WRITE (*,740)           ! OVERRIDE CONTROLLER TUNING ? [0=NO , 1=YES]
      READ (*,*) TUNE
      IF (TUNE .EQ. 0) GOTO 650
610   WRITE (*,745)           ! INDICATE CONTROLLER: 1=LEVEL, 2=TEMP, 3=COOLFLOW
      READ (*,*) TUNEI
C      IF ((TUNEI .GT. 4) .OR. (TUNEI .LT. 1)) GOTO 610
      IF ((TUNEI .GT. 3) .OR. (TUNEI .LT. 1)) GOTO 610
620   WRITE (*,750)           ! ENTER PARAMETER TO CHANGE: 1=GAIN 2=INTEGAL 3=DERIVATIVE
      READ (*,*) TUNEJ
      IF ((TUNEJ .GT. 3) .OR. (TUNEJ .LT. 1)) GOTO 620
C      write(*,*) 'TUNEI=', TUNEI, 'TUNEJ=', TUNEJ, 'B=', B
      WRITE (*,760) B(TUNEI,TUNEJ)  ! ENTER NEW VALUE [CURRENT VALUE IS
      READ (*,*) B(TUNEI,TUNEJ)
      GOTO 600
650   CONTINUE
C
C ENTER EXPONENTIAL FILTER CONSTANT (1=NO EWMA)
C
660   WRITE (*,770)           ! ENTER FILTER CONSTANT (0.0 - 1.0]
      READ (*,*) THETA
      IF ((THETA .LE. 0.0) .OR. (THETA .GT. 1.0)) GOTO 660
C
C ENTER PRINT FORMAT
C
      WRITE (*,772)           ! PRINT OUTPUT IN RANDOM ORDER? [0=NO , 1=YES]
      READ (*,*) S7
      WRITE (*,773)           ! VARIABLE OUTPUT FREQUENCY? [0=NO , 1=YES]
      READ (*,*) S8
C
C COMPUTE SEED FOR RANDOM NUMBER GENERATOR
C
      WRITE (*,790)           ! ENTER INTEGER SEED FOR RANDOM NUMBER GENERATOR
      READ (*,*) ISEED
      DSEED=DBLE(ISEED)
C      write(*,*) 'ISEED=',ISEED,'DSEED=',DSEED
C      stop
C
C ENTER FAULT INFORMATION
C
20    CONTINUE
      DO 30 I=1,23
        WRITE(*,810) I,FTYPE(I)
30    CONTINUE
C
C35    S0=S0 + 1         ! counter of the fault (one or more)
      S0=S0 + 1         ! counter of the fault (one or more)
      WRITE (*,820)                 ! ENTER DESIRED FAULT TYPE
      READ (*,*) F(S0)
      DO 37 I=1,(S0-1)
      IF (F(I) .EQ. F(S0)) THEN
        WRITE (*,822)               ! INDICATED FAULT HAS ALREADY BEEN SELECTED!
        WRITE (*,823)               ! CONTINUE? [0=NO , 1=YES)
        READ (*,*) S6
        IF (S6 .EQ. 1) GOTO 37
        S0=S0 - 1
        GOTO 40
      ENDIF
37    CONTINUE
      IF (F(S0) .EQ. 1) GOTO 40
      IF (F(S0) .EQ. 23) THEN
        WRITE (*,821) (I,SENSORS(I),I=1,14)  ! PROCESS SENSORS ARE : 
        WRITE (*,825)               ! ENTER NUMBER OF FAILED SENSOR
        READ (*,*) S1(S0)
        WRITE (*,830)         ! ENTER TYPE OF SENSOR FAILURE: O=FIXED BIAS 1=FIXED VALUE
        READ (*,*) S2(S0)
        WRITE (*,*) 'NOMINAL VALUE IS : ' ,MEAS(S1(S0),1)
      ELSE
      IF (F(S0) .EQ. 2) WRITE (*,835) R(1) ! NOMINAL VALUE IS
      IF (F(S0) .EQ. 3) WRITE (*,835) R(9) ! 'BLOCKAGE IN JACKET', ! Fault 3
      IF (F(S0) .EQ. 4) WRITE (*,835) R(8)
      IF (F(S0) .EQ. 5) WRITE (*,835) R(7)
      IF (F(S0) .EQ. 6) WRITE (*,835) R(2)
      IF (F(S0) .EQ. 7) WRITE (*,835) PP
      IF (F(S0) .EQ. 8) WRITE (*,835) UA
      IF (F(S0) .EQ. 9) WRITE (*,835) QEXT
      IF (F(S0) .EQ. 10) WRITE (*,835) BETA1
      IF (F(S0) .EQ. 11) WRITE (*,835) BETA2
      IF (F(S0) .EQ. 12) WRITE (*,835) FLOW(1)
      IF (F(S0) .EQ. 13) WRITE (*,835) T(1)
      IF (F(S0) .EQ. 14) WRITE (*,835) CA0
      IF (F(S0) .EQ. 15) WRITE (*,835) T(3)
      IF (F(S0) .EQ. 16) WRITE (*,835) PCW
      IF (F(S0) .EQ. 17) WRITE (*,835) JEP
      IF (F(S0) .EQ. 18) WRITE (*,835) REP
      IF (F(S0) .EQ. 19) WRITE (*,835) SP(1)
      IF (F(S0) .EQ. 20) WRITE (*,835) SP(2)
      IF (F(S0) .EQ. 21) WRITE (*,835) 100.0D0 - V(1)
      IF (F(S0) .EQ. 22) WRITE (*,835) 100.0D0 - V(2)
      ENDIF
      WRITE (*,840)           ! ENTER EXTENT OF FAULT (IN APPROPRIATE UNITS)
      READ (*,*) EXTENT0(S0)
      WRITE (*,850)           ! ENTER FAULT DELAY [min]
      READ (*,*) DELAY(S0)
C     IF (F(S0) .GE. 20) GOTO 40
      WRITE (*,860)           ! ENTER TIME CONSTANT [(min)-1]
      READ (*,*) TC(S0)
40    IF (S0 .GE. 3) GOTO 42
      WRITE (*,865)           ! ADD ANOTHER FAULT? [0=NO , 1=YES]
      READ (*,*) S6
      IF (S6 .EQ. 1) GOTO 20
42    M=S0  ! if condition is only normal, M will have value 1
      WRITE (*,870)           ! 'ENTER TIME HORIZON FOR SIMULATION [min]
      READ (*,*) TH
      WRITE (*,873)           ! 'PRINT INTERMEDIATE RESULTS? [1=YES/0=NO]
      READ (*,*) ZLIM
      IF (ZLIM .EQ. 0.0) THEN
        ZLIM=TH
        GOTO 45
      ENDIF
      WRITE (*,875)           ! 'ENTER PRINT INTERVAL [min]
      READ (*,*) ZLIM
45    ZLIM=ZLIM/DT
      ZCOUNT=IDNINT(ZLIM)
      DO 47 S0=1,M
      IF (F(S0) .EQ. 23) THEN
        WRITE (*,880) FTYPE(23),SENSORS(S1(S0)),FMODE(S2(S0)+1),
C
C     PROCESS CONDITION IS -- SENSOR FAULT
C     IN (which sensor)
C     SENSOR FAULT IS TYPE -- (fixed or bias) 
C         FAULT EXTENT IS ...
C         FAULT DELAY IS ... (MIN)
C         FILTER CONSTANT IS  --
C
     1  EXTENT0(S0),DELAY(S0)
      ELSE
        WRITE (*,890) FTYPE(F(S0)),EXTENT0(S0),DELAY(S0),TC(S0)
      ENDIF
C     write(*,*) 'F(S0)=',F(S0),'S1(S0)=',S1(S0),'S2(S0)=',S2(S0)
47    CONTINUE

      WRITE (*,895) THETA
      WRITE (*,950)   !  RUNNING.....

C     write(*,*) 'S0=', S0, 'M=', M
C     stop
C
C
C
C======================================================================
C TOP OF MAIN ITERATION LOOP
C======================================================================
C
      iter=0   ! Number of 50 Hz iterations
      nsmp = 0 ! Number of data samples

      CALL CLASSLABEL(S1, S2, F, M, TIME, DELAY, labelstr)

      CALL MEASOUT( dataout, MEAS, csvsep, labelstr, iter, nsmp )
C      if(.TRUE.) THEN
      if(.FALSE.) THEN
        write(*,*) 'INITIAL VALUES'
        write(*, '(A11,18(A15,A1))')
     1         '',(SENSORACRONYMS(I),csvsep,I=1,18)
        CALL PEEPVARS(MEAS, labelstr, iter)
C      stop
      END IF

50    CONTINUE
C
C INTRODUCE FAULT (OR NOT)


C
      DO 60 S0=1,M  ! M = NUMBER OF FAULTS TO BE SIMULATED

      IF ((F(S0) .EQ. 1).OR.(TIME .LT. DELAY(S0))) GOTO 60 ! Normal or delay not yet reached
C
C CALCULATE EXTENT DIFFERENTIAL , S0 is fault number, cannot be 1
C
      IF (DEXT0(S0) .EQ. 0.0) THEN
        IF (F(S0) .EQ. 2) DEXT0(S0)=EXTENT0(S0) - R(1) ! Fault BLOCKAGE AT TANK OUTLET, R(1) nominal=100
        IF (F(S0) .EQ. 3) DEXT0(S0)=EXTENT0(S0) - R(9)
        IF (F(S0) .EQ. 4) DEXT0(S0)=EXTENT0(S0) - R(8)
        IF (F(S0) .EQ. 5) DEXT0(S0)=EXTENT0(S0) - R(7)
        IF (F(S0) .EQ. 6) DEXT0(S0)=EXTENT0(S0) - R(2)
        IF (F(S0) .EQ. 7) DEXT0(S0)=EXTENT0(S0) - PP
        IF (F(S0) .EQ. 8) DEXT0(S0)=EXTENT0(S0) - UA ! Fault JACKET EXCHANGE SURFACE FOULING
        IF (F(S0) .EQ. 9) DEXT0(S0)=EXTENT0(S0) - QEXT
        IF (F(S0) .EQ. 10) DEXT0(S0)=EXTENT0(S0) - BETA1
        IF (F(S0) .EQ. 11) DEXT0(S0)=EXTENT0(S0) - BETA2
        IF (F(S0) .EQ. 12) DEXT0(S0)=EXTENT0(S0) - FLOW(1)
        IF (F(S0) .EQ. 13) DEXT0(S0)=EXTENT0(S0) - T(1)
        IF (F(S0) .EQ. 14) DEXT0(S0)=EXTENT0(S0) - CA0
        IF (F(S0) .EQ. 15) DEXT0(S0)=EXTENT0(S0) - T(3)
        IF (F(S0) .EQ. 16) DEXT0(S0)=EXTENT0(S0) - PCW
        IF (F(S0) .EQ. 17) DEXT0(S0)=EXTENT0(S0) - JEP
        IF (F(S0) .EQ. 18) DEXT0(S0)=EXTENT0(S0) - REP
        IF (F(S0) .EQ. 19) DEXT0(S0)=EXTENT0(S0) - SP(1) ! Fault 'ABNORMAL LEVEL CONTROLLER SETPOINT'
        IF (F(S0) .EQ. 20) DEXT0(S0)=EXTENT0(S0) - SP(2)
      ENDIF
C
C CALCULATE FAULT EXTENT EXPONENTIAL GROWTH FUNCTION
C
      DEXT(S0)=EXP(TC(S0)*(DELAY(S0)-TIME))
C     write(*,*)'S0=',S0,'F(S0)=',F(S0),'DEXT(S0)=',DEXT(S0),
C    1          'TC(S0)=',TC(S0),'DELAY(S0)=',DELAY(S0),'TIME=',TIME
C
C CALCULATE FAULT EXTENT
C
      IF (F(S0) .EQ. 2) R(1)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0) ! Fault BLOCKAGE AT TANK OUTLET
C      write(*,*) 'R(1)=',R(1)
C      stop
      IF (F(S0) .EQ. 3) R(9)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 4) R(8)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0) ! Fault 'JACKET LEAK TO ENVIRONMENT', diminish resistance R(8)
      IF (F(S0) .EQ. 5) R(7)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0) ! Fault 'JACKET LEAK TO TANK', diminish resistance R(7)
      IF (F(S0) .EQ. 6) R(2)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 7) PP=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 8) UA=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 9) QEXT=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 10) BETA1=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 11) BETA2=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 12) FLOW(1)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 13) T(1)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 14) CA0=EXTENT0(S0) - DEXT0(S0)*DEXT(S0) ! ABNORMAL FEED CONCENTRATION
      IF (F(S0) .EQ. 15) T(3)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 16) PCW=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 17) JEP=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 18) REP=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 19) SP(1)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
      IF (F(S0) .EQ. 20) SP(2)=EXTENT0(S0) - DEXT0(S0)*DEXT(S0)
C
60    CONTINUE
C
C CALCULATE CONTROLLER OUTPUTS
C https://en.wikipedia.org/wiki/PID_controller
C https://controls.engin.umich.edu/wiki/index.php/CascadeControl
C
C      DATA SP  /2.00,80.0,0.9/  ! Set points of the controllers: REACTOR LEVEL, REACTOR TEMP, COOLING WATER FLOW RATE
C
C REG1(3),REG2(3)       - STORAGE REGISTERS FOR PAST CONTROLLER ERRORS
C ERROR(3)              - VECTOR OF CONTROLLER ERRORS
C DERROR(3)             - VECTOR OF DIFFERENCES BETWEEN ERROR AT SUCCESSIVE TIME STEPS
C B(3,3)                - ARRAY OF CONTROLLER CONSTANTS FOR PID CONTROL K_P=(1,1), K_I=(1,2), K_D=(1,3)
C                         (IF B(1,2)=0 & B(1,3)=0 THEN P CONTROL, IF B(1,3)=0 THEN PI CONTROL)
C      DATA (B(I,1),I=1,3) /35.0,-0.040,-25.0/  K_P for Level controller, temperature controller, flow controller
C      DATA (B(I,2),I=1,3) /5.00,-0.020,-75.0/  K_I for Level controller, temperature controller, flow controller
C      DATA (B(I,3),I=1,3) /0.00,0.00,0.00/     K_D for Level controller, temperature controller, flow controller
C
C LEVEL CONTROLLER
C
      REG2(1)=REG1(1)
      REG1(1)=DERROR(1)
      DERROR(1)=ERROR(1) - (SP(1) - MEAS(4,1))  ! SP(1)=setpoint Level; MEAS(4,1)=LEVEL
      ERROR(1)=ERROR(1) - DERROR(1)
      CNT(1)=MIN(100.0D0,MAX(0.0D0,CNT(1)-B(1,1)*DERROR(1)
     1 +B(1,2)*(0.5D0*DERROR(1)+ERROR(1))*DT
     2 +B(1,3)*(2.0D0*REG1(1)-0.5D0*REG2(1)-1.5D0*DERROR(1))/DT))
C      write(*,*) 'iter=',iter,'dif=',SP(1) - MEAS(4,1),
C     1            'DERROR=',DERROR(1),'REG2(1)=',REG2(1),
C     1            'REG1(1)=',REG1(1),'ERROR(1)=',ERROR(1)
C      if (iter.eq.4) stop
C
C REACTOR TEMPERATURE CONTROLLER
C
      REG2(2)=REG1(2)
      REG1(2)=DERROR(2)
      DERROR(2)=ERROR(2) - (SP(2) - MEAS(7,1))  ! SP(2)=setpoint CSTR temperature; MEAS(7,1)=CSTR temperature
      ERROR(2)=ERROR(2) - DERROR(2)
      CNT(2)=MAX(0.0D0,CNT(2)-B(2,1)*DERROR(2)
     1 +B(2,2)*(0.5D0*DERROR(2)+ERROR(2))*DT
     2 +B(2,3)*(2.0D0*REG1(2)-0.5D0*REG2(2)-1.5D0*DERROR(2))/DT)
C
C COOLING WATER FLOW CONTROLLER
C
      REG2(3)=REG1(3)
      REG1(3)=DERROR(3)
      DERROR(3)=ERROR(3) - (SP(3) - MEAS(8,1))  ! SP(3)=cooling water flowrate; MEAS(8,1)=inflow cooling water (FLOW(5))
      ERROR(3)=ERROR(3) - DERROR(3)
      CNT(3)=MIN(100.0D0,MAX(0.0D0,CNT(3)-B(3,1)*DERROR(3)
     1 +B(3,2)*(0.5D0*DERROR(3)+ERROR(3))*DT
     2 +B(3,3)*(2.0D0*REG1(3)-0.5D0*REG2(3)-1.5D0*DERROR(3))/DT))

C      do i=1,3
C        write(*,*) REG2(i), REG1(i), DERROR(i), ERROR(i), CNT(i)
C      enddo
C      stop
C
C EVALUATE SAFETY SYSTEMS
C
C Attention: Only a message is issued and a flag is set, simulation continues util reaching timeout, no simulation abort
C
C
      IF ((MEAS(4,1) .GE. 2.75).OR.(MEAS(7,1) .GE. 130.0)) THEN ! MEAS(4,1)=LEVEL, MEAS(7,1)=T(2)=Reactor temperature
        IF (S3 .EQ. 1) GOTO 150  ! If inflow FLOW(1) is already closed, continue normally
        DO 110 S0=1,M
        IF (F(S0) .EQ. 12) THEN  ! Eliminate Fault 12, ABNORMAL FEED FLOWRATE, because inflow has been set to zero
             EXTENT0(S0)=0.0D0
             DEXT0(S0)=0.0D0
        ENDIF
110     CONTINUE
        FLOW(1)=0.0D0          ! Close inflow FLOW(1)
        S3=1                   ! Set closed inflow flag true
        WRITE (*,970) TIME     ! '***** EMERGENCY SHUTDOWN INITIATED AT '
        WRITE (prot,970) TIME
      ELSEIF (MEAS(4,1) .LE. 1.2) THEN ! MEAS(4,1)=LEVEL, reactor below level
        IF (S4 .EQ. 1) GOTO 150 ! If pump is already shut down, continue normally
        DO 120 S0=1,M
        IF (F(S0) .EQ. 7) THEN ! Eliminate Fault 7, LOSS OF PUMP PRESSURE, because pump head has been set to zero
             EXTENT0(S0)=0.0D0
             DEXT0(S0)=0.0D0
        ENDIF
120     CONTINUE
        PP=0.0D0               ! Shut down pump, no more pump head
        S4=1                   ! Set shut down pump flag true
        WRITE (*,980) TIME     ! '***** LOW LEVEL FORCES PUMP SHUTDOWN AT '
        WRITE (prot,980) TIME
      ENDIF
C
C CALCULATE VALVE POSITIONS
C
150   DO 160 S0=1,M
        IF ((F(S0) .EQ. 21) .AND. (TIME .GE. DELAY(S0))) THEN
          V(1)=100.0D0 - EXTENT0(S0)*(1.0D0-DEXT(S0))
        ELSE
          V(1)=MIN(100.0D0,MAX(0.0D0,100.0D0-MEAS(12,2))) ! MEAS(12,2)=CNT(1)
C          write(*,*) 'V(1)=',V(1)
        ENDIF
        IF ((F(S0) .EQ. 22) .AND. (TIME .GE. DELAY(S0))) THEN
          V(2)=100.0D0 - EXTENT0(S0)*(1.0D0-DEXT(S0))
        ELSE
          V(2)=MIN(100.0D0,MAX(0.0D0,100.0D0-MEAS(13,2))) ! MEAS(13,2)=CNT(3)
C          write(*,*) 'V(2)=',V(2)
        ENDIF
C        write(*,*) 'Fault ',S0,'V(1)=',V(1),'V(2)=',V(2)
160   CONTINUE
C
C CALCULATE FLOWRATES
C
C Valve resistances from the valve positions obtained from valve controllers
C
      R(3)=5.0D0 * EXP(0.0545D0*V(1))
      R(6)=5.0D0 * EXP(0.0545D0*V(2))
C
C Combined pipe, valve and leak resistances: Finch p. 319, p. 322
C RE = R_effluent
C RC = R_coolant
      RE=(1/((1/R(2))+(1/(R(3)+R(4)))))+R(1)
      RC=(1/((1/R(7))+(1/R(8))+(1/(R(9)+R(10)))))+R(5)+R(6)
C
C NOTE: THESE FORMULAS WILL NOT WORK WELL IF BOTH AN
C ABNORMAL EFFLUENT PRESSURE AND A LEAK ARE SIMULATED
C
C JEP = Jacket effluent pressure (fault 17, nominal zero)
C REP = Reactor effluent pressure (fault 18, nominal zero)
C PCW = COOLING WATER SUPPLY PRESSURE, nominal 56250.0 KG/M^2 
C Simulates probably some obstruction
C
      aux=(PP+PB-REP)
      if (aux .le. 0) THEN
        FLOW(2)=0
      else
        FLOW(2)=(1D0/RE)*aux**0.5D0
      endif
      aux = ((PP+PB-REP)-(FLOW(2)*R(1))**2.0D0)
      if (aux .le. 0) THEN
        FLOW(3)=0
      else
        FLOW(3)=(1D0/R(2))*aux**0.5D0
      endif
      FLOW(4)=FLOW(2)-FLOW(3)
C
      FLOW(5)=(1D0/RC)*(PCW-JEP)**0.5D0
      aux=(PCW-JEP-(FLOW(5)*(R(5)+R(6)))**2.0D0)
      if (aux .le. 0) THEN
        FLOW(6)=0
        FLOW(7)=0
      else
        FLOW(6)=(1D0/R(7))*aux**0.5D0
        FLOW(7)=(1D0/R(8))*aux**0.5D0
      endif
      FLOW(8)=FLOW(5)-FLOW(6)-FLOW(7)
C      if (.TRUE.) then
      if (.FALSE.) then
C      write(prot,
      write(*,
C     1    '(1x,A,8(A,e12.5),A,3(f7.2),2(A,f6.2,A),4(A,f12.2))')
     1    '(1x,A,e25.15,A,8(A,e15.5),A,3(e15.5),
     1      2(A,f9.4,A),4(A,e15.5))')
     1    'aux=',aux,'  FLOWRATES Q=',
     1    ' Q1=',FLOW(1),' Q2=',FLOW(2),' Q3=',FLOW(3),' Q4=',FLOW(4),
     1    ' Q5=',FLOW(5),' Q6=',FLOW(6),' Q7=',FLOW(7),' Q8=',FLOW(8),
     1    '  CONTROLLERS=',CNT(1),CNT(2),CNT(3),
     1    '  VALVE POS: V(1)=',V(1),'%',' V(2)=',V(2),'%',
     1    '  R: R(3)=',R(3),' R(6)=',R(6),' RE=',RE,' RC=',RC
C      stop
      endif

C
C CALCULATE JACKET TEMPERATURE
C
      T(4)=(UA*T(2)+RHO2*CP2*FLOW(8)*T(3)) / (UA+RHO2*CP2*FLOW(8))
C
C CALCULATE HEAT FLUX
C
      QJAC=UA*(T(2)-T(4))
C
C CALCULATE CSTR VARIABLES
C
C LEVEL/VOLUME/DENSITY/HEAT CAPACITY
C
C220   VOLD=VOL
      VOLD=VOL
      RHOLD=RHO
      CPOLD=CP
C
      VOL=VOLD+(FLOW(1)+FLOW(6)-FLOW(2))*DT
      RHO=(1/VOL)*(VOLD*RHO)+
     1 (1/VOL)*(DT*(FLOW(1)*RHO1+FLOW(6)*RHO2-FLOW(2)*RHO))
C
      CP=(1/VOL)*(VOLD*CP)+
     1 (1/VOL)*(DT*(FLOW(1)*CP1+FLOW(6)*CP2-FLOW(2)*CP))
C
      LEVEL=VOL/TAREA
C      if (.TRUE.) then
      if (.FALSE.) then
      write(prot,'(1x,6(A,e12.5))') 'Level=', LEVEL, ' RHO=',RHO,
     1    ' RHO1=',RHO1, ' RHO2=',RHO2, ' CP=',CP, ' VOL=',VOL
      endif

      IF (LEVEL .LE. 0.0) THEN
        WRITE (*,*) 'FAILURE DUE TO LOW LEVEL'
        STOP
      ENDIF
      PB=RHO*LEVEL
C     write(*,*) 'k=',iter,'FLOW(1)=',FLOW(1),'FLOW(6)=',FLOW(6),
C    1           'FLOW(2)=',FLOW(2),'LEVEL=',LEVEL,'VOL=',VOL
C
C CONCENTRATIONS p. 235, https://en.wikipedia.org/wiki/Arrhenius_equation
C
C http://www.nyu.edu/classes/tuckerman/pchem/lectures/lecture_21.pdf
C
C
      RRATE1=ALPHA1*EXP(-BETA1/(8.314D0*(273.15D0+T(2))))*CA
      RRATE2=ALPHA2*EXP(-BETA2/(8.314D0*(273.15D0+T(2))))*CA
C
      CA=(1/VOL)*(VOLD*CA)+
     1 (1/VOL)*(FLOW(1)*CA0-FLOW(2)*CA-RRATE1*VOLD-RRATE2*VOLD)*DT
C
      CB=(1/VOL)*(VOLD*CB)+
     1 (1/VOL)*(RRATE1*VOLD-FLOW(2)*CB)*DT
C
      CC=(1/VOL)*(VOLD*CC)+
     1 (1/VOL)*(RRATE2*VOLD-FLOW(2)*CC)*DT
C
C TEMPERATURE
C
      T(2)=(1/(VOL*RHO*CP))*(VOLD*RHOLD*CPOLD*T(2))+
     1 (1/(VOL*RHO*CP))*(((QREAC1*RRATE1+QREAC2*RRATE2)*VOLD)*DT)+
     1 (1/(VOL*RHO*CP))*((QEXT-QJAC)*DT)+
     1 (1/(VOL*RHO*CP))*(FLOW(1)*RHO1*CP1*T(1)*DT)+
     1 (1/(VOL*RHO*CP))*(FLOW(6)*RHO2*CP2*T(4)*DT)-
     1 (1/(VOL*RHO*CP))*(FLOW(2)*RHOLD*CPOLD*T(2)*DT)

C      if (.TRUE.) then
      if (.FALSE.) then
      write(*,*) 'T(4)=',T(4),'QJAC=',QJAC
      write(*,*) 'VOL=',VOL,'RHO=',RHO
      write(*,*) 'T(2)=',T(2),'RRATE1=',RRATE1,'RRATE2=',RRATE2
      write(*,*) 'CA=',CA,'CB=',CB,'CC=',CC
      endif
C
C MEASURE SENSED VARIABLES
C
      MEAS(1,2)=CA0
      MEAS(2,2)=FLOW(1)
      MEAS(3,2)=T(1)
      MEAS(4,2)=LEVEL
      MEAS(5,2)=CA
      MEAS(6,2)=CB
      MEAS(7,2)=T(2)
      MEAS(8,2)=FLOW(5)
      MEAS(9,2)=FLOW(4)
      MEAS(10,2)=T(3)
      MEAS(11,2)=PCW ! COOLING WATER SUPPLY PRESSURE, nominal 56250.0 KG/M^2  
      MEAS(12,2)=100.0D0 - CNT(1)
      MEAS(13,2)=100.0D0 - CNT(3)
      MEAS(14,2)=CNT(2)
C
C MODIFY SENSOR READINGS
C
C CALL GGNML(A,B,C)
C
C GGNML IS A GAUSIAN RANDOM NUMBER GENERATOR PRODUCING A VECTOR C
C OF NORMAL (0,1) RANDOM NUMBERS OF DIMENSION B. THE SEED (A) MUST
C BE DOUBLE PRECISION.
C
C     write(*,*)'S0:',S0,'S1:',(S1(I),I=1,3),
C    1          'S2:',(S2(I),I=1,3),'I=',I,'S1(S0)=',S1(S0)
      measold = MEAS(1,2)
      DO I=1,14
C       measold = MEAS(I,2)
C       write(*,*) 'BEFORE: Sensor=', I, 'seed=', DSEED
        CALL GGNML(DSEED,3,RAND)
C       write(*,*) 'AFTER: Sensor=', I, 'seed=', DSEED, 'RAND=', RAND
C       stop
        DO S0=1,M
C         SENSOR FAULTS
          aux=RAND(S0)*SDEV(I)
          IF ((F(S0) .EQ. 23) .AND. (TIME .GE. DELAY(S0))) THEN ! IF SENSOR FAULT ACTIVE
            IF ((I .EQ. S1(S0)) .AND. (S2(S0) .EQ. 0)) THEN ! 0=BIAS
              MEAS(I,2)=MEAS(I,2)+RAND(S0)*SDEV(I)+
     1             EXTENT0(S0)*(1.0D0-DEXT(S0))
            ELSE IF ((I .EQ. S1(S0)) .AND. (S2(S0) .EQ. 1)) THEN ! 1=VALUE
              MEAS(I,2)=EXTENT0(S0)*(1.0D0-DEXT(S0))
            ELSE
              MEAS(I,2)=MEAS(I,2)+RAND(S0)*SDEV(I) ! <><><><><><><><>
            END IF
C           write(*,*) 'Sensor fault:','S0=',S0,'F(S0)=',F(S0),'I=',I,
C    1                'S1(S0)=',S1(S0),'S2(S0)=',S2(S0),'measold=',
C    1                 measold,'MEAS(I,2)=',MEAS(I,2),'iter=',iter,
C    1                 'EXTENT0(S0)=',EXTENT0(S0),'DEXT(S0)=',DEXT(S0),
C    1                 'DEXT0(S0)=',DEXT0(S0),'aux=',RAND(S0)*SDEV(I)
C           if (i.eq.14) stop
C         NON SENSOR FAULTS
C           stop
          ELSE
C           write(*,*) 'S0=',S0,'F(S0)=',F(S0)
            MEAS(I,2)=MEAS(I,2)+RAND(S0)*SDEV(I)
          END IF
C         write(*,'(1x,A,I4,A,e12.5,A,e12.5,A)') 'NOISE:MEAS #',
C    1           I,' ',measold,' --> ',MEAS(I,2),' <==='
        END DO
      END DO
C     write(*,*) '>>>>iter=',iter,'aux=',aux,'measold=', measold,
C    1           'measCA0new=MEAS(1,2)=',MEAS(1,2),'TIME=',TIME
C
C CALCULATE EXPONENTIAL WEIGHTED MOVING AVERAGE FOR USE IN CONTROLLERS
C
      DO I=1,14
        MEAS(I,1)=THETA*MEAS(I,2) + (1.0D0 - THETA)*MEAS(I,1)
      ENDDO
C
C EVALUATE QUANTITATIVE CONSTRAINTS:
C
C INVENTORY, COOLING WATER PRESSURE DROP, EFFLUENT PRESSURE DROP
C
C INVENTORY CONSTRAINT
C
C <<<<<<<<<<<< CORRECTION: Flow 6 is increasing the mass, Flow 3 is reducing
      FINTEG=FINTEG+(MEAS(2,2)-MEAS(9,2)+FLOW(6)-FLOW(3))*DT ! MEAS(2,2)=FLOW(1), MEAS(9,2)=FLOW(4)
C      FINTEG=FINTEG+(MEAS(2,2)-MEAS(9,2))*DT ! MEAS(2,2)=FLOW(1), MEAS(9,2)=FLOW(4)
      MASSBAL=TAREA*MEAS(4,2)-FINTEG-3.00D0 ! MEAS(4,2)=L, L(t=0)=0.2 ==> A*L_0=3.0
C
C EFFLUENT FLOW CONSTRAINT
C
C     RCOMP(10) - COMPUTED VALUES FOR FLOW RESISTANCE BASED ON SENSOR MEASUREMENT DATA
C
      PBCOMP=RHO1*MEAS(4,2) ! MEAS(4,2)=L
      RCOMP(3)=5.0D0*EXP(0.0545D0*(100.0D0-MEAS(12,2))) ! MEAS(12,2)=100-CNT(1), MEAS(9,2)=FLOW(4)
      EPD=MEAS(9,2)-
     1 ((1D0/(RCOMP(3)+R0(1)+R0(4)))*(PBCOMP+PP0)**0.5D0)
C
C COOLING WATER FLOW CONSTRAINT
C
      RCOMP(6)=5.0D0*EXP(0.0545D0*(100.0D0-MEAS(13,2))) ! MEAS(13,2)=100-CNT(2), MEAS(8,2)=FLOW(5), MEAS(11,2)=PCW
      CWPD=MEAS(8,2)-
     1 ((1/(RCOMP(6)+R0(5)+R0(9)+R0(10)))*MEAS(11,2)**0.5D0)
C
C MOL BALANCE CONSTRAINT: MEAS: 1=CA0; 2=FLOW(1); 3=T(1); 4=L; 5=CA; 6=CB, 9=FLOW(4)
C
      sumConcABC = MEAS(5,2)+MEAS(6,2)+CCNOM  ! CA + CB + CC Nominal
      dmolin=MEAS(1,2)*MEAS(2,2)*DT  ! CA0 * FLOW(1) * DT
C <<<<<<<<<<<< CORRECTION: Flow 3 is reducing mol
      dmolout=sumConcABC*(MEAS(9,2)+FLOW(3))*DT ! (CA + CB + CC Nominal) * (FLOW(4)+FLOW(3)) * DT
C      dmolout=sumConcABC*MEAS(9,2)*DT ! (CA + CB + CC Nominal) * FLOW(4) * DT
      MOLIN=MOLIN+dmolin
      MOLOUT=MOLOUT+dmolout
      mol=sumConcABC*TAREA*MEAS(4,2)
      MOLBAL=mol-60.0D0-MOLIN+MOLOUT ! Mol balance, cA(t=0)*V(t=0)=20*3=60
C
C DETERMINE VALUES OF CONSTRAINTS FROM SENSORS
C
      MEAS(15,2)=MASSBAL
      MEAS(16,2)=CWPD
      MEAS(17,2)=EPD
      MEAS(18,2)=MOLBAL
C
C CONVERT TIME TO HOURS, MINUTES, AND SECONDS (DFRAC WAS NOT IN THE CODE)
C
      SEC=IDNINT(DFRAC(TIME)*60D0)
      TEMP=(TIME-DFRAC(TIME))/60D0
      HOUR=IDNINT(TEMP-DFRAC(TEMP))
      MINUTE=IDNINT((TIME-DFRAC(TIME))-(TEMP-DFRAC(TEMP))*60D0)
C
      IF (SEC .EQ. 60) THEN
        SEC=0
        MINUTE=MINUTE + 1
      ENDIF
C
C PRINT UPDATED STATUS
C
      IF ((ZCOUNT .EQ. IDNINT(ZLIM/3.0D0)).AND.(S8 .EQ. 1)) THEN
        SS=3
      ELSEIF ((ZCOUNT .EQ. IDNINT(ZLIM/1.5D0)).AND.(S8 .EQ. 1)) THEN
        SS=3
      ELSEIF ((ZCOUNT .EQ. IDNINT(ZLIM/2.0D0)).AND.(S8 .EQ. 1)) THEN
        SS=2
      ELSEIF (ZCOUNT .EQ. IDNINT(ZLIM)) THEN
        ZCOUNT=0
        SS=1
      ELSE
        GOTO 500
      ENDIF
C
C PRINT TO RAW DATA OUTPUT FILE
C
C LIST SENSOR READINGS IN RANDOM ORDER
C
C CALL GGUBS(A,B,C)
C
C GGUBS IS A UNIFORM RANDOM NUMBER GENERATOR PRODUCING A VECTOR C
C OF UNIFORM (0,1) RANDOM NUMBERS OF DIMENSION B. THE SEED (A) MUST
C BE DOUBLE PRECISION.
C
      DO 261 J=1,18
      RO(J,1)=SAMP(SS,J)
      RO(J,2)=SAMP(SS,J)
261   CONTINUE
      IF (S7 .EQ. 0) GOTO 266
      CALL GGUBS(DSEED,18,RV)
      DO 265 J=1,SAMP(SS,19)
      K=IDNINT(RV(J)*18.0D0+0.5D0)
262   IF (RO(K,1) .NE. 0) THEN
      RO(J,2)=RO(K,1)
      RO(K,1)=0
      ELSE
      K=K+1
      IF (K .GT. 18) K=1
      GOTO 262
      ENDIF
265   CONTINUE
C
C WRITE SENSOR READING TO FILE, char(9)=TAB
C
266   DO 270 J=1,SAMP(SS,19)
      WRITE (14,999) SENSORS(RO(J,2)),YEAR,MON,DAY,HOUR,MINUTE,SEC,
     1 MEAS(RO(J,2),2),UNITS(RO(J,2))
270   CONTINUE
C
C CHECK FOR TERMINATION AND ITERATE
C
C NOTE: SIMPLOT IS A SUBPROGRAM USED CREATING DATA FILES FOR PPLOT AND
C    CAN BE REMOVED WITHOUT OTHER MODIFICATIONS.
C
500   TIME=TIME + DT
      ZCOUNT=ZCOUNT + 1
      SP(3)=MEAS(14,2)  ! MEAS(14,2)=CNT(2)
C     write(*,*)'TIME=',TIME,'ZCOUNT=',ZCOUNT,'iter=',iter
      IF (TIME .GE. TH) THEN
C      IF (TIME .GT. TH) THEN
C      IF (TIME .GT. (TH + DT)) THEN
        WRITE (*,960)   ! 'END OF RUN'
        CLOSE (UNIT=14)
        WRITE (*,990)   ! 'PERFORM ANOTHER RUN? [0=NO , 1=YES]'
        READ (*,*) S5
        IF (S5 .EQ. 1) GOTO 5
        WRITE (*,995)   ! PLOT RESULTS? [0=NO , 1=YES]
        READ (*,*) S6
        IF (S6 .EQ. 1) THEN
C550       CONTINUE
          CALL SIMPLOT(DATAFILE,N,NORMVAL,MEAS,SDEV,SENSORS,UNITS)
        ENDIF

C ZCOUNT                - COUNTER FOR PRINTING OUTPUT
        write(prot,'(1x,A,f10.2,A,i6,A,i6,A,f10.2,A,f7.3)')
     1    'Closing at TIME=',TIME-DT,'  Iteration=',iter,
     1    '  Print output counter=',ZCOUNT,
     1    '  Time Horizon=',TH,'  Delta t=',DT
        write(*,*) 'Terminating with', nsmp,
     1             'data samples at iteration =', iter, 'TIME=',TIME
        write(*,*) 'Closing protocol file ''output/log.txt''...'
        close(prot)
        write(*,*) 'Closing data output file ''output/X.csv''...'
        close(dataout)
C       write(*,*)'M=',M,'fcnt=',fcnt

        STOP
      ENDIF
C
C ITERATE FOR NEXT TIME STEP
C
      iter=iter+1
C     CALL PEEPVARS(MEAS,classstr,iter)
C     write(*,*) 'iter=', iter, 'labelstr=', labelstr
C     write(*,*) 'RE=',RE,' RC=',RC
      if (MOD(iter,50).EQ.0) then
        CALL CLASSLABEL(S1, S2, F, M, TIME, DELAY, labelstr)
        CALL MEASOUT( dataout, MEAS, csvsep, labelstr, iter, nsmp )
C       write(*,*)'>>> Data dump at iter=',iter
C       stop
      endif

C      if (iter.eq.1) then
C      if (.FALSE.) then
C        stop
C      endif
      GOTO 50
C END OF MAIN LOOP
C
C==========================================================================
C
C FORMAT STATEMENTS
C
700   FORMAT (/,5X,'ENTER OUTPUT FILE NAME [DEFAULT=DUMP.RAW]')
710   FORMAT (A32)
C720   FORMAT (/,5X,'ENTER DATE STAMP [MM/DD/YYYY] ')
C730   FORMAT (I2,1X,I2,1X,I4)
735   FORMAT (/,5X,'PLOT RESULTS OF PREVIOUS RUN? [0=NO , 1=YES] ')
740   FORMAT (/,5X,'OVERRIDE CONTROLLER TUNING ? [0=NO , 1=YES] ')
C745   FORMAT (/,5X,'INDICATE CONTROLLER:',3X,'1=CONTROLLER',
C     1 /,28X,'2=TEMPERATURE CONTROLLER',
C     1 /,28X,'3=RECYCLE FLOW CONTROLLER',
C     1 /,28X,'4=COOLING WATER FLOW CONTROLLER ')
745   FORMAT (/,5X,'INDICATE CONTROLLER:',3X,'1=LEVEL CONTROLLER',
     1 /,28X,'2=TEMPERATURE CONTROLLER',
     1 /,28X,'3=COOLING WATER FLOW CONTROLLER ')
750   FORMAT (/,5X,'ENTER PARAMETER TO CHANGE:',3X,'1=GAIN',
     1 /,34X,'2=INTEGAL'
     1 /,34X,'3=DERIVATIVE ')
760   FORMAT (/,5X,'ENTER NEW VALUE [CURRENT VALUE IS ',F10.4,'] ')
770   FORMAT (/,5X,'ENTER FILTER CONSTANT [0.0 - 1.0] ')
772   FORMAT (/,5X,'PRINT OUTPUT IN RANDOM ORDER? [0=NO , 1=YES] ')
773   FORMAT (/,5X,'VARIABLE OUTPUT FREQUENCY? [0=NO , 1=YES] ')
C775   FORMAT (5X,F4.2)
790   FORMAT (/,5X,'ENTER INTEGER SEED FOR RANDOM NUMBER GENERATOR ')
800   FORMAT (//,5X,'*** JACKETED CSTR DYNAMIC SIMULATION ***',/)
810   FORMAT (5X,I2,') ',A40)
820   FORMAT (5X,'ENTER DESIRED FAULT TYPE ')
821   FORMAT (//,15X,'PROCESS SENSORS ARE : ',//,15(5X,I4,' ',A,/)
     1 ,//)
822   FORMAT (/,5X,'INDICATED FAULT HAS ALREADY BEEN SELECTED! ')
823   FORMAT (/,5X,'CONTINUE? [0=NO , 1=YES) ')
825   FORMAT (/,5X,'ENTER NUMBER OF FAILED SENSOR ')
830   FORMAT (/,5X,'ENTER TYPE OF SENSOR FAILURE: O=FIXED BIAS',
     1 /,35X,'1=FIXED VALUE ')
835   FORMAT (/,5X,'NOMINAL VALUE IS: ',3X,F14.4)
840   FORMAT (/,5X,'ENTER EXTENT OF FAULT (IN APPROPRIATE UNITS) ')
850   FORMAT (/,5X,'ENTER FAULT DELAY [min] ')
860   FORMAT (/,5X,'ENTER TIME CONSTANT [(min)-1]: ')
865   FORMAT (/,5X,'ADD ANOTHER FAULT? [0=NO , 1=YES] ')
870   FORMAT (/,5X,'ENTER TIME HORIZON FOR SIMULATION [min] ')
873   FORMAT (/,5X,'PRINT INTERMEDIATE RESULTS? [1=YES/0=NO] ')
875   FORMAT (/,5X,'ENTER PRINT INTERVAL [min] ')
880   FORMAT (//,10X,'PROCESS CONDITION IS -- ',A40,
     1   /,28X,'IN ',A32,
     2   //,10X,'SENSOR FAULT IS TYPE -- ',A12,
     3   /,10X,'FAULT EXTENT IS  -- ',F8.2,
     4   /,10X,'FAULT DELAY IS   -- ',F8.2,' (MIN)')
890   FORMAT (//,10X,'PROCESS CONDITION IS -- ',A40,
     1   //,10X,'FAULT EXTENT IS  --',F8.2,
     2   /,10X,'FAULT DELAY IS   --',F8.2,' (MIN)',
     3   /,10X,'TIME CONSTANT IS  --',F8.2)
895   FORMAT (/,10X,'FILTER CONSTANT IS  --',F8.2)
C930   FORMAT (/)
C940   FORMAT (5X,A36,' IS ',A16)
950   FORMAT (/,5X,'RUNNING.....',/)
960   FORMAT (/,5X,'END OF RUN',/)
970   FORMAT (//,5X,'***** EMERGENCY SHUTDOWN INITIATED AT ',F6.2,
     1   ' MIN ',//)
980   FORMAT (//,5X,'***** LOW LEVEL FORCES PUMP SHUTDOWN AT ',F6.2,
     1   ' MIN ',//)
990   FORMAT (/,5X,'PERFORM ANOTHER RUN? [0=NO , 1=YES] ')
995   FORMAT (/,5X,'PLOT RESULTS? [0=NO , 1=YES] ')
999   FORMAT (A32,I4,I2,I2,I2,I2,I2,F14.4,A7)
C
      END

C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C+++++++++++++  ROUTINES NOT IN THE ORIGINAL CODE  ++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


C----------------------------------------------------------------------
C Dump the measurements to the output file
      SUBROUTINE MEASOUT( dataout, MEAS, csvsep, labelstr, iter, nsmp )
      integer dataout, iter, nsmp
      REAL*8 MEAS(18,2)
      character*1 csvsep
      character *40 labelstr

      write(dataout,'(18(e15.5,A1),A1,A10)')
     1               (MEAS(I,2),csvsep,I=1,18),' ',labelstr
      nsmp=nsmp+1
      RETURN
      END
C----------------------------------------------------------------------


C----------------------------------------------------------------------
C Produce condition class label from current process condition
      SUBROUTINE CLASSLABEL(S1, S2, F, M, TIME, DELAY, labelstr)

      INTEGER *2 S1(3), S2(3), F(3), M
      REAL *8 DELAY(3), TIME
      INTEGER *2 sensfaultnr, faultnr
      INTEGER S0, I, POS, active_faults
      character *40 format_string, labelstr, condbuf

C      write(*,'(A15,I3,A10,3(I3),A5, 3(I3),A5, 3(I3),
C     1         A7,F15.3,A7,3(F15.3))')
C     1              'CLASSLABEL> M=',M, 'FAULTS=', (F(S0),S0=1,3),
C     1              'S1=', (S1(S0),S0=1,3),
C     1              'S2=', (S2(S0),S0=1,3),
C     1              'TIME=',TIME,'DELAY=', (DELAY(S0),S0=1,3)
C      stop

      labelstr='normal'
      IF (M .EQ. 1 .AND. F(1) .EQ. 1) RETURN  ! No fault

      active_faults = 0
      DO S0=1,M
        IF (TIME .GE. DELAY(S0)) active_faults = active_faults + 1
      END DO

      IF (active_faults .EQ. 0) RETURN

      labelstr='                    '
      POS = 1

      DO S0=1,active_faults
        faultnr = F(S0)
        condbuf='                    '
        IF (faultnr .EQ. 23) THEN
          sensfaultnr = S1(S0)
          IF (sensfaultnr .lt. 10) THEN
            format_string = '(A1,I1,A1)'
          ELSE
            format_string = '(A1,I2,A1)'
          END IF
          write(condbuf,format_string) 'S', sensfaultnr, ' '
        ELSE ! No sensor fault
          IF (faultnr .lt. 10) THEN
            format_string = '(I1,A1)'
          ELSE
            format_string = '(I2,A1)'
          END IF
          write(condbuf,format_string) faultnr, ' '
        END IF
C        write(*,*) 'S0=', S0, 'condbuf=>>>', condbuf, '<<<'
        I = 1
        DO WHILE ( condbuf(I:I) .NE. ' ' )
          labelstr(POS:POS) = condbuf(I:I)
C           write(*,*) '>>>>> POS=', POS, 'I=', I, 'labelstr=',labelstr
          I = I + 1
          POS = POS + 1
        END DO
        IF (S0 .LT. active_faults) THEN
          labelstr(POS:POS) = '+'
          POS = POS + 1
        END IF
      END DO

      RETURN
      END
C-----------------------------------------------------------------------


C----------------------------------------------------------------------
C Peep measured variables
      SUBROUTINE PEEPVARS(MEAS,labelstr,iter)

      REAL *8 MEAS(18,2)
      character *40 labelstr
      INTEGER *4 iter
C
      character *16 dummy, csvsep*1
      csvsep = ';'
C      WRITE(*,*) 'Peeping variables ...'
      write(*,'(A5,I6,18(e15.5,A1),A1,A10)')
     1               'k=',iter,(MEAS(I,2),csvsep,I=1,18),' ',labelstr
C      write(*, *) '...'
C      read(*,*) dummy
      RETURN
      END

C-----------------------------------------------------------------------



C----------------------------------------------------------------------
C Fractional part of double precision real value
      REAL*8 FUNCTION DFRAC( X )
      REAL*8 X
      DFRAC = X-INT(X)
      RETURN
      END
C----------------------------------------------------------------------



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C DUMMY IMPLEMENTATION

      SUBROUTINE SIMPLOT(DATAFILE,N,NORMVAL,MEAS,SDEV,SENSORS,UNITS)

      CHARACTER*8 UNITS(18)
      CHARACTER *32 SENSORS(18),DATAFILE,gnufile*40
      INTEGER *4 N, i
      REAL *8 NORMVAL(18),SDEV(18),MEAS(18,2)
C
      WRITE(*,*) 'Executing SIMPLOT ...'
C      WRITE(*,*) 'dummy function. Not yet implemented ...'
C      gnufile=adjustr(trim(DATAFILE)//'.gnu')
C      OPEN (UNIT=15,FILE=gnufile)
C      write(*,*) 'Opening plot file ''',gnufile,''''
C      write(*,*) 'n=',N
C      write(*,*) (i,SENSORS(i),i=1,N)
C      write(*,*) (i,':',NORMVAL(i),' ',MEAS(i,1),' ',SDEV(i),i=1,N)
C
C
C      write(*,*) 'Closing plot file ''',gnufile,''''
C      CLOSE(UNIT=15)
      RETURN
      END

C-----------------------------------------------------------------------
C
C
C p. 319:
C ... Note that the Jacketed CSTR Simulation program makes
C three external function calls to GGNML, GGUBS, and SIMPLOT. GGNML and GGUBS
C are IMSL random number generator routines. SIMPLOT is a custom plotting routine. The
C SIMPLOT call can be removed with no effect on the performance of the simulator
C
C
C http://www.netlib.org/list/imsl
C imsl/GGNML   NORMAL OR GAUSSIAN RANDOM DEVIATE GENERATOR
C imsl/GGUBS   BASIC UNIFORM (0, 1) PSEUDO-RANDOM NUMBER GENERATOR
C
C
C DSEED: inital seed (DSEED is automatically updated)
C NR: dimension of the double precision array R
C R: REAL*8 R(NR)
C GGUBS generates NR independent uniform random variables
C GGNML generates NR independent normal variables


      SUBROUTINE GGUBS(DSEED,NR,R)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   GGUBS GENERATES NR SINGLE PRECISION RANDOM VARIATES UNIFORM C
C ON (0,1) BY A LINEAR CONGRUENTIAL SCHEME.  THIS ROUTINE IS    C
C DEPENDENT ON MACHINE WORD SIZE.                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C   ON ENTRY                                                    C
C       DSEED   DOUBLE PRECISION                                C
C               SEED FOR GENERATOR                              C
C       NR      INTEGER                                         C
C               NUMBER OF VARIATES TO GENERATE                  C
C   ON RETURN                                                   C
C       R       REAL (NR)                                       C
C               SINGLE PRECISION ARRAY CONTAINING THE VARIATES  C
C   GGUBS CALLS                                                 C
C               DMOD                                            C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
      INTEGER            NR
      REAL*8             R(NR)
      DOUBLE PRECISION   DSEED
C
C                              LOCAL
C
      INTEGER            I
      DOUBLE PRECISION   D2P31M,D2P31,DMULTX
C
C                              MACHINE CONSTANTS
C                              D2P31M=(2**31) - 1
C                              D2P31 =(2**31)(OR AN ADJUSTED VALUE)
C
      DATA               D2P31M/2147483647.D0/
      DATA               D2P31/2147483711.D0/
      DATA               DMULTX/16807.0D+00/
C
      DO 5 I=1,NR
         DSEED=DMOD(DMULTX*DSEED,D2P31M)
         R(I) =DSEED / D2P31
  5   CONTINUE
C
C                               END OF GGUBS
C
      RETURN
      END
C===============================================================C
C===============================================================C
      SUBROUTINE GGNML(DSEED,NR,R)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   GGNML GENERATES NR SINGLE PRECISION N(0,1) RANDOM VARIATES  C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C   ON ENTRY                                                    C
C       DSEED   DOUBLE PRECISION                                C
C               SEED FOR GENERATOR                              C
C       NR      INTEGER                                         C
C               NUMBER OF VARIATES TO GENERATE                  C
C   ON RETURN                                                   C
C       R       REAL (NR)                                       C
C               SINGLE PRECISION ARRAY CONTAINING THE VARIATES  C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
      INTEGER            NR
      REAL*8             R(NR)
      DOUBLE PRECISION   DSEED
C                              LOCAL
      INTEGER             IER
C
C                              GET NR RANDOM NUMBERS
C                              UNIFORM (0,1)
C
      CALL GGUBS(DSEED,NR,R)
C
C                              TRANSFORM EACH UNIFORM DEVIATE
C
      DO 5 I=1,NR
         CALL MDNRIS(R(I),R(I),IER)
    5 CONTINUE
C
C                               END OF GGNML
C
      RETURN
      END


      SUBROUTINE MDNRIS (P,Y,IER)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL*8             P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL*8             EPS,G0,G1,G2,G3,H0,H1,H2,A,W,WI,SN,SD
      REAL*8             SIGMA,SQRT2,X,XINF
      DATA               XINF/1.7014E+38/
      DATA               SQRT2/1.414214/
      DATA               EPS/1.1921E-07/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      IF (P .GT. 0.0 .AND. P .LT. 1.0) GO TO 5
      IER = 129

      if (p .lt. 0.0D+00) then
         sigma = -1.0D+00
      else
         sigma = 1.0D+00
      endif

C      SIGMA = SIGN(1.0,P)
C      write(6,666) p,sigma
C666   format(' Sign #1: ',2f15.8)
C      pause

      Y = SIGMA * XINF
      GO TO 20
    5 IF(P.LE.EPS) GO TO 10
      X = 1.0D0 -(P + P)
      CALL MERFI (X,Y,IER)
      Y = -SQRT2 * Y
      GO TO 20
C                                  P TOO SMALL, COMPUTE Y DIRECTLY
   10 A = P+P
      W = SQRT(-dLOG(A+(A-A*A)))
C                                  USE A RATIONAL FUNCTION IN 1./W
      WI = 1.0D0/W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      Y = W + W*(G0+SN/SD)
      Y = -Y*SQRT2
C                               END OF MDNRIS
  20  RETURN
      END
C===============================================================C
      SUBROUTINE MERFI (P,Y,IER)
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C                                                               C
C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
C   LAST MODIFIED                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                                               C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$C
C                                  SPECIFICATIONS FOR ARGUMENTS
      REAL*8             P,Y
      INTEGER            IER
C                                  SPECIFICATIONS FOR LOCAL VARIABLES
      REAL*8             A,B,X,Z,W,WI,SN,SD,F,Z2,RINFM,A1,A2,A3,B0,B1,
     *                   B2,B3,C0,C1,C2,C3,D0,D1,D2,E0,E1,E2,E3,F0,F1,
     *                   F2,G0,G1,G2,G3,H0,H1,H2,SIGMA
      DATA               A1/-.5751703/,A2/-1.896513/,A3/-.5496261E-1/
      DATA               B0/-.1137730/,B1/-3.293474/,B2/-2.374996/
      DATA               B3/-1.187515/
      DATA               C0/-.1146666/,C1/-.1314774/,C2/-.2368201/
      DATA               C3/.5073975E-1/
      DATA               D0/-44.27977/,D1/21.98546/,D2/-7.586103/
      DATA               E0/-.5668422E-1/,E1/.3937021/,E2/-.3166501/
      DATA               E3/.6208963E-1/
      DATA               F0/-6.266786/,F1/4.666263/,F2/-2.962883/
      DATA               G0/.1851159E-3/,G1/-.2028152E-2/
      DATA               G2/-.1498384/,G3/.1078639E-1/
      DATA               H0/.9952975E-1/,H1/.5211733/
      DATA               H2/-.6888301E-1/
      DATA               RINFM/1.7014E+38/
C                                  FIRST EXECUTABLE STATEMENT
      IER = 0
      X = P

      if (x .lt. 0.0D+00) then
         sigma = -1.0D+00
      else
         sigma = 1.0D+00
      endif

C      SIGMA = SIGN(1.0,X)
C      write(6,666) x,sigma
C666   format(' Sign #2: ',2f15.8)
C      pause

C                                  TEST FOR INVALID ARGUMENT
      IF (.NOT.(X.GT.-1. .AND. X.LT.1.)) GO TO 30
      Z = ABS(X)
      IF (Z.LE. .85D0) GO TO 20
      A = 1.0D0-Z
      B = Z
C                                  REDUCED ARGUMENT IS IN (.85,1.),
C                                     OBTAIN THE TRANSFORMED VARIABLE
C   5 W = SQRT(-dLOG(A+A*B))
      W = SQRT(-dLOG(A+A*B))
      IF (W.LT.2.5D0) GO TO 15
      IF (W.LT.4.0D0) GO TO 10
C                                  W GREATER THAN 4., APPROX. F BY A
C                                     RATIONAL FUNCTION IN 1./W
      WI = 1.0D0/W
      SN = ((G3*WI+G2)*WI+G1)*WI
      SD = ((WI+H2)*WI+H1)*WI+H0
      F = W + W*(G0+SN/SD)
      GO TO 25
C                                  W BETWEEN 2.5 AND 4., APPROX. F
C                                     BY A RATIONAL FUNCTION IN W
   10 SN = ((E3*W+E2)*W+E1)*W
      SD = ((W+F2)*W+F1)*W+F0
      F = W + W*(E0+SN/SD)
      GO TO 25
C                                  W BETWEEN 1.13222 AND 2.5, APPROX.
C                                     F BY A RATIONAL FUNCTION IN W
   15 SN = ((C3*W+C2)*W+C1)*W
      SD = ((W+D2)*W+D1)*W+D0
      F = W + W*(C0+SN/SD)
      GO TO 25
C                                  Z BETWEEN 0. AND .85, APPROX. F
C                                     BY A RATIONAL FUNCTION IN Z
   20 Z2 = Z*Z
      F = Z+Z*(B0+A1*Z2/(B1+Z2+A2/(B2+Z2+A3/(B3+Z2))))
C                                  FORM THE SOLUTION BY MULT. F BY
C                                     THE PROPER SIGN
   25 Y = SIGMA*F
      IER = 0
      GO TO 40
C                                  ERROR EXIT. SET SOLUTION TO PLUS
C                                     (OR MINUS) INFINITY
   30 IER = 129
      Y = SIGMA * RINFM
C                               END OF MERFI
   40 RETURN
      END
C===============================================================C
