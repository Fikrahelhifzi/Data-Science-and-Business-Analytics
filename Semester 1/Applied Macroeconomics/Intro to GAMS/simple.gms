* The model with one household and one producer
* 1. Fill in the gaps
* 2. Check if model is homogeneous
* 3. Check if walras law holds
* 4. Run a simulation of an increase of one factor's endowment

$ONTEXT
            X     Y     L     K     W     CONS  Total
X                             100                     100
Y                             100                     100
W                                         200         200
L           40    60                                  100
K           60    40                                  100
CONS              100   100                     200
Total      100    100   100   100   200   200
$OFFTEXT

PARAMETERS
      BETA_KX
      BETA_KY
      BETA_LX
      BETA_LY
      ALPHA_X
      ALPHA_Y
      W0
      OMEGA_K
      OMEGA_L
      A
      B_X
      B_Y
      ;

* Assigns numbers from SAM to parameters:
* Parameters of the Cobb Douglas production function
BETA_KX = 60/100;
BETA_LX = 40/100;
BETA_KY = 40/100;
BETA_LY = 60/100;


* Parameters of the Cobb Dougas utility function
ALPHA_X = 0.5;
ALPHA_Y = 0.5;

* Endowments of factors
OMEGA_K = 100;
OMEGA_L = 100;

A    = 200 / (100**ALPHA_X * 100**ALPHA_Y);
* For Production
B_X =  100 / ( 40 ** BETA_LX * 60**BETA_KX);
B_Y =  100 / ( 40 ** BETA_KY * 60**BETA_LY);

POSITIVE VARIABLES
      XS                Output X
      YS                Output Y
      PX                Price of X
      PY                Price of Y
      PK                Wage of K
      PL                Wage of L

      PW                Consumer price index
      INC               Consumer income

;

VARIABLE
      TRICK             Trick variable
      W                 Welfare index
;


EQUATIONS
      PRF_XS       Zero profits on goods production production
      PRF_YS       Zero profits on goods production production
      PRF_W        Determination of the consumer price index

      MKT_XD       Market clearing for X
      MKT_YD       Market clearing for Y
      MKT_W        Demand for aggregate consumption

      MKT_K        Labor market clearing
      MKT_L        Capital market clearing

      I_INC        Income balance
      TRCK         Trick equation
;


PRF_XS..
      ((PK/BETA_KX)**BETA_KX) * ((PL/BETA_LX)**BETA_LX) / B_X =E= FILL IN;
PRF_YS..
      ((PK/BETA_KY)**BETA_KY) * ((PL/BETA_LY)**BETA_LY) / B_Y =E= FILL IN;


PRF_W..
      ((PX/ALPHA_X)**ALPHA_X) *  ((PY/ALPHA_Y)**ALPHA_Y) / A =E= PW;

MKT_XD..
      XS  =E= ALPHA_X * W * PW / PX;

MKT_YD..
      YS  =E= ALPHA_Y * W * PW / PY;


MKT_W..
      PW * W =E= INC ;

MKT_K..
      OMEGA_K =E= BETA_KX * PX * XS / PK + BETA_KY * PY * YS / PK;

MKT_L..
      FILL_IN

I_INC..
      INC =E= OMEGA_K * PK + OMEGA_L * PL;

TRCK..
      TRICK=E=1;

MODEL SIMPLE /
PRF_XS,
PRF_YS,
PRF_W,
MKT_XD,
MKT_YD,
MKT_W,
MKT_K,
MKT_L,
I_INC,
TRCK
/;


*Numeraire
PK.FX=1;

PK.L=1;
PL.L=1;
PX.L=1;
PY.L=1;
PW.L=1;
W.L=200;
XS.L=100;
YS.L=100;
INC.L=200;

*OMEGA_K=1.2*OMEGA_K;

option nlp = pathnlp ;
SOLVE  SIMPLE USING NLP maximizing Trick;

** An alternative formulation with explicit maximization of utility
$ontext
EQUATIONS
Utility
;

Utility..
     W =E= A*XS**ALPHA_X*YS**ALPHA_Y;

MODEL SIMPLE2 /
PRF_XS,
PRF_YS,
Utility,
MKT_W,
MKT_K,
MKT_L,
I_INC
/;

SOLVE  SIMPLE2 USING NLP maximizing W;
$offtext
