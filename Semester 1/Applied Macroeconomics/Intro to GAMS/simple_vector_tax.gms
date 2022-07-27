* The model with one household and one producer
* 1. Fill in the gaps
* 2. Check if model is homogeneous
* 3. Check if walras law holds
* 4. Run a simulation of an increase of one factor's endowment
* 5. Add a uniform factor use tax on all sectors/factors
* 6. How are the results different from when you introduced the tax on only one
* good and factor?
* 7. When taxes are distortionary?

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

SETS
      SEC /X,Y/
      FAC /K,L/
      ;

PARAMETERS
      BETA(FAC,SEC)
      ALPHA(SEC)
      XS0(SEC)
      W0
      OMEGA(FAC)
      USE0(FAC,SEC)
      TAX(FAC,SEC)
      ;

* USE0 is the factor use by sector
USE0("K","X") = 60;
USE0("L","X") = 40;
USE0("K","Y") = 40;
USE0("L","Y") = 60;

* Benchmark outputs (from SAM)
XS0(SEC) = 100;


* Assigns numbers from SAM to parameters:
* Parameters of the Cobb Douglas production function
BETA(FAC,SEC) = USE0(FAC,SEC) / XS0(SEC);

* Parameters of the Cobb Dougas utility function
ALPHA(SEC) = 0.5;

* Benchmark utility
W0 = 200;

* Endowments of factors
OMEGA(FAC) = 100;

* Important: we have to introduce the tax in all sectors/factors but we set it to zero.
TAX(FAC,SEC) = 0;

* The Cobb Douglas scaling factor
* You need that in order for the output and total consumption aggregate
* to be exactly as much as you want.

PARAMETERS
A
B(SEC);

* For consumption:
A    =  W0 / (PROD(SEC, XS0(SEC)**ALPHA(SEC)));
* For Production
B(SEC) =  XS0(SEC) / ( PROD(FAC, USE0(FAC,SEC) ** BETA(FAC,SEC)));


VARIABLES
      XS(SEC)                 Output
      P(SEC)                  Prices of final goods
      W                       Welfare index
      PW                      Consumer price index
      PFAC(FAC)         Factor wages
      INC                      Consumer income
      TRICK
;



EQUATIONS
      PRF_XS(SEC)       Zero profits on goods production production
      PRF_W             Determination of the consumer price index

      MKT_XS(SEC)       Market clearing for goods
      MKT_W             Demand for aggregate consumption

      MKT_FAC(FAC)      Market clearing for factors

      I_INC             Income balance
      TRCK
;

PRF_XS(SEC)..
      PROD(FAC, (PFAC(FAC)/BETA(FAC,SEC))**BETA(FAC,SEC)) / B(SEC)   =E= P(SEC);

PRF_W..
      PROD(SEC, (P(SEC) / ALPHA(SEC))**ALPHA(SEC))  / A =E= PW;

MKT_XS(SEC)..
      XS(SEC) =E= ALPHA(SEC) * W * PW / P(SEC);

MKT_W..
      W =E= INC / PW  ;

MKT_FAC(FAC)..    OMEGA(FAC) =E=  FILL_IN ;

I_INC..                  INC =E=  FILL_IN ;


TRCK..  TRICK=E=1;

MODEL SIMPLE
      /PRF_XS
     PRF_W
     MKT_XS
     MKT_W
     MKT_FAC
     I_INC
     TRCK
/
;

* Provide starting values for the solver

      XS.L(SEC)= XS0(SEC);
      W.L= W0;
      P.L(SEC)= 1;
      PW.L= 1;
      PFAC.L(FAC)= 1;
      INC.L= 200;

* Provide lower bounds for the solver
      XS.LO(SEC)= 0.00001;
      W.LO= 0.00001;
      P.LO(SEC)= 0.00001;
      PW.LO= 0.00001;
      PFAC.LO(FAC)= 0.00001;
      INC.LO= 0.00001;

* Provide a numeraire:
PFAC.FX("K")=1;

*simple.iterlim=0;

SOLVE  SIMPLE USING NLP maximizing TRICK;

