* The model with one household and one producer
* 1. Fill in the gaps
* 2. Check if model is homogeneous
* 3. Check if walras law holds
* 4. Introduce a welfare function that will be a sum of welfare
* * of households and government
* 5. Introduce a tax on capital in production of good 4. Interpret the results
* 6. Run a tax-replacement scenario: cut the income taxes by 5pp at the same time increasing
* * a tax on capital in production of good 4. Make sure that government consumption remains unchanged.
* * Hint: you need to change one parameter to become a variable and this variable needs to be endogeneous.

$ONTEXT
* -----------------------------------------------------
        X1        X2        X3        X4        L        K        HOU1        GOV        INV
X1                                                                68,8          20         17.2
X2                                                                32            30         8
X3                                                                64            14         16
X4                                                                52            10         13
L        46        15        64        10
K        60        55        30        65
HOU1                                        135        210
GOV                                                               74
SAV                                                               54.2          0

* -----------------------------------------------------

$OFFTEXT

SETS
      SEC /X1, X2, X3, X4/
      FAC /K,L/
      HOU /HOU1/
      ;

ALIAS(FACC,FAC);
ALIAS(SECC,SEC);
ALIAS(HOUSS,HOU);

PARAMETERS
      XS0(SEC)
      W0(HOU)
      OMEGA0(HOU,FAC)      Endowment of the household
      USE0(FAC,SEC)         Input-output matrix
      INC0(HOU)           Gross income
      TRANSFER(HOU,HOU)
      BETA(FAC,SEC)         Share parameter of the prod. f-n
      ALPHA(HOU,SEC)      Share parameter of the utility f-n
      B(SEC)                Shift parameter of the C-D production f-n
      A(HOU)              Scaling parameter of the utility function
      SAV0(HOU)           Initial savings
      SR(HOU)             Savings rate
      ;

TABLE FACTORUSE(FAC,SEC)
      X1   X2   X3   X4
L     46   15   64   10
K     60   55   30   65
;

USE0(FAC,SEC) = FACTORUSE(FAC,SEC);

XS0(SEC)= SUM(FAC, USE0(FAC,SEC));

BETA(FAC,SEC) = USE0(FAC,SEC) / SUM(FACC, USE0(FACC,SEC));

B(SEC) =  XS0(SEC) / ( PROD(FAC, USE0(FAC,SEC) ** BETA(FAC,SEC)));


* Calibration of the government
PARAMETER XDG0(SEC) /
X1    20
X2    30
X3    14
X4    10
/;

* Parameters of govt C-D function
PARAMETERS
ALPHAG(SEC)
AG
WG0;

WG0 = SUM(SEC, XDG0(SEC));
ALPHAG(SEC) =  XDG0(SEC) / SUM(SECC, XDG0(SECC));
AG = WG0 /  (PROD(SEC, XDG0(SEC)**ALPHAG(SEC)));

PARAMETER
TAXREV0(HOU)  Tax revenues /HOU1 74/
ITAX(HOU)     Income tax rate
FTAX(FAC,SEC)   Tax on factor use;




* Endowments of factors
TABLE ENDOWMENT(HOU,FAC)
         L    K
HOU1   135  210
;
OMEGA0(HOU,FAC) = ENDOWMENT(HOU,FAC);
INC0(HOU) = SUM(FAC, OMEGA0(HOU,FAC));
SAV0("HOU1")  =   54.2;


TABLE XD0(SEC,HOU)
     HOU1
X1   68.8
X2   32
X3   64
X4   52
;

* Initial transfers will be zero
TRANSFER(HOU,HOU) = 0;
* Benchmark utility
W0(HOU) = SUM(SEC, XD0(SEC,HOU));

Display XD0, W0;
* Parameters of the Cobb Dougas utility function
ALPHA(HOU,SEC) =  XD0(SEC,HOU) / W0(HOU);

A(HOU)    =  W0(HOU) /  (PROD(SEC, XD0(SEC,HOU)**ALPHA(HOU,SEC)));

ITAX(HOU) = TAXREV0(HOU) / INC0(HOU);

**** savings rate is based on disposable (net) income
SR(HOU)   =  SAV0(HOU) / (INC0(HOU) - TAXREV0(HOU));


FTAX(FAC,SEC) = 0;

*** Calibration of investment
PARAMETER INVD0(SEC)
/
X1 17.2
X2 8
X3 16
X4 13
/;

PARAMETERS
INV0             Total investmestment
ALPHAI(SEC)      Parameter of the investment production function
AI               Scale parameter of the production function
;
INV0 = SUM(SEC, INVD0(SEC));
ALPHAI(SEC) =  INVD0(SEC) / SUM(SECC, INVD0(SECC));
AI = INV0 /  (PROD(SEC, INVD0(SEC)**ALPHAI(SEC)));

display itax;
VARIABLES
      XS(SEC)           Output
      P(SEC)            Prices of final goods
      INVD(SEC)         Investment demand
      INV               Aggregate investment
      W(HOU)          Welfare index
      WG                Index of govt consumption
      PW(HOU)         Consumer price index
      PWG               Govt price index
      PINV              Investment price index
      PFAC(FAC)         Factor wages
      INC(HOU)        Consumer income
      INCG              Government budget
      TRICK             Trick var
      GSAV              Govt balance
      INC_SHIFT         Income tax rate shifter

;



EQUATIONS
      PRF_XS(SEC)       Zero profits on goods production production
      PRF_W(HOU)      Determination of the consumer price index
      PRF_WG            Determination of the govt price index
      PRF_INV           Determination of the investment price index
      MKT_XS(SEC)       Market clearing for goods
      MKT_W (HOU)     Demand for aggregate consumption
      MKT_WG            Demand for govt consumption
      MKT_INV           Demand for investment
      MKT_FAC(FAC)      Market clearing for factors
      I_INC(HOU)      Income balance (households)
      I_INCG            Income balance (govt)
      TRCK              THe trick equation
;


* Here we have to introduce the tax
PRF_XS(SEC)..
      (PROD(FAC, ( PFAC(FAC) * (1 + FTAX(FAC,SEC)) / BETA(FAC,SEC) )**BETA(FAC,SEC))) / B(SEC)   =E= P(SEC);

PRF_W(HOU)..
      (PROD(SEC, (P(SEC) / ALPHA(HOU,SEC) )**ALPHA(HOU,SEC)))  / A(HOU) =E= PW(HOU)
;

PRF_WG..
       (PROD(SEC, (P(SEC) / ALPHAG(SEC) )**ALPHAG(SEC)))  / AG =E= PWG
;

PRF_INV..
      (PROD(SEC, (P(SEC) / ALPHAI(SEC) )**ALPHAI(SEC)))  / AI =E= PINV
;


MKT_XS(SEC)..
      XS(SEC) =E= SUM(HOU,ALPHA(HOU,SEC) * W(HOU) * PW(HOU) / P(SEC))
                         +  ALPHAG(SEC) * WG * PWG / P(SEC)
                         +  ALPHAI(SEC) * INV * PINV / P(SEC)
;

MKT_W(HOU)..
      PW(HOU) * W(HOU) =E= INC(HOU) *  (1-ITAX(HOU)*INC_SHIFT)
                                                       * (1- SR(HOU));  ;
;

MKT_WG..
      PWG * WG =E= INCG - GSAV;
;

MKT_FAC(FAC)..     SUM(HOU, OMEGA0(HOU,FAC)) =E=
                   SUM(SEC, BETA(FAC,SEC) * P(SEC) * XS(SEC) / (PFAC(FAC)
                                                                 * (1 + FTAX(FAC,SEC)))) ;

;

I_INC(HOU)..      INC(HOU) =E= SUM(FAC, PFAC(FAC) * OMEGA0(HOU,FAC)) + SUM(HOUSS, TRANSFER(HOUSS,HOU))
                                                                           - SUM(HOUSS, TRANSFER(HOU,HOUSS));


I_INCG..            INCG =E=  SUM(HOU, INC(HOU) *  ITAX(HOU) * INC_SHIFT)
                         +   SUM(FAC,  SUM(SEC, FTAX(FAC,SEC) * PFAC(FAC) * BETA(FAC,SEC) * P(SEC) * XS(SEC) / (PFAC(FAC)
                                                                 * (1 + FTAX(FAC,SEC)))));

TRCK..  TRICK=E=1;

MODEL SIMPLE
      / PRF_XS
        PRF_W
        PRF_WG
        MKT_XS
        MKT_W
        MKT_WG
        MKT_FAC
        I_INC
        I_INCG
        PRF_INV
        TRCK

        /
;


* Provide starting values for the solver

      XS.L(SEC)= XS0(SEC);
      W.L(HOU)= W0(HOU);
      WG.L = WG0;

      P.L(SEC)= 1;
      PW.L(HOU)= 1;
      PWG.L = 1;
      PFAC.L(FAC)= 1;
      PINV.L=1;
      INV.L = INV0;
      INVD.L(SEC) = INVD0(SEC);

      INC.L(HOU)= INC0(HOU) ;
      INCG.L = WG0;


* Provide lower bounds for the solver
      XS.LO(SEC)= 0.00001;
      W.LO(HOU)= 0.00001;
      P.LO(SEC)= 0.00001;
      PW.LO(HOU)= 0.00001;
      PFAC.LO(FAC)= 0.00001;
      INC.LO(HOU)= 0.00001;
      TRICK.L=1;




PARAMETER
NGDP
GDP
GDP0
CHANGE_GDP
;

* Provide a numeraire:
PFAC.FX("L")=1;
INC_SHIFT.FX = 1;

* Set the GOVT balance
*GSAV.FX = 0;

WG.FX = 1.5*WG0;

option nlp = pathnlp ;
SOLVE  SIMPLE USING NLP maximizing TRICK;
NGDP = SUM(SEC, P.L(SEC) * XS.L(SEC));
GDP  = SUM(SEC,  XS.L(SEC));
GDP0 = SUM(SEC, XS0(SEC));
CHANGE_GDP = 100 *(GDP / GDP0 -1);
Display NGDP, GDP, GDP0, CHANGE_GDP;

*$ontext
* Calculate steady state values of depr, growth and interest rate
PARAMETERS
r        rate of return on capital
depr     depreciation                    /0.04/
gp       population growth               /0.00/;
r = SUM(HOU, OMEGA0(HOU, "K")) * (gp + depr) / INV0;

FTAX(FAC,SEC) = 0;
;

*Second period
* make sure capital accumulates by investment and labor grows as well
* Note that we take into account PREVIOUS period investment
OMEGA0(HOU, "L") = (1 + gp) * OMEGA0(HOU, "L");
OMEGA0(HOU, "K") = (1 - depr) *OMEGA0(HOU, "K")
+ r * INV.L * OMEGA0(HOU, "K") / SUM(HOUSS,OMEGA0(HOUSS, "K"));

option nlp = pathnlp ;
SOLVE  SIMPLE USING NLP maximizing TRICK;

NGDP = SUM(SEC, P.L(SEC) * XS.L(SEC));
GDP  = SUM(SEC,  XS.L(SEC));
GDP0 = SUM(SEC, XS0(SEC));
CHANGE_GDP = 100 *(GDP / GDP0 -1);
Display NGDP, GDP, GDP0, CHANGE_GDP;
*$offtext

*Third period

* make sure capital accumulates by investment and labor grows as well
* Note that we take into account PREVIOUS period investment
OMEGA0(HOU, "L") = (1 + gp) * OMEGA0(HOU, "L");
OMEGA0(HOU, "K") = (1 - depr) *OMEGA0(HOU, "K")
+ r * INV.L * OMEGA0(HOU, "K") / SUM(HOUSS,OMEGA0(HOUSS, "K"));

option nlp = pathnlp ;
SOLVE  SIMPLE USING NLP maximizing TRICK;

NGDP = SUM(SEC, P.L(SEC) * XS.L(SEC));
GDP  = SUM(SEC,  XS.L(SEC));
GDP0 = SUM(SEC, XS0(SEC));
CHANGE_GDP = 100 *(GDP / GDP0 -1);
Display NGDP, GDP, GDP0, CHANGE_GDP;