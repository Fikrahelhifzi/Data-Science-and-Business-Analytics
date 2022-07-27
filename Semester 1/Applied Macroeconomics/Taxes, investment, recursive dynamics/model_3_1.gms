$ONTEXT
* -----------------------------------------------------
      X1   X2   X3   X4   L    K    W1   W2   WG   CONS1 CONS2 CONSG
X1                                  60   26   20
X2                                  15   25   30
X3                                  40   40   14
X4                                  10   55   10
L     46   15   64   10
K     60   55   30   65
W1                                                 125
W2                                                       146
WG                                                             74
CONS1                     100  55
CONS2                     35   155
CONSG                                              30    44
* -----------------------------------------------------

$OFFTEXT

SETS
      SEC /X1, X2, X3, X4/
      FAC /K,L/
      HOU /HOU1,HOU2/
      ;

ALIAS(FACC,FAC);
ALIAS(SECC,SEC);
ALIAS(HOUSS,HOU);

PARAMETERS
      XD0(SEC)
      W0(HOU)
      OMEGA0(HOU,FAC)
      USE0(FAC,SEC)
      CONS0(HOU,SEC)
      TAX(FAC,SEC)
      TRANSFER(HOU,HOU)
      BETA(FAC,SEC)      Share parameter of the prod. f-n
      ALPHA(HOU,SEC)   Share parameter of the utility f-n
      B(SEC)      Shift parameter of the C-D production f-n
      A(HOU)     Scaling parameter of the utility function
      ;

TABLE FACTORUSE(FAC,SEC)
      X1   X2   X3   X4
L     46   15   64   10
K     60   55   30   65
;

USE0(FAC,SEC) = FACTORUSE(FAC,SEC);

* Benchmark outputs (from SAM)
XD0(SEC)= SUM(FAC, USE0(FAC,SEC));


BETA(FAC,SEC) = USE0(FAC,SEC) / SUM(FACC, USE0(FACC,SEC));

B(SEC) =  XD0(SEC) / ( PROD(FAC, USE0(FAC,SEC) ** BETA(FAC,SEC)));

* Endowments of factors
TABLE ENDOWMENT(HOU,FAC)
         L    K
HOU1     100  55
HOU2     35   155
;
OMEGA0(HOU,FAC) = ENDOWMENT(HOU,FAC);

TABLE CONSUMPTION(SEC,HOU)
    HOU1   HOU2
X1  60       26
X2  15       25
X3  40       40
X4  10       55
;
CONS0(HOU,SEC) = CONSUMPTION(SEC,HOU);

display CONS0;

TRANSFER(HOU,HOU) = 0;

* Parameters of the Cobb Dougas utility function
ALPHA(HOU,SEC) =  CONS0(HOU,SEC) / SUM(SECC, CONS0(HOU,SECC));

* Benchmark utility
W0(HOU) = SUM(SEC, CONS0(HOU,SEC));

A(HOU)    =  W0(HOU) /  (PROD(SEC, CONS0(HOU,SEC)**ALPHA(HOU,SEC)));


* Calibration of the government
PARAMETER CONSG0(SEC) /
X1    20
X2    30
X3    14
X4    10
/;
display consg0;


* Parameters of govt C-D function
PARAMETERS
ALPHAG(SEC)
AG
WG0;

WG0 = SUM(SEC, CONSG0(SEC));
ALPHAG(SEC) =  CONSG0(SEC) / SUM(SECC, CONSG0(SECC));
AG = WG0 /  (PROD(SEC, CONSG0(SEC)**ALPHAG(SEC)));

PARAMETER TAXREV0(HOU)  Tax revenues /HOU1 30, HOU2 44/;
PARAMETER ITAX(HOU)     Income tax rate;
PARAMETER GINC0(HOU)    Household gross income;

GINC0(HOU) = SUM(FAC, OMEGA0(HOU,FAC));
ITAX(HOU) = TAXREV0(HOU) / GINC0(HOU);

display itax;
*$ontext

VARIABLES
      XD(SEC)           Output
      P(SEC)            Prices of final goods
      W(HOU)          Welfare index
      WG                Index of govt consumption
      PW(HOU)         Consumer price index
      PWG               Govt price index
      PFAC(FAC)         Factor wages
      INC(HOU)        Consumer income
      INCG              Government budget
      FTAX(FAC,SEC)
      FTAXK
      TRICK             Trick var
;



EQUATIONS
      PRF_XD(SEC)       Zero profits on goods production production
      PRF_W(HOU)      Determination of the consumer price index
      PRF_WG            Determination of the govt price index
      MKT_XD(SEC)       Market clearing for goods
      MKT_W (HOU)     Demand for aggregate consumption
      MKT_WG            Demand for govt consumption
      MKT_FAC(FAC)      Market clearing for factors
     I_INC(HOU)      Income balance (households)
      I_INCG            Income balance (govt)
      TRCK              THe trick equation
      FTAXEQ(SEC)      Equation for uniform tax rate
;

PRF_XD(SEC)..
      (PROD(FAC, ( PFAC(FAC) * (1 + FTAX(FAC,SEC)) / BETA(FAC,SEC) )**BETA(FAC,SEC))) / B(SEC)   =E= P(SEC);

PRF_W(HOU)..
      (PROD(SEC, (P(SEC) / ALPHA(HOU,SEC) )**ALPHA(HOU,SEC)))  / A(HOU) =E= PW(HOU)
;

PRF_WG..
         (PROD(SEC, (P(SEC) / ALPHAG(SEC) )**ALPHAG(SEC)))  / AG =E= PWG
;

MKT_XD(SEC)..
      XD(SEC) =E= SUM(HOU,ALPHA(HOU,SEC) * W(HOU) * PW(HOU) / P(SEC))
                         +  ALPHAG(SEC) * WG * PWG / P(SEC)
;

MKT_W(HOU)..
      PW(HOU) * W(HOU) =E= INC(HOU)
;

MKT_WG..

        WG*PWG  =E= INCG;

;

MKT_FAC(FAC)..     SUM(HOU, OMEGA0(HOU,FAC)) =E=
                   SUM(SEC, BETA(FAC,SEC) * P(SEC) * XD(SEC) / (PFAC(FAC)
                                                                 * (1 + FTAX(FAC,SEC)) ))

;

I_INC(HOU)..      INC(HOU) =E= (SUM(FAC, PFAC(FAC) * OMEGA0(HOU,FAC)) + SUM(HOUSS, TRANSFER(HOUSS,HOU))
                                                                          - SUM(HOUSS, TRANSFER(HOU,HOUSS))

                                         ) * (1-ITAX(HOU));


I_INCG..            INCG =E=  SUM(HOU, (SUM(FAC, PFAC(FAC) * OMEGA0(HOU,FAC)) + SUM(HOUSS, TRANSFER(HOUSS,HOU))
                                                                          - SUM(HOUSS, TRANSFER(HOU,HOUSS))
                                                                           ) * ITAX(HOU))

                         + SUM(FAC,  SUM(SEC, FTAX(FAC,SEC) * PFAC(FAC) * BETA(FAC,SEC) * P(SEC) * XD(SEC) / (PFAC(FAC)
                                                                 * (1 + FTAX(FAC,SEC)))));
;

FTAXEQ(SEC)..
                    FTAXK =E= FTAX("K",SEC);

TRCK..  TRICK=E=1;

MODEL SIMPLE
      /ALL/
;

* Provide starting values for the solver

      XD.L(SEC)= XD0(SEC);
      W.L(HOU)= W0(HOU);
      WG.L = WG0;

      P.L(SEC)= 1;
      PW.L(HOU)= 1;
      PWG.L = 1;
      PFAC.L(FAC)= 1;

      INC.L(HOU)= W0(HOU);
      INCG.L = WG0;


* Provide lower bounds for the solver
      XD.LO(SEC)= 0.00001;
      W.LO(HOU)= 0.00001;
      P.LO(SEC)= 0.00001;
      PW.LO(HOU)= 0.00001;
      PFAC.LO(FAC)= 0.00001;
      INC.LO(HOU)= 0.00001;
      TRICK.L=1;
* Provide a numeraire:
*FTAX.FX(FAC,SEC) = 0;
FTAX.FX("L",SEC)=0;

PFAC.FX("L")=1;
*WG.FX=WG0;
WG.FX=83.3037;
option nlp = pathnlp ;

*simple.iterlim=0;
*FTAX.FX("K",SEC) = 0.0605;
* We want WG to be 83.3037;
* So the FTAX.FX("K",SEC) = 0.0605;
SOLVE  SIMPLE USING NLP maximizing TRICK;

display XD.L,W.L,P.L,PW.L,PFAC.L,INC.L;

PARAMETERS
W_IDX(HOU)
WG_IDX
XD_IDX(SEC)
CONS_IDX(HOU,SEC)
CONSG_IDX(SEC)
WELFARE_IDX;


*** Recompute the display variables
W_IDX(HOU) = W.L(HOU)/W0(HOU);
WG_IDX = WG.L / WG0;
XD_IDX(SEC) = XD.L(SEC) / XD0(SEC);
CONS_IDX(HOU,SEC) = ALPHA(HOU,SEC) * W.L(HOU) * PW.L(HOU) / P.L(SEC) / CONS0(HOU,SEC);
CONSG_IDX(SEC) = ALPHAG(SEC) * WG.L * PWG.L / P.L(SEC) / CONSG0(SEC);
WELFARE_IDX = (SUM(HOU, W.L(HOU))  + WG.L) /  (SUM(HOU, W0(HOU))  + WG0);

display XD.L,W.L,P.L,PW.L,PFAC.L,INC.L, WG.L;
display W_IDX, WG_IDX, XD_IDX, CONS_IDX, CONSG_IDX, WELFARE_IDX;
display ftax.l;
*$offtext