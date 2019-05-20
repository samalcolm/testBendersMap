$ontext

   An example of Benders Decomposition on fixed charge transportation
   problem bk4x3.

   Optimal objective in reference : 350.

   Erwin Kalvelagen, December 2002

   See:
   http://www.in.tu-clausthal.de/~gottlieb/benchmarks/fctp/
   http://www.gams.com/~erwin/benders/benders.pdf

$offtext

set i 'sources' /i1*i4/;
set j 'demands' /j1*j3/;

parameter supply(i) /
   i1 10
   i2 30
   i3 40
   i4 20
/;

parameter demand(j) /
   j1 20
   j2 50
   j3 30
/;

* modify nodelist to include category of node
file nodelist /"D:\RShiny\testBendersMap\nodes.ben"/
put nodelist;
put "id        city        lon            lat          Category"  /
put "i1        Seattle    -122.3320708    47.6062095   Source"   /
put "i2        San_Diego  -117.1610838    32.715738    Source"   /
put "i3        Dallas     -96.7969879     32.7766642   Source"     /
put "i4        Memphis    -90.0489801     35.1495343   Source"   /
put "j1        New_York   -74.0059728     40.7127753   Demand"  /
put "j2        Chicago    -87.6297982     41.8781136   Demand"   /
put "j3        Topeka     -95.6751576     39.0473451   Demand"   /

table c(i,j) 'variable cost'
       j1    j2    j3
i1    2.0   3.0   4.0
i2    3.0   2.0   1.0
i3    1.0   4.0   3.0
i4    4.0   5.0   2.0
;

table f(i,j) 'fixed cost'
        j1     j2     j3
i1    10.0   30.0   20.0
i2    10.0   30.0   20.0
i3    10.0   30.0   20.0
i4    10.0   30.0   20.0
;

*
* check supply-demand balance
*
scalar totdemand, totsupply;
totdemand = sum(j, demand(j));
totsupply = sum(i, supply(i));
abort$(abs(totdemand-totsupply)>0.001) "Supply does not equal demand.";

*
* for big-M formulation we need tightest possible upperbounds on x
*
parameter xup(i,j) 'tight upperbounds for x(i,j)';
xup(i,j) = min(supply(i),demand(j));

variables
   cost   'objective variable'
   x(i,j) 'shipments'
   y(i,j) 'on-off indicator for link'
;
positive variable x;
binary variable y;

*---------------------------------------------------------------------
* Benders Decomposition Initialization
*---------------------------------------------------------------------

display "--------------------- BENDERS ALGORITHM ----------------------------";

scalar UB 'upperbound' /INF/;
scalar LB 'lowerbound' /-INF/;

y.l(i,j) = 1;


*---------------------------------------------------------------------
* Benders Subproblem
*---------------------------------------------------------------------

variable z 'objective variable';

positive variables
   u(i) 'duals for capacity constraint'
   v(j) 'duals for demand constraint'
   w(i,j) 'duals for xy constraint'
;

equations
   subobj          'objective'
   subconstr(i,j)  'dual constraint'
;

* to detect unbounded subproblem
scalar unbounded /1.0e6/;
z.up = unbounded;

subobj..  z =e= sum(i, -supply(i)*u(i)) + sum(j, demand(j)*v(j))
                + sum((i,j), -xup(i,j)*y.l(i,j)*w(i,j))
                 ;

subconstr(i,j)..  -u(i) + v(j) - w(i,j) =l= c(i,j);

model subproblem /subobj, subconstr/;
* reduce output to listing file:
subproblem.solprint=2;
* speed up by keeping GAMS in memory:
subproblem.solvelink=2;

*---------------------------------------------------------------------
* Benders Modified Subproblem to find unbounded ray
*---------------------------------------------------------------------

variable dummy 'dummy objective variable';

equations
   modifiedsubobj          'objective'
   modifiedsubconstr(i,j)  'dual constraint'
   edummy;
;

modifiedsubobj..
    sum(i, -supply(i)*u(i)) + sum(j, demand(j)*v(j))
           + sum((i,j), -xup(i,j)*y.l(i,j)*w(i,j)) =e= 1;

modifiedsubconstr(i,j)..
    -u(i) + v(j) - w(i,j) =l= 0;

edummy.. dummy =e= 0;

model modifiedsubproblem /modifiedsubobj, modifiedsubconstr, edummy/;
* reduce output to listing file:
modifiedsubproblem.solprint=2;
* speed up by keeping GAMS in memory:
modifiedsubproblem.solvelink=2;


*---------------------------------------------------------------------
* Benders Relaxed Master Problem
*---------------------------------------------------------------------

set iter /iter1*iter50/;

set cutset(iter) 'dynamic set';
cutset(iter)=no;
set unbcutset(iter) 'dynamic set';
unbcutset(iter)=no;


variable z0 'relaxed master objective variable';
equations
   cut(iter)           'Benders cut for optimal subproblem'
   unboundedcut(iter)  'Benders cut for unbounded subproblem'
;

parameters
   cutconst(iter)     'constant term in cuts'
   cutcoeff(iter,i,j)
;

cut(cutset).. z0 =g= sum((i,j), f(i,j)*y(i,j))
                      + cutconst(cutset)
                      + sum((i,j), cutcoeff(cutset,i,j)*y(i,j));
unboundedcut(unbcutset)..
                cutconst(unbcutset)
                + sum((i,j), cutcoeff(unbcutset,i,j)*y(i,j)) =l= 0;

model master /cut,unboundedcut/;
* reduce output to listing file:
master.solprint=2;
* speed up by keeping GAMS in memory:
master.solvelink=2;
* solve to optimality
master.optcr=0;

*---------------------------------------------------------------------
* Benders Algorithm
*---------------------------------------------------------------------

scalar newincumbent;
scalar converged /0/;
scalar iteration;
scalar bound;
parameter ybest(i,j);
parameter log(iter,*) 'logging info';
parameter curr_solution;
file out /"D:\RShiny\testBendersMap\iterations.out"/;
out.ap=0;
put out;
put ""  /;
putclose out;
out.ap=1;

file conv /"D:\RShiny\testBendersMap\convergence.out"/;
conv.ap=0;
put conv;
put ""  /;
putclose conv;
conv.ap=1;

loop(iter$(not converged),

*
* solve Benders subproblem
*
   solve subproblem maximizing z using lp;
newincumbent = 0;

*
* check results.
*
if (subproblem.modelstat>=2, put conv);
   abort$(subproblem.modelstat>=2) "Subproblem not solved to optimality";

*
* was subproblem unbounded?
*

   if (z.l+1 < unbounded,

*
* no, so update upperbound
*
      bound = sum((i,j), f(i,j)*y.l(i,j)) + z.l;
      if (bound < UB,
          UB = bound;
          ybest(i,j) = y.l(i,j);
          display ybest;
* Hooray, a new integer solution has been found!
newincumbent = 1;
* Write solution information
LOOP((i,j)$(subconstr.M(i,j) > 0),
put out;
PUT  iter.TL ","  i.TL ","  j.TL ","  subconstr.M(i,j)   /;
);
putclose out;

      );

*
* and add Benders' cut to Relaxed Master
*
      cutset(iter) = yes;

   else

*
* solve modified subproblem
*

     solve modifiedsubproblem maximizing dummy using lp;

*
* check results.
*

     abort$(modifiedsubproblem.modelstat>=2)
            "Modified subproblem not solved to optimality";


*
* and add Benders' cut to Relaxed Master
*
      unbcutset(iter) = yes;
   );


*
* cut data
*
   cutconst(iter) = sum(i, -supply(i)*u.l(i)) + sum(j, demand(j)*v.l(j));
   cutcoeff(iter,i,j) = -xup(i,j)*w.l(i,j);

*
* solve Relaxed Master Problem
*

   option optcr=0;
   solve master minimizing z0 using mip;

*
* check results.
*

   abort$(master.modelstat=4) "Relaxed Master is infeasible";
   abort$(master.modelstat>=2) "Masterproblem not solved to optimality";

*
* update lowerbound
*

   LB = z0.l;

   log(iter,'LB') = LB;
   log(iter,'UB') = UB;
   curr_solution(iter,i,j)=ybest(i,j);
   iteration = ord(iter);
* file to save convergence information
put conv;
put ord(iter) "," UB "," LB "," newincumbent /;
putclose conv;
display iteration,LB,UB;

   converged$( (UB-LB) < 0.1 ) = 1;
   display$converged "Converged";

);

display log;

abort$(not converged) "No convergence";

