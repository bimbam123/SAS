%macro bowker(data=_last_, var1=, var2=, weight=, where=, k=10);

/*
Exact (permutation) Bowker symmetry test for K x K frequency table.

Parameters:

data name of input dataset (default=_last_)
var1,var2 names of row & column variables (required)
weight name of a case weighting variable (default=none)
where where clause as %str() to be applied to PROC FREQ
(default=none)
k max. floating point value is 2**(2**k)-1 (default=10)

Restrictions:

No. of rows I & no. of columns J of the frequency table resulting
from crosstabulating var1 with var2 should satisfy 1<I=J<10.
The special case of I=J=2 corresponds to McNemar's test, for which
PROC FREQ supplies an exact test.

PROC FREQ must produce a square table - use a weighting variable
with zero weights if necessary to achieve this. If all you have are
cell frequencies, you can use a weighting variable to quickly set
up an input dataset - only off-diagonal cells need be entered.
Macro will fail if any weight is other than a non-negative integer.
Note that SAS releases prior to 9.1 do not support zero weights,
but if you're able to make a square table without zero weights the
macro will work for earlier releases if you remove the "/zeros"
code.

The k parameter is used to avoid floating point overflows. Reduce
the
value if 'Invalid argument(s) to the exponential operator "**"'
message appears in log. Due to the issue of underflows, increasing k
could increase accuracy of the p-value if your system will support
larger numbers. Default value is for SAS Release 9.1.3 for
Windows/Intel.

Output listing generated:

PROC FREQ output including AGREE option material
Number of permutations, weighted permutations
Exact (2-tailed) p-value for symmetry test

Output dataset generated:

_out cell frequencies for K x K table

Example of use:

data v;
input x y f; datalines;
1 2 8
1 3 0
2 1 0
2 3 1
3 1 0
3 2 0
run;
%bowker(data=v, var1=x, var2=y, weight=f, where=%str(x<3 and y<3))

Exact Bowker symmetry test
Number of permutations computed = 18
Number of weighted permutations = 2**9
Exact p-value = 0.0078125

Written by Jay Weedon
Last revised November 22, 2005
Please forward comments/bug reports to .

Rough outline of algorithm:

For each pair of cells Cij & Cji (i<j) having observed frequencies
Fij & Fji, compute Tij=Fij+Fji, then generate all possible
permutations of frequencies in Cij & Cji that sum to Tij, and
compute the Bowker chi-square statistic for each permutation.
Since the permutations are not equiprobable, a weight is assigned
to each permutation corresponding to the product of relevant
binomial coefficients. The sum of weights of chi-square values
greater than or equal to the observed value is used to generate
p-value.

The number of permutations tested is equal to the product of terms
(Tij+1) for all i<j. The number of weighted permutations that this
corresponds to is equal to the product of terms 2**(Tij) for all
i<j.
*/

*Generate cell frequencies;
proc freq data=&data;
tables &var1*&var2 / agree out=_out sparse;
%if &weight^= %then %do;
weight &weight / zeros;
%end;
%if &where^= %then %do;
where &where;
%end;
run;

*Copy off-diagonal cell frequencies into macro vars fij.
Size of table also copied into macro var size;
data _null_;
set _out end=last;
by &var1 &var2;
if first.&var1 then do;
_i+1; _j=1;
end;
else _j+1;
if _i^=_j then call
symput('f'||put(_i,1.)||put(_j,1.),put(count,best16.));
if last then call symput('size',put(_i,1.));
run;

*Compute tij=fij+fji, then Bowker statistic as well
as number of permutations needed and base-2 log of
sum of weighted permutations;
%let bowker=0; %let npermute=1; %let log2sumweights=0;
%do i=1 %to &size-1;
%do j=&i+1 %to &size;
%let t&i&j=%eval(&&f&i&j+&&f&j&i);
%if &&t&i&j>0 %then %let
bowker=%sysevalf(&bowker+(&&f&i&j-&&f&j&i)**2/&&t&i&j);
%let npermute=%eval(&npermute*(&&t&i&j+1));
%let log2sumweights=%eval(&log2sumweights+&&t&i&j);
%end;
%end;

*Check whether overflow would occur - if so, use normalization
factor;
%if &log2sumweights<2**&k %then %let adjust=0;
%else %let
adjust=%sysfunc(floor(%sysfunc(log2(&log2sumweights))-&k+1));
%let factor=%sysevalf(&log2sumweights*(1-2**-&adjust));
%let factorln2=%sysevalf(%sysfunc(log(2))*&factor);
%put Computing &npermute permutations...;

*Generate permutations and calculate chi-square & weight;
data _null_;

*Write data step DO statements;
%do i=1 %to &size-1;
%do j=&i+1 %to &size;
%if &&t&i&j>0 %then %do;
do f&i&j=0 to &&t&i&j;
%end;
%end;
%end;

*Calculate chi-square for current permutation;
chisquare=
%do m=1 %to &size-1;
%do n=&m+1 %to &size;
%if &&t&m&n>0 %then %do;
+(f&m&n+f&m&n-&&t&m&n)**2/&&t&m&n
%end;
%end;
%end;;

*If chi-square for permutation is greater than or equal
to observed value, compute weight, which is product
of binomial(Tmn,Fmn) for all m<n. Logs of weights are
used to prevent overflows;
if chisquare>=&bowker then do;
logeweight=
%do m=1 %to &size-1;
%do n=&m+1 %to &size;
%if &&t&m&n>0 %then %do;
+lgamma(&&t&m&n+1)
-lgamma(f&m&n+1)
-lgamma(&&t&m&n-f&m&n+1)
%end;
%end;
%end;;
cumweight+exp(logeweight-&factorln2);
end;

*Write data step END statements;
%do i=1 %to &size-1;
%do j=&i+1 %to &size;
%if &&t&i&j>0 %then %do;
end;
%end;
%end;
%end;

p=cumweight/(2**(&log2sumweights-&factor));
if fuzz(p)=1 then p=1;
file print;
put / @22 "Exact Bowker symmetry test"
/ @22 "Number of permutations computed = &npermute"
/ @22 "Number of weighted permutations = 2**&log2sumweights";
if p>0 then put @22 "Exact p-value = " p;
else put @22 "The Bowker exact test cannot be computed with
sufficient precision for this sample size";
run;
%mend;
 
