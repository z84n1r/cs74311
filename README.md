java c
Exercises for Stats305a
December 7, 2022
0    Matrix review questions
Question 0.1:   Let R be a square right-triangular (upper triangular) matrix, that is,

Assume diag(R) is all non-zero. Give a (short!) algorithm to solve for x in Rx = b.
Question 0.2 (Projections):    Let A ∈ Rm×n  where m < n and A has full row rank (that is, rank(A) = m), and let b ∈ Rm. Consider the affine space (subspace plus an ofset)
S := {x ∈ Rn  | Ax = b}.
(a) Let y ∈ Rn. Give the (Euclidean) projection


(b) Draw a picture of your result above.
Question 0.3:   A Householder  refiection (or transformation) is Hu  = I — 2uuT  where u is a unit vector and I is the identity.
(a) Show that Hu  is symmetric and unitary, meaning Hu(T) = Hu  and Hu(T)Hu  = I.
(b) Draw a picture of the mapping x }→ Hux exhibiting why this is called a refiection.
(c) We show how to refiect a vector x about the line between the direction of x and the ﬁrst standard basis vector e1 , so that the transform. Hux is on the line {te1   |  t  ∈ R}.  Let x ∈ Rn  be an arbitrary vector, and deﬁne

(taking v = 0 if x = —te1  for some t > 0). Show that for Hvx = — ⅡxⅡ e1 .
(d) Let u1 , . . . , uk  be arbitrary unit vectors. Show that the following matrix is unitary: k

(e) Let A ∈ Rn×n  be full rank with ﬁrst column a1 . Give a Householder transformation (by specifying the unit vector u) so that

where *s are arbitrary numbers, 0n-1  ∈ Rn-1  is all zeros, and B ∈ Rn-1×n-1  is a matrix.
(f) Given a matrix A with the block structurewhere R is square and upper triangular, B ∈ Rk×k  is square and has ﬁrst column b1 , and 0 is an all-zeros matrix of appropriate size, give a symmetric unitary matrix S so that
(g) Describe (in one or two sentences) how one may use Householder transformation to con- struct a QR factorization of any full rank matrix A.Question 0.4 (The power method):    Let A ∈ Rn×n be diagonalizable, meaning A = SΛS-1 for a diagonal matrix Λ = diag(λ1 , . . . , λn), and let S have columns S = [s1   · · ·  sn]. Assume the largest eigenvalue is unique, with jλ1 j  > jλjj for each j  1, and consider the iteration
where Ⅱ·Ⅱ is the usual l2-norm.Recall that a diagonalizable matrix has left and right eigenvectors, where left eigenvectors satisfy vT A = λvT , that is, ATv = λv, while right eigenvectors satisfy Av = λv.  (Note that the left and right eigenvalues are necessarily identical.) Let v1 , . . . , vn  be the left eigenvectors of A, and assume that v 1(T)x0  ≠  0, that is, the initial iterate x0  is not  orthogonal to the left eigenvector v1  corresponding to the largest eigenvalue λ 1 .
(a) To what does xk  converge? (b) Prove it.Hint.  Convince yourself that the left eigenvectors are the rows of the matrix S-1 . Because S is full rank, we can write any x ∈ Rn  as x = SQ for some Q ∈ Rn. Show that if x0 = SQ then Q1   0 as v 1(T)x0   0, then develop the iteration for xk .Question 0.5 (Block matrix inversion, linear algebra, and the Sherman-Morrison-Woodbury formula):     In this question, we will work out a few formulae for inverting block-structured matrices, using them to develop various inversion formulas for structured matrices.
(a) Consider the matrix equation
where A and D are square matrices and B, C have appropriate sizes (this is not important for this question).  Give a formula for x in terms of A,B, C, D, a, and b; your formula, if correct, will involve inverses of some of these.  You may assume that A  and D  are invertible and that A — BD-1 C  is invertible.  (Note:  B and C may be rectangular, so don’t try to invert them alone.)
(b) We now consider inverting a matrix plus a (typically) low rank matrix. We wish to solve
(A + UCVT)x = z
for x, where A ∈ Rn×n , U  ∈ Rn×k , C  ∈ Rk×k , and V  ∈ Rn×k .  Introducing the variable y = CV Tx, or VT x — C-1y = 0, solve


for x. Using this, show that
x = (A-1 — A-1 U(C-1 + VT A-1 U)-1 VT A-1) z,
i.e.,
(A + UCVT)-1 = A-1 — A-1 U(C-1 + VT A-1 U)-1 VT A-1 .
(As an aside, if you know A-1  and C-1  already, the largest matrix you must invert to compute (A + UCVT)-1  is then k × k, which is smaller.)Question 0.6 (Majorization inequalities):     A matrix P  ∈ Rn×n   is a permutation  matrix if its entries are all {0, 1}-valued and PT P = I , that is, P has a single 1 in each row and column. Let Pn be the collection of n × n permutation matrices. The Birkhof polytope is the convex hull of the permutation matrices and coincides with the doubly  stochastic  matrices, where we recall a matrix S ∈ Rn×n  is doubly stochastic if S1 = ST 1 = 1 and its entries are nonnegative. That is, letting Sn  be the doubly stochastic matrices, we have
so each S ∈ Sn  is a convex combination of permutation matrices.
(a) Let a1  ≥ a2  and b1  ≥ b2 . Show that
a1 b1 + a2 b2  ≥ a2 b1  + a1 b2 .
(b) Let u, v ∈ Rn  and assume v1  ≥ v2  ≥ · · · ≥ vn  and u1  ≥ u2  ≥ · · · ≥ un. Show that uTv ≥ uTP v  for all P ∈ Pn.You may use that if σ : [n] → [n] is any permutation, there is a sequence of transpositions (i.e., swaps of two elements, so that if σI  and σ are identical except that σ(i) = σI(i + 1) and σ(i + 1) = σI(i) for some index i, then they are transpositions) that transform. σ into the identity permutation.
(c) Assume u and v are as in part (a). Show that

Question 0.7:   Let A, B ∈ Rn×n. The von  Neumann  trace  inequality  states that

where σ 1 (A) ≥ · · · ≥ σn(A) ≥ 0 denote the singular values of A (and similarly for σi(B)). In this question, you will demonstrate this inequality.
(a) Show that it is no loss of generality to assume that A is diagonal, that is, to show that inequality (vNti) holds when A = diag(a1 , . . . , an) with a1  ≥ · · · ≥ an  ≥ 0.
(b) Let B  = U ΣVT   where U  =  [u1    · · · un] and V  = [v1    · · ·  v n] are unitary, and A  = diag(a1 , . . . , an) where a i  ≥ 0. Show that

(c) Under the same conditions as in part (b), show that

Hint. For any x, y ∈ R, we have xy ≤ 2/1x2 + 2/1y2 because 0 ≤ 2/1(x − y)2 = 2/1x2 − xy + 2/1y2.
(d) Using the preceding parts and Question 0.6, show the trace inequality (vNti).
1    A few preliminary questions
Question 1.1 (A rank-one update to a linear regression solution):    Consider a least-squares problem with data X ∈ Rn×d , Y ∈ Rn, where n ≥ d and X has rank d, and let

a(G)iv(v)o(s)o,m(r)g(u),ll2(-)X,)o(,-1+o,1-.ie(})t.do(ha)e(t)s(m)co(u)lm(ti)p(p)lu(y)ting(ing)yo(H)u(b)r(y)
update take? (Note: you can simply say “A few multiples of d,” or “A few multiples of d2 ,” or “A few multiples of d3 ,” depending on which is true.)Question 1.2 (The most negatively correlated distribution):     A ﬁnancial analyst tells you that he has a great stock tip that will allow you to short stocks based on others that are doing well. He assures youthathe has found a correlation matrix between the n stocks, C ∈ Rn×n , with entries

(a) Is his correlation matrix possible? Would you trust him with your graduate stipend? (b) More generally, consider a correlation matrix Cρ  with entries of the following form.	
where P  ≥  0.   How large can P  be while  Cρ   remains  a potentially valid correlation matrix?  Hint:   A correlation matrix for a random vector X  ∈  Rn   has entries Cij   = Cov(Xi, Xj)/√Var(Xi)Var(Xj), and so may be written
C = diag(v)-1/2Cov(X) diag(v)-1/2 ,
where v is a vector with entries vi   = Var(Xi) and diag(v) is the diagonal matrix with diagonal v.Question  1.3  (Some basic plotting and data processing):      The UCI Machine Learning repository has a collection of useful datasets for experiments and data exploration. This is a question that simply serves as a forcing function for you to pick a computer language, read in data, and plot it appropriately. Using the data in the UCI Wine quality dataset (https: //archive.ics.uci.edu/ml/datasets/Wine+Quality, and see the winequality-red .csv ﬁle in the data folder there), plot a scatterplot matrix showing the ﬁve variables density, alcohol, pH, volatile.acidity, and the target variable quality, which is a measure of wine quality. Such pairwise scatterplots can be a useful tool for data exploration and summarization (see, e.g., Fig. 1.1 of [5]).
In your plots, what do you notice about density, alcohol, and quality?  (Just a sentence su伍ces here.)Question  1.4  (Predicting high temperatures at SFO):      In this question, we use linear regression to predict the high temperature at San Francisco International Airport (SFO). A natural model of the temperature is the following: let x be the day of the year (that is, x = 1 corresponds to January 1, while x = 365 corresponds to December 31, except in leap years); we assume that the temperature
(Admittedly this ignores the issue of leap years, but we will punt on that.)  Let φ(x) =  - T  be the vector feature representation above.The data ﬁle simplified-sfo-weather .csv contains weather data for SFO since 1960, including precipitation (in inches), low, and high temperatures (in degrees Fahrenheit). Note that in May 2018, the temperature sensor broke and consequently a few days report NA as the high and low temperatures.  You should simply omit those from any averages or model ﬁtting.1    Using this data, ﬁt the model (1.1) to predict high temperature (this is column "temphigh" in the ﬁle) from the date x of the year for years prior to 1990. For each decade 1961–1970, 1971–1980, 1981–1990, 1991–2000, 2001–2010, and 2011–2020, print the mean of the actual high temperature minus the predicted high temperature for the decade.  That is,
if β = [β0  β1  β2]T  denotes your ﬁt model, compute the average diference y - y- = y - φ(x)T β-
over each of those decades. Include your code and a printout of the results.Question 1.5:   We consider monitoring changes in rainfall/precipitation over the years at San Francisco International Airport (SFO) using the data in simplified-sfo-weather .csv.代 写Exercises for Stats305aR
代做程序编程语言 To do so, we will set up a standard linear model with d = 3 dimensions, where for dates (times) t ∈ {1, 2, 3, . . . , 366} (we have 366 for leap years) we setx = [1  sin  (t - 1))   cos  (t - 1))]T  ∈ Rd                                   (1.2)
where d = 3. Under the standard linear model
yi = xi(T)β +  N(0, σ2 ),    i = 1, 2, 3, . . . ,
we would like to test whether future data follows a similar distribution to the past data. We begin with a few mathematical generalities. Consider two datasets modeled by
Y = Xβ + ε,     Ynew  = Zβ + εnew ,where X ∈ Rm×d  and Z ∈ Rn×d  are the given covariates, and we model ε ~ N(0, σ2 Im) and εnew ~ N(0, σ2 In) independently; we will think of (X, Y) as the initial data pair and (Z, Ynew) as the new data. (Their particular form. is immaterial; we assume both X and Z are rank d.)
Let β = (XTX)-1XTY , be the usual least-squares estimate on the “initial” data pair (X, Y), let HX  = X(XTX)-1XT  ∈ Rm×m  be the usual hat matrix, and deﬁne the predicted values
Y := Xβ = HXY  and  Ynew  := Zβ .
1 In R, you may do this automatically in the lm methods using the keyword na.action  =  na.omit, and in computing a mean using na.rm  =  TRUE.
((b(a)))  G(Sh)ive(ow)ath(symm(at the)etric(resi)d, p(u)o(a)ls(s)itiv(Y)e—ﬁn(an)ite(d)ma(new)tri(—) R(re)n(i)n so(dep)et(n)h(d)a(e)t(n)t.
where In  denotes the n × n identity matrix.
(c) Give the distribution of the ratiounder the null hypothesis                                         (1.3)
H0  : { ,   ε 丄 εnew.We now consider implementing a series of hypothesis tests about whether the precipitation at SFO is remaining consistent over the years or whether it is changing in some meaningful way.
(d) For each of the years 1966, 1967, . . ., 2020, repeat the following. Deﬁne a data matrix X using the featurization (1.2) consisting of all  dates prior to that year (so that for 1966, X will be a data matrix for the years 1960–1965, for 1967, X will be the data for years 1960–1966, and so on). Deﬁne the responses Y to consist of precipitation (column precip in the csv ﬁle) for the given years.  Deﬁne the new data matrix Z  ∈ Rn×d  to consist of the n = 365 (or 366 in a leap year) rows in the given year and the responses Ynew  to be the precipitation in the given year.  For this data, compute the statistic A in Eq. (1.2) and its p-value, that is, conditional on A = a, report
p := P(A ≥ a)
where A follows the null above.Give a printed list of your p values for each of the years, and provide a plot of the p values for each of the years as well.  Provide a few sentences of discussion about the observed p-values.
Please include your code in your solution.
(e) For a ﬁxed year, suppose we perform. the procedure in part (d) for only that year, getting a single p-value.  Suppose that this p-value is very small, for example, something like p  =  10-5   or p  =  10-6 .  In one or two sentences, does rejecting the null hypothesis necessarily mean that the distribution of precipitation is changing over time?
2    Sampling error and testingQuestion 2.1 (Regression with random samples and best linear approximations):    Say that (X, Y) ∈ Rd  × R come from a joint probability distribution P , where E[XXT] = C > 0 and X, Y both have ﬁnite second moments. Assume that the ﬁrst coordinate of X  is constant, with X1  = 1. For x ∈ Rd, deﬁne the regression function
f大 (x) := E[Y j X = x].
(a) Show that f minimizes E[(Y − f (X))2] over all functions f : Rd  → R.
(b) Instead of ﬁtting a model of Y  j  X  over the space of all functions, consider ﬁtting one over all linear predictors, and choosing β大 to minimize the expected squared loss
over b  ∈ Rd.  Characterize the solution β大   (i.e., give a formula for it), and show that the linear function ϕ(x) = β大Tx is the best linear approximation to f  (in mean-squared distance). Note:  in case you are worried about it, it is ﬁne to exchange expectation and diferentiation in this case; you deﬁnitely don’t need to show that, though it is  true (see, for example, Bertsekas [2]).
(c) Say that we have an i.i.d. sample (xi, yi), i = 1, . . . , n from P , with yi  = f (xi) + εi, and

X\j  = [x(1)   · · · x(j-1)  x(j+1)   · · · x(d)] ,
which is a submodel as we have discussed in class. The F-statistic for coordinate j is then

where H = X(XTX)-1XT  is the usual hat matrix (projection onto range(X)) and Hj  is the projection matrix onto range(X\j). Show that tj(2) = Fj .Hint.  Assume without loss of generality that j = d, the dth component.  (One can do so by permutating the columns of X.) Consider the QR factorization of X , i.e., X = QR where Q ∈ Rn×n  is orthogonal and R ∈ Rn×d  has the form.	

for an upper triangular (invertible) matrix T with entries Tij  = Rij  for all 1 ≤ i, j ≤ d.Question 2.3 (Non-independent noise and testing challenges):      A subtle but problematic situation occurs in linear models when noise is correlated instead of independent—indeed, this is often much  worse than non-normality of noise, which the central limit theorem more
or less addresses. To make this a bit more concrete, we consider a 2-group ANOVA model,   y1j  = μ + Q1 + ε 1j,    y2j  = μ + Q2 + ε2j ,                                 (2.1)where we assume we observe a sample of size n for each group (i.e. j = 1, . . . , n). The standard assumption is that εi  ~ N(0, σ2 In), where we use εi  = [εi1   · · ·  εin]T  ∈ Rn  for shorthand, and we have the null model
H0  : Q1  = Q2 N(0, σ2 ). In this case, for yi  =  Σj(n)=1 yij  and standard error estimate

the usual t-statistic is

the t-distribution with 2n — 2 degrees of freedom, or equivalently,
We will show that a test that re may reject
unrealistically frequently when the errors are correlated.
To that end, consider the situation that ε 1 , ε2  are independent, butεi ~ N (0, σ2 (1 — P)In + σ2P11T ) ,                                       (2.2)where P ∈ [0, 1] indicates correlation within the group.  Such correlation may be reasonable, e.g., when (hidden) confounding relates members of a group. Through the remainder of this question, let Cρ  = (1 — P)In + P11T  be a shorthand for the covariance.
(a) Show that if Z ~ N(0, Cρ), then (I — 11T)Z = Z — 1Zn  ~ N(0, (1 — P)(In  — 11T)). (b) Show that - Sn(2) ~ -  · χ2(2)n-2  under the correlation structure (2.2).
(c) Show that if Q1  = Q2  and the correlation (2.2) holds,

(d) Argue that y1  — y2  is independent of Sn(2) even with correlation (2.2). (e) Show that under correlation (2.2),

H0  too frequently when the correlation (2.2) holds.Question  2.4  (Intuition for correlated rejections via simulation):      Here we revisit ques- tion 2.3, except that we perform. some simulations and corrections.  First, we describe a general strategy for eliminating correlations in the noise, making the prima facie  ridiculous assumption that we know the noise covariance.
(a) Let y = Xβ + ε where ε ~ (0, Σ). Show that Σ-1/2y = Σ-1/2Xβ + ξ, where ξ ~ (0, In).
By part (a), if we knew Σ, we could make the substitutions
= Σ-1/2y,    X(~) = Σ-1/2X
~ 
have developed would then work.
(b) Repeat the following experiment several (say, 100) times for values of n = 2, 4, 8, 16, 32, 64, 128, 256, 512. Generate data from the ANOVA model (2.1), except that the noise is correlated (2.2) with ρ = .1, with µ = α 1   = α2   = 0.  Perform. an F-test of signiﬁcance
It is of interest to correct  an estimate for possible correlations, thereby achieving a test whose nominal level is closer to accurate. In general, one never has enough data to estimate Cov(ε) in a linear regression model except under assumptions on the noise model.  In the ANOVA model (2.1), it may be reasonable to assume that within  a group, the noises are all equally correlated, that is, the noise model (2.2) holds, and we can approximate Cov(εi). Note that E[(y1j  — y2l)2] = 2σ2  and that y1  丄 y2  under model (2.2) and the null α 1  = α2 . Deﬁne the estimates

Question  2.5 (Clumpy testing errors):     In the data ﬁle abalone .data we have data on abalone (a type of mollusc) age, where the dataset is explained in ﬁle abalone .names.  The goal is to predict the age of an abalone (given by the count of rings in its shell) from other characteristics. Here we use this dataset to investigate false discoveries and whether they come alone or in groups by adding additional complete noise variables to the data, then regressing a linear model including these noise variables.
Write code to perform. the following: ﬁrst, load the abalone data. Theni. Add two columns (call them x(1)  and x(2), say) to the data, where their entries are i.i.d.

that is, i.i.d. normal random variables with correlation P ∈ (—1, 1).
ii. Fit a linear model for the response y = Rings against all other variables (including the noise variables x(1)  and x(2))
iii. Perform. at-test for association of variable x(1)  (adjusting for all other variables) and x(2) (again, adjusting for all other variables) with y, rejecting at the level Q = .05.For the values P ∈ {— .9, — .8, — .4, 0, .4, .8, .9}, repeat the experiment instepsi–iii N = 1000 times. Across the experiments, record the number of times there is a false discovery of x(1) , a false discovery for x(2), and a false discovery of both simultaneously.  Report your false discovery counts and describe them. (Include your code in your solution.)Hints  and  pointers.   You will want to represent the abalone’s sex as a factor, that is, instead of the raw character M, F, or I  (infant), represent it in a 0-1 encoding over 3 levels. That is, if S ∈ {M, F, I} represents the sex of the abolone, transform. it into
In R this is achieved by using the method factor.  Also, the t-test in step iii is simply the standard t-test we have developed in class and is that performed by R’s summary method of a linear model.

         
加QQ：99515681  WX：codinghelp
