\\ Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science and Technology,
\\ Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambery, France)
\\
\\ Margins of Euclidean minimality for every reduced P/Q with Q <= 100 or P + Q <= 160.
\\ For each occurrence: the label m (Cohn trace of the lower Christoffel word), the occurrence pair (the
\\ column He of the central Gram factor), the complete factorization of m, the whole primitive fibre
\\ (the Gaussian sign cube), the Euclidean cost of every pair, and the margin: the cost of the cheapest
\\ other pair minus P + Q - 2, or "-" when the fibre has a single element.  The run stops on any
\\ inconsistency (an occurrence pair of the wrong norm or cost, a factor that is not a pseudoprime, a
\\ fibre of the wrong size, a margin that is not positive).
\\
\\ Usage, from the supplement directory (about ten minutes of processor time, part 1 of 1):
\\   gp -q -D parisizemax=2000000000 code/margins_census.gp
\\ or in s parallel parts, part r = 0, ..., s-1, then concatenate and sort by Q and P:
\\   echo 'PART=r; PARTS=s; read("code/margins_census.gp")' | gp -q -D parisizemax=2000000000
\\ Output: data/margins-census.txt, one line "P Q margin factorization", factors written p^e*q^f.

if(type(PART) == "t_POL", PART = 0);
if(type(PARTS) == "t_POL", PARTS = 1);
OUT = if(PARTS == 1, "data/margins-census.txt", Str("data/margins-census-part", PART, ".txt"));
MA = [2,1;1,1]; MB = [5,2;2,1]; MP1 = [1,1;1,0]; MP2 = [2,1;1,0];
word(P,Q) = vector(Q, i, if((i*P)\Q - ((i-1)*P)\Q, 2, 1));
mu(w) = {my(M = matid(2)); for(i = 1, #w, M = M*if(w[i] == 2, MB, MA)); M};
occpair(P,Q) =
{
  my(w = word(P,Q), z, n, H, v);
  z = vector(Q-2, i, w[i+1]); n = #z;
  H = mu(vector(n - (n+1)\2, i, z[(n+1)\2 + i]));
  if(n % 2, H = if(z[n\2 + 1] == 1, MP1, MP2)*H);
  v = H*[2;1];
  vecsort([v[1,1], v[2,1]]);
}
cost(x,y) = vecsum(contfrac(max(x,y)/min(x,y))) - 1;
fibre(m, f) =
{
  my(eps = 0, pr = List(), out = List(), t, z, g, s);
  for(j = 1, #f[,1],
    my(p = f[j,1], e = f[j,2]);
    if(p == 2,
      if(e > 1, error("4 divides a label"));
      eps = 1
    ,
      if(p % 4 != 1, error("a prime 3 mod 4 divides a label"));
      my(ab = qfbcornacchia(1, p)); listput(pr, [ab[1] + I*ab[2], e])
    )
  );
  t = #pr;
  for(mask = 0, 2^(t-1) - 1,
    z = if(eps, 1 + I, 1);
    for(j = 1, t,
      g = pr[j][1];
      if(j > 1 && bittest(mask, j-2), g = conj(g));
      z *= g^pr[j][2]
    );
    s = vecsort([abs(real(z)), abs(imag(z))]);
    if(s[1]^2 + s[2]^2 != m || gcd(s[1], s[2]) != 1 || s[1] == s[2], error("bad pair"));
    listput(out, s)
  );
  if(#Set(Vec(out)) != 2^(t-1), error("wrong fibre size"));
  Vec(out);
}
fstring(f) = {my(s = ""); for(j = 1, #f[,1], s = Str(s, if(j > 1, "*", ""), f[j,1], if(f[j,2] > 1, Str("^", f[j,2]), ""))); s};
{
  for(Q = 2, 159, if(Q % PARTS != PART, next);
    for(P = 1, if(Q <= 100, Q-1, min(Q-1, 160-Q)),
      if(gcd(P,Q) != 1, next);
      my(m = mu(word(P,Q))[1,2], occ = occpair(P,Q), L = P+Q-2, f, fb, others, margin);
      if(occ[1]^2 + occ[2]^2 != m || cost(occ[1], occ[2]) != L, error("occurrence pair"));
      f = factor(m);
      for(j = 1, #f[,1], if(!ispseudoprime(f[j,1]), error("composite factor")));
      fb = fibre(m, f);
      others = [cost(s[1], s[2]) | s <- fb, s != occ];
      if(#others + 1 != #fb, error("occurrence pair not in its fibre"));
      margin = if(#others, vecmin(others) - L, "-");
      if(#others && margin <= 0, error(Str("margin not positive at ", P, "/", Q)));
      write(OUT, Str(P, " ", Q, " ", margin, " ", fstring(f)))));
}
