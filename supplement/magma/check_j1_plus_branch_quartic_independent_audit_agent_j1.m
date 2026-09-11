Qplus := [18, 72, 132, 120, 49];
Qminus := [18, 72, 84, 24, 1];

major, minor, patch := GetVersion();
printf "VERSION V%o.%o-%o\n", major, minor, patch;

SetSeed(42424243);
seed_plus, step_plus_before := GetSeed();
A := IntegralQuarticPoints(Qplus);
seed_plus_after, step_plus_after := GetSeed();

SetSeed(98765431);
seed_minus, step_minus_before := GetSeed();
B := IntegralQuarticPoints(Qminus);
seed_minus_after, step_minus_after := GetSeed();

Anorm := Sort([[row[1], Abs(row[2])] : row in A]);
Bnorm := Sort([[row[1], Abs(row[2])] : row in B]);
expected_plus := [
    [-4, 41],
    [-2, 7],
    [0, 7],
    [2, 41]
];
expected_minus := [
    [-2, 1],
    [0, 1]
];

assert Anorm eq expected_plus;
assert Bnorm eq expected_minus;

print "PLUS_SEED", seed_plus, step_plus_before, step_plus_after;
print "Qplus_normalized", Anorm;
print "MINUS_SEED", seed_minus, step_minus_before, step_minus_after;
print "Qminus_normalized", Bnorm;

Qx<x> := PolynomialRing(Rationals());
fplus := 18*x^4 + 72*x^3 + 132*x^2 + 120*x + 49;
fminus := 18*x^4 + 72*x^3 + 84*x^2 + 24*x + 1;
assert Discriminant(fplus) ne 0;
assert Discriminant(fminus) ne 0;

Eplus, map_plus := AssociatedEllipticCurve(fplus);
Eminus, map_minus := AssociatedEllipticCurve(fminus);

rank_plus, generators_plus, sha_plus := MordellWeilShaInformation(Eplus);
Gplus, mw_plus, rank_plus_proved, full_plus_proved :=
    MordellWeilGroup(Eplus);
assert rank_plus_proved;
assert full_plus_proved;

rank_minus, generators_minus, sha_minus := MordellWeilShaInformation(Eminus);
Gminus, mw_minus, rank_minus_proved, full_minus_proved :=
    MordellWeilGroup(Eminus);
assert rank_minus_proved;
assert full_minus_proved;

print "DISCRIMINANTS", Discriminant(fplus), Discriminant(fminus);
print "ASSOCIATED_PLUS", Eplus;
print "PLUS_RANK_INFO", rank_plus;
print "PLUS_MW_STATUS", rank_plus_proved, full_plus_proved;
print "ASSOCIATED_MINUS", Eminus;
print "MINUS_RANK_INFO", rank_minus;
print "MINUS_MW_STATUS", rank_minus_proved, full_minus_proved;
print "PASS: both complete shifted-quartic lists and Mordell-Weil statuses were independently reproduced.";
