// Complete integral-point superset for C0068's plus-branch j=1 equation.
// If y = z + 1, the two quartics have square constant terms.

SetSeed(8675309);

Qeven := [18, 72, 132, 120, 49];
Qodd := [18, 72, 84, 24, 1];

print "Qeven:";
for row in IntegralQuarticPoints(Qeven) do
    print row;
end for;

print "Qodd:";
for row in IntegralQuarticPoints(Qodd) do
    print row;
end for;
