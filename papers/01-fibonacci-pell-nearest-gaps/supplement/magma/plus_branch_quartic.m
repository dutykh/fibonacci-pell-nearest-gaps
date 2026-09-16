// Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
// and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
// Blanc, CNRS, LAMA, Chambery, France).
// Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
// Distributed under the GNU Lesser General Public License, version 2.1;
// see the LICENSE file of this package.
// Complete integral-point superset for the nearest-gap criterion plus-branch j=1 equation.
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
