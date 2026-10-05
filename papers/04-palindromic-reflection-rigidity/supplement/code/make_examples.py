#!/usr/bin/env python3
# Authors: Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
# and Technology, Abu Dhabi, UAE) and Prof. Laurent Vuillon (Univ. Savoie Mont
# Blanc, CNRS, LAMA, Chambery, France).
# Copyright (C) 2026 Dr. Denys Dutykh and Prof. Laurent Vuillon.
# Distributed under the GNU Lesser General Public License, version 2.1;
# see the LICENSE file of this package.
"""Regenerate the five prescribed table rows; no enumeration is performed."""
import json
import sys
from pathlib import Path
from palindromic_words import (pal, paired_preimages, standard_image, gram_data,
                               completed_counts)

if not __debug__ or len(sys.argv) != 1:
    raise SystemExit("Run without -O; this fixed generator accepts no arguments.")

rows = []
for directive in ("", "a", "b", "ab", "ba"):
    words = [standard_image(directive, z) for z in paired_preimages("ab")]
    rows.append({"directive": directive, "ancestor": pal(directive),
                 "ancestor_gram": gram_data(pal(directive)),
                 "images": [{"word": z, "gram": gram_data(z),
                             "completed_counts": completed_counts(z)} for z in words]})
path = Path(__file__).resolve().parents[1] / "data" / "midpoint_examples.json"
path.write_text(json.dumps(rows, indent=2) + "\n")
print("Wrote five fixed exact examples to " + str(path))
