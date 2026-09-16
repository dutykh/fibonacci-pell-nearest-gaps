<!--
Around the Markoff uniqueness conjecture: manuscripts and reproduction
material.

Authors:
  Dr. Denys Dutykh (Mathematics Department, Khalifa University of Science
  and Technology, Abu Dhabi, UAE)
  Prof. Laurent Vuillon (Univ. Savoie Mont Blanc, CNRS, LAMA, Chambéry,
  France)
-->

# Project rules

These rules govern every manuscript under `papers/` in this repository. They
sit on top of the global rules and never override them.

## The repository URL is displayed, centred, on its own line

In the code-availability, data-availability or supplementary-material section
of every manuscript, the GitHub URL of this repository is never left inline in
the running text. The sentence that introduces it ends with a colon, and the
address is then displayed on a line of its own, centred:

```latex
\section*{Code and data availability}
The source code and the data supporting this article, with the instructions
that repeat every calculation, are contained in the supplementary material of
this article, which is available, together with that of the other articles of
this series:

{\centering
\url{https://github.com/dutykh/fibonacci-pell-nearest-gaps/}
\par}
```

Three details matter. The group `{\centering ... \par}` is used rather than
`\begin{center}` because it centres the address without adding the vertical
glue of a `center` environment, which is what keeps the block from being pushed
onto a page of its own. The blank line before it closes the introducing
sentence, so that `\centering` applies to the address alone. The address keeps
its trailing slash and carries no full stop, since a final period is easily
mistaken for part of the link.

A displayed address is one line taller than an inline one, which is enough to
leave the closing text page short and to make the strict build gate fail on an
underfull `\vbox`. Rather than shortening the sentence to fit, put

```latex
\raggedbottom
```

immediately before the first back-matter section, that is, before
`\section*{Supplementary material}` when there is one and before
`\section*{Code and data availability}` otherwise. The closing pages are then
allowed to end short instead of being stretched, the body of the article is
untouched, and the gate passes.

## Every manuscript directory answers the same five targets

The `Makefile` at the top of the repository discovers the directories under
`papers/` and forwards to each. A manuscript directory therefore provides
`all`, `check`, `checks`, `release` and `clean`, and a `help` besides. Targets
beyond those are its own business and are reached with
`make -C papers/<directory> <target>`.

## Manuscripts are recompiled after every source change

A tracked PDF that no longer matches its sources is worse than no PDF at all,
because it is the file a reader opens. Any change to a `.tex` file, a section,
a figure or a bibliography is followed by `make check` in that manuscript
directory, or `make check` at the top of the repository for all of them, before
the work is reported as done.

## Working material stays out of this repository

Plans, acceptance notes, proof maps, internal reviews, cover letters and
journal correspondence live in the research tree. What lands here is the
manuscript, the material a reader needs to check it, and the record of what
that material establishes.
