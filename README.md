# How many solutions are there?

## Overview

 My goal with this project is in the repo description, but to state it more carefully:

 For a fixed prime q, when is there more than one solution for a prime p = 2^m + q^n, for m, n >= 1. The conjecture here is that (p,q) = (11,3) is the only pair that survives.

 Perhaps the most significant value is A, which is the co-primality between the two partitions. See the roadmap or concise notes for a quick primer. Note the A = 1 regime is quickly dispensed with and has a parallel to Mihăilescu's theorem. Clearly, the A = 1 having just the one example doesn't need that result to be established. But the relevance of the expression is worth noting before moving to the A > 1 case, which is where we spend the preponderance of our time.

 The exponent differences are typically labeled M and N in the repo. Often you'll find r = ord_q(2), and E = lcm(2, r), and various labels for the exponents vary based on context. In time this point will get some clean up.

## Repository contents

```text
.
├── src/tex/                  roadmap, concise notes, long notes and BibTeX
├── docs/
│   ├── notes/                compiled notes PDFs
│   │   └── old/              original TeX manuscripts
│   ├── pf/                   reserved for the proof
│   └── ref/                  local reference copies
├── lean/DiophantineBinaryExp/ Lean project (Git submodule)
├── scripts/                  Sage, PARI/GP and native experiments
│   ├── exp-gcd/              exponent-gap checks and multiplier search
│   └── mers-tab/             GAP calculations on Mersenne factors
├── labs/                     interactive Sage notebook
├── data/                     local experiment outputs
├── build/                    generated TeX files and native binaries
├── makefile                  document and native build targets
└── pixi.toml                 Sage and Jupyter environment
```

- Start with the [proof roadmap](docs/notes/prime_power_roadmap.pdf). The [concise notes](docs/notes/prime_power_notes_concise.pdf) give the arithmetic reductions, and the [long notes](docs/notes/prime_power_notes.pdf) include the supporting calculations and alternative arguments.
- The Lean project maps the written argument to prospective Lean statements and library dependencies. Its outline and proof plan are in `lean/DiophantineBinaryExp/README.md`. [Lean project](https://github.com/genamethi/DiophantineBinaryExp#readme)
- The [notebook](labs/prime-power.ipynb) collects interactive experiments. The [q-comparison guide](scripts/qcomp_readme.md) documents the scans and their outputs.

## Usage

| Command | Output |
|---|---|
| `make -j` | All three notes PDFs and native tools |
| `make notes` | PDFs in `docs/notes/`, with LaTeX files in `build/tex/` |
| `make prime_power_roadmap` | The roadmap PDF only. Use `prime_power_notes` or `prime_power_notes_concise` for the other PDFs. |
| `make native` | The three executables in `build/bin/` |
| `make scan_q_native` | Just the scanner. The other tool targets are `multiplier_orders` and `explore_2adic_lifts`. |
| `make display` | Open the main notes with `PDFVIEWER` (defaults to `pdfterm`) |
| `make clean` | Remove `build/` and the three generated notes PDFs |

The notes use latexmk, pdfLaTeX, and BibTeX. The native tools use
C17 or C++17, with FLINT and GMP for `scan_q_native`, POSIX threads for
`scan_q_native` and `multiplier_orders`, and OpenMP and GMP for
`explore_2adic_lifts`. Run them from the repository root, for example
`./build/bin/scan_q_native --help`.

## AI Assistance

These readme comments I've written myself. A lot of the drafting, type-setting and various stages of brainstorming and derivation have been done by some combination of Gemini 3.1 Pro, Claude Opus and Fable (Ranging from 4.5 to 4.6 for the former. 5.0 for the latter.), and more recently some combination of Astra, Sol and Terra (6.2, 6.1, 6.1).

The ideas I would say are mostly my own in terms of research direction, however, Astra and Sol have been helpful for specific insights and the discovery of better bounds. A few topics, such as the more abstract parts of the power-residue reciprocity argument in section 8 are a bit outside of my current knowledge. My goal is to rely on these arguments as little as possible, while extending that knowledge.1

## About Me

I figure since the AI Assistance part of thing tends to go hand in hand with credibility questions, it's worth asking what my backgrond is.
In short I have a BS in Neuroscience and Cognitive Science (Emphasis in Language Sciences) & Mathematics (Comprehensive Emphasis) from University of Arizona (2019). While I've done some graduate level coursework in algebra, my knowledge of algebraic number theory largely comes from projects such as these. This problem stems from various inquiries into this particular expression, which you'll find littered in various states of disarray in my repos.

I likely will just keep to amateur mathematics, but at some point may extend my education as I do have the intent to get a PhD in Linguistics.

## Related Repos

This particular question slowly got my attention while studying the partitions of prime powers, length two while building https://github.com/genamethi/dioph.pp That application is a more open ended tool and attempt at trying different data science tools out with experimental mathematics.
