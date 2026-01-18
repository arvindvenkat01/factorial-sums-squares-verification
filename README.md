# Factorial Sums Equal to Perfect Squares: Computational Verification

Exhaustive verification of OEIS sequence [A014597](https://oeis.org/A014597): 
Numbers k such that k² is a sum of distinct factorials.

## Problem Statement

Find all integers k where:
- **General case**: k² = Σ(distinct factorials)
- **Variant**: k² = Σ(distinct factorials) + 1

## Computational Results

Exhaustively searched all 2^40 subsets of {1!, 2!, ..., 40!}:

- **Search space**: 1,099,511,627,776 combinations
- **Nodes evaluated**: 2,201,375,574 
- **Computation time**: ~84 minutes
- **Solutions found**: 
  - General case: 15 (all known in A014597)
  - +1 variant: 13

**Conclusion**: No new solutions exist up to N=40.

## Algorithm

Uses advanced pruning strategies:
- Quadratic residue filtering (95 primes)
- Binary tree traversal with early termination
- Tree/Leaf prime splitting
- Trigger-point checking

See [docs/algorithm_explanation.md](docs/algorithm_explanation.md) for details.

## Usage


### Requires Python 3.8+ and gmpy2
```bash
pip install gmpy2
```

### Run general case (A014597)
```bash
python src/a014597_general.py
```

### Run +1 variant
```
python src/plus_one_variant.py
```

## Results
General Case (A014597)
All 15 known solutions found:

```bash 
1² = 1!
3² = 1! + 2! + 3!
5² = 1! + 4!
11² = 1! + 5!
12² = 4! + 5!
...
1183893² = 1! + 2! + 3! + 7! + 8! + 9! + 10! + 11! + 12! + 13! + 14! + 15!
```
See results/N40_general_output.txt for full output.


## Related Work
This verification extends the computational frontier from N≈20-21 to N=40.


Citation
text
@misc{venkat2026factorial,
  author = {Venkat, Arvind N.},
  title = {Computational Verification of OEIS A014597 to N=40},
  year = {2026},
  url = {https://github.com/arvindvenkat01/factorial-sums-squares-verification}
}

License
MIT License

Contact
Arvind N. Venkat - arvind.venkat01@gmail.com
