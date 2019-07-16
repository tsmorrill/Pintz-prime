
# An elementary bound on Siegel zeroes

The code in this repository verifies a calculation for L-functions associated to real Dirichlet characters, based on a method of Pintz. The code is written for Sage 8.4, in order to take advantage of its built in zetaderiv function.

## Getting Started

In the Sage terminal, navigate to the code folder and attach Pintz.sage.
```
sage: attach("Pintz.sage")
Siegel zeroes for primitive real characters.
```

## Running the verification

In the proof, we wish to show that F, a of the variables c, x, and q, is negative for some values of c and x. This is calculated by F(c, q0, q1, x, parity). Here, parity must be one of the strings 'even' or 'odd', which have been aliased as even and odd for convinience.

```
sage: F(1.02, 400000.000000000, 500000.000000000, 10^5.49, odd)
-1.1183346064235717
```

### Verifying Table 1

The commands to verify Table 1 are contained in data.txt in the data folder. Simply copy these into the terminal and execute.

### Generating a new table

The function cq_table(q_list) generates a table of c and x values for which F(c, q0, q1, x) is negative between the values of q given in q_list. The first block of the output is the table, formatted for the tabular environment in LaTeX. The second and third blocks are the list of Sage commands to verify the table for even and odd characters, respecitvely.


## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Thanks to Billie Thompson for the README.md template.