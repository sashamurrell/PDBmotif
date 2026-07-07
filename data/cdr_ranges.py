"""Chothia CDR residue-number boundaries.

Ranges are inclusive (start, end) integer residue numbers in the Chothia numbering
scheme that SAbDab applies to the variable region of antibody chains. Insertion-coded
residues (e.g. "100A") share the integer base of their preceding residue, so they fall
into the same range automatically -- no special-casing needed.

Source: the standard Chothia CDR definition table (Al-Lazikani, Lesk & Chothia 1997;
reproduced at bioinf.org.uk/abs). Validated directly against PDB 1AHW chain B (heavy):
residues 95-102 = DNSYYFDY, flanked by the conserved Cys92/Trp103 expected for CDR-H3.
"""

CDR_RANGES = {
    "h1": (26, 32),
    "h2": (52, 56),
    "h3": (93, 102),
    "l1": (24, 34),
    "l2": (50, 56),
    "l3": (89, 97),
}


def cdr_chain(cdr):
    """Return 'H' or 'L' for a CDR name like 'h1'/'H3'."""
    return cdr[0].upper()
