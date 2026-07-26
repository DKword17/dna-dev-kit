#!/usr/bin/env python3
"""
dna-dev-kit/README.md

DNA Development Kit — Synthetic Biology Design Toolkit
=======================================================

A comprehensive toolkit for DNA sequence design, manipulation, and
assembly for synthetic biology applications.

Features:
    - DNA sequence I/O (FASTA, GenBank)
    - Restriction enzyme digestion (100+ Type II and Type IIS)
    - PCR primer design (nearest-neighbor thermodynamics)
    - Codon optimization for 12 expression hosts (E. coli, yeast, human, etc.)
    - DNA assembly planning (Golden Gate, Gibson, CPEC)
    - Sequence translation (standard genetic code)
    - Melting temperature calculation (SantaLucia 1998)
    - Reverse complement, GC content, restriction site mapping

Architecture:
    genome_engine.py
        ├── DNASequence          (core sequence object)
        ├── codon_optimize       (host-specific optimization)
        ├── design_primers       (Primer3-style primer design)
        ├── plan_golden_gate_assembly (multipart assembly)
        ├── RestrictionEnzyme    (database)
        └── AssemblyPart         (DNA assembly part)

References:
    - Gibson et al. (2009) Nat Methods 6:343 — Gibson assembly
    - Engler et al. (2008) PLoS ONE 3:e3647 — Golden Gate
    - SantaLucia (1998) PNAS 95:1460 — DNA thermodynamics
    - Sambrook & Russell (2001) Molecular Cloning

Example:
    >>> from genome_engine import DNASequence, codon_optimize
    >>> seq = DNASequence("ATGCGTACGTCG")
    >>> print(seq.gc_content(), seq.translate())

Quick Start:
    $ python genome_engine.py
