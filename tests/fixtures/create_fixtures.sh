#!/usr/bin/env bash
# Generates synthetic FASTA and FASTQ fixtures for integration tests.
# Sequences are non-repetitive and non-palindromic so reads align unambiguously.
# R2 in the FASTQ is the Illumina-convention reverse-complement of the reference
# at the R2 insert position.
#
# Run once from the project root: bash tests/fixtures/create_fixtures.sh

set -euo pipefail
FIXTURES_DIR="$(cd "$(dirname "$0")" && pwd)"

# ── Reference FASTA (2 sequences, 200 bp each) ────────────────────────────────
# Deliberately non-repetitive, non-palindromic sequences.
cat > "$FIXTURES_DIR/ref.fasta" <<'EOF'
>ref_A
AGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGATCGATCGAATCGATCGATCGAATCGATC
GATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCG
ATCGATCGAATCGATCGATCGAAT
>ref_B
TTAGCCTTAGCCTTAGCCGAATCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCT
TAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGC
CTTAGCCTTAGCCTTAGCC
EOF

python3 - "$FIXTURES_DIR" << 'PYEOF'
import sys, os

fixtures = sys.argv[1]

def rc(s):
    comp = str.maketrans('ACGTacgt', 'TGCAtgca')
    return s.translate(comp)[::-1]

# Use a long single string without newlines
ref_a = ("AGCTAGCTAGCTAGCTTACGATCGATCGAATCGAATCGATCGATCGATCGATCGATCGAATCGATCGATCGAATCGATC"
         "GATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCGATCGATCGAATCG"
         "ATCGATCGAATCGATCGATCGAAT")

ref_b = ("TTAGCCTTAGCCTTAGCCGAATCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCT"
         "TAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGCCTTAGC"
         "CTTAGCCTTAGCCTTAGCC")

# Remove newlines (they're formatting artifacts from the heredoc above)
ref_a = ref_a.replace('\n', '').replace(' ', '')
ref_b = ref_b.replace('\n', '').replace(' ', '')

READ_LEN = 50

# Build read pairs: R1=forward strand at pos O, R2=rc(forward at O+READ_LEN)
def make_pair(ref, o):
    r1 = ref[o : o + READ_LEN]
    r2_fwd = ref[o + READ_LEN : o + 2 * READ_LEN]
    r2_fastq = rc(r2_fwd)  # illumina R2 is rc of forward insert
    return r1, r2_fastq

pairs_a = [make_pair(ref_a, o) for o in [0, 10, 20]]
pairs_b = [make_pair(ref_b, o) for o in [0, 10]]

all_pairs = pairs_a + pairs_b
qual = 'I' * READ_LEN

r1_lines = []
r2_lines = []
for i, (r1, r2) in enumerate(all_pairs, 1):
    r1_lines += [f'@read{i}/1', r1, '+', qual]
    r2_lines += [f'@read{i}/2', r2, '+', qual]

with open(os.path.join(fixtures, 'r1.fastq'), 'w') as f:
    f.write('\n'.join(r1_lines) + '\n')

with open(os.path.join(fixtures, 'r2.fastq'), 'w') as f:
    f.write('\n'.join(r2_lines) + '\n')

print("Read pairs generated:")
for i, (r1, r2) in enumerate(all_pairs, 1):
    print(f"  read{i}: R1={r1} R2={r2} (R2 in FASTQ is rc of ref at pos {([0,10,20]+[0,10])[i-1]+50})")
    assert len(r1) == READ_LEN, f"R1 too short: {len(r1)}"
    assert len(r2) == READ_LEN, f"R2 too short: {len(r2)}"

print("Done.")
PYEOF

echo "Fixtures written to $FIXTURES_DIR"
