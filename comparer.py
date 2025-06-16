# part3 
from sequence import MitochondrialDNA
from tools import SequenceAligner, MotifFinder
from typing import List


class SequenceAlignWrapper:
    def __init__(self):
        self.aligner = SequenceAligner()

    def align(self, seq1: str, seq2: str, method: str = 'global'):
        self.aligner.run(seq1, seq2, method=method)
        return self.aligner.get_alignment_data()

    def report(self):
        self.aligner.report()


class AlignmentVisualizer:
    def __init__(self, aligner: SequenceAlignWrapper):
        self.aligner = aligner

    def display(self, idx1: int, idx2: int, sequences: List[MitochondrialDNA], method: str = 'global', width: int = 80):
        result = self.aligner.align(sequences[idx1].sequence, sequences[idx2].sequence, method=method)
        a1, a2, m = result['aligned_seq1'], result['aligned_seq2'], result['matches']

        print(f"Alignment between sequence {idx1} and {idx2} ({method}):\n")
        for i in range(0, len(a1), width):
            print(''.join(a1[i:i + width]))
            print(''.join(m[i:i + width]))
            print(''.join(a2[i:i + width]))
            print()

        print(f"Score: {result['score']}")
        print(f"Matches: {result['match_count']}, Mismatches: {result['mismatch_count']}, Gaps: {result['gap_count']}")


class SequenceComparer:
    def __init__(self, sequences: List[MitochondrialDNA]):
        self.sequences = sequences
        self.wrapper = SequenceAlignWrapper()

    def compare_pair(self, idx1: int, idx2: int, method: str = 'global'):
        return self.wrapper.align(self.sequences[idx1].sequence, self.sequences[idx2].sequence, method)

    def compare_all(self):
        results = []
        for i in range(len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                data = self.compare_pair(i, j)
                results.append({
                    'pair': (i, j),
                    'score': data['score'],
                    'matches': data['match_count'],
                    'mismatches': data['mismatch_count'],
                    'gaps': data['gap_count']
                })
        return results

    def compare_to_reference(self, ref_index: int = 0, method: str = 'global'):
        results = []
        for i in range(len(self.sequences)):
            if i == ref_index:
                continue
            data = self.compare_pair(ref_index, i, method=method)
            results.append({
                'reference_vs': i,
                'score': data['score'],
                'matches': data['match_count'],
                'mismatches': data['mismatch_count'],
                'gaps': data['gap_count']
            })
        return results


class ConservedMotifAnalyzer:
    def __init__(self, sequences: List[MitochondrialDNA]):
        self.sequences = sequences
        self.finder = MotifFinder()

    def find_conserved(self, k=5, threshold=2):
        motif_counts: dict[str, int] = {}
        for seq in self.sequences:
            result = self.finder.run(seq.sequence, k=k, threshold=threshold)
            for motif in result:
                motif_counts[motif] = motif_counts.get(motif, 0) + 1

        conserved = {m: c for m, c in motif_counts.items() if c == len(self.sequences)}
        return conserved

    def motif_occurrences(self, motif: str):
        positions = {}
        for i, seq in enumerate(self.sequences):
            result = self.finder.run(seq.sequence, motif=motif)
            positions[i] = result.get('positions', [])
        return positions


class MultiAligner:
    def __init__(self, sequences: List[MitochondrialDNA]):
        self.sequences = sequences
        self.wrapper = SequenceAlignWrapper()

    def pairwise_scores(self, method='global'):
        score_map = {}
        for i in range(len(self.sequences)):
            for j in range(i + 1, len(self.sequences)):
                result = self.wrapper.align(self.sequences[i].sequence, self.sequences[j].sequence, method)
                score_map[(i, j)] = result['score']
        return score_map

    def align_pair(self, idx1: int, idx2: int, method='global'):
        return self.wrapper.align(self.sequences[idx1].sequence, self.sequences[idx2].sequence, method)

    def report_alignment(self, idx1: int, idx2: int, method='global'):
        self.wrapper.align(self.sequences[idx1].sequence, self.sequences[idx2].sequence, method)
        self.wrapper.report()


# test usage
if __name__ == "__main__":
    from tools import Parser

    parser = Parser()
    df = parser.run("synthetic_mtDNA_dataset.fasta.txt")
    mito_seqs = [MitochondrialDNA(df.loc[i]) for i in range(3)]  # using first 3 sequences

    comparer = SequenceComparer(mito_seqs)
    visualizer = AlignmentVisualizer(comparer.wrapper)

    print("Pairwise Comparison Result (0 vs 1):")
    print(comparer.compare_pair(0, 1))

    print("\nAll Pairwise Comparisons:")
    for result in comparer.compare_all():
        print(result)

    print("\nComparison to Reference (index 0):")
    for result in comparer.compare_to_reference(ref_index=0):
        print(result)

    print("\nVisualizing alignment between 0 and 1:")
    visualizer.display(0, 1, mito_seqs, method='global', width=60)
