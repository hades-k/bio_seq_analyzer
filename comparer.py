from sequence import MitochondrialDNA
from tools import SequenceAligner, MotifFinder
from typing import List
import multiprocessing

class SequenceAlignWrapper:
    def __init__(self, match: int = 2, mismatch: int = -1, gap: int = -2):
        self.aligner = SequenceAligner(match=match, mismatch=mismatch, gap=gap)

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
    def __init__(self, sequences: List[MitochondrialDNA], match: int = 2, mismatch: int = -1, gap: int = -2):
        self.sequences = sequences
        self.wrapper = SequenceAlignWrapper(match=match, mismatch=mismatch, gap=gap)

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
