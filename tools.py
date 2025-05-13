from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from Bio import SeqIO

class Tool(ABC):
    @abstractmethod
    def run(self, *args, **kwargs):
        pass


class Parser (Tool):
    def __init__(self, format="fasta"):
        self.format = format

    def run(self, file_path):
        self.file_path = file_path
        try:
            with open(self.file_path) as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self.file_path}")

        self.__data = {}
        self._df = pd.DataFrame()

        self.__records = list(SeqIO.parse(self.file_path, self.format))

        if not self.__records:
            raise ValueError(f"No valid records found in file: {self.file_path}")
        for record in self.__records:
            for key, value in record.__dict__.items():
                if key not in self.__data:
                    self.__data[key] = []
                if hasattr(value, '_data'):
                    value = str(value)
                self.__data[key].append(value)

        for key, values in self.__data.items():
            if values and len(values) == len(self.__records):
                self._df[key.strip('_')] = values

        if 'seq' in self._df.columns:
            self._df['length'] = self._df['seq'].str.len()
        return self._df

    def save_to_csv(self, output_path=None):
        if output_path is None:
            output_path = self.file_path.rsplit('.', 1)[0] + '.csv'
        try:
            self._df.to_csv(output_path, index=False)
            print(f"Successfully saved records to {output_path}")
        except Exception as e:
            print(f"Error saving CSV file: {str(e)}, run parser method first")


class SequenceAligner(Tool):
    def __init__(self, match=2, mismatch=-1, gap=-2, show_matrix: bool = False):
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.show_matrix = show_matrix
        self.result = {}

    def run(self, seq1, seq2, method='global'):
        if method == 'global':
            self._global_align(seq1, seq2)
        elif method == 'local':
            self._local_align(seq1, seq2)
        else:
            raise ValueError("Invalid method. Use 'global' or 'local'.")

    def _print_matrix(self, matrix):
        for row in matrix:
            print("\t".join(str(cell) for cell in row))

    def _global_align(self, seq1, seq2):
        rows, cols = len(seq1) + 1, len(seq2) + 1
        score_matrix = np.zeros((rows, cols), dtype=int)

        score_matrix[:, 0] = np.arange(0, rows) * self.gap
        score_matrix[0, :] = np.arange(0, cols) * self.gap

        for i in range(1, rows):
            for j in range(1, cols):
                match_score = score_matrix[i - 1][j - 1] + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch)
                del_score = score_matrix[i - 1][j] + self.gap
                ins_score = score_matrix[i][j - 1] + self.gap
                score_matrix[i][j] = max(match_score, del_score, ins_score)

        if self.show_matrix:
            print("Global Alignment Score Matrix:")
            self._print_matrix(score_matrix)

        aligned_seq1, aligned_seq2, matches, match_count, mismatch_count, gap_count = self._traceback(seq1, seq2, score_matrix, method="global")

        self.result = {
            "type": "global",
            "score": score_matrix[-1][-1],
            "aligned_seq1": aligned_seq1,
            "aligned_seq2": aligned_seq2,
            "matches": matches,
            "match_count": match_count,
            "mismatch_count": mismatch_count,
            "gap_count": gap_count
        }

    def _local_align(self, seq1, seq2):
        rows, cols = len(seq1) + 1, len(seq2) + 1
        score_matrix = np.zeros((rows, cols), dtype=int)

        max_score, max_pos = 0, (0, 0)

        for i in range(1, rows):
            for j in range(1, cols):
                match_score = score_matrix[i - 1][j - 1] + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch)
                del_score = score_matrix[i - 1][j] + self.gap
                ins_score = score_matrix[i][j - 1] + self.gap
                score_matrix[i][j] = max(0, match_score, del_score, ins_score)

                if score_matrix[i][j] > max_score:
                    max_score = score_matrix[i][j]
                    max_pos = (i, j)

        if self.show_matrix:
            print("Local Alignment Score Matrix:")
            self._print_matrix(score_matrix)

        aligned_seq1, aligned_seq2, matches, match_count, mismatch_count, gap_count = self._traceback(seq1, seq2, score_matrix, max_pos, method="local")

        self.result = {
            "type": "local",
            "score": max_score,
            "aligned_seq1": aligned_seq1,
            "aligned_seq2": aligned_seq2,
            "matches": matches,
            "match_count": match_count,
            "mismatch_count": mismatch_count,
            "gap_count": gap_count
        }

    def _traceback(self, seq1, seq2, score_matrix, max_pos=None, method="global"):
        aligned_seq1, aligned_seq2, matches = [], [], []
        match_count = mismatch_count = gap_count = 0  # Initialize counts
        if method == "global":
            i, j = len(seq1), len(seq2)
        else:  # Local
            i, j = max_pos

        while i > 0 and j > 0 and (score_matrix[i][j] > 0 or method == "global"):
            current = score_matrix[i][j]
            diag = score_matrix[i - 1][j - 1]
            up = score_matrix[i][j - 1]
            left = score_matrix[i - 1][j]

            if current == diag + (self.match if seq1[i - 1] == seq2[j - 1] else self.mismatch):
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append(seq2[j - 1])
                if seq1[i - 1] == seq2[j - 1]:
                    matches.append('|')
                    match_count += 1  # Increment match count
                else:
                    matches.append(' ')
                    mismatch_count += 1  # Increment mismatch count
                i -= 1
                j -= 1
            elif current == up + self.gap:
                aligned_seq1.append('-')
                aligned_seq2.append(seq2[j - 1])
                matches.append(' ')
                gap_count += 1  # Increment gap count
                j -= 1
            else:  # left
                aligned_seq1.append(seq1[i - 1])
                aligned_seq2.append('-')
                matches.append(' ')
                gap_count += 1  # Increment gap count
                i -= 1

        while i > 0:
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append('-')
            matches.append(' ')
            gap_count += 1  # Increment gap count
            i -= 1
        while j > 0:
            aligned_seq1.append('-')
            aligned_seq2.append(seq2[j - 1])
            matches.append(' ')
            gap_count += 1  # Increment gap count
            j -= 1

        aligned_seq1.reverse()
        aligned_seq2.reverse()
        matches.reverse()

        return aligned_seq1, aligned_seq2, matches, match_count, mismatch_count, gap_count

    def report(self, width=50, print_alignment=True):
        if not self.result:
            print("No alignment has been run.")
            return

        print(f"{self.result['type'].capitalize()} alignment:")
        print(f"Score: {self.result['score']}")
        print(f"Matches: {self.result['match_count']}")
        print(f"Mismatches: {self.result['mismatch_count']}")
        print(f"Gaps: {self.result['gap_count']}\n")

        a1 = self.result['aligned_seq1']
        a2 = self.result['aligned_seq2']
        m = self.result['matches']

        if print_alignment:
            for i in range(0, len(a1), width):
                print(' '.join(a1[i:i + width]))
                print(' '.join(m[i:i + width]))
                print(' '.join(a2[i:i + width]))
                print()

    def get_alignment_data(self):
        if not self.result:
            raise ValueError("No alignment has been run.")
        return self.result


class MotifFinder(Tool):
    def __init__(self, motif: str):
        self._motif = motif.upper()
        self.__last_result = None

    def run(self, sequence: str):
        sequence = sequence.upper()
        positions = self._find_motif_occurrences(sequence)
        self.__last_result = {
            "motif": self._motif,
            "count": len(positions),
            "positions": positions
        }
        return self.__last_result

    def _find_motif_occurrences(self, sequence: str):
        positions = []
        motif_len = len(self._motif)
        for i in range(len(sequence) - motif_len + 1):
            if sequence[i:i + motif_len] == self._motif:
                positions.append(i)
        return positions

    def get_result(self):
        if not self.__last_result:
            return {"error": "No motif search has been run yet."}
        return self.__last_result

    def report(self):
        if not self.__last_result:
            print("No motif search has been run yet.")
            return
        print(f"Motif: {self.__last_result['motif']}")
        print(f"Count: {self.__last_result['count']}")
        print(f"Positions: {', '.join(map(str, self.__last_result['positions']))}")



if __name__ == '__main__':
    from sequence import *

    parser = Parser('fasta')
    records = parser.run("synthetic_mtDNA_dataset.fasta")
    NC10seq = MitochondrialDNA(records.loc[10]).sequence
    NC11seq = MitochondrialDNA(records.loc[11]).sequence

    aligner = SequenceAligner()
    aligner.run(NC10seq, NC11seq)
    aligner.report(print_alignment=False)
    #print (aligner.get_alignment_data())
    aligner.run(NC10seq, NC11seq, method="local")
    aligner.report(print_alignment=False)

    motif = MotifFinder("GATC")
    motif.run(NC10seq)
    motif.report()
    #print(motif.get_result())
