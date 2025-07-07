from abc import ABC, abstractmethod
import numpy as np
import pandas as pd
from Bio import SeqIO
from collections import Counter
from sequence import MitochondrialDNA

class Tool(ABC):
    @abstractmethod
    def run(self, *args, **kwargs):
        '''
        abstract method for running the tool
        '''
        pass

    @abstractmethod
    def report(self, *args, **kwargs):
        '''
        abstract method for graphical representation
        '''
        pass

class Parser (Tool):
    def __init__(self, format="fasta"):
        '''
        :param format: fasta format for now, can be anything SeqIO supports
        '''
        self.format = format
        self.__data = {}
        self._df = pd.DataFrame()
        self.__records = None

    def run(self, file_path, return_objects=False):
        '''
        :param file_path: path to the file to parse
        :param return_objects: if True, returns a list of MitochondrialDNA objects
        :return: a pandas DataFrame or a list of MitochondrialDNA objects
        '''
        self._file_path = file_path
        try:
            with open(self._file_path) as f:
                pass
        except FileNotFoundError:
            raise FileNotFoundError(f"File not found: {self._file_path}")

        self.__data = {}
        self._df = pd.DataFrame()

        self.__records = list(SeqIO.parse(self._file_path, self.format))

        if not self.__records:
            raise ValueError(f"No valid records found in file: {self._file_path}")
        for record in self.__records:
            for key, value in record.__dict__.items():
                if key not in self.__data:
                    self.__data[key] = []
                if hasattr(value, '_data'):
                    value = str(value)
                self.__data[key].append(value)

        for key, values in self.__data.items():
            self._df[key.strip('_')] = values

        if 'seq' in self._df.columns:
            self._df['length'] = self._df['seq'].str.len()

        if return_objects:
            return [MitochondrialDNA(self._df.loc[i]) for i in range(len(self._df))]
        else:
            return self._df

    def save_to_csv(self, output_path=None):
        if output_path is None:
            output_path = self._file_path.rsplit('.', 1)[0] + '.csv'
        try:
            self._df.to_csv(output_path, index=False)
            print(f"Successfully saved records to {output_path} \n")
        except Exception as e:
            print(f"Error saving CSV file: {str(e)}, run parser method first \n")

    def report(self, print_header: bool = True):
        print(f'Successfully parsed {self._file_path}')
        print(f'Found {len(self._df)} records \n')
        if print_header:
            print(f'{self._df.head()} \n')


class SequenceAligner(Tool):
    def __init__(self, match=2, mismatch=-1, gap=-2, show_matrix: bool = False):
        '''
        :param match: score for match
        :param mismatch: score for mismatch
        :param gap: score for gap
        :param show_matrix: boolean, do we want to print the score matrix
        '''
        self.match = match
        self.mismatch = mismatch
        self.gap = gap
        self.show_matrix = show_matrix
        self.result = {}

    def run(self, seq1, seq2, method='global'):
        '''
        :param seq1: sequence to align
        :param seq2: sequence to align
        :param method: global or local alignment, default is global
        '''
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
        '''
        :param seq1: sequence to align
        :param seq2: sequence to align
        Updates:
        self.result (dict): Dictionary containing:
            - 'type': Alignment type ("global")
            - 'score': Final alignment score
            - 'aligned_seq1': Aligned version of sequence 1 (with gaps)
            - 'aligned_seq2': Aligned version of sequence 2 (with gaps)
            - 'matches': List of match/mismatch indicators ('|' or ' ')
            - 'match_count': Number of exact matches
            - 'mismatch_count': Number of mismatches
            - 'gap_count': Number of gaps introduced in either sequence
        '''
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
        '''
        :param seq1: sequence to align
        :param seq2: sequence to align
        Updates:
        self.result (dict): Dictionary containing:
            - 'type': Alignment type ("global")
            - 'score': Final alignment score
            - 'aligned_seq1': Aligned version of sequence 1 (with gaps)
            - 'aligned_seq2': Aligned version of sequence 2 (with gaps)
            - 'matches': List of match/mismatch indicators ('|' or ' ')
            - 'match_count': Number of exact matches
            - 'mismatch_count': Number of mismatches
            - 'gap_count': Number of gaps introduced in either sequence
        '''
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
        '''Trace back through the score matrix to reconstruct the optimal alignment.
        :param seq1: sequence 1
        :param seq2: sequence 2
        :param score_matrix: score matrix from the alignment
        :param max_pos: the maximum position of the local alignment
        :param method: global or local alignment, default is global
        '''
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
        '''Graphical report of the statistics and the alignment to the given width.
        :param width: the width of the alignment
        :param print_alignment: boolean, print the graphical representation of the alignment
        '''
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
        '''
        :return: dictionary of the entire alignment data
        '''
        if not self.result:
            raise ValueError("No alignment has been run.")
        return self.result


from typing import List
from sequence import MitochondrialDNA

class MotifFinder(Tool):
    def __init__(self):
        self.__last_result = {}
        self.__discovered = False

    def run(self, sequences: List[MitochondrialDNA], motif: str = None, k: int = 5, threshold: int = 2):
        '''
        :param sequences: A list of MitochondrialDNA objects to analyze.
        :param motif: Specific motif to find across all sequences. If None, discovers conserved motifs.
        :param k: Length of k-mers for motif discovery.
        :param threshold: Minimum number of sequences a motif must appear in to be considered conserved.
        :return: dict containing results.
        '''
        if not sequences:
            return {'error': 'No sequences provided.'}

        if motif:
            self.__last_result = self._find_specific_motif_across_sequences(sequences, motif)
            self.__discovered = False
        else:
            self.__last_result = self._discover_conserved_motifs(sequences, k, threshold)
            self.__discovered = True
            if not self.__last_result:
                print("No conserved motifs found.")
        return self.__last_result

    def _find_specific_motif_across_sequences(self, sequences: List[MitochondrialDNA], motif: str):
        results = []
        for i, seq_obj in enumerate(sequences):
            positions = self._find_motif_occurrences(seq_obj.sequence.upper(), motif.upper())
            if positions:
                results.append({
                    'sequence_index': i,
                    'sequence_name': seq_obj.name,
                    'motif': motif,
                    'count': len(positions),
                    'positions': positions
                })
        return results

    def _find_motif_occurrences(self, sequence: str, motif: str):
        positions = []
        motif_len = len(motif)
        for i in range(len(sequence) - motif_len + 1):
            if sequence[i:i + motif_len] == motif:
                positions.append(i)
        return positions

    def _discover_conserved_motifs(self, sequences: List[MitochondrialDNA], k: int, threshold: int):
        motif_occurrences_details = {}
        for seq_idx, seq_obj in enumerate(sequences):
            kmer_counts = Counter()
            sequence = seq_obj.sequence.upper()
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i + k]
                kmer_counts[kmer] += 1
    
            found_motifs_in_seq = {motif: count for motif, count in kmer_counts.items() if count >= 1}
    
            for motif in found_motifs_in_seq:
                if motif not in motif_occurrences_details:
                    motif_occurrences_details[motif] = {}
                motif_occurrences_details[motif][seq_idx] = self._find_motif_occurrences(sequence, motif)
    
    
        conserved_motifs = []
        for motif, seq_details in motif_occurrences_details.items():
            if len(seq_details) >= threshold:
                conserved_motifs.append({
                    'motif': motif,
                    'sequences': [
                        {
                            'sequence_index': i,
                            'sequence_name': sequences[i].name,
                            'positions': positions
                        }
                        for i, positions in seq_details.items()
                    ]
                })
    
        return conserved_motifs

    def get_result(self):
        if not self.__last_result:
            return {"error": "No motif search has been run yet."}
        return self.__last_result

    def report(self):
        if not self.__last_result:
            print("No motif search has been run yet.")
            return
        if self.__discovered:
            print(f'Discovered the following motifs:')
            for motif, seq_details in self.__last_result.items():
                print(f'\t{motif}:')
                for seq_idx, positions in seq_details.items():
                    print(f'\t\tSequence {seq_idx}: {len(positions)} occurrences at {positions}')
            print()
        else:
            print(f"Specific Motif Search Results:")
            if not self.__last_result:
                print("\tMotif not found in any sequence.")
            for result in self.__last_result:
                print(f"\tSequence {result['sequence_index']} ({result['sequence_name']}):")
                print(f"\t\tMotif: {result['motif']}")
                print(f"\t\tCount: {result['count']}")
                print(f"\t\tPositions: {', '.join(map(str, result['positions']))}")
            print()
