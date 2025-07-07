from dataclasses import dataclass

@dataclass
class AlignmentWindowResult:
    name1: str
    name2: str
    score: int
    matches: int
    mismatches: int
    gaps: int
    window_text: str
    max_slider: int

class AlignmentWindow:
    def __init__(self, seq1_name, seq2_name, aligned_seq1, aligned_seq2, matches, score, match_count, mismatch_count, gap_count):
        self.name1 = seq1_name
        self.name2 = seq2_name
        self.seq1 = aligned_seq1
        self.seq2 = aligned_seq2
        self.match_line = matches
        self.score = score
        self.match_count = match_count
        self.mismatch_count = mismatch_count
        self.gap_count = gap_count

    def get_window(self, slider: int, window_size: int = 80) -> AlignmentWindowResult:
        start = slider * window_size
        end = start + window_size
        aligned_a = ''.join(self.seq1[start:end])
        aligned_b = ''.join(self.seq2[start:end])
        match_line = ''.join(self.match_line[start:end])
        combined = f"{self.name1}:\n{aligned_a}\n{match_line}\n{aligned_b}\n{self.name2}"
        max_slider = len(self.seq1) // window_size
        return AlignmentWindowResult(
            name1=self.name1,
            name2=self.name2,
            score=self.score,
            matches=self.match_count,
            mismatches=self.mismatch_count,
            gaps=self.gap_count,
            window_text=combined,
            max_slider=max_slider
        )
