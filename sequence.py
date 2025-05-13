from abc import ABC, abstractmethod

class Sequence(ABC):
    def __init__(self, df):
        self.__info = df

    @property
    @abstractmethod
    def sequence(self):
        pass

    @property
    @abstractmethod
    def length(self):
        pass

class MitochondrialDNA(Sequence):
    def __init__(self, df):
        super().__init__(df)
        self.__sequence = df['seq']
        self.__length = df['length']
        self.__id = df['id']
        self.__name = df['name']
        self.__description = df['description']

    @property
    def sequence(self):
        return self.__sequence

    @property
    def length(self):
        return self.__length

    @property
    def gc_content(self):
        return ((self.sequence.count('G') + self.sequence.count('C')) / self.length) * 100

    def get_subsequence(self, start: int, end: int):
        if start < 0 or end > len(self.sequence):
            raise ValueError(f"Subsequence indices out of range: start={start}, end={end}, length={self.length}")
        else:
            return self.sequence[start:end]

    def find_irregular_bases(self):
        valid_bases = {'A', 'T', 'G', 'C'}
        irregular = []
        for base in self.sequence.upper():
            if base not in valid_bases:
                irregular.append(base)
        return irregular


if __name__ == "__main__":
    from tools import Parser
    parser = Parser()
    records = parser.run("synthetic_mtDNA_dataset.fasta")
    NC10 = MitochondrialDNA(records.loc[10])
    print(NC10.gc_content)
    print(NC10.find_irregular_bases())