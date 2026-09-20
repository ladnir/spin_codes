"""Construct binomial composition tables only when a box requests them."""
import occupation_refresh_v1 as transfer


class BinomialTable:
    def __init__(self, probability):
        self.probability = probability
        self.cache = {}

    def __getitem__(self, count):
        if count not in self.cache:
            self.cache[count] = transfer.boxes.binomial_logs(count, self.probability)
        return self.cache[count]


class CompositionBoxes(transfer.CompositionBoxes):
    tables = {}

    def __init__(self, counts, block, bands, probabilities, maximum):
        # The density measure and envelope do not depend on the ceiling.
        # Avoid eagerly building O(maximum^2) coefficients for every bank.
        super().__init__(counts, block, bands, probabilities, 1)
        self.maximum = maximum
        self.binomials = []
        for p in probabilities:
            if p not in self.tables:
                self.tables[p] = BinomialTable(p)
            self.binomials.append(self.tables[p])
