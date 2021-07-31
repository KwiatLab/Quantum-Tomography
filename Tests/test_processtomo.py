import unittest
from QuantumTomography.ProcessTomo import ProcessTomography
import numpy as np

coinc = np.random.binomial(1000, 0.5, (6,6))



print(ProcessTomography(coinc))

