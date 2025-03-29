from numpy import random

from aliases.common import *

def random_biases(size: int) -> Vector:
    return random.rand(size)

def random_weights(rows: int, cols: int) -> Vector:
    return random.rand(rows, cols)
