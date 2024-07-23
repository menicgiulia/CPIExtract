'''Util functions and metaclasses.'''

# import pandas as pd
# import numpy as np

from typing import Generator


def generate_subsets(set_size: int, n: int) -> Generator[tuple[int, int], None, None]:
    """
    Generate starting and ending subset indices, diving set_size into n equal subsets
    
    Parameters
    ----------
    set_size : int
        size of the set to divide
    n : int
        number of subset to generate

    Yields
    -------
    int
        index of the first element of current subset
    int
        index of the last element of current subset
    """
    start = 0
    while start < set_size:
        end = min(start + n, set_size)
        yield start, end
        start = end

class Singleton(type):
    '''The metaclass for singleton design pattern.'''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]
