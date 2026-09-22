'''Util functions and metaclasses.'''

# import pandas as pd
# import numpy as np

import time
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


def call_with_retry(func, *args, max_retries: int = 3, backoff_base: float = 1.0,
                     verbose: bool = False, **kwargs):
    """
    Calls func with retry-with-backoff. Originally written for PubChem's PUG-REST
    service specifically (confirmed against a real ServerBusyError), but the retry
    mechanism itself has nothing PubChem-specific about it - also used for ChEMBL API
    calls (see Database.py's _pubchem_search_chembl). Without this, a single transient
    failure inside a per-row/per-batch loop is indistinguishable from "genuinely not
    found", silently leaving that row's data unresolved with no indication anything
    went wrong, rather than surfacing as a retryable, temporary issue. Retries
    max_retries times total, waiting backoff_base * 2^attempt seconds between attempts
    (1s, 2s, 4s by default) before letting the final exception propagate to the
    caller's own handling.

    verbose=True prints each retry attempt as it happens (used by
    utils/identifiers.py's compound_identifiers(); Database.py's own callers don't
    currently pass this, so their retries stay silent by default, matching prior
    behavior there).

    Shared here (rather than duplicated in Database.py and utils/identifiers.py, which
    both use it) specifically to avoid the two copies drifting out of sync if either
    is tuned later.
    """
    last_exception = None
    for attempt in range(max_retries):
        try:
            return func(*args, **kwargs)
        except Exception as e:
            last_exception = e
            if attempt < max_retries - 1:
                wait = backoff_base * (2 ** attempt)
                if verbose:
                    print(f"  Call failed (attempt {attempt + 1}/{max_retries}): "
                          f"{e} - retrying in {wait:.1f}s...")
                time.sleep(wait)
    raise last_exception

class Singleton(type):
    '''The metaclass for singleton design pattern.'''
    _instances = {}
    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super(Singleton, cls).__call__(*args, **kwargs)
        return cls._instances[cls]