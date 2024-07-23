'''Type alias for other modules.'''

from typing import Sequence, TypeVar, Generic, Generator, Callable
from mysql.connector.abstracts import MySQLConnectionAbstract
from mysql.connector.pooling import PooledMySQLConnection

Connection = MySQLConnectionAbstract | PooledMySQLConnection
'''SQL connection type.'''

T = TypeVar('T')
'''Immutable type.'''

T_co = TypeVar('T_co', covariant=True)
'''Any type covariant containers.'''

T_contra = TypeVar('T_contra', contravariant=True)
'''Any type contravariant containers.'''

Strs = str | list[str]
'''Single str or a list of strs.'''

CIDs = Strs | int | list[int]
'''Compound IDs,may be int(s) also.'''
