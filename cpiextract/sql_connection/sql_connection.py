'''Connect to mysql.'''

import time
import mysql.connector
from ..utils.typing import Connection

def connect_to_mysql(config:dict,
                     attempts: int = 3,
                     delay: int = 2) -> None | Connection:
    '''Connect to mysql with configs and attempts.
    
    Parameters
    ----------
    config: dist
        Config arguments passed into mysql's connect API.
    attempts: int
        Retrying times(3) for mysql connection.
    delay: int
        Wait time(2s) for reconnection.
    
    Returns
    -------
    None
        When connections fail.
    MySQLConnectionAbstract | PooledMySQLConnection
        A connection object when connections sussess.
    '''
    attempt = 0
    # Implement a reconnection routine
    while attempt < attempts:
        try:
            return mysql.connector.connect(**config)
        except (mysql.connector.Error, IOError) as err:
            if (attempts is attempt):
                # Attempts to reconnect failed; returning None
                print("Failed to connect, exiting without a connection: %s", err)
                return None
            print("Connection failed: %s. Retrying (%d/%d)...",
                err,
                attempt,
                attempts-1)
            # progressive reconnect delay
            time.sleep(delay ** attempt)
            attempt += 1
    return None