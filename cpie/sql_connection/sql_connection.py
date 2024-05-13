import time
import mysql.connector

def connect_to_mysql(config, attempts=3, delay=2):
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