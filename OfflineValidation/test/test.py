import os
import time
while True:
    env = os.getenv("testvar")
    if env is not None:
        print env
    time.sleep(5)
