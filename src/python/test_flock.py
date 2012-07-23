from filelock import FileLock
import time

with FileLock("blade_surf_copy.dat"):
    print("Lock acquired")
    time.sleep(5)

print("Lock released")
    
