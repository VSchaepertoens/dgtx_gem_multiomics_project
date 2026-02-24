import os
import hashlib

# Get the directory where this script is located
base_dir = os.path.dirname(os.path.abspath(__file__))

def sha256sum(filename):
    h = hashlib.sha256()
    with open(filename, 'rb') as f:
        for chunk in iter(lambda: f.read(8192), b''):
            h.update(chunk)
    return h.hexdigest()

files = [
    os.path.join(base_dir, "..", "models", "original", "iCHOv1.xml"),
    os.path.join(base_dir, "..", "models", "original", "iCHO2441.xml"),
    os.path.join(base_dir, "..", "models", "original", "iCHO3K_cho_prod_generic_unblocked.xml")
]

for f in files:
    print(f"{f}: {sha256sum(f)}")