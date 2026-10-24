# Script to extract data and setup datasets
import re
from pathlib import Path

# Create data directory
Path('data').mkdir(exist_ok=True)
Path('output').mkdir(exist_ok=True)

print("Data directory created.")
