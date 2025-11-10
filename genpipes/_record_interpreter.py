import sys
import os

target_file = os.path.join(os.path.dirname(__file__), "_interpreter.py")
with open(target_file, "w") as f:
    f.write(f'INTERPRETER = "{sys.executable}"\n')
