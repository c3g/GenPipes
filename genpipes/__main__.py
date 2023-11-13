"""
Stub function and module used as a setuptools entry point.
"""
import sys

import genpipes

# Entry point for setuptools pipelines.
def main():
    return genpipes.main( sys.argv[1:] )

# Run when called as `python -m augur`, here for good measure.
if __name__ == "__main__":
    sys.exit( main() )