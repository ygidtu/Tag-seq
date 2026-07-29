#!/usr/bin/env python3
"""Tag-seq pipeline entry point.

Usage:
    python main.py -c config.txt all
    python main.py -c config.txt --resume align
    python main.py -c config.txt report
"""

import sys
import os

# Ensure src/ is on the path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from tagseq.cli import main

if __name__ == "__main__":
    main()
