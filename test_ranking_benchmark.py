"""Release gate wrapper for the frozen ranking-truth benchmark."""
import os, sys
sys.path.insert(0, os.path.dirname(__file__))
from run_ranking_benchmark import main
raise SystemExit(main())
