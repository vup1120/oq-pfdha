#!/usr/bin/env python
"""Run norcia test and compare with reference values."""
import json
import sys
import os
import importlib
from pathlib import Path

import numpy as np

# Add parent directories to path
sys.path.insert(0, str(Path(__file__).parent.parent.parent.parent.parent))

# Force reload modules to avoid Python cache issues
if 'openquake.fdha.calc.hazard' in sys.modules:
    importlib.reload(sys.modules['openquake.fdha.calc.hazard'])
if 'openquake.fdha.calc.calculators' in sys.modules:
    importlib.reload(sys.modules['openquake.fdha.calc.calculators'])

from openquake.fdha.calc.calculators import FaultRuptureProbabilityCalculator

def main():
    test_dir = Path(__file__).parent
    config_path = test_dir / "config_youngs2003_AD_85.toml"
    source_model_path = test_dir / "source_model_norcia.xml"
    results_path = test_dir / "results.json"
    
    print(f"Running calculation with config: {config_path}")
    print(f"Source model: {source_model_path}")
    
    # Run calculation
    calculator = FaultRuptureProbabilityCalculator(
        str(config_path),
        [str(source_model_path)],
        rupture_mesh_spacing=0.1,
        width_of_mfd_bin=0.1
    )
    
    results = calculator.run()
    
    # Save results
    with open(results_path, 'w') as f:
        json.dump({
            'imls': np.asarray(results['imls']).tolist(),
            'poes': np.asarray(results['poes']).tolist()
        }, f, indent=2)
    
    print(f"\nResults saved to: {results_path}")
    print("\n" + "="*80)
    print("To compare with reference values, run:")
    print(f"  python {test_dir / 'compare_results.py'}")
    print(f"  python {test_dir / 'plot_results.py'}")
    print("="*80)
    
    return 0

if __name__ == "__main__":
    sys.exit(main())
