# Norcia Sensitivity Test - Youngs2003 AD 85

This test verifies the implementation against reference values from the paper (Youngs2003 AD 85).

## Running the Test

### Option 1: Using the test script

```bash
cd /home/ychen/GIT/pfdha/openquake/fdha/test/benchmark/norcia_sensitivity_youngs2003
python run_test.py
```

This will:
1. Run the calculation using the current implementation
2. Save results to `results.json`

### Option 2: Compare with reference values

After running the calculation, compare with reference values:

```bash
python compare_results.py
```

This will compare the current results with reference values from the paper and report statistics.

### Option 3: Generate comparison plot

```bash
python plot_results.py
```

This will generate `hazard_curve_comparison.png` showing:
- Hazard curve comparison (Current vs Reference)
- Ratio plot (Current / Reference)
- Statistics

## Files

- `config_youngs2003_AD_85.toml` - Configuration file
- `source_model_norcia.xml` - Source model
- `results.json` - Current calculation results
- `run_test.py` - Run calculation script
- `compare_results.py` - Compare with reference values
- `plot_results.py` - Generate comparison plot
- `hazard_curve_comparison.png` - Generated comparison plot

## Reference Values

The reference values are from the paper (Youngs2003 AD 85) and are hardcoded in the comparison scripts.

## Expected Results

The current implementation should match the reference values within reasonable tolerance (rtol=1e-2, atol=1e-6).
