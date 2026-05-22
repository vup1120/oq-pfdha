import logging
from .calculators import (
    BaseFaultRuptureCalculator,
    FaultRuptureProbabilityCalculator,
)
from .hazard_map import compute_hazard_map
from .io_utils import (
    plot_hazard_map,
    plot_fault_displacement_hazard,
    save_results_to_json,
    save_map_to_json,
)

logger = logging.getLogger(__name__)


def run_calculation(
    config_path: str,
    source_model_path: str,
    plot: bool = False,
    output_file: str = None,
    plot_file: str = None,
    return_period: float = None,
    calculation_type: str = None,
):
    logger.info(f"Starting calculation with config={config_path}, source_model={source_model_path}, "
                f"plot={plot}, output={output_file}, plot_file={plot_file}, "
                f"return_period={return_period}, calculation_type={calculation_type}")

    config = BaseFaultRuptureCalculator._load_configuration(config_path)
    site_config = config.get("site_location", {})

    # Define the calculator type: hazard map or curve
    if calculation_type is None:
        if "geometry" in config:
            calculation_type = "hazard_map"
        else:
            calculation_type = "hazard_curve"

    if calculation_type == "hazard_map":
        hazard_map, lons, lats, fault_lons, fault_lats, fault_displacements = compute_hazard_map(
            config_path,
            source_model_path,
        )

        if plot:
            plot_hazard_map(
                hazard_map,
                lons,
                lats,
                fault_lons,
                fault_lats,
                fault_displacements,
                plot_file=plot_file,
            )
        if output_file:
            save_map_to_json(
                hazard_map,
                lons,
                lats,
                fault_lons,
                fault_lats,
                fault_displacements,
                output_file=output_file,
            )
        return hazard_map

    elif calculation_type == "hazard_curve":
        has_single_site = "latitude" in site_config and "longitude" in site_config
        has_multi_site = "sites_list" in site_config
        if has_single_site or has_multi_site:
            calculator = FaultRuptureProbabilityCalculator(
                config_path, source_model_path
            )
        else:
            raise ValueError(
                "Hazard curve requires 'latitude' and 'longitude', "
                "or a multi-site 'sites'/'sites_csv' entry in [geometry]"
            )

        results = calculator.run()

        if plot:
            plot_fault_displacement_hazard(results, plot_file=plot_file)
        if output_file:
            save_results_to_json(results, output_file)
        return results

    else:
        raise ValueError(f"Unsupported calculation_type: {calculation_type}")
