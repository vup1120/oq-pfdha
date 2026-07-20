# -*- coding: utf-8 -*-
# vim: tabstop=4 shiftwidth=4 softtabstop=4
#
# Copyright (C) 2024-2026 Yen-Shin Chen, OGS
#
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <https://www.gnu.org/licenses/>.

"""Cross-tool parity harness (oq-engine integration Phase 1, plan §8.2).

While the FDHA model library exists both here and inside the oq-engine
fork (``openquake/fdha`` in the engine tree), this module is the
regression net keeping the two implementations in sync:

* **model parity** — every model class ported to the engine must return
  bit-identical ``get_prob`` values over a magnitude grid;
* **scalerel parity** — the FDHA scaling relations added to
  ``openquake.hazardlib.scalerel`` (WC1994 AD/MD, Leonard2010 interplate,
  Thingbaijam2017 average slip, Leonard2014 widths) must match this
  repository's ``openquake.fdha.scalerel`` implementations;
* **distance parity** — the engine's ``rtor``/``x_l`` surface methods
  must match :class:`RuptureDistanceCalculator` exactly when both consume
  the same trace, and stay within the documented §4.4 envelope when the
  calculator uses the original NRML trace while the engine uses the
  resampled mesh top edge.

Import mechanics: both trees define ``openquake.fdha``, and this
repository shadows the engine's copy on ``sys.path``. The engine's
hazardlib is imported normally (it only exists in the engine tree); the
engine's ``openquake/fdha`` model modules are loaded by file path under
alias module names. When the engine checkout does not contain
``openquake/fdha`` (e.g. running against upstream master), the model
parity tests are skipped.
"""

import importlib.util
from pathlib import Path

import numpy as np
import pytest

import openquake.hazardlib as _hazardlib
from openquake.hazardlib.contexts import RuptureContext

# local (reference) implementations
from openquake.fdha.primary_surf_rup.youngs2003 import Youngs2003PrimarySR
from openquake.fdha.primary_surf_rup.wells_coppersmith1993 import (
    WC1993PrimarySR)
from openquake.fdha.primary_surf_rup.moss_ross2011 import MossRoss2011PrimarySR
from openquake.fdha.primary_surf_rup.takao2013 import Takao2013PrimarySR
from openquake.fdha.primary_surf_rup.moss2013 import Moss2013PrimarySR
from openquake.fdha.primary_surf_rup.yang2021 import Yang2021PrimarySR
from openquake.fdha.primary_surf_rup.pizza2023 import Pizza2023PrimarySR
from openquake.fdha.primary_surf_rup.fixed import FixedPrimarySR
from openquake.fdha.scalerel.wc1994 import WellsCoppersmith1994
from openquake.fdha.scalerel.leonard2010 import Leonard2010
from openquake.fdha.scalerel.thingbaijam2017 import ThingbaijamInterface
from openquake.fdha.primary_surf_rup.mammarella2024 import TAB1
from openquake.fdha.calc.utils.rupture_distance import (
    RuptureDistanceCalculator)

ENGINE_ROOT = Path(_hazardlib.__file__).resolve().parent.parent.parent
ENGINE_FDHA = ENGINE_ROOT / 'openquake' / 'fdha'

MAGS = np.round(np.arange(4.5, 8.51, 0.1), 2)

# Norcia MVFS trace, the shared fixture of the engine-side tests
NORCIA = [(13.1015, 43.0131), (13.1332, 42.9931), (13.1523, 42.9766),
          (13.1685, 42.9544), (13.1627, 42.9394), (13.1750, 42.9201),
          (13.1970, 42.9097), (13.2264, 42.8846), (13.2476, 42.8449),
          (13.2482, 42.8211), (13.2616, 42.8094), (13.2725, 42.7865),
          (13.2813, 42.7550)]


def _load_engine_module(stem):
    """Load ``openquake/fdha/primary_surf_rup/<stem>.py`` from the engine
    checkout under an alias name (the regular import path is shadowed by
    this repository's own ``openquake.fdha``)."""
    path = ENGINE_FDHA / 'primary_surf_rup' / (stem + '.py')
    if not path.exists():
        pytest.skip(f'engine checkout has no openquake/fdha '
                    f'({path} missing)')
    spec = importlib.util.spec_from_file_location(
        'engine_fdha_primary_surf_rup_' + stem, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def _ctx(**pairs):
    return RuptureContext(list(pairs.items()))


# --------------------------------------------------------------------------
# Model parity: engine openquake.fdha classes vs local reference classes
# --------------------------------------------------------------------------
class TestPrimarySurfRupParity:
    """get_prob parity over the magnitude grid, engine vs reference."""

    def test_youngs2003(self):
        mod = _load_engine_module('youngs2003')
        # engine _GB == local style='normal' (Great Basin subset);
        # the 'all' set is the WC1993 regression, tested below.
        # _ExC and _nBR have no local counterpart (regional additions).
        got = mod.Youngs2003PrimarySR_GB().get_prob(_ctx(mag=MAGS))
        expected = Youngs2003PrimarySR().get_prob(MAGS, style='normal')
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_wells_coppersmith1993(self):
        mod = _load_engine_module('wells_coppersmith1993')
        got = mod.WC1993PrimarySR().get_prob(_ctx(mag=MAGS))
        expected = WC1993PrimarySR().get_prob(MAGS)
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_moss_ross2011(self):
        mod = _load_engine_module('moss_ross2011')
        got = mod.MossRoss2011PrimarySR().get_prob(_ctx(mag=MAGS))
        expected = MossRoss2011PrimarySR().get_prob(MAGS)
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_takao2013(self):
        mod = _load_engine_module('takao2013')
        got = mod.Takao2013PrimarySR().get_prob(_ctx(mag=MAGS))
        expected = Takao2013PrimarySR().get_prob(MAGS)
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_yang2021(self):
        mod = _load_engine_module('yang2021')
        got = mod.Yang2021PrimarySR().get_prob(_ctx(mag=MAGS))
        expected = Yang2021PrimarySR().get_prob(MAGS)
        np.testing.assert_allclose(got, expected, rtol=1e-12)

    def test_pizza2023(self):
        mod = _load_engine_module('pizza2023')
        pairs = [(mod.Pizza2023PrimarySR, 'all'),
                 (mod.Pizza2023PrimarySR_Normal, 'normal'),
                 (mod.Pizza2023PrimarySR_Reverse, 'reverse'),
                 (mod.Pizza2023PrimarySR_SS, 'strike-slip')]
        ref = Pizza2023PrimarySR()
        for cls, style in pairs:
            got = cls().get_prob(_ctx(mag=MAGS))
            expected = ref.get_prob(MAGS, style=style)
            np.testing.assert_allclose(got, expected, rtol=1e-12,
                                       err_msg=style)

    def test_moss2013(self):
        mod = _load_engine_module('moss2013')
        ref = Moss2013PrimarySR()
        for cls, style in [(mod.Moss2013PrimarySR_Reverse, 'reverse'),
                           (mod.Moss2013PrimarySR_SS, 'strike-slip')]:
            for vs30 in (500.0, 800.0):
                got = cls().get_prob(_ctx(mag=MAGS, vs30=vs30))
                expected = ref.get_prob(MAGS, style=style, vs30=vs30)
                np.testing.assert_allclose(got, expected, rtol=1e-12,
                                           err_msg=f'{style} vs30={vs30}')

    def test_fixed(self):
        mod = _load_engine_module('fixed')
        got = mod.FixedPrimarySR(value=0.42).get_prob(_ctx(mag=MAGS))
        expected = FixedPrimarySR(value=0.42).get_prob(MAGS)
        np.testing.assert_allclose(got, np.full_like(MAGS, expected),
                                   rtol=1e-12)


# --------------------------------------------------------------------------
# Scalerel parity: engine hazardlib.scalerel additions vs local fdha.scalerel
# --------------------------------------------------------------------------
RAKE_BY_STYLE = {'strike-slip': 0.0, 'normal': -90.0, 'reverse': 90.0,
                 'all': None}


class TestScalerelParity:

    def test_wc1994_displacement(self):
        from openquake.hazardlib.scalerel.wc1994 import WC1994
        if not hasattr(WC1994, 'get_average_displacement'):
            pytest.skip('engine WC1994 has no FDHA displacement relations')
        eng, ref = WC1994(), WellsCoppersmith1994()
        for style, rake in RAKE_BY_STYLE.items():
            for mag in MAGS:
                ad, sad = eng.get_average_displacement(
                    mag, rake, return_sigma=True)
                ad_ref, sad_ref = ref.get_average_displacement(
                    mag, style, return_sigma=True)
                np.testing.assert_allclose(ad, ad_ref, rtol=1e-12)
                assert sad == sad_ref
                md, smd = eng.get_maximum_displacement(
                    mag, rake, return_sigma=True)
                md_ref, smd_ref = ref.get_maximum_displacement(
                    mag, style, return_sigma=True)
                np.testing.assert_allclose(md, md_ref, rtol=1e-12)
                assert smd == smd_ref

    def test_leonard2010_interplate(self):
        try:
            from openquake.hazardlib.scalerel.leonard2010 import (
                Leonard2010_Interplate)
        except ImportError:
            pytest.skip('engine has no Leonard2010_Interplate')
        eng, ref = Leonard2010_Interplate(), Leonard2010()
        np.testing.assert_allclose(
            eng.get_median_length(MAGS), ref.get_rupture_length(MAGS),
            rtol=1e-12)
        np.testing.assert_allclose(
            eng.get_average_displacement(MAGS),
            ref.get_average_displacement(MAGS), rtol=1e-12)
        np.testing.assert_allclose(
            eng.get_median_area(MAGS, 90.0),
            ref.get_median_area(MAGS, 90.0), rtol=1e-12)
        areas = np.array([10.0, 100.0, 1000.0, 5000.0, 20000.0])
        np.testing.assert_allclose(
            eng.get_median_mag(areas, 90.0),
            ref.get_median_mag(areas, 90.0), rtol=1e-12)

    def test_thingbaijam_displacement(self):
        from openquake.hazardlib.scalerel import thingbaijam2017 as eng_mod
        if not hasattr(eng_mod.ThingbaijamStrikeSlip,
                       'get_average_displacement'):
            pytest.skip('engine Thingbaijam has no displacement relations')
        ref = ThingbaijamInterface()
        pairs = [(eng_mod.ThingbaijamStrikeSlip, 'strike-slip'),
                 (eng_mod.ThingbaijamNormalFault, 'normal'),
                 (eng_mod.ThingbaijamReverseFault, 'reverse')]
        for cls, style in pairs:
            np.testing.assert_allclose(
                cls().get_average_displacement(MAGS),
                ref.get_average_displacement(MAGS, style), rtol=1e-12,
                err_msg=style)

    def test_leonard2014_width(self):
        from openquake.hazardlib.scalerel.leonard2014 import (
            Leonard2014_Interplate, Leonard2014_SCR)
        if not hasattr(Leonard2014_Interplate, 'get_median_width'):
            pytest.skip('engine Leonard2014 has no width relations')
        # reference: TAB1 of Mammarella et al. (2024) as implemented in
        # openquake.fdha.primary_surf_rup.mammarella2024 (MSR codes 0/1,
        # m = a + 2.5 log10 W); SoF codes 3/4 dip-slip, 5 strike-slip
        for msr_code, cls in ((0, Leonard2014_Interplate),
                              (1, Leonard2014_SCR)):
            for sof_code, rake in ((3, -90.0), (4, 90.0), (5, 0.0)):
                mask = (TAB1[:, 0] == msr_code) & (TAB1[:, 1] == sof_code)
                a = TAB1[mask][0, 2]
                expected = 10.0 ** ((MAGS - a) / 2.5)
                np.testing.assert_allclose(
                    cls().get_median_width(MAGS, rake), expected,
                    rtol=1e-12, err_msg=f'MSR={msr_code} SoF={sof_code}')
                sigma = TAB1[mask][0, 4]
                assert cls().get_std_dev_width(7.0, rake) == sigma


# --------------------------------------------------------------------------
# Distance parity: engine rtor / x_l vs RuptureDistanceCalculator
# --------------------------------------------------------------------------
class _Sites:
    """Duck-typed site collection over a lon/lat grid around Norcia."""

    def __init__(self):
        lons = np.linspace(13.00, 13.40, 21)
        lats = np.linspace(42.70, 43.10, 21)
        LO, LA = np.meshgrid(lons, lats)
        self.lons = LO.flatten()
        self.lats = LA.flatten()


def _norcia_surface():
    from openquake.hazardlib.geo import Line, Point
    from openquake.hazardlib.geo.surface.simple_fault import (
        SimpleFaultSurface)
    trace = Line([Point(lo, la) for lo, la in NORCIA])
    return SimpleFaultSurface.from_fault_data(
        trace, upper_seismogenic_depth=0.0, lower_seismogenic_depth=11.0,
        dip=47.0, mesh_spacing=1.0)


class TestDistanceParity:

    def setup_method(self, method):
        surf = _norcia_surface()
        if not hasattr(surf, 'get_tor_distance'):
            pytest.skip('engine surfaces have no get_tor_distance')
        from openquake.hazardlib.geo.mesh import Mesh
        self.surf = surf
        self.sites = _Sites()
        self.mesh = Mesh(self.sites.lons, self.sites.lats, depths=None)
        self.rtor = surf.get_tor_distance(self.mesh)
        self.x_l, self.l_km = surf.get_x_l_ratio(self.mesh)

    def test_same_trace_parity(self):
        """With both tools consuming the engine mesh top edge the
        projections must agree to numerical noise (the only difference
        is the anchoring of the orthographic projection)."""
        calc = RuptureDistanceCalculator(
            self.sites, self.surf, reference_line_method='segments')
        assert not calc.trace_is_original
        r_ref = calc.calculate_site_to_trace_distances()
        xl_ref, l_ref = calc.calculate_x_l_ratios()
        assert np.max(np.abs(self.rtor - r_ref)) < 5e-3      # km
        assert np.max(np.abs(self.x_l - xl_ref)) < 1e-4
        assert abs(self.l_km - l_ref) < 1e-3                 # km

    def test_original_trace_envelope(self):
        """Engine mesh-top-edge metrics vs the calculator on the original
        NRML trace: the deviation is the documented §4.4 corner-cutting
        effect of the resampled mesh (at 1 km spacing). The envelope is
        checked on distribution quantiles; single sites near trace kinks
        may project onto a different segment (worst observed |dx/L| ~0.2
        at one grid site 8 km on the footwall side)."""
        self.surf.original_trace = np.array(NORCIA)
        try:
            calc = RuptureDistanceCalculator(
                self.sites, self.surf, reference_line_method='segments')
            assert calc.trace_is_original
            r_ref = calc.calculate_site_to_trace_distances()
            xl_ref, l_ref = calc.calculate_x_l_ratios()
        finally:
            del self.surf.original_trace
        dr = np.abs(self.rtor - r_ref)
        dx = np.abs(self.x_l - xl_ref)
        assert np.median(dr) < 0.05          # km
        assert np.percentile(dr, 95) < 0.5   # km
        assert dr.max() < 1.0                # km
        assert np.median(dx) < 0.02
        assert np.percentile(dx, 95) < 0.05
        assert dx.max() < 0.25
        # trace length: mesh resampling shortens the wiggly trace a little
        assert abs(self.l_km - l_ref) < 0.5  # km

    def test_length_consistency(self):
        assert abs(self.surf.get_tor_length() - self.l_km) < 1e-9
