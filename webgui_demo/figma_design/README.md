# Figma Make design — extracted source (for review)

This folder holds the **design you created in Figma Make**, extracted
from `Web_UI_for_Hazard_Analysis.make`. A `.make` file is a ZIP bundle of
a generated **React + Tailwind** app plus its chat history and image
assets — it is *not* a `.fig` vector file, and there is nothing to import
into Canva. The design *is* this code.

These files are kept here **for reference only**; they are not wired into
the build. The live app is still the Streamlit GUI in `webgui_demo/app.py`.

## Renders

- `render_dashboard.png` — the polished dashboard concept (navy nav rail,
  top bar, map + hazard-curve + displacement-summary cards).
- `render_configure_page.png` — a restyle of the current Configure page.

## Design tokens (extracted)

| Token | Value |
|---|---|
| Nav rail / primary navy | `#0f2846` (hover `#1a3d66`) |
| Accent (links, active, focus) | blue-600 `#2563eb` / blue-500 `#3b82f6` |
| App background | `#f3f4f6` |
| Card surface / borders | white / `slate-200` `#e2e8f0` |
| Section header tint | `slate-50` |
| Text | `slate-800` body, `slate-500/600` labels |
| Chart series | Petersen `#0f2846`, Youngs `#3b82f6`, Mean `#f97316` (dashed) |
| Sans font | **Inter** (300–700) |
| Serif / Mono | Source Serif 4 / JetBrains Mono |
| Aesthetic | OGS institutional — deep navy, white, cool grays |

(Source: `src/styles/fonts.css`, and the Tailwind classes in each
component. The Figma Make theme prompt explicitly targeted the OGS
institutional palette.)

## Component map

```
src/app/App.tsx                       router root
src/app/routes.tsx                     single dashboard route (rest are stubs)
src/app/components/RootLayout.tsx      sidebar + header + <Outlet/> shell
src/app/components/Sidebar.tsx         navy nav rail (Dashboard/Map/Data/...)
src/app/components/Header.tsx          top bar: title, search, bell, avatar
src/app/components/Dashboard.tsx       grid: params | (map+chart) / table
src/app/components/ParameterPanel.tsx  accordion: geometry/site/models + Run
src/app/components/MapView.tsx         map card with site marker + fault line
src/app/components/HazardCurveChart.tsx Recharts log-y hazard curves
src/app/components/ResultsTable.tsx    displacement summary table
```

## ⚠️ Important: this is a visual prototype, not connected to the engine

Every number in the design is **hardcoded mock data**, e.g.:

- `ResultsTable.tsx` — `tableData` (Return Period / Exceedance / Petersen /
  Youngs / Mean) is a static array.
- `HazardCurveChart.tsx` — `mockData` is a static array.
- `ParameterPanel.tsx` — inputs have default values but no wiring; the
  "Run Analysis" button does nothing.

To make this design *actually run PFDHA*, it needs a backend (e.g. FastAPI)
that calls `FdhaLogicTree.from_ini(...).run(...)` and feeds real results
into these components in place of the mock arrays. That is "Path B"
discussed with the project owner.

## Two ways to use this design

- **Path A — restyle Streamlit** to match these tokens (navy sidebar,
  Inter, white cards). Faster; keeps the working engine; ~80% visual
  fidelity (Streamlit cannot reproduce the fixed rail + top-bar + grid
  exactly).
- **Path B — adopt this React frontend for real**, adding a FastAPI layer
  that replaces the mock arrays with live engine output. Pixel-perfect to
  the Figma design; a larger build and a new deployment story.
