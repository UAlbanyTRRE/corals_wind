# =============================================================================
# CMIP6 projected mean-wind change (delta U) at every reef cell
# Part of: Wind effects on coral bleaching severity (Lapenis)
# Language: Python 3 (designed for Google Colab; runs locally with minor edits)
#
# WHAT IT DOES
#   Reads the equal-area reef coordinates, selects the FULL available CMIP6
#   ensemble from the Pangeo Google-Cloud catalogue (one realisation per model;
#   ~43 models for ssp585, ~41 for ssp245), computes delta U = future - historical
#   per model, samples it to every reef cell, writes one delta-U column per model
#   beside the coordinates, computes one-model-one-vote and one-centre-one-vote
#   ensemble means, the sign-agreement statistics, and the GLO/SAF bleaching-odds
#   fields, and produces the wind-change map, histogram, and odds maps.
#
# INPUT  (CONFIG below):  data/input/reef_coordinates_5km.xlsx  (cols LATITUDE, LONGITUDE)
# OUTPUT (written to cwd / data/output): delU_per_reef.{parquet,csv,xlsx} + figures
#
# HOW TO RUN
#   Google Colab : New notebook -> paste -> Run all (you will be prompted to
#                  upload the coordinates file).
#   Locally      : set IN_COORDS to your file and run `python 05_cmip6_wind_delU_per_reef.py`.
# =============================================================================
# Deposited code does not auto-install. Verify required packages are present and
# point to requirements.txt if not. (Install with: pip install -r requirements.txt)
# In Google Colab you may still need to `pip install` cartopy/gcsfs in a cell first.
import importlib.util, sys
_required = {  # module name -> pip name (shown to the user if missing)
    "xarray": "xarray", "gcsfs": "gcsfs", "zarr": "zarr", "cftime": "cftime",
    "pandas": "pandas", "numpy": "numpy", "matplotlib": "matplotlib",
    "openpyxl": "openpyxl", "pyarrow": "pyarrow",
}
_missing = [pip for mod, pip in _required.items() if importlib.util.find_spec(mod) is None]
if _missing:
    sys.exit("Missing required package(s): " + ", ".join(_missing) +
             "\nInstall them with:  pip install -r requirements.txt")
# Note: cartopy is OPTIONAL (maps fall back to no coastlines if it is absent).

import os, re, warnings
from collections import defaultdict
import pandas as pd, numpy as np, xarray as xr, gcsfs
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")

try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    HAVE_CARTOPY = True
except Exception:
    HAVE_CARTOPY = False
    print("   (cartopy unavailable -> maps drawn without coastlines)")

# ---- 0. I/O CONFIG ---------------------------------------------------
# Equal-area reef coordinates (columns LATITUDE, LONGITUDE), e.g. 5x5 km cells.
IN_COORDS = "data/input/reef_coordinates_5km.xlsx"
OUT_DIR   = "data/output"
os.makedirs(OUT_DIR, exist_ok=True)

# ---- 1. settings -----------------------------------------------------
VARIABLE   = 'sfcWind'
TABLE      = 'Amon'
HIST_EXP   = 'historical'
SSPS       = ['ssp585', 'ssp245']
HIST_YEARS = ('1995', '2014')
FUT_YEARS  = ('2080', '2099')

REQUIRE_R1F1         = False           # True -> reproduces the ~8-model run
ONE_MEMBER_PER_MODEL = True
PREFER_GRID          = ['gn', 'gr', 'gr1']

# Sign-agreement thresholds reported below
AGREE_LO, AGREE_HI = 0.66, 0.80

# ---- bleaching-odds coefficients ------------------------------------
# These are NOT re-estimated here; they are the published wind odds ratios fed
# into the projection. Provenance:
#   OR_GLO = 0.70  per s.d.  -> GLO within-region 6-month mean-wind OR, from
#            06_glo_shapley_mundlak.R / 11_glo_selection_robustness.R
#            (6-mo OR ~0.70, p~2e-7); Results / Supplementary robustness table.
#   OR_SAF = 0.45  per s.d.  -> SAF 12-month mean-wind OR (n=81), from
#            02_saf_mediation_residual_wind.R; main SAF results figure.
#   SD     = 1.20  m/s       -> s.d. of the wind metric used to scale the OR.
# Update these if the corresponding regression outputs change.
# GLO: within-region 6-month mean wind (conceptual match to a sustained
#      multidecadal change in CLIMATOLOGICAL mean wind) -> primary.
OR_GLO, SD_GLO = 0.70, 1.20            # odds ratio per s.d.; s.d. in m/s
# SAF: 12-month mean wind preceding events (steeper slope, n=81) -> sensitivity.
OR_SAF, SD_SAF = 0.45, 1.20            # set SD_SAF to the SAF metric's s.d. (m/s)

# ---- 2. reef polygons (equal-area) ----------------------------------
try:
    from google.colab import files
    print(">> Upload the reef coordinates file (cols LATITUDE, LONGITUDE) or a lon,lat CSV:")
    up = files.upload(); fname = list(up.keys())[0]
    IN_USE = fname
except Exception:
    IN_USE = IN_COORDS
    print(f">> Reading reef coordinates from {IN_USE}")

reef = pd.read_csv(IN_USE) if IN_USE.lower().endswith('.csv') else pd.read_excel(IN_USE)
reef.columns = [c.strip().lower() for c in reef.columns]
reef = reef.rename(columns={'latitude': 'lat', 'longitude': 'lon'})[['lon', 'lat']].dropna()
reef = reef.reset_index(drop=True)
reef_lat = reef['lat'].to_numpy()
reef_lon360 = reef['lon'].to_numpy() % 360.0
print(f"   {len(reef):,} reef cells "
      f"(lat {reef.lat.min():.1f}..{reef.lat.max():.1f}, "
      f"lon {reef.lon.min():.1f}..{reef.lon.max():.1f})")

# ---- 3. catalog + cloud filesystem ----------------------------------
CAT = 'https://storage.googleapis.com/cmip6/cmip6-zarr-consolidated-stores.csv'
fs  = gcsfs.GCSFileSystem(token='anon')
cat = pd.read_csv(CAT)
cat = cat[(cat.variable_id == VARIABLE) & (cat.table_id == TABLE)]


def member_rank(m):
    mm = re.match(r'r(\d+)i(\d+)p(\d+)f(\d+)', str(m))
    if not mm:
        return (9, 99, 99, 99)
    r, i, p, f = map(int, mm.groups())
    return (f, r, i, p)


def grid_rank(g):
    return PREFER_GRID.index(g) if g in PREFER_GRID else len(PREFER_GRID)


def pick_member_grid(hist_rows, ssp_rows):
    common_mem = set(hist_rows.member_id) & set(ssp_rows.member_id)
    if REQUIRE_R1F1:
        common_mem &= {'r1i1p1f1'}
    if not common_mem:
        return None, "no shared member" + (" (r1i1p1f1 required)" if REQUIRE_R1F1 else "")
    mem = sorted(common_mem, key=member_rank)[0]
    h = hist_rows[hist_rows.member_id == mem]
    s = ssp_rows[ssp_rows.member_id == mem]
    common_grid = set(h.grid_label) & set(s.grid_label)
    if not common_grid:
        return None, f"member {mem} has no shared grid_label"
    grid = sorted(common_grid, key=grid_rank)[0]
    hz = h[h.grid_label == grid].sort_values('version').zstore.iloc[-1]
    sz = s[s.grid_label == grid].sort_values('version').zstore.iloc[-1]
    inst = h.institution_id.iloc[0]
    return (mem, grid, inst, hz, sz), None


def open_clim(zstore, years):
    ds = xr.open_zarr(fs.get_mapper(zstore), consolidated=True, use_cftime=True)
    da = ds[VARIABLE]
    for a, b in [('latitude', 'lat'), ('longitude', 'lon')]:
        if a in da.coords:
            da = da.rename({a: b})
    da = da.sel(time=slice(years[0], years[1])).mean('time').squeeze(drop=True)
    return da.load()


def sample_to_reefs(dU):
    return dU.sel(lat=xr.DataArray(reef_lat, dims='reef'),
                  lon=xr.DataArray(reef_lon360, dims='reef'),
                  method='nearest').values.astype('float32')


def add_coast(ax):
    if HAVE_CARTOPY:
        ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='0.35')
        ax.add_feature(cfeature.LAND, facecolor='0.93', zorder=0)


# ---- 4. process every SSP, every model ------------------------------
result = reef.copy()
provenance = []

for ssp in SSPS:
    hist = cat[cat.experiment_id == HIST_EXP]
    futr = cat[cat.experiment_id == ssp]
    models = sorted(set(hist.source_id) & set(futr.source_id))
    print(f"\n=== {ssp}: {len(models)} models have both {HIST_EXP} and {ssp} ===")
    kept = {}
    for src in models:
        sel, reason = pick_member_grid(hist[hist.source_id == src], futr[futr.source_id == src])
        if sel is None:
            print(f"   [drop] {src:22s} {reason}");  continue
        mem, grid, inst, hz, sz = sel
        try:
            dU = open_clim(sz, FUT_YEARS) - open_clim(hz, HIST_YEARS)
            col = f"dU_{ssp}_{src}"
            result[col] = sample_to_reefs(dU)
            kept[col] = inst
            provenance.append({'experiment': ssp, 'source_id': src, 'institution_id': inst,
                               'member_id': mem, 'grid_label': grid, 'hist_zstore': hz, 'ssp_zstore': sz})
            print(f"   [keep] {src:22s} {inst:14s} {mem:10s} {grid}")
        except Exception as e:
            print(f"   [fail] {src:22s} {type(e).__name__}: {str(e)[:60]}")

    cols = list(kept.keys())
    if not cols:
        continue
    block = result[cols].to_numpy()

    # one-model-one-vote ensemble
    ens = np.nanmean(block, axis=1)
    result[f"ens_mean_{ssp}"]   = ens.astype('float32')
    result[f"ens_median_{ssp}"] = np.nanmedian(block, axis=1).astype('float32')
    result[f"ens_std_{ssp}"]    = np.nanstd(block, axis=1).astype('float32')
    result[f"n_models_{ssp}"]   = np.sum(~np.isnan(block), axis=1).astype('int16')
    frac_neg = np.mean(block < 0, axis=1)
    result[f"frac_neg_{ssp}"]   = frac_neg.astype('float32')

    # one-center-one-vote ensemble (average within institution, then across)
    groups = defaultdict(list)
    for c, inst in kept.items():
        groups[inst].append(c)
    inst_means = np.column_stack([result[c2].to_numpy() if len(c2) == 1
                                  else result[c2].mean(axis=1).to_numpy()
                                  for c2 in groups.values()])
    ens_center = np.nanmean(inst_means, axis=1)
    result[f"ens_mean_center_{ssp}"] = ens_center.astype('float32')

    # bleaching-odds fields (per cell, from one-model-one-vote mean)
    beta_glo = np.log(OR_GLO) / SD_GLO
    beta_saf = np.log(OR_SAF) / SD_SAF
    result[f"OR_GLO_{ssp}"] = np.exp(beta_glo * ens).astype('float32')
    result[f"OR_SAF_{ssp}"] = np.exp(beta_saf * ens).astype('float32')

    # ---- headline statistics (area-weighted, since cells are equal-area) ----
    p5, p50, p95 = np.percentile(ens, [5, 50, 95])
    agree_decline = np.mean(frac_neg >= AGREE_LO)
    agree_sign    = np.mean(np.maximum(frac_neg, 1 - frac_neg) >= AGREE_HI)
    print(f"\n   --- {ssp} headline (n={len(cols)} models, {len(groups)} centers) ---")
    print(f"   area-mean  delU  (1 model/vote)  = {ens.mean():+.4f} m/s")
    print(f"   area-mean  delU  (1 center/vote) = {ens_center.mean():+.4f} m/s")
    print(f"   area-weighted median delU        = {p50:+.4f} m/s")
    print(f"   5th / 95th percentile delU       = {p5:+.3f} / {p95:+.3f} m/s")
    print(f"   reef AREA with declining winds   = {np.mean(ens < 0)*100:.1f} %")
    print(f"   AREA where >={int(AGREE_LO*100)}% models decline = {agree_decline*100:.1f} %")
    print(f"   AREA where >={int(AGREE_HI*100)}% models agree on sign = {agree_sign*100:.1f} %")
    print(f"   GLO odds at area-mean = x{result[f'OR_GLO_{ssp}'].mean():.3f} "
          f"|  SAF odds at area-mean = x{result[f'OR_SAF_{ssp}'].mean():.3f}")

# ---- 5. save ---------------------------------------------------------
prov = pd.DataFrame(provenance)
result.to_parquet(os.path.join(OUT_DIR, "delU_per_reef.parquet"), index=False)
result.to_csv(os.path.join(OUT_DIR, "delU_per_reef.csv"), index=False)
with pd.ExcelWriter(os.path.join(OUT_DIR, "delU_per_reef.xlsx"), engine="openpyxl") as xw:
    result.to_excel(xw, sheet_name="delU_per_reef", index=False)
    prov.to_excel(xw,   sheet_name="models_used",  index=False)
print(f"\nSaved {len(result):,} reef cells x {result.shape[1]} columns to {OUT_DIR}/")

# ---- 6. FIGURES ------------------------------------------------------
ssp_plot = 'ssp585' if 'ens_mean_ssp585' in result else SSPS[0]
proj = ccrs.PlateCarree() if HAVE_CARTOPY else None

# (a) wind-change map with coastlines
c = result[f"ens_mean_{ssp_plot}"]
vmax = np.nanpercentile(np.abs(c), 99)
fig = plt.figure(figsize=(11, 4.4))
ax = fig.add_subplot(111, projection=proj) if HAVE_CARTOPY else fig.add_subplot(111)
add_coast(ax)
kw = dict(transform=proj) if HAVE_CARTOPY else {}
sc = ax.scatter(reef['lon'], reef['lat'], c=c, s=2, cmap='RdBu', vmin=-vmax, vmax=vmax, **kw)
ax.set_xlim(-180, 180); ax.set_ylim(-40, 40)
ax.set_xlabel('Longitude (E)'); ax.set_ylabel('Latitude (N)')
ax.set_title(f'Projected mean-wind change over reefs, {ssp_plot} '
             f'({FUT_YEARS[0]}-{FUT_YEARS[1]} vs {HIST_YEARS[0]}-{HIST_YEARS[1]})')
plt.colorbar(sc, ax=ax, label=r'$\Delta U$ (m s$^{-1}$)', shrink=0.85)
plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "map_delU_coast.png"), dpi=160); plt.show()

# (b) histogram
fig, ax = plt.subplots(figsize=(6, 3))
ax.hist(c, bins=60, color='steelblue'); ax.axvline(0, ls='--', c='grey')
ax.set_xlabel(fr'$\Delta U$ per reef (m s$^{{-1}}$), {ssp_plot}'); ax.set_ylabel('reef polygons')
plt.tight_layout(); plt.savefig(os.path.join(OUT_DIR, "hist_delU.png"), dpi=160); plt.show()

# (c) bleaching-odds maps: GLO (primary) and SAF (sensitivity)
for tag, col in [('GLO 6-month coefficient', f'OR_GLO_{ssp_plot}'),
                 ('SAF 12-month coefficient', f'OR_SAF_{ssp_plot}')]:
    pct = (result[col] - 1.0) * 100.0
    vlim = np.nanpercentile(np.abs(pct), 99)
    fig = plt.figure(figsize=(11, 4.4))
    ax = fig.add_subplot(111, projection=proj) if HAVE_CARTOPY else fig.add_subplot(111)
    add_coast(ax)
    sc = ax.scatter(reef['lon'], reef['lat'], c=pct, s=2, cmap='RdBu_r',
                    vmin=-vlim, vmax=vlim, **kw)
    ax.set_xlim(-180, 180); ax.set_ylim(-40, 40)
    ax.set_xlabel('Longitude (E)'); ax.set_ylabel('Latitude (N)')
    ax.set_title(f'Implied change in severe-bleaching odds, {ssp_plot}\n({tag})')
    plt.colorbar(sc, ax=ax, label='change in odds (%)', shrink=0.85)
    plt.tight_layout()
    plt.savefig(os.path.join(OUT_DIR, f"map_odds_{'GLO' if 'GLO' in tag else 'SAF'}.png"), dpi=160); plt.show()
