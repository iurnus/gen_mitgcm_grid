#!/usr/bin/env python3
"""
Make an MITgcm bathymetry file (bathy.bin) and grid settings for a regional box.

HOW TO USE
  1. Edit the five settings below.
  2. Run:  python3 gen_bathy.py
  3. Copy bathy.bin and the lines in data_grid.txt into your MITgcm input/data.

Works with ETOPO 2022 / GEBCO (NetCDF4, needs: pip install netCDF4) and with
older NetCDF3 files such as ETOPO2v2 (no extra package needed).
"""
import numpy as np
from scipy import ndimage

# ============================== SETTINGS ===================================
ELEVATION_FILE = "ETOPO_2022_v1_60s_N90W180_surface.nc"

LON_WEST, LON_EAST = 30, 50      # Red Sea example
LAT_SOUTH, LAT_NORTH = 10, 30
# LON_WEST, LON_EAST = 170, -170 # date-line example: just give west, then east

RESOLUTION = 0.08                # degrees

OPEN_SIDES = "E"                 # sides that stay open to the ocean (for OBCS):
                                 # any of "N", "S", "E", "W", or "" for a closed box
# ===========================================================================

MIN_DEPTH = 10.0                 # m; shallower ocean becomes land or is deepened
NUM_LEVELS = 50                  # vertical levels, 5 m at the top


# ---- 1. model grid ---------------------------------------------------------
west, east = LON_WEST, LON_EAST
if east <= west:                 # box crosses the date line, e.g. 170 -> -170
    east += 360
nx = round((east - west) / RESOLUTION)
ny = round((LAT_NORTH - LAT_SOUTH) / RESOLUTION)
print(f"grid: {nx} x {ny} cells, longitude {west} to {east}, "
      f"latitude {LAT_SOUTH} to {LAT_NORTH}")


# ---- 2. read the elevation data in the box ---------------------------------
def open_netcdf(path):
    if open(path, "rb").read(3) == b"CDF":           # old NetCDF3 file
        from scipy.io import netcdf_file
        return netcdf_file(path, "r", mmap=False)
    from netCDF4 import Dataset                       # NetCDF4 file
    return Dataset(path)

def find(nc, names):
    return next(nc.variables[n] for n in names if n in nc.variables)

nc = open_netcdf(ELEVATION_FILE)
lon = np.array(find(nc, ["lon", "longitude", "x"])[:], dtype=float)
lat = np.array(find(nc, ["lat", "latitude", "y"])[:], dtype=float)
elev_var = find(nc, ["z", "elevation", "Band1"])

# write every longitude as degrees east of the western edge (0..360), so a box
# across the date line becomes one continuous range
lon_e = west + (lon - west) % 360
ilon = np.where(lon_e < east)[0]
ilon = ilon[np.argsort(lon_e[ilon])]
ilon = ilon[np.r_[True, np.diff(lon_e[ilon]) > 1e-9]]   # drop the duplicate 180/-180
ilat = np.where((lat >= LAT_SOUTH) & (lat < LAT_NORTH))[0]

lo, hi = ilon.min(), ilon.max() + 1                      # read only what is needed
rows = elev_var[ilat.min():ilat.max() + 1, lo:hi]        # stored as (lat, lon)
elev = np.array(rows, dtype=np.float32)[:, ilon - lo]
lon_e, lat = lon_e[ilon], lat[ilat.min():ilat.max() + 1]
if lat[0] > lat[-1]:                                     # some files run north->south
    lat, elev = lat[::-1], elev[::-1]
print(f"read {elev.shape[1]} x {elev.shape[0]} elevation points")


# ---- 3. average the elevation points inside each model cell ----------------
i = np.minimum(((lon_e - west) / RESOLUTION).astype(int), nx - 1)
j = np.minimum(((lat - LAT_SOUTH) / RESOLUTION).astype(int), ny - 1)
cell = (j[:, None] * nx + i[None, :]).ravel()
total = np.bincount(cell, weights=elev.ravel().astype(float), minlength=nx * ny)
count = np.bincount(cell, minlength=nx * ny)
mean_elev = (total / count).reshape(ny, nx)


# ---- 4. turn elevation into an MITgcm depth field --------------------------
depth = np.minimum(mean_elev, 0.0)                    # negative = ocean, 0 = land
depth[depth > -MIN_DEPTH / 2] = 0.0                   # very shallow -> land
depth[(depth < 0) & (depth > -MIN_DEPTH)] = -MIN_DEPTH

# MITgcm wraps around (periodic) by default, so put walls on closed sides
if "S" not in OPEN_SIDES: depth[0, :] = 0
if "N" not in OPEN_SIDES: depth[-1, :] = 0
if "W" not in OPEN_SIDES: depth[:, 0] = 0
if "E" not in OPEN_SIDES: depth[:, -1] = 0

def keep_largest_ocean(d):
    """remove lakes and seas that are not connected to the main ocean"""
    labels, n = ndimage.label(d < 0)
    if n > 1:
        largest = np.argmax(np.bincount(labels.ravel())[1:]) + 1
        d[labels != largest] = 0
    return d

def remove_dead_ends(d):
    """remove one-cell-wide inlets (ocean cells with only one ocean neighbour)"""
    while True:
        wet = np.pad(d < 0, 1).astype(int)
        neighbours = wet[:-2, 1:-1] + wet[2:, 1:-1] + wet[1:-1, :-2] + wet[1:-1, 2:]
        dead = (d < 0) & (neighbours <= 1)
        # cells on an open side are fine
        if "S" in OPEN_SIDES: dead[0, :] = False
        if "N" in OPEN_SIDES: dead[-1, :] = False
        if "W" in OPEN_SIDES: dead[:, 0] = False
        if "E" in OPEN_SIDES: dead[:, -1] = False
        if not dead.any():
            return d
        d[dead] = 0

depth = keep_largest_ocean(depth)
depth = remove_dead_ends(depth)
depth = keep_largest_ocean(depth)
print(f"ocean cells: {(depth < 0).sum()} of {nx * ny}, "
      f"deepest: {-depth.min():.0f} m")


# ---- 5. vertical levels: thin at the surface, thicker at depth -------------
s = np.linspace(0, 1, NUM_LEVELS) ** 2
delR = 5.0 + s * (-depth.min() * 1.02 - 5.0 * NUM_LEVELS) / s.sum()
delR = np.round(delR, 1)


# ---- 6. write the files ----------------------------------------------------
depth.astype(">f4").tofile("bathy.bin")      # big-endian float32, as MITgcm reads

delR_text = ",\n          ".join(
    ", ".join(f"{v:.1f}" for v in delR[k:k + 8]) for k in range(0, NUM_LEVELS, 8))
with open("data_grid.txt", "w") as f:
    f.write(f"""# --- put these in input/data ---
 &PARM04
 usingSphericalPolarGrid = .TRUE.,
 xgOrigin = {west},
 ygOrigin = {LAT_SOUTH},
 delX = {nx}*{RESOLUTION},
 delY = {ny}*{RESOLUTION},
 delZ = {delR_text},
 &
 &PARM05
 bathyFile = 'bathy.bin',
 &

# --- and in code/SIZE.h: Nx = {nx}, Ny = {ny}, Nr = {NUM_LEVELS}
#     (sNx*nPx must equal {nx}, sNy*nPy must equal {ny})
""")
    if OPEN_SIDES:
        f.write("\n# --- open boundaries, in input/data.obcs ---\n &OBCS_PARM01\n")
        if "N" in OPEN_SIDES: f.write(f" OB_Jnorth = {nx}*-1,\n")
        if "S" in OPEN_SIDES: f.write(f" OB_Jsouth = {nx}*1,\n")
        if "E" in OPEN_SIDES: f.write(f" OB_Ieast = {ny}*-1,\n")
        if "W" in OPEN_SIDES: f.write(f" OB_Iwest = {ny}*1,\n")
        f.write(" &\n")


# ---- 7. picture to check the result ----------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

x_edges = west + RESOLUTION * np.arange(nx + 1)
y_edges = LAT_SOUTH + RESOLUTION * np.arange(ny + 1)
fig, ax = plt.subplots(figsize=(8, 7))
pic = ax.pcolormesh(x_edges, y_edges, np.ma.masked_equal(-depth, 0), cmap="viridis_r")
ax.set_facecolor("tan")                                # land
fig.colorbar(pic, label="depth (m)")
ax.set_xlabel("longitude (degrees east)")
ax.set_ylabel("latitude")
ax.set_title("bathy.bin")
fig.savefig("bathy_check.png", dpi=120)
print("wrote bathy.bin, data_grid.txt, bathy_check.png")
