"""
parquet_to_cells_zarr.py
------------------------
Convert a cell-boundaries parquet file into a cells.zarr.zip that is fully
compatible with XeniumCells reader.

Required parquet columns
    cell_id    – str, unique cell identifier
    vertex_x   – float32, x coordinate in Xenium micron space
    vertex_y   – float32, y coordinate in Xenium micron space
    (label_id  – ignored)

Usage
-----
    from parquet_to_cells_zarr import parquet_to_cells_zarr

    parquet_to_cells_zarr(
        parquet_path  = "cell_boundaries.parquet",
        output_path   = "cells.zarr.zip",
        boundary_set  = 1,   # 1 = cell boundaries (reader default), 0 = nucleus,
        scalefactor   = 0.2125,  # scaling applied to vertex coordinates (match Xenium mpp), or 0.3250 for Cell DIVE
    )

Reader compatibility notes
--------------------------
- Row index i in cell_summary  ==  row index i in polygon_sets/1/vertices
  The reader uses KDTree row indices as cell_id_int, so ordering is identity.
- Vertices are stored interleaved: [x0,y0, x1,y1, …, xN,yN, NaN,NaN, …]
  Slots beyond the actual vertex count are NaN-padded.
  The reader filters with np.isfinite() so padding is transparent.
- Only polygon_sets/1/vertices and cell_summary[:,0:2] are read by XeniumCells.
  All other arrays (cell_id, masks/, polygon_sets/N/cell_index, …) are written
  for spec-completeness but are not accessed by the reader.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd
import zarr
import zarr.storage
import numcodecs


# ── defaults matching the original Xenium cells.zarr.zip ─────────────────────

_BLOSC_DEFAULT = numcodecs.Blosc(cname="lz4", clevel=5, shuffle=numcodecs.Blosc.SHUFFLE)

# root-level metadata written verbatim (reader does not check these)
_ROOT_ATTRS = {
    "major_version": 6,
    "minor_version": 2,
    "name": "CellSegmentationDataset",
    "polygon_set_descriptions": ["DAPI-based nuclei segmentation", "Cell Segmentation"],
    "polygon_set_display_names": ["Nucleus boundaries", "Cell boundaries"],
    "polygon_set_names": ["nucleus", "cell"],
    "segmentation_methods": [
        "Segmented by boundary stain (None)",
        "Segmented by interior stain (None)",
        "Segmented by nucleus expansion of 5.0µm",
        "Segmented by nuclear stain (DAPI)",
        "Imported Cell Segmentation",
        "NA",
    ],
    "spatial_units": "microns",
}

_CELL_SUMMARY_COLS = [
    "cell_centroid_x", "cell_centroid_y", "cell_area",
    "nucleus_centroid_x", "nucleus_centroid_y", "nucleus_area",
    "z_level", "nucleus_count",
]
_CELL_SUMMARY_DESCS = [
    "Cell centroid in X", "Cell centroid in Y", "Cell area",
    "Nucleus centroid in X", "Nucleus centroid in Y", "Nucleus area",
    "z_level", "Nucleus count",
]


# ── helpers ───────────────────────────────────────────────────────────────────

def _chunk_rows(n: int, target_chunks: int = 16) -> int:
    """Row chunk size that divides n into ~target_chunks pieces."""
    return max(1, math.ceil(n / target_chunks))


def _create(group, name, data, chunks, compressor=_BLOSC_DEFAULT, attrs=None):
    """Create a zarr array and optionally attach attributes (v2 + v3 API)."""
    kw = dict(shape=data.shape, chunks=chunks, dtype=data.dtype,
              compressor=compressor, overwrite=True)
    # zarr 3.x prefers create_array; zarr 2.x prefers create_dataset
    creator = getattr(group, "create_array", None) or group.create_dataset
    arr = creator(name, **kw)
    arr[:] = data
    if attrs:
        arr.attrs.update(attrs)
    return arr


def _build_vertices(df: pd.DataFrame,
                    cell_order: np.ndarray,
                    codes: np.ndarray,
                    max_verts: int) -> np.ndarray:
    """
    Pack vertex coordinates into a [N, max_verts*2] float32 array.

    Layout per row: [x0,y0, x1,y1, …, xK,yK, NaN,NaN, …]
    Slots beyond the cell's actual vertex count are NaN-padded.
    The reader filters with np.isfinite() so padding is transparent.

    Vectorised: no Python loop over cells — runs in a few seconds even for
    600k+ cells with ~15M rows.
    """
    n_cells = len(cell_order)
    n_coords = max_verts * 2

    # Per-row vertex index within its cell (0, 1, 2, … per group)
    vertex_idx = df.groupby(codes, sort=False).cumcount().values

    # Drop any vertices beyond the cap
    mask       = vertex_idx < max_verts
    codes_m    = codes[mask]
    vidx_m     = vertex_idx[mask]
    vx         = df["vertex_x"].values[mask].astype(np.float32)
    vy         = df["vertex_y"].values[mask].astype(np.float32)

    # Write x and y directly into flat index of the output array
    vertices = np.full((n_cells, n_coords), np.nan, dtype=np.float32)
    ix = codes_m * n_coords + vidx_m * 2       # x slot
    iy = ix + 1                                 # y slot
    vertices.ravel()[ix] = vx
    vertices.ravel()[iy] = vy

    return vertices


# ── public API ────────────────────────────────────────────────────────────────

def parquet_to_cells_zarr(
    parquet_path: str | Path,
    output_path: str | Path = "cells.zarr.zip",
    boundary_set: int = 1,
    scalefactor: float = 0.2125,
    max_vertices: int | None = None,
    compressor: numcodecs.abc.Codec = _BLOSC_DEFAULT,
    extra_root_attrs: dict | None = None,
    cell_summary_extra: np.ndarray | None = None,
) -> None:
    """
    Convert a cell-boundaries parquet into a cells.zarr.zip.

    Parameters
    ----------
    parquet_path : str | Path
        Input parquet with columns: cell_id, vertex_x, vertex_y.
    output_path : str | Path
        Destination .zarr.zip (overwritten if exists).
    boundary_set : int
        Which polygon_set slot to write boundaries into.
        Must be 1 (cell) for the reader's default path; use 0 for nucleus.
    max_vertices : int or None
        Cap on vertices per cell. If None, uses the maximum found in data.
        Must be ≥ actual max. Original Xenium files use 25.
    compressor : numcodecs codec
        Blosc (default) matches original files. Pass None for no compression.
    extra_root_attrs : dict or None
        Merged into root-level attrs (overrides defaults).
    cell_summary_extra : np.ndarray or None
        Shape [N, 6] float64 for columns 2-7 of cell_summary
        (cell_area, nucleus_x, nucleus_y, nucleus_area, z_level, nucleus_count).
        If None, those columns are zero-filled.
    """
    parquet_path = Path(parquet_path)
    output_path  = Path(output_path)

    # ── 1. Load and validate ──────────────────────────────────────────────────
    print(f"Reading {parquet_path} …")
    df = pd.read_parquet(parquet_path)

    required = {"cell_id", "vertex_x", "vertex_y"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"Parquet is missing required columns: {missing}")

    df["vertex_x"] = df["vertex_x"].astype(np.float32) * scalefactor
    df["vertex_y"] = df["vertex_y"].astype(np.float32) * scalefactor

    print(scalefactor, "scaling applied to vertex coordinates")

    # ── 2. Factorize cell_id → integer codes in one vectorised pass ───────────
    # codes[i] = integer index of df["cell_id"].iloc[i] in cell_order.
    # Row index i in zarr == position i in cell_order.
    # The reader uses KDTree row indices as cell_id_int, so this ordering
    # becomes the ground truth for cell identity.
    codes, cell_order = pd.factorize(df["cell_id"], sort=False)
    n_cells = len(cell_order)
    print(f"  {n_cells:,} unique cells")

    # ── 3. Vertex cap (vectorised via bincount on integer codes) ──────────────
    counts   = np.bincount(codes, minlength=n_cells).astype(np.int32)
    data_max = int(counts.max())
    print(f"  Vertices per cell — min: {counts.min()}  "
          f"max: {data_max}  mean: {counts.mean():.1f}")

    if max_vertices is None:
        max_vertices = data_max
    elif max_vertices < data_max:
        raise ValueError(
            f"max_vertices={max_vertices} < data maximum {data_max}. "
            "Increase max_vertices or set to None to auto-detect."
        )
    print(f"  Vertex slots per cell: {max_vertices}  "
          f"(zarr dim1 = {max_vertices * 2} floats)")

    # ── 4. Build core arrays (fully vectorised — no Python loops) ─────────────
    print("Building vertex array …")
    vertices = _build_vertices(df, cell_order, codes, max_vertices)

    # Centroids via bincount weighted sum — no Python loop
    print("Computing centroids …")
    vx_vals = df["vertex_x"].values.astype(np.float64)
    vy_vals = df["vertex_y"].values.astype(np.float64)
    cx = np.bincount(codes, weights=vx_vals, minlength=n_cells) / counts
    cy = np.bincount(codes, weights=vy_vals, minlength=n_cells) / counts

    # cell_summary: [N, 8]  cols = [cx, cy, area, ncx, ncy, n_area, z, n_count]
    cell_summary = np.zeros((n_cells, 8), dtype=np.float64)
    cell_summary[:, 0] = cx
    cell_summary[:, 1] = cy
    if cell_summary_extra is not None:
        extra = np.asarray(cell_summary_extra, dtype=np.float64)
        if extra.shape != (n_cells, 6):
            raise ValueError(
                f"cell_summary_extra must be shape ({n_cells}, 6), got {extra.shape}"
            )
        cell_summary[:, 2:] = extra

    # cell_id array: [N, 2] uint32 — not read by reader
    cell_id_arr = np.zeros((n_cells, 2), dtype=np.uint32)
    cell_id_arr[:, 0] = np.arange(n_cells, dtype=np.uint32)

    # num_vertices: clipped to max_vertices cap
    num_vertices_arr = np.minimum(counts, max_vertices)

    # method: zero-filled (not read by reader)
    method_arr = np.zeros(n_cells, dtype=np.uint32)

    # cell_index: 0..N-1 (not read by reader)
    cell_index_arr = np.arange(n_cells, dtype=np.uint32)

    # ── 5. Chunk sizes matching original Xenium file pattern ─────────────────
    row_chunk    = _chunk_rows(n_cells, target_chunks=16)
    cs_col_chunk = 1                         # cell_summary: one chunk per col
    vx_col_chunk = min(7, max_vertices * 2)  # vertices: 7 cols (original pattern)

    # ── 6. Write zarr store ───────────────────────────────────────────────────
    print(f"Writing {output_path} …")
    store = zarr.storage.ZipStore(str(output_path), mode="w")
    try:
        root = zarr.open_group(store, mode="w", zarr_format=2)
    except TypeError:
        root = zarr.open_group(store, mode="w")   # zarr 2.x: v2 is the default

    # root attrs
    root_attrs = dict(_ROOT_ATTRS)
    root_attrs["number_cells"] = n_cells
    if extra_root_attrs:
        root_attrs.update(extra_root_attrs)
    root.attrs.update(root_attrs)

    # cell_id  [N, 2]
    _create(root, "cell_id", cell_id_arr,
            chunks=(row_chunk, 1))

    # cell_summary  [N, 8]
    _create(root, "cell_summary", cell_summary,
            chunks=(n_cells, cs_col_chunk),
            compressor=compressor,
            attrs={
                "column_names":        _CELL_SUMMARY_COLS,
                "column_descriptions": _CELL_SUMMARY_DESCS,
            })

    # masks/  (not read by reader — written as empty stubs for spec compliance)
    masks = root.require_group("masks")
    stub_shape  = (4, 4)          # minimal non-zero shape
    stub_chunks = (4, 4)
    for slot in (0, 1):
        _create(masks, str(slot),
                np.zeros(stub_shape, dtype=np.uint32),
                chunks=stub_chunks, compressor=compressor)
    _create(masks, "homogeneous_transform",
            np.eye(4, dtype=np.float32),
            chunks=(4, 4), compressor=compressor)

    # polygon_sets/
    psets = root.require_group("polygon_sets")

    # Write both slots 0 and 1; the reader only ever reads slot 1 (boundary_set=1).
    # Slot 0 (nucleus) gets zero-filled vertices as a stub.
    for slot in (0, 1):
        grp = psets.require_group(str(slot))

        if slot == boundary_set:
            verts_data = vertices
        else:
            # stub: zeros (reader never reads this slot)
            verts_data = np.zeros((n_cells, max_vertices * 2), dtype=np.float32)

        _create(grp, "vertices",    verts_data,
                chunks=(row_chunk, vx_col_chunk), compressor=compressor)
        _create(grp, "num_vertices", num_vertices_arr,
                chunks=(row_chunk,), compressor=compressor)
        _create(grp, "cell_index",   cell_index_arr,
                chunks=(row_chunk,), compressor=compressor)
        _create(grp, "method",       method_arr,
                chunks=(row_chunk,), compressor=compressor)

    store.close()
    print(f"Done — {output_path.stat().st_size / 1e6:.1f} MB")
