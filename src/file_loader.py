"""NetCDF input loading and grid/PFT normalization utilities."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, Sequence, Tuple

import netCDF4 as nc
import numpy as np


@dataclass
class GridMeta:
    lat: np.ndarray
    lon: np.ndarray
    lat_desc: bool
    lon_0_360: bool
    lat_name: str
    lon_name: str
    time_name: Optional[str]
    time_values: Optional[np.ndarray]
    internal_lat: np.ndarray
    internal_lon: np.ndarray


@dataclass
class DatasetInfo:
    path: Path
    variables: tuple[str, ...]
    dimensions: dict[str, int]
    lat_name: str
    lon_name: str
    time_name: Optional[str]
    lat: np.ndarray
    lon: np.ndarray
    time_values: Optional[np.ndarray]


@dataclass
class PFTData:
    """Static or dynamic PFT fractions aligned to the model grid."""

    values: np.ndarray  # (time, lat, lon, pft), time length is 1 for static
    years: np.ndarray   # calendar years corresponding to values[:, ...]
    dynamic: bool
    variable_name: str
    source_mode: str

    def get_cell(self, calendar_year: int, i: int, j: int) -> np.ndarray:
        if not self.dynamic:
            return np.asarray(self.values[0, i, j], dtype=np.float64)
        pos = int(np.searchsorted(self.years, int(calendar_year)))
        if pos < len(self.years) and int(self.years[pos]) == int(calendar_year):
            idx = pos
        elif pos <= 0:
            idx = 0
        elif pos >= len(self.years):
            idx = len(self.years) - 1
        else:
            left = pos - 1
            right = pos
            idx = left if abs(self.years[left] - calendar_year) <= abs(self.years[right] - calendar_year) else right
        return np.asarray(self.values[idx, i, j], dtype=np.float64)


class _MemoryDim:
    def __init__(self, size: int):
        self.size = int(size)


class _MemoryDataset:
    """Small in-memory object exposing netCDF4-like variables/dimensions."""

    def __init__(self, variables: dict, dimensions: dict):
        self.variables = variables
        self.dimensions = {key: _MemoryDim(value) for key, value in dimensions.items()}


def _guess_lat_lon_names(ds: nc.Dataset) -> Tuple[str, str]:
    lat_name = next((n for n in ("lat", "latitude", "y") if n in ds.variables), None)
    lon_name = next((n for n in ("lon", "longitude", "x") if n in ds.variables), None)
    if lat_name is None or lon_name is None:
        raise KeyError("Cannot find lat/lon coordinate variables.")
    return lat_name, lon_name


def _guess_time_name(ds: nc.Dataset) -> Optional[str]:
    for name in ("time", "year", "time_counter"):
        if name in ds.variables:
            return name
    for name, variable in ds.variables.items():
        if getattr(variable, "standard_name", "").lower() == "time":
            return name
    return None


def _as_1d_coord(variable, coord_type: str) -> np.ndarray:
    arr = np.asarray(variable[:], dtype=float)
    if arr.ndim == 1:
        out = arr
    elif arr.ndim == 2:
        if coord_type == "lat":
            col0 = arr[:, 0]
            row0 = arr[0, :]
            if np.allclose(arr, col0[:, None], equal_nan=True):
                out = col0
            elif np.allclose(arr, row0[None, :], equal_nan=True):
                out = row0
            else:
                raise ValueError(f"Irregular two-dimensional latitude coordinate: {arr.shape}")
        else:
            row0 = arr[0, :]
            col0 = arr[:, 0]
            if np.allclose(arr, row0[None, :], equal_nan=True):
                out = row0
            elif np.allclose(arr, col0[:, None], equal_nan=True):
                out = col0
            else:
                raise ValueError(f"Irregular two-dimensional longitude coordinate: {arr.shape}")
    else:
        raise ValueError(f"Unsupported {coord_type} coordinate ndim={arr.ndim}")
    return np.asarray(out, dtype=float).ravel()


def _is_descending(values: np.ndarray) -> bool:
    return bool(values.size >= 2 and np.all(np.diff(values) < 0))


def _is_lon_0_360(values: np.ndarray) -> bool:
    finite = values[np.isfinite(values)]
    return bool(finite.size and finite.min() >= 0 and finite.max() > 180)


def _wrap_lon_180(values: np.ndarray) -> np.ndarray:
    return ((np.asarray(values, dtype=float) + 180.0) % 360.0) - 180.0


def _optimize_indices_to_slice(indices: np.ndarray):
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size == 0:
        return slice(0, 0)
    if indices.size == 1:
        return slice(int(indices[0]), int(indices[0]) + 1)
    step = int(indices[1] - indices[0])
    if np.all(np.diff(indices) == step):
        if step > 0:
            return slice(int(indices[0]), int(indices[-1]) + 1, step)
        stop = int(indices[-1]) - 1
        return slice(int(indices[0]), None if stop < 0 else stop, step)
    return indices


def _coord_resolution(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float).ravel()
    unique = np.unique(np.round(values[np.isfinite(values)], 10))
    if unique.size < 2:
        return np.nan
    return abs(float(np.median(np.diff(np.sort(unique)))))


def _nearest_indices(source: np.ndarray, target: np.ndarray, *, periodic: bool) -> np.ndarray:
    source = np.asarray(source, dtype=float).ravel()
    target = np.asarray(target, dtype=float).ravel()
    result = np.empty(target.size, dtype=np.int64)
    for k, value in enumerate(target):
        if periodic:
            distance = np.abs(((source - value + 180.0) % 360.0) - 180.0)
        else:
            distance = np.abs(source - value)
        result[k] = int(np.nanargmin(distance))
    return result


def calendarize_time_values(
    raw_values: Optional[np.ndarray],
    *,
    length: int,
    input_base_year: int,
    encoding: str,
) -> np.ndarray:
    """Convert a NetCDF time axis or positional indices to calendar years."""
    mode = str(encoding).lower()
    if raw_values is None:
        raw = np.arange(length, dtype=int)
    else:
        raw = np.asarray(raw_values).reshape(-1)
        if raw.size != length:
            raise ValueError(f"Time coordinate length {raw.size} != dimension length {length}")

    if mode == "calendar_year":
        return np.rint(raw).astype(int)
    if mode == "index":
        return int(input_base_year) + np.rint(raw).astype(int)
    if mode != "auto":
        raise ValueError("time_encoding must be index, calendar_year, or auto")

    finite = raw[np.isfinite(raw)]
    if finite.size and np.nanmin(finite) >= 500:
        return np.rint(raw).astype(int)
    return int(input_base_year) + np.rint(raw).astype(int)


def resolve_year_index(
    info: DatasetInfo,
    calendar_year: int,
    *,
    input_base_year: int,
    time_encoding: str,
) -> int:
    time_len = info.dimensions.get(info.time_name, 1) if info.time_name else 1
    years = calendarize_time_values(
        info.time_values,
        length=time_len,
        input_base_year=input_base_year,
        encoding=time_encoding,
    )
    matches = np.flatnonzero(years == int(calendar_year))
    if matches.size == 0:
        raise IndexError(
            f"Calendar year {calendar_year} is not present in {info.path}; "
            f"available={int(years.min())}..{int(years.max())}."
        )
    return int(matches[0])


class FileLoader:
    def inspect_dataset(self, path: str | Path) -> DatasetInfo:
        resolved = Path(path).resolve()
        with nc.Dataset(resolved, "r") as ds:
            lat_name, lon_name = _guess_lat_lon_names(ds)
            time_name = _guess_time_name(ds)
            time_values = (
                np.asarray(ds.variables[time_name][:])
                if time_name is not None
                else None
            )
            return DatasetInfo(
                path=resolved,
                variables=tuple(ds.variables.keys()),
                dimensions={name: len(dim) for name, dim in ds.dimensions.items()},
                lat_name=lat_name,
                lon_name=lon_name,
                time_name=time_name,
                lat=_as_1d_coord(ds.variables[lat_name], "lat"),
                lon=_as_1d_coord(ds.variables[lon_name], "lon"),
                time_values=time_values,
            )

    def load_luh2_dataset(
        self,
        path: str | Path,
        *,
        lat_slice: slice | None = None,
        lon_slice: slice | None = None,
        time_slice: slice | None = None,
        variable_names: Optional[Iterable[str]] = None,
    ) -> Tuple[_MemoryDataset, GridMeta]:
        """Load only requested variables/slices and close the NetCDF immediately."""
        requested = None if variable_names is None else set(variable_names)
        resolved = Path(path).resolve()

        ds = nc.Dataset(resolved, "r")
        try:
            lat_name, lon_name = _guess_lat_lon_names(ds)
            time_name = _guess_time_name(ds)
            lat_full = _as_1d_coord(ds.variables[lat_name], "lat")
            lon_full = _as_1d_coord(ds.variables[lon_name], "lon")
            is_desc = _is_descending(lat_full)
            lon_0_360 = _is_lon_0_360(lon_full)
            nlat = lat_full.size
            nlon = lon_full.size

            lat_slice = lat_slice or slice(0, nlat)
            lon_slice = lon_slice or slice(0, nlon)

            if is_desc:
                i0 = 0 if lat_slice.start is None else int(lat_slice.start)
                i1 = nlat if lat_slice.stop is None else int(lat_slice.stop)
                lat_indices = np.arange(nlat)[slice(nlat - i1, nlat - i0)]
            else:
                lat_indices = np.arange(nlat)[lat_slice]
            lat_raw = lat_full[lat_indices]
            lat_meta = lat_raw[::-1] if is_desc else lat_raw

            lon_res = abs(float(np.median(np.diff(lon_full))))
            nlon_internal = int(round(360.0 / lon_res))
            lon_internal_full = np.linspace(
                -180.0 + 0.5 * lon_res,
                180.0 - 0.5 * lon_res,
                nlon_internal,
            )
            target_lon = lon_internal_full[lon_slice]
            lon_wrapped = _wrap_lon_180(lon_full)
            lon_indices = _nearest_indices(lon_wrapped, target_lon, periodic=True)
            lon_raw = lon_full[lon_indices]

            lat_res = abs(float(np.median(np.diff(lat_full))))
            internal_lat_full = np.linspace(
                -90.0 + 0.5 * lat_res,
                90.0 - 0.5 * lat_res,
                nlat,
            )

            raw_time_values = (
                np.asarray(ds.variables[time_name][:])
                if time_name is not None
                else None
            )
            if time_name is not None and time_slice is not None:
                selected_time_values = raw_time_values[time_slice]
            else:
                selected_time_values = raw_time_values

            meta = GridMeta(
                lat=lat_meta,
                lon=lon_raw,
                lat_desc=is_desc,
                lon_0_360=lon_0_360,
                lat_name=lat_name,
                lon_name=lon_name,
                time_name=time_name,
                time_values=selected_time_values,
                internal_lat=internal_lat_full[lat_slice],
                internal_lon=target_lon,
            )

            dimensions: dict[str, int] = {}
            for dim_name, dim in ds.dimensions.items():
                if dim_name == lat_name:
                    dimensions[dim_name] = len(lat_indices)
                elif dim_name == lon_name:
                    dimensions[dim_name] = len(lon_indices)
                elif dim_name == time_name and time_slice is not None:
                    dimensions[dim_name] = len(np.arange(len(dim))[time_slice])
                else:
                    dimensions[dim_name] = len(dim)

            variables: dict[str, np.ndarray] = {
                lat_name: np.asarray(lat_raw, dtype=np.float32),
                lon_name: np.asarray(lon_raw, dtype=np.float32),
            }
            if time_name is not None and selected_time_values is not None:
                variables[time_name] = np.asarray(selected_time_values)

            lat_indexer = _optimize_indices_to_slice(lat_indices)
            lon_indexer = _optimize_indices_to_slice(lon_indices)
            for name, variable in ds.variables.items():
                if name in (lat_name, lon_name, time_name):
                    continue
                if requested is not None and name not in requested:
                    continue
                index = []
                for dim_name in variable.dimensions:
                    if dim_name == lat_name:
                        index.append(lat_indexer)
                    elif dim_name == lon_name:
                        index.append(lon_indexer)
                    elif dim_name == time_name and time_slice is not None:
                        index.append(time_slice)
                    else:
                        index.append(slice(None))
                try:
                    arr = variable[tuple(index)]
                except Exception as exc:
                    raise RuntimeError(
                        f"Failed reading {name!r} from {resolved} with index={index}"
                    ) from exc
                if hasattr(arr, "filled"):
                    arr = arr.filled(np.nan)
                variables[name] = np.asarray(arr, dtype=np.float32)

            print(
                f"[LOAD] {resolved}: lat={len(lat_indices)}, lon={len(lon_indices)}, "
                f"time={dimensions.get(time_name, 1) if time_name else 1}, "
                f"data_vars={len(variables) - 2 - int(time_name is not None)}",
                flush=True,
            )
            return _MemoryDataset(variables, dimensions), meta
        finally:
            ds.close()

    def load_pft_map(
        self,
        pft_path: str | Path,
        *,
        pft_var: Optional[str],
        target_lat_asc: np.ndarray,
        target_lon_180: np.ndarray,
        target_lon_res: float,
        expected_n_pft: int,
        model_years: Sequence[int],
        pft_base_year: int,
        time_encoding: str = "auto",
        min_year_policy: str = "clip",
        max_year_policy: str = "clip",
    ) -> PFTData:
        """Load static/dynamic dominant or fractional PFT data.

        Supported layouts are ``lat,lon``; ``pft,lat,lon``;
        ``time,lat,lon``; and ``time,pft,lat,lon`` in any dimension order.
        """
        resolved = Path(pft_path).resolve()
        target_lat = np.asarray(target_lat_asc, dtype=float).ravel()
        target_lon = _wrap_lon_180(np.asarray(target_lon_180, dtype=float).ravel())
        model_years = np.asarray(model_years, dtype=int)

        pft_aliases = {"pft", "pfts", "npft", "veget", "veg", "vegtype", "vegetation", "class", "classes", "type", "types"}
        time_aliases = {"time", "year", "time_counter"}
        lat_aliases = {"lat", "latitude", "y"}
        lon_aliases = {"lon", "longitude", "x"}

        ds = nc.Dataset(resolved, "r")
        try:
            lat_name, lon_name = _guess_lat_lon_names(ds)
            time_name = _guess_time_name(ds)

            if pft_var is None:
                preferred = (
                    "pft_fraction", "pft_frac", "pft_fractions",
                    "vegetation_fraction", "vegetfrac", "maxvegetfrac",
                    "pft", "PFT", "dominant_pft", "pft_map", "pft_dominant",
                )
                chosen = next((name for name in preferred if name in ds.variables), None)
                if chosen is None:
                    candidates = []
                    for name, candidate in ds.variables.items():
                        dims = {dim.lower() for dim in candidate.dimensions}
                        if dims & lat_aliases and dims & lon_aliases:
                            candidates.append(name)
                    if len(candidates) != 1:
                        raise KeyError(
                            "Cannot infer PFT variable. Set inputs.pft_variable in the run YAML. "
                            f"Candidates={candidates}"
                        )
                    chosen = candidates[0]
            else:
                chosen = str(pft_var)
                if chosen not in ds.variables:
                    raise KeyError(
                        f"PFT variable {chosen!r} not found in {resolved}; "
                        f"available={list(ds.variables)}"
                    )

            variable = ds.variables[chosen]
            dims = list(variable.dimensions)
            lower = [dim.lower() for dim in dims]
            lat_axis = next((i for i, dim in enumerate(lower) if dim in lat_aliases), None)
            lon_axis = next((i for i, dim in enumerate(lower) if dim in lon_aliases), None)
            time_axis = next((i for i, dim in enumerate(lower) if dim in time_aliases), None)
            if lat_axis is None or lon_axis is None:
                raise ValueError(f"PFT variable {chosen!r} lacks lat/lon dimensions: {dims}")

            src_lat = _as_1d_coord(ds.variables[lat_name], "lat")
            src_lon = _wrap_lon_180(_as_1d_coord(ds.variables[lon_name], "lon"))
            src_lat_res = _coord_resolution(src_lat)
            src_lon_res = _coord_resolution(src_lon)
            tgt_lat_res = _coord_resolution(target_lat)
            tgt_lon_res = float(target_lon_res) if target_lon_res > 0 else _coord_resolution(target_lon)
            if (
                np.isfinite(src_lat_res) and np.isfinite(tgt_lat_res)
                and src_lat_res < tgt_lat_res * (1 - 1e-6)
            ) or (
                np.isfinite(src_lon_res) and np.isfinite(tgt_lon_res)
                and src_lon_res < tgt_lon_res * (1 - 1e-6)
            ):
                raise ValueError(
                    "PFT map is finer than the model grid. Pre-aggregate fractions "
                    "with area weighting before simulation."
                )

            lat_indices = _nearest_indices(src_lat, target_lat, periodic=False)
            lon_indices = _nearest_indices(src_lon, target_lon, periodic=True)
            lat_indexer = _optimize_indices_to_slice(lat_indices)
            lon_indexer = _optimize_indices_to_slice(lon_indices)

            pft_axis = next((i for i, dim in enumerate(lower) if dim in pft_aliases), None)
            if pft_axis is None:
                candidates = [
                    axis for axis, size in enumerate(variable.shape)
                    if axis not in (lat_axis, lon_axis, time_axis)
                    and int(size) == int(expected_n_pft)
                ]
                if len(candidates) == 1:
                    pft_axis = candidates[0]

            dynamic = (
                time_axis is not None and int(variable.shape[time_axis]) > 1
            )
            if dynamic:
                raw_time = np.asarray(ds.variables[time_name][:]) if time_name else None
                source_years = calendarize_time_values(
                    raw_time,
                    length=variable.shape[time_axis],
                    input_base_year=pft_base_year,
                    encoding=time_encoding,
                )
                requested_indices = []
                selected_years = []
                for year in model_years:
                    exact = np.flatnonzero(source_years == year)
                    if exact.size:
                        idx = int(exact[0])
                    elif year < source_years.min() and min_year_policy == "clip":
                        idx = int(np.argmin(source_years))
                    elif year > source_years.max() and max_year_policy == "clip":
                        idx = int(np.argmax(source_years))
                    else:
                        raise IndexError(
                            f"PFT year {year} unavailable; source range "
                            f"{source_years.min()}..{source_years.max()}"
                        )
                    requested_indices.append(idx)
                    selected_years.append(int(year))
                time_indexer = np.asarray(requested_indices, dtype=np.int64)
            else:
                time_indexer = 0 if time_axis is not None else None
                selected_years = [int(model_years[0])]

            index = []
            retained_dims = []
            for axis, dim in enumerate(dims):
                if axis == lat_axis:
                    index.append(lat_indexer)
                    retained_dims.append(dim)
                elif axis == lon_axis:
                    index.append(lon_indexer)
                    retained_dims.append(dim)
                elif axis == time_axis:
                    index.append(time_indexer)
                    if dynamic:
                        retained_dims.append(dim)
                else:
                    index.append(slice(None))
                    retained_dims.append(dim)

            arr = variable[tuple(index)]
            if hasattr(arr, "filled"):
                arr = arr.filled(np.nan)
            arr = np.asarray(arr)
        finally:
            ds.close()

        retained_lower = [dim.lower() for dim in retained_dims]
        lat_axis_now = next(i for i, dim in enumerate(retained_lower) if dim in lat_aliases)
        lon_axis_now = next(i for i, dim in enumerate(retained_lower) if dim in lon_aliases)
        time_axis_now = next((i for i, dim in enumerate(retained_lower) if dim in time_aliases), None)
        pft_axis_now = next((i for i, dim in enumerate(retained_lower) if dim in pft_aliases), None)
        if pft_axis_now is None:
            candidates = [
                axis for axis, size in enumerate(arr.shape)
                if axis not in (lat_axis_now, lon_axis_now, time_axis_now)
                and int(size) == int(expected_n_pft)
            ]
            if len(candidates) == 1:
                pft_axis_now = candidates[0]

        nlat = target_lat.size
        nlon = target_lon.size
        if pft_axis_now is None:
            # Dominant map, static 2-D or dynamic 3-D.
            if time_axis_now is None:
                dominant = np.moveaxis(arr, (lat_axis_now, lon_axis_now), (0, 1))[None, ...]
            else:
                dominant = np.moveaxis(
                    arr,
                    (time_axis_now, lat_axis_now, lon_axis_now),
                    (0, 1, 2),
                )
            if dominant.shape[1:] != (nlat, nlon):
                raise ValueError(f"Dominant PFT shape mismatch: {dominant.shape}")
            rounded = np.rint(dominant)
            valid = (
                np.isfinite(dominant)
                & (np.abs(dominant - rounded) <= 1e-5)
                & (rounded >= 1)
                & (rounded <= int(expected_n_pft))
            )
            fractions = np.zeros(
                (dominant.shape[0], nlat, nlon, int(expected_n_pft)),
                dtype=np.float32,
            )
            tt, ii, jj = np.nonzero(valid)
            fractions[tt, ii, jj, rounded[tt, ii, jj].astype(int) - 1] = 1.0
            source_mode = "dynamic_dominant" if dynamic else "static_dominant"
        else:
            if time_axis_now is None:
                fractions = np.moveaxis(
                    arr,
                    (lat_axis_now, lon_axis_now, pft_axis_now),
                    (0, 1, 2),
                )[None, ...]
            else:
                fractions = np.moveaxis(
                    arr,
                    (time_axis_now, lat_axis_now, lon_axis_now, pft_axis_now),
                    (0, 1, 2, 3),
                )
            fractions = np.asarray(fractions, dtype=np.float32)
            if fractions.shape[1:] != (nlat, nlon, int(expected_n_pft)):
                raise ValueError(
                    f"Fractional PFT shape mismatch: {fractions.shape}; expected "
                    f"(time,{nlat},{nlon},{expected_n_pft})"
                )
            fractions = np.where(np.isfinite(fractions), fractions, 0.0)
            fractions[fractions < 0.0] = 0.0
            sums = fractions.sum(axis=3, keepdims=True, dtype=np.float64)
            fractions = np.divide(
                fractions,
                sums,
                out=np.zeros_like(fractions),
                where=sums > 0.0,
            )
            source_mode = "dynamic_fraction" if dynamic else "static_fraction"

        zero_cells = int(np.count_nonzero(fractions.sum(axis=3) <= 0.0))
        print(
            f"[PFT] {resolved}: var={chosen}, mode={source_mode}, "
            f"shape={fractions.shape}, zero_cells={zero_cells}",
            flush=True,
        )
        return PFTData(
            values=fractions,
            years=np.asarray(selected_years, dtype=int),
            dynamic=bool(dynamic),
            variable_name=chosen,
            source_mode=source_mode,
        )
