#!/usr/bin/env python3
# ******************************************************************************
#
# TOMATO - TOMogram Annotation TOol
# Copyright (C) 2026 Simon J. Alkema
#
# Author: Simon J. Alkema with contributions from Euan W. Pyle, Higor Rosa and Anastasiia Babenko
# EMBL Imaging Centre, Heidelberg, Germany - Mattei Lab
# contact alkemasimon@gmail.com for questions
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
#
# This complete copyright notice must be included in any revised version of the
# source code. Additional authorship citations may be added, but existing
# author citations must be preserved.
#
# ******************************************************************************
"""
Tomogram annotator that inventories tomogram contents, with an in-GUI orthogonal XYZ viewer.

Each .mrc is loaded with mrcfile and shown as a quad view: XY (main,
bottom-left), XZ (top-left), YZ (bottom-right), and a position readout
in the remaining corner.

Annotation:
    Mito / ER / Golgi are always present. Extra flags can be added
    with --flags, or live in the GUI via the "Add flag" box.
    Every flag starts as an on/off toggle and can be promoted to an
    integer counter from its corner control, or demoted back again. That
    choice is a session-wide UI preference, not per-tomogram data.
    If the output CSV already exists it is loaded on startup, so a run
    can be resumed. --previous_csv continues one specific CSV instead:
    its flags and annotations are loaded and annotation is written back
    to it, unless --output names somewhere else.
Outputs:
    <output>.csv  full table: a count per flag, plus quality and comment
                  (a plain toggle is just 0/1)
    <output>.txt  compact summary, one line per tomogram:
                  TS_0001.mrc,Good,Mito,ER,Shell:3
                  (bare flag name when count == 1, name:N when higher)
Usage:
    python TOMATO.py annotate --input /path/to/tomograms
                              [--flags Lysosome Nucleus ...]
                              [--output annotations.csv]
                              [--suffix .mrc]
                              [--prefix TS_]
                              [--previous_csv prior_annotations.csv]
                              [--pixel_size 3.42]
    python TOMATO.py link --csv annotations.csv
                          --input /path/to/tomograms
                          --object Mito (one or more options available)
Navigation:
    Click a panel to move the crosshair there. The wheel steps through
    the axis normal to that panel (XY -> Z, XZ -> Y, YZ -> X); the slider
    below XY also moves Z. The crosshair resets to the volume centre
    whenever a new tomogram loads.
Zoom / pan:
    Ctrl + wheel zooms about the pointer; right- or middle-button drag
    pans. Both are synced across all three panels, since every pair
    shares an axis (XY/XZ share X, XY/YZ share Y, XZ/YZ share Z).
Measure:
    Left-button drag on any panel draws a ruler and reports its length in
    the readout corner, in both angstrom and nanometres. It uses the Å/px
    value shown there, which defaults to the MRC header's voxel size and
    can be edited. Display-only: nothing is saved.
Window sizing:
    The window is fitted to the screen, so nothing is hidden or scrolled
    off. The quad view takes whatever the annotation and list columns
    don't need, down to a floor of MIN_XY pixels on the XY panel's
    longest side; the flag list scrolls rather than growing the window.
    Resizing or maximizing afterwards grows the quad view. The window can
    be shrunk again, but not past the size it was originally fitted to.
Performance:
    Tomograms ahead of the current one are loaded and pre-rendered on
    background threads, so "Next" is instant. Running out of memory
    reduces how far ahead that reaches — down to fully synchronous
    loading — instead of crashing.
"""

# ══════════════════════════════════════════════════════════════════════════════
# READING THIS FILE
# The file is in 6 numbered SECTIONs, listed here in the order they appear.
# Search for "SECTION " to jump between them.
#
#   SECTION 1  Configuration
#              Flag names, quality grades, colours, panel size limits.
#              Safe to edit, and the most likely thing you want.
#   SECTION 2  Reading and writing the annotation files
#              The output CSV and TXT. Edit for a different format.
#   SECTION 3  Turning a tomogram into images
#              Contrast, downsampling, pre-rendering. Display only —
#              nothing here changes what gets saved.
#   SECTION 4  The annotator GUI  (class Annotator)
#              Long, but split into 8 PARTs — see the map just below
#              "class Annotator:". You rarely need to go in here.
#   SECTION 5  The 'link' subcommand
#   SECTION 6  Command line and entry point
#
# Headings get lighter as you go deeper, so you can tell at a glance how
# far in you are:  ═ box = SECTION,  ─ box = PART,  "# ── … ──" = a group
# of related functions,  plain "#" = a note about the lines beneath it.
# ══════════════════════════════════════════════════════════════════════════════

import argparse
import csv
import math
import os
import re
import sys
import threading
from pathlib import Path
try:
    import tkinter as tk
except ImportError:
    print("This script requires the 'tkinter' package.")
    print("tkinter not installed - for EMBL ml Tkinter/3.13.5-GCCcore-14.3.0")
    sys.exit(1)


try:
    import mrcfile
except ImportError:
    print("This script requires the 'mrcfile' package.")
    print("mrcfile not installed - for EMBL ml mrcfile")
    sys.exit(1)


try:
    import numpy as np
except ImportError:
    print("This script requires numpy.")
    print("numpy not installed - for EMBL ml SciPy-bundle/2025.07-gfbf-2025b")
    sys.exit(1)


try:
    from PIL import Image, ImageTk
except ImportError:
    print("This script requires the 'Pillow' package.")
    print("Pillow not installed - for EMBL ml Pillow/11.3.0-GCCcore-14.3.0")
    sys.exit(1)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 1 — Configuration
# Flag names, quality grades, colours and panel size limits.
# The part you are most likely to want to change.
# ══════════════════════════════════════════════════════════════════════════════
BASE_FLAGS      = ["Mito", "ER", "Golgi"]
GRADES          = ["FigureQuality", "Good", "Ok", "Bad"]
DEFAULT_SUFFIX  = ".mrc"
DEFAULT_OUTPUT  = "annotations.csv"
# Bounds on the XY panel's longest side, in pixels — not on the window.
# The quad view is fitted to whatever the other columns leave, then
# clamped to this range (see _fit_scale).
# MIN_XY: below this you can't make anything out. Applied last, so it can
# overflow the budget rather than go under it.
# MAX_XY: every Z slice is prerendered and cached at the panel's size, so
# an uncapped panel risks running out of memory on a big monitor. It also
# sets how large the window opens, and is only the startup ceiling —
# _apply_quad_resize raises it when the window is grown.
MAX_XY          = 750
MIN_XY          = 320
# How many tomograms ahead of the current one are prefetched into memory.
# At runtime this only ever steps DOWN (2 -> 1 -> 0), in reaction to a
# MemoryError. There is deliberately no upfront free-memory check: with
# several threads allocating at once the reading is outdated the moment
# you take it, so reacting to the real failure is simpler and safer.
PREFETCH_AHEAD_DEFAULT = 2
# Padding around the columns (the tomato "skin" border), in pixels.
WINDOW_PAD      = 24
FLAG_COLORS = [
    ("#791F1F", "#F7C1C1"),
    ("#1B4F8A", "#CFE3F7"),
    ("#0F6E56", "#CDEEE1"),
    ("#854F0B", "#FAEEDA"),
    ("#3C3489", "#E2DFFB"),
    ("#A3460C", "#FBDFCB"),
    ("#1D6E78", "#D4EEF0"),
    ("#7A2E58", "#F6DCEA"),
]
GRADE_COLORS = {
    "FigureQuality": ("#3C3489", "#E2DFFB"),
    "Good":          ("#0F6E56", "#CDEEE1"),
    "Ok":            ("#854F0B", "#FAEEDA"),
    "Bad":           ("#791F1F", "#F7C1C1"),
}


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 2 — Reading and writing the annotation files
# Everything that touches the output CSV and TXT. Edit here if you
# want a different output format.
# ══════════════════════════════════════════════════════════════════════════════
RESERVED_COLS = {"tomo_name", "quality", "comment"}


def csv_error(csv_path, problem, headers=None):
    """Explain why a CSV can't be used, show what one looks like, and stop.

    Same idea as the import checks at the top of the file: say what is
    wrong, say what was expected, and exit — rather than carrying on with
    nothing loaded and overwriting the file at the first save.

    Never returns.
    """
    print(f"Could not read the annotation CSV: {csv_path}")
    print(f"  Problem: {problem}")
    if headers is not None:
        print(f"  Columns found: {', '.join(headers) if headers else '(none)'}")
    print()
    print("  A TOMATO CSV has this header line, in this order:")
    print("      tomo_name,quality,<one column per flag>,comment")
    print()
    print("  for example:")
    print("      tomo_name,quality,Mito,ER,Golgi,comment")
    print("      TS_0001.mrc,Good,1,1,0,nice ribosomes")
    print("      TS_0002.mrc,Bad,0,0,0,")
    print()
    print("  tomo_name is the tomogram's filename, quality is one of")
    print(f"  {'/'.join(GRADES)} (or NA), and every other column is a flag")
    print("  holding a count: 0 for absent, 1 or more for present.")
    print()
    print("  This file is normally written by TOMATO itself — check that it")
    print("  is an annotation file and not some other CSV.")
    sys.exit(1)


def parse_count(value):
    """Read one CSV cell as a flag count.

    Deliberately forgiving: anything unrecognised becomes 0 rather than
    raising, so one malformed cell can't stop a whole CSV from loading and
    cost someone a session's annotation. Accepts a plain integer, and also
    the words True and False, as 1 and 0.
    """
    value = (value or "").strip()
    if value == "":
        return 0
    if value == "True":
        return 1
    if value == "False":
        return 0
    try:
        return max(0, int(float(value)))
    except ValueError:
        return 0


def _row_to_entry(row, flags):
    """Turn one CSV DictReader row into {flag: count, 'quality', 'comment'}.

    Called once per row by load_csv_with_flags.
    """
    entry = {flag: parse_count(row.get(flag)) for flag in flags}
    entry["quality"] = row.get("quality", "")
    entry["comment"] = row.get("comment", "")
    return entry


def load_csv_with_flags(csv_path):
    """Read a previously written TOMATO CSV.

    Returns (flags_list, data_dict), where flags_list is the ordered flag
    column names — everything except tomo_name, quality and comment.
    """
    if not csv_path.exists():
        csv_error(csv_path, "the file does not exist")
    data = {}
    try:
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            headers = list(reader.fieldnames or [])
            if "tomo_name" not in headers:
                # Without this check every column becomes a flag, so an
                # unrelated CSV would fill the GUI with buttons named after
                # somebody else's spreadsheet.
                csv_error(csv_path, "there is no 'tomo_name' column", headers)
            flags_from_csv = [h for h in headers if h not in RESERVED_COLS]
            for row in reader:
                name = row.get("tomo_name")
                if not name:
                    continue
                data[name] = _row_to_entry(row, flags_from_csv)
    except OSError as e:
        csv_error(csv_path, f"the file could not be opened ({e.strerror})")
    except UnicodeDecodeError:
        csv_error(csv_path, "the file is not text — a binary file named .csv?")
    except csv.Error as e:
        csv_error(csv_path, f"the file is not valid CSV ({e})")
    return flags_from_csv, data


def write_csv(csv_path, flags, all_data):
    """Write the full annotation table: one row per tomogram.

        tomo_name,quality,Mito,ER,Golgi,comment
        TS_0001.mrc,Good,1,1,0,nice ribosomes
        TS_0002.mrc,Bad,0,0,0,

    A simple flag is 0 or 1; a counter flag can be any whole number.
    Quality is written as "NA" if it was never set. The same data also
    goes to a compact .txt summary — see write_txt.
    """
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tomo_name", "quality"] + flags + ["comment"])
        for name in sorted(all_data.keys()):
            quality = all_data[name].get("quality", "") or "NA"
            row_vals = [str(all_data[name].get(flag, 0)) for flag in flags]
            comment = all_data[name].get("comment", "")
            writer.writerow([name, quality] + row_vals + [comment])


def write_txt(txt_path, flags, all_data):
    """Write the compact summary, one line per tomogram.

        TS_0001.mrc,Good,Mito,ER,Shell:3

    Only flags with a count above 0 appear — as a bare name when the count
    is 1, as name:N when it is higher. Easier to grep or paste into a
    figure legend than the full CSV; write_csv has the same data in full.
    """
    with open(txt_path, "w") as f:
        for name in sorted(all_data.keys()):
            entry = all_data[name]
            quality = entry.get("quality", "") or "NA"
            active = []
            for flag in flags:
                count = entry.get(flag, 0)
                if count <= 0:
                    continue
                active.append(flag if count == 1 else f"{flag}:{count}")
            parts = [name, quality] + active
            f.write(",".join(parts) + "\n")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 3 — Turning a tomogram into images
# Loading an .mrc, choosing its contrast, downsampling it and
# pre-rendering the slices the GUI displays.
# ══════════════════════════════════════════════════════════════════════════════
def _read_pixel_size(mrc):
    """Read the Å per voxel value out of an MRC header.

    Assumes isotropic voxels — the X value is used for all three axes.
    Falls back to 1.0 if the field is missing, zero or unreadable; the
    Å/px box in the readout corner lets the user correct a wrong header
    by hand.
    """
    try:
        val = float(mrc.voxel_size.x)
        if val > 0:
            return val
    except (AttributeError, TypeError, ValueError):
        pass
    return 1.0


def load_volume(path):
    """Load one tomogram from disk into a volume.

    Returns (volume, Å per voxel). permissive=True so a slightly
    non-standard header is a warning rather than a refusal to open.
    """
    with mrcfile.open(str(path), mode="r", permissive=True) as mrc:
        data = mrc.data.copy()
        pixel_size = _read_pixel_size(mrc)
    return data, pixel_size


def compute_contrast_window(data, n_sample_slices=5, std_factor=3.0):
    """Choose which range of intensities is shown as black-to-white.

    Tomogram intensities have no fixed scale, so "black" and "white" have
    to be picked per tomogram. Everything from mean - 3*std to mean + 3*std
    is spread across the greyscale; anything outside is clipped to pure
    black or pure white. A larger std_factor gives a flatter, greyer image,
    a smaller one more contrast.

    DISPLAY ONLY. The volume is never modified and none of this reaches the
    CSV — two people using different contrast annotate the same numbers.

    Measured on 5 evenly spaced Z-slices rather than the whole volume:
    enough for a stable window, at a small fraction of the cost.
    """
    n_slices = data.shape[0]
    if n_slices <= n_sample_slices:
        sample_indices = range(n_slices)
    else:
        sample_indices = np.linspace(0, n_slices - 1, n_sample_slices, dtype=int)
    samples = [np.asarray(data[i]).astype(np.float32) for i in sample_indices]
    stacked = np.concatenate([s.ravel() for s in samples])
    mean = float(np.mean(stacked))
    std  = float(np.std(stacked))
    if std == 0:
        std = 1.0
    return mean - std_factor * std, mean + std_factor * std


def bin_downsample(arr, target_size):
    """Shrink a slice by averaging square blocks of pixels together.

    Used when a tomogram is larger than the panel it has to fit in.
    Averaging first is what keeps noisy tomogram data readable: simply
    throwing pixels away (which is what plain interpolation does) makes
    noise alias into what looks like structure — actively misleading on a
    tomogram.

    The array is cropped to a whole number of blocks first, so up to
    bin_factor-1 rows/columns at the far edge are dropped.
    """
    h, w = arr.shape
    if h <= target_size or w <= target_size:
        return arr
    bin_factor = max(1, min(h // target_size, w // target_size))
    if bin_factor <= 1:
        return arr
    new_h = (h // bin_factor) * bin_factor
    new_w = (w // bin_factor) * bin_factor
    cropped = arr[:new_h, :new_w]
    reshaped = cropped.reshape(new_h // bin_factor, bin_factor,
                               new_w // bin_factor, bin_factor)
    return reshaped.mean(axis=(1, 3))


def array_to_image(arr2d, contrast_window, out_w, out_h):
    """Turn one 2D slice of raw intensities into a greyscale image.

    Three steps: average down if the slice is bigger than the panel
    (bin_downsample), map the contrast window onto black-to-white
    (compute_contrast_window), then resize to exactly (out_w, out_h).
    """
    arr2d = np.asarray(arr2d).astype(np.float32)
    arr2d = bin_downsample(arr2d, max(int(out_w), int(out_h)))
    lo, hi = contrast_window
    if hi <= lo:
        hi = lo + 1.0
    arr2d = np.clip((arr2d - lo) / (hi - lo), 0, 1)
    arr2d = (arr2d * 255).astype(np.uint8)
    img = Image.fromarray(arr2d, mode="L")
    return img.resize((max(1, int(out_w)), max(1, int(out_h))), Image.BILINEAR)


def prerender_z_stack(data, contrast_window, out_w, out_h, progress_callback=None):
    """Render every Z slice up front, so scrolling through Z is instant.

    Z is by far the most scrolled axis, so it earns the cost; XZ and YZ are
    rendered on demand instead (see Annotator._render_xz).

    Costs memory: one image per Z slice, held for as long as the tomogram
    is open. That is the main reason MAX_XY exists, and why prefetching
    steps itself down under memory pressure (see PART 8).
    """
    n_slices = data.shape[0]
    images = []
    for z in range(n_slices):
        img = array_to_image(data[z], contrast_window, out_w, out_h)
        images.append(ImageTk.PhotoImage(img))
        if progress_callback is not None and (z % 5 == 0 or z == n_slices - 1):
            progress_callback(z + 1, n_slices)
    return images


def panel_dims(nx, ny, nz, scale):
    """Canvas pixel sizes for the three panels at a given voxel->pixel scale.

    One shared scale for all three keeps the real aspect ratio, so XZ and
    YZ come out as true thin strips rather than stretched squares. Used by
    both _apply_panel_geometry and _prefetch_worker, so the current and
    prefetched tomograms can never be sized differently.
    """
    xy_w, xy_h = max(1, round(nx * scale)), max(1, round(ny * scale))
    xz_w, xz_h = xy_w, max(1, round(nz * scale))
    yz_w, yz_h = xz_h, xy_h
    return xy_w, xy_h, xz_w, xz_h, yz_w, yz_h


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 4 — The annotator GUI
# The whole window: layout, controls, rendering and navigation.
# Long, but split into 8 PARTs — see the map below.
# ══════════════════════════════════════════════════════════════════════════════
# Filenames are split into words and the separators between them, so that
# trimming happens on word boundaries. Without this, TS_0001 and TS_0002
# would be found to share the characters "TS_000" and get displayed as "1"
# and "2".
_NAME_TOKENS = re.compile(r"[^_.\-]+|[_.\-]")


def shorten_display_names(filenames):
    """Shorten filenames for the GUI by dropping the parts they all share.

    Tomogram names are long and mostly identical — the informative part is
    usually a few characters in the middle. Given

        TS_20240115_0001_rec.mrc
        TS_20240115_0002_rec.mrc

    the shared start (TS_20240115_) and end (_rec.mrc) carry no information
    within this run, so the list shows 0001 and 0002.

    DISPLAY ONLY. The CSV always records the full filename. Returns
    {filename: label}; a filename with no entry is shown in full.

    Backs off to full names if trimming would make two tomograms look the
    same — better long than ambiguous.
    """
    names = list(filenames)
    if len(names) < 2:
        return {n: n for n in names}
    tokens = [_NAME_TOKENS.findall(n) for n in names]
    # Count the leading, then the trailing, tokens every name agrees on.
    # The "len(t) > lead + tail + 1" guard keeps at least one token, so a
    # set of identically-named files can never trim down to nothing.
    lead = 0
    while (all(len(t) > lead + 1 for t in tokens)
           and len({t[lead] for t in tokens}) == 1):
        lead += 1
    tail = 0
    while (all(len(t) > lead + tail + 1 for t in tokens)
           and len({t[-1 - tail] for t in tokens}) == 1):
        tail += 1
    short = {}
    for name, t in zip(names, tokens):
        label = "".join(t[lead:len(t) - tail]).strip("_-.")
        short[name] = label or name
    if len(set(short.values())) != len(names):
        return {n: n for n in names}
    return short


class Annotator:
    """The annotation window: everything TOMATO does interactively.

    Everything below lives in one class because it all shares state (the
    loaded volume, the crosshair position, the current annotation...).
    That makes it long, so here's a map of the 8 parts it's organised
    into, in the order they appear in the file:

      PART 1 — Setup & Lifecycle
               __init__, small path/save helpers, and _build_gui (which
               constructs every widget in the window).
      PART 2 — Window Sizing & Layout Fitting
               Works out how big the window and each column should be,
               both at startup and when the user resizes it by hand.
      PART 3 — Flag & Quality Controls
               Builds and updates the organelle flag buttons and the
               quality (Good/Ok/Bad/FigureQuality) segmented control.
      PART 4 — Crosshair Position & Panel Rendering
               Moves the crosshair and redraws the XY/XZ/YZ panels to
               match it.
      PART 5 — Pixel Size (Å/voxel)
               Tracks and validates the Å/px value shown in the readout
               corner, used only by the ruler tool.
      PART 6 — Mouse & Wheel Input, Zoom / Pan
               Turns raw mouse/wheel events into crosshair moves, zoom,
               pan, or a ruler measurement.
      PART 7 — Tomogram Navigation & Loading
               Saving the current tomogram's annotation, moving to the
               next/previous one, and loading a volume from disk.
      PART 8 — Background Prefetch
               Loads upcoming tomograms on worker threads while the
               current one is being annotated, so "Next" is instant.

    Search for "PART " to jump between them.
    """

    # ──────────────────────────────────────────────────────────────────────
    # PART 1 — Setup & Lifecycle
    # __init__, small path/save helpers, and _build_gui (which
    # constructs every widget in the window).
    # ──────────────────────────────────────────────────────────────────────

    def __init__(self, tomo_paths, flags, output_path, txt_path, initial_data,
                 initial_flag_mode=None, pixel_size=None):
        """Set up all the state one annotation session needs, then open it.

        Two kinds of state live here, and mixing them up is the easiest
        mistake to make in this class:

          per-tomogram   current_state, current_quality, volume, crosshair.
                         Reset by _load_current on every navigation.
          per-session    flags, flag_mode, prefetch_ahead,
                         _pixel_size_override. Set once and deliberately
                         survive navigation — promoting a flag to a counter
                         keeps it a counter for every tomogram after it.

        Ends by calling _build_gui, which does not return until the window
        is closed.
        """
        self.tomo_paths  = tomo_paths
        self.flags       = list(flags)
        self.output_path = output_path
        self.txt_path    = txt_path
        # Worked out once from the whole selection, since what is worth
        # showing depends on what the other filenames look like.
        self.display_names = shorten_display_names(p.name for p in tomo_paths)
        self.idx         = 0
        self.volume           = None
        self.xy_cache         = []     # prerendered PhotoImages, one per Z slice
        self.contrast_window  = (0.0, 1.0)
        # Volume dimensions and crosshair position, all in voxel coords
        self.nz = self.ny = self.nx = 1
        self.x = self.y = self.z    = 0
        self.scale = 1.0
        # The panel pixel sizes (xy_w ... yz_h) are set in _build_gui, so
        # they need no placeholder here.
        self.photo_xy = self.photo_xz = self.photo_yz = None

        # Synced zoom/pan state — see PART 6.
        # self.scale above stays the layout value that fits the whole volume
        # in the canvas; zoom is a separate render scale that starts equal
        # to it, and x0/y0/z0 are the voxel coordinate at each canvas's
        # top-left corner. All reset to fit on load and on resize.
        self.fit_zoom  = 1.0
        self.zoom      = 1.0
        self.MAX_ZOOM  = 40.0
        self.x0 = self.y0 = self.z0 = 0.0
        self._drag_start      = None   # (panel, start_x, start_y) for left-button click-vs-ruler-drag
        self._drag_is_measuring = False
        self._pan_start       = None   # (panel, start_x, start_y, x0, y0, z0) for right/middle-drag panning

        # Pixel size (Å/voxel), shown and editable in the readout corner.
        # _pixel_size_override stays None until the user types a value (or
        # passes --pixel_size), after which it wins over the MRC header for
        # every tomogram from then on: a wrong header is usually wrong for
        # the whole batch, not just one file.
        self._pixel_size_override = pixel_size if pixel_size else None
        self.pixel_size = pixel_size if pixel_size else 1.0
        self.data            = initial_data
        self.current_state   = {flag: 0 for flag in self.flags}
        self.current_quality = ""
        self._list_updating  = False   # guard against recursive listbox events
        # Which flags are plain toggles and which are counters. This is a
        # per-session UI choice, NOT per-tomogram state: _load_current
        # resets the annotation fields on every navigation but deliberately
        # never touches flag_mode, so a promoted flag stays a counter for
        # every tomogram afterwards. Seeded from initial_flag_mode (see
        # _main_annotate); anything unlisted defaults to "simple".
        self.flag_mode = dict(initial_flag_mode or {})
        for flag in self.flags:
            self.flag_mode.setdefault(flag, "simple")
        # Background prefetch of the next tomogram(s) while the current one
        # is being annotated, so "Next" is instant. See
        # PREFETCH_AHEAD_DEFAULT for how prefetch_ahead behaves at runtime.
        self.prefetch_ahead    = PREFETCH_AHEAD_DEFAULT
        self._prefetch_lock    = threading.Lock()
        self._prefetch_cache   = {}    # {idx: {...}} — up to prefetch_ahead entries
        self._prefetch_targets = set() # indices with a fetch thread currently in flight
        self._build_gui()

    # ── Small helpers used throughout the class ──
    def _tomo_name(self, path):
        """The label to show for a tomogram, from shorten_display_names.

        DISPLAY ONLY. What gets saved is always the full filename — see
        _commit_current.
        """
        return self.display_names.get(path.name, path.name)

    def _save_outputs(self):
        write_csv(self.output_path, self.flags, self.data)
        write_txt(self.txt_path,    self.flags, self.data)

    # ── GUI construction ──
    def _build_gui(self):
        """Build the whole window and hand control to Tk.

        Does not return until the user quits — root.mainloop() at the end
        blocks. Everything after that line runs on Tk's events, not here.
        """
        # This one method builds the entire window, in the order Tk needs
        # things to exist. It is long, so it is split into 8 numbered
        # STEPs — search for "STEP " to jump between them:
        #
        #   STEP 1  The window itself
        #   STEP 2  Starting sizes (placeholders, refined later)
        #   STEP 3  Bottom navigation row  (Previous / Quit / Z-slider / Next)
        #   STEP 4  The main three-column frame
        #   STEP 5  Column 0 — image viewer (the XY/XZ/YZ quad view)
        #   STEP 6  Column 1 — quality, flags, add-flag box, comment box
        #   STEP 7  Column 2 — tomogram list
        #   STEP 8  Final fit, first tomogram, and hand over to Tk
        #
        # Nothing appears on screen until STEP 8 calls root.mainloop().

        # ── STEP 1 — The window itself ──
        self.root = tk.Tk()
        self.root.title("🍅 TOMATO - TOMogram Annotation TOol 🍅")
        self.root.resizable(True, True)
        self.root.protocol("WM_DELETE_WINDOW", self._on_quit)
        # On some platforms a wheel event's Control bit arrives a tick late,
        # so Ctrl+scroll sometimes just steps Z instead of zooming. Track
        # Control from its own key events too, as a backstop. bind_all, so
        # it doesn't matter which widget has focus.
        self._ctrl_held = False
        self.root.bind_all("<KeyPress-Control_L>",   lambda e: setattr(self, "_ctrl_held", True))
        self.root.bind_all("<KeyPress-Control_R>",   lambda e: setattr(self, "_ctrl_held", True))
        self.root.bind_all("<KeyRelease-Control_L>", lambda e: setattr(self, "_ctrl_held", False))
        self.root.bind_all("<KeyRelease-Control_R>", lambda e: setattr(self, "_ctrl_held", False))
        # If focus leaves the app while Ctrl is held (alt-tabbing), the
        # KeyRelease above never fires and _ctrl_held would stick on
        # forever, so clear it whenever the window loses focus.
        self.root.bind_all("<FocusOut>", lambda e: setattr(self, "_ctrl_held", False), add="+")
        self.screen_w = self.root.winfo_screenwidth()
        self.screen_h = self.root.winfo_screenheight()

        # ── STEP 2 — Starting sizes ──
        # Placeholder numbers only. Real ones are measured in STEP 8,
        # once the widgets below actually exist and can be asked how
        # big they are.
        # The quad view gets whatever the annotation and list columns and
        # the nav row leave, so nothing is ever hidden off-screen. These
        # are just placeholders so the canvases can be created — real
        # budgets come from _measure_budgets, later in this method.
        self.min_xy = MIN_XY
        self.max_xy = MAX_XY
        self.quad_w_budget = self.max_xy
        self.quad_h_budget = self.max_xy
        self.list_height   = 24
        # Height of the scrollable flag list; refined in
        # _update_flag_area_budget once the widgets can be measured.
        self.flag_area_h   = max(200, self.screen_h - 420)
        # Usable space for the columns; set by _measure_budgets once the
        # widgets exist. None means "not measured yet".
        self.avail_h       = None
        self.avail_w       = None
        # Measuring before a volume is loaded is slightly optimistic (fonts,
        # borders and label text all settle later), so _calibrate_fit folds
        # any error into these after the first real load — otherwise the
        # window can end up taller than the screen with the slider clipped
        # off the bottom.
        self.quad_h_correction = 0
        self.quad_w_correction = 0
        self._calibrated       = False
        self.xy_w = self.xy_h = self.min_xy
        self.xz_w = self.xz_h = self.min_xy
        self.yz_w = self.yz_h = self.min_xy
        BORDER_RED  = "#8B1E1E"
        INTERIOR_RED = "#F4E2DE"
        self.root.configure(bg=BORDER_RED)

        # ── STEP 3 — Bottom navigation row ──
        # Previous / Quit / Z-slider / Next, across the bottom of the
        # window. Packed first and pinned to the bottom, so it stays
        # visible whatever happens to the columns above it.
        nav_frame = tk.Frame(self.root, bg="#F4E2DE")
        nav_frame.pack(side="bottom", fill="x", padx=6, pady=6)
        self.nav_frame_bottom = nav_frame
        self.prev_btn = tk.Button(nav_frame, text="← Previous",
                  font=("Helvetica", 12), width=12,
                  fg="#2C2C2A", bg="#F1EFE8",
                  relief="flat", bd=0, cursor="hand2",
                  command=self._go_previous)
        self.prev_btn.pack(side="left")
        tk.Button(nav_frame, text="Quit",
                  font=("Helvetica", 11), fg="#888780", bg="#F4E2DE",
                  relief="flat", bd=0, cursor="hand2",
                  command=self._on_quit).pack(side="left", padx=20)
        self.next_btn = tk.Button(nav_frame, text="Next →",
                  font=("Helvetica", 12, "bold"), width=12,
                  fg="#FAFAF8", bg="#2C2C2A",
                  relief="flat", bd=0, cursor="hand2",
                  command=self._go_next)
        self.next_btn.pack(side="right")
        # Z slider fills the gap between Quit and Next.
        slider_box = tk.Frame(nav_frame, bg="#F4E2DE")
        slider_box.pack(side="left", fill="x", expand=True, padx=16)
        self.slice_var = tk.StringVar()
        tk.Label(slider_box, textvariable=self.slice_var,
                 font=("Helvetica", 10), bg="#F4E2DE", fg="#888780",
                 width=16, anchor="e").pack(side="left", padx=(0, 8))
        self.slider = tk.Scale(slider_box, from_=0, to=1, orient="horizontal",
                               showvalue=False,
                               bg="#F4E2DE", troughcolor="#F1EFE8",
                               highlightthickness=0, command=self._on_slider)
        self.slider.pack(side="left", fill="x", expand=True)

        # ── STEP 4 — The main three-column frame ──
        # Everything above the navigation row sits in `main`, as three
        # side-by-side columns: 0 = image viewer, 1 = annotation panel,
        # 2 = tomogram list. Built in STEPs 5, 6 and 7 respectively.
        main = tk.Frame(self.root, bg=INTERIOR_RED)
        # fill+expand so `main` grows when the window is dragged bigger than
        # its natural size. At exactly root.minsize() there is no slack to
        # distribute, so this has no effect there.
        main.pack(side="top", padx=6, pady=6, fill="both", expand=True)
        self.main_frame = main
        # Only column 0 (the image-viewer column) gets weight, so any
        # slack `main` receives goes entirely to the quad view — the
        # annotation and list columns (1, 2) stay at their fitted size.
        main.columnconfigure(0, weight=1)
        main.rowconfigure(0, weight=1)

        # ── STEP 5 — Column 0: image viewer ──
        # The name/counter header, the 2x2 quad view (XY, XZ, YZ plus the
        # position readout), and the mouse and wheel bindings that make
        # those three panels navigable.
        left = tk.Frame(main, bg="#F4E2DE")
        left.grid(row=0, column=0, sticky="nsew")
        self.left_col = left
        # Counter and name share one line: two stacked labels cost ~60px of
        # height that the quad view can use instead.
        header = tk.Frame(left, bg="#F4E2DE")
        header.pack(fill="x", pady=(0, 4))
        # A StringVar is a small Tk-managed box holding a string. A widget
        # built with textvariable=some_var displays whatever's inside that
        # box, and updates on screen automatically whenever .set() is
        # called on the var — no need to reach back into the widget itself
        # later. Used throughout this method for every label/entry whose
        # text changes after creation (progress_var, name_var, pos_var...).
        self.progress_var = tk.StringVar()
        tk.Label(header, textvariable=self.progress_var,
                 font=("Helvetica", 11), bg="#F4E2DE",
                 fg="#888780", anchor="w").pack(side="left", padx=(0, 10))
        self.name_var = tk.StringVar()
        tk.Label(header, textvariable=self.name_var,
                 font=("Helvetica", 14, "bold"),
                 bg="#F4E2DE", fg="#2C2C2A", anchor="w").pack(side="left")
        # 2x2 quad view:
        #   XZ (top-left)      | position readout (top-right)
        #   XY (bottom-left, main) | YZ (bottom-right)
        quad = tk.Frame(left, bg="#F4E2DE")
        # fill+expand: quad claims any slack `left` has (the header stays
        # fill="x" only, so it never competes for it). This is what lets
        # _apply_quad_resize read the quad's size directly — see the note
        # above _on_quad_configure.
        quad.pack(fill="both", expand=True)
        self.quad_frame = quad
        self.canvas_xz = tk.Canvas(quad, width=self.xz_w, height=self.xz_h,
                                   bg="black", highlightthickness=1,
                                   highlightbackground="#8B1E1E")
        self.canvas_xz.grid(row=0, column=0, sticky="sw")
        self.nav_frame = tk.Frame(quad, bg="#2C2C2A",
                                  width=self.yz_w, height=self.xz_h)
        self.nav_frame.grid(row=0, column=1, sticky="nw", padx=(2, 0))
        self.nav_frame.grid_propagate(False)
        self.pos_var = tk.StringVar()
        tk.Label(self.nav_frame, textvariable=self.pos_var,
                 font=("Helvetica", 9), bg="#2C2C2A", fg="#E8E6E0",
                 justify="left", anchor="nw").pack(fill="x", padx=4, pady=(4, 6))
        # Å/px: the pixel size currently in use, from the MRC header or
        # from the user (see _pixel_size_override in __init__).
        px_row = tk.Frame(self.nav_frame, bg="#2C2C2A")
        px_row.pack(fill="x", padx=4, pady=(0, 6))
        tk.Label(px_row, text="Å/px:", font=("Helvetica", 9),
                 bg="#2C2C2A", fg="#E8E6E0").pack(side="left")
        self.pixel_size_var = tk.StringVar(value=f"{self.pixel_size:.4f}")
        self.pixel_size_entry = tk.Entry(px_row, textvariable=self.pixel_size_var,
                 font=("Helvetica", 9), width=8, bg="#F1EFE8", fg="#2C2C2A",
                 relief="flat", bd=0, justify="right")
        self.pixel_size_entry.pack(side="left", padx=(4, 0))
        self.pixel_size_entry.bind("<Return>", self._on_pixel_size_entry)
        self.pixel_size_entry.bind("<FocusOut>", self._on_pixel_size_entry)
        # Ruler readout: a drag on any panel shows its length here.
        # Display-only — nothing is saved, and the next drag replaces it.
        self.measured_var = tk.StringVar(value="Measured size: —")
        tk.Label(self.nav_frame, textvariable=self.measured_var,
                 font=("Helvetica", 9), bg="#2C2C2A", fg="#F5E642",
                 justify="left", anchor="nw", wraplength=max(80, self.yz_w - 8)
                 ).pack(fill="x", padx=4, pady=(0, 4))
        self.canvas_xy = tk.Canvas(quad, width=self.xy_w, height=self.xy_h,
                                   bg="black", highlightthickness=1,
                                   highlightbackground="#0F6E56")
        self.canvas_xy.grid(row=1, column=0, sticky="nw", pady=(2, 0))
        self.canvas_yz = tk.Canvas(quad, width=self.yz_w, height=self.yz_h,
                                   bg="black", highlightthickness=1,
                                   highlightbackground="#1B4F8A")
        self.canvas_yz.grid(row=1, column=1, sticky="nw", padx=(2, 0), pady=(2, 0))
        # Mouse map, the same on all three panels: left click moves the
        # crosshair, a left drag past DRAG_THRESHOLD draws the ruler
        # instead, right or middle drag pans, the wheel steps the axis the
        # panel doesn't show, and Ctrl+wheel zooms.
        panel_canvases = {"xy": self.canvas_xy, "xz": self.canvas_xz, "yz": self.canvas_yz}
        step_fns       = {"xy": self._step_z,   "xz": self._step_y,   "yz": self._step_x}
        for panel, canvas in panel_canvases.items():
            canvas.bind("<ButtonPress-1>",   lambda e, p=panel: self._on_panel_press(p, e))
            canvas.bind("<B1-Motion>",       lambda e, p=panel: self._on_panel_motion(p, e))
            canvas.bind("<ButtonRelease-1>", lambda e, p=panel: self._on_panel_release(p, e))
            wheel_handler = (lambda e, p=panel, fn=step_fns[panel]: self._panel_wheel(e, p, fn))
            canvas.bind("<MouseWheel>", wheel_handler)
            canvas.bind("<Button-4>",   wheel_handler)
            canvas.bind("<Button-5>",   wheel_handler)
            for b in (2, 3):
                canvas.bind(f"<ButtonPress-{b}>",   lambda e, p=panel: self._on_pan_press(p, e))
                canvas.bind(f"<B{b}-Motion>",       self._on_pan_drag)
                canvas.bind(f"<ButtonRelease-{b}>", self._on_pan_release)
        # The Z slider and its readout live in the bottom navigation row,
        # not in this column — that hands their height to the quad view.

        # ── STEP 6 — Column 1: quality + annotation panel ──
        # Stacked top to bottom: the quality grades, the scrollable flag
        # list, the "Add flag" box, the comment box, and the status line.
        right = tk.Frame(main, bg="#F4E2DE")
        right.grid(row=0, column=1, sticky="n", padx=(24, 0))
        self.right_col = right
        tk.Label(right, text="Quality",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(0, 6))
        # Segmented control: all four grade buttons in one bordered frame,
        # separated by hairlines.
        GRADE_LABELS = {
            "FigureQuality": "Figure\nQuality",
            "Good":          "Good",
            "Ok":            "Ok",
            "Bad":           "Bad",
        }
        SEGMENT_BORDER = "#2C2C2A"
        segment_outer = tk.Frame(right, bg=SEGMENT_BORDER)
        segment_outer.pack(pady=(0, 16))
        self.segment_outer = segment_outer
        segment_inner = tk.Frame(segment_outer, bg=SEGMENT_BORDER)
        segment_inner.pack(padx=1, pady=1, fill="x")
        self.grade_buttons = {}
        for i, grade in enumerate(GRADES):
            fg, bg_on = GRADE_COLORS[grade]
            # Equal weight everywhere, so the buttons stretch to fill
            # segment_outer's width (set by _match_quality_width_to_flags)
            # instead of each button sizing itself from its own text.
            segment_inner.columnconfigure(i, weight=1)
            b = tk.Button(segment_inner, text=GRADE_LABELS[grade],
                          font=("Helvetica", 12, "bold"),
                          fg=fg, bg="#F1EFE8",
                          activebackground=bg_on, activeforeground=fg,
                          relief="flat", bd=0, cursor="hand2",
                          height=2,
                          command=lambda g=grade: self._set_quality(g))
            b.grid(row=0, column=i, sticky="nsew",
                   padx=(0 if i == 0 else 1, 0))
            self.grade_buttons[grade] = (b, fg, bg_on)
        tk.Label(right, text="Annotation",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(0, 6))
        tk.Label(right, text="Click a flag to toggle it. Click its − / + "
                             "to turn it into a counter, or × simple on a "
                             "counter to turn it back — both stick for the "
                             "rest of the session.",
                 font=("Helvetica", 9), bg="#F4E2DE", fg="#888780",
                 anchor="w", justify="left",
                 wraplength=220).pack(fill="x", pady=(0, 6))
        # Flags live in a fixed-height scrollable list: flag_outer never
        # changes size once laid out, so adding a flag can never resize the
        # annotation column or the window. A scrollbar appears once there
        # are more flags than fit in flag_area_h.
        flag_outer = tk.Frame(right, bg="#F4E2DE", highlightthickness=1,
                              highlightbackground="#C9C6BC")
        flag_outer.pack(fill="x")
        self.flag_outer = flag_outer
        self.flag_canvas = tk.Canvas(flag_outer, bg="#F4E2DE",
                                     highlightthickness=0,
                                     height=self.flag_area_h)
        flag_scroll = tk.Scrollbar(flag_outer, orient="vertical",
                                   command=self.flag_canvas.yview)
        self.flag_canvas.configure(yscrollcommand=flag_scroll.set)
        flag_scroll.pack(side="right", fill="y")
        self.flag_canvas.pack(side="left", fill="both", expand=True)
        # Tkinter has no built-in "scrollable frame" widget. The standard
        # workaround: put the real content (buttons_frame, holding all the
        # flag rows) inside a Canvas as an embedded window, then scroll the
        # Canvas itself. flag_canvas.bbox("all") below tells the scrollbar
        # how tall buttons_frame currently is, growing as flags are added.
        self.buttons_frame = tk.Frame(self.flag_canvas, bg="#F4E2DE")
        self._flag_window = self.flag_canvas.create_window(
            (0, 0), window=self.buttons_frame, anchor="nw")
        self.buttons_frame.bind(
            "<Configure>",
            lambda e: self.flag_canvas.configure(scrollregion=self.flag_canvas.bbox("all")))
        self.flag_canvas.bind(
            "<Configure>",
            lambda e: self.flag_canvas.itemconfigure(self._flag_window, width=e.width))
        self._bind_flag_wheel(self.flag_canvas)
        self.flag_buttons = {}
        for flag in self.flags:
            self._add_flag_button(flag)
        # The flag block's real width is known now, so size the quality
        # control to match it. Only needed once — see the docstring.
        self._match_quality_width_to_flags()
        tk.Label(right, text="Add flag",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(14, 4))
        add_frame = tk.Frame(right, bg="#F4E2DE")
        add_frame.pack(fill="x")
        self.add_entry = tk.Entry(add_frame, font=("Helvetica", 11),
                                  bg="#F1EFE8", fg="#2C2C2A",
                                  relief="flat", bd=6, width=14,
                                  insertbackground="#2C2C2A")
        self.add_entry.pack(side="left")
        self.add_entry.bind("<Return>", lambda e: self._add_flag_from_entry())
        tk.Button(add_frame, text="+", font=("Helvetica", 12, "bold"),
                  fg="#FAFAF8", bg="#2C2C2A",
                  relief="flat", bd=0, cursor="hand2",
                  padx=10, pady=4,
                  command=self._add_flag_from_entry).pack(side="left", padx=(6, 0))
        tk.Label(right, text="Comment",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(16, 4))
        self.comment_box = tk.Text(right, width=24, height=6,
                                   font=("Helvetica", 11),
                                   bg="#F1EFE8", fg="#2C2C2A",
                                   relief="flat", bd=6, wrap="word",
                                   insertbackground="#2C2C2A")
        self.comment_box.pack(fill="x")
        self.status_var = tk.StringVar()
        tk.Label(right, textvariable=self.status_var,
                 font=("Helvetica", 10), bg="#F4E2DE",
                 fg="#888780", anchor="w", justify="left",
                 wraplength=220).pack(fill="x", pady=(10, 4))

        # ── STEP 7 — Column 2: tomogram list ──
        # Every tomogram found in --input, as a clickable list for
        # jumping straight to one instead of paging through with "Next".
        list_col = tk.Frame(main, bg="#F4E2DE")
        list_col.grid(row=0, column=2, sticky="ns", padx=(24, 0))
        self.list_col = list_col
        tk.Label(list_col, text="Tomograms",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(0, 6))
        list_frame = tk.Frame(list_col, bg="#F4E2DE")
        list_frame.pack(fill="both", expand=True)
        scrollbar = tk.Scrollbar(list_frame, orient="vertical")
        scrollbar.pack(side="right", fill="y")
        self.tomo_list = tk.Listbox(
            list_frame,
            yscrollcommand=scrollbar.set,
            font=("Helvetica", 11),
            bg="#F1EFE8", fg="#2C2C2A",
            selectbackground="#2C2C2A", selectforeground="#FAFAF8",
            relief="flat", bd=0,
            highlightthickness=0,
            activestyle="none",
            width=40,
            height=self.list_height,
        )
        self.tomo_list.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=self.tomo_list.yview)
        for path in self.tomo_paths:
            self.tomo_list.insert("end", self._tomo_name(path))
        self.tomo_list.bind("<<ListboxSelect>>", self._on_list_select)

        # ── STEP 8 — Final fit, first tomogram, hand over to Tk ──
        # Every widget now exists and can be measured, so the real sizes
        # from STEP 2 get worked out here. root.mainloop() on the last
        # line shows the window and waits for the user — it does not
        # return until the window is closed.
        self._measure_budgets()
        self._load_current()
        # The columns were fitted before any volume existed, i.e. to the
        # screen rather than to the quad view. Now that its real height is
        # known, re-fit them to that instead.
        self._tighten_columns_to_quad_height()
        # Centre the window on screen now that its actual size is known.
        self.root.update_idletasks()
        w = self.root.winfo_reqwidth()
        h = self.root.winfo_reqheight()
        x = max(0, (self.screen_w - w) // 2)
        y = max(0, (self.screen_h - h) // 2)
        self.root.geometry(f"+{x}+{y}")
        # Floor: don't let the user (or a window manager "restore") drag the
        # window below the size just fitted, where the quad view would
        # shrink past MIN_XY and the tomogram would barely be visible.
        self.root.minsize(w, h)
        self._resize_job = None
        # Baseline for the change detection in _on_quad_configure, read from
        # the widget itself rather than assumed.
        self._last_quad_size = (self.quad_frame.winfo_width(),
                                self.quad_frame.winfo_height())
        # Bound last, once layout has settled: <Configure> also fires
        # several times while the window is being built and positioned, and
        # those early events don't reflect the final fit.
        self.quad_frame.bind("<Configure>", self._on_quad_configure)
        self.root.mainloop()


    # ──────────────────────────────────────────────────────────────────────
    # PART 2 — Window Sizing & Layout Fitting
    # Works out how big the window and each column should be, both at
    # startup and when the user resizes it by hand.
    # ──────────────────────────────────────────────────────────────────────
    # ── Reacting to manual resize / maximize ──
    #
    # The invariant, set up in _build_gui: quad_frame — not the root
    # window — absorbs whatever space the window gains or loses, so
    # quad_frame.winfo_width()/height() IS the quad-view budget, with no
    # "root size minus everything else" subtraction needed. It holds
    # because main expands and only its column 0 has weight, so slack
    # flows through left into quad while the annotation and list columns
    # stay fixed. root.minsize() keeps the window at or above the fitted
    # floor, so this only ever sees more room than that, never less.
    def _on_quad_configure(self, event):
        """Debounce quad_frame resizes: record the size, restart the timer.

        Fires for genuine window resizes and, harmlessly, for our own canvas
        resizes when a tomogram loads. The real work happens once, in
        _apply_quad_resize, 150ms after resizing settles.

        Tk fires a <Configure> event whenever a widget's size or position
        changes — including many times in a row during a manual window
        drag. Doing the expensive rescale work on every single one of
        those would be wasteful and janky, so this "debounces" it: each
        event just records the new size and (re)starts a 150ms timer via
        root.after(); root.after_cancel() below throws away the previous
        timer if another event arrives before it fires, so only the last
        one in a burst actually runs _apply_quad_resize.
        """
        size = (event.width, event.height)
        if size == self._last_quad_size:
            return
        self._last_quad_size = size
        if self._resize_job is not None:
            self.root.after_cancel(self._resize_job)
        self._resize_job = self.root.after(150, self._apply_quad_resize)

    def _apply_quad_resize(self):
        """Rescale the quad view to fill however much room it has been given.

        Reads quad_frame's own size directly — see the invariant above.
        """
        self._resize_job = None
        self.root.update_idletasks()
        self.quad_w_budget = max(self.min_xy, self.quad_frame.winfo_width())
        self.quad_h_budget = max(self.min_xy, self.quad_frame.winfo_height())
        # Raise the max_xy ceiling to match a deliberately grown window, or
        # growth would stop dead at MAX_XY however large the window gets.
        # The budget is always a generous enough ceiling: quad_w_budget
        # bounds (nx+nz)*scale, which is more than the nx*scale (or
        # ny*scale) that max_xy limits.
        self.max_xy = max(self.max_xy, self.quad_w_budget, self.quad_h_budget)
        # This can resize the canvases and so retrigger <Configure> here.
        # Harmless: the size will match the one just recorded, and
        # _rescale_current_view's own "scale barely changed" guard absorbs
        # any one-pixel rounding difference.
        self._rescale_current_view()

    # ── Fit-to-screen sizing ──
    def _usable_area(self):
        """Ask the window manager how large a window is allowed to be.

        Its own answer already excludes the title bar and any taskbar or
        dock, so we don't have to guess at those. With no window manager at
        all (a bare VNC session) Tk falls back to the screen minus a small
        margin, so this is safe everywhere.
        """
        try:
            max_w, max_h = self.root.wm_maxsize()
        except tk.TclError:
            max_w, max_h = 0, 0
        # Guard against a WM reporting something unusable or absurd.
        if not 200 < max_w <= self.screen_w * 2:
            max_w = self.screen_w - 60
        if not 200 < max_h <= self.screen_h * 2:
            max_h = self.screen_h - 80
        return min(max_w, self.screen_w), min(max_h, self.screen_h)

    def _measure_budgets(self):
        """Work out how much room each part of the window may occupy.

        The three columns sit side by side, so the window's height is set
        by the tallest — shrinking the quad view alone doesn't help if the
        annotation column is the long one. Each column is fitted to the
        available height in turn, then the quad view gets what is left.
        Runs once at startup.
        """
        self.root.update_idletasks()
        usable_w, usable_h = self._usable_area()
        nav_h = self.nav_frame_bottom.winfo_reqheight()
        # WINDOW_PAD covers the root frame's own padding (the tomato
        # "skin" border) above and below the columns.
        self.avail_h = avail_h = usable_h - nav_h - WINDOW_PAD
        self.avail_w = avail_w = usable_w - WINDOW_PAD
        # 1) Annotation column: give the flag list whatever height the rest
        #    of the column leaves, then shrink the comment box and tighten
        #    the padding if it still doesn't fit.
        self._update_flag_area_budget()
        self._fit_annotation_column(avail_h)
        # 2) Tomogram list: pick a row count that fills but doesn't exceed,
        #    and narrow it if the three columns don't fit side by side.
        self._fit_tomo_list(avail_h)
        self._fit_tomo_list_width(avail_w)
        # 3) Quad view gets whatever is left.
        self._recompute_quad_budgets()

    def _recompute_quad_budgets(self):
        """Give the quad view whatever space the other columns don't need.

        Re-run whenever those columns change size — for example when a flag
        was added and the buttons wrapped into another column.
        """
        if self.avail_h is None:
            return
        self.root.update_idletasks()
        left_overhead_h = max(0, self.left_col.winfo_reqheight()
                                 - self.quad_frame.winfo_reqheight())
        other_cols_w = max(0, self.main_frame.winfo_reqwidth()
                              - self.quad_frame.winfo_reqwidth())
        self.quad_w_budget = max(self.min_xy,
                                 self.avail_w - other_cols_w - self.quad_w_correction)
        self.quad_h_budget = max(self.min_xy,
                                 self.avail_h - left_overhead_h - self.quad_h_correction)

    def _calibrate_fit(self):
        """Check the laid-out window really fits, and correct it if not.

        Any overflow is folded into the quad budgets, so every later load
        and prefetch uses the corrected numbers.

        Runs once, after the first volume is laid out but before its slices
        are rendered — resizing canvases is cheap, re-rendering a few hundred
        slices is not.
        """
        self.root.update_idletasks()
        usable_w, usable_h = self._usable_area()
        excess_h = self.root.winfo_reqheight() - usable_h
        excess_w = self.root.winfo_reqwidth() - usable_w
        changed = False
        if excess_h > 0:
            self.quad_h_correction += excess_h
            changed = True
        if excess_w > 0:
            self.quad_w_correction += excess_w
            changed = True
        if changed:
            self._recompute_quad_budgets()
        return changed

    def _apply_panel_geometry(self):
        """Push the current scale out to the three canvases and slider."""
        (self.xy_w, self.xy_h, self.xz_w, self.xz_h,
         self.yz_w, self.yz_h) = panel_dims(self.nx, self.ny, self.nz, self.scale)
        for canvas, w, h in ((self.canvas_xy, self.xy_w, self.xy_h),
                             (self.canvas_xz, self.xz_w, self.xz_h),
                             (self.canvas_yz, self.yz_w, self.yz_h)):
            canvas.config(width=w, height=h)
        self.nav_frame.config(width=self.yz_w, height=self.xz_h)
        # The slider lives in the nav row and stretches to fill it, so only
        # its range depends on the volume.
        self.slider.config(from_=0, to=max(1, self.nz - 1))

    def _rescale_current_view(self):
        """Re-render the loaded volume at a new scale.

        Used when the quad view's budget changes while a tomogram is open,
        i.e. after the window has been resized. Works from the in-memory
        volume, so there is no disk read, and the crosshair and annotation
        state are left alone.
        """
        if self.volume is None or not self.xy_cache:
            return
        new_scale = self._fit_scale(self.nx, self.ny, self.nz)
        if abs(new_scale - self.scale) < 0.01:
            return
        self.scale = new_scale
        # The canvas geometry just changed, so a zoom/pan computed against
        # the old canvas size is no longer valid — reset to fit instead.
        self._reset_zoom_to_fit()
        self._apply_panel_geometry()
        self.xy_cache = prerender_z_stack(self.volume, self.contrast_window,
                                          self.xy_w, self.xy_h)
        self._render_xy()
        self._render_xz()
        self._render_yz()
        self._update_pos_label()
        # Anything prefetched was rendered at the old size — drop it and
        # refetch, or "Next" would show mis-sized images.
        with self._prefetch_lock:
            self._prefetch_cache = {}
        self._start_prefetch(self.idx + 1)

    def _fits(self, avail_h):
        # winfo_reqheight() asks Tkinter "how tall would this frame be if
        # it could have whatever size it wants?" — i.e. the natural height
        # of everything packed inside it, ignoring the window's actual
        # current size. It's used throughout this PART to check whether a
        # column's natural size still fits the room available.
        # update_idletasks() forces Tk to finish laying out any pending
        # widget changes first, since winfo_reqheight() is stale until
        # that's happened.
        self.root.update_idletasks()
        return self.right_col.winfo_reqheight() <= avail_h

    def _fit_annotation_column(self, avail_h):
        """Shrink the annotation column until it fits avail_h.

        The flag list is a fixed-height scrollable area, so it can never be
        why the column overflows. That leaves the comment box (the most
        forgiving thing to shrink) and, as a last resort, the padding
        between sections.
        """
        if self._fits(avail_h):
            return
        # 1) Comment box: 6 rows -> 3.
        self._shrink_comment(avail_h, floor=3)
        # Shrinking the comment box frees room the flag list can use, so
        # recompute its budget rather than leaving that space empty.
        self._update_flag_area_budget()
        if self._fits(avail_h):
            return
        # 2) Tighten the gaps between sections.
        for child in self.right_col.winfo_children():
            try:
                child.pack_configure(pady=1)
            except tk.TclError:
                pass
        if self._fits(avail_h):
            return
        # 3) Last resort: squeeze the comment box down to 2 rows.
        self._shrink_comment(avail_h, floor=2)
        self._update_flag_area_budget()

    def _update_flag_area_budget(self):
        """Give the flag list whatever height the rest of the column leaves.

        Sums the other sections rather than subtracting flag_outer from the
        column total — sizing it from its own height would feed back on
        itself. The floor keeps it usable on a short screen.
        """
        if self.avail_h is None:
            return
        self.root.update_idletasks()
        children  = self.right_col.winfo_children()
        padding_h = max(0, self.right_col.winfo_reqheight()
                           - sum(c.winfo_reqheight() for c in children))
        non_flag_h = padding_h + sum(c.winfo_reqheight() for c in children
                                     if c is not self.flag_outer)
        self.flag_area_h = max(180, self.avail_h - non_flag_h)
        self.flag_canvas.config(height=self.flag_area_h)

    def _tighten_columns_to_quad_height(self):
        """Re-fit the side columns to the quad view's real height.

        Runs once the first tomogram has loaded. _measure_budgets had to fit
        them to the screen instead, since the quad view's size wasn't known
        yet.
        """
        if not self.xy_cache or self.avail_h is None:
            return
        self.root.update_idletasks()
        quad_h = self.xy_h + self.xz_h + 4   # padding and borders between the rows
        new_avail_h = max(200, min(self.avail_h, quad_h))
        if new_avail_h >= self.avail_h - 1:
            return   # quad is already at least as tall as the fitted columns
        self.avail_h = new_avail_h
        self._update_flag_area_budget()
        self._fit_annotation_column(self.avail_h)
        self._fit_tomo_list(self.avail_h)

    def _shrink_comment(self, avail_h, floor):
        rows = int(self.comment_box.cget("height"))
        while rows > floor and not self._fits(avail_h):
            rows -= 1
            self.comment_box.config(height=rows)

    def _fit_tomo_list(self, avail_h):
        """Choose how many rows the tomogram list shows.

        Enough to fill the column, but not so many that the window gets
        pushed past the bottom of the screen.
        """
        self.root.update_idletasks()
        list_h  = self.tomo_list.winfo_reqheight()
        row_px  = max(12, list_h // max(1, self.list_height))
        # Everything in the list column that isn't the listbox itself.
        col_overhead = max(0, self.list_col.winfo_reqheight() - list_h)
        self.list_height = int(max(6, min(32, (avail_h - col_overhead) // row_px)))
        self.tomo_list.config(height=self.list_height)

    def _fit_tomo_list_width(self, avail_w):
        """Narrow the tomogram list until all three columns fit side by side.

        Measured with the quad view at its minimum size. Tomogram names are
        long but the tail is usually the informative part, so trimming this
        column is cheaper than shrinking the image further.
        """
        self.root.update_idletasks()
        other_w    = max(0, self.main_frame.winfo_reqwidth()
                             - self.quad_frame.winfo_reqwidth())
        min_quad_w = self.min_xy + 80          # XY at minimum + a thin YZ strip
        over = (other_w + min_quad_w) - avail_w
        if over <= 0:
            return
        chars   = int(self.tomo_list.cget("width"))
        char_px = max(6, self.tomo_list.winfo_reqwidth() // max(1, chars))
        self.tomo_list.config(width=max(14, chars - -(-over // char_px)))

    def _fit_scale(self, nx, ny, nz):
        """Voxel->pixel scale that makes the whole quad view fit its budget.

        Panel widths are  XY=nx, YZ=nz  and heights are  XZ=nz, XY=ny  (all
        times the scale), plus `gap` for the space between the panels.
        """
        gap = 4
        s_w = (self.quad_w_budget - gap) / max(1, nx + nz)
        s_h = (self.quad_h_budget - gap) / max(1, ny + nz)
        s = min(s_w, s_h)
        # Clamp the XY panel's longest side to a usable range. The floor is
        # applied last, so on a very tight budget it wins and the quad view
        # is allowed to overflow rather than become unreadable.
        longest = max(nx, ny)
        s = min(s, self.max_xy / longest)
        s = max(s, self.min_xy / longest)
        return s


    # ──────────────────────────────────────────────────────────────────────
    # PART 3 — Flag & Quality Controls
    # Builds and updates the organelle flag buttons and the quality
    # (Good/Ok/Bad/FigureQuality) segmented control.
    # ──────────────────────────────────────────────────────────────────────
    # ── Flag button creation ──
    def _add_flag_button(self, flag):
        i = len(self.flag_buttons)
        fg, bg_on = FLAG_COLORS[i % len(FLAG_COLORS)]
        self.flag_mode.setdefault(flag, "simple")
        # Each flag is a Frame, not a single Button: a main label button on
        # the left, plus a control cluster on the right that is either a
        # tiny "− +" affordance (simple mode) or a full "− N +" stepper
        # (counter mode). See _rebuild_flag_control.
        frame = tk.Frame(self.buttons_frame, bg="#F1EFE8")
        frame.columnconfigure(0, weight=1)
        main_btn = tk.Button(frame, text=flag, width=14,
                             font=("Helvetica", 13, "bold"),
                             fg=fg, bg="#F1EFE8",
                             activebackground=bg_on,
                             relief="flat", bd=0, cursor="hand2",
                             anchor="w", padx=12, pady=8,
                             command=lambda fl=flag: self._flag_primary_click(fl))
        main_btn.grid(row=0, column=0, sticky="nsew")
        control = tk.Frame(frame, bg="#F1EFE8")
        control.grid(row=0, column=1, sticky="ne", padx=(4, 6), pady=(4, 0))
        self.flag_buttons[flag] = {
            "frame": frame, "main_btn": main_btn, "control": control,
            "fg": fg, "bg_on": bg_on,
        }
        self._rebuild_flag_control(flag)
        # Always one roomy column, top to bottom; the scrollable
        # flag_canvas is what absorbs however many flags there are.
        frame.pack(fill="x", padx=2, pady=2)
        self._bind_flag_wheel(frame)

    def _rebuild_flag_control(self, flag):
        """(Re)draw one flag's corner control to match its current flag_mode.

        Called when the button is created and whenever its mode flips, but
        never on ordinary navigation — the mode doesn't change there.
        """
        info = self.flag_buttons[flag]
        control = info["control"]
        for child in control.winfo_children():
            child.destroy()
        info.pop("count_label", None)
        fg = info["fg"]
        if self.flag_mode.get(flag, "simple") == "simple":
            mini = tk.Frame(control, bg="#F1EFE8", highlightthickness=1,
                            highlightbackground="#C9C6BC")
            mini.pack()
            tk.Button(mini, text="−", font=("Helvetica", 8), fg="#888780",
                      bg="#F1EFE8", relief="flat", bd=0, width=2, cursor="hand2",
                      command=lambda fl=flag: self._enable_counter_mode(fl, -1)
                      ).pack(side="left")
            tk.Button(mini, text="+", font=("Helvetica", 8), fg="#888780",
                      bg="#F1EFE8", relief="flat", bd=0, width=2, cursor="hand2",
                      command=lambda fl=flag: self._enable_counter_mode(fl, 1)
                      ).pack(side="left")
        else:
            # The revert control is deliberately small and sits above the
            # stepper, so demoting a flag is an on-purpose click rather
            # than something that happens while adjusting its count.
            revert_row = tk.Frame(control, bg="#F1EFE8")
            revert_row.pack(fill="x")
            tk.Button(revert_row, text="× simple", font=("Helvetica", 7),
                      fg="#888780", bg="#F1EFE8", relief="flat", bd=0,
                      cursor="hand2",
                      command=lambda fl=flag: self._disable_counter_mode(fl)
                      ).pack(side="right")
            stepper_row = tk.Frame(control, bg="#F1EFE8")
            stepper_row.pack()
            tk.Button(stepper_row, text="−", font=("Helvetica", 11, "bold"),
                      fg=fg, bg="#F1EFE8", relief="flat", bd=0, width=2, cursor="hand2",
                      command=lambda fl=flag: self._step_count(fl, -1)
                      ).pack(side="left")
            count_label = tk.Label(stepper_row, text=str(self.current_state.get(flag, 0)),
                                   font=("Helvetica", 12, "bold"),
                                   fg=fg, bg="#F1EFE8", width=2)
            count_label.pack(side="left")
            info["count_label"] = count_label
            tk.Button(stepper_row, text="+", font=("Helvetica", 11, "bold"),
                      fg=fg, bg="#F1EFE8", relief="flat", bd=0, width=2, cursor="hand2",
                      command=lambda fl=flag: self._step_count(fl, 1)
                      ).pack(side="left")

    def _switch_flag_mode(self, flag, new_mode):
        """Flip one flag's mode and redraw its control.

        Shared by _enable_counter_mode and _disable_counter_mode. The flag
        list is a fixed-size scrollable area, so a control changing size
        never affects anything outside flag_canvas — nothing to refit.
        """
        self.flag_mode[flag] = new_mode
        self._rebuild_flag_control(flag)
        self._refresh_buttons()

    def _enable_counter_mode(self, flag, delta):
        """Promote a flag from toggle to counter, from its corner control.

        Sticky for the rest of the session (see flag_mode in __init__).
        Also applies the +1/-1 that was clicked, so the control does
        something immediately rather than needing a second click.
        """
        if self.flag_mode.get(flag) != "counter":
            self._switch_flag_mode(flag, "counter")
        self._step_count(flag, delta)

    def _disable_counter_mode(self, flag):
        """Demote a flag back to a plain on/off toggle, via '× simple'.

        A toggle can't hold a count above 1, so any existing count
        collapses to 0/1 — only whether it was zero is preserved. Sticky
        for the rest of the session, like promotion.
        """
        if self.flag_mode.get(flag) != "simple":
            self.current_state[flag] = 1 if self.current_state.get(flag, 0) > 0 else 0
            self._switch_flag_mode(flag, "simple")

    def _on_flag_wheel(self, event):
        num = getattr(event, "num", None)
        if num == 4:
            delta = -1
        elif num == 5:
            delta = 1
        else:
            delta = -1 if getattr(event, "delta", 0) > 0 else 1
        self.flag_canvas.yview_scroll(delta, "units")
        return "break"

    def _bind_flag_wheel(self, widget):
        """Make the wheel scroll the flag list while the pointer is over it.

        Bound per-widget and recursively into children, rather than
        globally, so it can't interfere with the wheel stepping Z/Y/X over
        the quad view panels.
        """
        for seq in ("<MouseWheel>", "<Button-4>", "<Button-5>"):
            try:
                widget.bind(seq, self._on_flag_wheel)
            except Exception:
                pass
        for child in widget.winfo_children():
            self._bind_flag_wheel(child)

    def _add_flag_from_entry(self):
        """Add a flag typed into the "Add flag" box.

        The new flag applies to every tomogram in the session and becomes a
        column in the CSV, recorded as 0 for tomograms annotated before it
        was added — indistinguishable from a flag that was checked and found
        absent. Duplicates are refused. There is no way to remove a flag from
        the GUI — edit the CSV, or start a run without it.
        """
        name = self.add_entry.get().strip()
        if not name:
            return
        if name in self.flags:
            self.status_var.set(f"'{name}' already exists.")
            return
        self.flags.append(name)
        self.current_state[name] = 0
        self._add_flag_button(name)
        self.add_entry.delete(0, "end")
        self.status_var.set(f"Added flag: {name}")

    # ── Toggle / counter / quality logic ──
    def _refresh_buttons(self):
        """Repaint the flag buttons to match current_state.

        Flag buttons only — the quality control is _refresh_quality_buttons.
        A flag with any count above 0 is drawn filled and outlined.
        """
        for flag, info in self.flag_buttons.items():
            count = self.current_state.get(flag, 0)
            if count > 0:
                info["main_btn"].config(bg=info["bg_on"], relief="solid", bd=1)
            else:
                info["main_btn"].config(bg="#F1EFE8", relief="flat", bd=0)
            count_label = info.get("count_label")
            if count_label is not None:
                count_label.config(text=str(count))

    def _flag_primary_click(self, flag):
        """Handle a click on the main body of a flag button.

        +1 for a counter flag, a plain 0/1 toggle otherwise.
        """
        if self.flag_mode.get(flag) == "counter":
            self._step_count(flag, 1)
        else:
            self.current_state[flag] = 0 if self.current_state.get(flag, 0) else 1
            self._refresh_buttons()

    def _step_count(self, flag, delta):
        """Add delta to one flag's count, floored at 0 (never negative)."""
        self.current_state[flag] = max(0, self.current_state.get(flag, 0) + delta)
        self._refresh_buttons()

    def _set_quality(self, grade):
        """Set the quality grade, or clear it by clicking the active one.

        Clicking the grade that is already selected sets quality back to
        empty. That is deliberate — it is the only way to undo a misclick —
        but it surprises people, so it is worth knowing. An empty quality is
        written to the CSV as "NA".
        """
        self.current_quality = "" if self.current_quality == grade else grade
        self._refresh_quality_buttons()

    def _refresh_quality_buttons(self):
        for grade, (btn, fg, bg_on) in self.grade_buttons.items():
            if self.current_quality == grade:
                btn.config(bg=bg_on, fg=fg)
            else:
                btn.config(bg="#F1EFE8", fg=fg)

    def _match_quality_width_to_flags(self):
        """Make the Quality control exactly as wide as the flag block below.

        flag_outer's width comes from the flag buttons' configured width,
        not from the data, so this only needs to run once.
        """
        self.root.update_idletasks()
        target_w = self.flag_outer.winfo_reqwidth()
        if target_w <= 1:
            return
        self.segment_outer.config(width=target_w,
                                  height=self.segment_outer.winfo_reqheight())
        # By default a Tk frame auto-sizes to fit whatever's packed inside
        # it — pack_propagate(False) turns that off, "locking" segment_outer
        # at the width/height just set, so it stays matched to the flag
        # block below even though its own contents (the grade buttons)
        # would otherwise size it differently.
        self.segment_outer.pack_propagate(False)


    # ──────────────────────────────────────────────────────────────────────
    # PART 4 — Crosshair Position & Panel Rendering
    # Moves the crosshair and redraws the XY/XZ/YZ panels to match it.
    # ──────────────────────────────────────────────────────────────────────
    # ── Crosshair navigation (XYZ quad view) ──
    def _clamp(self, v, lo, hi):
        return max(lo, min(v, hi))

    def _update_position(self, x=None, y=None, z=None):
        """Move the crosshair, redrawing only what has to change.

        A panel is re-rendered only if its slicing axis moved, or if its
        visible window had to shift to keep the crosshair in view (see
        _ensure_axis_visible). The others just get their crosshair lines
        redrawn, which is a cheap canvas op with no new image.
        """
        if not self.xy_cache:
            return
        old_x, old_y, old_z = self.x, self.y, self.z
        if x is not None:
            self.x = self._clamp(int(round(x)), 0, self.nx - 1)
        if y is not None:
            self.y = self._clamp(int(round(y)), 0, self.ny - 1)
        if z is not None:
            self.z = self._clamp(int(round(z)), 0, self.nz - 1)
        # While zoomed in, each axis shows only a fraction of its range. If
        # the crosshair has just moved outside that window, nudge the axis's
        # pan offset to keep it visible rather than letting it scroll
        # off-canvas. No-op at fit zoom.
        x_shifted = self._ensure_axis_visible("x")
        y_shifted = self._ensure_axis_visible("y")
        z_shifted = self._ensure_axis_visible("z")
        if self.z != old_z or x_shifted or y_shifted:
            self._render_xy()
        else:
            self._draw_xy_crosshair()
        if self.y != old_y or x_shifted or z_shifted:
            self._render_xz()
        else:
            self._draw_xz_crosshair()
        if self.x != old_x or y_shifted or z_shifted:
            self._render_yz()
        else:
            self._draw_yz_crosshair()
        self._update_pos_label()

    # ── XY panel (main): Z-slice ──
    # At fit zoom this uses the prerendered per-Z cache, so scrolling Z is
    # instant; zoomed in it renders on demand like XZ and YZ.
    def _render_xy(self):
        """Redraw the main XY panel at the current Z.

        Two routes: at fit zoom it reuses the prerendered image for this Z
        slice, which is why scrolling Z feels instant; zoomed in there is no
        cache to use, so it renders on demand like XZ and YZ.
        """
        # Tkinter Canvas images must be kept referenced somewhere in Python
        # (e.g. as self.photo_xy) for as long as they're displayed — if the
        # PhotoImage object gets garbage-collected, the canvas just goes
        # blank with no error. That's why every rendered image in this
        # class is stashed on self rather than a local variable, and why
        # xy_cache holds one PhotoImage per Z-slice rather than discarding
        # them after display.
        if self._at_fit_zoom():
            self.photo_xy = self.xy_cache[self.z]
        else:
            self.photo_xy = self._render_zoomed(
                self.volume[self.z], self.x0, self.nx, self.xy_w,
                self.y0, self.ny, self.xy_h)
        self.canvas_xy.delete("all")
        self.canvas_xy.create_image(0, 0, anchor="nw", image=self.photo_xy)
        self._draw_xy_crosshair()
        self.slice_var.set(f"Z-slice {self.z + 1} / {self.nz}")
        self.slider.set(self.z)

    def _draw_xy_crosshair(self):
        self.canvas_xy.delete("crosshair")
        px, py = (self.x - self.x0) * self.zoom, (self.y - self.y0) * self.zoom
        self.canvas_xy.create_line(px, 0, px, self.xy_h, fill="#F5E642", tags="crosshair")
        self.canvas_xy.create_line(0, py, self.xy_w, py, fill="#F5E642", tags="crosshair")

    # ── XZ panel (top-left): Y-slice, rendered on demand ──
    def _render_xz(self):
        plane = self.volume[:, self.y, :]     # (Z, X) — see AXIS ORDER in _load_current
        self.photo_xz = self._render_zoomed(
            plane, self.x0, self.nx, self.xz_w, self.z0, self.nz, self.xz_h)
        self.canvas_xz.delete("all")
        self.canvas_xz.create_image(0, 0, anchor="nw", image=self.photo_xz)
        self._draw_xz_crosshair()

    def _draw_xz_crosshair(self):
        self.canvas_xz.delete("crosshair")
        px, pz = (self.x - self.x0) * self.zoom, (self.z - self.z0) * self.zoom
        self.canvas_xz.create_line(px, 0, px, self.xz_h, fill="#F5E642", tags="crosshair")
        self.canvas_xz.create_line(0, pz, self.xz_w, pz, fill="#F5E642", tags="crosshair")

    # ── YZ panel (bottom-right): X-slice, rendered on demand ──
    def _render_yz(self):
        plane = self.volume[:, :, self.x].T   # (Y, Z) — see AXIS ORDER in _load_current
        self.photo_yz = self._render_zoomed(
            plane, self.z0, self.nz, self.yz_w, self.y0, self.ny, self.yz_h)
        self.canvas_yz.delete("all")
        self.canvas_yz.create_image(0, 0, anchor="nw", image=self.photo_yz)
        self._draw_yz_crosshair()

    def _draw_yz_crosshair(self):
        self.canvas_yz.delete("crosshair")
        pz, py = (self.z - self.z0) * self.zoom, (self.y - self.y0) * self.zoom
        self.canvas_yz.create_line(pz, 0, pz, self.yz_h, fill="#F5E642", tags="crosshair")
        self.canvas_yz.create_line(0, py, self.yz_w, py, fill="#F5E642", tags="crosshair")

    def _update_pos_label(self):
        self.pos_var.set(
            f"X {self.x:>4} / {self.nx}\nY {self.y:>4} / {self.ny}\nZ {self.z:>4} / {self.nz}"
        )


    # ──────────────────────────────────────────────────────────────────────
    # PART 5 — Pixel Size (Å/voxel)
    # Tracks and validates the Å/px value shown in the readout corner,
    # used only by the ruler tool.
    # ──────────────────────────────────────────────────────────────────────
    def _apply_pixel_size(self, header_px):
        """Set the Å/px value for a newly loaded tomogram.

        The user's override if one has been given, otherwise the MRC
        header's value. See _pixel_size_override in __init__.
        """
        self.pixel_size = self._pixel_size_override if self._pixel_size_override else header_px
        self.pixel_size_var.set(f"{self.pixel_size:.4f}")

    def _on_pixel_size_entry(self, _event=None):
        """Validate and apply a value typed into the Å/px box.

        The box is reverted to the previous value if what was typed cannot
        be read as a positive number.
        """
        text = self.pixel_size_var.get().strip()
        try:
            val = float(text)
            if val <= 0:
                raise ValueError
        except ValueError:
            self.pixel_size_var.set(f"{self.pixel_size:.4f}")
            self.status_var.set("Å/px must be a positive number — reverted.")
            return
        self.pixel_size = val
        self._pixel_size_override = val
        self.pixel_size_var.set(f"{val:.4f}")
        self.status_var.set(f"Å/px set to {val:.4f} — applied to this and future tomograms.")

    def _commit_pixel_size_entry(self):
        """Apply whatever is currently typed in the Å/px box.

        Does not wait for <Return> or <FocusOut> to fire first. Called before
        every navigation (see _commit_current): clicking a button doesn't
        reliably move focus out of the entry on every platform, so a value
        just typed would otherwise be lost. A box that isn't a positive
        number is left alone rather than reverted — _on_pixel_size_entry
        handles that when the user edits directly. Updates the status bar if
        the value changed.
        """
        text = self.pixel_size_var.get().strip()
        try:
            val = float(text)
            if val <= 0:
                raise ValueError
        except ValueError:
            return
        if val != self._pixel_size_override:
            self._pixel_size_override = val
            self.pixel_size = val
            self.status_var.set(f"Å/px set to {val:.4f} — applied to this and future tomograms.")


    # ──────────────────────────────────────────────────────────────────────
    # PART 6 — Mouse & Wheel Input, Zoom / Pan
    # Turns raw mouse/wheel events into crosshair moves, zoom, pan, or
    # a ruler measurement.
    # ──────────────────────────────────────────────────────────────────────

    # Covers the wheel handlers below, "Synced zoom/pan" after them, and
    # "Click handlers" further down. For the user-facing description of
    # what each gesture does, see the Navigation / Zoom / Measure sections
    # of the module docstring at the top of the file.
    # ── Wheel handlers: step through the axis each panel doesn't show ──
    def _step_z(self, delta):
        self._update_position(z=self.z + delta)

    def _step_y(self, delta):
        self._update_position(y=self.y + delta)

    def _step_x(self, delta):
        self._update_position(x=self.x + delta)

    def _on_slider(self, value):
        self._update_position(z=int(float(value)))

    def _panel_wheel(self, event, panel, step_fn):
        """Route a wheel event: plain wheel steps a slice, Ctrl+wheel zooms.

        The plain wheel steps the axis normal to the panel; Ctrl+wheel zooms
        about the pointer instead, synced across all three panels. Checks
        both the event's own Control bit and the separately tracked
        _ctrl_held, since that bit can arrive stale — see _build_gui.

        Tk detail: Tk reports which modifier keys were held during a mouse or
        keyboard event as a bitmask in event.state, where each bit corresponds to
        one key (Shift, Ctrl, Alt...). 0x0004 is the bit for Control, so
        `event.state & 0x0004` is nonzero exactly when Control was down.
        """
        if (event.state & 0x0004) or self._ctrl_held:    # Control held
            self._on_panel_zoom_wheel(panel, event)
        else:
            up = getattr(event, "num", None) == 4 or getattr(event, "delta", 0) > 0
            step_fn(1 if up else -1)
        return "break"

    # ── Synced zoom/pan across XY/XZ/YZ ──
    # zoom is canvas pixels per voxel, shared by all three panels. x0/y0/z0
    # are the voxel coordinate at each canvas's top-left corner along that
    # axis; XY uses (x0,y0), XZ uses (x0,z0), YZ uses (z0,y0) — so panning
    # or zooming in any one panel moves the other two panels that share an
    # axis with it too, which is what keeps the three views in sync.
    def _at_fit_zoom(self):
        """True when the whole volume fits the panels, i.e. not zoomed in.

        "Fit zoom" is the scale at which a tomogram first appears, and the
        minimum the user can zoom out to. Several shortcuts are only valid
        here — most importantly the prerendered Z-slice cache.
        """
        return self.zoom <= self.fit_zoom + 1e-6

    def _reset_zoom_to_fit(self):
        """Zoom back out to the whole volume and forget any pan.

        Called on every load and every resize, because a zoom or pan worked
        out against the old canvas size means nothing at the new one.
        """
        self.fit_zoom = self.scale
        self.zoom     = self.fit_zoom
        self.x0 = self.y0 = self.z0 = 0.0
        self.MAX_ZOOM = max(self.fit_zoom * 15.0, self.fit_zoom + 20.0)

    def _clamp_offsets(self):
        """Keep the pan offsets inside the volume, on all three axes.

        Without this, dragging past an edge would scroll empty black space
        into view and let the crosshair coordinates drift outside the data.
        """
        def clamp_axis(off, n, canvas_len):
            if self.zoom <= 0:
                return 0.0
            visible = canvas_len / self.zoom
            return self._clamp(off, 0.0, max(0.0, n - visible))
        self.x0 = clamp_axis(self.x0, self.nx, self.xy_w)
        self.y0 = clamp_axis(self.y0, self.ny, self.xy_h)
        self.z0 = clamp_axis(self.z0, self.nz, self.xz_h)

    def _axis_geometry(self, axis):
        """(voxel count, canvas pixels, crosshair value) for one shared axis.

        X uses xy_w (== xz_w), Y uses xy_h (== yz_h) and Z uses xz_h
        (== yz_w); those equalities hold because all three panels share one
        scale — see panel_dims.
        """
        if axis == "x":
            return self.nx, self.xy_w, self.x
        if axis == "y":
            return self.ny, self.xy_h, self.y
        return self.nz, self.xz_h, self.z

    def _ensure_axis_visible(self, axis):
        """Scroll one axis so the crosshair stays inside the visible window.

        Returns True if it had to scroll, False if the crosshair was already
        in view — always the case at fit zoom, which shows the whole axis.
        This matters because the crosshair is not just a marker: it selects
        what each panel displays. XY shows the Z slice it sits on, XZ the Y
        slice, YZ the X slice. Zoomed in, each panel shows only part of an
        axis, and the pan offset is per axis, shared by the two panels using
        it. So if Z scrolls out of view, XY shows Z=175 while XZ beside it
        shows Z 0-30 — two panels claiming to be the same place and showing
        different ones. Scrolling back onto the crosshair keeps all three
        intersecting at one point.
        """
        if self._at_fit_zoom():
            return False
        n, canvas_len, val = self._axis_geometry(axis)
        off = getattr(self, f"{axis}0")
        visible = canvas_len / self.zoom
        if off <= val <= off + visible:
            return False
        new_off = self._clamp(val - visible / 2.0, 0.0, max(0.0, n - visible))
        if abs(new_off - off) < 1e-9:
            return False
        setattr(self, f"{axis}0", new_off)
        return True

    def _refresh_zoomed_views(self):
        if not self.xy_cache:
            return
        self._render_xy()
        self._render_xz()
        self._render_yz()

    def _render_zoomed(self, plane2d, off_a, n_a, canvas_w, off_b, n_b, canvas_h):
        """Crop a plane to the current zoom/pan window and scale it to the canvas.

        plane2d has shape (n_b, n_a): rows are the vertical (b) axis,
        columns the horizontal (a) one. At fit zoom the crop is the whole
        plane.
        """
        src_w = max(1, min(n_a, int(math.ceil(canvas_w / self.zoom))))
        src_h = max(1, min(n_b, int(math.ceil(canvas_h / self.zoom))))
        start_a = int(self._clamp(round(off_a), 0, max(0, n_a - src_w)))
        start_b = int(self._clamp(round(off_b), 0, max(0, n_b - src_h)))
        sub = plane2d[start_b:start_b + src_h, start_a:start_a + src_w]
        img = array_to_image(sub, self.contrast_window, canvas_w, canvas_h)
        return ImageTk.PhotoImage(img)

    def _on_panel_zoom_wheel(self, panel, event):
        if not self.xy_cache:
            return
        factor = 1.25 if (getattr(event, "num", None) == 4 or getattr(event, "delta", 0) > 0) else 1 / 1.25
        self._zoom_at(panel, factor, event.x, event.y)
    # The axis each panel does NOT display: the one _zoom_at has to
    # recentre, since zooming via that panel never touches it directly.
    _THIRD_AXIS = {"xy": "z", "xz": "y", "yz": "x"}

    def _zoom_at(self, panel, factor, cx, cy):
        """Zoom about the pointer, keeping the voxel under it in place.

        The trick is to work out which voxel is under the cursor first, then
        choose the new pan offset so that voxel lands back under the cursor
        at the new zoom — that is what makes zoom feel anchored rather than
        jumping to the panel centre.

        All three panels share one zoom, so this also has to recentre the
        axis this panel does not show (see _THIRD_AXIS).
        """
        old_zoom = self.zoom
        a_before, b_before = self._panel_to_voxel(panel, cx, cy)
        new_zoom = self._clamp(old_zoom * factor, self.fit_zoom, self.MAX_ZOOM)
        if abs(new_zoom - old_zoom) < 1e-9:
            return
        self.zoom = new_zoom
        a_off = a_before - cx / self.zoom
        b_off = b_before - cy / self.zoom
        if panel == "xy":
            self.x0, self.y0 = a_off, b_off
        elif panel == "xz":
            self.x0, self.z0 = a_off, b_off
        else:
            self.z0, self.y0 = a_off, b_off
        self._clamp_offsets()
        # The axis this panel doesn't show keeps a now-stale offset, since
        # its window just shrank around that offset rather than around the
        # crosshair. Recentre it, so the two panels sharing that axis still
        # show where you are — see _ensure_axis_visible.
        self._ensure_axis_visible(self._THIRD_AXIS[panel])
        self._clamp_offsets()
        self._refresh_zoomed_views()

    def _on_pan_press(self, panel, event):
        if not self.xy_cache:
            return
        self._pan_start = (panel, event.x, event.y, self.x0, self.y0, self.z0)

    def _on_pan_drag(self, event):
        if self._pan_start is None:
            return
        panel, sx, sy, x0_0, y0_0, z0_0 = self._pan_start
        dx, dy = (event.x - sx) / self.zoom, (event.y - sy) / self.zoom
        if panel == "xy":
            self.x0, self.y0 = x0_0 - dx, y0_0 - dy
        elif panel == "xz":
            self.x0, self.z0 = x0_0 - dx, z0_0 - dy
        else:
            self.z0, self.y0 = z0_0 - dx, y0_0 - dy
        self._clamp_offsets()
        self._refresh_zoomed_views()

    def _on_pan_release(self, _event=None):
        self._pan_start = None

    # ── Click handlers ──
    # A click jumps the crosshair to that point; a left-button drag draws
    # the ruler instead (see _on_panel_motion).
    DRAG_THRESHOLD = 4   # canvas pixels below which a press+release is a click, not a drag

    def _panel_canvas(self, panel):
        return {"xy": self.canvas_xy, "xz": self.canvas_xz, "yz": self.canvas_yz}[panel]

    def _panel_to_voxel(self, panel, cx, cy):
        if panel == "xy":
            return self.x0 + cx / self.zoom, self.y0 + cy / self.zoom
        if panel == "xz":
            return self.x0 + cx / self.zoom, self.z0 + cy / self.zoom
        return self.z0 + cx / self.zoom, self.y0 + cy / self.zoom   # yz

    def _on_panel_press(self, panel, event):
        if not self.xy_cache:
            return
        self._commit_pixel_size_entry()
        self._drag_start = (panel, event.x, event.y)
        self._drag_is_measuring = False

    def _on_panel_motion(self, panel, event):
        if self._drag_start is None or self._drag_start[0] != panel:
            return
        _, sx, sy = self._drag_start
        if not self._drag_is_measuring:
            if abs(event.x - sx) < self.DRAG_THRESHOLD and abs(event.y - sy) < self.DRAG_THRESHOLD:
                return
            self._drag_is_measuring = True
        canvas = self._panel_canvas(panel)
        for c in (self.canvas_xy, self.canvas_xz, self.canvas_yz):
            c.delete("measure")
        canvas.create_line(sx, sy, event.x, event.y, fill="#6EC62F", width=4, tags="measure")
        dist_px  = math.hypot(event.x - sx, event.y - sy)
        dist_ang = (dist_px / self.zoom) * self.pixel_size
        self.measured_var.set(f"Measured size:\n{dist_ang:,.1f} Å  ({dist_ang / 10:,.2f} nm)")

    def _on_panel_release(self, panel, event):
        drag = self._drag_start
        was_measuring = self._drag_is_measuring
        self._drag_start = None
        self._drag_is_measuring = False
        if drag is None or drag[0] != panel:
            return
        if was_measuring:
            return   # line + "Measured size" readout stay as-is; nothing else happens
        if not self.xy_cache:
            return
        a, b = self._panel_to_voxel(panel, event.x, event.y)
        if panel == "xy":
            self._update_position(x=a, y=b)
        elif panel == "xz":
            self._update_position(x=a, z=b)
        else:
            self._update_position(z=a, y=b)


    # ──────────────────────────────────────────────────────────────────────
    # PART 7 — Tomogram Navigation & Loading
    # Saving the current tomogram's annotation, moving to the
    # next/previous one, and loading a volume from disk.
    # ──────────────────────────────────────────────────────────────────────
    # ── Persisting current tomo's state ──
    def _commit_current(self):
        """Store the current tomogram's annotation and rewrite the outputs."""
        # Catch a just-typed Å/px value before focus-event ordering can
        # lose it (see _commit_pixel_size_entry).
        self._commit_pixel_size_entry()
        # The full filename, not the shortened display name: the CSV has to
        # identify a file on disk, and the display name is relative to the
        # rest of the selection (see shorten_display_names).
        name = self.tomo_paths[self.idx].name
        entry = dict(self.current_state)
        entry["quality"] = self.current_quality
        entry["comment"] = self.comment_box.get("1.0", "end").strip()
        self.data[name]  = entry
        self._save_outputs()

    # ── Tomogram list interaction ──
    def _on_list_select(self, *_):
        """Jump to the tomogram clicked in the list, saving the current one.

        The _list_updating guard matters: _sync_list changes the selection
        programmatically, which makes Tk fire this same event, which would
        navigate again — an endless loop.
        """
        if self._list_updating:
            return
        sel = self.tomo_list.curselection()
        if not sel or sel[0] == self.idx:
            return
        self._commit_current()
        self.idx = sel[0]
        self._load_current()

    def _sync_list(self):
        """Highlight the current tomogram in the list and scroll it into view.

        Sets _list_updating while it works, so the selection change it makes
        does not come back as a user click — see _on_list_select.
        """
        self._list_updating = True
        self.tomo_list.selection_clear(0, "end")
        self.tomo_list.selection_set(self.idx)
        self.tomo_list.see(self.idx)
        self._list_updating = False

    # ── Loading a tomogram ──
    def _close_volume(self):
        """Drop the loaded volume and clear everything that depended on it."""
        self.volume   = None
        self.xy_cache = []
        for canvas in (self.canvas_xy, self.canvas_xz, self.canvas_yz):
            canvas.delete("all")
        self.pos_var.set("")
        # The ruler is a display-only overlay, so the last measurement
        # goes with the view it was taken on.
        self.measured_var.set("Measured size: —")
        self._drag_start = None
        self._drag_is_measuring = False
        self._pan_start = None

    def _load_current(self):
        """Show the tomogram at self.idx: its annotation, then its images.

        Two halves. First the annotation side — name, saved flags, quality
        and comment are restored from self.data, so revisiting a tomogram
        shows what you last recorded. Then the image side, by one of two
        routes:

          - a prefetch worker already loaded it (PART 8): only the final
            image conversion is left, and it appears immediately
          - it did not: load from disk, compute contrast, and render every
            Z slice, updating the status line as it goes

        Either way the crosshair starts at the centre of the volume and the
        next tomograms start prefetching. Running out of memory is handled
        rather than fatal: the tomogram is dropped, prefetching steps down,
        and the user is told what to do about it.
        """
        total = len(self.tomo_paths)
        self.progress_var.set(f"{self.idx + 1} / {total}")
        path  = self.tomo_paths[self.idx]
        # Two names, deliberately: the shortened one is what the user reads,
        # the filename is what the annotation is filed under.
        self.name_var.set(self._tomo_name(path))
        saved = self.data.get(path.name, {})
        self.current_state   = {flag: saved.get(flag, 0) for flag in self.flags}
        self.current_quality = saved.get("quality", "")
        self._refresh_buttons()
        self._refresh_quality_buttons()
        self.comment_box.delete("1.0", "end")
        if saved.get("comment"):
            self.comment_box.insert("1.0", saved["comment"])
        self.prev_btn.config(state="normal" if self.idx > 0 else "disabled")
        self.next_btn.config(text="Next →" if self.idx < total - 1 else "Finish ✓")
        self._sync_list()
        self._close_volume()
        if not path.exists():
            self.status_var.set(f"File not found: {path}")
            return
        cached = self._take_prefetched(self.idx)
        if cached is not None:
            self.volume = cached["volume"]
            self.nz, self.ny, self.nx = cached["nz"], cached["ny"], cached["nx"]   # see AXIS ORDER below
            self.scale = cached["scale"]
            self.contrast_window = cached["contrast_window"]
            self._apply_panel_geometry()
            self._reset_zoom_to_fit()
            self._apply_pixel_size(cached.get("pixel_size", 1.0))
            self.x, self.y, self.z = self.nx // 2, self.ny // 2, self.nz // 2
            self.status_var.set("Loading prefetched volume…")
            self.root.update()
            self.xy_cache = [ImageTk.PhotoImage(img) for img in cached["pil_slices"]]
            self.status_var.set(
                f"Ready — {self.nz} Z-slices (prefetched). Click/scroll to navigate."
            )
        else:
            self.status_var.set("Loading volume…")
            self.root.update()
            try:
                self.volume, header_px = load_volume(path)
            except MemoryError:
                # numpy/mrcfile raise a bare MemoryError with no message, so
                # the generic handler below would show "Failed to load: "
                # and give no clue what went wrong.
                self.status_var.set(
                    "Out of memory while loading this tomogram — close other "
                    "programs and try again, or annotate it on a machine "
                    "with more RAM."
                )
                return
            except Exception as e:
                self.status_var.set(f"Failed to load: {e}")
                return
            # numpy detail — AXIS ORDER. This catches everyone once. An
            # MRC volume is indexed volume[z, y, x], slowest axis first,
            # so .shape is (nz, ny, nx) — the REVERSE of the x, y, z
            # order used everywhere else in this class. Hence:
            #
            #     volume[z]           one XY slice  (the main panel)
            #     volume[:, y, :]     one XZ slice  (top-left panel)
            #     volume[:, :, x]     one YZ slice  (bottom-right panel)
            #
            # If a panel ever shows the wrong plane, check this first.
            #
            # One shared voxel->pixel scale for all three panels, so the
            # aspect ratio stays real — see panel_dims.
            self.nz, self.ny, self.nx = self.volume.shape
            self.scale = self._fit_scale(self.nx, self.ny, self.nz)
            self._apply_panel_geometry()
            # First real layout: check it fits and rescale if the pre-load
            # estimate was optimistic. Cheap here — nothing is rendered yet.
            if not self._calibrated:
                self._calibrated = True
                if self._calibrate_fit():
                    self.scale = self._fit_scale(self.nx, self.ny, self.nz)
                    self._apply_panel_geometry()
            self._reset_zoom_to_fit()
            self._apply_pixel_size(header_px)
            self.x, self.y, self.z = self.nx // 2, self.ny // 2, self.nz // 2
            self.status_var.set("Computing contrast…")
            self.root.update()
            try:
                self.contrast_window = compute_contrast_window(self.volume)

                def report_progress(done, total):
                    self.status_var.set(f"Processing tomogram… {done} / {total} Z-slices")
                    self.root.update()
                self.xy_cache = prerender_z_stack(
                    self.volume, self.contrast_window, self.xy_w, self.xy_h,
                    progress_callback=report_progress
                )
            except MemoryError:
                # Rendering is the memory peak: raw volume plus every
                # rendered Z-slice at once. This failed on the tomogram
                # being viewed, so drop it and stop prefetching too.
                self._close_volume()
                self._step_down_prefetch()
                self.status_var.set(
                    "Out of memory while rendering this tomogram — close other "
                    "programs and try again, or annotate it on a machine with "
                    "more RAM."
                )
                return
            self.status_var.set(
                f"Ready — {self.nz} Z-slices. Click/scroll the XY, XZ, YZ panels to navigate."
            )
        self._render_xy()
        self._render_xz()
        self._render_yz()
        self._update_pos_label()
        self._start_prefetch(self.idx + 1)

    # ── Navigation actions ──
    def _go_next(self):
        self._commit_current()
        if self.idx < len(self.tomo_paths) - 1:
            self.idx += 1
            self._load_current()
        else:
            self.status_var.set("Last tomogram — all annotations saved.")

    def _go_previous(self):
        self._commit_current()
        if self.idx > 0:
            self.idx -= 1
            self._load_current()

    def _on_quit(self):
        self._commit_current()
        self._close_volume()
        self.root.destroy()


    # ──────────────────────────────────────────────────────────────────────
    # PART 8 — Background Prefetch
    # Loads upcoming tomograms on worker threads while the current one
    # is being annotated, so "Next" is instant.
    # ──────────────────────────────────────────────────────────────────────

    # Threading: the load, contrast and per-slice resize are all numpy/PIL
    # and never touch Tk. Only the PIL -> PhotoImage step must happen on
    # the main thread, since it talks to Tcl, so each worker stops short of
    # it and _load_current finishes that quickly.
    #
    # Running out of memory here reduces prefetch_ahead instead of crashing
    # (see PREFETCH_AHEAD_DEFAULT); at 0 loading is fully synchronous.
    def _step_down_prefetch(self):
        """Permanently reduce how far ahead tomograms are prefetched.

        Steps down to a floor of 0, which is fully synchronous loading. It
        never steps back up, so a machine already under memory pressure does
        not keep re-triggering the same problem on every "Next".

        Threading: must run on the main thread — it touches Tk state and the
        shared cache.
        """
        if self.prefetch_ahead <= 0:
            return
        self.prefetch_ahead -= 1
        with self._prefetch_lock:
            keep = set(range(self.idx, self.idx + self.prefetch_ahead + 1))
            for idx in list(self._prefetch_cache):
                if idx not in keep:
                    del self._prefetch_cache[idx]
        if self.prefetch_ahead == 0:
            self.status_var.set("Low memory: prefetch disabled — loading tomograms synchronously.")
        else:
            self.status_var.set(f"Low memory: now prefetching only {self.prefetch_ahead} tomogram(s) ahead.")

    def _start_prefetch(self, base_idx):
        """Start background loads for the next few tomograms.

        Loads up to prefetch_ahead of them, starting at base_idx. Indices
        already cached or in flight are skipped, so this is safe to call on
        every navigation.
        """
        for offset in range(self.prefetch_ahead):
            idx = base_idx + offset
            if idx < 0 or idx >= len(self.tomo_paths):
                continue
            with self._prefetch_lock:
                already_have = idx in self._prefetch_targets or idx in self._prefetch_cache
            if already_have:
                continue
            with self._prefetch_lock:
                self._prefetch_targets.add(idx)
            threading.Thread(target=self._prefetch_worker, args=(idx,), daemon=True).start()

    def _prefetch_worker(self, idx):
        path = self.tomo_paths[idx]
        try:
            if not path.exists():
                return
            volume, header_px = load_volume(path)
            nz, ny, nx = volume.shape          # see AXIS ORDER in _load_current
            scale = self._fit_scale(nx, ny, nz)
            xy_w, xy_h, xz_w, xz_h, yz_w, yz_h = panel_dims(nx, ny, nz, scale)
            contrast_window = compute_contrast_window(volume)
            pil_slices = [array_to_image(volume[z], contrast_window, xy_w, xy_h)
                          for z in range(nz)]
            result = dict(volume=volume, contrast_window=contrast_window,
                          pil_slices=pil_slices, nz=nz, ny=ny, nx=nx, scale=scale,
                          xy_w=xy_w, xy_h=xy_h, xz_w=xz_w, xz_h=xz_h,
                          yz_w=yz_w, yz_h=yz_h, pixel_size=header_px)
            with self._prefetch_lock:
                self._prefetch_cache[idx] = result   # keep at most prefetch_ahead entries
        except MemoryError:
            # Discard any partial result for this index, and step the
            # prefetch depth down on the main thread — Tk state isn't safe
            # to touch from a worker.
            with self._prefetch_lock:
                self._prefetch_cache.pop(idx, None)
            try:
                self.root.after(0, self._step_down_prefetch)
            except Exception:
                # A worker can only schedule work like this once the main
                # thread has entered root.mainloop(), and the very first
                # prefetch starts before _build_gui gets there. A
                # MemoryError landing in that tiny window would raise
                # "main thread is not in main loop". Skipping it is
                # harmless: the next prefetch attempt still reports.
                pass
        except Exception:
            pass  # prefetch is a best-effort speedup; _load_current falls back to a normal load
        finally:
            with self._prefetch_lock:
                self._prefetch_targets.discard(idx)

    def _take_prefetched(self, idx):
        with self._prefetch_lock:
            return self._prefetch_cache.pop(idx, None)


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 5 — The 'link' subcommand
# Symlinks every tomogram matching a set of filters into one folder,
# so a subset can be handed to another program.
# ══════════════════════════════════════════════════════════════════════════════
def soft_link_object(csv_path, object_names, input_dir):
    """Read a TOMATO CSV and soft-link every tomogram that satisfies ALL filters.

    Each name in object_names is checked against:
      - the 'quality' column  (e.g. "Good", "Bad") if the name is a known grade
      - the flag column of that name  (count > 0) otherwise
    All conditions are combined with AND — a tomogram must match every filter.
    Links are placed in  ./<name1>_<name2>_..._tomograms/

    Existing links in that folder are replaced, so re-running with the same
    filters is safe. A missing source file is warned about and skipped, not
    fatal.

    tomo_name is the filename, so there is nothing to rebuild and no
    suffix option here. A CSV that stores a name without its suffix still
    works: the file is found by looking for a single filename starting with
    that name.
    """
    csv_path  = Path(csv_path)
    input_dir = Path(input_dir)
    if not csv_path.exists():
        print(f"CSV not found: {csv_path}")
        sys.exit(1)
    if not input_dir.exists():
        print(f"Input directory not found: {input_dir}")
        sys.exit(1)
    matching = []
    try:
        with open(csv_path, newline="") as f:
            reader = csv.DictReader(f)
            headers = list(reader.fieldnames or [])
            if "tomo_name" not in headers:
                csv_error(csv_path, "there is no 'tomo_name' column", headers)
            flag_cols = [h for h in headers if h not in RESERVED_COLS]
            for obj in object_names:
                if obj not in GRADES and obj not in headers:
                    print(f"'{obj}' not found in {csv_path}.")
                    print(f"Available flags   : {', '.join(flag_cols)}")
                    print(f"Available grades  : {', '.join(GRADES)}")
                    sys.exit(1)
            for row in reader:
                name = row.get("tomo_name", "").strip()
                if not name:
                    continue
                hit = True
                for obj in object_names:
                    if obj in GRADES:
                        if row.get("quality", "").strip() != obj:
                            hit = False
                            break
                    else:
                        if parse_count(row.get(obj)) <= 0:
                            hit = False
                            break
                if hit:
                    matching.append(name)
    except OSError as e:
        csv_error(csv_path, f"the file could not be opened ({e.strerror})")
    except UnicodeDecodeError:
        csv_error(csv_path, "the file is not text — a binary file named .csv?")
    except csv.Error as e:
        csv_error(csv_path, f"the file is not valid CSV ({e})")
    label   = "_".join(object_names)
    out_dir = Path(f"{label}_tomograms")
    if not matching:
        print(f"No tomograms found matching all of: {', '.join(object_names)}")
        sys.exit(0)
    out_dir.mkdir(exist_ok=True)
    created = 0
    for name in sorted(matching):
        src = input_dir / name
        if not src.exists():
            # The stored name may be missing its suffix. Accept a match
            # only when exactly one file starts with it, so an ambiguous
            # name is reported rather than guessed at.
            candidates = sorted(p for p in input_dir.glob(name + "*") if p.is_file())
            if len(candidates) == 1:
                src = candidates[0]
            elif len(candidates) > 1:
                print(f"  Warning: '{name}' matches {len(candidates)} files — skipped")
                continue
            else:
                print(f"  Warning: source not found — {src}")
                continue
        link = out_dir / src.name
        if link.is_symlink() or link.exists():
            link.unlink()
        os.symlink(src.resolve(), link)
        print(f"  {link.name}")
        created += 1
    print(f"\n{created} symlink(s) created in ./{out_dir}/")


# ══════════════════════════════════════════════════════════════════════════════
# SECTION 6 — Command line and entry point
# Argument parsing for 'annotate' and 'link', and the startup
# sequence: resume, seed from a previous CSV, open the GUI.
# ══════════════════════════════════════════════════════════════════════════════
def _main_annotate():
    """Handle `TOMATO.py annotate`: parse the options, then open the GUI.

    The startup order matters and is easy to misread:

      1. flags = the base set, plus anything from --flags
      2. --previous_csv, if given, adds its flag columns and its
         annotations, and becomes the file written back to unless
         --output names another one
      3. with no --previous_csv, an existing output file is loaded
         instead, so an interrupted run resumes where it stopped
      4. any flag found with a count above 1 is opened as a counter, since
         only a counter could have produced that number

    Nothing is written until the user navigates or quits.
    """
    parser = argparse.ArgumentParser(
        prog="TOMATO.py annotate",
        description="Open the annotation GUI — load tomograms one by one, toggle flags, grade quality.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--input", required=True,
                        help="Path to folder containing tomogram .mrc files (REQUIRED)")
    parser.add_argument("--flags", nargs="+", default=[],
                        help="Extra organelle button labels, added alongside the "
                             f"always-present base set {BASE_FLAGS}")
    parser.add_argument("--output", default=None,
                        help=f"Output path for the annotation CSV; a matching .txt "
                             f"compact summary is written alongside with the same stem. "
                             f"Defaults to --previous_csv when that is given, so "
                             f"annotation continues in the same file; otherwise to "
                             f"{DEFAULT_OUTPUT} (→ annotations.csv + annotations.txt)")
    parser.add_argument("--suffix", default=DEFAULT_SUFFIX,
                        help=f"Only annotate files ending with this "
                             f"(default: {DEFAULT_SUFFIX})")
    parser.add_argument("--prefix", default="",
                        help="Only annotate files starting with this, e.g. "
                             "--prefix TS_ . Combines with --suffix, so the two "
                             "together select <prefix>*<suffix>. Optional — by "
                             "default every file matching --suffix is opened. "
                             "Neither flag changes what is saved: the CSV always "
                             "records the full filename.")
    parser.add_argument("--previous_csv", default=None,
                        help="Continue a specific TOMATO CSV instead of the default "
                             "one: its flag columns and annotations are loaded, and "
                             "annotation is written back to it. Add --output to write "
                             "somewhere else and leave the original untouched.")
    parser.add_argument("--pixel_size", type=float, default=None,
                        help="Å/voxel to use instead of trusting each MRC header's "
                             "voxel size (which is sometimes wrong or gets copied "
                             "incorrectly between programs). Optional — by default "
                             "the header value is used and shown/editable in the "
                             "GUI's position readout. Only used for the ruler tool's "
                             "size output; navigation is unaffected either way.")
    args = parser.parse_args(sys.argv[2:])
    tomo_dir = Path(args.input)
    if not tomo_dir.exists():
        print(f"Directory not found: {tomo_dir}")
        sys.exit(1)
    flags = list(BASE_FLAGS)
    for f in args.flags:
        if f not in flags:
            flags.append(f)
    initial_data = {}
    if args.previous_csv:
        prev_path = Path(args.previous_csv)
        # load_csv_with_flags reports the CSV's own flag columns, so a flag
        # added in an earlier session comes back as a button instead of
        # being dropped and then erased on the next save.
        csv_flags, initial_data = load_csv_with_flags(prev_path)
        for f in csv_flags:
            if f not in flags:
                flags.append(f)
        print(f"Previous CSV   : {prev_path} ({len(initial_data)} tomograms loaded)")
    # Where annotation is written. --previous_csv means "continue this file",
    # so it is also the output unless --output names somewhere else.
    if args.output:
        output_path = Path(args.output)
    elif args.previous_csv:
        output_path = Path(args.previous_csv)
    else:
        output_path = Path(DEFAULT_OUTPUT)
    txt_path    = output_path.with_suffix(".txt")
    # Only when no --previous_csv was given: that flag already said which
    # file to continue, and loading the output on top of it would mix two
    # sources into one run.
    if not args.previous_csv and output_path.exists():
        csv_flags, resume_data = load_csv_with_flags(output_path)
        for f in csv_flags:
            if f not in flags:
                flags.append(f)
        if resume_data:
            initial_data.update(resume_data)
            print(f"Auto-resuming  : {output_path} ({len(resume_data)} tomograms)")
    pattern = f"{args.prefix}*{args.suffix}"
    tomo_paths = sorted(tomo_dir.glob(pattern))
    if not tomo_paths:
        print(f"No files matching {pattern} found in {tomo_dir}")
        sys.exit(1)
    # A CSV may name its tomograms without the suffix, where this writes
    # the full filename. Rename any such row that matches a file in this
    # run, so the annotations still resume instead of silently starting
    # over and leaving both spellings in the CSV.
    full_names = {p.name for p in tomo_paths}
    stripped_to_full = {}
    for p in tomo_paths:
        stem = (p.name[: -len(args.suffix)]
                if args.suffix and p.name.endswith(args.suffix) else p.name)
        stripped_to_full[stem] = p.name
    renamed = 0
    for old_name in list(initial_data):
        if old_name not in full_names and old_name in stripped_to_full:
            initial_data[stripped_to_full[old_name]] = initial_data.pop(old_name)
            renamed += 1
    if renamed:
        print(f"Updated names  : {renamed} row(s) from an older TOMATO version "
              f"now stored under their full filename")
    print(f"Tomo directory : {tomo_dir}")
    print(f"Selecting      : {pattern}")
    print(f"Flags          : {flags}")
    print(f"Output CSV     : {output_path}")
    print(f"Output TXT     : {txt_path}")
    print(f"Tomogram count : {len(tomo_paths)}")
    # A value above 1 can only have come from a counter (a plain toggle
    # stores 0/1), so promote that flag and open the GUI with its stepper
    # rather than silently clamping the existing counts back to 0/1.
    initial_flag_mode = {}
    for entry in initial_data.values():
        for flag in flags:
            if entry.get(flag, 0) > 1:
                initial_flag_mode[flag] = "counter"
    print()
    Annotator(tomo_paths, flags, output_path, txt_path, initial_data,
              initial_flag_mode=initial_flag_mode, pixel_size=args.pixel_size)


def _main_link():
    parser = argparse.ArgumentParser(
        prog="TOMATO.py link",
        description="Soft-link tomograms that match an object or quality grade from a TOMATO CSV.",
    )
    parser.add_argument("--csv", required=True,
                        help="TOMATO CSV to read annotations from (REQUIRED)")
    parser.add_argument("--object", required=True, nargs="+",
                        help="One or more object names / quality grades to filter by (AND logic). "
                             "E.g. --object Good Mito  → tomograms rated Good that also contain Mito. "
                             "Output directory is named after all values joined with '_'.")
    parser.add_argument("--input", required=True,
                        help="Directory containing the tomogram .mrc files (REQUIRED)")
    args = parser.parse_args(sys.argv[2:])
    soft_link_object(
        csv_path     = args.csv,
        object_names = args.object,
        input_dir    = args.input,
    )


SUBCOMMANDS = {
    "annotate": (_main_annotate, "Open the annotation GUI"),
    "link":     (_main_link,     "Soft-link tomograms matching an object or quality grade"),
}


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in SUBCOMMANDS:
        print("usage: TOMATO.py <subcommand> [options]\n")
        print("Subcommands:")
        for name, (_, desc) in SUBCOMMANDS.items():
            print(f"  {name:<12}{desc}")
        print("\nRun  python TOMATO.py <subcommand> --help  for options.")
        sys.exit(0 if len(sys.argv) == 1 else 1)
    SUBCOMMANDS[sys.argv[1]][0]()


if __name__ == "__main__":
    main()
