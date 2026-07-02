# ******************************************************************************
#
# TOMATO - TOMogram Annotation TOol
# Copyright (C) 2026 Simon J. Alkema
#
# Author: Simon J. Alkema with contributions from Euan W. Pyle and Anastasiia Babenko
# EMBL Imaging Centre, Heidelberg, Germany - Mattei Lab
# contact simon.alkema@embl.de for questions
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

#!/usr/bin/env python3
"""
TOMATO-style tomogram annotator — mrcfile + in-GUI slice scroller.

Loads each .mrc directly via mrcfile and lets you
scroll through Z-slices inside the same window as the annotation buttons.
Mito / ER / Golgi are always present as base organelles; extra organelles
can be added with --flags (without removing the base three), or added
live inside the GUI at any time using the "Add flag" entry box.

Auto-resume: if the output CSV already exists when the tool starts, prior
annotations are loaded automatically.  Pass --previous_csv to seed from a
different run's CSV (new output is still written to --output).

Outputs:
  <output>.csv  — full annotation table (True/False per flag + quality + comment)
  <output>.txt  — compact summary: TS_0001,Good,Mito,ER,Shell  (one line per tomo)

Usage:
    python TOMATO.py annotate --input /path/to/tomograms
                              [--flags Lysosome Nucleus ...]
                              [--output annotations.csv]
                              [--suffix .mrc]
                              [--previous_csv prior_annotations.csv]

    python TOMATO.py link --previous_csv annotations.csv
                          --input /path/to/tomograms
                          --object Mito (one or more options available)
                          [--suffix .mrc]

Scrolling:
    Mouse wheel over the image, or the slider below it, moves through
    Z-slices. Slice position is independent per tomogram (resets to
    the middle slice when a new tomogram loads).
"""

import argparse
import csv
import os
import sys
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

# ── Configuration ─────────────────────────────────────────────────────────────

BASE_FLAGS      = ["Mito", "ER", "Golgi"]
GRADES          = ["FigureQuality", "Good", "Ok", "Bad"]
DEFAULT_SUFFIX  = ".mrc"
DEFAULT_OUTPUT  = "annotations.csv"
CANVAS_SIZE     = 600

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

# ── CSV / TXT helpers ─────────────────────────────────────────────────────────

RESERVED_COLS = {"tomo_name", "quality", "comment"}


def load_csv(csv_path, flags):
    """Returns dict: tomo_name -> {flag: bool, 'quality': str, 'comment': str}"""
    data = {}
    if not csv_path.exists():
        return data
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row.get("tomo_name")
            if not name:
                continue
            entry = {flag: row.get(flag, "False") == "True" for flag in flags}
            entry["quality"] = row.get("quality", "")
            entry["comment"] = row.get("comment", "")
            data[name] = entry
    return data


def load_csv_with_flags(csv_path):
    """Read a previously written TOMATO CSV.

    Returns (flags_list, data_dict). flags_list is the ordered list of flag
    column names (tomo_name, quality, comment excluded).
    """
    if not csv_path.exists():
        print(f"Warning: --previous_csv path not found: {csv_path}")
        return [], {}
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        headers = list(reader.fieldnames or [])
        flags_from_csv = [h for h in headers if h not in RESERVED_COLS]
        data = {}
        for row in reader:
            name = row.get("tomo_name")
            if not name:
                continue
            entry = {flag: row.get(flag, "False") == "True" for flag in flags_from_csv}
            entry["quality"] = row.get("quality", "")
            entry["comment"] = row.get("comment", "")
            data[name] = entry
    return flags_from_csv, data


def write_csv(csv_path, flags, all_data):
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["tomo_name", "quality"] + flags + ["comment"])
        for name in sorted(all_data.keys()):
            quality = all_data[name].get("quality", "") or "NA"
            row_vals = [str(all_data[name].get(flag, False)) for flag in flags]
            comment = all_data[name].get("comment", "")
            writer.writerow([name, quality] + row_vals + [comment])


def write_txt(txt_path, flags, all_data):
    """Compact summary: TS_0001,Good,Mito,ER,Shell  (active flags only)."""
    with open(txt_path, "w") as f:
        for name in sorted(all_data.keys()):
            entry = all_data[name]
            quality = entry.get("quality", "") or "NA"
            active = [flag for flag in flags if entry.get(flag, False)]
            parts = [name, quality] + active
            f.write(",".join(parts) + "\n")


# ── MRC slice loading ─────────────────────────────────────────────────────────

def load_volume(path):
    with mrcfile.open(str(path), mode="r", permissive=True) as mrc:
        data = mrc.data.copy()
    return data


def compute_contrast_window(data, n_sample_slices=5, std_factor=3.0):
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


def slice_to_image(data, z_index, contrast_window, size=CANVAS_SIZE):
    z_index = max(0, min(z_index, data.shape[0] - 1))
    sl = np.asarray(data[z_index]).astype(np.float32)
    sl = bin_downsample(sl, size)
    lo, hi = contrast_window
    if hi <= lo:
        hi = lo + 1.0
    sl = np.clip((sl - lo) / (hi - lo), 0, 1)
    sl = (sl * 255).astype(np.uint8)
    img = Image.fromarray(sl, mode="L")
    return img.resize((size, size), Image.BILINEAR)


def prerender_all_slices(data, contrast_window, size=CANVAS_SIZE, progress_callback=None):
    n_slices = data.shape[0]
    images = []
    for z in range(n_slices):
        img = slice_to_image(data, z, contrast_window, size=size)
        images.append(ImageTk.PhotoImage(img))
        if progress_callback is not None and (z % 5 == 0 or z == n_slices - 1):
            progress_callback(z + 1, n_slices)
    return images


# ── Main annotator ────────────────────────────────────────────────────────────

class Annotator:
    def __init__(self, tomo_paths, flags, output_path, txt_path, suffix, initial_data):
        self.tomo_paths  = tomo_paths
        self.flags       = list(flags)
        self.output_path = output_path
        self.txt_path    = txt_path
        self.suffix      = suffix
        self.idx         = 0

        self.volume          = None
        self.rendered_images = []
        self.z               = 0
        self.photo_img       = None

        self.data            = initial_data
        self.current_state   = {flag: False for flag in self.flags}
        self.current_quality = ""
        self._list_updating  = False   # guard against recursive listbox events

        self._build_gui()

    # ── Helpers ──

    def _tomo_name(self, path):
        name = path.name
        if self.suffix and name.endswith(self.suffix):
            return name[: -len(self.suffix)]
        return name

    def _save_outputs(self):
        write_csv(self.output_path, self.flags, self.data)
        write_txt(self.txt_path,    self.flags, self.data)

    # ── GUI construction ──

    def _build_gui(self):
        self.root = tk.Tk()
        self.root.title("🍅 TOMATO - TOMogram Annotation TOol 🍅")
        self.root.resizable(False, False)
        self.root.protocol("WM_DELETE_WINDOW", self._on_quit)

        # Tomato look: darker red "skin" border around a lighter red
        # interior, instead of a colored title bar.
        BORDER_RED  = "#8B1E1E"
        INTERIOR_RED = "#F4E2DE"
        self.root.configure(bg=BORDER_RED)

        main = tk.Frame(self.root, bg=INTERIOR_RED)
        main.pack(padx=6, pady=6)

        # ── Column 0: image viewer ──
        left = tk.Frame(main, bg="#F4E2DE")
        left.grid(row=0, column=0, sticky="n")

        self.progress_var = tk.StringVar()
        tk.Label(left, textvariable=self.progress_var,
                 font=("Helvetica", 11), bg="#F4E2DE",
                 fg="#888780", anchor="w").pack(fill="x")

        self.name_var = tk.StringVar()
        tk.Label(left, textvariable=self.name_var,
                 font=("Helvetica", 15, "bold"),
                 bg="#F4E2DE", fg="#2C2C2A", anchor="w").pack(fill="x", pady=(0, 8))

        self.canvas = tk.Canvas(left, width=CANVAS_SIZE, height=CANVAS_SIZE,
                                bg="black", highlightthickness=0)
        self.canvas.pack()
        self.canvas.bind("<MouseWheel>", self._on_scroll)
        self.canvas.bind("<Button-4>",   self._on_scroll_linux_up)
        self.canvas.bind("<Button-5>",   self._on_scroll_linux_down)

        self.slice_var = tk.StringVar()
        tk.Label(left, textvariable=self.slice_var,
                 font=("Helvetica", 10), bg="#F4E2DE", fg="#888780").pack(pady=(4, 0))

        self.slider = tk.Scale(left, from_=0, to=1, orient="horizontal",
                               length=CANVAS_SIZE, showvalue=False,
                               bg="#F4E2DE", troughcolor="#F1EFE8",
                               highlightthickness=0, command=self._on_slider)
        self.slider.pack(pady=(4, 0))

        # ── Column 1: quality + annotation panel (stacked) ──
        right = tk.Frame(main, bg="#F4E2DE")
        right.grid(row=0, column=1, sticky="n", padx=(24, 0))

        tk.Label(right, text="Quality",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(0, 6))

        # Segmented control: one bordered frame holding all grade buttons
        # in a square 2x2 grid, separated by thin hairlines, like a
        # connected pill block rather than a single long row.
        GRADE_LABELS = {
            "FigureQuality": "Figure\nQuality",
            "Good":          "Good",
            "Ok":            "Ok",
            "Bad":           "Bad",
        }
        SEGMENT_BORDER = "#2C2C2A"
        segment_outer = tk.Frame(right, bg=SEGMENT_BORDER)
        segment_outer.pack(pady=(0, 16))
        segment_inner = tk.Frame(segment_outer, bg=SEGMENT_BORDER)
        segment_inner.pack(padx=1, pady=1)

        self.grade_buttons = {}
        for i, grade in enumerate(GRADES):
            row, col = divmod(i, 2)
            fg, bg_on = GRADE_COLORS[grade]
            b = tk.Button(segment_inner, text=GRADE_LABELS[grade],
                          font=("Helvetica", 12, "bold"),
                          fg=fg, bg="#F1EFE8",
                          activebackground=bg_on, activeforeground=fg,
                          relief="flat", bd=0, cursor="hand2",
                          width=8, height=2,
                          command=lambda g=grade: self._set_quality(g))
            b.grid(row=row, column=col, sticky="nsew",
                   padx=(0 if col == 0 else 1, 0),
                   pady=(0 if row == 0 else 1, 0))
            self.grade_buttons[grade] = (b, fg, bg_on)

        tk.Label(right, text="Annotation",
                 font=("Helvetica", 12, "bold"), bg="#F4E2DE",
                 fg="#2C2C2A", anchor="w").pack(fill="x", pady=(0, 6))

        self.buttons_frame = tk.Frame(right, bg="#F4E2DE")
        self.buttons_frame.pack(fill="x")

        self.flag_buttons = {}
        for flag in self.flags:
            self._add_flag_button(flag)

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

        # ── Column 2: tomogram list ──
        list_col = tk.Frame(main, bg="#F4E2DE")
        list_col.grid(row=0, column=2, sticky="ns", padx=(24, 0))

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
            height=32,
        )
        self.tomo_list.pack(side="left", fill="both", expand=True)
        scrollbar.config(command=self.tomo_list.yview)

        for path in self.tomo_paths:
            self.tomo_list.insert("end", self._tomo_name(path))

        self.tomo_list.bind("<<ListboxSelect>>", self._on_list_select)

        # ── Navigation ──
        nav_frame = tk.Frame(self.root, bg="#F4E2DE")
        nav_frame.pack(fill="x", padx=6, pady=(0, 6))

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

        self._load_current()
        self.root.mainloop()

    # ── Flag button creation ──

    def _add_flag_button(self, flag):
        i = len(self.flag_buttons)
        fg, bg_on = FLAG_COLORS[i % len(FLAG_COLORS)]
        b = tk.Button(self.buttons_frame, text=flag, width=18,
                      font=("Helvetica", 13, "bold"),
                      fg=fg, bg="#F1EFE8",
                      activebackground=bg_on,
                      relief="flat", bd=0, cursor="hand2",
                      anchor="w", padx=16, pady=10,
                      command=lambda fl=flag: self._toggle(fl))
        b.pack(fill="x", pady=4)
        self.flag_buttons[flag] = (b, fg, bg_on)

    def _add_flag_from_entry(self):
        name = self.add_entry.get().strip()
        if not name:
            return
        if name in self.flags:
            self.status_var.set(f"'{name}' already exists.")
            return
        self.flags.append(name)
        self.current_state[name] = False
        self._add_flag_button(name)
        self.add_entry.delete(0, "end")
        self.status_var.set(f"Added flag: {name}")

    # ── Toggle / quality logic ──

    def _refresh_buttons(self):
        for flag, (btn, fg, bg_on) in self.flag_buttons.items():
            if self.current_state.get(flag, False):
                btn.config(bg=bg_on, relief="solid", bd=1)
            else:
                btn.config(bg="#F1EFE8", relief="flat", bd=0)

    def _toggle(self, flag):
        self.current_state[flag] = not self.current_state.get(flag, False)
        self._refresh_buttons()

    def _set_quality(self, grade):
        self.current_quality = "" if self.current_quality == grade else grade
        self._refresh_quality_buttons()

    def _refresh_quality_buttons(self):
        for grade, (btn, fg, bg_on) in self.grade_buttons.items():
            if self.current_quality == grade:
                btn.config(bg=bg_on, fg=fg)
            else:
                btn.config(bg="#F1EFE8", fg=fg)

    # ── Slice navigation ──

    def _render_slice(self):
        if not self.rendered_images:
            return
        self.photo_img = self.rendered_images[self.z]
        self.canvas.delete("all")
        self.canvas.create_image(0, 0, anchor="nw", image=self.photo_img)
        self.slice_var.set(f"Slice {self.z + 1} / {len(self.rendered_images)}")
        self.slider.set(self.z)

    def _on_slider(self, value):
        self.z = int(float(value))
        self._render_slice()

    def _on_scroll(self, event):
        self._step_slice(1 if event.delta > 0 else -1)

    def _on_scroll_linux_up(self, event):
        self._step_slice(1)

    def _on_scroll_linux_down(self, event):
        self._step_slice(-1)

    def _step_slice(self, delta):
        if not self.rendered_images:
            return
        self.z = max(0, min(self.z + delta, len(self.rendered_images) - 1))
        self._render_slice()

    # ── Persisting current tomo's state ──

    def _commit_current(self):
        name = self._tomo_name(self.tomo_paths[self.idx])
        entry = dict(self.current_state)
        entry["quality"] = self.current_quality
        entry["comment"] = self.comment_box.get("1.0", "end").strip()
        self.data[name]  = entry
        self._save_outputs()

    # ── Tomogram list interaction ──

    def _on_list_select(self, *_):
        if self._list_updating:
            return
        sel = self.tomo_list.curselection()
        if not sel or sel[0] == self.idx:
            return
        self._commit_current()
        self.idx = sel[0]
        self._load_current()

    def _sync_list(self):
        self._list_updating = True
        self.tomo_list.selection_clear(0, "end")
        self.tomo_list.selection_set(self.idx)
        self.tomo_list.see(self.idx)
        self._list_updating = False

    # ── Loading a tomogram ──

    def _close_volume(self):
        self.volume = None
        self.rendered_images = []

    def _load_current(self):
        total = len(self.tomo_paths)
        self.progress_var.set(f"{self.idx + 1} / {total}")
        path  = self.tomo_paths[self.idx]
        name  = self._tomo_name(path)
        self.name_var.set(name)

        saved = self.data.get(name, {})
        self.current_state   = {flag: saved.get(flag, False) for flag in self.flags}
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
        self.status_var.set("Loading volume…")
        self.root.update()

        if not path.exists():
            self.status_var.set(f"File not found: {path}")
            self.canvas.delete("all")
            return

        try:
            self.volume = load_volume(path)
        except Exception as e:
            self.status_var.set(f"Failed to load: {e}")
            self.canvas.delete("all")
            return

        n_slices = self.volume.shape[0]
        self.slider.config(from_=0, to=n_slices - 1)
        self.z = n_slices // 2

        self.status_var.set("Computing contrast…")
        self.root.update()
        contrast_window = compute_contrast_window(self.volume)

        def report_progress(done, total):
            self.status_var.set(f"Processing tomogram… {done} / {total} slices")
            self.root.update()

        self.rendered_images = prerender_all_slices(
            self.volume, contrast_window, progress_callback=report_progress
        )
        self.status_var.set(f"Ready — {n_slices} slices, scroll or use slider to navigate")
        self._render_slice()

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


# ── Soft-link helper ──────────────────────────────────────────────────────────

def soft_link_object(csv_path, object_names, input_dir, suffix):
    """Read a TOMATO CSV and soft-link every tomogram that satisfies ALL filters.

    Each name in object_names is checked against:
      - the 'quality' column  (e.g. "Good", "Bad") if the name is a known grade
      - the flag column of that name  (value == "True") otherwise
    All conditions are combined with AND — a tomogram must match every filter.
    Links are placed in  ./<name1>_<name2>_..._tomograms/
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
    with open(csv_path, newline="") as f:
        reader = csv.DictReader(f)
        headers = list(reader.fieldnames or [])
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
                    if row.get(obj, "False").strip() != "True":
                        hit = False
                        break
            if hit:
                matching.append(name)

    label   = "_".join(object_names)
    out_dir = Path(f"{label}_tomograms")

    if not matching:
        print(f"No tomograms found matching all of: {', '.join(object_names)}")
        sys.exit(0)

    out_dir.mkdir(exist_ok=True)

    created = 0
    for name in sorted(matching):
        src  = input_dir / (name + suffix)
        link = out_dir   / (name + suffix)
        if not src.exists():
            print(f"  Warning: source not found — {src}")
            continue
        if link.is_symlink() or link.exists():
            link.unlink()
        os.symlink(src.resolve(), link)
        print(f"  {link.name}")
        created += 1

    print(f"\n{created} symlink(s) created in ./{out_dir}/")


# ── Entry point ───────────────────────────────────────────────────────────────

def _main_annotate():
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
    parser.add_argument("--output", default=DEFAULT_OUTPUT,
                        help=f"Output path for the annotation CSV; a matching .txt "
                             f"compact summary is written alongside with the same stem "
                             f"(default: {DEFAULT_OUTPUT} → annotations.csv + annotations.txt)")
    parser.add_argument("--suffix", default=DEFAULT_SUFFIX,
                        help=f"Tomogram filename suffix stripped from tomo_name "
                             f"(default: {DEFAULT_SUFFIX})")
    parser.add_argument("--previous_csv", default=None,
                        help="Path to a previously written TOMATO CSV. Flags and "
                             "annotations are loaded as a starting point. Output is "
                             "still written to --output.")
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
        csv_flags, initial_data = load_csv_with_flags(prev_path)
        for f in csv_flags:
            if f not in flags:
                flags.append(f)
        print(f"Previous CSV   : {prev_path} ({len(initial_data)} tomograms loaded)")

    output_path = Path(args.output)
    txt_path    = output_path.with_suffix(".txt")
    if output_path.exists():
        resume_data = load_csv(output_path, flags)
        if resume_data:
            initial_data.update(resume_data)
            print(f"Auto-resuming  : {output_path} ({len(resume_data)} tomograms)")

    tomo_paths = sorted(tomo_dir.glob(f"*{args.suffix}"))
    if not tomo_paths:
        print(f"No files matching *{args.suffix} found in {tomo_dir}")
        sys.exit(1)

    print(f"Tomo directory : {tomo_dir}")
    print(f"Flags          : {flags}")
    print(f"Output CSV     : {output_path}")
    print(f"Output TXT     : {txt_path}")
    print(f"Tomogram count : {len(tomo_paths)}")
    print()

    Annotator(tomo_paths, flags, output_path, txt_path, args.suffix, initial_data)


def _main_link():
    parser = argparse.ArgumentParser(
        prog="TOMATO.py link",
        description="Soft-link tomograms that match an object or quality grade from a TOMATO CSV.",
    )
    parser.add_argument("--previous_csv", required=True,
                        help="TOMATO CSV to read annotations from (REQUIRED)")
    parser.add_argument("--object", required=True, nargs="+",
                        help="One or more object names / quality grades to filter by (AND logic). "
                             "E.g. --object Good Mito  → tomograms rated Good that also contain Mito. "
                             "Output directory is named after all values joined with '_'.")
    parser.add_argument("--input", required=True,
                        help="Directory containing the tomogram .mrc files (REQUIRED)")
    parser.add_argument("--suffix", default=DEFAULT_SUFFIX,
                        help=f"Tomogram filename suffix (default: {DEFAULT_SUFFIX})")
    args = parser.parse_args(sys.argv[2:])

    soft_link_object(
        csv_path     = args.previous_csv,
        object_names = args.object,
        input_dir    = args.input,
        suffix       = args.suffix,
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
