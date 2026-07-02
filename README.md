# TOMATO — TOMogram Annotation TOol

TOMATO (TOMogram Annotation TOol) is a lightweight Python GUI for browsing cryo-ET tomograms and annotating them with organelle flags and quality grades. It loads .mrc files directly, autosaves progress, and can soft-link tomograms matching any annotation for downstream processing.

It has two subcommands:
- `annotate` — open the GUI and annotate a folder of tomograms
- `link` — soft-link tomograms that match a saved annotation (e.g. all "Good" tomograms containing Mito)

---

## License and Copyright

TOMATO is open-source software licensed under the GNU GPLv2 (or later). For full copyright, author attributions, and license details, please refer to the header at the top of the main script file.

---

## 1. Dependencies

TOMATO requires Python 3 with the following packages:

```bash
pip install mrcfile numpy Pillow
```

Tkinter is included with most Python installations. If it is missing:
- **Debian/Ubuntu:** `sudo apt install python3-tk`
- **Fedora:** `sudo dnf install python3-tkinter`
- **macOS (Homebrew):** `brew install python-tk`
- **conda:** `conda install tk`

**EMBL users:** load the following modules  instead of using pip
(module versions current as of 2026 — check with `ml avail` if these fail):
```bash
ml mrcfile
ml SciPy-bundle/2025.07-gfbf-2025b
ml Tkinter/3.13.5-GCCcore-14.3.0
ml Pillow/11.3.0-GCCcore-14.3.0
```

---

## 2. Annotating tomograms

### Basic usage

```bash
python TOMATO.py annotate --input /path/to/reconstruction
```

This opens the GUI, listing every .mrc file in `--input` that matches the default suffix (`.mrc`), and writes results to `annotations.csv` / `annotations.txt` in the current directory.

### Full argument list

<table>
<tr>
  <th width="160">Argument</th>
  <th>Required</th>
  <th>Default</th>
  <th>Description</th>
</tr>
<tr>
  <td><code>--input</code></td>
  <td>Yes</td>
  <td>—</td>
  <td>Folder containing the tomogram <code>.mrc</code> files to annotate.</td>
</tr>
<tr>
  <td><code>--flags</code></td>
  <td>No</td>
  <td><em>(none)</em></td>
  <td>Extra button labels, added on top of the always-present base set Mito, ER, Golgi. Space-separated, e.g. <code>--flags Lysosome Nucleus Vesicle</code>. You can also add flags live from inside the GUI ("Add flag" box) without restarting.</td>
</tr>
<tr>
  <td><code>--output</code></td>
  <td>No</td>
  <td><code>annotations.csv</code></td>
  <td>Where the full annotation table is written. A compact <code>.txt</code> summary is written alongside it automatically, using the same filename stem (e.g. <code>--output myrun.csv</code> → <code>myrun.csv</code> + <code>myrun.txt</code>).</td>
</tr>
<tr>
  <td><code>--suffix</code></td>
  <td>No</td>
  <td><code>.mrc</code></td>
  <td>The filename suffix used to select which tomograms to annotate, and stripped when generating the <code>tomo_name</code> key in the CSV. Set this to match your actual filenames (e.g. <code>_10.00Apx.mrc</code>, or <code>.mrc</code> if there's no pixel-size suffix).</td>
</tr>
<tr>
  <td><code>--previous_csv</code></td>
  <td>No</td>
  <td><em>(none)</em></td>
  <td>Path to a previously written TOMATO CSV to seed annotations from. Useful for continuing someone else's annotation run or branching off an existing CSV without overwriting it (new output still goes to <code>--output</code>). Any extra flag columns found in that CSV are automatically added to the flag list.</td>
</tr>
</table>

### Auto-resume

If the file at `--output` already exists when you start, TOMATO automatically loads it and resumes — so you can simply re-run the same command to continue an interrupted session. This is separate from `--previous_csv`, which is for explicitly seeding from a different file.

### Example: annotating denoised oocyte tomograms with extra flags

```bash
python TOMATO.py annotate \
    --input /path/to/tomograms \
    --flags Lysosome Vesicle Shell \
    --output oocyte_annotations.csv \
    --suffix _19.12Apx.mrc
```

### Using the GUI

- **Image viewer (left):** shows the current tomogram's mid-slice by default. Scroll the mouse wheel over the image, or drag the slider underneath it, to move through Z-slices. All slices are pre-rendered when a tomogram loads, so scrolling is smooth once the "Ready" status appears.
- **Annotation (top right):** click a flag button (Mito, ER, Golgi, or any custom flag) to toggle it on/off for the current tomogram. Active flags are highlighted in their assigned color.
- **Add flag:** type a new flag name and press Enter or click "+" to add a button for it on the fly. It applies to all tomograms in the session from that point on.
- **Comment:** free-text notes for the current tomogram.
- **Quality:** a single segmented control — click one of `FigureQuality`, `Good`, `Ok`, `Bad` to grade the tomogram. Clicking the same grade again clears it.
- **Tomograms (right list):** click any entry to jump directly to that tomogram. The current tomogram's annotations are saved automatically before switching.
- **Previous / Next:** step through tomograms in order. The button reads "Finish ✓" on the last tomogram. Annotations are saved on every navigation action, so there's no separate "save" button.
- **Quit:** saves the current tomogram's state and closes the GUI.

Progress is written to disk after every tomogram you move away from, so it's safe to close the window at any point — restarting with the same `--output` path resumes where you left off.

---

## 3. Soft-linking tomograms by annotation (`link`)

Once you've annotated a batch, use `link` to create a folder of symlinks pointing at tomograms matching one or more conditions (combined with AND logic).

```bash
python TOMATO.py link \
    --previous_csv oocyte_annotations.csv \
    --object Good Mito \
    --input /path/to/tomograms
```

This creates a folder `Good_Mito_tomograms/` containing symlinks to every tomogram graded `Good` and flagged `Mito`.

### Argument list

<table>
<tr>
  <th width="160">Argument</th>
  <th>Required</th>
  <th>Default</th>
  <th>Description</th>
</tr>
<tr>
  <td><code>--previous_csv</code></td>
  <td>Yes</td>
  <td>—</td>
  <td>The TOMATO CSV to read annotations from.</td>
</tr>
<tr>
  <td><code>--object</code></td>
  <td>Yes</td>
  <td>—</td>
  <td>One or more values to filter by, AND-combined. Each value is matched against either the <code>quality</code> column (if it's one of <code>FigureQuality</code>, <code>Good</code>, <code>Ok</code>, <code>Bad</code>) or a flag column. Output folder is named after all values joined with <code>_</code>.</td>
</tr>
<tr>
  <td><code>--input</code></td>
  <td>Yes</td>
  <td>—</td>
  <td>Folder containing the original tomogram <code>.mrc</code> files (the symlink targets).</td>
</tr>
<tr>
  <td><code>--suffix</code></td>
  <td>No</td>
  <td><code>.mrc</code></td>
  <td>Must match the suffix used when the CSV was generated, so filenames can be reconstructed from <code>tomo_name</code>.</td>
</tr>
</table>

### More examples

```bash
# All tomograms flagged Mito, regardless of quality
python TOMATO.py link --previous_csv oocyte_annotations.csv --object Mito --input /path/to/reconstruction

# All tomograms rated FigureQuality (for making figures)
python TOMATO.py link --previous_csv oocyte_annotations.csv --object FigureQuality --input /path/to/tomograms

# Good tomograms with both ER and Golgi
python TOMATO.py link --previous_csv oocyte_annotations.csv --object Good ER Golgi --input /path/to/tomograms
```

If a value in `--object` doesn't match any known flag or quality grade, TOMATO prints the available options and exits without creating anything.

---

## 4. Output file formats

**CSV** (`<output>.csv`): one row per tomogram, with columns `tomo_name`, `quality`, one column per flag (`True`/`False`), and `comment`. This is the canonical record and what `--previous_csv` / `link` read from.

**TXT** (`<output>.txt`): a compact summary, one line per tomogram, listing only the active flags:

```
TS_0001,Good,Mito,ER,Shell
TS_0002,Bad
TS_0003,FigureQuality,Mito,Golgi
```

Useful for quick `grep`/`awk` filtering without parsing the full CSV.
