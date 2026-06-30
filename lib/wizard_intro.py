#!/usr/bin/env python3
"""CLIPittyCLIP wizard-mode intro animation.

Plays a terminal animation: a paperclip appears, a wizard hat sweeps in
from the left and absorbs it, then the pipeline title appears to the right.

Requires a terminal of at least 96 cols × 28 rows.
Falls back to a one-liner if the terminal is too small or non-interactive.

Usage:
    python3 wizard_intro.py          # run directly
    from wizard_intro import intro   # import and call
"""

import os
import sys
import signal
import time

# ── sprites ──────────────────────────────────────────────────────────────────

HAT = [
    "                         +++",
    "                        ++++",
    "                      .##++#",
    "                     .++++++",
    "                    .+#++#++",
    "                    +#+++#++",
    "                   -++++++++",
    "                  .+#++++#++",
    "                  +#+++++++#",
    "                 .+++#++####+",
    "                 ++++##+#+#++",
    "                 +#+++###+++#",
    "                ++#+++#+##+##+.",
    "                ++###+#+#++#++#",
    "               ####++++++###++#.",
    "               +#++#+++++#+#+##+.",
    "               ##+##+##+++#++++##",
    "              +#+###+#+#+++#++####.",
    "       .++++++++++#######+#+-+-##+-++++++++++",
    "     -+++#+++-#+--+####+###---++---##+#+++##+++",
    "  .++++++++++++#+##++##+##+-+#+##++-----#+++++++",
    ".+++++#++++++##+-----#++##+##+##+#++++#+++#++",
    ".+++++++-#++#++++#####+++++##+##+++++##++++#",
    "  .++++++#+#+++++++###+++#++###++++#++#+#",
    "      ++++#++++++++#+++#++++++++#+#++#+",
    "         .  +#+++++###+##+#++++## .",
]

CLIP = [
    "         .-.",
    "       -.  .+",
    "      -.    .+",
    "      -.      #.",
    "       .+.-.   ...",
    "            ...  .-. .##.",
    "              .#.   #-####",
    "     .+.        .#-.###-.###.#-",
    "   ..   --        -#.#####+##+.",
    " .+        ..+##++-##########.",
    "  -...-+++-.    +-#####+####+",
    "   ..       ..+#######+## ##",
    "              .##########++",
    "                ..#######+",
    "                      ###",
]

# ── layout constants ─────────────────────────────────────────────────────────

MIN_W, MIN_H = 96, 28
GRID_W, GRID_H = 96, 28
BUF = 3                       # edge-wipe buffer (chars between hat edge & clip)

HAT_W = max(len(l) for l in HAT)
HAT_H = len(HAT)
CLIP_W = max(len(l) for l in CLIP)
CLIP_H = len(CLIP)

HAT_Y = 1
HAT_FINAL_X = 1
CLIP_X = 12                   # positioned so hat fully covers clip at rest
CLIP_Y = HAT_Y + HAT_H - CLIP_H  # bottom-aligned with hat

TEXT_X = HAT_FINAL_X + HAT_W + 4
TITLE_Y = HAT_Y + (HAT_H - 8) // 2
BANNER_Y = TITLE_Y + 5

TAGLINE = "End-to-End CLIP-seq Analysis Pipeline"

# ── grid compositing ────────────────────────────────────────────────────────

def _grid():
    return [[" "] * GRID_W for _ in range(GRID_H)]


def _place(g, rows, x, y):
    for j, line in enumerate(rows):
        ry = y + j
        if ry < 0 or ry >= GRID_H:
            continue
        for i, ch in enumerate(line):
            if ch == " ":
                continue
            rx = x + i
            if 0 <= rx < GRID_W:
                g[ry][rx] = ch


def _place_edge_wiped(g, rows, x, y, row_min_x):
    for j, line in enumerate(rows):
        ry = y + j
        if ry < 0 or ry >= GRID_H:
            continue
        mn = row_min_x[j] if j < len(row_min_x) and row_min_x[j] >= 0 else -1
        for i, ch in enumerate(line):
            if ch == " ":
                continue
            rx = x + i
            if rx < 0 or rx >= GRID_W:
                continue
            if mn >= 0 and rx < mn:
                continue
            g[ry][rx] = ch


def _hat_edge(hat_x):
    """Per-row right-edge wipe position for the hat at hat_x."""
    edges = []
    for line in HAT:
        right = -1
        for i in range(len(line) - 1, -1, -1):
            if line[i] != " ":
                right = i
                break
        edges.append(-1 if right < 0 else hat_x + right + 1 + BUF)
    return edges


def _place_box(g, label, y, x):
    bar = "+" + "-" * (len(label) + 2) + "+"
    mid = "| " + label + " |"
    _place(g, [bar], x, y)
    _place(g, [mid], x, y + 1)
    _place(g, [bar], x, y + 2)


def _render(g):
    return "\n".join("".join(row).rstrip() for row in g)


# ── frame builder ────────────────────────────────────────────────────────────

def _compose(hat_x=None, show_clip=False, show_title=False):
    g = _grid()

    if show_clip:
        if hat_x is not None:
            all_edges = _hat_edge(hat_x)
            y_off = CLIP_Y - HAT_Y
            clip_edges = []
            for j in range(CLIP_H):
                hr = y_off + j
                clip_edges.append(all_edges[hr] if 0 <= hr < len(all_edges) else -1)
            _place_edge_wiped(g, CLIP, CLIP_X, CLIP_Y, clip_edges)
        else:
            _place(g, CLIP, CLIP_X, CLIP_Y)

    if hat_x is not None:
        _place(g, HAT, hat_x, HAT_Y)

    if show_title:
        _place_box(g, "CLIPittyCLIP", TITLE_Y, TEXT_X)
        _place_box(g, TAGLINE, BANNER_Y, TEXT_X)

    return _render(g)


# ── animation timeline ──────────────────────────────────────────────────────

def _build_frames():
    frames = []

    # phase 1: hold paperclip
    clip_frame = _compose(show_clip=True)
    for _ in range(8):
        frames.append(clip_frame)

    # phase 2: hat slides in, edge-wiping clip
    for hx in range(-HAT_W, HAT_FINAL_X + 1, 3):
        frames.append(_compose(hat_x=hx, show_clip=True))
    # snap to final
    frames.append(_compose(hat_x=HAT_FINAL_X, show_clip=True))

    # phase 3: hold hat (clip absorbed)
    hat_frame = _compose(hat_x=HAT_FINAL_X)
    for _ in range(5):
        frames.append(hat_frame)

    # phase 4: title appears + hold
    final_frame = _compose(hat_x=HAT_FINAL_X, show_title=True)
    for _ in range(12):
        frames.append(final_frame)

    return frames


# ── terminal helpers ─────────────────────────────────────────────────────────

def _hide_cursor():
    sys.stdout.write("\033[?25l")
    sys.stdout.flush()


def _show_cursor():
    sys.stdout.write("\033[?25h")
    sys.stdout.flush()


def _clear_and_home():
    sys.stdout.write("\033[2J\033[H")
    sys.stdout.flush()


def _move_home():
    sys.stdout.write("\033[H")
    sys.stdout.flush()


# ── public API ───────────────────────────────────────────────────────────────

def intro():
    """Play the wizard intro animation.

    Gates on terminal size (96×28 minimum).
    Falls back to a text one-liner if terminal is too small or not a TTY.
    """
    if not sys.stdout.isatty():
        print(f"\n  CLIPittyCLIP — {TAGLINE}\n")
        return

    try:
        cols, rows = os.get_terminal_size()
    except OSError:
        print(f"\n  CLIPittyCLIP — {TAGLINE}\n")
        return

    if cols < MIN_W or rows < MIN_H:
        print(f"\n  CLIPittyCLIP — {TAGLINE}\n")
        return

    # center the animation in the terminal
    pad_left = max(0, (cols - GRID_W) // 2)
    pad_top = max(0, (rows - GRID_H) // 2)
    left_margin = " " * pad_left

    frames = _build_frames()

    # clean exit on Ctrl+C
    interrupted = False

    def _on_sigint(sig, frame):
        nonlocal interrupted
        interrupted = True

    prev_handler = signal.signal(signal.SIGINT, _on_sigint)

    try:
        _hide_cursor()
        _clear_and_home()

        for frame_str in frames:
            if interrupted:
                break
            _move_home()
            # top padding
            sys.stdout.write("\n" * pad_top)
            # frame lines with left padding
            for line in frame_str.split("\n"):
                sys.stdout.write(left_margin + line + "\n")
            sys.stdout.flush()
            time.sleep(0.12)

    finally:
        _show_cursor()
        _clear_and_home()
        signal.signal(signal.SIGINT, prev_handler)


# ── entry point ──────────────────────────────────────────────────────────────

if __name__ == "__main__":
    intro()
