import shlex
import subprocess
import sys
from pathlib import Path
from typing import List, Tuple

from PIL import Image


FRAME_DIR = Path("bin_1280_6400/save/temp")
FRAME_PATTERN = "temp_%05d.png"
OUTPUT_VIDEO = Path("bin_1280_6400/save/temp/temperature.mp4")
FPS = 8
START_NUMBER = None  # auto-detect from frames
FFMPEG_PRESET = "ultrafast"
FFMPEG_CRF = 10  # higher CRF lowers quality but speeds encoding
MAX_WIDTH = 4096  # aggressive downscale to reduce work
BLACK_BACKGROUND = True


def ensure_frames_exist(frame_dir: Path, frame_pattern: str) -> List[Path]:
    if not frame_dir.exists():
        sys.exit(f"Frame directory {frame_dir} does not exist")

    matches = sorted(frame_dir.glob(frame_pattern.replace("%05d", "*")))
    if not matches:
        sys.exit(
            f"No frames matching pattern {frame_pattern} found in {frame_dir}. "
            "Generate frames before running this script."
        )
    return matches


def get_frame_size(frame: Path) -> Tuple[int, int]:
    with Image.open(frame) as img:
        return img.width, img.height


def build_filter_complex(frame_size: Tuple[int, int]) -> Tuple[str, str]:
    width, height = frame_size
    filter_nodes = []
    current_label = "[0:v]"

    if BLACK_BACKGROUND:
        filter_nodes.append("[0:v]format=rgba[fg]")
        filter_nodes.append(f"color=color=black:size={width}x{height}[bg]")
        filter_nodes.append("[bg][fg]overlay=shortest=1[comp]")
        current_label = "[comp]"
    else:
        filter_nodes.append("[0:v]format=rgba[comp]")
        current_label = "[comp]"

    if MAX_WIDTH:
        mw = MAX_WIDTH
        filter_nodes.append(
            f"{current_label}scale='if(gt(iw,{mw}),{mw},iw)':"
            f"'if(gt(iw,{mw}),-2,ih)'[scaled]"
        )
        current_label = "[scaled]"

    filter_nodes.append(f"{current_label}format=yuv420p[outv]")
    return ";".join(filter_nodes), "[outv]"


def run_ffmpeg(
    frame_dir: Path,
    frame_pattern: str,
    output: Path,
    fps: int,
    preset: str,
    crf: int,
    frame_size: Tuple[int, int],
    start_number: int,
) -> None:
    pattern_path = frame_dir / frame_pattern
    filter_complex, map_label = build_filter_complex(frame_size)

    ffmpeg_cmd = (
        "module load ffmpeg && "
        f"ffmpeg -hide_banner -loglevel warning -y "
        f"-framerate {fps} -start_number {start_number} "
        f"-i {shlex.quote(str(pattern_path))} "
        f"-filter_complex {shlex.quote(filter_complex)} "
        f"-map {shlex.quote(map_label)} "
        "-c:v libx264 "
        f"-preset {shlex.quote(preset)} "
        f"-crf {crf} "
        f"{shlex.quote(str(output))}"
    )

    print(f"Running: {ffmpeg_cmd}")
    subprocess.run(["bash", "-lc", ffmpeg_cmd], check=True)
    print(f"Video written to {output}")


def main() -> None:
    frames = ensure_frames_exist(FRAME_DIR, FRAME_PATTERN)
    frame_size = get_frame_size(frames[0])
    if START_NUMBER is None:
        try:
            start_number = min(
                int(frame.stem.split("_")[-1]) for frame in frames if frame.stem
            )
        except ValueError:
            start_number = 0
    else:
        start_number = START_NUMBER
    run_ffmpeg(
        frame_dir=FRAME_DIR,
        frame_pattern=FRAME_PATTERN,
        output=OUTPUT_VIDEO,
        fps=FPS,
        start_number=start_number,
        preset=FFMPEG_PRESET,
        crf=FFMPEG_CRF,
        frame_size=frame_size,
    )


if __name__ == "__main__":
    main()
