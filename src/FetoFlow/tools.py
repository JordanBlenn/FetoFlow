import subprocess
import sys
import pathlib


def _run(cmd: list[str]) -> None:
    print(" ".join(cmd))
    if (code := subprocess.call(cmd)) != 0:
        sys.exit(code)


def lint() -> None:  # this becomes the `lint` console command
    root = pathlib.Path(__file__).resolve().parents[2]  # project root
    _run(["ruff", "format", str(root / "src"), str(root / "tests")])
    _run(["ruff", "check", str(root / "src"), str(root / "tests"), "--fix"])
    _run(["pylint", "--rcfile", str(root / ".pylintrc"), "src", "tests"])
