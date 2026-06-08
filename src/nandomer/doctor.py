import platform
import shutil
import subprocess
import sys


def _check_python():
    version = platform.python_version()
    return True, f"Python {version}"


def _check_tool(tool_name):
    tool_path = shutil.which(tool_name)
    if not tool_path:
        return False, f"{tool_name}: not found in PATH"

    try:
        result = subprocess.run(
            [tool_name, "--version"],
            check=True,
            capture_output=True,
            text=True,
        )
    except subprocess.CalledProcessError as exc:
        details = (exc.stderr or exc.stdout or "").strip()
        return False, f"{tool_name}: found at {tool_path}, but '--version' failed ({details})"

    first_line = (result.stdout or result.stderr or "").splitlines()
    version_line = first_line[0].strip() if first_line else "version unknown"
    return True, f"{tool_name}: {version_line} ({tool_path})"


def run_doctor():
    checks = []

    py_ok, py_info = _check_python()
    checks.append(("python", py_ok, py_info))

    sam_ok, sam_info = _check_tool("samtools")
    checks.append(("samtools", sam_ok, sam_info))

    ok = all(item[1] for item in checks)

    lines = ["nandomer doctor report"]
    for name, passed, details in checks:
        status = "OK" if passed else "MISSING"
        lines.append(f"- {name}: {status} - {details}")

    if not sam_ok:
        lines.append("- hint: install samtools, e.g. 'brew install samtools' or 'conda install -c bioconda samtools'.")

    lines.append(f"Overall: {'PASS' if ok else 'FAIL'}")
    return "\n".join(lines), ok
