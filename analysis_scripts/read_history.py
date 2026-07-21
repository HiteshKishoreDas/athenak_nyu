"""Read an Athena++ history file into a dictionary.

Athena++ history files usually store the column names in a commented header
line like ``#  [1]=time [2]=dt ...`` followed by numeric rows. This parser
extracts those names when present, while still falling back to the first
non-comment row for simpler whitespace-separated tables.
"""

import re
import numpy as np
from pathlib import Path
from typing import Dict, List, Union


Number = Union[int, float]
_HEADER_TOKEN_RE = re.compile(r"\[\d+\]=([^\s]+)")


def _convert(value: str) -> Union[str, Number]:
	"""Convert a string token to int/float when possible."""
	try:
		if any(ch in value for ch in (".", "e", "E")):
			return float(value)
		return int(value)
	except ValueError:
		try:
			return float(value)
		except ValueError:
			return value


def read_history_file(
	filepath: Union[str, Path] = Path(__file__).with_name("Turb.hydro.hst"),
) -> Dict[str, List[Union[str, Number]]]:
	"""Read a Turb.hydro.hst file into a dictionary.

	Parameters
	----------
	filepath
		Path to the history file.

	Returns
	-------
	dict
		Dictionary mapping each column name to a list of parsed values.
	"""
	path = Path(filepath)
	if not path.is_file():
		raise FileNotFoundError(f"History file not found: {path}")

	header = None
	rows = []

	with path.open("r", encoding="utf-8") as f:
		for line in f:
			line = line.strip()
			if not line:
				continue
			if line.startswith("#"):
				if header is None:
					comment_header = _HEADER_TOKEN_RE.findall(line)
					if comment_header:
						header = comment_header
				continue
			tokens = line.split()
			if header is None:
				header = tokens
				continue
			if len(tokens) != len(header):
				raise ValueError(
					f"Row has {len(tokens)} columns, expected {len(header)}: {line}"
				)
			rows.append([_convert(token) for token in tokens])

	if header is None:
		raise ValueError(f"No data header found in file: {path}")
	if not rows:
		raise ValueError(f"No data rows found in file: {path}")

	history = {name: [] for name in header}
	for row in rows:
		for name, value in zip(header, row):
			history[name].append(value)

	return history


history = read_history_file()
print(history.keys())

usr_hist = read_history_file(filepath = "Turb.user.hst") 
print(usr_hist.keys())



skip = 0
vol = 1.0

time = np.array(history["time"][skip:])
mass = np.array(history["mass"][skip:])
KE1 = np.array(history["1-KE"][skip:])
KE2 = np.array(history["2-KE"][skip:])
KE3 = np.array(history["3-KE"][skip:])
Etot = np.array(history["tot-E"][skip:])

time_usr = np.array(usr_hist["time"][skip:])
Tavg = np.array(usr_hist["Tsumvol"][skip:])/vol

KE = KE1 + KE2 + KE3
TE = Etot - KE

import matplotlib.pyplot as plt

# plt.plot(time)

plt.plot(time, TE, label="TE")
plt.plot(time, KE, label="KE")
plt.plot(time, Etot, label="Etot")

plt.legend()
plt.xlabel("time")
plt.ylabel("E")
plt.show()

plt.plot(time_usr, Tavg, label="Tavg")
plt.legend()
plt.xlabel("time")
plt.ylabel("Tavg")
plt.show()
