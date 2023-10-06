import csv
import sys
from typing import Any

sys.stderr = open(snakemake.log[0], "w")


class UnknownIvarHeaderFormatError(Exception):
    """Raised when an unknown ivar header format is encountered."""


class InvalidValueError(Exception):
    """Raised when an unknown ivar header format is encountered."""


class MixedPositionDeterminator:
    alt_depth: int
    min_alt_freq: float
    max_alt_freq: float
    total_depth: int

    def __init__(self, *, alt_depth: int, min_alt_freq: float, max_alt_freq: float, total_depth: int):
        if alt_depth < 0:
            raise InvalidValueError(f"{alt_depth=} value should be greater or equal than 0")

        if total_depth < 0:
            raise InvalidValueError(f"{total_depth=} value should be greater or equal than 0")

        if min_alt_freq < 0 and min_alt_freq > 1:
            raise InvalidValueError(f"{min_alt_freq=} value should be between 0 and 1")

        if max_alt_freq < 0 and max_alt_freq > 1:
            raise InvalidValueError(f"{max_alt_freq=} value should be between 0 and 1")

        if min_alt_freq > max_alt_freq:
            raise InvalidValueError(f"{min_alt_freq=} value should be greater than {max_alt_freq=}")

        self.alt_depth = alt_depth
        self.min_alt_freq = min_alt_freq
        self.max_alt_freq = max_alt_freq
        self.total_depth = total_depth

    def _is_row_a_mixed_position(self, row: dict[str | Any, str | Any]):
        return (
            int(row["ALT_DP"]) >= self.alt_depth
            and float(row["ALT_FREQ"]) >= self.min_alt_freq
            and float(row["ALT_FREQ"]) < self.max_alt_freq
            and int(row["TOTAL_DP"]) >= self.total_depth
        )

    def process_rows(self, rows: list[dict[str | Any, str | Any]]):
        return [row for row in rows if self._is_row_a_mixed_position(row)]


def load_ivar_variants(ivar_tsv: str):
    with open(ivar_tsv, "r") as ivar_file:
        ivar_reader = csv.DictReader(ivar_file, delimiter="\t")
        if not ivar_reader.fieldnames:
            raise ValueError("No header in ivar .tsv file")
        ivar_header: list[str] = list(ivar_reader.fieldnames)
        ivar_rows = list(ivar_reader)

        if any(
            required_column not in ivar_header
            for required_column in ["REGION", "POS", "ALT_DP", "ALT_FREQ", "TOTAL_DP"]
        ):
            raise UnknownIvarHeaderFormatError("Found unknown header: %s" % ivar_header)
    return ivar_rows, ivar_header


def compute_mixed_positions(
    *,
    ivar_tsv: str,
    out_mixed_positions_tsv: str,
    out_count_file: str,
    mixed_position_determinator: MixedPositionDeterminator,
):
    ivar_rows, header = load_ivar_variants(ivar_tsv)
    mixed_positions = mixed_position_determinator.process_rows(ivar_rows)
    count = len(set([(i["REGION"], i["POS"]) for i in mixed_positions]))

    with open(out_mixed_positions_tsv, "w") as summary_file:
        summary_writer = csv.DictWriter(summary_file, delimiter="\t", fieldnames=header)
        summary_writer.writeheader()
        summary_writer.writerows(mixed_positions)

    with open(out_count_file, "w") as count_file:
        count_file.write(f"{count}")


if __name__ == "__main__":
    mixed_position_determinator = MixedPositionDeterminator(
        alt_depth=int(snakemake.params["alt_depth"]),
        min_alt_freq=float(snakemake.params["min_alt_freq"]),
        max_alt_freq=float(snakemake.params["max_alt_freq"]),
        total_depth=int(snakemake.params["total_depth"]),
    )
    compute_mixed_positions(
        ivar_tsv=snakemake.input[0],
        out_mixed_positions_tsv=snakemake.output.mixed_positions,
        out_count_file=snakemake.output.readcount,
        mixed_position_determinator=mixed_position_determinator,
    )
