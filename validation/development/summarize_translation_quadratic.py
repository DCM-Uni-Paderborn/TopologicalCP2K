"""Verify the archive and re-extract numerical bounds without rerunning CP2K."""

import argparse
import hashlib
import json
import re
import tarfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("archive")
    args = parser.parse_args()
    with tarfile.open(args.archive, "r:gz") as archive:
        members = archive.getmembers()
        names = [member.name for member in members]
        assert len(set(names)) == len(names), "Duplicate archive member"
        manifest = json.load(archive.extractfile("manifest.json"))
        assert set(names) == set(manifest["files"]) | {"manifest.json"}
        logs = {}
        for name, record in manifest["files"].items():
            data = archive.extractfile(name).read()
            assert len(data) == record["bytes"], name
            assert hashlib.sha256(data).hexdigest() == record["sha256"], name
            if name.endswith((".log", ".out")):
                logs[name] = data.decode()

    number = r"[-+0-9.Ee]+"
    gram = []
    folded = []
    distributed = []
    regressions = {}
    for name, content in logs.items():
        for line in content.splitlines():
            if line.startswith("Translation Gram/"):
                values = list(map(float, line.split(":", 1)[1].split()))
                assert len(values) == 4
                gram.append({"log": name, "values": values})
            match = re.fullmatch(
                rf"Trimerized full/Bloch spectra: (\d+) eigenvalues; maximum error\s+({number})", line)
            if match:
                assert int(match[1]) == 1170
                folded.append({"log": name, "eigenvalues_per_scan": int(match[1]),
                               "maximum_error": float(match[2])})
            match = re.fullmatch(r"Quadratic DBCSR metric mode ([01]) forward/adjoint/solve/spectrum/Ritz:\s+(.+)", line)
            if match:
                values = list(map(float, match[2].split()))
                assert len(values) == 5
                distributed.append({"log": name, "metric_mode": int(match[1]), "values": values})
        if name.endswith("translation-quadratic-final-regtests.log"):
            match = re.search(r"Summary: correct: (\d+) / (\d+);", content)
            assert match and match[1] == match[2] and "Status: OK" in content, name
            regressions[name] = int(match[1])
    assert gram and folded and distributed and len(regressions) == 2
    assert regressions["build-serial/translation-quadratic-final-regtests.log"] == 53
    assert regressions["build-mpi/translation-quadratic-final-regtests.log"] == 88
    bounds = {
        "gram_metric_complex_mean_split_maxima": [max(r["values"][j] for r in gram) for j in range(4)],
        "folded_spectrum_maximum_error": max(r["maximum_error"] for r in folded),
        "dbcsr_forward_adjoint_metric_spectrum_ritz_maxima": [
            max(r["values"][j] for r in distributed) for j in range(5)],
    }
    assert bounds["gram_metric_complex_mean_split_maxima"][0] < 4e-15
    assert bounds["gram_metric_complex_mean_split_maxima"][1] < 6e-15
    assert bounds["gram_metric_complex_mean_split_maxima"][2] < 3e-14
    assert bounds["gram_metric_complex_mean_split_maxima"][3] > 0.08
    assert bounds["folded_spectrum_maximum_error"] < 1.61e-14
    assert bounds["dbcsr_forward_adjoint_metric_spectrum_ritz_maxima"][3] < 9e-15
    assert bounds["dbcsr_forward_adjoint_metric_spectrum_ritz_maxima"][4] < 1e-10
    print(json.dumps({"source_commit": manifest["source_commit"],
                      "verified_files": len(manifest["files"]), "bounds": bounds,
                      "regressions": regressions, "native_records": gram,
                      "bloch_records": folded, "distributed_records": distributed,
                      "units": "model units with t=1; eigenvalues and Ritz residuals in t^2",
                      "scope": "numerical translation kernel, not Gaussian DFT unfolding"}, indent=2))


if __name__ == "__main__":
    main()
