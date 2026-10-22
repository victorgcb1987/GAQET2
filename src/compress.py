import os
import shutil
import tarfile

from pathlib import Path


ARCHIVE_NAME = "GAQET_intermediate_results.tar.gz"


def compress_intermediate(basedir, keep, protect=()):
    """Pack every intermediate result of a GAQET run into one tar.gz and delete the originals.

    basedir: GAQET output directory.
    keep: files that stay where they are, uncompressed (stats table, log, DeTEnGA summary).
    protect: files that are never archived nor deleted (e.g. the user's input files if they
             live inside basedir).
    Originals are only removed once the archive has been written and checked. If anything
    fails the run is left untouched and the error is raised.
    Returns (archive path, number of files, size before, size after) in bytes.
    """
    basedir = Path(basedir).resolve()
    skip = {Path(path).resolve() for path in list(keep) + list(protect)}
    archive = basedir / ARCHIVE_NAME
    partial = basedir / (ARCHIVE_NAME + ".part")
    skip.update({archive, partial})

    files, directories = [], []
    for root, dirs, names in os.walk(basedir):
        root = Path(root)
        for name in dirs:
            if (root / name).is_symlink():
                files.append(root / name)  # link to a directory, stored as a link
            else:
                directories.append(root / name)
        for name in names:
            if root.resolve() / name not in skip:
                files.append(root / name)
    if not files:
        return None, 0, 0, 0

    size_before = sum(path.lstat().st_size for path in files)
    try:
        with tarfile.open(partial, "w:gz") as tar:
            for path in sorted(files):
                tar.add(path, arcname=str(path.relative_to(basedir)), recursive=False)
        with tarfile.open(partial, "r:gz") as tar:
            archived = {member.name for member in tar.getmembers()}
        expected = {str(path.relative_to(basedir)) for path in files}
        if archived != expected:
            raise RuntimeError("archive check failed: {} files expected, {} found".format(len(expected), len(archived)))
        partial.replace(archive)
    except Exception:
        partial.unlink(missing_ok=True)
        raise

    for path in files:
        path.unlink()
    for path in sorted(directories, key=lambda p: len(p.parts), reverse=True):
        try:
            path.rmdir()
        except OSError:
            pass  # not empty: holds a kept or protected file
    return archive, len(files), size_before, archive.stat().st_size
