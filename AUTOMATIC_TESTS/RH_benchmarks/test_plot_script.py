'''Guard: every figure the report includes must still be generatable.

Two edits in the course of this work silently deleted figure functions by
replacing a slice of the file that reached further than intended.  Nothing
complained, because the PDFs were already on disk and LaTeX kept building.
This test fails the moment that happens again.
'''
import ast
import pathlib
import re

HERE = pathlib.Path(__file__).resolve().parent
REPORT = HERE.parent.parent / 'DOCUMENTATION' / 'stella_RH_report' / 'stella_RH_report.tex'


def test_whether_every_registered_figure_function_exists():
    source = (HERE / 'plot_rh_report.py').read_text()
    defined = {n.name for n in ast.parse(source).body if isinstance(n, ast.FunctionDef)}
    registered = set(re.findall(r"':\s*'(figure_\w+)'", source))
    assert registered <= defined, f'registered but undefined: {sorted(registered - defined)}'


def test_whether_every_figure_the_report_includes_is_on_disk():
    if not REPORT.exists():
        return
    figures = set(re.findall(r'\\includegraphics\[[^\]]*\]\{figures/([^}]+)\}', REPORT.read_text()))
    missing = sorted(f for f in figures if not (REPORT.parent / 'figures' / f).exists())
    assert not missing, f'included by the report but not on disk: {missing}'


def test_whether_the_figures_directory_has_no_orphans():
    if not REPORT.exists():
        return
    included = set(re.findall(r'\\includegraphics\[[^\]]*\]\{figures/([^}]+)\}', REPORT.read_text()))
    on_disk = {p.name for p in (REPORT.parent / 'figures').glob('*.pdf')}
    orphans = sorted(on_disk - included)
    assert not orphans, f'on disk but no longer included: {orphans}'
